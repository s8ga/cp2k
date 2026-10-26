!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for Geometry optimization using  Conjugate Gradients
!> \author Teodoro Laino [teo]
!>      10.2005
! **************************************************************************************************
MODULE cg_optimizer

   USE cg_utils,                        ONLY: cg_linmin,&
                                              get_conjugate_direction
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_iterate
   USE cp_subsys_types,                 ONLY: cp_subsys_type
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
   USE global_types,                    ONLY: global_environment_type
   USE gopt_driver,                     ONLY: gopt_backend_type,&
                                              gopt_close_geolog,&
                                              gopt_open_geolog,&
                                              gopt_optimize_loop,&
                                              gopt_report_nc,&
                                              gopt_setup_space_group,&
                                              gopt_show_space_group
   USE gopt_f_methods,                  ONLY: cp_eval_at,&
                                              gopt_f_io_finalize,&
                                              gopt_f_io_init,&
                                              print_geo_opt_header
   USE gopt_f_types,                    ONLY: gopt_f_type
   USE gopt_param_types,                ONLY: gopt_param_type
   USE input_constants,                 ONLY: default_cell_method_id,&
                                              default_minimization_method_id,&
                                              default_ts_method_id
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_walltime
   USE space_groups,                    ONLY: spgr_apply_rotations_coord,&
                                              spgr_apply_rotations_force
   USE space_groups_types,              ONLY: spgr_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: geoopt_cg
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cg_optimizer'

! **************************************************************************************************
!> \brief CG backend state driven by gopt_optimize_loop: search direction and
!>        gradient history. The shared gradient g aliases xi (the gradient at
!>        the current point), gf keeps the forces.
! **************************************************************************************************
   TYPE, EXTENDS(gopt_backend_type) :: cg_opt_type
      PRIVATE
      INTEGER                                            :: iter_nr = 0, max_steep_steps = 0
      LOGICAL                                            :: evaluate_before_io = .TRUE., Fletcher_Reeves = .FALSE.
      REAL(KIND=dp)                                      :: res_lim = 0.0_dp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: xold
      REAL(KIND=dp), DIMENSION(:), POINTER               :: gf => NULL(), h => NULL(), xi => NULL(), x0 => NULL()
      TYPE(gopt_f_type), POINTER                         :: gopt_env => NULL()
   CONTAINS
      PROCEDURE, PUBLIC :: finalize => cg_finalize
      PROCEDURE, PUBLIC :: init => cg_init
      PROCEDURE, PUBLIC :: next => cg_next
      PROCEDURE, PUBLIC :: on_stop => cg_on_stop
      PROCEDURE, PUBLIC :: post_step => cg_post_step
   END TYPE cg_opt_type

CONTAINS

! **************************************************************************************************
!> \brief Driver for conjugate gradient optimization technique
!> \param force_env ...
!> \param gopt_param ...
!> \param globenv ...
!> \param geo_section ...
!> \param gopt_env ...
!> \param x0 ...
!> \param do_update ...
!> \par History
!>      10.2005 created [tlaino]
!> \author Teodoro Laino
! **************************************************************************************************
   RECURSIVE SUBROUTINE geoopt_cg(force_env, gopt_param, globenv, geo_section, &
                                  gopt_env, x0, do_update)

      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(gopt_param_type), POINTER                     :: gopt_param
      TYPE(global_environment_type), POINTER             :: globenv
      TYPE(section_vals_type), POINTER                   :: geo_section
      TYPE(gopt_f_type), POINTER                         :: gopt_env
      REAL(KIND=dp), DIMENSION(:), POINTER               :: x0
      LOGICAL, INTENT(OUT), OPTIONAL                     :: do_update

      CHARACTER(len=*), PARAMETER                        :: routineN = 'geoopt_cg'

      INTEGER                                            :: handle, iter_nr, its, output_unit
      LOGICAL                                            :: my_do_update
      TYPE(cg_opt_type)                                  :: cg
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(spgr_type), POINTER                           :: spgr

      CALL timeset(routineN, handle)

      NULLIFY (logger, subsys, spgr)
      logger => cp_get_default_logger()
      spgr => gopt_env%spgr

      CALL gopt_open_geolog(logger, geo_section, output_unit)
      CALL print_geo_opt_header(gopt_env, output_unit, "CONJUGATE GRADIENTS")

      ! find space_group
      CALL force_env_get(force_env, subsys=subsys)
      IF (spgr%keep_space_group) THEN
         SELECT CASE (gopt_env%type_id)
         CASE (default_minimization_method_id, default_ts_method_id, default_cell_method_id)
            CALL gopt_setup_space_group(spgr, subsys, geo_section, gopt_env, x0, output_unit)
         CASE DEFAULT
            spgr%keep_space_group = .FALSE.
         END SELECT
      END IF

      CALL cg%init(force_env, gopt_param, globenv, geo_section, gopt_env, x0, output_unit, iter_nr)
      CALL gopt_optimize_loop(cg, force_env, gopt_param, globenv, geo_section, gopt_env, x0, &
                              gopt_param%max_iter, iter_nr, output_unit, its)

      my_do_update = cg%conv
      CALL gopt_report_nc(gopt_env, output_unit, its, gopt_param%max_iter, cg%conv)
      CALL cg%finalize(force_env, gopt_env, x0, cg%conv, its, force_env%root_section, output_unit)

      ! show space_group
      CALL gopt_show_space_group(spgr, subsys, geo_section, gopt_env, output_unit)

      CALL gopt_close_geolog(output_unit, logger, geo_section)
      IF (PRESENT(do_update)) do_update = my_do_update

      CALL timestop(handle)

   END SUBROUTINE geoopt_cg

! **************************************************************************************************
!> \brief Sets up the search direction arrays and evaluates the start point
!>        with the CG force-consistency setting.
!> \param backend ...
!> \param force_env ...
!> \param gopt_param ...
!> \param globenv ...
!> \param geo_section ...
!> \param gopt_env ...
!> \param x0 ...
!> \param output_unit ...
!> \param iter_nr ...
! **************************************************************************************************
   RECURSIVE SUBROUTINE cg_init(backend, force_env, gopt_param, globenv, geo_section, gopt_env, &
                                x0, output_unit, iter_nr)
      CLASS(cg_opt_type), INTENT(INOUT)                  :: backend
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(gopt_param_type), POINTER                     :: gopt_param
      TYPE(global_environment_type), POINTER             :: globenv
      TYPE(section_vals_type), POINTER                   :: geo_section
      TYPE(gopt_f_type), POINTER                         :: gopt_env
      REAL(KIND=dp), DIMENSION(:), POINTER               :: x0
      INTEGER, INTENT(IN)                                :: output_unit
      INTEGER, INTENT(OUT)                               :: iter_nr

      LOGICAL                                            :: save_consistent_energy_force
      REAL(KIND=dp)                                      :: t_diff, t_now, t_old
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(spgr_type), POINTER                           :: spgr

      MARK_USED(force_env)
      MARK_USED(globenv)
      MARK_USED(geo_section)

      NULLIFY (logger, spgr)
      logger => cp_get_default_logger()
      spgr => gopt_env%spgr
      t_old = m_walltime()

      backend%conv = .FALSE.
      backend%max_steep_steps = gopt_param%max_steep_steps
      backend%Fletcher_Reeves = gopt_param%Fletcher_Reeves
      ! Dimer evaluations update the rotational state, so retain their original timing.
      backend%evaluate_before_io = gopt_env%type_id /= default_ts_method_id
      backend%res_lim = gopt_param%restart_limit
      backend%output_unit = output_unit
      backend%gopt_env => gopt_env
      backend%x0 => x0

      ALLOCATE (backend%xi(SIZE(x0)))
      ALLOCATE (backend%gf(SIZE(x0)))
      ALLOCATE (backend%h(SIZE(x0)))
      ALLOCATE (backend%xold(SIZE(x0)))
      ALLOCATE (backend%dr(SIZE(x0)))
      backend%dr(:) = 0.0_dp
      backend%g => backend%xi

      ! applies rotation matrices to coordinates
      IF (spgr%keep_space_group) THEN
         CALL spgr_apply_rotations_coord(spgr, x0)
      END IF

      ! Evaluate energy and forces at the first step
      ![NB] consistent energies and forces not required for CG, but some line minimizers might set it
      save_consistent_energy_force = gopt_env%require_consistent_energy_force
      gopt_env%require_consistent_energy_force = .FALSE.

      CALL cp_eval_at(gopt_env, x0, backend%etot, backend%xi, master=gopt_env%force_env%para_env%mepos, &
                      para_env=gopt_env%force_env%para_env)

      gopt_env%require_consistent_energy_force = save_consistent_energy_force

      ! Symmetrize coordinates and forces
      IF (spgr%keep_space_group) THEN
         CALL spgr_apply_rotations_coord(spgr, x0)
         CALL spgr_apply_rotations_force(spgr, backend%xi)
      END IF

      backend%gf = -backend%xi
      backend%h = backend%gf
      backend%xi = backend%h
      backend%emin = HUGE(0.0_dp)
      CALL cp_iterate(logger%iter_info, increment=0, iter_nr_out=iter_nr)
      backend%iter_nr = iter_nr
      backend%wildcard = "   SD"
      t_now = m_walltime()
      t_diff = t_now - t_old
      CALL gopt_f_io_init(gopt_env, output_unit, backend%etot, backend%wildcard, used_time=t_diff, its=iter_nr)
      backend%eold = backend%etot

   END SUBROUTINE cg_init

! **************************************************************************************************
!> \brief Runs the line minimization along the current search direction and
!>        evaluates the accepted point.
!> \param backend ...
!> \param its ...
!> \param force_env ...
!> \param gopt_param ...
!> \param globenv ...
!> \param geo_section ...
!> \param gopt_env ...
!> \param x0 ...
!> \return ...
! **************************************************************************************************
   RECURSIVE FUNCTION cg_next(backend, its, force_env, gopt_param, globenv, geo_section, &
                              gopt_env, x0) RESULT(keep)
      CLASS(cg_opt_type), INTENT(INOUT)                  :: backend
      INTEGER, INTENT(IN)                                :: its
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(gopt_param_type), POINTER                     :: gopt_param
      TYPE(global_environment_type), POINTER             :: globenv
      TYPE(section_vals_type), POINTER                   :: geo_section
      TYPE(gopt_f_type), POINTER                         :: gopt_env
      REAL(KIND=dp), DIMENSION(:), POINTER               :: x0
      LOGICAL                                            :: keep

      TYPE(spgr_type), POINTER                           :: spgr

      MARK_USED(force_env)
      MARK_USED(geo_section)
      MARK_USED(its)

      NULLIFY (spgr)
      spgr => gopt_env%spgr

      ! Symmetrize coordinates and forces
      IF (spgr%keep_space_group) THEN
         CALL spgr_apply_rotations_coord(spgr, x0)
         CALL spgr_apply_rotations_force(spgr, backend%gf)
         CALL spgr_apply_rotations_force(spgr, backend%xi)
      END IF

      backend%xold(:) = x0

      ! Line minimization
      CALL cg_linmin(gopt_env, x0, backend%xi, backend%gf, backend%etot, backend%output_unit, gopt_param, globenv)

      ! Applies rotation matrices to coordinates
      IF (spgr%keep_space_group) THEN
         CALL spgr_apply_rotations_coord(spgr, x0)
      END IF

      IF (backend%evaluate_before_io) THEN
         ! Keep the accepted point and the physical state consistent even on external stop.
         CALL cg_eval_current(gopt_env, x0, backend%etot, backend%xi, spgr)
      END IF

      backend%dr(:) = x0(:) - backend%xold(:)

      keep = .TRUE.

   END FUNCTION cg_next

! **************************************************************************************************
!> \brief Prepares the conjugate (or reset steepest-descent) direction for the
!>        next line minimization.
!> \param backend ...
!> \param its ...
! **************************************************************************************************
   RECURSIVE SUBROUTINE cg_post_step(backend, its)
      CLASS(cg_opt_type), INTENT(INOUT)                  :: backend
      INTEGER, INTENT(IN)                                :: its

      TYPE(spgr_type), POINTER                           :: spgr

      NULLIFY (spgr)
      spgr => backend%gopt_env%spgr

      IF (.NOT. backend%evaluate_before_io) THEN
         CALL cg_eval_current(backend%gopt_env, backend%x0, backend%etot, backend%xi, spgr)
      END IF

      ! Get Conjugate Directions:  updates the searching direction (h)
      backend%wildcard = "   CG"
      CALL get_conjugate_direction(backend%gopt_env, backend%Fletcher_Reeves, backend%gf, backend%xi, backend%h)

      ! Symmetrize coordinates and forces
      IF (spgr%keep_space_group) THEN
         CALL spgr_apply_rotations_force(spgr, backend%gf)
         CALL spgr_apply_rotations_force(spgr, backend%h)
      END IF

      ! Reset Condition or Steepest Descent Requested
      ! ABS(DOT_PRODUCT(g, h))/SQRT((DOT_PRODUCT(g, g)*DOT_PRODUCT(h, h))) > res_lim ...
      IF ((DOT_PRODUCT(backend%gf, backend%h)*DOT_PRODUCT(backend%gf, backend%h)) > &
          (backend%res_lim*backend%res_lim*DOT_PRODUCT(backend%gf, backend%gf)*DOT_PRODUCT(backend%h, backend%h)) &
          .OR. its + 1 <= backend%max_steep_steps) THEN
         ! Steepest Descent
         backend%wildcard = "   SD"
         backend%h = -backend%xi
      END IF
      backend%gf = -backend%xi
      backend%xi = backend%h

   END SUBROUTINE cg_post_step

! **************************************************************************************************
!> \brief Nothing to stop for CG on an external exit request.
!> \param backend ...
! **************************************************************************************************
   SUBROUTINE cg_on_stop(backend)
      CLASS(cg_opt_type), INTENT(INOUT)                  :: backend

      MARK_USED(backend)

   END SUBROUTINE cg_on_stop

! **************************************************************************************************
!> \brief Writes the final information and releases the direction arrays.
!> \param backend ...
!> \param force_env ...
!> \param gopt_env ...
!> \param x0 ...
!> \param conv ...
!> \param its ...
!> \param root_section ...
!> \param output_unit ...
! **************************************************************************************************
   RECURSIVE SUBROUTINE cg_finalize(backend, force_env, gopt_env, x0, conv, its, root_section, &
                                    output_unit)
      CLASS(cg_opt_type), INTENT(INOUT)                  :: backend
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(gopt_f_type), POINTER                         :: gopt_env
      REAL(KIND=dp), DIMENSION(:), POINTER               :: x0
      LOGICAL, INTENT(IN)                                :: conv
      INTEGER, INTENT(IN)                                :: its
      TYPE(section_vals_type), POINTER                   :: root_section
      INTEGER, INTENT(IN)                                :: output_unit

      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (logger)
      logger => cp_get_default_logger()

      ! Write final particle information and restart, if converged
      CALL cp_iterate(logger%iter_info, last=.TRUE., increment=0)
      CALL gopt_f_io_finalize(gopt_env, force_env, x0, conv, its, root_section, &
                              gopt_env%force_env%para_env, gopt_env%force_env%para_env%mepos, output_unit)

      NULLIFY (backend%g)
      DEALLOCATE (backend%gf)
      DEALLOCATE (backend%h)
      DEALLOCATE (backend%xi)
      DEALLOCATE (backend%xold)
      DEALLOCATE (backend%dr)

   END SUBROUTINE cg_finalize

! **************************************************************************************************
!> \brief Evaluate the current CG point with the force-consistency setting used by CG.
!> \param gopt_env Optimization environment
!> \param x Current coordinates, possibly projected by cp_eval_at
!> \param energy Energy at the evaluated coordinates
!> \param gradient Gradient at the evaluated coordinates
!> \param spgr Space-group data
! **************************************************************************************************
   SUBROUTINE cg_eval_current(gopt_env, x, energy, gradient, spgr)
      TYPE(gopt_f_type), POINTER                         :: gopt_env
      REAL(KIND=dp), DIMENSION(:), POINTER               :: x
      REAL(KIND=dp), INTENT(OUT)                         :: energy
      REAL(KIND=dp), DIMENSION(:), POINTER               :: gradient
      TYPE(spgr_type), POINTER                           :: spgr

      LOGICAL                                            :: save_consistent_energy_force

      save_consistent_energy_force = gopt_env%require_consistent_energy_force
      gopt_env%require_consistent_energy_force = .FALSE.
      CALL cp_eval_at(gopt_env, x, energy, gradient, master=gopt_env%force_env%para_env%mepos, &
                      para_env=gopt_env%force_env%para_env)
      gopt_env%require_consistent_energy_force = save_consistent_energy_force

      IF (spgr%keep_space_group) CALL spgr_apply_rotations_force(spgr, gradient)
   END SUBROUTINE cg_eval_current

END MODULE cg_optimizer
