!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Main driver for L-BFGS optimizer
!> \par History
!>      01.2020 Space Group Symmetry introduced by Pierre-André Cazade [pcazade]
! **************************************************************************************************
MODULE cp_lbfgs_geo
   USE cp_lbfgs_optimizer_gopt,         ONLY: cp_lbfgs_opt_gopt_type,&
                                              cp_opt_gopt_create,&
                                              cp_opt_gopt_next,&
                                              cp_opt_gopt_release,&
                                              cp_opt_gopt_stop
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
   USE gopt_f_methods,                  ONLY: gopt_f_io_finalize,&
                                              print_geo_opt_header
   USE gopt_f_types,                    ONLY: gopt_f_type
   USE gopt_param_types,                ONLY: gopt_param_type
   USE input_constants,                 ONLY: default_ts_method_id
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE space_groups_types,              ONLY: spgr_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_lbfgs_geo'

   PUBLIC :: geoopt_lbfgs

   INTEGER, PARAMETER :: lbfgs_abnormal_exit_status = 5

! **************************************************************************************************
!> \brief L-BFGS backend: wraps the master/minion setulb state machine, which
!>        already reports its steps inside cp_opt_gopt_next (io_inside_next).
! **************************************************************************************************
   TYPE, EXTENDS(gopt_backend_type) :: lbfgs_opt_type
      PRIVATE
      TYPE(cp_lbfgs_opt_gopt_type), POINTER             :: optimizer => NULL()
   CONTAINS
      PROCEDURE, PUBLIC :: finalize => lbfgs_finalize
      PROCEDURE, PUBLIC :: init => lbfgs_init
      PROCEDURE, PUBLIC :: next => lbfgs_next
      PROCEDURE, PUBLIC :: on_stop => lbfgs_on_stop
      PROCEDURE, PUBLIC :: post_step => lbfgs_post_step
   END TYPE lbfgs_opt_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param force_env ...
!> \param gopt_param ...
!> \param globenv ...
!> \param geo_section ...
!> \param gopt_env ...
!> \param x0 ...
!> \par History
!>      08.2003 created [fawzi]
!>      01.2020 modified [pcazade]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE geoopt_lbfgs(force_env, gopt_param, globenv, geo_section, gopt_env, x0)
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(gopt_param_type), POINTER                     :: gopt_param
      TYPE(global_environment_type), POINTER             :: globenv
      TYPE(section_vals_type), POINTER                   :: geo_section
      TYPE(gopt_f_type), POINTER                         :: gopt_env
      REAL(KIND=dp), DIMENSION(:), POINTER               :: x0

      CHARACTER(len=*), PARAMETER                        :: routineN = 'geoopt_lbfgs'

      INTEGER                                            :: handle, iter_nr, its, output_unit
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(lbfgs_opt_type)                               :: lbfgs
      TYPE(spgr_type), POINTER                           :: spgr

      CALL timeset(routineN, handle)

      NULLIFY (logger, subsys, spgr)
      logger => cp_get_default_logger()
      spgr => gopt_env%spgr
      CPASSERT(ASSOCIATED(force_env))
      CPASSERT(ASSOCIATED(gopt_param))

      ! collecting subsys
      CALL force_env_get(force_env, subsys=subsys)

      ! Geometry optimization starts now
      CALL gopt_open_geolog(logger, geo_section, output_unit)
      CALL print_geo_opt_header(gopt_env, output_unit, "L-BFGS")

      ! finds space group
      CALL gopt_setup_space_group(spgr, subsys, geo_section, gopt_env, x0, output_unit)

      CALL lbfgs%init(force_env, gopt_param, globenv, geo_section, gopt_env, x0, output_unit, iter_nr)
      CALL gopt_optimize_loop(lbfgs, force_env, gopt_param, globenv, geo_section, gopt_env, x0, &
                              gopt_param%max_iter, iter_nr, output_unit, its)

      IF (lbfgs%optimizer%status == lbfgs_abnormal_exit_status) THEN
         IF (output_unit > 0) THEN
            WRITE (output_unit, FMT="(/,T2,A)") "*** L-BFGS ABNORMAL EXIT: OPTIMIZATION NOT CONVERGED ***"
            WRITE (output_unit, FMT="(T2,A)") "L-BFGS task: "//TRIM(lbfgs%optimizer%task)
         END IF
      ELSE
         CALL gopt_report_nc(gopt_env, output_unit, its, gopt_param%max_iter, lbfgs%conv)
      END IF

      ! show space_group
      CALL gopt_show_space_group(spgr, subsys, geo_section, gopt_env, output_unit)

      ! Write final output information, if converged
      CALL lbfgs%finalize(force_env, gopt_env, x0, lbfgs%conv, its, force_env%root_section, output_unit)

      CALL gopt_close_geolog(output_unit, logger, geo_section)

      CALL timestop(handle)

   END SUBROUTINE geoopt_lbfgs

! **************************************************************************************************
!> \brief Creates the setulb state machine and opens the iteration counters.
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
   RECURSIVE SUBROUTINE lbfgs_init(backend, force_env, gopt_param, globenv, geo_section, gopt_env, &
                                   x0, output_unit, iter_nr)
      CLASS(lbfgs_opt_type), INTENT(INOUT)               :: backend
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(gopt_param_type), POINTER                     :: gopt_param
      TYPE(global_environment_type), POINTER             :: globenv
      TYPE(section_vals_type), POINTER                   :: geo_section
      TYPE(gopt_f_type), POINTER                         :: gopt_env
      REAL(KIND=dp), DIMENSION(:), POINTER               :: x0
      INTEGER, INTENT(IN)                                :: output_unit
      INTEGER, INTENT(OUT)                               :: iter_nr

      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mp_para_env_type), POINTER                    :: para_env

      MARK_USED(globenv)
      MARK_USED(geo_section)
      MARK_USED(output_unit)

      NULLIFY (logger, para_env)
      logger => cp_get_default_logger()
      backend%io_inside_next = .TRUE.
      backend%conv = .FALSE.

      ! Stop if not implemented
      IF (gopt_env%type_id == default_ts_method_id) THEN
         CPABORT("L-BFGS method not yet working with DIMER")
      END IF

      CALL force_env_get(force_env, para_env=para_env)
      ALLOCATE (backend%optimizer)
      CALL cp_opt_gopt_create(backend%optimizer, para_env=para_env, obj_funct=gopt_env, &
                              x0=x0, m=gopt_param%max_h_rank, &
                              print_every=gopt_param%print_control, &
                              wanted_relative_f_delta=gopt_param%wanted_rel_f_error, &
                              wanted_projected_gradient=gopt_param%wanted_proj_gradient, &
                              max_f_per_iter=gopt_param%max_f_per_iter, &
                              trust_radius=gopt_param%trust_radius)
      CALL cp_iterate(logger%iter_info, increment=0, iter_nr_out=iter_nr)

   END SUBROUTINE lbfgs_init

! **************************************************************************************************
!> \brief Advances the state machine, which evaluates and reports internally.
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
   RECURSIVE FUNCTION lbfgs_next(backend, its, force_env, gopt_param, globenv, geo_section, &
                                 gopt_env, x0) RESULT(keep)
      CLASS(lbfgs_opt_type), INTENT(INOUT)               :: backend
      INTEGER, INTENT(IN)                                :: its
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(gopt_param_type), POINTER                     :: gopt_param
      TYPE(global_environment_type), POINTER             :: globenv
      TYPE(section_vals_type), POINTER                   :: geo_section
      TYPE(gopt_f_type), POINTER                         :: gopt_env
      REAL(KIND=dp), DIMENSION(:), POINTER               :: x0
      LOGICAL                                            :: keep

      MARK_USED(globenv)
      MARK_USED(its)
      MARK_USED(x0)

      ! Real optimization step..
      keep = cp_opt_gopt_next(backend%optimizer, geo_section=geo_section, &
                              force_env=force_env, gopt_param=gopt_param, &
                              converged=backend%conv, spgr=gopt_env%spgr)

   END FUNCTION lbfgs_next

! **************************************************************************************************
!> \brief No after-IO updates for the L-BFGS state machine.
!> \param backend ...
!> \param its ...
! **************************************************************************************************
   SUBROUTINE lbfgs_post_step(backend, its)
      CLASS(lbfgs_opt_type), INTENT(INOUT)               :: backend
      INTEGER, INTENT(IN)                                :: its

      MARK_USED(backend)
      MARK_USED(its)

   END SUBROUTINE lbfgs_post_step

! **************************************************************************************************
!> \brief Stops the setulb state machine on an external exit request.
!> \param backend ...
! **************************************************************************************************
   SUBROUTINE lbfgs_on_stop(backend)
      CLASS(lbfgs_opt_type), INTENT(INOUT)               :: backend

      CALL cp_opt_gopt_stop(backend%optimizer)

   END SUBROUTINE lbfgs_on_stop

! **************************************************************************************************
!> \brief Writes the final information and releases the state machine.
!> \param backend ...
!> \param force_env ...
!> \param gopt_env ...
!> \param x0 ...
!> \param conv ...
!> \param its ...
!> \param root_section ...
!> \param output_unit ...
! **************************************************************************************************
   RECURSIVE SUBROUTINE lbfgs_finalize(backend, force_env, gopt_env, x0, conv, its, root_section, &
                                       output_unit)
      CLASS(lbfgs_opt_type), INTENT(INOUT)               :: backend
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(gopt_f_type), POINTER                         :: gopt_env
      REAL(KIND=dp), DIMENSION(:), POINTER               :: x0
      LOGICAL, INTENT(IN)                                :: conv
      INTEGER, INTENT(IN)                                :: its
      TYPE(section_vals_type), POINTER                   :: root_section
      INTEGER, INTENT(IN)                                :: output_unit

      TYPE(cp_logger_type), POINTER                      :: logger

      MARK_USED(x0)

      NULLIFY (logger)
      logger => cp_get_default_logger()

      ! Write final output information, if converged
      CALL cp_iterate(logger%iter_info, last=.TRUE., increment=0)
      CALL gopt_f_io_finalize(gopt_env, force_env, backend%optimizer%x, conv, its, root_section, &
                              backend%optimizer%para_env, backend%optimizer%master, output_unit)

      CALL cp_opt_gopt_release(backend%optimizer)
      DEALLOCATE (backend%optimizer)

   END SUBROUTINE lbfgs_finalize

END MODULE cp_lbfgs_geo
