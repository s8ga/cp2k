!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Shared scaffolding of the geometry-optimizer drivers (BFGS, CG, L-BFGS):
!>        .geoLog plumbing, space-group setup and epilogue, the per-step iteration
!>        header, and the not-converged report. All pieces are lifted verbatim from
!>        the three drivers. Their call order inside each driver is unchanged.
! **************************************************************************************************
MODULE gopt_driver
   USE cp_external_control,             ONLY: external_control
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_iterate,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_subsys_types,                 ONLY: cp_subsys_type
   USE force_env_types,                 ONLY: force_env_type
   USE global_types,                    ONLY: global_environment_type
   USE gopt_f_methods,                  ONLY: gopt_f_ii,&
                                              gopt_f_io,&
                                              print_geo_opt_nc
   USE gopt_f_types,                    ONLY: gopt_f_type
   USE gopt_param_types,                ONLY: gopt_param_type
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_set
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_walltime
   USE space_groups,                    ONLY: identify_space_group,&
                                              print_spgr,&
                                              spgr_apply_rotations_coord
   USE space_groups_types,              ONLY: spgr_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gopt_driver'

   PUBLIC :: gopt_backend_type, gopt_close_geolog, gopt_io_plain, gopt_io_pred, gopt_io_trust, &
             gopt_open_geolog, gopt_optimize_loop, gopt_report_nc, gopt_setup_space_group, &
             gopt_show_space_group

   ! Observables gopt_optimize_loop hands to gopt_f_io: plain (CG), with the
   ! trust-radius quantities (BFGS inverse), or with the energy prediction on
   ! top (BFGS direct).
   INTEGER, PARAMETER :: gopt_io_plain = 0, gopt_io_pred = 2, gopt_io_trust = 1

! **************************************************************************************************
!> \brief Contract every geometry-optimizer backend fulfils for the shared loop:
!>        init sets up state and evaluates the start point, next advances to the
!>        next accepted point, post_step performs after-IO updates, finalize
!>        writes the final information and tears down.
!>        The loop owns the iteration header, the external-exit check, the
!>        gopt_f_io report and the eold/emin bookkeeping. Backends that
!>        print inside next (the L-BFGS state machine) set io_inside_next.
! **************************************************************************************************
   TYPE, ABSTRACT :: gopt_backend_type
      ! Per-step observables produced by next and consumed by gopt_f_io.
      REAL(KIND=dp)                                      :: ediff = 0.0_dp, emin = 0.0_dp, &
                                                            eold = 0.0_dp, etot = 0.0_dp, &
                                                            pred = 0.0_dp, rad = 0.0_dp, &
                                                            rat = 0.0_dp, step = 0.0_dp
      REAL(KIND=dp), DIMENSION(:), POINTER              :: dr => NULL(), g => NULL()
      CHARACTER(LEN=5)                                   :: wildcard = ""
      LOGICAL                                            :: conv = .FALSE., io_inside_next = .FALSE.
      INTEGER                                            :: io_variant = gopt_io_plain, output_unit = -1
   CONTAINS
      PROCEDURE(gopt_backend_init_fn), DEFERRED         :: init
      PROCEDURE(gopt_backend_next_fn), DEFERRED         :: next
      PROCEDURE(gopt_backend_post_fn), DEFERRED         :: post_step
      PROCEDURE(gopt_backend_stop_fn), DEFERRED         :: on_stop
      PROCEDURE(gopt_backend_finalize_fn), DEFERRED     :: finalize
   END TYPE gopt_backend_type

   ABSTRACT INTERFACE

! **************************************************************************************************
!> \brief Prepares the backend and evaluates the start point (gopt_f_io_init).
!> \param backend ...
!> \param force_env ...
!> \param gopt_param ...
!> \param globenv ...
!> \param geo_section ...
!> \param gopt_env ...
!> \param x0 ...
!> \param output_unit ...
!> \param iter_nr iteration counter of the start point
! **************************************************************************************************
      SUBROUTINE gopt_backend_init_fn(backend, force_env, gopt_param, globenv, geo_section, &
                                      gopt_env, x0, output_unit, iter_nr)
         IMPORT :: dp, force_env_type, global_environment_type, gopt_backend_type, gopt_f_type, &
            gopt_param_type, section_vals_type
         CLASS(gopt_backend_type), INTENT(INOUT)        :: backend
         TYPE(force_env_type), POINTER                  :: force_env
         TYPE(gopt_param_type), POINTER                 :: gopt_param
         TYPE(global_environment_type), POINTER         :: globenv
         TYPE(section_vals_type), POINTER               :: geo_section
         TYPE(gopt_f_type), POINTER                     :: gopt_env
         REAL(KIND=dp), DIMENSION(:), POINTER           :: x0
         INTEGER, INTENT(IN)                            :: output_unit
         INTEGER, INTENT(OUT)                           :: iter_nr
      END SUBROUTINE gopt_backend_init_fn

! **************************************************************************************************
!> \brief Advances to the next accepted point, filling the observables of the
!>        base type. Returns .FALSE. when the backend ended the optimization
!>        on its own (e.g. the L-BFGS state machine converged or failed).
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
      FUNCTION gopt_backend_next_fn(backend, its, force_env, gopt_param, globenv, geo_section, &
                                    gopt_env, x0) RESULT(keep)
         IMPORT :: dp, force_env_type, global_environment_type, gopt_backend_type, gopt_f_type, &
            gopt_param_type, section_vals_type
         CLASS(gopt_backend_type), INTENT(INOUT)        :: backend
         INTEGER, INTENT(IN)                            :: its
         TYPE(force_env_type), POINTER                  :: force_env
         TYPE(gopt_param_type), POINTER                 :: gopt_param
         TYPE(global_environment_type), POINTER         :: globenv
         TYPE(section_vals_type), POINTER               :: geo_section
         TYPE(gopt_f_type), POINTER                     :: gopt_env
         REAL(KIND=dp), DIMENSION(:), POINTER           :: x0
         LOGICAL                                        :: keep
      END FUNCTION gopt_backend_next_fn

! **************************************************************************************************
!> \brief After-IO updates while the optimization continues (default: no-op).
!> \param backend ...
!> \param its ...
! **************************************************************************************************
      SUBROUTINE gopt_backend_post_fn(backend, its)
         IMPORT :: gopt_backend_type
         CLASS(gopt_backend_type), INTENT(INOUT)        :: backend
         INTEGER, INTENT(IN)                            :: its
      END SUBROUTINE gopt_backend_post_fn

! **************************************************************************************************
!> \brief Called when the external-exit check stops the optimization. Most
!>        backends have no state to stop. The L-BFGS state machine does.
!> \param backend ...
! **************************************************************************************************
      SUBROUTINE gopt_backend_stop_fn(backend)
         IMPORT :: gopt_backend_type
         CLASS(gopt_backend_type), INTENT(INOUT)        :: backend
      END SUBROUTINE gopt_backend_stop_fn

! **************************************************************************************************
!> \brief Writes the final information (cp_iterate last, backend restart data,
!>        gopt_f_io_finalize) and releases backend resources.
!> \param backend ...
!> \param force_env ...
!> \param gopt_env ...
!> \param x0 ...
!> \param conv ...
!> \param its ...
!> \param root_section ...
!> \param output_unit ...
! **************************************************************************************************
      SUBROUTINE gopt_backend_finalize_fn(backend, force_env, gopt_env, x0, conv, its, &
                                          root_section, output_unit)
         IMPORT :: dp, force_env_type, gopt_backend_type, gopt_f_type, section_vals_type
         CLASS(gopt_backend_type), INTENT(INOUT)        :: backend
         TYPE(force_env_type), POINTER                  :: force_env
         TYPE(gopt_f_type), POINTER                     :: gopt_env
         REAL(KIND=dp), DIMENSION(:), POINTER           :: x0
         LOGICAL, INTENT(IN)                            :: conv
         INTEGER, INTENT(IN)                            :: its
         TYPE(section_vals_type), POINTER               :: root_section
         INTEGER, INTENT(IN)                            :: output_unit
      END SUBROUTINE gopt_backend_finalize_fn
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Opens the geometry-optimization log unit (PRINT%PROGRAM_RUN_INFO, .geoLog).
!> \param logger ...
!> \param geo_section ...
!> \param output_unit ...
! **************************************************************************************************
   SUBROUTINE gopt_open_geolog(logger, geo_section, output_unit)
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: geo_section
      INTEGER, INTENT(OUT)                               :: output_unit

      output_unit = cp_print_key_unit_nr(logger, geo_section, "PRINT%PROGRAM_RUN_INFO", &
                                         extension=".geoLog")

   END SUBROUTINE gopt_open_geolog

! **************************************************************************************************
!> \brief Closes the geometry-optimization log unit.
!> \param output_unit ...
!> \param logger ...
!> \param geo_section ...
! **************************************************************************************************
   SUBROUTINE gopt_close_geolog(output_unit, logger, geo_section)
      INTEGER, INTENT(INOUT)                             :: output_unit
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: geo_section

      CALL cp_print_key_finished_output(output_unit, logger, geo_section, &
                                        "PRINT%PROGRAM_RUN_INFO")

   END SUBROUTINE gopt_close_geolog

! **************************************************************************************************
!> \brief Identifies the space group and symmetrizes the start coordinates, so that
!>        the optimization keeps it. Callers must invoke this before the first
!>        evaluation. BFGS does so only after its Hessian restart read.
!>        That read may re-express x0 in the stored reference frame.
!> \param spgr ...
!> \param subsys ...
!> \param geo_section ...
!> \param gopt_env ...
!> \param x0 ...
!> \param output_unit ...
! **************************************************************************************************
   SUBROUTINE gopt_setup_space_group(spgr, subsys, geo_section, gopt_env, x0, output_unit)
      TYPE(spgr_type), POINTER                           :: spgr
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(section_vals_type), POINTER                   :: geo_section
      TYPE(gopt_f_type), POINTER                         :: gopt_env
      REAL(KIND=dp), DIMENSION(:), POINTER               :: x0
      INTEGER, INTENT(IN)                                :: output_unit

      IF (spgr%keep_space_group) THEN
         CALL identify_space_group(subsys, geo_section, gopt_env, output_unit)
         CALL spgr_apply_rotations_coord(spgr, x0)
         CALL print_spgr(spgr)
      END IF

   END SUBROUTINE gopt_setup_space_group

! **************************************************************************************************
!> \brief Re-identifies and prints the space group at the end of the optimization.
!> \param spgr ...
!> \param subsys ...
!> \param geo_section ...
!> \param gopt_env ...
!> \param output_unit ...
! **************************************************************************************************
   SUBROUTINE gopt_show_space_group(spgr, subsys, geo_section, gopt_env, output_unit)
      TYPE(spgr_type), POINTER                           :: spgr
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(section_vals_type), POINTER                   :: geo_section
      TYPE(gopt_f_type), POINTER                         :: gopt_env
      INTEGER, INTENT(IN)                                :: output_unit

      IF (spgr%show_space_group) THEN
         CALL identify_space_group(subsys, geo_section, gopt_env, output_unit)
         CALL print_spgr(spgr)
      END IF

   END SUBROUTINE gopt_show_space_group

! **************************************************************************************************
!> \brief Prints the not-converged report when the iteration budget was exhausted.
!> \param gopt_env ...
!> \param output_unit ...
!> \param its ...
!> \param maxiter ...
!> \param conv ...
! **************************************************************************************************
   SUBROUTINE gopt_report_nc(gopt_env, output_unit, its, maxiter, conv)
      TYPE(gopt_f_type), POINTER                         :: gopt_env
      INTEGER, INTENT(IN)                                :: output_unit, its, maxiter
      LOGICAL, INTENT(IN)                                :: conv

      IF (its == maxiter .AND. (.NOT. conv)) THEN
         CALL print_geo_opt_nc(gopt_env, output_unit)
      END IF

   END SUBROUTINE gopt_report_nc

! **************************************************************************************************
!> \brief Shared geometry-optimization loop driving any gopt_backend_type:
!>        iteration header, backend advance, external-exit check, the gopt_f_io
!>        report (skipped for backends that print inside next), exit conditions
!>        and the eold/emin bookkeeping. its_final reports the last iteration.
!> \param backend ...
!> \param force_env ...
!> \param gopt_param ...
!> \param globenv ...
!> \param geo_section ...
!> \param gopt_env ...
!> \param x0 ...
!> \param maxiter ...
!> \param iter_nr ...
!> \param output_unit ...
!> \param its_final ...
! **************************************************************************************************
   SUBROUTINE gopt_optimize_loop(backend, force_env, gopt_param, globenv, geo_section, gopt_env, &
                                 x0, maxiter, iter_nr, output_unit, its_final)
      CLASS(gopt_backend_type), INTENT(INOUT)           :: backend
      TYPE(force_env_type), POINTER                     :: force_env
      TYPE(gopt_param_type), POINTER                    :: gopt_param
      TYPE(global_environment_type), POINTER            :: globenv
      TYPE(section_vals_type), POINTER                  :: geo_section
      TYPE(gopt_f_type), POINTER                        :: gopt_env
      REAL(KIND=dp), DIMENSION(:), POINTER              :: x0
      INTEGER, INTENT(IN)                               :: maxiter, iter_nr, output_unit
      INTEGER, INTENT(OUT)                              :: its_final

      CHARACTER(len=*), PARAMETER                       :: routineN = 'gopt_optimize_loop'

      INTEGER                                           :: handle, its, ndf
      LOGICAL                                           :: should_stop
      REAL(KIND=dp)                                     :: t_diff, t_now, t_old
      TYPE(cp_logger_type), POINTER                     :: logger
      TYPE(section_vals_type), POINTER                  :: root_section

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()
      root_section => force_env%root_section
      ndf = SIZE(x0)
      t_old = m_walltime()

      DO its = iter_nr + 1, maxiter
         CALL cp_iterate(logger%iter_info, last=(its == maxiter))
         CALL section_vals_val_set(geo_section, "STEP_START_VAL", i_val=its)
         CALL gopt_f_ii(its, output_unit)

         IF (.NOT. backend%next(its, force_env, gopt_param, globenv, geo_section, gopt_env, x0)) EXIT

         ! check for an external exit command
         CALL external_control(should_stop, "GEO", globenv=globenv)
         IF (should_stop) THEN
            CALL backend%on_stop()
            EXIT
         END IF

         IF (backend%io_inside_next) THEN
            ! gopt_f_io ran inside next already. The backend itself
            ! reports convergence exits.
            IF (its == maxiter) EXIT
         ELSE
            t_now = m_walltime()
            t_diff = t_now - t_old
            t_old = t_now
            SELECT CASE (backend%io_variant)
            CASE (gopt_io_pred)
               CALL gopt_f_io(gopt_env, force_env, root_section, its, backend%etot, output_unit, &
                              backend%eold, backend%emin, backend%wildcard, gopt_param, ndf, backend%dr, &
                              backend%g, backend%conv, backend%pred, backend%rat, backend%step, &
                              backend%rad, used_time=t_diff)
            CASE (gopt_io_trust)
               CALL gopt_f_io(gopt_env, force_env, root_section, its, backend%etot, output_unit, &
                              backend%eold, backend%emin, backend%wildcard, gopt_param, ndf, backend%dr, &
                              backend%g, backend%conv, rat=backend%rat, step=backend%step, &
                              rad=backend%rad, used_time=t_diff)
            CASE DEFAULT
               CALL gopt_f_io(gopt_env, force_env, root_section, its, backend%etot, output_unit, &
                              backend%eold, backend%emin, backend%wildcard, gopt_param, ndf, backend%dr, &
                              backend%g, backend%conv, used_time=t_diff)
            END SELECT

            IF (backend%conv .OR. its == maxiter) EXIT
            backend%eold = backend%etot
            IF (backend%etot < backend%emin) backend%emin = backend%etot
            CALL backend%post_step(its)
         END IF
      END DO
      its_final = its

      CALL timestop(handle)

   END SUBROUTINE gopt_optimize_loop

END MODULE gopt_driver
