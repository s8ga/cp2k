!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Contains type used for a Simulation Cell Optimization
!> \par History
!>      none
!> \author Teodoro Laino - created [tlaino] - 03.2008 - Zurich University
! **************************************************************************************************
MODULE cell_opt_types

   USE cell_manifold_utils,             ONLY: &
        cell_manifold_init, cell_manifold_type, comp_bx, comp_cx, comp_cy, gfix_alpha, gfix_beta, &
        gfix_gamma, gfix_v, info_ok, n_constraint_bits, n_geom, n_slot
   USE cell_methods,                    ONLY: canonicalize_cell_matrix,&
                                              cell_create
   USE cell_opt_utils,                  ONLY: cell_opt_constraint_label,&
                                              parse_cell_constraint,&
                                              read_external_press_tensor,&
                                              symmetry_constraint_rows
   USE cell_types,                      ONLY: cell_clone,&
                                              cell_release,&
                                              cell_sym_none,&
                                              cell_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
   USE input_constants,                 ONLY: fix_none
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE string_utilities,                ONLY: uppercase
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cell_opt_types'

   PUBLIC :: cell_opt_env_type, &
             cell_opt_env_create, &
             cell_opt_env_release, &
             cell_opt_setup_manifold, &
             cell_opt_fold_symmetry_rows

! **************************************************************************************************
!> \brief Type containing all informations abour the simulation cell optimization
!> \par History
!>      none
!> \author Teodoro Laino - created [tlaino] - 03.2008 - Zurich University
! **************************************************************************************************
   TYPE cell_opt_env_type
      ! Simulation cell optimization parameters
      ! constraint_id is a mask of fixed cell components and geometric
      ! invariants (comp_* / gfix_*): a set bit means fixed, 0 (= fix_none)
      ! leaves all degrees of freedom relaxable
      INTEGER                                    :: constraint_id = fix_none
      LOGICAL                                    :: keep_angles = .FALSE., &
                                                    keep_symmetry = .FALSE., &
                                                    keep_volume = .FALSE.
      REAL(KIND=dp)                              :: pres_ext = 0.0_dp, pres_int = 0.0_dp, pres_tol = 0.0_dp, pres_constr = 0.0_dp
      REAL(KIND=dp), DIMENSION(3, 3)             :: mtrx = 0.0_dp
      ! Fixed rotations between the user's Cartesian frame and the canonical
      ! frame used for the six independent cell optimization variables.
      REAL(KIND=dp), DIMENSION(3, 3)             :: input_to_opt = 0.0_dp, &
                                                    opt_to_input = 0.0_dp
      TYPE(cell_type), POINTER                   :: opt_cell => NULL()
      TYPE(cell_manifold_type)                   :: manifold
   END TYPE cell_opt_env_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param cell_env ...
!> \param force_env ...
!> \param geo_section ...
!> \par History
!>      none
!> \author Teodoro Laino - created [tlaino] - 03.2008 - Zurich University
! **************************************************************************************************
   SUBROUTINE cell_opt_env_create(cell_env, force_env, geo_section)
      TYPE(cell_opt_env_type), INTENT(OUT)               :: cell_env
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(section_vals_type), POINTER                   :: geo_section

      CHARACTER(LEN=default_string_length)               :: label
      CHARACTER(LEN=default_string_length), &
         DIMENSION(:), POINTER                           :: canonical_tokens, constraint_tokens
      INTEGER                                            :: output_unit
      LOGICAL                                            :: constraint_explicit
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (cell_env%opt_cell, cell, canonical_tokens, constraint_tokens)
      CALL force_env_get(force_env, cell=cell)
      CALL cell_create(cell_env%opt_cell)
      CALL cell_clone(cell, cell_env%opt_cell, tag="CELL_OPT_CANONICAL")
      CALL canonicalize_cell_matrix(cell_env%opt_cell)
      ! Keep the physical cell in the input orientation while the optimizer
      ! works with the canonical cell and its six non-redundant components.
      cell_env%input_to_opt = MATMUL(cell_env%opt_cell%hmat, cell%h_inv)
      cell_env%opt_to_input = MATMUL(cell%hmat, cell_env%opt_cell%h_inv)
      CALL section_vals_val_get(geo_section, "KEEP_VOLUME", l_val=cell_env%keep_volume)
      CALL section_vals_val_get(geo_section, "KEEP_ANGLES", l_val=cell_env%keep_angles)
      CALL section_vals_val_get(geo_section, "KEEP_SYMMETRY", l_val=cell_env%keep_symmetry)
      CALL section_vals_val_get(geo_section, "PRESSURE_TOLERANCE", r_val=cell_env%pres_tol)
      CALL section_vals_val_get(geo_section, "CONSTRAINT", &
                                c_vals=constraint_tokens, explicit=constraint_explicit)
      CALL parse_cell_constraint(constraint_tokens, cell, constraint_explicit, &
                                 cell_env%constraint_id)

      IF (cell_env%keep_symmetry .AND. cell_env%opt_cell%symmetry_id == cell_sym_none) THEN
         ! The crystal system is only ever set by the SYMMETRY keyword of
         ! &SUBSYS/&CELL (default NONE, no detection): without it there is
         ! nothing to keep and the run silently proceeds unconstrained.
         CALL cp_warn(__LOCATION__, &
                      "CELL_OPT: KEEP_SYMMETRY is on, but &SUBSYS/&CELL SYMMETRY is not "// &
                      "set (default NONE). No crystal system is known, so no symmetry "// &
                      "constraint is applied. Set SYMMETRY to the crystal system of the "// &
                      "cell to keep it during the optimization.")
      END IF

      ! Only for explicitly set values: an unset keyword must stay unparsed so
      ! the restart echo keeps omitting the line, as it always has.
      IF (constraint_explicit) THEN
         ALLOCATE (canonical_tokens(SIZE(constraint_tokens)))
         canonical_tokens = constraint_tokens
         CALL uppercase(canonical_tokens)
         CALL section_vals_val_set(geo_section, "CONSTRAINT", c_vals_ptr=canonical_tokens)
      END IF

      CALL cell_opt_setup_manifold(cell_env)

      ! Get the external pressure
      CALL read_external_press_tensor(geo_section, cell_env%opt_cell, cell_env%pres_ext, &
                                      cell_env%mtrx, cell_env%input_to_opt)

      ! Print cell optimisation setup
      NULLIFY (logger)
      logger => cp_get_default_logger()
      output_unit = cp_print_key_unit_nr(logger, geo_section, "PRINT%CELL", extension=".Log")
      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T2,A,T61,F20.1)") &
            "CELL_OPT| Pressure tolerance [bar]: ", cp_unit_from_cp2k(cell_env%pres_tol, "bar")
         IF (cell_env%keep_volume) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,T78,A3)") &
               "CELL_OPT| Keep volume of cell: ", "YES"
         ELSE
            WRITE (UNIT=output_unit, FMT="(T2,A,T78,A3)") &
               "CELL_OPT| Keep volume of cell: ", " NO"
         END IF
         IF (cell_env%keep_angles) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,T78,A3)") &
               "CELL_OPT| Keep angles between the cell vectors: ", "YES"
         ELSE
            WRITE (UNIT=output_unit, FMT="(T2,A,T78,A3)") &
               "CELL_OPT| Keep angles between the cell vectors: ", " NO"
         END IF
         IF (cell_env%keep_symmetry) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,T78,A3)") &
               "CELL_OPT| Keep cell symmetry: ", "YES"
         ELSE
            WRITE (UNIT=output_unit, FMT="(T2,A,T78,A3)") &
               "CELL_OPT| Keep cell symmetry: ", " NO"
         END IF
         label = cell_opt_constraint_label(cell_env%constraint_id)
         IF (COUNT(cell%perd /= 0) == 2 .AND. .NOT. constraint_explicit) THEN
            label = TRIM(label)//" (auto)"
         END IF
         WRITE (UNIT=output_unit, FMT="(T2,A,T77,A)") &
            "CELL_OPT| Constraint: ", TRIM(label)
      END IF
      CALL cp_print_key_finished_output(output_unit, logger, geo_section, "PRINT%CELL")

   END SUBROUTINE cell_opt_env_create

! **************************************************************************************************
!> \brief Initialize the cell constraint manifold from CONSTRAINT, KEEP_SYMMETRY,
!>        KEEP_ANGLES, and KEEP_VOLUME.
!> \param cell_env ...
! **************************************************************************************************
   SUBROUTINE cell_opt_setup_manifold(cell_env)
      TYPE(cell_opt_env_type), INTENT(INOUT)             :: cell_env

      INTEGER                                            :: angle_bits, b, info
      REAL(KIND=dp)                                      :: q0(n_slot)
      TYPE(cell_manifold_type)                           :: m_probe

      q0 = [cell_env%opt_cell%hmat(1, 1), cell_env%opt_cell%hmat(1, 2), &
            cell_env%opt_cell%hmat(2, 2), cell_env%opt_cell%hmat(1, 3), &
            cell_env%opt_cell%hmat(2, 3), cell_env%opt_cell%hmat(3, 3)]

      ! Validate explicit CONSTRAINT tokens first so redundant explicit tokens
      ! (such as CONSTRAINT A LENGTH_A) fail LICQ immediately.
      CALL cell_manifold_init(cell_env%manifold, q0, cell_env%constraint_id, info)
      IF (info /= info_ok) THEN
         CALL cp_abort(__LOCATION__, &
                       "CELL_OPT: degenerate constraint set (info="// &
                       cp_to_string(info)//"): the active constraints are linearly "// &
                       "dependent and do not define a manifold.")
      END IF

      IF (cell_env%keep_symmetry) THEN
         CALL cell_opt_fold_symmetry_rows(cell_env, cell_env%opt_cell%symmetry_id)
      END IF

      IF (cell_env%keep_angles) THEN
         ! On an orthorhombic cell, pinning BX = CX = CY = 0 is identical to
         ! pinning alpha = beta = gamma = 90 deg and keeps off-diagonal hmat
         ! entries identically 0.0_dp so cell%orthorhombic stays TRUE.
         IF (cell_env%opt_cell%orthorhombic) THEN
            angle_bits = comp_bx + comp_cx + comp_cy
         ELSE
            angle_bits = gfix_alpha + gfix_beta + gfix_gamma
         END IF
         DO b = 0, n_constraint_bits - 1
            IF (BTEST(angle_bits, b) .AND. .NOT. BTEST(cell_env%manifold%mask, b)) THEN
               CALL cell_manifold_init(m_probe, q0, IBSET(cell_env%manifold%mask, b), info, &
                                       lin_w=cell_env%manifold%lin_w, &
                                       n_lin_rows=cell_env%manifold%n_lin_rows, &
                                       eq_pair=cell_env%manifold%eq_pair, &
                                       n_eq_rows=cell_env%manifold%n_eq_rows)
               IF (info == info_ok) cell_env%manifold = m_probe
            END IF
         END DO
      END IF

      IF (cell_env%keep_volume) THEN
         IF (IAND(cell_env%manifold%mask, gfix_v) == 0) THEN
            CALL cell_manifold_init(m_probe, q0, IOR(cell_env%manifold%mask, gfix_v), info, &
                                    lin_w=cell_env%manifold%lin_w, &
                                    n_lin_rows=cell_env%manifold%n_lin_rows, &
                                    eq_pair=cell_env%manifold%eq_pair, &
                                    n_eq_rows=cell_env%manifold%n_eq_rows)
            IF (info == info_ok) cell_env%manifold = m_probe
         END IF
      END IF

   END SUBROUTINE cell_opt_setup_manifold

! **************************************************************************************************
!> \brief Greedily add the crystal-system rows of symmetry_id to the manifold.
!> \param cell_env ...
!> \param symmetry_id ...
!> \note Each row joins only when the extended set still passes the rank
!>       check, so rows folded earlier (KEEP_SYMMETRY with &CELL SYMMETRY)
!>       are skipped. KEEP_SPACE_GROUP calls this with the point group of
!>       its identified space group: the point group pins the cell metric
!>       exactly like the crystal system of that group.
! **************************************************************************************************
   SUBROUTINE cell_opt_fold_symmetry_rows(cell_env, symmetry_id)
      TYPE(cell_opt_env_type), INTENT(INOUT)             :: cell_env
      INTEGER, INTENT(IN)                                :: symmetry_id

      INTEGER                                            :: b, cand_eq_pair(n_geom, 2), cand_n_eq, &
                                                            cand_n_lin, i, info, sym_bits
      REAL(KIND=dp)                                      :: cand_lin_w(n_slot, n_slot)
      TYPE(cell_manifold_type)                           :: m_probe

      CALL symmetry_constraint_rows(symmetry_id, sym_bits, cand_lin_w, cand_n_lin, &
                                    cand_eq_pair, cand_n_eq)
      DO b = 0, n_slot - 1
         IF (BTEST(sym_bits, b) .AND. .NOT. BTEST(cell_env%manifold%mask, b)) THEN
            CALL cell_manifold_init(m_probe, cell_env%manifold%q0, &
                                    IBSET(cell_env%manifold%mask, b), info, &
                                    lin_w=cell_env%manifold%lin_w, &
                                    n_lin_rows=cell_env%manifold%n_lin_rows, &
                                    eq_pair=cell_env%manifold%eq_pair, &
                                    n_eq_rows=cell_env%manifold%n_eq_rows)
            IF (info == info_ok) cell_env%manifold = m_probe
         END IF
      END DO
      DO i = 1, cand_n_lin
         IF (cell_env%manifold%n_lin_rows >= n_slot) EXIT
         cell_env%manifold%lin_w(:, cell_env%manifold%n_lin_rows + 1) = cand_lin_w(:, i)
         CALL cell_manifold_init(m_probe, cell_env%manifold%q0, cell_env%manifold%mask, info, &
                                 lin_w=cell_env%manifold%lin_w, &
                                 n_lin_rows=cell_env%manifold%n_lin_rows + 1, &
                                 eq_pair=cell_env%manifold%eq_pair, &
                                 n_eq_rows=cell_env%manifold%n_eq_rows)
         IF (info == info_ok) THEN
            cell_env%manifold = m_probe
         ELSE
            cell_env%manifold%lin_w(:, cell_env%manifold%n_lin_rows + 1) = 0.0_dp
         END IF
      END DO
      DO i = 1, cand_n_eq
         IF (cell_env%manifold%n_eq_rows >= n_geom) EXIT
         cell_env%manifold%eq_pair(cell_env%manifold%n_eq_rows + 1, :) = cand_eq_pair(i, :)
         CALL cell_manifold_init(m_probe, cell_env%manifold%q0, cell_env%manifold%mask, info, &
                                 lin_w=cell_env%manifold%lin_w, &
                                 n_lin_rows=cell_env%manifold%n_lin_rows, &
                                 eq_pair=cell_env%manifold%eq_pair, &
                                 n_eq_rows=cell_env%manifold%n_eq_rows + 1)
         IF (info == info_ok) THEN
            cell_env%manifold = m_probe
         ELSE
            cell_env%manifold%eq_pair(cell_env%manifold%n_eq_rows + 1, :) = 0
         END IF
      END DO

   END SUBROUTINE cell_opt_fold_symmetry_rows

! **************************************************************************************************
!> \brief ...
!> \param cell_env ...
!> \par History
!>      none
!> \author Teodoro Laino - created [tlaino] - 03.2008 - Zurich University
! **************************************************************************************************
   SUBROUTINE cell_opt_env_release(cell_env)
      TYPE(cell_opt_env_type), INTENT(INOUT)             :: cell_env

      CALL cell_release(cell_env%opt_cell)

   END SUBROUTINE cell_opt_env_release

END MODULE cell_opt_types
