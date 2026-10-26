!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Unit test for cell_manifold_utils (CELL_OPT constraint manifold, #6031)
!> \par Testing principle
!>      The analytic Jacobian is verified against central finite differences
!>      on five cell shapes, plus golden values cross-validated to 2e-16
!>      against an independent Python implementation of the same math. The
!>      null-space basis is checked through its defining properties
!>      (A Z = 0, Z^T Z = I) rather than through the QR code path. The
!>      retraction is checked through its defining property (first-order
!>      consistency: the normal correction is O(t^2)) and through exact
!>      constraint restoration. The bit-identity of the component masks with
!>      cell_opt_utils is asserted over all 64 legacy masks, so the two
!>      modules cannot drift apart.
!> \par Coverage
!>      Geometry and Jacobian (FD + golden), null-space over representative
!>      constraint sets, LICQ rejection of redundant pairs and mid-run rank
!>      loss (both exact Gram-Cholesky failure and pivot-ratio failure),
!>      degenerate/coplanar cell rejection on equality rows, gradient
!>      projection idempotency, curved-manifold retraction, the flat
!>      all-angles theorem (the {alpha,beta,gamma} manifold is the per-vector
!>      stretch cone: zero Gauss-Newton iterations by construction), mixed
!>      component/geometric sets, volume retraction, and first-order
!>      retraction consistency.
!> \par Failure behavior
!>      The first mismatch prints the offending quantity and stops with exit
!>      code 1. do_regtest runs this program from the UNIT_TESTS registry and
!>      reads only the exit code.
! **************************************************************************************************
PROGRAM cell_manifold_utils_TEST
   USE cell_manifold_utils,             ONLY: &
        cell_manifold_active_jacobian, cell_manifold_expand_step, cell_manifold_geometry, &
        cell_manifold_init, cell_manifold_jacobian, cell_manifold_proj_pressure, &
        cell_manifold_project_grad, cell_manifold_residual, cell_manifold_retract, &
        cell_manifold_type, cell_manifold_update_basis, comp_ax_cm => comp_ax, &
        comp_bx_cm => comp_bx, comp_by_cm => comp_by, comp_cx_cm => comp_cx, &
        comp_cy_cm => comp_cy, comp_cz_cm => comp_cz, gfix_a, gfix_alpha, gfix_b, gfix_beta, &
        gfix_c, gfix_gamma, gfix_v, info_bad_input, info_degenerate, info_licq, info_no_converge, &
        info_ok, mask_a_cm => mask_a, mask_b_cm => mask_b, mask_c_cm => mask_c, &
        mask_x_cm => mask_x, mask_y_cm => mask_y, mask_z_cm => mask_z, n_geom, n_slot
   USE cell_opt_utils,                  ONLY: &
        comp_ax_cou => comp_ax, comp_bx_cou => comp_bx, comp_by_cou => comp_by, &
        comp_cx_cou => comp_cx, comp_cy_cou => comp_cy, comp_cz_cou => comp_cz, &
        mask_a_cou => mask_a, mask_b_cou => mask_b, mask_c_cou => mask_c, mask_x_cou => mask_x, &
        mask_y_cou => mask_y, mask_z_cou => mask_z
   USE ieee_arithmetic,                 ONLY: ieee_positive_inf,&
                                              ieee_value
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi

   IMPLICIT NONE

   CALL check_bit_compat_64()
   CALL check_geometry()
   CALL check_jacobian_fd()
   CALL check_golden_jacobian()
   CALL check_nullspace()
   CALL check_licq_redundant()
   CALL check_update_basis_rank()
   CALL check_projection()
   CALL check_proj_pressure()
   CALL check_retraction_curved()
   CALL check_flat_angles()
   CALL check_symmetry_rows()
   CALL check_kernel_robustness()
   CALL check_mixed_and_volume()
   CALL check_first_order()

   PRINT *, "cell_manifold_utils_unittest: all checks passed"

CONTAINS

! **************************************************************************************************
!> \brief The 64 legacy component masks of cell_manifold_utils must be
!>        bit-identical to cell_opt_utils (single source of semantics for
!>        parse_cell_constraint masks).
! **************************************************************************************************
   SUBROUTINE check_bit_compat_64()
      INTEGER                                            :: mask_value

      DO mask_value = 0, 63
         IF (IAND(cell_manifold_bit(mask_value), cell_opt_bit(mask_value)) /= &
             cell_opt_bit(mask_value) .OR. &
             cell_manifold_bit(mask_value) /= cell_opt_bit(mask_value)) THEN
            PRINT *, "check_bit_compat_64: mask mismatch at", mask_value
            ERROR STOP 1
         END IF
      END DO
      ! union shorthands must agree as well
      IF (mask_a_cm /= mask_a_cou .OR. mask_b_cm /= mask_b_cou .OR. &
          mask_c_cm /= mask_c_cou .OR. mask_x_cm /= mask_x_cou .OR. &
          mask_y_cm /= mask_y_cou .OR. mask_z_cm /= mask_z_cou) THEN
         PRINT *, "check_bit_compat_64: union mask mismatch"
         ERROR STOP 1
      END IF
      PRINT *, "check_bit_compat_64: 64 masks bit-identical to cell_opt_utils"
   END SUBROUTINE check_bit_compat_64

! **************************************************************************************************
!> \brief ...
!> \param m ...
!> \return ...
! **************************************************************************************************
   FUNCTION cell_manifold_bit(m) RESULT(b)
      INTEGER, INTENT(IN)                                :: m
      INTEGER                                            :: b

      INTEGER                                            :: cm(n_slot)

      cm = [comp_ax_cm, comp_bx_cm, comp_by_cm, comp_cx_cm, comp_cy_cm, comp_cz_cm]
      b = 0
      IF (BTEST(m, 0)) b = b + cm(1)
      IF (BTEST(m, 1)) b = b + cm(2)
      IF (BTEST(m, 2)) b = b + cm(3)
      IF (BTEST(m, 3)) b = b + cm(4)
      IF (BTEST(m, 4)) b = b + cm(5)
      IF (BTEST(m, 5)) b = b + cm(6)
   END FUNCTION cell_manifold_bit

! **************************************************************************************************
!> \brief ...
!> \param m ...
!> \return ...
! **************************************************************************************************
   FUNCTION cell_opt_bit(m) RESULT(b)
      INTEGER, INTENT(IN)                                :: m
      INTEGER                                            :: b

      INTEGER                                            :: co(n_slot)

      co = [comp_ax_cou, comp_bx_cou, comp_by_cou, comp_cx_cou, comp_cy_cou, comp_cz_cou]
      b = 0
      IF (BTEST(m, 0)) b = b + co(1)
      IF (BTEST(m, 1)) b = b + co(2)
      IF (BTEST(m, 2)) b = b + co(3)
      IF (BTEST(m, 3)) b = b + co(4)
      IF (BTEST(m, 4)) b = b + co(5)
      IF (BTEST(m, 5)) b = b + co(6)
   END FUNCTION cell_opt_bit

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE check_geometry()
      INTEGER                                            :: info
      REAL(KIND=dp)                                      :: q(n_slot), s(n_geom)

      q = monoclinic_q0()
      CALL cell_manifold_geometry(q, s, info)
      IF (info /= info_ok) ERROR STOP 1
      IF (ABS(s(4) - 0.5_dp*pi) > 1.0E-12_dp .OR. &
          ABS(s(5) - 105.0_dp*pi/180.0_dp) > 1.0E-12_dp .OR. &
          ABS(s(6) - 0.5_dp*pi) > 1.0E-12_dp .OR. &
          ABS(s(7) - q(1)*q(3)*q(6)) > 1.0E-12_dp) THEN
         PRINT *, "check_geometry: monoclinic angles/volume wrong"
         ERROR STOP 1
      END IF
      q = [0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
      CALL cell_manifold_geometry(q, s, info)
      IF (info /= info_degenerate) THEN
         PRINT *, "check_geometry: degenerate cell not detected"
         ERROR STOP 1
      END IF
      ! Coplanar cell (CZ = 0 with valid lengths and non-collinear angles) must also be rejected
      q = [1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp]
      CALL cell_manifold_geometry(q, s, info)
      IF (info /= info_degenerate) THEN
         PRINT *, "check_geometry: coplanar zero-volume cell not detected"
         ERROR STOP 1
      END IF
      PRINT *, "check_geometry: angles, volume, degeneracy OK"
   END SUBROUTINE check_geometry

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE check_jacobian_fd()
      INTEGER                                            :: c, i, info, jj
      REAL(KIND=dp)                                      :: J(n_geom, n_slot), Jfd(n_geom, n_slot), &
                                                            q(n_slot), rel
      REAL(KIND=dp), ALLOCATABLE                         :: cells(:, :)

      ALLOCATE (cells(n_slot, 5))
      cells(:, 1) = monoclinic_q0()
      cells(:, 2) = triclinic_q0()
      cells(:, 3) = [4.0_dp, -2.0_dp, 2.0_dp*SQRT(3.0_dp), 0.0_dp, 0.0_dp, 6.5_dp]
      cells(:, 4) = [3.5_dp, 0.0_dp, 3.5_dp, 0.0_dp, 0.0_dp, 3.5_dp]
      cells(:, 5) = [4.7_dp, -1.1_dp, 5.3_dp, 0.6_dp, 2.2_dp, 4.1_dp]
      rel = 0.0_dp
      DO c = 1, 5
         q = cells(:, c)
         CALL cell_manifold_jacobian(q, J, info)
         IF (info /= info_ok) ERROR STOP 1
         CALL fd_jacobian(q, 1.0E-6_dp, Jfd)
         DO i = 1, n_geom
            DO jj = 1, n_slot
               rel = MAX(rel, ABS(J(i, jj) - Jfd(i, jj))/MAX(1.0_dp, ABS(Jfd(i, jj))))
            END DO
         END DO
      END DO
      DEALLOCATE (cells)
      IF (rel > 1.0E-7_dp) THEN
         PRINT *, "check_jacobian_fd: rel err", rel
         ERROR STOP 1
      END IF
      PRINT *, "check_jacobian_fd: max rel err", rel
   END SUBROUTINE check_jacobian_fd

! **************************************************************************************************
!> \brief Golden values cross-validated against an independent Python
!>        implementation (agreement 2e-16 on this cell).
! **************************************************************************************************
   SUBROUTINE check_golden_jacobian()
      INTEGER                                            :: info
      REAL(KIND=dp)                                      :: J(n_geom, n_slot), q(n_slot)

      q = triclinic_q0()
      CALL cell_manifold_jacobian(q, J, info)
      IF (info /= info_ok) ERROR STOP 1
      IF (ABS(J(4, 2) - 3.9793083023372110E-02_dp) > 1.0E-12_dp .OR. &
          ABS(J(4, 5) - (-1.5623896481013405E-01_dp)) > 1.0E-12_dp .OR. &
          ABS(J(5, 4) - (-1.6940129870132348E-01_dp)) > 1.0E-12_dp .OR. &
          ABS(J(6, 2) - (-1.9607843137254902E-01_dp)) > 1.0E-12_dp .OR. &
          ABS(J(7, 1) - 2.64E+01_dp) > 1.0E-12_dp .OR. &
          ABS(J(7, 3) - 2.75E+01_dp) > 1.0E-12_dp) THEN
         PRINT *, "check_golden_jacobian: golden values drifted"
         ERROR STOP 1
      END IF
      PRINT *, "check_golden_jacobian: golden values OK"
   END SUBROUTINE check_golden_jacobian

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE check_nullspace()
      INTEGER                                            :: i, info, jj, masks(6), mi, r
      REAL(KIND=dp)                                      :: az, q(n_slot), ztz
      REAL(KIND=dp), ALLOCATABLE                         :: A(:, :)
      TYPE(cell_manifold_type)                           :: m

      masks(1) = gfix_gamma
      masks(2) = gfix_alpha + gfix_beta + gfix_gamma
      masks(3) = mask_x_cm
      masks(4) = mask_a_cm + gfix_gamma
      masks(5) = gfix_a + gfix_b + gfix_c + gfix_v
      masks(6) = 0
      q = triclinic_q0()
      DO mi = 1, 6
         CALL cell_manifold_init(m, q, masks(mi), info)
         IF (info /= info_ok) THEN
            PRINT *, "check_nullspace: init failed at spec", mi
            ERROR STOP 1
         END IF
         ALLOCATE (A(m%k, n_slot))
         CALL cell_manifold_active_jacobian(m, q, A, info)
         IF (info /= info_ok) ERROR STOP 1
         az = 0.0_dp
         DO r = 1, m%k
            DO jj = 1, m%ndof
               az = MAX(az, ABS(SUM(A(r, :)*m%Z(:, jj))))
            END DO
         END DO
         ztz = 0.0_dp
         DO i = 1, m%ndof
            DO jj = 1, m%ndof
               ztz = MAX(ztz, ABS(DOT_PRODUCT(m%Z(:, i), m%Z(:, jj)) - &
                                  MERGE(1.0_dp, 0.0_dp, i == jj)))
            END DO
         END DO
         DEALLOCATE (A)
         IF (m%ndof /= n_slot - m%k .OR. az > 1.0E-12_dp .OR. ztz > 1.0E-12_dp) THEN
            PRINT *, "check_nullspace: spec", mi, "az", az, "ztz", ztz
            ERROR STOP 1
         END IF
      END DO
      PRINT *, "check_nullspace: 6 specs AZ=0, ZtZ=I"
   END SUBROUTINE check_nullspace

! **************************************************************************************************
!> \brief In the canonical frame fixing length a equals freezing AX: the
!> pair is linearly dependent and LICQ must reject it.
! **************************************************************************************************
   SUBROUTINE check_licq_redundant()
      INTEGER                                            :: info
      TYPE(cell_manifold_type)                           :: m

      CALL cell_manifold_init(m, triclinic_q0(), comp_ax_cm + gfix_a, info)
      IF (info /= info_licq) THEN
         PRINT *, "check_licq_redundant: A + gfix(a) not rejected, info =", info
         ERROR STOP 1
      END IF
      PRINT *, "check_licq_redundant: redundant pair rejected"
   END SUBROUTINE check_licq_redundant

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE check_projection()
      INTEGER                                            :: i, info
      REAL(KIND=dp)                                      :: diff, g(n_slot), gred(4), gtan(n_slot), &
                                                            pp(n_slot, n_slot), ptan(n_slot), &
                                                            q(n_slot)
      TYPE(cell_manifold_type)                           :: m

      q = triclinic_q0()
      CALL cell_manifold_init(m, q, gfix_alpha + gfix_beta, info)
      g = [1.7_dp, -0.4_dp, 2.2_dp, 0.9_dp, -1.1_dp, 0.3_dp]
      CALL cell_manifold_project_grad(m, g, gred)
      CALL cell_manifold_expand_step(m, gred, gtan)
      CALL cell_manifold_project_grad(m, gtan, gred)
      CALL cell_manifold_expand_step(m, gred, ptan)
      diff = MAXVAL(ABS(ptan - gtan))
      pp = MATMUL(m%Z(:, 1:m%ndof), TRANSPOSE(m%Z(:, 1:m%ndof)))
      DO i = 1, n_slot
         diff = MAX(diff, MAXVAL(ABS(pp(:, i) - MATMUL(pp, pp(:, i)))))
      END DO
      IF (diff > 1.0E-13_dp) THEN
         PRINT *, "check_projection: idempotency broken", diff
         ERROR STOP 1
      END IF
      PRINT *, "check_projection: expand/project round trip, P^2=P"
   END SUBROUTINE check_projection

! **************************************************************************************************
!> \brief Projected-stress pressure residual: exact zero at constrained
!> stationary points (the KKT property the legacy masked-diagonal average
!> lacks), legacy free-slot equivalence on orthogonal component
!> constraints, and the two bypass cases.
! **************************************************************************************************
   SUBROUTINE check_proj_pressure()
      INTEGER                                            :: info
      LOGICAL                                            :: projected
      REAL(KIND=dp)                                      :: A(2, n_slot), g(n_slot), lambda(2), &
                                                            pdiff, q(n_slot), ref, uv(n_slot)
      TYPE(cell_manifold_type)                           :: m

      ! (1) KKT exactness: a pure constraint-reaction gradient (g = A^T
      ! lambda, the general shape of the enthalpy gradient at a constrained
      ! stationary point) must give exactly zero pressure residual - here
      ! for curved angle constraints, where the legacy criterion is blind
      q = triclinic_q0()
      CALL cell_manifold_init(m, q, gfix_beta + gfix_gamma, info)
      CALL cell_manifold_active_jacobian(m, q, A, info)
      IF (info /= info_ok) ERROR STOP 1
      lambda = [0.3_dp, -1.7_dp]
      g = MATMUL(TRANSPOSE(A), lambda)
      CALL cell_manifold_proj_pressure(m, q, g, pdiff, projected)
      IF (.NOT. projected .OR. ABS(pdiff) > 1.0E-11_dp) THEN
         PRINT *, "check_proj_pressure: KKT zero broken", projected, pdiff
         ERROR STOP 1
      END IF

      ! (2) orthogonal component constraints reduce to the legacy
      ! free-diagonal-slot average: BY and CZ pinned, slot AX free
      q = [5.0_dp, 0.0_dp, 6.0_dp, 0.0_dp, 0.0_dp, 7.0_dp]
      CALL cell_manifold_init(m, q, comp_by_cm + comp_cz_cm, info)
      g = 0.0_dp
      g(2) = 0.31_dp; g(4) = -0.22_dp; g(5) = 0.17_dp ! shear slots: vt has no support there
      g(1) = -q(3)*q(6)*2.5_dp ! = -V*(sigma_xx - p_ext)/AX with excess 2.5
      CALL cell_manifold_proj_pressure(m, q, g, pdiff, projected)
      IF (.NOT. projected .OR. ABS(pdiff - 2.5_dp) > 1.0E-12_dp) THEN
         PRINT *, "check_proj_pressure: legacy equivalence broken", projected, pdiff
         ERROR STOP 1
      END IF

      ! (3) pinned volume: bypass (the pressure direction is a constraint row)
      q = monoclinic_q0()
      CALL cell_manifold_init(m, q, gfix_v, info)
      CALL cell_manifold_proj_pressure(m, q, fixed_w(), pdiff, projected)
      IF (projected) THEN
         PRINT *, "check_proj_pressure: volume must bypass"
         ERROR STOP 1
      END IF

      ! (4) all diagonal slots pinned: bypass, a fully fixed cell has no
      ! pressure criterion
      CALL cell_manifold_init(m, q, comp_ax_cm + comp_bx_cm + comp_by_cm + &
                              comp_cx_cm + comp_cy_cm + comp_cz_cm, info)
      CALL cell_manifold_proj_pressure(m, q, fixed_w(), pdiff, projected)
      IF (projected) THEN
         PRINT *, "check_proj_pressure: all-diagonal must bypass"
         ERROR STOP 1
      END IF

      ! (5) unconstrained: plain least-squares pressure along dV/dq
      CALL cell_manifold_init(m, q, 0, info)
      g = fixed_w()
      uv = 0.0_dp
      uv(1) = q(3)*q(6)
      uv(3) = q(1)*q(6)
      uv(6) = q(1)*q(3)
      ref = -DOT_PRODUCT(g, uv)/DOT_PRODUCT(uv, uv)
      CALL cell_manifold_proj_pressure(m, q, g, pdiff, projected)
      IF (.NOT. projected .OR. ABS(pdiff - ref) > 1.0E-13_dp) THEN
         PRINT *, "check_proj_pressure: fix_none formula broken", pdiff, ref
         ERROR STOP 1
      END IF
      PRINT *, "check_proj_pressure: KKT zero, legacy equivalence, bypasses OK"
   END SUBROUTINE check_proj_pressure

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE check_retraction_curved()
      INTEGER                                            :: info
      REAL(KIND=dp)                                      :: dqtan(n_slot), du(4), q(n_slot), &
                                                            qout(n_slot), rel, rpre(2), w(n_slot)
      TYPE(cell_manifold_type)                           :: m

      q = triclinic_q0()
      CALL cell_manifold_init(m, q, gfix_beta + gfix_gamma, info)
      w = fixed_w()
      CALL cell_manifold_project_grad(m, w, du)
      CALL cell_manifold_expand_step(m, du, dqtan)
      CALL cell_manifold_residual(m, q + dqtan, rpre, info)
      IF (info /= info_ok) ERROR STOP 1
      CALL cell_manifold_retract(m, q, du, qout, info)
      rel = normal_rel(m, qout, q, dqtan)
      IF (NORM2(rpre) <= 1.0E-8_dp) THEN
         PRINT *, "check_retraction_curved: no visible curvature drift"
         ERROR STOP 1
      END IF
      IF (info /= info_ok .OR. m%last_resid > 1.0E-12_dp .OR. m%last_iter < 1) THEN
         PRINT *, "check_retraction_curved: not restored, info", info
         ERROR STOP 1
      END IF
      IF (rel > 0.10_dp) THEN
         PRINT *, "check_retraction_curved: drift", rel
         ERROR STOP 1
      END IF
      PRINT *, "check_retraction_curved: restored in", m%last_iter, "GN iters"
   END SUBROUTINE check_retraction_curved

! **************************************************************************************************
!> \brief Fixing all three angles is exactly integrable: the manifold is the
!> per-vector stretch cone, any tangent line stays on it, so the
!> Gauss-Newton correction must terminate in zero iterations.
! **************************************************************************************************
   SUBROUTINE check_flat_angles()
      INTEGER                                            :: c, info
      REAL(KIND=dp)                                      :: du(3), q(n_slot), qout(n_slot)
      TYPE(cell_manifold_type)                           :: m

      DO c = 1, 2
         IF (c == 1) q = triclinic_q0()
         IF (c == 2) q = monoclinic_q0()
         CALL cell_manifold_init(m, q, gfix_alpha + gfix_beta + gfix_gamma, info)
         CALL cell_manifold_project_grad(m, 2.0_dp*fixed_w(), du)
         CALL cell_manifold_retract(m, q, du, qout, info)
         IF (info /= info_ok .OR. m%last_iter /= 0 .OR. m%last_resid > 1.0E-12_dp) THEN
            PRINT *, "check_flat_angles: cell", c, "iters", m%last_iter
            ERROR STOP 1
         END IF
      END DO
      PRINT *, "check_flat_angles: stretch-cone flatness holds (0 GN iters)"
   END SUBROUTINE check_flat_angles

! **************************************************************************************************
!> \brief Crystal-symmetry rows (linear couplings and geometric equalities):
!> slot pins and linear coupling rows are exactly flat (zero Gauss-Newton
!> iterations by construction, like the stretch cone); geometric equality
!> rows are curved and must be restored by the retraction (rhombohedral
!> a=b=c, alpha=beta=gamma asserted here; the full crystal-system matrix
!> is asserted in cell_opt_utils_unittest). Redundant equality sets hit
!> LICQ, malformed rows are rejected.
! **************************************************************************************************
   SUBROUTINE check_symmetry_rows()
      INTEGER                                            :: eqp(n_geom, 2), info
      REAL(KIND=dp)                                      :: ct, cx, cy, cz, du(n_slot), &
                                                            lin(n_slot, n_slot), q(n_slot), &
                                                            qout(n_slot), rh(n_slot), s(n_geom), &
                                                            st, worst
      TYPE(cell_manifold_type)                           :: m

      ! only 1:ndof is written by project_grad; zero-fill so the 0.05*du
      ! argument expression never touches uninitialized elements (traps
      ! under -ffpe-trap=invalid)
      du = 0.0_dp

      ! --- flat systems: slots + linear coupling rows, 0 GN iterations ---
      ! (1) tetragonal_ab: BX,CX,CY = 0 and AX = BY
      q = [5.1_dp, 0.0_dp, 5.1_dp, 0.0_dp, 0.0_dp, 6.3_dp]
      lin = 0.0_dp
      lin(1, 1) = 1.0_dp
      lin(3, 1) = -1.0_dp
      CALL cell_manifold_init(m, q, comp_bx_cm + comp_cx_cm + comp_cy_cm, info, &
                              lin_w=lin, n_lin_rows=1)
      IF (info /= info_ok .OR. m%ndof /= 2) THEN
         PRINT *, "check_symmetry_rows: tetragonal init", info, m%ndof
         ERROR STOP 1
      END IF
      CALL cell_manifold_project_grad(m, fixed_w(), du)
      CALL cell_manifold_retract(m, q, 0.05_dp*du, qout, info)
      IF (info /= info_ok .OR. m%last_iter /= 0) THEN
         PRINT *, "check_symmetry_rows: tetragonal not flat", info, m%last_iter
         ERROR STOP 1
      END IF

      ! (2) cubic: slots + AX=BY and BY=CZ
      q = [5.1_dp, 0.0_dp, 5.1_dp, 0.0_dp, 0.0_dp, 5.1_dp]
      lin = 0.0_dp
      lin(1, 1) = 1.0_dp
      lin(3, 1) = -1.0_dp
      lin(3, 2) = 1.0_dp
      lin(6, 2) = -1.0_dp
      CALL cell_manifold_init(m, q, comp_bx_cm + comp_cx_cm + comp_cy_cm, info, &
                              lin_w=lin, n_lin_rows=2)
      IF (info /= info_ok .OR. m%ndof /= 1) THEN
         PRINT *, "check_symmetry_rows: cubic init", info, m%ndof
         ERROR STOP 1
      END IF
      CALL cell_manifold_project_grad(m, fixed_w(), du)
      CALL cell_manifold_retract(m, q, 0.05_dp*du, qout, info)
      IF (info /= info_ok .OR. m%last_iter /= 0 .OR. ABS(qout(1) - qout(3)) > 1.0E-12_dp) THEN
         PRINT *, "check_symmetry_rows: cubic not flat/equal", info, m%last_iter
         ERROR STOP 1
      END IF

      ! (3) hexagonal gamma=120: CX,CY = 0, BX + AX/2 = 0, BY - sqrt(3)/2 AX = 0
      q = [5.1_dp, -2.55_dp, 0.5_dp*SQRT(3.0_dp)*5.1_dp, 0.0_dp, 0.0_dp, 6.3_dp]
      lin = 0.0_dp
      lin(1, 1) = 0.5_dp
      lin(2, 1) = 1.0_dp
      lin(1, 2) = -0.5_dp*SQRT(3.0_dp)
      lin(3, 2) = 1.0_dp
      CALL cell_manifold_init(m, q, comp_cx_cm + comp_cy_cm, info, lin_w=lin, n_lin_rows=2)
      IF (info /= info_ok .OR. m%ndof /= 2) THEN
         PRINT *, "check_symmetry_rows: hexagonal init", info, m%ndof
         ERROR STOP 1
      END IF
      CALL cell_manifold_project_grad(m, fixed_w(), du)
      CALL cell_manifold_retract(m, q, 0.05_dp*du, qout, info)
      IF (info /= info_ok .OR. m%last_iter /= 0) THEN
         PRINT *, "check_symmetry_rows: hexagonal not flat", info, m%last_iter
         ERROR STOP 1
      END IF

      ! --- rhombohedral: 4 equality rows, curved, 2 dof (length AND angle) ---
      ct = COS(70.0_dp*pi/180.0_dp)
      st = SIN(70.0_dp*pi/180.0_dp)
      cx = ct
      cy = (ct - ct*ct)/st
      cz = SQRT(MAX(1.0_dp - cx*cx - cy*cy, 0.0_dp))
      rh = 5.1_dp*[1.0_dp, ct, st, cx, cy, cz]
      eqp = 0
      eqp(1, :) = [1, 2]
      eqp(2, :) = [2, 3]
      eqp(3, :) = [4, 5]
      eqp(4, :) = [5, 6]
      CALL cell_manifold_init(m, rh, 0, info, eq_pair=eqp, n_eq_rows=4)
      IF (info /= info_ok .OR. m%ndof /= 2) THEN
         PRINT *, "check_symmetry_rows: rhombohedral init", info, m%ndof
         ERROR STOP 1
      END IF
      CALL cell_manifold_project_grad(m, fixed_w(), du)
      CALL cell_manifold_retract(m, rh, 0.05_dp*du, qout, info)
      CALL cell_manifold_geometry(qout, s, info)
      worst = MAX(MAX(ABS(s(1) - s(2)), ABS(s(2) - s(3))), &
                  MAX(ABS(s(4) - s(5)), ABS(s(5) - s(6))))
      IF (info /= info_ok .OR. m%last_iter < 1 .OR. worst > 1.0E-11_dp*5.1_dp) THEN
         PRINT *, "check_symmetry_rows: rhombohedral failed", info, m%last_iter, worst
         ERROR STOP 1
      END IF

      ! --- LICQ: a fully redundant equality chain a=b, b=c, c=a ---
      eqp(5, :) = [3, 1]
      CALL cell_manifold_init(m, rh, 0, info, eq_pair=eqp, n_eq_rows=5)
      IF (info /= info_licq) THEN
         PRINT *, "check_symmetry_rows: redundant chain not rejected", info
         ERROR STOP 1
      END IF

      ! --- bad input: V (index 7) in an equality, zero linear weight ---
      eqp = 0
      eqp(1, :) = [1, 7]
      CALL cell_manifold_init(m, rh, 0, info, eq_pair=eqp, n_eq_rows=1)
      IF (info /= info_bad_input) THEN
         PRINT *, "check_symmetry_rows: V equality not rejected", info
         ERROR STOP 1
      END IF
      lin = 0.0_dp
      CALL cell_manifold_init(m, rh, 0, info, lin_w=lin, n_lin_rows=1)
      IF (info /= info_bad_input) THEN
         PRINT *, "check_symmetry_rows: zero weight not rejected", info
         ERROR STOP 1
      END IF
      PRINT *, "check_symmetry_rows: flat systems 0 iters, rhombohedral restored, LICQ OK"
   END SUBROUTINE check_symmetry_rows

! **************************************************************************************************
!> \brief Mid-run LICQ and degeneracy gates on update_basis: fixing all three
!>        side lengths plus volume passes on a skewed triclinic cell (ndof=2),
!>        collapses in rank on a healthy orthorhombic cell (grad V lies in the
!>        span of grad a, grad b, grad c -> exact Cholesky failure), trips the
!>        pivot-ratio gate near orthorhombicity, reports info_degenerate on a
!>        collapsed cell, and recovers on a valid skewed cell.
! **************************************************************************************************
   SUBROUTINE check_update_basis_rank()
      INTEGER                                            :: info
      REAL(KIND=dp)                                      :: q(n_slot)
      TYPE(cell_manifold_type)                           :: m

      q = triclinic_q0()
      CALL cell_manifold_init(m, q, gfix_a + gfix_b + gfix_c + gfix_v, info)
      IF (info /= info_ok .OR. m%ndof /= 2) THEN
         PRINT *, "check_update_basis_rank: init", info, m%ndof
         ERROR STOP 1
      END IF

      ! exact Jacobian rank collapse on a healthy orthorhombic cell:
      ! grad V = BY*CZ e_1 + AX*CZ e_3 + AX*BY e_6 lies in span(grad a, grad b, grad c)
      q = [5.0_dp, 0.0_dp, 5.0_dp, 0.0_dp, 0.0_dp, 5.0_dp]
      CALL cell_manifold_update_basis(m, q, info)
      IF (info /= info_licq) THEN
         PRINT *, "check_update_basis_rank: exact collapse not rejected", info
         ERROR STOP 1
      END IF

      ! near rank collapse on a healthy cell: Cholesky succeeds, pivot ratio < licq_pivot_tol
      q = [5.0_dp, 1.0E-6_dp, 5.0_dp, 0.0_dp, 0.0_dp, 5.0_dp]
      CALL cell_manifold_update_basis(m, q, info)
      IF (info /= info_licq) THEN
         PRINT *, "check_update_basis_rank: near collapse not rejected", info
         ERROR STOP 1
      END IF

      ! degenerate cell (BY = 0) must report info_degenerate rather than info_licq
      q = [5.0_dp, 3.7_dp, 0.0_dp, 0.0_dp, 0.0_dp, 5.0_dp]
      CALL cell_manifold_update_basis(m, q, info)
      IF (info /= info_degenerate) THEN
         PRINT *, "check_update_basis_rank: degenerate cell not reported", info
         ERROR STOP 1
      END IF

      ! recovery on a valid skewed cell
      q = triclinic_q0()
      CALL cell_manifold_update_basis(m, q, info)
      IF (info /= info_ok) THEN
         PRINT *, "check_update_basis_rank: recovery rejected", info
         ERROR STOP 1
      END IF
      PRINT *, "check_update_basis_rank: exact/near LICQ, degeneracy, recovery OK"
   END SUBROUTINE check_update_basis_rank

! **************************************************************************************************
!> \brief Robustness checks: scale invariance of linear coupling rows,
!>      non-finite steps rejected, degenerate probe points rejected even on
!>      pure equality constraints (rhombohedral), and weighted-mean projected
!>      pressure across multiple free diagonal slots.
! **************************************************************************************************
   SUBROUTINE check_kernel_robustness()
      INTEGER                                            :: eqp(n_geom, 2), info
      LOGICAL                                            :: projected
      REAL(KIND=dp)                                      :: ct, cx, cy, cz, du(2), g(n_slot), &
                                                            lin(n_slot, n_slot), pdiff, q(n_slot), &
                                                            qout(n_slot), ref, rh(n_slot), st
      TYPE(cell_manifold_type)                           :: m

      ! (1) scale invariance of linear coupling rows
      q = [5.1_dp, 0.0_dp, 5.1_dp, 0.0_dp, 0.0_dp, 6.3_dp]
      lin = 0.0_dp
      lin(1, 1) = 1.0E9_dp
      lin(3, 1) = -1.0E9_dp
      CALL cell_manifold_init(m, q, comp_bx_cm + comp_cx_cm + comp_cy_cm, info, &
                              lin_w=lin, n_lin_rows=1)
      IF (info /= info_ok) THEN
         PRINT *, "check_kernel_robustness: scaled lin row rejected", info
         ERROR STOP 1
      END IF
      CALL cell_manifold_project_grad(m, fixed_w(), du)
      CALL cell_manifold_retract(m, q, 27.0_dp*du, qout, info)
      IF (info /= info_ok .OR. m%last_resid > 1.0E-12_dp) THEN
         PRINT *, "check_kernel_robustness: scaled lin retraction failed", info, m%last_resid
         ERROR STOP 1
      END IF

      ! (2) non-finite step: reported, never validated as converged
      CALL cell_manifold_init(m, q, gfix_beta, info)
      du = 0.0_dp
      du(1) = ieee_value(1.0_dp, ieee_positive_inf)
      CALL cell_manifold_retract(m, q, du, qout, info)
      IF (info /= info_no_converge) THEN
         PRINT *, "check_kernel_robustness: non-finite step not reported", info
         ERROR STOP 1
      END IF

      ! (3) degenerate probe point on pure equality constraints (rhombohedral)
      ! must return info_degenerate, never cancel to zero residual
      ct = COS(70.0_dp*pi/180.0_dp)
      st = SIN(70.0_dp*pi/180.0_dp)
      cx = ct
      cy = (ct - ct*ct)/st
      cz = SQRT(MAX(1.0_dp - cx*cx - cy*cy, 0.0_dp))
      rh = 5.1_dp*[1.0_dp, ct, st, cx, cy, cz]
      eqp = 0
      eqp(1, :) = [1, 2]
      eqp(2, :) = [2, 3]
      eqp(3, :) = [4, 5]
      eqp(4, :) = [5, 6]
      CALL cell_manifold_init(m, rh, 0, info, eq_pair=eqp, n_eq_rows=4)
      IF (info /= info_ok) ERROR STOP 1
      du = 0.0_dp
      CALL cell_manifold_retract(m, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], du, qout, info)
      IF (info /= info_degenerate) THEN
         PRINT *, "check_kernel_robustness: degenerate rhombohedral point accepted, info =", info
         ERROR STOP 1
      END IF

      ! (4) several free diagonal slots: (dV/dq)-weighted mean
      q = [5.0_dp, 0.0_dp, 6.0_dp, 0.0_dp, 0.0_dp, 7.0_dp]
      CALL cell_manifold_init(m, q, 0, info)
      g = 0.0_dp
      g(1) = -2.0_dp*q(3)*q(6)
      g(6) = -5.0_dp*q(1)*q(3)
      ref = (2.0_dp*(q(3)*q(6))**2 + 5.0_dp*(q(1)*q(3))**2)/ &
            ((q(3)*q(6))**2 + (q(1)*q(6))**2 + (q(1)*q(3))**2)
      CALL cell_manifold_proj_pressure(m, q, g, pdiff, projected)
      IF (.NOT. projected .OR. ABS(pdiff - ref) > 1.0E-12_dp*ref) THEN
         PRINT *, "check_kernel_robustness: weighted mean broken", pdiff, ref
         ERROR STOP 1
      END IF
      PRINT *, "check_kernel_robustness: scale invariance, NaN/degeneracy guards, weighted mean OK"
   END SUBROUTINE check_kernel_robustness

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE check_mixed_and_volume()
      INTEGER                                            :: info
      REAL(KIND=dp)                                      :: du(5), q(n_slot), qout(n_slot)
      TYPE(cell_manifold_type)                           :: m

      q = triclinic_q0()
      CALL cell_manifold_init(m, q, comp_ax_cm + gfix_b, info)
      CALL cell_manifold_project_grad(m, fixed_w(), du)
      CALL cell_manifold_retract(m, q, du, qout, info)
      IF (info /= info_ok .OR. ABS(qout(1) - q(1)) > 1.0E-12_dp .OR. &
          m%last_resid > 1.0E-12_dp) THEN
         PRINT *, "check_mixed_and_volume: mixed A + b failed"
         ERROR STOP 1
      END IF
      q = monoclinic_q0()
      CALL cell_manifold_init(m, q, gfix_v, info)
      CALL cell_manifold_project_grad(m, fixed_w(), du)
      CALL cell_manifold_retract(m, q, du, qout, info)
      ! volume targets are large: the residual converges to the relative
      ! tolerance used by the retraction, not to an absolute 1e-12
      IF (info /= info_ok .OR. m%last_resid > 1.0E-10_dp*ABS(q(1)*q(3)*q(6)) .OR. &
          m%last_iter < 1) THEN
         PRINT *, "check_mixed_and_volume: volume failed"
         ERROR STOP 1
      END IF
      PRINT *, "check_mixed_and_volume: mixed set and volume restored"
   END SUBROUTINE check_mixed_and_volume

! **************************************************************************************************
!> \brief Retraction defining property: |R(q, t*du) - (q + t*P du)| = O(t^2),
!> i.e. the error shrinks ~4x per halving of t.
! **************************************************************************************************
   SUBROUTINE check_first_order()
      INTEGER                                            :: info, it
      REAL(KIND=dp)                                      :: d(n_slot), du(4), e(4), p(n_slot), &
                                                            q(n_slot), qout(n_slot), t(4), &
                                                            w(n_slot)
      TYPE(cell_manifold_type)                           :: m

      q = triclinic_q0()
      CALL cell_manifold_init(m, q, gfix_beta + gfix_gamma, info)
      w = fixed_w()
      CALL cell_manifold_project_grad(m, w, du)
      CALL cell_manifold_expand_step(m, du, p)
      t = [0.2_dp, 0.1_dp, 0.05_dp, 0.025_dp]
      DO it = 1, 4
         CALL cell_manifold_retract(m, q, t(it)*du, qout, info)
         d = qout - q - t(it)*p
         e(it) = NORM2(d)
      END DO
      IF (e(3)/e(2) < 0.2_dp .OR. e(3)/e(2) > 0.3_dp .OR. &
          e(4)/e(3) < 0.2_dp .OR. e(4)/e(3) > 0.3_dp) THEN
         PRINT *, "check_first_order: ratios", e(3)/e(2), e(4)/e(3)
         ERROR STOP 1
      END IF
      PRINT *, "check_first_order: quadratic shrink ratios OK"
   END SUBROUTINE check_first_order

! **************************************************************************************************
!  helpers
! **************************************************************************************************
!> \brief Central-difference Jacobian for verifying cell_manifold_jacobian.
!> \param q ...
!> \param rel_step ...
!> \param J ...
! **************************************************************************************************
   SUBROUTINE fd_jacobian(q, rel_step, J)
      REAL(KIND=dp), INTENT(IN)                          :: q(n_slot), rel_step
      REAL(KIND=dp), INTENT(OUT)                         :: J(n_geom, n_slot)

      INTEGER                                            :: i, info_chk
      REAL(KIND=dp)                                      :: h, qm(n_slot), qp(n_slot), sm(n_geom), &
                                                            sp(n_geom)

      DO i = 1, n_slot
         h = rel_step*MAX(1.0_dp, ABS(q(i)))
         qp = q
         qm = q
         qp(i) = qp(i) + h
         qm(i) = qm(i) - h
         CALL cell_manifold_geometry(qp, sp, info_chk)
         IF (info_chk /= info_ok) ERROR STOP 1
         CALL cell_manifold_geometry(qm, sm, info_chk)
         IF (info_chk /= info_ok) ERROR STOP 1
         J(:, i) = (sp - sm)/(2.0_dp*h)
      END DO
   END SUBROUTINE fd_jacobian

! **************************************************************************************************
!> \brief ...
!> \return ...
! **************************************************************************************************
   FUNCTION monoclinic_q0() RESULT(q)
      REAL(KIND=dp)                                      :: q(n_slot)

      q = [5.1_dp, 0.0_dp, 4.9_dp, 6.2_dp*COS(105.0_dp*pi/180.0_dp), &
           0.0_dp, 6.2_dp*SIN(105.0_dp*pi/180.0_dp)]
   END FUNCTION monoclinic_q0

! **************************************************************************************************
!> \brief ...
!> \return ...
! **************************************************************************************************
   FUNCTION triclinic_q0() RESULT(q)
      REAL(KIND=dp)                                      :: q(n_slot)

      q = [5.0_dp, 1.2_dp, 4.8_dp, -0.7_dp, 1.9_dp, 5.5_dp]
   END FUNCTION triclinic_q0

! **************************************************************************************************
!> \brief ...
!> \return ...
! **************************************************************************************************
   FUNCTION fixed_w() RESULT(w)
      REAL(KIND=dp)                                      :: w(n_slot)

      w = [0.13_dp, -0.27_dp, 0.11_dp, 0.31_dp, -0.07_dp, 0.19_dp]
   END FUNCTION fixed_w

! **************************************************************************************************
!> \brief ...
!> \param m ...
!> \param qout ...
!> \param q ...
!> \param dqtan ...
!> \return ...
! **************************************************************************************************
   FUNCTION normal_rel(m, qout, q, dqtan) RESULT(rel)
      TYPE(cell_manifold_type), INTENT(IN)               :: m
      REAL(KIND=dp), INTENT(IN)                          :: qout(n_slot), q(n_slot), dqtan(n_slot)
      REAL(KIND=dp)                                      :: rel

      REAL(KIND=dp)                                      :: d(n_slot), ptan(n_slot), du(m%ndof)

      d = qout - q
      CALL cell_manifold_project_grad(m, d, du)
      CALL cell_manifold_expand_step(m, du, ptan)
      rel = SQRT(MAX(0.0_dp, DOT_PRODUCT(d, d) - DOT_PRODUCT(ptan, ptan)))/NORM2(dqtan)
   END FUNCTION normal_rel

END PROGRAM cell_manifold_utils_TEST
