!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Cell constraint-manifold kernel for CELL_OPT (#6031).
!>
!> Pure-math layer: depends only on kinds, ieee_arithmetic, and LAPACK
!> (DGEQRF/DORGQR, DPOTRF/DPOTRS) - no force_env, cell or input dependencies,
!> so it stays unit-testable as a standalone program.
!>
!> Representation (identical to the CELL_OPT optimization vector):
!>   canonical cell (CELL_OPT_CANONICAL): a along x, b in the xy plane,
!>   hmat(:,1)=a=(AX,0,0), hmat(:,2)=b=(BX,BY,0), hmat(:,3)=c=(CX,CY,CZ);
!>   q = (AX,BX,BY,CX,CY,CZ) is the last-six cell block of the gopt vector
!>   (gopt_f_create_x0 / apply_cell_change) and of the get_dg_dh gradient.
!>   Angles: alpha = angle(b,c), beta = angle(a,c), gamma = angle(a,b),
!>   in radians; V = det(h) = AX*BY*CZ.
!>
!> Constraint model (Boumal 2023, sec. 7.7; issue #6031): the active set
!> defines M = {q : h(q) = 0}; component constraints contribute identity
!> rows, geometric constraints Jacobian rows; crystal symmetry adds
!> constant linear coupling rows (w . q = w . q0, exactly flat) and
!> geometric equality rows (s_i(q) = s_j(q), e.g. a = b of rhombohedral
!> symmetry). All rows live in one k x 6 matrix A(q); the tangent space is
!> ker A(q) with orthonormal basis Z (complete QR); the retraction is
!> Gauss-Newton on the residual along the constraint row space (Absil &
!> Malick 2012, projection-like retraction, second order).
!>
!> Gradient sign convention: get_dg_dh returns -dE/dq (force-like); every
!> routine here is linear in the gradient and sign-agnostic.
! **************************************************************************************************
MODULE cell_manifold_utils

   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp

   IMPLICIT NONE
   PRIVATE

   ! Component bits: bit-identical to comp_ax..comp_cz of cell_opt_utils
   ! (kept duplicated to preserve the zero-dependency property; the unit
   ! test asserts bit-by-bit equality over all 64 masks of cell_opt_utils).
   INTEGER, PARAMETER, PUBLIC :: comp_ax = 2**0, & ! a(x)
                                 comp_bx = 2**1, & ! b(x)
                                 comp_by = 2**2, & ! b(y)
                                 comp_cx = 2**3, & ! c(x)
                                 comp_cy = 2**4, & ! c(y)
                                 comp_cz = 2**5  ! c(z)
   INTEGER, PARAMETER, PUBLIC :: mask_a = comp_ax, &
                                 mask_b = comp_bx + comp_by, &
                                 mask_c = comp_cx + comp_cy + comp_cz, &
                                 mask_x = comp_ax + comp_bx + comp_cx, &
                                 mask_y = comp_bx + comp_by + comp_cy, &
                                 mask_z = mask_c, &
                                 mask_all = comp_ax + comp_bx + comp_by + comp_cx + comp_cy + comp_cz

   ! Geometric constraint bits (fix a length, an angle, or the volume to the
   ! starting value). parse_token of cell_opt_utils maps the CONSTRAINT
   ! tokens LENGTH_A..LENGTH_C, ALPHA, BETA, GAMMA and VOLUME onto them.
   ! Note: in the canonical frame fixing length a is the same constraint as
   ! freezing AX (a = |AX|); the LICQ check rejects the redundant pair.
   INTEGER, PARAMETER, PUBLIC :: gfix_a = 2**6, &
                                 gfix_b = 2**7, &
                                 gfix_c = 2**8, &
                                 gfix_alpha = 2**9, &
                                 gfix_beta = 2**10, &
                                 gfix_gamma = 2**11, &
                                 gfix_v = 2**12
   INTEGER, PARAMETER, PUBLIC :: gfix_all = gfix_a + gfix_b + gfix_c + gfix_alpha + &
                                 gfix_beta + gfix_gamma + gfix_v

   INTEGER, PARAMETER, PUBLIC :: n_slot = 6, n_geom = 7
   INTEGER, PARAMETER, PUBLIC :: n_constraint_bits = n_slot + n_geom

   ! Numerical thresholds for degeneracy, rank detection, and retraction
   REAL(KIND=dp), PARAMETER, PUBLIC  :: min_cell_length = 1.0E-8_dp
   REAL(KIND=dp), PARAMETER, PRIVATE :: collinear_cos_margin = 1.0E-8_dp
   REAL(KIND=dp), PARAMETER, PRIVATE :: licq_pivot_tol = 1.0E-6_dp
   REAL(KIND=dp), PARAMETER, PRIVATE :: default_retract_tol = 1.0E-12_dp
   INTEGER, PARAMETER, PRIVATE       :: default_retract_max_iter = 100
   REAL(KIND=dp), PARAMETER, PRIVATE :: proj_pressure_rel_tol = 1.0E-12_dp

   ! info codes (negative = error)
   INTEGER, PARAMETER, PUBLIC :: info_ok = 0
   INTEGER, PARAMETER, PUBLIC :: info_degenerate = -1 ! collapsed cell or 0/180 deg angle
   INTEGER, PARAMETER, PUBLIC :: info_licq = -2        ! A(q) row-rank deficient
   INTEGER, PARAMETER, PUBLIC :: info_no_converge = -3 ! retraction hit max_iter or non-finite step
   INTEGER, PARAMETER, PUBLIC :: info_bad_input = -4   ! invalid mask bit or row shape

   !> Active constraint set on the canonical cell q = (AX,BX,BY,CX,CY,CZ).
   TYPE, PUBLIC :: cell_manifold_type
      INTEGER              :: mask = 0        ! union of comp_*/gfix_* bits
      REAL(KIND=dp)        :: q0(n_slot) = 0.0_dp
      REAL(KIND=dp)        :: s0(n_geom) = 0.0_dp  ! targets: s(q0)
      INTEGER              :: k = 0           ! number of active constraints
      INTEGER              :: ndof = n_slot   ! 6 - k
      INTEGER              :: slot_rows(n_slot) = 0  ! frozen slot indices (1-based)
      INTEGER              :: n_slot_rows = 0
      INTEGER              :: geom_rows(n_geom) = 0  ! s-indices (1=a..7=V)
      INTEGER              :: n_geom_rows = 0
      ! Linear coupling rows w . q = w . q0 (crystal-system couplings such as
      ! AX=BY or the hexagonal gamma pattern; constant Jacobian, so their
      ! portion of the manifold is exactly flat). Column i holds the unit
      ! weights of row i.
      REAL(KIND=dp)        :: lin_w(n_slot, n_slot) = 0.0_dp
      INTEGER              :: n_lin_rows = 0
      ! Geometric equality rows s_i(q) = s_j(q) with i,j in 1..6 (a=b,
      ! alpha=beta; the natural description of rhombohedral and monoclinic
      ! a=b symmetry; no target needed - the residual is the difference).
      INTEGER              :: eq_pair(n_geom, 2) = 0
      INTEGER              :: n_eq_rows = 0
      REAL(KIND=dp)        :: Z(n_slot, n_slot) = 0.0_dp ! orthonormal, cols 1..ndof
      REAL(KIND=dp)        :: tol = default_retract_tol
      INTEGER              :: max_iter = default_retract_max_iter
      ! diagnostics of the last retract call
      INTEGER              :: last_iter = 0
      REAL(KIND=dp)        :: last_resid = 0.0_dp
   END TYPE cell_manifold_type

   PUBLIC :: cell_manifold_init, cell_manifold_geometry, cell_manifold_jacobian, &
             cell_manifold_residual, cell_manifold_active_jacobian, &
             cell_manifold_update_basis, cell_manifold_retract, &
             cell_manifold_project_grad, cell_manifold_expand_step, &
             cell_manifold_proj_pressure

   EXTERNAL :: dgeqrf, dorgqr, dpotrf, dpotrs

CONTAINS

! **************************************************************************************************
!> \brief Snapshot targets at q0, decompose the mask, build the tangent basis.
!> \param m manifold state
!> \param q0 canonical start cell, source of all constraint targets
!> \param mask union of comp_* (frozen slots) and gfix_* (geometric) bits
!> \param info info_ok, or info_degenerate / info_licq / info_bad_input
!> \param lin_w optional linear coupling rows, weights in column i of row i
!> \param n_lin_rows number of linear coupling rows (SIZE(lin_w, 2) used if absent)
!> \param eq_pair optional geometric equality rows, (s_i, s_j) in row i,
!>        indices 1..6 (a,b,c,alpha,beta,gamma; a V equality is meaningless
!>        and rejected)
!> \param n_eq_rows number of equality rows (SIZE(eq_pair, 1) used if absent)
! **************************************************************************************************
   SUBROUTINE cell_manifold_init(m, q0, mask, info, lin_w, n_lin_rows, eq_pair, n_eq_rows)
      TYPE(cell_manifold_type), INTENT(OUT)              :: m
      REAL(KIND=dp), INTENT(IN)                          :: q0(n_slot)
      INTEGER, INTENT(IN)                                :: mask
      INTEGER, INTENT(OUT)                               :: info
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: lin_w(:, :)
      INTEGER, INTENT(IN), OPTIONAL                      :: n_lin_rows, eq_pair(:, :), n_eq_rows

      INTEGER                                            :: i, ne, nl
      REAL(KIND=dp)                                      :: wnorm

      m%mask = mask
      m%q0 = q0

      CALL cell_manifold_geometry(q0, m%s0, info)
      IF (info /= info_ok) RETURN

      IF (IAND(mask, NOT(mask_all + gfix_all)) /= 0) THEN
         info = info_bad_input
         RETURN
      END IF

      m%n_slot_rows = 0
      DO i = 1, n_slot
         IF (BTEST(mask, i - 1)) THEN
            m%n_slot_rows = m%n_slot_rows + 1
            m%slot_rows(m%n_slot_rows) = i
         END IF
      END DO
      m%n_geom_rows = 0
      DO i = 1, n_geom
         IF (BTEST(mask, i - 1 + n_slot)) THEN
            m%n_geom_rows = m%n_geom_rows + 1
            m%geom_rows(m%n_geom_rows) = i
         END IF
      END DO

      m%n_lin_rows = 0
      IF (PRESENT(lin_w)) THEN
         IF (SIZE(lin_w, 1) /= n_slot) THEN
            info = info_bad_input
            RETURN
         END IF
         nl = SIZE(lin_w, 2)
         IF (PRESENT(n_lin_rows)) nl = n_lin_rows
         IF (nl < 0 .OR. nl > n_slot .OR. nl > SIZE(lin_w, 2)) THEN
            info = info_bad_input
            RETURN
         END IF
         DO i = 1, nl
            ! The constraint w . q = w . q0 is invariant under the scale of
            ! w, while the residual roundoff floor grows with |w|. Unit-norm
            ! rows keep the retraction tolerance scale-free.
            wnorm = NORM2(lin_w(:, i))
            IF (wnorm <= TINY(1.0_dp)) THEN
               info = info_bad_input
               RETURN
            END IF
            m%lin_w(:, i) = lin_w(:, i)/wnorm
         END DO
         m%n_lin_rows = nl
      ELSE IF (PRESENT(n_lin_rows)) THEN
         IF (n_lin_rows /= 0) THEN
            info = info_bad_input
            RETURN
         END IF
      END IF

      m%n_eq_rows = 0
      IF (PRESENT(eq_pair)) THEN
         IF (SIZE(eq_pair, 2) /= 2) THEN
            info = info_bad_input
            RETURN
         END IF
         ne = SIZE(eq_pair, 1)
         IF (PRESENT(n_eq_rows)) ne = n_eq_rows
         IF (ne < 0 .OR. ne > n_geom - 1 .OR. ne > SIZE(eq_pair, 1)) THEN
            info = info_bad_input
            RETURN
         END IF
         DO i = 1, ne
            IF (eq_pair(i, 1) < 1 .OR. eq_pair(i, 1) > n_geom - 1 .OR. &
                eq_pair(i, 2) < 1 .OR. eq_pair(i, 2) > n_geom - 1 .OR. &
                eq_pair(i, 1) == eq_pair(i, 2)) THEN
               info = info_bad_input
               RETURN
            END IF
            m%eq_pair(i, :) = eq_pair(i, :)
         END DO
         m%n_eq_rows = ne
      ELSE IF (PRESENT(n_eq_rows)) THEN
         IF (n_eq_rows /= 0) THEN
            info = info_bad_input
            RETURN
         END IF
      END IF

      m%k = m%n_slot_rows + m%n_geom_rows + m%n_lin_rows + m%n_eq_rows
      IF (m%k > n_slot) THEN
         info = info_licq
         RETURN
      END IF
      m%ndof = n_slot - m%k

      CALL cell_manifold_update_basis(m, q0, info)

   END SUBROUTINE cell_manifold_init

! **************************************************************************************************
!> \brief s(q) = (a, b, c, alpha, beta, gamma, V), angles in radians.
!> \param q ...
!> \param s ...
!> \param info ...
! **************************************************************************************************
   SUBROUTINE cell_manifold_geometry(q, s, info)
      REAL(KIND=dp), INTENT(IN)                          :: q(n_slot)
      REAL(KIND=dp), INTENT(OUT)                         :: s(n_geom)
      INTEGER, INTENT(OUT)                               :: info

      REAL(KIND=dp)                                      :: cosa(3), v1(3), v2(3), v3(3)

      s = 0.0_dp
      IF (.NOT. ALL(ieee_is_finite(q))) THEN
         info = info_degenerate
         RETURN
      END IF
      ! Positive upper-triangular diagonal ensures positive orientation and
      ! non-coplanar canonical vectors (V = AX*BY*CZ > 0).
      IF (MIN(q(1), q(3), q(6)) < min_cell_length) THEN
         info = info_degenerate
         RETURN
      END IF

      v1 = [q(1), 0.0_dp, 0.0_dp]
      v2 = [q(2), q(3), 0.0_dp]
      v3 = [q(4), q(5), q(6)]
      s(1) = NORM2(v1)
      s(2) = NORM2(v2)
      s(3) = NORM2(v3)
      cosa(1) = DOT_PRODUCT(v2, v3)/(s(2)*s(3))
      cosa(2) = DOT_PRODUCT(v1, v3)/(s(1)*s(3))
      cosa(3) = DOT_PRODUCT(v1, v2)/(s(1)*s(2))
      IF (ANY(ABS(cosa) > 1.0_dp - collinear_cos_margin)) THEN
         info = info_degenerate
         s = 0.0_dp
         RETURN
      END IF
      s(4:6) = ACOS(MIN(1.0_dp, MAX(-1.0_dp, cosa)))
      s(7) = q(1)*q(3)*q(6)
      info = info_ok

   END SUBROUTINE cell_manifold_geometry

! **************************************************************************************************
!> \brief Analytic Jacobian ds/dq (n_geom x n_slot); the rows are the
!>        constraint gradients (grad h_i in the manifold model).
!> \param q ...
!> \param J ...
!> \param info ...
! **************************************************************************************************
   SUBROUTINE cell_manifold_jacobian(q, J, info)
      REAL(KIND=dp), INTENT(IN)                          :: q(n_slot)
      REAL(KIND=dp), INTENT(OUT)                         :: J(n_geom, n_slot)
      INTEGER, INTENT(OUT)                               :: info

      REAL(KIND=dp)                                      :: cosang(3), db(n_slot), dc(n_slot), &
                                                            du(n_slot), s(n_geom), sina(3)

      J = 0.0_dp
      CALL cell_manifold_geometry(q, s, info)
      IF (info /= info_ok) RETURN

      ! d a / dq
      J(1, 1) = SIGN(1.0_dp, q(1))
      ! d b / dq
      db = 0.0_dp
      db(2) = q(2)/s(2)
      db(3) = q(3)/s(2)
      J(2, :) = db
      ! d c / dq
      dc = 0.0_dp
      dc(4) = q(4)/s(3)
      dc(5) = q(5)/s(3)
      dc(6) = q(6)/s(3)
      J(3, :) = dc

      cosang = [DOT_PRODUCT([q(2), q(3), 0.0_dp], [q(4), q(5), q(6)])/(s(2)*s(3)), &
                q(1)*q(4)/(s(1)*s(3)), &
                q(1)*q(2)/(s(1)*s(2))]
      sina = SIN(s(4:6))

      ! alpha = angle(b,c): d(b.c)/dq = (0, CX, CY, BX, BY, 0)
      du = [0.0_dp, q(4), q(5), q(2), q(3), 0.0_dp]
      J(4, :) = dangle(cosang(1), sina(1), du, s(2), s(3), db, dc)
      ! beta = angle(a,c): d(a.c)/dq = (CX, 0, 0, AX, 0, 0)
      du = [q(4), 0.0_dp, 0.0_dp, q(1), 0.0_dp, 0.0_dp]
      J(5, :) = dangle(cosang(2), sina(2), du, s(1), s(3), J(1, :), dc)
      ! gamma = angle(a,b): d(a.b)/dq = (BX, AX, 0, 0, 0, 0)
      du = [q(2), q(1), 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      J(6, :) = dangle(cosang(3), sina(3), du, s(1), s(2), J(1, :), db)

      ! d V / dq for V = AX*BY*CZ
      J(7, :) = [q(3)*q(6), 0.0_dp, q(1)*q(6), 0.0_dp, 0.0_dp, q(1)*q(3)]

   END SUBROUTINE cell_manifold_jacobian

! **************************************************************************************************
!> \brief r(q) = [q_i - q0_i ; s_lab(q) - s_lab(q0) ; s_a - s_b ; w . (q - q0)] - zero on M.
!> \param m ...
!> \param q ...
!> \param r ...
!> \param info ...
! **************************************************************************************************
   SUBROUTINE cell_manifold_residual(m, q, r, info)
      TYPE(cell_manifold_type), INTENT(IN)               :: m
      REAL(KIND=dp), INTENT(IN)                          :: q(n_slot)
      REAL(KIND=dp), INTENT(OUT)                         :: r(m%k)
      INTEGER, INTENT(OUT)                               :: info

      INTEGER                                            :: i, row
      REAL(KIND=dp)                                      :: s(n_geom)

      r = 0.0_dp
      info = info_ok
      IF (.NOT. ALL(ieee_is_finite(q))) THEN
         info = info_degenerate
         RETURN
      END IF
      IF (m%k == 0) RETURN

      ! Check cell validity even on pure slot/linear constraints so collapsed
      ! or inverted cells cannot pass feasibility retraction silently.
      CALL cell_manifold_geometry(q, s, info)
      IF (info /= info_ok) RETURN

      DO i = 1, m%n_slot_rows
         r(i) = q(m%slot_rows(i)) - m%q0(m%slot_rows(i))
      END DO
      DO i = 1, m%n_geom_rows
         row = m%n_slot_rows + i
         r(row) = s(m%geom_rows(i)) - m%s0(m%geom_rows(i))
      END DO
      DO i = 1, m%n_eq_rows
         row = m%n_slot_rows + m%n_geom_rows + i
         r(row) = s(m%eq_pair(i, 1)) - s(m%eq_pair(i, 2))
      END DO
      DO i = 1, m%n_lin_rows
         row = m%n_slot_rows + m%n_geom_rows + m%n_eq_rows + i
         r(row) = DOT_PRODUCT(m%lin_w(:, i), q) - DOT_PRODUCT(m%lin_w(:, i), m%q0)
      END DO

   END SUBROUTINE cell_manifold_residual

! **************************************************************************************************
!> \brief A(q) (k x 6): identity rows for frozen slots, Jacobian rows for
!>        geometric labels, Jacobian differences for geometric equalities, and
!>        constant weight rows for linear couplings (in that order).
!> \param m ...
!> \param q ...
!> \param A ...
!> \param info ...
! **************************************************************************************************
   SUBROUTINE cell_manifold_active_jacobian(m, q, A, info)
      TYPE(cell_manifold_type), INTENT(IN)               :: m
      REAL(KIND=dp), INTENT(IN)                          :: q(n_slot)
      REAL(KIND=dp), INTENT(OUT)                         :: A(m%k, n_slot)
      INTEGER, INTENT(OUT)                               :: info

      INTEGER                                            :: i, row
      REAL(KIND=dp)                                      :: J(n_geom, n_slot)

      A = 0.0_dp
      info = info_ok
      IF (m%k == 0) RETURN

      CALL cell_manifold_jacobian(q, J, info)
      IF (info /= info_ok) RETURN

      DO i = 1, m%n_slot_rows
         A(i, m%slot_rows(i)) = 1.0_dp
      END DO
      DO i = 1, m%n_geom_rows
         row = m%n_slot_rows + i
         A(row, :) = J(m%geom_rows(i), :)
      END DO
      DO i = 1, m%n_eq_rows
         row = m%n_slot_rows + m%n_geom_rows + i
         A(row, :) = J(m%eq_pair(i, 1), :) - J(m%eq_pair(i, 2), :)
      END DO
      DO i = 1, m%n_lin_rows
         row = m%n_slot_rows + m%n_geom_rows + m%n_eq_rows + i
         A(row, :) = m%lin_w(:, i)
      END DO

   END SUBROUTINE cell_manifold_active_jacobian

! **************************************************************************************************
!> \brief Orthonormal null-space basis Z (6 x ndof) with A Z = 0 by a
!>        complete QR of A^T (DGEQRF + DORGQR). LICQ gate: the Gram matrix
!>        A A^T must pass a Cholesky factorization with a bounded diagonal
!>        ratio. The Cholesky fails on exact row dependence; the diagonal
!>        ratio flags near dependence that QR alone would pass (Boumal 2023,
!>        sec. 7.7).
!> \param m ...
!> \param q ...
!> \param info ...
! **************************************************************************************************
   SUBROUTINE cell_manifold_update_basis(m, q, info)
      TYPE(cell_manifold_type), INTENT(INOUT)            :: m
      REAL(KIND=dp), INTENT(IN)                          :: q(n_slot)
      INTEGER, INTENT(OUT)                               :: info

      INTEGER                                            :: i, info_lapack, lwork
      REAL(KIND=dp)                                      :: Qfull(n_slot, n_slot), rnorm, &
                                                            s_chk(n_geom), work_query(1)
      REAL(KIND=dp), ALLOCATABLE                         :: work(:)
      REAL(KIND=dp)                                      :: tau(m%k), gram(m%k, m%k), gdiag(m%k), &
                                                            AT(n_slot, m%k), A(m%k, n_slot)

      CALL cell_manifold_geometry(q, s_chk, info)
      IF (info /= info_ok) RETURN

      IF (m%k == 0) THEN
         m%ndof = n_slot
         m%Z = 0.0_dp
         DO i = 1, n_slot
            m%Z(i, i) = 1.0_dp
         END DO
         RETURN
      END IF

      CALL cell_manifold_active_jacobian(m, q, A, info)
      IF (info /= info_ok) RETURN

      ! Row scaling does not change the row space or rank: normalize before
      ! the Gram-Cholesky gate so mixed row norms (e.g. a volume gradient
      ! next to unit slot rows) can neither mask nor fake rank deficiency.
      DO i = 1, m%k
         rnorm = NORM2(A(i, :))
         IF (rnorm <= TINY(1.0_dp)) THEN
            info = info_licq
            RETURN
         END IF
         A(i, :) = A(i, :)/rnorm
      END DO

      gram = MATMUL(A, TRANSPOSE(A))
      CALL dpotrf('U', m%k, gram, m%k, info_lapack)
      IF (info_lapack /= 0) THEN
         info = info_licq
         RETURN
      END IF
      gdiag = [(gram(i, i), i=1, m%k)]
      IF (MINVAL(gdiag) <= licq_pivot_tol*MAXVAL(gdiag)) THEN
         info = info_licq
         RETURN
      END IF

      AT = TRANSPOSE(A)

      CALL dgeqrf(n_slot, m%k, AT, n_slot, tau, work_query, -1, info_lapack)
      lwork = INT(work_query(1))
      ALLOCATE (work(MAX(1, lwork)))

      CALL dgeqrf(n_slot, m%k, AT, n_slot, tau, work, SIZE(work), info_lapack)
      IF (info_lapack /= 0) THEN
         info = info_licq
         DEALLOCATE (work)
         RETURN
      END IF
      Qfull = 0.0_dp
      Qfull(1:n_slot, 1:m%k) = AT(1:n_slot, 1:m%k)
      CALL dorgqr(n_slot, n_slot, m%k, Qfull, n_slot, tau, work, SIZE(work), info_lapack)
      DEALLOCATE (work)
      IF (info_lapack /= 0) THEN
         info = info_licq
         RETURN
      END IF

      m%ndof = n_slot - m%k
      m%Z = 0.0_dp
      m%Z(:, 1:m%ndof) = Qfull(:, m%k + 1:n_slot)

   END SUBROUTINE cell_manifold_update_basis

! **************************************************************************************************
!> \brief Retraction R_q(Z du) (Absil & Malick 2012 projection-like
!>        retraction): tangent step q~ = q + Z du, then Gauss-Newton
!>        dq = -A^T (A A^T)^{-1} r(q~) until |r| <= resid_tol.
!> \param m ...
!> \param q_in ...
!> \param du ...
!> \param q_out ...
!> \param info ...
! **************************************************************************************************
   SUBROUTINE cell_manifold_retract(m, q_in, du, q_out, info)
      TYPE(cell_manifold_type), INTENT(INOUT)            :: m
      REAL(KIND=dp), INTENT(IN)                          :: q_in(n_slot), du(m%ndof)
      REAL(KIND=dp), INTENT(OUT)                         :: q_out(n_slot)
      INTEGER, INTENT(OUT)                               :: info

      INTEGER                                            :: it, lapack_info
      REAL(KIND=dp)                                      :: qt(n_slot), resid_tol, rn, &
                                                            s_chk(n_geom), zdu(n_slot), r(m%k), &
                                                            lam(m%k), chol(m%k, m%k), &
                                                            A(m%k, n_slot)

      ! Scale the residual tolerance by the geometric target magnitudes so
      ! volume and equality residuals remain compatible with double-precision
      ! roundoff on large cells.
      resid_tol = m%tol
      IF (m%n_geom_rows > 0) THEN
         resid_tol = m%tol*MAX(1.0_dp, MAXVAL(ABS(m%s0(m%geom_rows(1:m%n_geom_rows)))))
      END IF
      IF (m%n_eq_rows > 0) THEN
         DO it = 1, m%n_eq_rows
            resid_tol = MAX(resid_tol, m%tol*ABS(m%s0(m%eq_pair(it, 1))), &
                            m%tol*ABS(m%s0(m%eq_pair(it, 2))))
         END DO
      END IF

      IF (.NOT. ALL(ieee_is_finite(q_in)) .OR. .NOT. ALL(ieee_is_finite(du)) .OR. &
          ANY(ABS(du) > SQRT(HUGE(1.0_dp)))) THEN
         ! Reject non-finite steps before Z*du so 0*Inf cannot trap under
         ! -ffpe-trap=invalid.
         q_out = q_in
         info = info_no_converge
         m%last_iter = 0
         m%last_resid = SQRT(HUGE(1.0_dp))
         RETURN
      END IF

      CALL cell_manifold_expand_step(m, du, zdu)
      qt = q_in + zdu
      m%last_iter = 0

      IF (m%k == 0) THEN
         CALL cell_manifold_geometry(qt, s_chk, info)
         q_out = qt
         m%last_resid = 0.0_dp
         RETURN
      END IF

      CALL cell_manifold_residual(m, qt, r, info)
      IF (info /= info_ok) THEN
         q_out = qt
         m%last_resid = SQRT(HUGE(1.0_dp))
         RETURN
      END IF
      rn = NORM2(r)
      m%last_resid = rn
      IF (.NOT. ieee_is_finite(rn)) THEN
         info = info_no_converge
         q_out = qt
         RETURN
      END IF

      DO it = 1, m%max_iter
         IF (rn <= resid_tol) EXIT
         CALL cell_manifold_active_jacobian(m, qt, A, info)
         IF (info /= info_ok) THEN
            q_out = qt
            RETURN
         END IF
         chol = MATMUL(A, TRANSPOSE(A))
         CALL dpotrf('U', m%k, chol, m%k, lapack_info)
         IF (lapack_info /= 0) THEN
            info = info_licq
            q_out = qt
            RETURN
         END IF
         lam = r
         CALL dpotrs('U', m%k, 1, chol, m%k, lam, m%k, lapack_info)
         qt = qt - MATMUL(TRANSPOSE(A), lam)
         m%last_iter = it
         CALL cell_manifold_residual(m, qt, r, info)
         IF (info /= info_ok) THEN
            q_out = qt
            m%last_resid = SQRT(HUGE(1.0_dp))
            RETURN
         END IF
         rn = NORM2(r)
         m%last_resid = rn
         IF (.NOT. ieee_is_finite(rn)) THEN
            info = info_no_converge
            q_out = qt
            RETURN
         END IF
      END DO
      q_out = qt
      IF (rn > resid_tol) info = info_no_converge

   END SUBROUTINE cell_manifold_retract

! **************************************************************************************************
!> \brief Constraint-conjugate (projected-stress) pressure residual (#6031).
!> \param m manifold state with the basis Z at the current cell q
!> \param q current canonical cell; only the volume gradient is taken from it
!> \param g_phys raw enthalpy gradient in q-space (get_dg_dh output, i.e.
!>        -d(E + pres_ext*V)/dq; sign convention of that routine)
!> \param pdiff P_constr - pres_ext, the pressure excess the optimizer can
!>        actually relax along the constraint manifold
!> \param projected .FALSE. when the tangent space carries no volume
!>        component (volume or all diagonal slots pinned): the pressure
!>        criterion must be bypassed, like the legacy bypass for a fully
!>        fixed cell
!> \note The isotropic volume direction uv = dV/dq is projected onto the
!>       tangent space, vt = Z Z^T uv, and pdiff = -(g_phys . vt)/(vt . vt):
!>       the enthalpy derivative along the best available isotropic scaling,
!>       per unit of volume change. At a constrained stationary point the
!>       tangent gradient vanishes, so pdiff = 0 for ANY active constraint
!>       set - angle and volume constraints included - which is what the
!>       legacy masked-diagonal average lacks (it is blind to geometric
!>       constraint rows). For orthogonal cells with pure component
!>       constraints the formula reduces to the legacy free-slot average
!>       when a single diagonal slot is free; with several free slots it is
!>       their (dV/dq_i)-weighted mean, where legacy averages unweighted.
! **************************************************************************************************
   SUBROUTINE cell_manifold_proj_pressure(m, q, g_phys, pdiff, projected)
      TYPE(cell_manifold_type), INTENT(IN)               :: m
      REAL(KIND=dp), INTENT(IN)                          :: q(n_slot), g_phys(n_slot)
      REAL(KIND=dp), INTENT(OUT)                         :: pdiff
      LOGICAL, INTENT(OUT)                               :: projected

      REAL(KIND=dp)                                      :: denom, uvn2, vt(n_slot), vred(m%ndof)

      pdiff = 0.0_dp
      projected = .FALSE.
      IF (m%ndof == 0) RETURN

      ! uv = dV/dq at the current cell (V = q(1)*q(3)*q(6))
      vt = 0.0_dp
      vt(1) = q(3)*q(6)
      vt(3) = q(1)*q(6)
      vt(6) = q(1)*q(3)
      uvn2 = DOT_PRODUCT(vt, vt)

      ! vt is tangent by construction, so g_phys . vt already picks up the
      ! tangent component of the gradient
      CALL cell_manifold_project_grad(m, vt, vred)
      CALL cell_manifold_expand_step(m, vred, vt)
      denom = DOT_PRODUCT(vt, vt)
      IF (denom <= proj_pressure_rel_tol*uvn2) RETURN

      projected = .TRUE.
      ! sign flip: g_phys = -dH/dq, and positive pdiff means internal
      ! pressure above external
      pdiff = -DOT_PRODUCT(g_phys, vt)/denom

   END SUBROUTINE cell_manifold_proj_pressure

! **************************************************************************************************
!> \brief Reduced (Riemannian) gradient g_red = Z^T g_phys.
!> \param m ...
!> \param g_phys ...
!> \param g_red ...
! **************************************************************************************************
   SUBROUTINE cell_manifold_project_grad(m, g_phys, g_red)
      TYPE(cell_manifold_type), INTENT(IN)               :: m
      REAL(KIND=dp), INTENT(IN)                          :: g_phys(n_slot)
      REAL(KIND=dp), INTENT(OUT)                         :: g_red(m%ndof)

      IF (m%ndof == 0) RETURN
      g_red = MATMUL(TRANSPOSE(m%Z(:, 1:m%ndof)), g_phys)

   END SUBROUTINE cell_manifold_project_grad

! **************************************************************************************************
!> \brief Expand a reduced step to physical space: dq = Z du.
!> \param m ...
!> \param du ...
!> \param dq ...
! **************************************************************************************************
   SUBROUTINE cell_manifold_expand_step(m, du, dq)
      TYPE(cell_manifold_type), INTENT(IN)               :: m
      REAL(KIND=dp), INTENT(IN)                          :: du(m%ndof)
      REAL(KIND=dp), INTENT(OUT)                         :: dq(n_slot)

      dq = 0.0_dp
      IF (m%ndof == 0) RETURN
      dq = MATMUL(m%Z(:, 1:m%ndof), du)

   END SUBROUTINE cell_manifold_expand_step

! **************************************************************************************************
!> \brief d(angle)/dq from d cos(angle):
!>        dcos = du/(|u||v|) - cosang*(du_n/|u| + dv_n/|v|).
!> \param cosang ...
!> \param sinang ...
!> \param du ...
!> \param un ...
!> \param vn ...
!> \param dun ...
!> \param dvn ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION dangle(cosang, sinang, du, un, vn, dun, dvn) RESULT(row)
      REAL(KIND=dp), INTENT(IN)                          :: cosang, sinang, du(n_slot), un, vn, &
                                                            dun(n_slot), dvn(n_slot)
      REAL(KIND=dp)                                      :: row(n_slot)

      REAL(KIND=dp)                                      :: dcos(n_slot)

      dcos = du/(un*vn) - cosang*(dun/un + dvn/vn)
      row = -dcos/sinang

   END FUNCTION dangle

END MODULE cell_manifold_utils
