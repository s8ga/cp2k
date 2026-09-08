!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Complex block-Krylov refinement kernel for the K-point SCF diagonalization.
!>        Faithful cfm port of lanczos_refinement_2v (qs_scf_lanczos.F): the settings
!>        and the result statistics are passed explicitly instead of living in a
!>        krylov_space_type, and all scratch matrices are local to the call, so the
!>        kernel stays free of CP2K state and independently testable.
!> \par History
!>      09.2026 created mirroring lanczos_refinement_2v [QS]
! **************************************************************************************************
MODULE qs_scf_lanczos_c
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_column_scale,&
                                              cp_cfm_gemm,&
                                              cp_cfm_qr_factorization,&
                                              cp_cfm_scale_and_add,&
                                              cp_cfm_triangular_multiply,&
                                              cp_cfm_vectorsnorm
   USE cp_cfm_cholesky,                 ONLY: cp_cfm_cholesky_decompose
   USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_submatrix,&
                                              cp_cfm_release,&
                                              cp_cfm_set_all,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: z_mone,&
                                              z_one,&
                                              z_zero
   USE message_passing,                 ONLY: mp_para_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_lanczos_c'

   PUBLIC :: lanczos_refinement_2v_c

CONTAINS

! **************************************************************************************************
!> \brief One block-Krylov refinement sweep on the standardized operator ks:
!>        Rayleigh-Ritz on the incoming c0, then nkrylov-step block-Krylov
!>        expansion of the non-converged columns, closed by a full
!>        Cholesky re-orthonormalization of c0 (which is updated in place
!>        through its pointer components, mirroring the Gamma kernel).
!> \param ks standardized operator A = U^-dag H U^-1, nao x nao
!> \param c0 current MOs in the standardized basis, nao x nmo, refined in place
!> \param c1 scratch of the same shape as c0
!> \param eval eigenvalues of c0, updated in place by the entry Rayleigh-Ritz
!> \param nao ...
!> \param nmo ...
!> \param eps_iter convergence threshold on the residual column norms
!> \param nkrylov number of Krylov expansion steps per block
!> \param nblock block size of the expansion
!> \param max_res_norm largest residual column norm on exit
!> \param min_res_norm smallest residual column norm on exit
!> \param nmo_conv number of leading converged columns on exit
!> \param nmo_nc number of remaining non-converged columns on exit
!> \param check_moconv_only only run the entry Rayleigh-Ritz and the
!>        convergence count, skip the expansion
! **************************************************************************************************
   SUBROUTINE lanczos_refinement_2v_c(ks, c0, c1, eval, nao, nmo, eps_iter, nkrylov, nblock, &
                                      max_res_norm, min_res_norm, nmo_conv, nmo_nc, check_moconv_only)

      TYPE(cp_cfm_type), INTENT(IN)                      :: ks, c0, c1
      REAL(dp), DIMENSION(:), POINTER                    :: eval
      INTEGER, INTENT(IN)                                :: nao, nmo
      REAL(dp), INTENT(IN)                               :: eps_iter
      INTEGER, INTENT(IN)                                :: nkrylov, nblock
      REAL(dp), INTENT(OUT)                              :: max_res_norm, min_res_norm
      INTEGER, INTENT(OUT)                               :: nmo_conv, nmo_nc
      LOGICAL, INTENT(IN), OPTIONAL                      :: check_moconv_only

      CHARACTER(LEN=*), PARAMETER :: routineN = 'lanczos_refinement_2v_c'

      COMPLEX(dp), ALLOCATABLE, DIMENSION(:)             :: c_eval_c, work
      COMPLEX(dp), ALLOCATABLE, DIMENSION(:, :)          :: a_block, b_block, q_mat, t_mat
      COMPLEX(dp), DIMENSION(:, :), POINTER              :: a_loc, b_loc
      INTEGER :: hand1, hand2, hand3, hand4, hand5, handle, i, ia, ib, ik, imo, imo_low, imo_up, &
         info, it, j, jt, liwork, lrwork, lwork, ndim, nmo_converged, nmob, num_blocks
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: itaken, iwork
      LOGICAL                                            :: my_check_moconv_only
      REAL(dp)                                           :: max_norm, min_norm, vmax
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: c_res, rwork, t_eval
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_cfm_type)                                  :: b_mat, c2_tmp, c_tmp, chc, evec, hc, &
                                                            v_tmp
      TYPE(cp_cfm_type), ALLOCATABLE, DIMENSION(:)       :: v_mat
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      NULLIFY (fm_struct_tmp, a_loc, b_loc, blacs_env, para_env)

      para_env => c0%matrix_struct%para_env
      blacs_env => c0%matrix_struct%context

      nmo_conv = 0
      nmo_nc = 0
      max_res_norm = 0.0_dp
      min_res_norm = HUGE(1.0_dp)
      my_check_moconv_only = .FALSE.
      IF (PRESENT(check_moconv_only)) my_check_moconv_only = check_moconv_only

      NULLIFY (fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nao, ncol_global=nmo, &
                               para_env=para_env, &
                               context=blacs_env)
      CALL cp_cfm_create(c_tmp, matrix_struct=fm_struct_tmp, name="c_tmp")
      CALL cp_cfm_create(hc, matrix_struct=fm_struct_tmp, name="hc")
      CALL cp_fm_struct_release(fm_struct_tmp)

      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nmo, ncol_global=nmo, &
                               para_env=para_env, &
                               context=blacs_env)
      CALL cp_cfm_create(chc, matrix_struct=fm_struct_tmp, name="chc")
      CALL cp_cfm_create(evec, matrix_struct=fm_struct_tmp, name="evec")
      CALL cp_fm_struct_release(fm_struct_tmp)

      ! Compute (C^dag)HC
      CALL cp_cfm_gemm('N', 'N', nao, nmo, nao, z_one, ks, c0, z_zero, hc)
      CALL cp_cfm_gemm('C', 'N', nmo, nmo, nao, z_one, c0, hc, z_zero, chc)

      ! Diagonalize (C^dag)HC
      CALL timeset(routineN//"diag_chc", hand1)
      CALL cp_cfm_heevd(chc, evec, eval)
      CALL timestop(hand1)

      CALL timeset(routineN//"check_conv", hand5)
      ! Rotate the C vectors
      CALL cp_cfm_gemm('N', 'N', nao, nmo, nmo, z_one, c0, evec, z_zero, c1)

      ! Check for converged states
      CALL cp_cfm_gemm('N', 'N', nao, nmo, nmo, z_one, hc, evec, z_zero, c_tmp)
      CALL cp_cfm_to_cfm(c1, c0, nmo, 1, 1)
      ALLOCATE (c_eval_c(nmo))
      c_eval_c = CMPLX(eval(1:nmo), 0.0_dp, KIND=dp)
      CALL cp_cfm_column_scale(c1, c_eval_c)
      DEALLOCATE (c_eval_c)
      CALL cp_cfm_scale_and_add(z_one, c_tmp, z_mone, c1)
      ALLOCATE (c_res(nmo))
      CALL cp_cfm_vectorsnorm(c_tmp, c_res)

      nmo_converged = 0
      nmo_nc = 0
      max_norm = 0.0_dp
      min_norm = 1.e10_dp
      CALL cp_cfm_set_all(c1, z_zero)
      DO imo = 1, nmo
         max_norm = MAX(max_norm, c_res(imo))
         min_norm = MIN(min_norm, c_res(imo))
      END DO
      DO imo = 1, nmo
         IF (c_res(imo) <= eps_iter) THEN
            nmo_converged = nmo_converged + 1
         ELSE
            nmo_nc = nmo - nmo_converged
            EXIT
         END IF
      END DO
      CALL timestop(hand5)

      CALL cp_cfm_release(c_tmp)
      CALL cp_cfm_release(hc)

      nmo_conv = nmo_converged
      nmo_nc = nmo_nc
      max_res_norm = max_norm
      min_res_norm = min_norm

      IF (my_check_moconv_only) THEN
         ! Do nothing
      ELSE IF (nmo_nc > 0) THEN

         CALL cp_cfm_to_cfm(c0, c1, nmo_nc, nmo_converged + 1, 1)

         IF (MODULO(nmo_nc, nblock) > 0) THEN
            num_blocks = nmo_nc/nblock + 1
         ELSE
            num_blocks = nmo_nc/nblock
         END IF

         DO ib = 1, num_blocks

            imo_low = (ib - 1)*nblock + 1
            imo_up = MIN(ib*nblock, nmo_nc)
            nmob = imo_up - imo_low + 1
            ndim = nkrylov*nmob

            ! Allocation
            CALL timeset(routineN//"alloc", hand5)
            ALLOCATE (a_block(nmob, nmob))
            ALLOCATE (b_block(nmob, nmob))
            ALLOCATE (t_mat(ndim, ndim))
            ALLOCATE (t_eval(ndim))

            NULLIFY (fm_struct_tmp)
            ! by forcing ncol_block=nmob, the needed part of the matrix is distributed on a subset of processes
            ! this is due to the use of two-dimensional grids of processes
            ! nrow_global is distributed over num_pe(1)
            ! a local_data array is anyway allocated for the processes non included
            ! this should have a minimum size
            ! with ncol_local=1
            CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nao, ncol_global=nmob, &
                                     ncol_block=nmob, &
                                     para_env=para_env, &
                                     context=blacs_env, &
                                     force_block=.TRUE.)
            CALL cp_cfm_create(c_tmp, matrix_struct=fm_struct_tmp, name="c_tmp")
            CALL cp_cfm_set_all(c_tmp, z_zero)
            CALL cp_cfm_create(v_tmp, matrix_struct=fm_struct_tmp, name="v_tmp")
            CALL cp_cfm_set_all(v_tmp, z_zero)
            CALL cp_fm_struct_release(fm_struct_tmp)
            NULLIFY (fm_struct_tmp)
            ALLOCATE (v_mat(nkrylov))
            CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nao, ncol_global=nmob, &
                                     ncol_block=nmob, &
                                     para_env=para_env, &
                                     context=blacs_env, &
                                     force_block=.TRUE.)
            DO ik = 1, nkrylov
               CALL cp_cfm_create(v_mat(ik), matrix_struct=fm_struct_tmp, name="v_mat")
            END DO
            CALL cp_cfm_create(hc, matrix_struct=fm_struct_tmp, name="hc")
            CALL cp_fm_struct_release(fm_struct_tmp)
            NULLIFY (fm_struct_tmp)
            CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nao, ncol_global=ndim, &
                                     ncol_block=ndim, &
                                     para_env=para_env, &
                                     context=blacs_env, &
                                     force_block=.TRUE.)
            CALL cp_cfm_create(c2_tmp, matrix_struct=fm_struct_tmp, name="c2_tmp")
            CALL cp_fm_struct_release(fm_struct_tmp)

            NULLIFY (fm_struct_tmp)
            CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nmob, ncol_global=nmob, &
                                     para_env=para_env, &
                                     context=blacs_env)
            CALL cp_cfm_create(b_mat, matrix_struct=fm_struct_tmp, name="b_mat")
            CALL cp_fm_struct_release(fm_struct_tmp)
            CALL timestop(hand5)
            ! End allocation

            CALL cp_cfm_set_all(b_mat, z_zero)
            CALL cp_cfm_to_cfm(c1, v_mat(1), nmob, imo_low, 1)

            ! Here starts the construction of krylov space
            CALL timeset(routineN//"_kry_loop", hand2)
            ! Compute A = (V^dag)HV
            CALL cp_cfm_gemm('N', 'N', nao, nmob, nao, z_one, ks, v_mat(1), z_zero, hc)

            a_block = z_zero
            a_loc => v_mat(1)%local_data
            b_loc => hc%local_data

            IF (SIZE(hc%local_data, 2) == nmob) THEN
               ! this is a work around to avoid problems due to the two dimensional grid of processes
               CALL zgemm('C', 'N', nmob, nmob, SIZE(hc%local_data, 1), z_one, a_loc(1, 1), &
                          SIZE(hc%local_data, 1), b_loc(1, 1), SIZE(hc%local_data, 1), z_zero, a_block(1, 1), nmob)
            END IF
            CALL para_env%sum(a_block)

            ! Compute the residual matrix R for next factorisation
            c_tmp%local_data = z_zero
            IF (SIZE(c_tmp%local_data, 2) == nmob) THEN
               b_loc => c_tmp%local_data
               CALL zgemm('N', 'N', SIZE(c_tmp%local_data, 1), nmob, nmob, z_one, a_loc(1, 1), &
                          SIZE(c_tmp%local_data, 1), a_block(1, 1), nmob, z_zero, &
                          b_loc(1, 1), SIZE(c_tmp%local_data, 1))
            END IF
            CALL cp_cfm_scale_and_add(z_mone, c_tmp, z_one, hc)

            ! deflate the search direction against the locked columns of c0
            IF (nmo_converged > 0) THEN
               CALL deflate_locked(c0, nmo_converged, c_tmp, nao, nmob)
            END IF

            ! Build the block tridiagonal matrix
            t_mat = z_zero
            DO i = 1, nmob
               t_mat(1:nmob, i) = a_block(1:nmob, i)
            END DO

            DO ik = 2, nkrylov
               ! Call lapack for QR factorization
               CALL cp_cfm_set_all(b_mat, z_zero)
               CALL cp_cfm_to_cfm(c_tmp, v_mat(ik), nmob, 1, 1)
               CALL cp_cfm_qr_factorization(c_tmp, b_mat, nao, nmob, 1, 1)

               CALL cp_cfm_triangular_multiply(b_mat, v_mat(ik), side="R", invert_tr=.TRUE., &
                                               n_rows=nao, n_cols=nmob)

               CALL cp_cfm_get_submatrix(b_mat, b_block, 1, 1, nmob, nmob)

               ! Compute A = (V^dag)HV
               CALL cp_cfm_gemm('N', 'N', nao, nmob, nao, z_one, ks, v_mat(ik), z_zero, hc)

               a_block = z_zero
               IF (SIZE(hc%local_data, 2) == nmob) THEN
                  a_loc => v_mat(ik)%local_data
                  b_loc => hc%local_data
                  CALL zgemm('C', 'N', nmob, nmob, SIZE(hc%local_data, 1), z_one, a_loc(1, 1), &
                             SIZE(hc%local_data, 1), b_loc(1, 1), SIZE(hc%local_data, 1), z_zero, a_block(1, 1), nmob)
               END IF
               CALL para_env%sum(a_block)

               ! Compute the residual matrix R for next factorisation
               c_tmp%local_data = z_zero
               IF (SIZE(c_tmp%local_data, 2) == nmob) THEN
                  a_loc => v_mat(ik)%local_data
                  b_loc => c_tmp%local_data
                  CALL zgemm('N', 'N', SIZE(c_tmp%local_data, 1), nmob, nmob, z_one, a_loc(1, 1), &
                             SIZE(c_tmp%local_data, 1), a_block(1, 1), nmob, z_zero, &
                             b_loc(1, 1), SIZE(c_tmp%local_data, 1))
               END IF
               CALL cp_cfm_scale_and_add(z_mone, c_tmp, z_one, hc)

               IF (SIZE(c_tmp%local_data, 2) == nmob) THEN
                  a_loc => v_mat(ik - 1)%local_data
                  DO j = 1, nmob
                     DO i = 1, j
                        DO ia = 1, SIZE(c_tmp%local_data, 1)
                           ! conjugated coupling, see the block tridiagonal
                           ! assembly below
                           b_loc(ia, i) = b_loc(ia, i) - a_loc(ia, j)*CONJG(b_block(i, j))
                        END DO
                     END DO
                  END DO
               END IF

               ! deflate the search direction against the locked columns of c0
               IF (nmo_converged > 0) THEN
                  CALL deflate_locked(c0, nmo_converged, c_tmp, nao, nmob)
               END IF

               ! Build the block tridiagonal matrix
               it = (ik - 2)*nmob
               jt = (ik - 1)*nmob
               DO j = 1, nmob
                  t_mat(jt + 1:jt + nmob, jt + j) = a_block(1:nmob, j)
                  DO i = 1, nmob
                     ! the coupling between V_ik-1 and V_ik enters the Hermitian
                     ! block tridiagonal as C^dag in the upper block; the Gamma
                     ! kernel mirrors b_block without conjugation, which is only
                     ! equivalent for real data
                     t_mat(it + i, jt + j) = CONJG(b_block(j, i))
                     t_mat(jt + j, it + i) = b_block(j, i)
                  END DO
               END DO
            END DO ! ik
            CALL timestop(hand2)

            CALL timeset(routineN//"_diag_tri", hand3)
            ! Diagonalize the block-tridiagonal matrix, workspace by the LAPACK query
            lwork = -1
            lrwork = -1
            liwork = -1
            ALLOCATE (work(2*ndim), rwork(2*ndim), iwork(2*ndim))
            CALL zheevd('V', 'U', ndim, t_mat(1, 1), ndim, t_eval(1), &
                        work(1), lwork, rwork(1), lrwork, iwork(1), liwork, info)
            lwork = MAX(INT(REAL(work(1))), 2*ndim)
            lrwork = MAX(INT(rwork(1)), 2*ndim)
            liwork = MAX(iwork(1), 2*ndim)
            DEALLOCATE (work, rwork, iwork)
            ALLOCATE (work(lwork), rwork(lrwork), iwork(liwork))
            CALL zheevd('V', 'U', ndim, t_mat(1, 1), ndim, t_eval(1), &
                        work(1), lwork, rwork(1), lrwork, iwork(1), liwork, info)
            DEALLOCATE (work, rwork, iwork)
            CALL timestop(hand3)

            CALL timeset(routineN//"_build_cnew", hand4)
            ! Compute the refined vectors

            c2_tmp%local_data = z_zero
            ALLOCATE (q_mat(nmob, ndim))
            q_mat = z_zero
            b_loc => c2_tmp%local_data
            DO ik = 1, nkrylov
               CALL cp_cfm_to_cfm(v_mat(ik), v_tmp, nmob, 1, 1)
               IF (SIZE(c2_tmp%local_data, 2) == ndim) THEN
                  a_loc => v_tmp%local_data
                  it = (ik - 1)*nmob
                  CALL zgemm('N', 'N', SIZE(c2_tmp%local_data, 1), ndim, nmob, z_one, a_loc(1, 1), &
                             SIZE(c2_tmp%local_data, 1), t_mat(it + 1, 1), ndim, z_one, &
                             b_loc(1, 1), SIZE(c2_tmp%local_data, 1))
               END IF
            END DO !ik

            ! Try to avoid linear dependencies
            CALL cp_cfm_to_cfm(v_mat(1), v_tmp, nmob, 1, 1)
            IF (SIZE(c2_tmp%local_data, 2) == ndim) THEN
               a_loc => v_tmp%local_data
               b_loc => c2_tmp%local_data
               CALL zgemm('C', 'N', nmob, ndim, SIZE(v_tmp%local_data, 1), z_one, a_loc(1, 1), &
                          SIZE(v_tmp%local_data, 1), b_loc(1, 1), SIZE(v_tmp%local_data, 1), &
                          z_zero, q_mat(1, 1), nmob)
            END IF
            CALL para_env%sum(q_mat)

            ALLOCATE (itaken(ndim))
            itaken = 0
            DO it = 1, nmob
               vmax = 0.0_dp
               ! select index ik
               DO jt = 1, ndim
                  IF (itaken(jt) == 0 .AND. ABS(q_mat(it, jt)) > vmax) THEN
                     vmax = ABS(q_mat(it, jt))
                     ik = jt
                  END IF
               END DO
               itaken(ik) = 1

               CALL cp_cfm_to_cfm(c2_tmp, v_mat(1), 1, ik, it)
            END DO
            DEALLOCATE (itaken)
            DEALLOCATE (q_mat)

            ! Copy in the converged set to enlarge the converged subspace
            CALL cp_cfm_to_cfm(v_mat(1), c0, nmob, 1, (nmo_converged + imo_low))
            CALL timestop(hand4)

            CALL cp_cfm_release(c2_tmp)
            CALL cp_cfm_release(c_tmp)
            CALL cp_cfm_release(hc)
            CALL cp_cfm_release(v_tmp)
            CALL cp_cfm_release(b_mat)

            DEALLOCATE (t_mat)
            DEALLOCATE (a_block)
            DEALLOCATE (b_block)
            DEALLOCATE (t_eval)

            DO ik = 1, nkrylov
               CALL cp_cfm_release(v_mat(ik))
            END DO
            DEALLOCATE (v_mat)

         END DO ! ib

         CALL timeset(routineN//"_ortho", hand1)
         CALL cp_cfm_gemm('C', 'N', nmo, nmo, nao, z_one, c0, c0, z_zero, chc)

         CALL cp_cfm_cholesky_decompose(chc, nmo)
         CALL cp_cfm_triangular_multiply(chc, c0, 'R', invert_tr=.TRUE.)
         CALL timestop(hand1)
      ELSE
         ! Do nothing
      END IF

      CALL cp_cfm_release(chc)
      CALL cp_cfm_release(evec)
      DEALLOCATE (c_res)

      CALL timestop(handle)
   END SUBROUTINE lanczos_refinement_2v_c

! **************************************************************************************************
!> \brief Deflate a search direction against the locked converged columns of
!>        c0: direction = (1 - Q Q^dag) direction with Q = c0(:, 1:nlocked).
!>        The block recurrence re-acquires components along converged
!>        eigenvectors, the direction into which Lanczos orthogonality is
!>        lost first, and the closing Cholesky rotation would smear those
!>        components back into the refined tail. The projection mirrors the
!>        locked-space deflation of blocked production solvers and costs two
!>        small gemms per expansion step.
!> \param c0 standardized MOS, orthonormal, its first nlocked columns converged
!> \param nlocked number of leading converged columns of c0
!> \param direction search direction, nao x nmob, deflated in place
!> \param nao ...
!> \param nmob ...
! **************************************************************************************************
   SUBROUTINE deflate_locked(c0, nlocked, direction, nao, nmob)

      TYPE(cp_cfm_type), INTENT(IN)                      :: c0
      INTEGER, INTENT(IN)                                :: nlocked
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: direction
      INTEGER, INTENT(IN)                                :: nao, nmob

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'deflate_locked'

      INTEGER                                            :: handle
      TYPE(cp_cfm_type)                                  :: proj
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp

      CALL timeset(routineN, handle)

      NULLIFY (fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nlocked, ncol_global=nmob, &
                               ncol_block=nmob, &
                               para_env=direction%matrix_struct%para_env, &
                               context=direction%matrix_struct%context, &
                               force_block=.TRUE.)
      CALL cp_cfm_create(proj, matrix_struct=fm_struct_tmp, name="defl_proj")
      CALL cp_fm_struct_release(fm_struct_tmp)

      CALL cp_cfm_gemm('C', 'N', nlocked, nmob, nao, z_one, c0, direction, z_zero, proj)
      CALL cp_cfm_gemm('N', 'N', nao, nmob, nlocked, z_mone, c0, proj, z_one, direction)

      CALL cp_cfm_release(proj)

      CALL timestop(handle)

   END SUBROUTINE deflate_locked

END MODULE qs_scf_lanczos_c
