!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Pair-domain PARI contraction for k-point RI-HFX (two-center localized RI,
!>        JPCL 2020, 10.1021/acs.jpclett.0c00481).
!>        Growing twin of hfx_ri_kp.F: this module hosts the per-pair machinery
!>        (pair domains, per-pair metric pseudo-inverse, pooled fit coefficients),
!>        while the legacy global path in hfx_ri_kp.F remains the default.
!>
!>        Storage layout of the pre-SCF layer (geometry-dependent, density-independent):
!>        - pari_mat_2c(b): replicated per-image (P_i^0 | Q_j^b), block (i,j) of image b
!>        - pari_c3: pooled C3(i,j,b) coefficients on aux(i) rows for products of
!>          own atom i and partner atom j^b, i.e. phi_mu(i^0) phi_nu(j^b) fitted over
!>          the pair domain aux(i^0) u aux(j^b), keeping only the aux(i^0) rows.
!>          The raw 3c integrals are folded into C3 and discarded.
!> \history 10.2026 created
! **************************************************************************************************
MODULE hfx_ri_kp_pari

   USE basis_set_types,                 ONLY: gto_basis_set_p_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_copy, dbcsr_create, dbcsr_distribution_release, dbcsr_distribution_type, &
        dbcsr_get_block_p, dbcsr_get_stored_coordinates, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_p_type, dbcsr_put_block, dbcsr_release, dbcsr_replicate_all, dbcsr_type, &
        dbcsr_type_no_symmetry
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_dist2d_to_dist
   USE dbt_api,                         ONLY: &
        dbt_create, dbt_destroy, dbt_distribution_destroy, dbt_distribution_new, &
        dbt_distribution_type, dbt_filter, dbt_get_block, dbt_mp_environ_pgrid, dbt_pgrid_create, &
        dbt_pgrid_destroy, dbt_pgrid_type, dbt_type
   USE distribution_2d_types,           ONLY: distribution_2d_type
   USE hfx_ri_kp,                       ONLY: get_kp_and_ri_images
   USE hfx_types,                       ONLY: hfx_ri_type
   USE input_cp2k_hfx,                  ONLY: ri_pmat
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE message_passing,                 ONLY: mp_cart_type,&
                                              mp_para_env_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
   USE particle_types, ONLY: particle_type
   USE qs_environment_types, ONLY: get_qs_env, qs_environment_type
   USE qs_integral_utils, ONLY: basis_set_list_setup
   USE qs_kind_types, ONLY: qs_kind_type
   USE qs_neighbor_list_types, ONLY: neighbor_list_set_p_type, release_neighbor_list_sets
   USE qs_tensors, ONLY: build_2c_integrals, build_2c_neighbor_lists, &
                         build_3c_integrals, build_3c_neighbor_lists, &
                         neighbor_list_3c_destroy
   USE qs_tensors_types, ONLY: create_3c_tensor, distribution_3d_create, &
                               distribution_3d_type, neighbor_list_3c_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_ri_kp_pari'

   ! relative regularization strength of the per-pair LM pseudo-inverse
   ! (mu = eps_eig * w_max). Provenance: the appendix-10 insensitivity
   ! calibration (at 1e-10, then in truncated form) predates the
   ! storage-parity fix and was overturned by the regularization-path A/B;
   ! 3e-3 is that A/B's choice -- a dedicated sensitivity scan is pending
   REAL(KIND=dp), PRIVATE :: pari_eps_eig = 3.0E-3_dp
   ! TODO(dev experiment, remove before PR): per-term ablation switch (bit 1=T1, 2=T2, 4=T3, 8=T4)
   INTEGER, PRIVATE :: pari_term_mask = 15
   ! TODO(dev experiment, remove before PR): image-index L1 cutoff knobs. The
   ! round-16 verdict: gating by image INDEX (not by physical content) breaks
   ! the Gram structure of the exchange operator -- truncating V blocks or
   ! gating T4 alone loses positive semidefiniteness and the T1<->T4 adjoint
   ! pair, which amplified a 3% step-1 error into a 6000x step-2 explosion on
   ! hBN4 (2D diagonal images at 3.5 A carry live content at L1=2). Defaults
   ! below are OFF (=99, no image is cut); writing 2/3/3 into the knob files
   ! restores the legacy gated behavior for A/B attribution and cost runs
   INTEGER, PRIVATE :: pari_remote_cut = 99
   INTEGER, PRIVATE :: pari_vgate = 99
   INTEGER, PRIVATE :: pari_hard_gate = 99
   REAL(KIND=dp), PRIVATE :: pari_sac_f = 1.0_dp
   ! TODO(dev experiment, remove before PR): aux-cell placement mode. 1 = extended
   ! (per-image RI blocks from get_kp_and_ri_images; the home block index stays
   ! katom, so aliased neighbor-aux contributions land elsewhere and are never
   ! read). 0 = legacy, every aux image folded into the main-cell block, which
   ! build_3c_integrals SUMS into that block (summation=.TRUE.)
   INTEGER, PRIVATE :: pari_ri_mode = 1
   ! TODO(dev diagnostic, remove before PR): per-term F images and far-chain
   ! source accumulators; off by default -- their GEMMs double every scatter

   PUBLIC :: hfx_ri_update_ks_kp_pari, pari_pair_pinv

CONTAINS

! **************************************************************************************************
!> \brief Levenberg-Marquardt regularized inverse of a dense symmetric
!>        pair-domain RI metric: w_inv = w/(w^2 + mu^2) with mu = eps_eig*w_max.
!>        Genuine modes (w >> mu) stay unbiased, near-null and noise modes damp
!>        smoothly to zero, and the global amplification is capped at 1/(2 mu);
!>        the response is differentiable, so the SCF mixer never sees
!>        eigenvalue-threshold jumps (unlike a hard cutoff). The price: modes
!>        near w ~ mu are attenuated ~50 percent, the fit is biased at that
!>        scale, and the exchange-energy error is first order in that bias
!>        (measured fit residuals 1e-3..1e-1 at the default -- see the probe in
!>        pari_solve_pairs); recalibrate eps_eig for a new system class.
!>        A direct Cholesky inverse is deliberately not attempted: the PARI
!>        union basis is intentionally redundant (pair-domain condition numbers
!>        of order 1e10 measured for Si), where a direct inverse deviates O(1)
!>        from the regularized solve in the near-null directions, while the
!>        contraction coefficients stay stable only under regularization.
!> \param vmat dense symmetric pair-domain metric (P|Q), untouched on exit
!> \param eps_eig relative regularization strength (default 3.0E-3, knob file PARI_EPS_EIG)
!> \param vinv regularized inverse, allocated by the caller (shape of vmat)
!> \param rcond reciprocal condition estimate 1/(||V||_1 ||V^+||_1), for logging
! **************************************************************************************************
   SUBROUTINE pari_pair_pinv(vmat, eps_eig, vinv, rcond)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: vmat
      REAL(KIND=dp), INTENT(IN)                          :: eps_eig
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: vinv
      REAL(KIND=dp), INTENT(OUT)                         :: rcond

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'pari_pair_pinv'

      INTEGER                                            :: handle, i, info, lwork, n
      REAL(KIND=dp)                                      :: a_norm, ainv_norm, mu, w_inv
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: w, work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: u_scaled, vwork
      REAL(KIND=dp), EXTERNAL                            :: dlange

      CALL timeset(routineN, handle)

      n = SIZE(vmat, 1)
      CPASSERT(SIZE(vmat, 2) == n)
      CPASSERT(SIZE(vinv, 1) == n .AND. SIZE(vinv, 2) == n)
      CPASSERT(eps_eig > 0.0_dp)

      ! Levenberg-Marquardt soft filtering of the pair metric inverse:
      ! w_inv = w / (w^2 + mu^2) with mu = eps_eig * w_max. Genuine modes
      ! (w >> mu) stay unbiased, near-null and noise modes damp smoothly to
      ! zero, and the global amplification is capped at 1/(2 mu); unlike a
      ! hard cutoff the response is differentiable, so the SCF mixer never
      ! sees eigenvalue threshold jumps
      ALLOCATE (vwork(n, n), work(3*n - 1), w(n), u_scaled(n, n))
      vwork(:, :) = vmat
      lwork = SIZE(work)
      CALL dsyev('V', 'U', n, vwork, n, w, work, lwork, info)
      CPASSERT(info == 0)

      mu = eps_eig*MAXVAL(w)
      DO i = 1, n
         IF (w(i) > 0.0_dp) THEN
            w_inv = w(i)/(w(i)*w(i) + mu*mu)
         ELSE
            w_inv = 0.0_dp
         END IF
         u_scaled(:, i) = w_inv*vwork(:, i)
      END DO
      CALL dgemm('N', 'T', n, n, n, 1.0_dp, u_scaled, n, vwork, n, 0.0_dp, vinv, n)

      ! condition estimate from the 1-norms of V and V^+ (the filtered spectrum
      ! makes dpocon-style bounds meaningless, so a norm product is used instead)
      a_norm = dlange('1', n, n, vmat, n, work)
      ainv_norm = dlange('1', n, n, vinv, n, work)
      rcond = 1.0_dp/(a_norm*ainv_norm)

      DEALLOCATE (vwork, work, w, u_scaled)

      CALL timestop(handle)

   END SUBROUTINE pari_pair_pinv

! **************************************************************************************************
!> \brief K-point RI-HFX exchange entry with pair atomic domain RI (PARI).
!>        Drop-in twin of hfx_ri_update_ks_kp (hfx_ri_kp.F), dispatched from
!>        hfx_admm_utils on RI%KP_METHOD PARI. The pre-SCF layer builds and caches
!>        the density-independent pair fit coefficients; the SCF-side four-term
!>        KS assembly lands in the follow-up (currently zero exchange + warning).
!> \param qs_env ...
!> \param ri_data ...
!> \param ks_matrix ...
!> \param ehfx ...
!> \param rho_ao ...
!> \param geometry_did_change ...
!> \param nspins ...
!> \param hf_fraction ...
! **************************************************************************************************
   SUBROUTINE hfx_ri_update_ks_kp_pari(qs_env, ri_data, ks_matrix, ehfx, rho_ao, &
                                       geometry_did_change, nspins, hf_fraction)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      REAL(KIND=dp), INTENT(OUT)                         :: ehfx
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao
      LOGICAL, INTENT(IN)                                :: geometry_did_change
      INTEGER, INTENT(IN)                                :: nspins
      REAL(KIND=dp), INTENT(IN)                          :: hf_fraction

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ri_update_ks_kp_pari'

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: fac

      CPASSERT(nspins > 0)
      CPASSERT(hf_fraction >= 0.0_dp)

      CALL timeset(routineN, handle)

      IF (geometry_did_change) THEN
         CALL hfx_ri_pre_scf_kp_pari(ks_matrix(1, 1)%matrix, ri_data, qs_env)
      END IF

      IF (nspins == 1) THEN
         fac = 0.5_dp*hf_fraction
      ELSE
         fac = 1.0_dp*hf_fraction
      END IF

      CALL pari_assemble_ks(ri_data, qs_env, rho_ao, ks_matrix, ehfx, fac)

      CALL timestop(handle)

   END SUBROUTINE hfx_ri_update_ks_kp_pari

! **************************************************************************************************
!> \brief Four-term PARI KS assembly (k_paper / JPCL eq 16 structure), dense
!>        image-general version: produces the KS matrix of EVERY periodic image
!>        F(mu^0, nu^b), b = 1..nimg; the k-point phases are applied downstream
!>        by the k-point machinery (the image basis carries them implicitly).
!>        The same-site 0.5 weight is folded into the cached C3 coefficients
!>        (pair (i,i,main)).
!>
!>        Periodic kernel per output (mu on atom i, cell 0; nu on atom j, cell b),
!>        with the density P(s^x, l^y) and the pair fits C3[X,Y,y] (fit of the
!>        product mu_X^0 nu_Y^y, rows on X; partner rows reachable as
!>        C3[Y,X,-y] by the (mu nu|P) pair symmetry):
!>          Ma=i:  T1  F(i,j,b) += BD_i(mu, la(l,y)) . V^y[i,l] . C3(l, j, b-y)
!>                 T2  F(i,j,b) += BD_i(mu, la(l,y)) . V^b[i,j] . C3(j, l, y-b)
!>          Ma=k:  T3  F(i,j,b) += BD2_k^{(x)}(la(l,y), mu) . V^{y-x}[k,l] . C3(l, j, b-y)
!>                 T4  F(i,j,b) += BD2_k^{(x)}(la(l,y), mu) . V^{b-x}[k,j] . C3(j, l, y-b)
!>        where BD_i(p, mu, la^y) = sum_{k,x,s} C3[i,k,x](p,mu,s) P(s^x, la^y)
!>              BD2_k^{(x)}(p, la^y, mu) = sum_s C3[k,i,-x](p,s,mu) P(s^x, la^y)
!>        and V^z[A,B] = (P_A^0|Q_B^z) block of the per-image 2c matrix.
!>        All four reduce to the molecular main-image recipe when only the main
!>        image is populated (verified against the pyscf oracle, gate G1/G2).
!>        Structural identity: F^b = (F^{-b})^T (asserted, warning if violated).
!> \param ri_data ...
!> \param qs_env ...
!> \param rho_ao density matrix images (symmetric-typed, P_ab^i = P_ba^-i)
!> \param ks_matrix KS matrix images (accumulated with -fac, legacy convention)
!> \param ehfx exchange energy 0.5 sum_b Tr((-fac F^b) P^b)
!> \param fac spin-reduced exchange prefactor
! **************************************************************************************************
   SUBROUTINE pari_assemble_ks(ri_data, qs_env, rho_ao, ks_matrix, ehfx, fac)
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao, ks_matrix
      REAL(KIND=dp), INTENT(OUT)                         :: ehfx
      REAL(KIND=dp), INTENT(IN)                          :: fac

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'pari_assemble_ks'

      CHARACTER(LEN=64)                                  :: progbuf
      INTEGER                                            :: ienv
      LOGICAL                                            :: pari_diag, qexist
      REAL(KIND=dp)                                      :: t_stage

      EXTERNAL :: dgemm
      INTEGER :: Aat, Bat, b_img, handle, htot, i, ie, imu, inu, ip, ithr, iu, iy, iz, j, jb, k, l, lam, &
                 lc, mu, nAOi, nAOj, nAOk, nAOl, gstep, nimg, naomax, nlive, npd, nRIi, &
                 ki_tot, nRIj, nRIk, nRIl, nrimax, nao, natom, npct, nthr, nu, nyl, o0, src, &
                 offI, offJ, offK, offL, lpct, unit_dump, u_img, v_img, w_img, x_img, &
                 y_img, z_img, z1, z2
      INTEGER, ALLOCATABLE, DIMENSION(:)                :: glen, slg, live_ip, live_y, mimg, pd_list, &
                                                           pd_rev, y_list
      INTEGER, ALLOCATABLE, DIMENSION(:, :)             :: hoff, imgcell, live_fill, live_off, reldiff
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)          :: pair3, pair3s
      INTEGER, ALLOCATABLE, DIMENSION(:, :)             :: live_cnt
      LOGICAL, ALLOCATABLE, DIMENSION(:)                :: ylive
      LOGICAL, ALLOCATABLE, DIMENSION(:, :, :)          :: liveH
      INTEGER, DIMENSION(3)                             :: dc
      INTEGER, DIMENSION(:, :), POINTER                 :: index_to_cell
      INTEGER, DIMENSION(:, :, :), POINTER              :: cell_to_index
      TYPE(mp_para_env_type), POINTER                   :: para_env
      REAL(KIND=dp)                                     :: d1y, d1z, d4e, d4v, d4w, fd_delta, fd_x, half, trb, worst_sym
      LOGICAL                                            :: lpair
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)          :: c3norm, c3norms, Gflat, Gj, Hflat, pdiff_norm
      REAL(KIND=dp)                                     :: cn_thr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)       :: ptmp, vblk
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)    :: bd2_blk, fimg, pdiff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :) :: fimg_t
      ! TODO(dev diagnostic, remove before PR): far-chain source accumulators
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :) :: f1src_y, f1src_z, f4src_e, f4src_v, f4src_w
      ! MPI: replicated density images + KS block owner
      TYPE(dbcsr_p_type), ALLOCATABLE, DIMENSION(:) :: rho_rep
      INTEGER :: owner
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :) :: bd_all
      ! bd2: (p, lam_global, mu, iz) -- the per-l planes are one nao-wide slab so
      ! the T4 fill is ONE dgemm per (k, x, z, mu) against the whole P row block
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :) :: bd2buf
      ! b7c stacks (T4 iz / T2 y): per-group live lists, packed A, stacked W
      INTEGER :: iw, lwi, nwl, nzl, izl, rr, nyl_l
      INTEGER, ALLOCATABLE, DIMENSION(:) :: w_iz, w_im, w_off, izl_of, yl_y, yl_off
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: Ap, Wt, Ap2, Wt2
      REAL(KIND=dp), DIMENSION(:, :), POINTER           :: blk
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(kpoint_type), POINTER                         :: kpoints

      NULLIFY (kpoints, cell_to_index, index_to_cell, blk, para_env)
      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, kpoints=kpoints, para_env=para_env)
      CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index, index_to_cell=index_to_cell)

      ! TODO(dev experiment, remove before PR): per-term ablation switch,
      ! selected by a PARI_TERM_MASK file in the run directory (bit 1=T1, 2=T2, 4=T3, 8=T4)
      INQUIRE (FILE="PARI_TERM_MASK", EXIST=qexist)
      IF (qexist) THEN
         OPEN (NEWUNIT=ienv, FILE="PARI_TERM_MASK", status="old", action="read")
         READ (ienv, *) pari_term_mask
         CLOSE (ienv)
      END IF
      INQUIRE (FILE="PARI_VGATE", EXIST=qexist)
      IF (qexist) THEN
         OPEN (NEWUNIT=ienv, FILE="PARI_VGATE", status="old", action="read")
         READ (ienv, *) pari_vgate
         CLOSE (ienv)
      END IF
      ! TODO(dev experiment, remove before PR): the five legacy hardcoded
      ! image-index gates of the four-term assembly (defaults off, see above)
      INQUIRE (FILE="PARI_HARD_GATE", EXIST=qexist)
      IF (qexist) THEN
         OPEN (NEWUNIT=ienv, FILE="PARI_HARD_GATE", status="old", action="read")
         READ (ienv, *) pari_hard_gate
         CLOSE (ienv)
      END IF
      ! TODO(dev diagnostic, remove before PR): opt-in diagnostics
      ! (single-rank only: the per-term accumulators live inside the
      ! rank-split atom loop and would be partial otherwise)
      ! per-call evaluation (a module variable would stay sticky across qs_envs)
      INQUIRE (FILE="PARI_KEEP_DIAG", EXIST=qexist)
      pari_diag = qexist .AND. para_env%num_pe == 1

      CPASSERT(SIZE(rho_ao, 1) == 1) ! spin-generalisation is a follow-up

      natom = SIZE(ri_data%bsizes_AO)
      nao = SUM(ri_data%bsizes_AO)
      nimg = ri_data%pari_nimg

      ! image cells, opposite-image and relative-image tables
      ALLOCATE (mimg(nimg), reldiff(nimg, nimg))
      mimg = 0
      reldiff = 0
      DO b_img = 1, nimg
         dc = -index_to_cell(:, b_img)
         IF (dc(1) >= LBOUND(cell_to_index, 1) .AND. dc(1) <= UBOUND(cell_to_index, 1) .AND. &
             dc(2) >= LBOUND(cell_to_index, 2) .AND. dc(2) <= UBOUND(cell_to_index, 2) .AND. &
             dc(3) >= LBOUND(cell_to_index, 3) .AND. dc(3) <= UBOUND(cell_to_index, 3)) THEN
            ! the k-point cell list may extend beyond the HFX image range
            IF (cell_to_index(dc(1), dc(2), dc(3)) <= nimg) THEN
               mimg(b_img) = cell_to_index(dc(1), dc(2), dc(3))
            END IF
         END IF
      END DO
      DO b_img = 1, nimg
         DO y_img = 1, nimg
            dc = index_to_cell(:, b_img) - index_to_cell(:, y_img)
            IF (dc(1) >= LBOUND(cell_to_index, 1) .AND. dc(1) <= UBOUND(cell_to_index, 1) .AND. &
                dc(2) >= LBOUND(cell_to_index, 2) .AND. dc(2) <= UBOUND(cell_to_index, 2) .AND. &
                dc(3) >= LBOUND(cell_to_index, 3) .AND. dc(3) <= UBOUND(cell_to_index, 3)) THEN
               IF (cell_to_index(dc(1), dc(2), dc(3)) <= nimg) THEN
                  reldiff(b_img, y_img) = cell_to_index(dc(1), dc(2), dc(3))
               END IF
            END IF
         END DO
      END DO

      ! dense density images, mirroring the legacy storage convention exactly:
      ! the symmetric neighbor list stores each unordered atom pair {a, b} in
      ! ONE orientation chosen by the parity of a+b (qs_neighbor_lists), so the
      ! stored block (a, b) of image z holds M(a@0, b@s*c_z) with s = +1 for
      ! even pairs (a+b even, incl. same atom) and s = -1 for odd pairs. Each
      ! stored block deposits into (a, b) of the s*z image and (b, a)^T of its
      ! opposite (P_ab^z = P_ba^-z); the same-atom blocks of a mirror pair are
      ! stored twice and both deposits land on the image diagonal, hence 0.5.
      ALLOCATE (pdiff(nao, nao, nimg), ptmp(nao, nao))
      ! MPI: replicated density images (nao^2, small) so the parity-rule fold
      ! sees every block on every rank; the caller's matrices stay untouched
      ALLOCATE (rho_rep(nimg))
      DO b_img = 1, nimg
         ALLOCATE (rho_rep(b_img)%matrix)
         CALL dbcsr_create(rho_rep(b_img)%matrix, template=rho_ao(1, b_img)%matrix)
         CALL dbcsr_copy(rho_rep(b_img)%matrix, rho_ao(1, b_img)%matrix)
         CALL dbcsr_replicate_all(rho_rep(b_img)%matrix)
      END DO
      pdiff = 0.0_dp
      DO b_img = 1, nimg
         CALL dbcsr_iterator_start(iter, rho_rep(b_img)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, Aat, Bat, blk)
            offI = SUM(ri_data%bsizes_AO(1:Aat - 1)); nAOi = ri_data%bsizes_AO(Aat)
            offJ = SUM(ri_data%bsizes_AO(1:Bat - 1)); nAOj = ri_data%bsizes_AO(Bat)
            lpair = MOD(Aat + Bat, 2) == 0
            IF (lpair) THEN
               half = 1.0_dp
               IF (Aat == Bat .AND. mimg(b_img) /= 0) half = 0.5_dp
               pdiff(offI + 1:offI + nAOi, offJ + 1:offJ + nAOj, b_img) = &
                  pdiff(offI + 1:offI + nAOi, offJ + 1:offJ + nAOj, b_img) + half*blk(:, :)
               IF (mimg(b_img) > 0) THEN
                  pdiff(offJ + 1:offJ + nAOj, offI + 1:offI + nAOi, mimg(b_img)) = &
                     pdiff(offJ + 1:offJ + nAOj, offI + 1:offI + nAOi, mimg(b_img)) &
                     + half*TRANSPOSE(blk(:, :))
               END IF
            ELSE
               ! odd pair: the stored block lives at the negated cell, so the
               ! deposits swap image and transpose
               pdiff(offJ + 1:offJ + nAOj, offI + 1:offI + nAOi, b_img) = &
                  pdiff(offJ + 1:offJ + nAOj, offI + 1:offI + nAOi, b_img) + TRANSPOSE(blk(:, :))
               IF (mimg(b_img) > 0) THEN
                  pdiff(offI + 1:offI + nAOi, offJ + 1:offJ + nAOj, mimg(b_img)) = &
                     pdiff(offI + 1:offI + nAOi, offJ + 1:offJ + nAOj, mimg(b_img)) + blk(:, :)
               END IF
            END IF
         END DO
         CALL dbcsr_iterator_stop(iter)
      END DO
      ! S1 (symmetric closure): images without a mirror participate one-sidedly
      ! only; their density is voided so that no term references them
      DO b_img = 1, nimg
         IF (mimg(b_img) == 0) pdiff(:, :, b_img) = 0.0_dp
      END DO
      ! TODO(dev experiment, remove before PR): finite-difference gate -- perturb
      ! pdiff by fd_delta * X with a deterministic symmetric pattern (X^b and its
      ! mirror transpose-paired, self-mirror images symmetrized; mirrorless stay
      ! void, matching the S1 convention). E(P+/-dX) from two runs vs
      ! tr(F_sym.X) tests whether the assembled operator is dE/dP
      INQUIRE (FILE="PARI_FD_DELTA", EXIST=qexist)
      IF (qexist) THEN
         OPEN (NEWUNIT=ienv, FILE="PARI_FD_DELTA", status="old", action="read")
         READ (ienv, *) fd_delta
         CLOSE (ienv)
         DO b_img = 1, nimg
            IF (mimg(b_img) > b_img) THEN
               DO imu = 1, nao
                  DO inu = 1, nao
                     fd_x = 0.5_dp*(SIN(0.1_dp*imu + 0.37_dp*inu + 0.71_dp*b_img) &
                                    + SIN(0.1_dp*inu + 0.37_dp*imu + 0.71_dp*mimg(b_img)))
                     pdiff(imu, inu, b_img) = pdiff(imu, inu, b_img) + fd_delta*fd_x
                     pdiff(inu, imu, mimg(b_img)) = pdiff(inu, imu, mimg(b_img)) + fd_delta*fd_x
                  END DO
               END DO
            ELSE IF (mimg(b_img) == b_img) THEN
               DO imu = 1, nao
                  DO inu = 1, imu
                     fd_x = 0.5_dp*(SIN(0.1_dp*imu + 0.37_dp*inu + 0.71_dp*b_img) &
                                    + SIN(0.1_dp*inu + 0.37_dp*imu + 0.71_dp*b_img))
                     pdiff(imu, inu, b_img) = pdiff(imu, inu, b_img) + fd_delta*fd_x
                     IF (inu /= imu) pdiff(inu, imu, b_img) = pdiff(inu, imu, b_img) + fd_delta*fd_x
                  END DO
               END DO
            END IF
         END DO
      END IF
      DO b_img = 1, nimg
         CALL dbcsr_release(rho_rep(b_img)%matrix)
      END DO
      DEALLOCATE (rho_rep)
      ALLOCATE (pdiff_norm(nimg))
      DO b_img = 1, nimg
         pdiff_norm(b_img) = SQRT(SUM(pdiff(:, :, b_img)**2))
      END DO

      ! pool lookup per (X, Y, y), pooled C3 norms, and the flat live-pair
      ! lists (dead image pairs have exactly-zero C3; the Coulomb 2c tail keeps
      ! them in the pair enumeration, so presence alone is not a screen)
      ALLOCATE (pair3(natom, natom, nimg), c3norm(ri_data%pari_npair))
      ALLOCATE (live_cnt(natom, natom), live_off(natom, natom))
      pair3 = 0
      c3norm = 0.0_dp
      live_cnt = 0
      live_off = 1
      DO ip = 1, ri_data%pari_npair
         c3norm(ip) = SQRT(SUM(ri_data%pari_c3(ri_data%pari_c3_offset(ip): &
                                               ri_data%pari_c3_offset(ip + 1) - 1)**2))
         IF (ri_data%pari_pair_img(ip) >= 1 .AND. ri_data%pari_pair_img(ip) <= nimg) THEN
            pair3(ri_data%pari_pair_i(ip), ri_data%pari_pair_j(ip), ri_data%pari_pair_img(ip)) = ip
         END IF
      END DO
      ! S3 (symmetric closure): a pool pair participates only if its mirror pair
      ! (j, i, -img) is also present (one-sided absence is dropped on both sides);
      ! liveness uses the union of the mirror-pair norms (both sides carry their
      ! own solve's values). Mirrorless images cannot have mirror pairs.
      ALLOCATE (pair3s(natom, natom, nimg), c3norms(ri_data%pari_npair))
      pair3s = pair3
      DO ip = 1, ri_data%pari_npair
         c3norms(ip) = c3norm(ip)
         IF (ri_data%pari_pair_img(ip) >= 1 .AND. ri_data%pari_pair_img(ip) <= nimg) THEN
            Aat = ri_data%pari_pair_i(ip); Bat = ri_data%pari_pair_j(ip)
            y_img = ri_data%pari_pair_img(ip)
            z_img = mimg(y_img)
            IF (z_img == 0) THEN
               pair3s(Aat, Bat, y_img) = 0
            ELSE IF (pair3(Bat, Aat, z_img) == 0) THEN
               pair3s(Aat, Bat, y_img) = 0
            ELSE
               c3norms(ip) = MAX(c3norm(ip), c3norm(pair3(Bat, Aat, z_img)))
            END IF
         END IF
      END DO
      ! live-pair gate relative to the pool's largest C3 norm: scale-free
      ! across basis sizes (the measured pool is bimodal -- exactly dead vs
      ! live, no gray band -- so this selects the same set as any threshold
      ! many orders below the smallest live norm)
      cn_thr = 1.0E-12_dp*MAXVAL(c3norms)
      DO ip = 1, ri_data%pari_npair
         IF (ri_data%pari_pair_img(ip) >= 1 .AND. ri_data%pari_pair_img(ip) <= nimg) THEN
            IF (pair3s(ri_data%pari_pair_i(ip), ri_data%pari_pair_j(ip), ri_data%pari_pair_img(ip)) /= 0 .AND. &
                c3norms(ip) > cn_thr) live_cnt(ri_data%pari_pair_i(ip), ri_data%pari_pair_j(ip)) = &
               live_cnt(ri_data%pari_pair_i(ip), ri_data%pari_pair_j(ip)) + 1
         END IF
      END DO
      nlive = SUM(live_cnt)
      ALLOCATE (live_y(nlive), live_ip(nlive))
      live_off(1, 1) = 1
      DO Aat = 1, natom
         DO Bat = 1, natom
            IF (Aat == 1 .AND. Bat == 1) THEN
               CYCLE
            ELSE IF (Bat == 1) THEN
               live_off(Aat, Bat) = live_off(Aat - 1, natom) + live_cnt(Aat - 1, natom)
            ELSE
               live_off(Aat, Bat) = live_off(Aat, Bat - 1) + live_cnt(Aat, Bat - 1)
            END IF
         END DO
      END DO
      ALLOCATE (live_fill(natom, natom))
      live_fill = live_off
      DO ip = 1, ri_data%pari_npair
         IF (ri_data%pari_pair_img(ip) >= 1 .AND. ri_data%pari_pair_img(ip) <= nimg) THEN
            Aat = ri_data%pari_pair_i(ip); Bat = ri_data%pari_pair_j(ip)
            IF (c3norms(ip) > cn_thr .AND. &
                pair3s(Aat, Bat, ri_data%pari_pair_img(ip)) /= 0) THEN
               live_y(live_fill(Aat, Bat)) = ri_data%pari_pair_img(ip)
               live_ip(live_fill(Aat, Bat)) = ip
               live_fill(Aat, Bat) = live_fill(Aat, Bat) + 1
            END IF
         END IF
      END DO
      DEALLOCATE (live_fill)

      ALLOCATE (imgcell(3, nimg))
      imgcell(:, :) = index_to_cell(:, 1:nimg)

      ! density-live images
      ALLOCATE (pd_list(nimg))
      npd = 0
      DO b_img = 1, nimg
         IF (pdiff_norm(b_img) > 1.0E-14_dp) THEN
            npd = npd + 1
            pd_list(npd) = b_img
         END IF
      END DO

      ALLOCATE (fimg(nao, nao, nimg))
      fimg = 0.0_dp
      ! TODO(dev diagnostic, remove before PR): per-term F images and far-chain
      ! sources (T1 by |y| / |z1|, T4 by |e| / |v| / |w|; cell norm >= 3 is far)
      IF (pari_diag) THEN
         ALLOCATE (fimg_t(nao, nao, nimg, 4))
         fimg_t = 0.0_dp
         ALLOCATE (f1src_y(nao, nao, nimg), f1src_z(nao, nao, nimg))
         ALLOCATE (f4src_e(nao, nao, nimg), f4src_v(nao, nao, nimg), f4src_w(nao, nao, nimg))
         f1src_y = 0.0_dp; f1src_z = 0.0_dp
         f4src_e = 0.0_dp; f4src_v = 0.0_dp; f4src_w = 0.0_dp
      END IF
      nrimax = MAXVAL(ri_data%bsizes_RI)
      naomax = MAXVAL(ri_data%bsizes_AO)
      ALLOCATE (vblk(nrimax, nrimax))
      ALLOCATE (bd2_blk(nrimax, naomax, naomax))
      ALLOCATE (Gj(naomax*nrimax*naomax))
      ALLOCATE (y_list(nimg), ylive(nimg), pd_rev(nimg))
      ALLOCATE (bd2buf(nrimax, nao, naomax, MAX(1, npd)))
      pd_rev = 0
      DO iz = 1, npd
         pd_rev(pd_list(iz)) = iz
      END DO

      DO i = 1, natom
         ! MPI: atoms interleaved over ranks; fimg row blocks are disjoint
         ! (the four terms of atom i write only rows offI+1:offI+nAOi)
         IF (MOD(i - 1, para_env%num_pe) /= para_env%mepos) CYCLE
         nRIi = ri_data%bsizes_RI(i); nAOi = ri_data%bsizes_AO(i)
         offI = SUM(ri_data%bsizes_AO(1:i - 1))
         WRITE (progbuf, '(I0,A,I0)') i, "/", natom

         ! ---- BD for atom i over all (l, y): live x/z enumeration; threads
         ! own the density-image index z, which maps one-to-one to y = x + z
         ! for a fixed x (disjoint writes); a barrier per live x protects the
         ! accumulation across different x into the same y. One joint GEMM per
         ! (k, x, z): the pool segment is read in place as (nRIi*nAOi, nAOk)
         ! and the y-slab of bd_all as (nRIi*nAOi, nao)
         ALLOCATE (bd_all(nRIi, nAOi, nao, nimg))
         bd_all = 0.0_dp
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(i, ri_data, natom, npd, nimg, nao, bd_all, cell_to_index, imgcell, live_off, live_cnt, &
!$OMP        live_y, live_ip, pd_list, pdiff, nRIi, nAOi) &
!$OMP PRIVATE(k, ie, x_img, iz, z_img, y_img, nAOk, offK, ithr, nthr)
         ithr = 0
         nthr = 1
!$       ithr = omp_get_thread_num()
!$       nthr = omp_get_max_threads()
         DO k = 1, natom
            nAOk = ri_data%bsizes_AO(k)
            offK = SUM(ri_data%bsizes_AO(1:k - 1))
            DO ie = live_off(i, k), live_off(i, k) + live_cnt(i, k) - 1
               x_img = live_y(ie)
               DO iz = 1, npd
                  IF (MOD(iz - 1, nthr) /= ithr) CYCLE
                  z_img = pd_list(iz)
                  y_img = pari_imgsum(cell_to_index, imgcell, nimg, imgcell(:, x_img), imgcell(:, z_img))
                  IF (y_img == 0) CYCLE
                  CALL dgemm('N', 'N', nRIi*nAOi, nao, nAOk, 1.0_dp, &
                             ri_data%pari_c3(ri_data%pari_c3_offset(live_ip(ie))), nRIi*nAOi, &
                             pdiff(offK + 1, 1, z_img), nao, 1.0_dp, bd_all(1, 1, 1, y_img), nRIi*nAOi)
               END DO
!$OMP BARRIER
            END DO
         END DO
!$OMP END PARALLEL
         nyl = 0
         ylive = .FALSE.
         DO b_img = 1, nimg
            IF (mimg(b_img) > 0 .AND. SQRT(SUM(bd_all(:, :, :, b_img)**2)) > 0.0_dp) THEN
               nyl = nyl + 1
               y_list(nyl) = b_img
               ylive(b_img) = .TRUE.
            END IF
         END DO
         CALL pari_stage_banner(ri_data, "KS assembly: BD done, four-term stages atom "//TRIM(progbuf))
         t_stage = m_walltime()

         ! ================= T1: hoist G(l,y) = BD_i(l,y)^T . V^y[i,l] once per
         ! (l, y) into the shared Gflat table (joint (q, lam) layout, block per
         ! (l, y), (l, iy) ownership, disjoint writes); the (j, b) blocks are
         ! then scattered under (j, b) ownership with ONE joint GEMM per combo
         ! (the pool segment is read in place as (nRIl*nAOl, nAOj))
         ALLOCATE (glen(natom))
         DO l = 1, natom
            glen(l) = nAOi*ri_data%bsizes_RI(l)*ri_data%bsizes_AO(l)
         END DO
         gstep = SUM(glen)
         ALLOCATE (Gflat(nyl*gstep))
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(pari_term_mask, i, ri_data, natom, nimg, nao, nyl, y_list, bd_all, fimg, reldiff, pair3s, c3norms, &
!$OMP        cn_thr, nRIi, nAOi, &
!$OMP        offI, glen, gstep, Gflat, nrimax, fimg_t, f1src_y, f1src_z, imgcell, pari_hard_gate, &
!$OMP        pari_diag) &
!$OMP PRIVATE(l, nRIl, nAOl, offL, iy, y_img, lam, lc, jb, j, nAOj, offJ, b_img, z1, vblk, ithr, nthr)
         ithr = 0
         nthr = 1
!$       ithr = omp_get_thread_num()
!$       nthr = omp_get_max_threads()
         IF (ALLOCATED(vblk)) DEALLOCATE (vblk)
         ALLOCATE (vblk(nrimax, nrimax))
         DO lc = 1, natom*nyl
            IF (MOD(lc - 1, nthr) /= ithr) CYCLE
            l = (lc - 1)/nyl + 1
            iy = lc - (l - 1)*nyl
            y_img = y_list(iy)
            nRIl = ri_data%bsizes_RI(l); nAOl = ri_data%bsizes_AO(l)
            offL = SUM(ri_data%bsizes_AO(1:l - 1))
            CALL pari_vblock(vblk, ri_data%pari_mat_2c(y_img), i, l, nRIi, nRIl)
            DO lam = 1, nAOl
               CALL dgemm('T', 'N', nAOi, nRIl, nRIi, 1.0_dp, bd_all(1, 1, offL + lam, y_img), nRIi, &
                          vblk, nrimax, 0.0_dp, &
                          Gflat((iy - 1)*gstep + SUM(glen(1:l - 1)) + (lam - 1)*nAOi*nRIl + 1), nAOi)
            END DO
         END DO
!$OMP BARRIER
         DO jb = 1, MERGE(natom*nimg, 0, IAND(pari_term_mask, 1) /= 0)
            IF (MOD(jb - 1, nthr) /= ithr) CYCLE
            j = (jb - 1)/nimg + 1
            b_img = jb - (j - 1)*nimg
            nAOj = ri_data%bsizes_AO(j)
            offJ = SUM(ri_data%bsizes_AO(1:j - 1))
            DO l = 1, natom
               nRIl = ri_data%bsizes_RI(l); nAOl = ri_data%bsizes_AO(l)
               DO iy = 1, nyl
                  y_img = y_list(iy)
                  z1 = reldiff(b_img, y_img)
                  IF (z1 <= 0) CYCLE
                  IF (pair3s(l, j, z1) == 0) CYCLE
                  IF (c3norms(pair3s(l, j, z1)) <= cn_thr) CYCLE
                  CALL dgemm('N', 'N', nAOi, nAOj, nRIl*nAOl, 1.0_dp, &
                             Gflat((iy - 1)*gstep + SUM(glen(1:l - 1)) + 1), nAOi, &
                             ri_data%pari_c3(ri_data%pari_c3_offset(pair3s(l, j, z1))), nRIl*nAOl, &
                             1.0_dp, fimg(offI + 1, offJ + 1, b_img), nao)
                  ! TODO(dev diagnostic, remove before PR): per-term + far-chain
                  IF (pari_diag) THEN
                     CALL dgemm('N', 'N', nAOi, nAOj, nRIl*nAOl, 1.0_dp, &
                                Gflat((iy - 1)*gstep + SUM(glen(1:l - 1)) + 1), nAOi, &
                                ri_data%pari_c3(ri_data%pari_c3_offset(pair3s(l, j, z1))), nRIl*nAOl, &
                                1.0_dp, fimg_t(offI + 1, offJ + 1, b_img, 1), nao)
                     IF (SUM(ABS(imgcell(:, y_img))) >= pari_hard_gate) THEN
                        CALL dgemm('N', 'N', nAOi, nAOj, nRIl*nAOl, 1.0_dp, &
                                   Gflat((iy - 1)*gstep + SUM(glen(1:l - 1)) + 1), nAOi, &
                                   ri_data%pari_c3(ri_data%pari_c3_offset(pair3s(l, j, z1))), nRIl*nAOl, &
                                   1.0_dp, f1src_y(offI + 1, offJ + 1, b_img), nao)
                     END IF
                     IF (SUM(ABS(imgcell(:, z1))) >= pari_hard_gate) THEN
                        CALL dgemm('N', 'N', nAOi, nAOj, nRIl*nAOl, 1.0_dp, &
                                   Gflat((iy - 1)*gstep + SUM(glen(1:l - 1)) + 1), nAOi, &
                                   ri_data%pari_c3(ri_data%pari_c3_offset(pair3s(l, j, z1))), nRIl*nAOl, &
                                   1.0_dp, f1src_z(offI + 1, offJ + 1, b_img), nao)
                     END IF
                  END IF
               END DO
            END DO
         END DO
         DEALLOCATE (vblk)
!$OMP END PARALLEL
         DEALLOCATE (Gflat, glen)
         WRITE (progbuf, '(A,I0,A,I0,A,F8.1,A)') "KS assembly: T1 done atom ", i, "/", natom, &
            " (wall ", m_walltime() - t_stage, " s)"
         CALL pari_stage_banner(ri_data, TRIM(progbuf))

         ! ================= T3: stream BD2(k,x,z,l) over the live (k, x, z)
         ! with per-mu GEMMs off the pool, contract with V^z[k,l] and fold
         ! directly into the shared Hflat table (joint (q, lam) layout, block
         ! per (k, l), image u = x+z outermost, beta=1 accumulation); the
         ! (j, b) scatter is ONE joint GEMM per (k, u, l, j, b). The
         ! y-in-y_list restriction of the original walk is preserved through
         ! the ylive gate on u.
         ALLOCATE (slg(natom), hoff(natom, natom))
         DO l = 1, natom
            slg(l) = nAOi*ri_data%bsizes_RI(l)*ri_data%bsizes_AO(l)
         END DO
         htot = nimg*SUM(slg)
         DO k = 1, natom
            DO l = 1, natom
               hoff(k, l) = (k - 1)*htot + nimg*SUM(slg(1:l - 1))
            END DO
         END DO
         ALLOCATE (Hflat(natom*htot))
         ALLOCATE (liveH(natom, nimg, natom))
         Hflat = 0.0_dp
         liveH = .FALSE.
         ki_tot = live_off(natom, i) + live_cnt(natom, i) - live_off(1, i)
         t_stage = m_walltime()
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(i, ri_data, natom, nimg, nao, npd, nrimax, naomax, cell_to_index, imgcell, live_off, live_cnt, &
!$OMP        live_y, live_ip, pd_list, pdiff, reldiff, mimg, ylive, pair3s, Hflat, hoff, slg, liveH, &
!$OMP        nAOi, ki_tot) &
!$OMP PRIVATE(k, nRIk, nAOk, offK, ie, x_img, iz, z_img, u_img, l, nRIl, nAOl, offL, lam, mu, bd2_blk, &
!$OMP        vblk, ithr, nthr, npct, lpct, progbuf)
         ithr = 0
         nthr = 1
         lpct = 0
!$       ithr = omp_get_thread_num()
!$       nthr = omp_get_max_threads()
         IF (ALLOCATED(bd2_blk)) DEALLOCATE (bd2_blk)
         IF (ALLOCATED(vblk)) DEALLOCATE (vblk)
         ALLOCATE (bd2_blk(nrimax, naomax, naomax))
         ALLOCATE (vblk(nrimax, nrimax))
         DO k = 1, natom
            nRIk = ri_data%bsizes_RI(k); nAOk = ri_data%bsizes_AO(k)
            offK = SUM(ri_data%bsizes_AO(1:k - 1))
            DO ie = live_off(k, i), live_off(k, i) + live_cnt(k, i) - 1
               IF (mimg(live_y(ie)) == 0) CYCLE
               x_img = mimg(live_y(ie))
               DO iz = 1, npd
                  IF (MOD(iz - 1, nthr) /= ithr) CYCLE
                  z_img = pd_list(iz)
                  u_img = pari_imgsum(cell_to_index, imgcell, nimg, imgcell(:, x_img), imgcell(:, z_img))
                  IF (u_img == 0) CYCLE
                  IF (.NOT. ylive(u_img)) CYCLE
                  DO l = 1, natom
                     nRIl = ri_data%bsizes_RI(l); nAOl = ri_data%bsizes_AO(l)
                     offL = SUM(ri_data%bsizes_AO(1:l - 1))
                     IF (pair3s(k, l, z_img) == 0) CYCLE
                     CALL pari_vblock(vblk, ri_data%pari_mat_2c(z_img), k, l, nRIk, nRIl)
                     ! BD2(p, lam, mu) = sum_s C3(k,i,-x)(p,s,mu) . P(s^z, lam^l),
                     ! one GEMM per mu writing the (p, lam) face in place
                     DO mu = 1, nAOi
                        CALL dgemm('N', 'N', nRIk, nAOl, nAOk, 1.0_dp, &
                                   ri_data%pari_c3(ri_data%pari_c3_offset(live_ip(ie)) + &
                                                   (mu - 1)*nRIk*nAOk), nRIk, &
                                   pdiff(offK + 1, offL + 1, z_img), nao, 0.0_dp, &
                                   bd2_blk(1, 1, mu), nrimax)
                     END DO
                     DO lam = 1, nAOl
                        CALL dgemm('T', 'N', nAOi, nRIl, nRIk, 1.0_dp, bd2_blk(1, lam, 1), nrimax*naomax, &
                                   vblk, nrimax, 1.0_dp, &
                                   Hflat(hoff(k, l) + (u_img - 1)*slg(l) + (lam - 1)*nAOi*nRIl + 1), nAOi)
                     END DO
                     liveH(k, u_img, l) = .TRUE.
                  END DO
               END DO
!$OMP BARRIER
               IF (ithr == 0) THEN
                  npct = 10*(ie - live_off(1, i) + 1)/MAX(1, ki_tot)
                  IF (npct > lpct) THEN
                     lpct = npct
                     WRITE (progbuf, '(A,I0,A,I0,A,I0,A,I0)') "KS assembly: T3 build atom ", i, "/", natom, ": ", &
                        ie - live_off(1, i) + 1, "/", ki_tot
                     CALL pari_stage_banner(ri_data, TRIM(progbuf))
                  END IF
               END IF
            END DO
         END DO
         DEALLOCATE (bd2_blk, vblk)
!$OMP END PARALLEL
         WRITE (progbuf, '(A,I0,A,I0,A,F8.1,A)') "KS assembly: T3 build done atom ", i, "/", natom, &
            " (wall ", m_walltime() - t_stage, " s)"
         CALL pari_stage_banner(ri_data, TRIM(progbuf))
         t_stage = m_walltime()
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(pari_term_mask, i, ri_data, natom, nimg, nao, fimg, reldiff, pair3s, c3norms, cn_thr, Hflat, hoff, slg, liveH, &
!$OMP        nAOi, offI, fimg_t, pari_diag) &
!$OMP PRIVATE(jb, j, nAOj, offJ, b_img, k, l, nRIl, nAOl, iu, z1, ithr, nthr)
         ithr = 0
         nthr = 1
!$       ithr = omp_get_thread_num()
!$       nthr = omp_get_max_threads()
         DO jb = 1, MERGE(natom*nimg, 0, IAND(pari_term_mask, 4) /= 0)
            IF (MOD(jb - 1, nthr) /= ithr) CYCLE
            j = (jb - 1)/nimg + 1
            b_img = jb - (j - 1)*nimg
            nAOj = ri_data%bsizes_AO(j)
            offJ = SUM(ri_data%bsizes_AO(1:j - 1))
            DO k = 1, natom
               DO l = 1, natom
                  nRIl = ri_data%bsizes_RI(l); nAOl = ri_data%bsizes_AO(l)
                  DO iu = 1, nimg
                     IF (.NOT. liveH(k, iu, l)) CYCLE
                     z1 = reldiff(b_img, iu)
                     IF (z1 <= 0) CYCLE
                     IF (pair3s(l, j, z1) == 0) CYCLE
                     IF (c3norms(pair3s(l, j, z1)) <= cn_thr) CYCLE
                     CALL dgemm('N', 'N', nAOi, nAOj, nRIl*nAOl, 1.0_dp, &
                                Hflat(hoff(k, l) + (iu - 1)*slg(l) + 1), nAOi, &
                                ri_data%pari_c3(ri_data%pari_c3_offset(pair3s(l, j, z1))), nRIl*nAOl, &
                                1.0_dp, fimg(offI + 1, offJ + 1, b_img), nao)
                     ! TODO(dev diagnostic, remove before PR): per-term F
                     IF (pari_diag) THEN
                        CALL dgemm('N', 'N', nAOi, nAOj, nRIl*nAOl, 1.0_dp, &
                                   Hflat(hoff(k, l) + (iu - 1)*slg(l) + 1), nAOi, &
                                   ri_data%pari_c3(ri_data%pari_c3_offset(pair3s(l, j, z1))), nRIl*nAOl, &
                                   1.0_dp, fimg_t(offI + 1, offJ + 1, b_img, 3), nao)
                     END IF
                  END DO
               END DO
            END DO
         END DO
!$OMP END PARALLEL
         DEALLOCATE (Hflat, liveH, slg, hoff)
         WRITE (progbuf, '(A,I0,A,I0,A,F8.1,A)') "KS assembly: T3 scatter done atom ", i, "/", natom, &
            " (wall ", m_walltime() - t_stage, " s)"
         CALL pari_stage_banner(ri_data, TRIM(progbuf))

         ! ================= T4: per (k, x) the bd2 buffer over (z, l) is
         ! filled once (per-mu GEMMs off the pool) and read as a shifted view:
         ! V is indexed by v = b - x and the right fit by w = y - b, so
         ! z = w + v with the implicit y = w + b required to be a valid
         ! y_list image (the original y walk); the per-(w) contraction folds
         ki_tot = live_off(natom, i) + live_cnt(natom, i) - live_off(1, i)
         t_stage = m_walltime()
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(pari_term_mask, i, ri_data, natom, nimg, nao, npd, nrimax, naomax, cell_to_index, imgcell, live_off, live_cnt, &
!$OMP        pari_vgate, pari_hard_gate, live_y, live_ip, pd_list, pd_rev, pdiff, reldiff, mimg, ylive, pair3s, &
!$OMP        c3norms, cn_thr, bd2buf, &
!$OMP        fimg, nAOi, offI, ki_tot, fimg_t, f4src_e, f4src_v, f4src_w, pari_diag) &
!$OMP PRIVATE(k, nRIk, nAOk, offK, ie, x_img, iz, z_img, l, nAOl, offL, mu, jb, j, nRIj, nAOj, offJ, &
!$OMP        b_img, v_img, w_img, y_img, lam, Gj, vblk, ithr, nthr, npct, lpct, progbuf, iw, nwl, nzl, &
!$OMP        izl, rr, w_iz, w_im, w_off, izl_of, Ap, Wt)
         ithr = 0
         nthr = 1
         lpct = 0
!$       ithr = omp_get_thread_num()
!$       nthr = omp_get_max_threads()
         IF (ALLOCATED(Gj)) DEALLOCATE (Gj)
         IF (ALLOCATED(vblk)) DEALLOCATE (vblk)
         ALLOCATE (Gj(naomax*nrimax*naomax))
         ALLOCATE (vblk(nrimax, nrimax))
         ALLOCATE (w_iz(nimg), w_im(nimg), w_off(nimg), izl_of(MAX(1, npd)))
         DO k = 1, natom
            nRIk = ri_data%bsizes_RI(k); nAOk = ri_data%bsizes_AO(k)
            offK = SUM(ri_data%bsizes_AO(1:k - 1))
            DO ie = live_off(k, i), live_off(k, i) + live_cnt(k, i) - 1
               IF (mimg(live_y(ie)) == 0) CYCLE
               x_img = mimg(live_y(ie))
               DO iz = 1, npd
                  IF (MOD(iz - 1, nthr) /= ithr) CYCLE
                  z_img = pd_list(iz)
                  ! one nao-wide GEMM per mu: the C3 mu-slab against the whole
                  ! density row block of atom k (every partner column at once);
                  ! the consumers read their lam range at offL + lam
                  DO mu = 1, nAOi
                     CALL dgemm('N', 'N', nRIk, nao, nAOk, 1.0_dp, &
                                ri_data%pari_c3(ri_data%pari_c3_offset(live_ip(ie)) + &
                                                (mu - 1)*nRIk*nAOk), nRIk, &
                                pdiff(offK + 1, 1, z_img), nao, 0.0_dp, &
                                bd2buf(1, 1, mu, iz), nrimax)
                  END DO
               END DO
!$OMP BARRIER
               DO jb = 1, MERGE(natom*nimg, 0, IAND(pari_term_mask, 8) /= 0)
                  IF (MOD(jb - 1, nthr) /= ithr) CYCLE
                  j = (jb - 1)/nimg + 1
                  b_img = jb - (j - 1)*nimg
                  nRIj = ri_data%bsizes_RI(j); nAOj = ri_data%bsizes_AO(j)
                  offJ = SUM(ri_data%bsizes_AO(1:j - 1))
                  v_img = reldiff(b_img, x_img)
                  IF (v_img <= 0) CYCLE
                  ! TODO(dev experiment): vacuum-V guard, mirrors the pool prune
                  IF (SUM(ABS(imgcell(:, v_img))) >= pari_vgate) CYCLE
                  IF (pair3s(k, j, v_img) == 0) CYCLE
                  CALL pari_vblock(vblk, ri_data%pari_mat_2c(v_img), k, j, nRIk, nRIj)
                  IF (SQRT(SUM(vblk(:nRIk, :nRIj)**2)) == 0.0_dp) CYCLE
                  DO l = 1, natom
                     nAOl = ri_data%bsizes_AO(l)
                     offL = SUM(ri_data%bsizes_AO(1:l - 1))
                     ! collect the live (w -> z) set with the gates of the
                     ! original walk, in the original w order
                     nwl = 0
                     nzl = 0
                     DO iz = 1, npd
                        izl_of(iz) = 0
                     END DO
                     ! hoist (oracle-certified, slice2b/b12_hoist.py): the two
                     ! pool gates below are exactly the live-table predicate
                     ! (same pool order, ascending w), so walk the (j,l) live-w
                     ! list instead of scanning all nimg per combo
                     DO lwi = live_off(j, l), live_off(j, l) + live_cnt(j, l) - 1
                        w_img = live_y(lwi)
                        ! y = w + b must be a valid image AND a y_list image
                        y_img = pari_imgsum(cell_to_index, imgcell, nimg, imgcell(:, w_img), imgcell(:, b_img))
                        IF (y_img == 0) CYCLE
                        IF (.NOT. ylive(y_img)) CYCLE
                        ! density image z = w + v must be density-live
                        z_img = pari_imgsum(cell_to_index, imgcell, nimg, imgcell(:, w_img), imgcell(:, v_img))
                        IF (z_img == 0) CYCLE
                        iz = pd_rev(z_img)
                        IF (iz == 0) CYCLE
                        nwl = nwl + 1
                        w_iz(nwl) = iz
                        w_im(nwl) = w_img
                        w_off(nwl) = ri_data%pari_c3_offset(pair3s(j, l, w_img))
                        IF (izl_of(iz) == 0) THEN
                           nzl = nzl + 1
                           izl_of(iz) = nzl
                        END IF
                     END DO
                     IF (nwl == 0) CYCLE
                     ALLOCATE (Ap(nrimax, nAOi*nzl), Wt(nRIj, nAOi*nzl))
                     ! b7c iz-stack one lam at a time (certified, slice2b/
                     ! b9_joint.py): ONE stacked GEMM per lam covers every
                     ! needed z image, and each (w, lam) scatter reads the
                     ! own=j pool segment's (p, mu) plane of THIS lam directly
                     ! (contiguous, lda = nRIj) -- no transposed pool copy and
                     ! no joint G buffer
                     DO lam = 1, nAOl
                        DO iz = 1, npd
                           izl = izl_of(iz)
                           IF (izl == 0) CYCLE
                           DO mu = 1, nAOi
                              Ap(1:nRIk, (izl - 1)*nAOi + mu) = bd2buf(1:nRIk, offL + lam, mu, iz)
                           END DO
                        END DO
                        CALL dgemm('T', 'N', nRIj, nzl*nAOi, nRIk, 1.0_dp, vblk, nrimax, &
                                   Ap, nrimax, 0.0_dp, Wt, nRIj)
                        DO iw = 1, nwl
                           izl = izl_of(w_iz(iw))
                           DO rr = 1, nRIj
                              Gj((rr - 1)*nAOi + 1:rr*nAOi) = Wt(rr, (izl - 1)*nAOi + 1:izl*nAOi)
                           END DO
                           CALL dgemm('N', 'N', nAOi, nAOj, nRIj, 1.0_dp, Gj, nAOi, &
                                      ri_data%pari_c3(w_off(iw) + (lam - 1)*nRIj*nAOj), nRIj, &
                                      1.0_dp, fimg(offI + 1, offJ + 1, b_img), nao)
                           ! TODO(dev diagnostic, remove before PR): per-term + far-chain
                           IF (pari_diag) THEN
                              CALL dgemm('N', 'N', nAOi, nAOj, nRIj, 1.0_dp, Gj, nAOi, &
                                         ri_data%pari_c3(w_off(iw) + (lam - 1)*nRIj*nAOj), nRIj, &
                                         1.0_dp, fimg_t(offI + 1, offJ + 1, b_img, 4), nao)
                              IF (SUM(ABS(imgcell(:, live_y(ie)))) >= pari_hard_gate) THEN
                                 CALL dgemm('N', 'N', nAOi, nAOj, nRIj, 1.0_dp, Gj, nAOi, &
                                            ri_data%pari_c3(w_off(iw) + (lam - 1)*nRIj*nAOj), nRIj, &
                                            1.0_dp, f4src_e(offI + 1, offJ + 1, b_img), nao)
                              END IF
                              IF (SUM(ABS(imgcell(:, v_img))) >= pari_hard_gate) THEN
                                 CALL dgemm('N', 'N', nAOi, nAOj, nRIj, 1.0_dp, Gj, nAOi, &
                                            ri_data%pari_c3(w_off(iw) + (lam - 1)*nRIj*nAOj), nRIj, &
                                            1.0_dp, f4src_v(offI + 1, offJ + 1, b_img), nao)
                              END IF
                              IF (SUM(ABS(imgcell(:, w_im(iw)))) >= pari_hard_gate) THEN
                                 CALL dgemm('N', 'N', nAOi, nAOj, nRIj, 1.0_dp, Gj, nAOi, &
                                            ri_data%pari_c3(w_off(iw) + (lam - 1)*nRIj*nAOj), nRIj, &
                                            1.0_dp, f4src_w(offI + 1, offJ + 1, b_img), nao)
                              END IF
                           END IF
                        END DO
                     END DO
                     DEALLOCATE (Ap, Wt)
                  END DO
               END DO
!$OMP BARRIER
               IF (ithr == 0) THEN
                  npct = 10*(ie - live_off(1, i) + 1)/MAX(1, ki_tot)
                  IF (npct > lpct) THEN
                     lpct = npct
                     WRITE (progbuf, '(A,I0,A,I0,A,I0,A,I0)') "KS assembly: T4 atom ", i, "/", natom, ": ", &
                        ie - live_off(1, i) + 1, "/", ki_tot
                     CALL pari_stage_banner(ri_data, TRIM(progbuf))
                  END IF
               END IF
            END DO
         END DO
         DEALLOCATE (Gj, vblk, w_iz, w_im, w_off, izl_of)
!$OMP END PARALLEL
         WRITE (progbuf, '(A,I0,A,I0,A,F8.1,A)') "KS assembly: T4 done atom ", i, "/", natom, &
            " (wall ", m_walltime() - t_stage, " s)"
         CALL pari_stage_banner(ri_data, TRIM(progbuf))

         ! ================= T2: BD_i . V^b[i,j] . C3(j, l, y-b); V depends on
         ! the output image b, so no hoist is available (per-(i,j,b) groups);
         ! per lam ONE stacked GEMM covers the live y set, then per (y, lam)
         ! scatters read the pool planes directly (b9_joint.py)
         t_stage = m_walltime()
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(pari_term_mask, i, ri_data, natom, nimg, nao, nrimax, naomax, nyl, y_list, bd_all, fimg, reldiff, pair3s, &
!$OMP        c3norms, cn_thr, nRIi, nAOi, offI, fimg_t, pari_diag) &
!$OMP PRIVATE(jb, j, nAOj, nRIj, offJ, b_img, l, nAOl, offL, iy, y_img, z2, lam, vblk, Gj, ithr, nthr, &
!$OMP        nyl_l, yl_y, yl_off, Ap2, Wt2, mu, rr)
         ithr = 0
         nthr = 1
!$       ithr = omp_get_thread_num()
!$       nthr = omp_get_max_threads()
         IF (ALLOCATED(vblk)) DEALLOCATE (vblk)
         IF (ALLOCATED(Gj)) DEALLOCATE (Gj)
         ALLOCATE (vblk(nrimax, nrimax))
         ALLOCATE (Gj(naomax*nrimax*naomax))
         ALLOCATE (yl_y(nyl), yl_off(nyl))
         DO jb = 1, MERGE(natom*nimg, 0, IAND(pari_term_mask, 2) /= 0)
            IF (MOD(jb - 1, nthr) /= ithr) CYCLE
            j = (jb - 1)/nimg + 1
            b_img = jb - (j - 1)*nimg
            nAOj = ri_data%bsizes_AO(j); nRIj = ri_data%bsizes_RI(j)
            offJ = SUM(ri_data%bsizes_AO(1:j - 1))
            IF (pair3s(i, j, b_img) == 0) CYCLE
            CALL pari_vblock(vblk, ri_data%pari_mat_2c(b_img), i, j, nRIi, nRIj)
            IF (SQRT(SUM(vblk(:nRIi, :nRIj)**2)) == 0.0_dp) CYCLE
            DO l = 1, natom
               nAOl = ri_data%bsizes_AO(l)
               offL = SUM(ri_data%bsizes_AO(1:l - 1))
               ! collect the live (y -> z2) set in the original walk order
               nyl_l = 0
               DO iy = 1, nyl
                  y_img = y_list(iy)
                  z2 = reldiff(y_img, b_img)
                  IF (z2 <= 0) CYCLE
                  IF (pair3s(j, l, z2) == 0) CYCLE
                  IF (c3norms(pair3s(j, l, z2)) <= cn_thr) CYCLE
                  nyl_l = nyl_l + 1
                  yl_y(nyl_l) = y_img
                  yl_off(nyl_l) = ri_data%pari_c3_offset(pair3s(j, l, z2))
               END DO
               IF (nyl_l == 0) CYCLE
               ALLOCATE (Ap2(nrimax, nAOi*nyl_l), Wt2(nRIj, nAOi*nyl_l))
               ! b7c y-stack one lam at a time (certified, slice2b/b9_joint.py):
               ! ONE stacked GEMM per lam covers every live y, and each (y, lam)
               ! scatter reads the own=j pool segment's (p, mu) plane of THIS
               ! lam directly (contiguous, lda = nRIj) -- no transposed copy
               DO lam = 1, nAOl
                  DO iy = 1, nyl_l
                     DO mu = 1, nAOi
                        Ap2(1:nRIi, (iy - 1)*nAOi + mu) = bd_all(1:nRIi, mu, offL + lam, yl_y(iy))
                     END DO
                  END DO
                  CALL dgemm('T', 'N', nRIj, nyl_l*nAOi, nRIi, 1.0_dp, vblk, nrimax, &
                             Ap2, nrimax, 0.0_dp, Wt2, nRIj)
                  DO iy = 1, nyl_l
                     DO rr = 1, nRIj
                        Gj((rr - 1)*nAOi + 1:rr*nAOi) = Wt2(rr, (iy - 1)*nAOi + 1:iy*nAOi)
                     END DO
                     CALL dgemm('N', 'N', nAOi, nAOj, nRIj, 1.0_dp, Gj, nAOi, &
                                ri_data%pari_c3(yl_off(iy) + (lam - 1)*nRIj*nAOj), nRIj, &
                                1.0_dp, fimg(offI + 1, offJ + 1, b_img), nao)
                     ! TODO(dev diagnostic, remove before PR): per-term F
                     IF (pari_diag) THEN
                        CALL dgemm('N', 'N', nAOi, nAOj, nRIj, 1.0_dp, Gj, nAOi, &
                                   ri_data%pari_c3(yl_off(iy) + (lam - 1)*nRIj*nAOj), nRIj, &
                                   1.0_dp, fimg_t(offI + 1, offJ + 1, b_img, 2), nao)
                     END IF
                  END DO
               END DO
               DEALLOCATE (Ap2, Wt2)
            END DO
         END DO
         DEALLOCATE (vblk, Gj, yl_y, yl_off)
!$OMP END PARALLEL
         WRITE (progbuf, '(A,I0,A,I0,A,F8.1,A)') "KS assembly: T2 done atom ", i, "/", natom, &
            " (wall ", m_walltime() - t_stage, " s)"
         CALL pari_stage_banner(ri_data, TRIM(progbuf))

         DEALLOCATE (bd_all)
      END DO
      DEALLOCATE (vblk, bd2_blk, Gj, y_list, ylive, pd_list, pd_rev, bd2buf)

      ! MPI: complete the atom-split partial images on every rank (dense
      ! sum over disjoint row blocks); everything downstream (symmetrization,
      ! energy, KS puts, dumps) then runs on identical complete data
      CALL para_env%sum(fimg)

      ! structural identity F^b == (F^{-b})^T, measured on the RAW assembly
      ! BEFORE the symmetric projection below (which would mask any violation);
      ! this is the primary probe of the image algebra, printed every SCF step
      worst_sym = 0.0_dp
      DO b_img = 1, nimg
         IF (mimg(b_img) > 0) THEN
            worst_sym = MAX(worst_sym, &
                            SQRT(SUM((fimg(:, :, b_img) - TRANSPOSE(fimg(:, :, mimg(b_img))))**2))/ &
                            MAX(SQRT(SUM(fimg(:, :, b_img)**2)), 1.0E-30_dp))
         END IF
      END DO
      WRITE (*, '(A,ES10.3)') " PARI| image symmetry F^b=(F^-b)^T, worst rel deviation: ", worst_sym
      IF (worst_sym > 1.0E-10_dp) THEN
         CALL cp_warn(__LOCATION__, "PARI KS assembly: image symmetry F^b = (F^-b)^T violated "// &
                      "beyond 1e-10 (check the image algebra)")
      END IF

      ! symmetrize the image pairs: the exact exchange satisfies F^b = (F^-b)^T,
      ! so any deviation is error (from one-sided screening or from noise of the
      ! per-pair truncated solves amplified by pair-domain near-null spaces);
      ! projecting onto the constraint manifold halves it and restores a
      ! Hermitian k-space Fock. Output images whose mirror lies outside the HFX
      ! image range are voided: keeping them one-sided is energy-invisible (S1
      ! already voids the matching density images), but breaks the Hermiticity
      ! of F(k) = sum_b F^b exp(ikR_b) and poisons the eigensolver
      DO b_img = 1, nimg
         IF (mimg(b_img) > b_img) THEN
            ptmp(:, :) = 0.5_dp*(fimg(:, :, b_img) + TRANSPOSE(fimg(:, :, mimg(b_img))))
            fimg(:, :, mimg(b_img)) = 0.5_dp*(fimg(:, :, mimg(b_img)) + TRANSPOSE(fimg(:, :, b_img)))
            fimg(:, :, b_img) = ptmp(:, :)
         ELSE IF (mimg(b_img) == b_img) THEN
            fimg(:, :, b_img) = 0.5_dp*(fimg(:, :, b_img) + TRANSPOSE(fimg(:, :, b_img)))
         ELSE IF (mimg(b_img) == 0) THEN
            fimg(:, :, b_img) = 0.0_dp
         END IF
      END DO

      !TODO(dev dump, remove before PR): per-image F, per-term F + symmetry, full P, V blocks
      IF (para_env%mepos == para_env%num_pe - 1) THEN
         OPEN (NEWUNIT=unit_dump, FILE="PARI_K_DUMP")
         WRITE (unit_dump, *) nimg, nao, worst_sym
         WRITE (unit_dump, *) imgcell
         DO b_img = 1, nimg
            WRITE (unit_dump, *) b_img
            WRITE (unit_dump, '(ES25.16)') pdiff(:, :, b_img)
         END DO
         DO b_img = 1, nimg
            WRITE (unit_dump, *) b_img
            WRITE (unit_dump, '(ES25.16)') fimg(:, :, b_img)
         END DO
         ! TODO(dev diagnostic, remove before PR): per-term sections only with
         ! PARI_KEEP_DIAG (the pdiff/fimg sections above stay -- regression channel)
         IF (pari_diag) THEN
            DO ip = 1, 4
               worst_sym = 0.0_dp
               DO b_img = 1, nimg
                  IF (mimg(b_img) > 0) THEN
                     worst_sym = MAX(worst_sym, &
                                     SQRT(SUM((fimg_t(:, :, b_img, ip) - TRANSPOSE(fimg_t(:, :, mimg(b_img), ip)))**2))/ &
                                     MAX(SQRT(SUM(fimg_t(:, :, b_img, ip)**2)), 1.0E-30_dp))
                  END IF
               END DO
               WRITE (unit_dump, *) ip
               WRITE (unit_dump, '(ES25.16)') worst_sym
               DO b_img = 1, nimg
                  WRITE (unit_dump, *) b_img
                  WRITE (unit_dump, '(ES25.16)') fimg_t(:, :, b_img, ip)
               END DO
            END DO
         END IF
         CLOSE (unit_dump)
      END IF

      ! energy, legacy convention: 0.5 sum_b Tr((-fac F^b) P^b)
      ehfx = 0.0_dp
      DO b_img = 1, nimg
         trb = SUM(fimg(:, :, b_img)*pdiff(:, :, b_img))
         ehfx = ehfx - 0.5_dp*fac*trb
      END DO

      ! accumulate -fac F^b into the symmetric-typed KS image matrices, upper
      ! blocks only, following the storage parity convention (same rule the
      ! legacy fold implements via checker_tr): even atom pairs store the
      ! block of image z at z; odd pairs store it at the NEGATED cell, i.e.
      ! the (a, b) block of image z carries F(a@0, b@-c_z), which by
      ! F^b = (F^-b)^T equals the (a, b) block of the opposite image
      DO b_img = 1, nimg
         DO Aat = 1, natom
            DO Bat = Aat, natom
               offI = SUM(ri_data%bsizes_AO(1:Aat - 1)); nAOi = ri_data%bsizes_AO(Aat)
               offJ = SUM(ri_data%bsizes_AO(1:Bat - 1)); nAOj = ri_data%bsizes_AO(Bat)
               src = b_img
               IF (MOD(Aat + Bat, 2) /= 0) src = mimg(b_img)
               IF (src == 0) CYCLE
               CALL dbcsr_get_stored_coordinates(ks_matrix(1, b_img)%matrix, Aat, Bat, owner)
               IF (owner == para_env%mepos) THEN
                  ! summation=.TRUE. so the contract does not depend on the
                  ! caller zeroing ks_matrix first (it does today)
                  CALL dbcsr_put_block(ks_matrix(1, b_img)%matrix, Aat, Bat, &
                                       -fac*fimg(offI + 1:offI + nAOi, offJ + 1:offJ + nAOj, src), summation=.TRUE.)
               END IF
            END DO
         END DO
      END DO

      ! TODO(dev diagnostic, remove before PR): far-output (cell norm >= 3)
      ! content split by chain source: T1 via |y| (BD/V side) or |z1| (right
      ! factor), T4 via |e| (BD2), |v| (V) or |w| (right factor)
      IF (pari_diag) THEN
         d1y = 0.0_dp; d1z = 0.0_dp; d4e = 0.0_dp; d4v = 0.0_dp; d4w = 0.0_dp
         DO b_img = 1, nimg
            IF (SUM(ABS(imgcell(:, b_img))) < 3) CYCLE
            d1y = d1y + SUM(f1src_y(:, :, b_img)**2)
            d1z = d1z + SUM(f1src_z(:, :, b_img)**2)
            d4e = d4e + SUM(f4src_e(:, :, b_img)**2)
            d4v = d4v + SUM(f4src_v(:, :, b_img)**2)
            d4w = d4w + SUM(f4src_w(:, :, b_img)**2)
         END DO
         WRITE (*, '(A,5ES13.4)') " PARI-DIAG| far-img source sq-norms T1(y),T1(z1),T4(e),T4(v),T4(w): ", &
            SQRT(d1y), SQRT(d1z), SQRT(d4e), SQRT(d4v), SQRT(d4w)
         DEALLOCATE (f1src_y, f1src_z, f4src_e, f4src_v, f4src_w)
      END IF

      IF (ALLOCATED(fimg_t)) DEALLOCATE (fimg_t)
      DEALLOCATE (fimg, pdiff, pdiff_norm, pair3, pair3s, mimg, reldiff, ptmp, live_cnt, live_off, &
                  live_y, live_ip, imgcell, c3norm, c3norms)

      CALL timestop(handle)

   END SUBROUTINE pari_assemble_ks

! **************************************************************************************************
!> \brief stage banner for the PARI k-point twin: goes to unit_nr when active,
!>        otherwise to stdout only when the KP_RI_PROGRESS_BAR print flag is on
!> \param ri_data ...
!> \param msg ...
! **************************************************************************************************
   SUBROUTINE pari_stage_banner(ri_data, msg)
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      CHARACTER(LEN=*), INTENT(IN)                       :: msg

      ! unit_nr comes from the HF_INFO print key (default level MEDIUM):
      ! silent at LOW, prints from MEDIUM upward, like the legacy KP banners
      IF (ri_data%unit_nr > 0) THEN
         WRITE (ri_data%unit_nr, '(T3,A,A)') "PARI| ", msg
         CALL m_flush(ri_data%unit_nr)
      END IF
   END SUBROUTINE pari_stage_banner

! **************************************************************************************************
!> \brief index of the image at cell ca+cb, 0 when outside the k-point cell
!>        list or beyond the HFX image range
!> \param cell_to_index ...
!> \param imgcell ...
!> \param nimg ...
!> \param ca ...
!> \param cb ...
!> \return ...
! **************************************************************************************************
   FUNCTION pari_imgsum(cell_to_index, imgcell, nimg, ca, cb) RESULT(img)
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: imgcell
      INTEGER, INTENT(IN)                                :: nimg
      INTEGER, DIMENSION(3), INTENT(IN)                  :: ca, cb
      INTEGER                                            :: img

      INTEGER, DIMENSION(3)                              :: cc

      cc = ca + cb
      IF (cc(1) >= LBOUND(cell_to_index, 1) .AND. cc(1) <= UBOUND(cell_to_index, 1) .AND. &
          cc(2) >= LBOUND(cell_to_index, 2) .AND. cc(2) <= UBOUND(cell_to_index, 2) .AND. &
          cc(3) >= LBOUND(cell_to_index, 3) .AND. cc(3) <= UBOUND(cell_to_index, 3)) THEN
         img = cell_to_index(cc(1), cc(2), cc(3))
         IF (img > nimg .OR. img < 1) img = 0
      ELSE
         img = 0
      END IF
   END FUNCTION pari_imgsum

! **************************************************************************************************
!> \brief copy the (rowat, colat) block of a per-image 2c matrix into vblk,
!>        zeroing it when the block is absent (screened image pair)
!> \param vblk ...
!> \param mat2c ...
!> \param rowat ...
!> \param colat ...
!> \param nr ...
!> \param nc ...
! **************************************************************************************************
   SUBROUTINE pari_vblock(vblk, mat2c, rowat, colat, nr, nc)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: vblk
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat2c
      INTEGER, INTENT(IN)                                :: rowat, colat, nr, nc

      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk

      NULLIFY (blk)
      CALL dbcsr_get_block_p(matrix=mat2c, row=rowat, col=colat, block=blk, found=found)
      IF (found) THEN
         vblk(1:nr, 1:nc) = blk(1:nr, 1:nc)
      ELSE
         vblk(1:nr, 1:nc) = 0.0_dp
      END IF
   END SUBROUTINE pari_vblock

! **************************************************************************************************
!> \brief PARI pre-SCF: build the per-image 2c metric, stream the 3c integrals
!>        through per-pair domain solves, and cache the pooled fit coefficients.
!>        Twin of hfx_ri_pre_scf_kp (hfx_ri_kp.F): the 3c tensor carries the
!>        extended per-image RI blocks of get_kp_and_ri_images, but PARI only
!>        ever reads the home-cell block (index katom), and the per-atom global
!>        (bumped) metric inverse is replaced by the per-pair truncated
!>        pseudo-inverse.
!> \param dbcsr_template ...
!> \param ri_data ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE hfx_ri_pre_scf_kp_pari(dbcsr_template, ri_data, qs_env)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: dbcsr_template
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ri_pre_scf_kp_pari'

      INTEGER                                            :: handle, i_img, i_RI, iatom, ienv, jatom, &
                                                            natom, ncell_RI_pari, nimg, nkind, &
                                                            npair, unit_dump
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bsizes_RI_ext, dist_1, dist_2, dist_3, &
                                                            edist_1, edist_2, edist_3, &
                                                            img_to_RI_cell_pari, pair_i, pair_img, &
                                                            pair_j
      INTEGER, DIMENSION(3)                              :: pcoord, pdims
      INTEGER, DIMENSION(:), POINTER                     :: col_bsize, row_bsize
      LOGICAL                                            :: qexist, qexist_ri
      REAL(KIND=dp)                                      :: ri_range_prv
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: vblk
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk_2c
      TYPE(dbcsr_distribution_type)                      :: dbcsr_dist
      TYPE(dbcsr_iterator_type)                          :: iter_2c
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: t_2c_pari
      TYPE(dbt_distribution_type)                        :: t_dist
      TYPE(dbt_pgrid_type)                               :: pgrid
      TYPE(dbt_type)                                     :: t_3c_dummy, t_3c_dummy_ext
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :)       :: t_3c_pari
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_2d_type), POINTER                :: dist_2d
      TYPE(distribution_3d_type)                         :: dist_3d
      TYPE(gto_basis_set_p_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: basis_set_AO, basis_set_RI
      TYPE(mp_cart_type)                                 :: mp_comm_t3c
      TYPE(mp_para_env_type)                             :: para_env_1
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_3c_type)                        :: nl_3c
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: nl_2c
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

! TODO(dev dump, remove before PR): RI_range runtime knob
! TODO(dev dump, remove before PR): V dump unit + view buffer
! MPI: size-1 sub-environment for the redundant per-rank 3c build

      NULLIFY (dft_control, dist_2d, para_env, particle_set, qs_kind_set, nl_2c, blk_2c, &
               col_bsize, row_bsize)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dft_control=dft_control, natom=natom, nkind=nkind, &
                      para_env=para_env, distribution_2d=dist_2d, particle_set=particle_set, &
                      qs_kind_set=qs_kind_set)

      IF (ri_data%flavor /= ri_pmat) CPABORT("K-points RI-HFX only with RHO flavor")
      IF (.NOT. ri_data%same_op) THEN
         CALL cp_abort(__LOCATION__, "RI%KP_METHOD PARI requires RI_METRIC HFX (metric "// &
                       "identical to the interaction potential) in this first version")
      END IF
      IF (ABS(ri_data%eps_pgf_orb - dft_control%qs_control%eps_pgf_orb) > 1.0E-16_dp) THEN
         CPABORT("RI%EPS_PGF_ORB and QS%EPS_PGF_ORB must be identical for RI-HFX k-points")
      END IF
      ! MPI: the 3c tensors are built redundantly per rank on a size-1
      ! sub-communicator (after the RI_range cut they hold home-aux blocks
      ! only, i.e. pool-scale); the pair solves below are rank-split
      CALL para_env_1%from_split(para_env, para_env%mepos)

      CALL pari_cleanup(ri_data)

      CALL get_kp_and_ri_images(ri_data, qs_env)
      nimg = ri_data%nimg

      ALLOCATE (basis_set_RI(nkind), basis_set_AO(nkind))
      CALL basis_set_list_setup(basis_set_RI, ri_data%ri_basis_type, qs_kind_set)
      CALL basis_set_list_setup(basis_set_AO, ri_data%orb_basis_type, qs_kind_set)

      ! ---------------------------------------------------------------- 2c: (P_i^0 | Q_j^b) per image
      CALL build_2c_neighbor_lists(nl_2c, basis_set_RI, basis_set_RI, ri_data%hfx_pot, &
                                   "PARI_2c_nl", qs_env, sym_ij=.FALSE., dist_2d=dist_2d)

      CALL cp_dbcsr_dist2d_to_dist(dist_2d, dbcsr_dist)
      ALLOCATE (row_bsize(natom), col_bsize(natom))
      row_bsize(:) = ri_data%bsizes_RI
      col_bsize(:) = ri_data%bsizes_RI

      ALLOCATE (t_2c_pari(nimg))
      CALL dbcsr_create(t_2c_pari(1), "(P|Q) PARI", dbcsr_dist, dbcsr_type_no_symmetry, &
                        row_bsize, col_bsize)
      DO i_img = 2, nimg
         CALL dbcsr_create(t_2c_pari(i_img), template=t_2c_pari(1))
      END DO
      DEALLOCATE (row_bsize, col_bsize)
      CALL dbcsr_distribution_release(dbcsr_dist)

      CALL pari_stage_banner(ri_data, "building per-image 2c metric integrals")
      CALL build_2c_integrals(t_2c_pari, ri_data%filter_eps_2c, qs_env, nl_2c, basis_set_RI, &
                              basis_set_RI, ri_data%hfx_pot, do_kpoints=.TRUE., do_hfx_kpoints=.TRUE.)
      CALL release_neighbor_list_sets(nl_2c)

      ! replicated so that every rank can harvest per-pair blocks (small natom^2 objects)
      DO i_img = 1, nimg
         CALL dbcsr_replicate_all(t_2c_pari(i_img))
      END DO

      ! ---------------------------------------------------------------- 3c: (mu^0 nu^b | P^0) per image
      ! aux-cell placement: get_kp_and_ri_images has already flagged the images
      ! hosting RI functions within kp_RI_range (ri_data%img_to_RI_cell, an
      ! extended basis of ri_data%ncell_RI cells). Folding every image into the
      ! main cell (the legacy hard-wire below) makes build_3c_integrals alias
      ! all aux images into the main-cell block and SUM them there
      ! (dbt_put_block summation=.TRUE.), so the tensor held
      ! sum_c (mu^0 nu^b | P_k^c) instead of the home-only (P_k^0) -- a
      ! cell-size dependent contamination of every pair solve. With the
      ! extended mapping the home block index stays (1-1)*natom + katom =
      ! katom, so all PARI reads are unchanged and clean
      ! TODO(dev experiment, remove before PR): PARI_RI_MODE file, 0 = legacy
      INQUIRE (FILE="PARI_RI_MODE", EXIST=qexist)
      IF (qexist) THEN
         OPEN (NEWUNIT=ienv, FILE="PARI_RI_MODE", status="old", action="read")
         READ (ienv, *) pari_ri_mode
         CLOSE (ienv)
      END IF
      ALLOCATE (img_to_RI_cell_pari(nimg))
      IF (pari_ri_mode == 0) THEN
         img_to_RI_cell_pari(:) = 1
         ncell_RI_pari = 1
      ELSE
         ! element-wise copy: ri_data%img_to_RI_cell is allocated with
         ! dft_control%nimages entries, which exceeds ri_data%nimg -- a whole-
         ! array assignment memcpy's the source length and overruns.
         ! Unflagged images (aux beyond kp_RI_range per the shared detection)
         ! get a private dead block: the 3c builder's looser, kind-radius based
         ! screening can still reach them, and a zero mapping there would make
         ! blk_idx(3) nonpositive; their out-of-range content lands in the dead
         ! block, is filtered away, and is never read (PARI reads home blocks)
         ncell_RI_pari = ri_data%ncell_RI
         DO i_img = 1, nimg
            img_to_RI_cell_pari(i_img) = ri_data%img_to_RI_cell(i_img)
            IF (img_to_RI_cell_pari(i_img) == 0) THEN
               ncell_RI_pari = ncell_RI_pari + 1
               img_to_RI_cell_pari(i_img) = ncell_RI_pari
            END IF
         END DO
      END IF

      pdims = 0
      ! tensor_dims should reflect the real block structure of the 3c tensor
      ! (GLOBAL practice, hfx_ri_kp.F); legacy keeps the historical value for
      ! bitwise reproduction
      IF (pari_ri_mode == 0) THEN
         CALL dbt_pgrid_create(para_env_1, pdims, pgrid, tensor_dims=[natom, natom, 1])
      ELSE
         CALL dbt_pgrid_create(para_env_1, pdims, pgrid, tensor_dims=[natom, natom, natom*ncell_RI_pari])
      END IF

      CALL create_3c_tensor(t_3c_dummy, dist_1, dist_2, dist_3, pgrid, &
                            ri_data%bsizes_AO, ri_data%bsizes_AO, ri_data%bsizes_RI, &
                            map1=[1, 2], map2=[3], name="(AO AO | RI) PARI")

      CALL dbt_mp_environ_pgrid(pgrid, pdims, pcoord)
      CALL mp_comm_t3c%create(pgrid%mp_comm_2d, 3, pdims)
      CALL distribution_3d_create(dist_3d, dist_1, dist_2, dist_3, nkind, particle_set, &
                                  mp_comm_t3c, own_comm=.TRUE.)

      CALL build_3c_neighbor_lists(nl_3c, basis_set_AO, basis_set_AO, basis_set_RI, dist_3d, &
                                   ri_data%ri_metric, "PARI_3c_nl", qs_env, op_pos=2, sym_ij=.FALSE., &
                                   own_dist=.TRUE.)

      ALLOCATE (t_3c_pari(nimg, 1))
      IF (pari_ri_mode == 0) THEN
         CALL dbt_distribution_new(t_dist, pgrid, dist_1, dist_2, dist_3)
         CALL dbt_create(t_3c_pari(1, 1), "PARI_3c_ints", t_dist, [1, 2], [3], &
                         ri_data%bsizes_AO, ri_data%bsizes_AO, ri_data%bsizes_RI)
      ELSE
         ! extended per-image RI blocks: the tensor's third dimension spans
         ! ncell_RI_pari cells of aux functions; neighbor-aux contributions
         ! land in their own blocks and are filtered away unread
         ALLOCATE (bsizes_RI_ext(natom*ncell_RI_pari))
         DO i_RI = 1, ncell_RI_pari
            bsizes_RI_ext((i_RI - 1)*natom + 1:i_RI*natom) = ri_data%bsizes_RI(:)
         END DO
         CALL create_3c_tensor(t_3c_dummy_ext, edist_1, edist_2, edist_3, pgrid, &
                               ri_data%bsizes_AO, ri_data%bsizes_AO, bsizes_RI_ext, &
                               map1=[1, 2], map2=[3], name="(AO AO | RI ext) PARI")
         CALL dbt_destroy(t_3c_dummy_ext)
         CALL dbt_distribution_new(t_dist, pgrid, edist_1, edist_2, edist_3)
         CALL dbt_create(t_3c_pari(1, 1), "PARI_3c_ints", t_dist, [1, 2], [3], &
                         ri_data%bsizes_AO, ri_data%bsizes_AO, bsizes_RI_ext)
         DEALLOCATE (edist_1, edist_2, edist_3, bsizes_RI_ext)
      END IF
      DO i_img = 2, nimg
         CALL dbt_create(t_3c_pari(1, 1), t_3c_pari(i_img, 1))
      END DO
      CALL dbt_distribution_destroy(t_dist)
      DEALLOCATE (dist_1, dist_2, dist_3)

      CALL pari_stage_banner(ri_data, "building per-image 3c integrals (the slow stage)")
      ! TODO(dev experiment, remove before PR): RI_range runtime knob
      ! (default eps = home-only; file restores the legacy kp_RI_range).
      ! The home-aux tuples carry positive roundoff in rik (~4 ulp, the
      ! iterator composes it), so the cut must be a small positive epsilon,
      ! not exact zero: the builder's gate is a strict dik > RI_range
      ri_range_prv = 1.0E-6_dp
      INQUIRE (FILE="PARI_RI_RANGE", EXIST=qexist_ri)
      IF (qexist_ri) THEN
         OPEN (NEWUNIT=ienv, FILE="PARI_RI_RANGE", status="old", action="read")
         READ (ienv, *) ri_range_prv
         CLOSE (ienv)
      END IF
      CALL build_3c_integrals(t_3c_pari, ri_data%filter_eps, qs_env, nl_3c, basis_set_AO, &
                              basis_set_AO, basis_set_RI, ri_data%ri_metric, &
                              int_eps=ri_data%eps_schwarz, op_pos=2, desymmetrize=.FALSE., &
                              do_kpoints=.TRUE., do_hfx_kpoints=.TRUE., &
                              ! every 3c read site contracts the home aux block only
                              ! (katom == first AO atom, dik within roundoff of zero),
                              ! so a tiny RI_range keeps every consumed block while the
                              ! builder skips far-aux images (~40x fewer 3c blocks in 2D/3D);
                              ! the knob above restores the legacy range for A/B
                              RI_range=ri_range_prv, img_to_RI_cell=img_to_RI_cell_pari)
      DO i_img = 1, nimg
         CALL dbt_filter(t_3c_pari(i_img, 1), ri_data%filter_eps)
      END DO

      CALL neighbor_list_3c_destroy(nl_3c)
      CALL dbt_destroy(t_3c_dummy)
      CALL dbt_pgrid_destroy(pgrid)
      DEALLOCATE (img_to_RI_cell_pari, basis_set_AO, basis_set_RI)

      ! ----------------------------------------------------- pairs and pooled C3 solve
      ! enumerate ordered pairs (i, j, b) from the 2c blocks: exactly the pairs whose
      ! pair-domain metric couples within range and above the 2c filter (norm
      ! screening: a pair may only be dropped when its whole cross metric block is
      ! negligible; distance cut-offs are not safe on the K axis)
      npair = 0
      DO i_img = 1, nimg
         CALL dbcsr_iterator_start(iter_2c, t_2c_pari(i_img))
         DO WHILE (dbcsr_iterator_blocks_left(iter_2c))
            CALL dbcsr_iterator_next_block(iter_2c, iatom, jatom, blk_2c)
            IF (SQRT(SUM(blk_2c**2)) <= 0.0_dp) CYCLE
            npair = npair + 1
         END DO
         CALL dbcsr_iterator_stop(iter_2c)
      END DO

      ALLOCATE (pair_i(npair), pair_j(npair), pair_img(npair))
      npair = 0
      DO i_img = 1, nimg
         CALL dbcsr_iterator_start(iter_2c, t_2c_pari(i_img))
         DO WHILE (dbcsr_iterator_blocks_left(iter_2c))
            CALL dbcsr_iterator_next_block(iter_2c, iatom, jatom, blk_2c)
            IF (SQRT(SUM(blk_2c**2)) <= 0.0_dp) CYCLE
            npair = npair + 1
            pair_i(npair) = iatom
            pair_j(npair) = jatom
            pair_img(npair) = i_img
         END DO
         CALL dbcsr_iterator_stop(iter_2c)
      END DO

      ! stream the 3c integrals through the per-pair solves and discard them
      CALL pari_stage_banner(ri_data, "solving pair domains")
      CALL pari_solve_pairs(ri_data, qs_env, t_2c_pari, t_3c_pari, pair_i, pair_j, pair_img)

      DO i_img = 1, nimg
         CALL dbt_destroy(t_3c_pari(i_img, 1))
      END DO
      DEALLOCATE (t_3c_pari)

      ALLOCATE (ri_data%pari_mat_2c(nimg))
      DO i_img = 1, nimg
         CALL dbcsr_create(ri_data%pari_mat_2c(i_img), template=t_2c_pari(i_img))
         CALL dbcsr_copy(ri_data%pari_mat_2c(i_img), t_2c_pari(i_img))
         CALL dbcsr_release(t_2c_pari(i_img))
      END DO
      DEALLOCATE (t_2c_pari)

      CALL MOVE_ALLOC(pair_i, ri_data%pari_pair_i)
      CALL MOVE_ALLOC(pair_j, ri_data%pari_pair_j)
      CALL MOVE_ALLOC(pair_img, ri_data%pari_pair_img)
      ri_data%pari_npair = npair
      ri_data%pari_nimg = nimg

      IF (ri_data%unit_nr > 0) THEN
         WRITE (ri_data%unit_nr, FMT="(T3,A,I0)") &
            "PARI| Number of ordered atom pairs with metric contact: ", npair
         WRITE (ri_data%unit_nr, FMT="(T3,A,ES12.3)") &
            "PARI| Worst pair-domain reciprocal condition number: ", ri_data%pari_rcond_max
         CALL m_flush(ri_data%unit_nr)
      END IF

      !TODO(dev dump, remove before PR): all nonzero V blocks (pari_vblock
      ! view), written once here -- the pair metric is fixed after pre-SCF
      ALLOCATE (vblk(MAXVAL(ri_data%bsizes_RI), MAXVAL(ri_data%bsizes_RI)))
      OPEN (NEWUNIT=unit_dump, FILE="PARI_V_DUMP")
      DO i_img = 1, nimg
         DO iatom = 1, natom
            DO jatom = 1, natom
               CALL pari_vblock(vblk, ri_data%pari_mat_2c(i_img), iatom, jatom, &
                                ri_data%bsizes_RI(iatom), ri_data%bsizes_RI(jatom))
               IF (SQRT(SUM(vblk(:ri_data%bsizes_RI(iatom), :ri_data%bsizes_RI(jatom))**2)) > 0.0_dp) THEN
                  WRITE (unit_dump, *) i_img, iatom, jatom
                  WRITE (unit_dump, '(ES25.16)') vblk(:ri_data%bsizes_RI(iatom), :ri_data%bsizes_RI(jatom))
               END IF
            END DO
         END DO
      END DO
      CLOSE (unit_dump)
      DEALLOCATE (vblk)
      CALL para_env_1%free()

      CALL timestop(handle)

   END SUBROUTINE hfx_ri_pre_scf_kp_pari

! **************************************************************************************************
!> \brief release the PARI pre-SCF cache and the image bookkeeping that
!>        get_kp_and_ri_images re-allocates unconditionally
!> \param ri_data ...
! **************************************************************************************************
   SUBROUTINE pari_cleanup(ri_data)
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data

      INTEGER                                            :: i

      IF (ALLOCATED(ri_data%pari_mat_2c)) THEN
         DO i = 1, SIZE(ri_data%pari_mat_2c)
            CALL dbcsr_release(ri_data%pari_mat_2c(i))
         END DO
         DEALLOCATE (ri_data%pari_mat_2c)
      END IF
      IF (ALLOCATED(ri_data%pari_c3)) DEALLOCATE (ri_data%pari_c3)
      IF (ALLOCATED(ri_data%pari_c3_offset)) DEALLOCATE (ri_data%pari_c3_offset)
      IF (ALLOCATED(ri_data%pari_pair_i)) DEALLOCATE (ri_data%pari_pair_i)
      IF (ALLOCATED(ri_data%pari_pair_j)) DEALLOCATE (ri_data%pari_pair_j)
      IF (ALLOCATED(ri_data%pari_pair_img)) DEALLOCATE (ri_data%pari_pair_img)
      IF (ALLOCATED(ri_data%pari_nrow_i)) DEALLOCATE (ri_data%pari_nrow_i)
      IF (ALLOCATED(ri_data%img_to_RI_cell)) DEALLOCATE (ri_data%img_to_RI_cell)
      IF (ALLOCATED(ri_data%RI_cell_to_img)) DEALLOCATE (ri_data%RI_cell_to_img)
      IF (ALLOCATED(ri_data%present_images)) DEALLOCATE (ri_data%present_images)
      ri_data%pari_npair = 0
      ri_data%pari_nimg = 0
      ri_data%pari_rcond_max = 0.0_dp

   END SUBROUTINE pari_cleanup

! **************************************************************************************************
!> \brief per-pair domain assembly and truncated pseudo-inverse solve; the pooled
!>        C3 coefficients and the dump consumed by the pyscf oracle gate are written
!>        into ri_data (TODO: remove the dev dump before the PR)
!> \param ri_data ...
!> \param qs_env ...
!> \param t_2c_pari per-image (P_i^0|Q_j^b)
!> \param t_3c_pari per-image (mu^0 nu^b | P^0)
!> \param pair_i ordered pair list: own atom
!> \param pair_j ordered pair list: partner atom
!> \param pair_img ordered pair list: image of the partner
! **************************************************************************************************
   SUBROUTINE pari_solve_pairs(ri_data, qs_env, t_2c_pari, t_3c_pari, pair_i, pair_j, pair_img)
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_type), DIMENSION(:), INTENT(INOUT)      :: t_2c_pari
      TYPE(dbt_type), DIMENSION(:, :), INTENT(INOUT)     :: t_3c_pari
      INTEGER, DIMENSION(:), INTENT(IN)                  :: pair_i, pair_j, pair_img

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'pari_solve_pairs'

      CHARACTER(LEN=64)                                  :: progbuf, progmsg
      INTEGER :: b_img, cell(3), handle, i, iatom, ienv, ipair, ithr, jatom, l, mb_img, nAOi, &
         nAOj, ncol, ndom, ni, nj, npair, nthr, offset, unit_dump, unit_pool
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: img_to_mimg
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: img_cells
      INTEGER, DIMENSION(3)                              :: ind
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: found, qexist, same_site
      REAL(KIND=dp)                                      :: fit_res, rcond, res_max
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: a2, c_full, vdom, vinv
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: blk_3c
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: vii, vij, vjj
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env

      NULLIFY (kpoints, cell_to_index, vij, vii, vjj, para_env)

      CALL timeset(routineN, handle)

      npair = SIZE(pair_i)
      CALL get_qs_env(qs_env, kpoints=kpoints, para_env=para_env)
      CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index)

      ! image index of the opposite cell -b for every image b (the k-cell box is symmetric)
      ALLOCATE (img_to_mimg(ri_data%nimg))
      img_to_mimg = 0
      ALLOCATE (img_cells(3, ri_data%nimg))
      img_cells = 0
      DO b_img = 1, ri_data%nimg
         found = .FALSE.
         cell = 0
         DO i = LBOUND(cell_to_index, 1), UBOUND(cell_to_index, 1)
            DO jatom = LBOUND(cell_to_index, 2), UBOUND(cell_to_index, 2)
               DO l = LBOUND(cell_to_index, 3), UBOUND(cell_to_index, 3)
                  IF (cell_to_index(i, jatom, l) == b_img) THEN
                     cell = [i, jatom, l]
                     found = .TRUE.
                  END IF
               END DO
            END DO
         END DO
         ! the main cell is [0, 0, 0], so emptiness must not be probed via
         ! ANY(cell /= 0) -- that silently zeroed the partner rows of every
         ! main-image pair (img_to_mimg(1) == 0) and blew up their C3 solves
         IF (found) THEN
            i = cell_to_index(-cell(1), -cell(2), -cell(3))
            IF (i >= 1 .AND. i <= ri_data%nimg) img_to_mimg(b_img) = i
         END IF
         ! TODO(dev experiment, remove before PR): store the cell vector for
         ! the distance prune of physically empty far pairs
         IF (found) img_cells(:, b_img) = cell
      END DO

      ALLOCATE (ri_data%pari_c3_offset(npair + 1))
      ri_data%pari_c3_offset(1) = 1
      DO ipair = 1, npair
         ni = ri_data%bsizes_RI(pair_i(ipair))
         ri_data%pari_c3_offset(ipair + 1) = ri_data%pari_c3_offset(ipair) + &
                                             ni*ri_data%bsizes_AO(pair_i(ipair))*ri_data%bsizes_AO(pair_j(ipair))
      END DO
      ALLOCATE (ri_data%pari_c3(ri_data%pari_c3_offset(npair + 1) - 1))
      ri_data%pari_c3 = 0.0_dp ! rank-split fills disjoint segments; sum below
      ALLOCATE (ri_data%pari_nrow_i(npair))
      ri_data%pari_nrow_i(:) = ri_data%bsizes_RI(pair_i(:))

      ! TODO(dev experiment, remove before PR): runtime remote-cut knob
      INQUIRE (FILE="PARI_REMOTE_CUT", EXIST=qexist)
      IF (qexist) THEN
         OPEN (NEWUNIT=ienv, FILE="PARI_REMOTE_CUT", status="old", action="read")
         READ (ienv, *) pari_remote_cut
         CLOSE (ienv)
      END IF
      ! TODO(dev experiment, remove before PR): eig truncation knob
      INQUIRE (FILE="PARI_EPS_EIG", EXIST=qexist)
      IF (qexist) THEN
         OPEN (NEWUNIT=ienv, FILE="PARI_EPS_EIG", status="old", action="read")
         READ (ienv, *) pari_eps_eig
         CLOSE (ienv)
      END IF
      ! TODO(dev experiment, remove before PR): same-atom cross-image pair scaling
      INQUIRE (FILE="PARI_SAC_F", EXIST=qexist)
      IF (qexist) THEN
         OPEN (NEWUNIT=ienv, FILE="PARI_SAC_F", status="old", action="read")
         READ (ienv, *) pari_sac_f
         CLOSE (ienv)
      END IF
      !TODO(dev dump, remove before PR): C3 pool dump for the pyscf oracle gate
      IF (para_env%num_pe == 1) OPEN (NEWUNIT=unit_dump, FILE="PARI_C3_DUMP")
      IF (para_env%num_pe == 1) OPEN (NEWUNIT=unit_pool, FILE="PARI_POOL_DUMP")
      progmsg = ""
      res_max = 0.0_dp
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(ri_data, t_2c_pari, t_3c_pari, pair_i, pair_j, pair_img, npair, img_to_mimg, img_cells, para_env, &
!$OMP        unit_dump, unit_pool, progmsg, pari_remote_cut, pari_sac_f, pari_eps_eig, res_max) &
!$OMP PRIVATE(ipair, iatom, jatom, b_img, mb_img, same_site, ni, nj, nAOi, nAOj, ncol, ndom, offset, &
!$OMP         ind, found, rcond, vdom, a2, c_full, vinv, blk_3c, vii, vij, vjj, progbuf, ithr, nthr, fit_res)
      ithr = 0
      nthr = 1
!$    ithr = omp_get_thread_num()
!$    nthr = omp_get_max_threads()
      DO ipair = 1, npair
         IF (MOD(ipair - 1, nthr*para_env%num_pe) /= para_env%mepos*nthr + ithr) CYCLE
         IF (MOD(ipair, MAX(1, npair/20)) == 0) THEN
!$OMP CRITICAL(pari_banner)
            WRITE (progbuf, '(I0,A,I0)') ipair, " / ", npair
            progmsg = "pair solves: "//TRIM(progbuf)
            CALL pari_stage_banner(ri_data, progmsg)
!$OMP END CRITICAL(pari_banner)
         END IF
         iatom = pair_i(ipair)
         jatom = pair_j(ipair)
         b_img = pair_img(ipair)
         same_site = (iatom == jatom .AND. b_img == 1)
         ni = ri_data%bsizes_RI(iatom)
         nj = ri_data%bsizes_RI(jatom)
         nAOi = ri_data%bsizes_AO(iatom)
         nAOj = ri_data%bsizes_AO(jatom)
         ncol = nAOi*nAOj
         offset = ri_data%pari_c3_offset(ipair) - 1

         IF (same_site) THEN
            ! single-atom domain
            ndom = ni
            ALLOCATE (vdom(ndom, ndom), a2(ndom, ncol), c_full(ndom, ncol), vinv(ndom, ndom))
            NULLIFY (vii)
            CALL dbcsr_get_block_p(matrix=t_2c_pari(1), row=iatom, col=iatom, block=vii, found=found)
            CPASSERT(found)
            vdom(:, :) = vii(:, :)
            ! (P_i^0 | mu_i^0 nu_i^0): a2(p, o + (q-1)*nAOi) = blk(o, q, p)
            ind = [iatom, iatom, iatom]
            CALL get_3c_or_zero(t_3c_pari(b_img, 1), ind, [nAOi, nAOi, ni], blk_3c)
            a2(:, :) = TRANSPOSE(RESHAPE(blk_3c, [nAOi*nAOi, ni]))
            DEALLOCATE (blk_3c)
         ELSE
            ! domain aux(i^0) u aux(j^b)
            ndom = ni + nj
            ALLOCATE (vdom(ndom, ndom), a2(ndom, ncol), c_full(ndom, ncol), vinv(ndom, ndom))
            NULLIFY (vii, vjj, vij)
            CALL dbcsr_get_block_p(matrix=t_2c_pari(1), row=iatom, col=iatom, block=vii, found=found)
            CPASSERT(found)
            CALL dbcsr_get_block_p(matrix=t_2c_pari(1), row=jatom, col=jatom, block=vjj, found=found)
            CPASSERT(found)
            CALL dbcsr_get_block_p(matrix=t_2c_pari(b_img), row=iatom, col=jatom, block=vij, found=found)
            CPASSERT(found)
            vdom(1:ni, 1:ni) = vii(:, :)
            vdom(ni + 1:ndom, ni + 1:ndom) = vjj(:, :)
            vdom(1:ni, ni + 1:ndom) = vij(:, :)
            vdom(ni + 1:ndom, 1:ni) = TRANSPOSE(vij(:, :))
            ! own-atom rows: a2(1:ni, p; o + (q-1)*nAOi) = (P_i^0 | mu_i^0 nu_j^b) = blk(o, q, p)
            ind = [iatom, jatom, iatom]
            CALL get_3c_or_zero(t_3c_pari(b_img, 1), ind, [nAOi, nAOj, ni], blk_3c)
            a2(1:ni, :) = TRANSPOSE(RESHAPE(blk_3c, [nAOi*nAOj, ni]))
            DEALLOCATE (blk_3c)
            ! partner rows: (P_j^b | mu_i^0 nu_j^b) = (P_j^0 | nu_j^0 mu_i^-b),
            ! fetched from the opposite image exploiting the (mu nu|P) pair symmetry;
            ! blk(nu, mu, p) with a2(ni+1:, p; mu + (nu-1)*nAOi) = blk(nu, mu, p)
            mb_img = img_to_mimg(b_img)
            ind = [jatom, iatom, jatom]
            IF (mb_img > 0) THEN
               CALL get_3c_or_zero(t_3c_pari(mb_img, 1), ind, [nAOj, nAOi, nj], blk_3c)
            ELSE
               ALLOCATE (blk_3c(nAOj, nAOi, nj))
               blk_3c = 0.0_dp
            END IF
            a2(ni + 1:ndom, :) = TRANSPOSE(RESHAPE( &
                                           RESHAPE(blk_3c, [nAOi, nAOj, nj], order=[2, 1, 3]), [nAOi*nAOj, nj]))
            DEALLOCATE (blk_3c)
         END IF

         ! dead-pair skip, tied to the tensor filter (two orders above it, the
         ! measured gap between overlap-dead and live 3c blocks; scale-aware
         ! when the user tightens EPS_FILTER)
         IF (SQRT(SUM(a2**2)) < 1.0E2_dp*ri_data%filter_eps .OR. SUM(ABS(img_cells(:, b_img))) >= pari_remote_cut) THEN
            ! dead or remote pair: the 3c products are overlap-dead (exactly
            ! zero, or pure integral-engine roundoff on a physically empty
            ! pair), or the pair cell lies beyond any orbital overlap at this
            ! range separation; feeding roundoff into the truncated inverse
            ! amplifies it by the near-null modes into spurious pool entries,
            ! so such pairs are skipped and the pooled C3 stays zero
            c_full = 0.0_dp
         ELSE
            CALL pari_pair_pinv(vdom, pari_eps_eig, vinv, rcond)
!$OMP CRITICAL(pari_rcond_max)
            ri_data%pari_rcond_max = MAX(ri_data%pari_rcond_max, rcond)
!$OMP END CRITICAL(pari_rcond_max)
            c_full = MATMUL(vinv, a2)
            ! fit-residual probe: ||A - V C||_F / ||A||_F of the regularized
            ! pair solve, before the same-site / same-atom scalings; grows when
            ! the pair domain cannot represent A (kernel too long-ranged for
            ! the home-only domains) -- cheap failure detector
            fit_res = SQRT(SUM((a2 - MATMUL(vdom, c_full))**2))/MAX(SQRT(SUM(a2**2)), 1.0E-30_dp)
!$OMP CRITICAL(pari_res_probe)
            IF (fit_res > res_max) res_max = fit_res
            IF (fit_res > 1.0E-2_dp) THEN
               WRITE (*, '(A,3I5,ES10.3)') " PARI| fit residual >1e-2, pair (i,j,b): ", &
                  iatom, jatom, b_img, fit_res
            END IF
!$OMP END CRITICAL(pari_res_probe)
            IF (same_site) c_full = 0.5_dp*c_full
            IF (iatom == jatom .AND. b_img /= 1) c_full = pari_sac_f*c_full
         END IF
         ri_data%pari_c3(offset + 1:offset + ni*ncol) = RESHAPE(c_full(1:ni, :), [ni*ncol])
         !TODO(dev dump, remove before PR): main-cell pairs only (the pyscf
         ! oracle gate); raw per-pair locals -- single-rank runs only
         IF (para_env%num_pe == 1) THEN
!$OMP CRITICAL(pari_c3_dump)
            IF (SUM(ABS(img_cells(:, b_img))) <= 1) THEN
               WRITE (unit_dump, *) iatom, jatom, b_img, ni, nAOi, nAOj, MERGE(1, 0, same_site)
               WRITE (unit_dump, '(ES25.16)') vdom
               WRITE (unit_dump, '(ES25.16)') c_full
               WRITE (unit_dump, '(ES25.16)') a2
            END IF
            WRITE (unit_pool, *) iatom, jatom, b_img, ni, nAOi, nAOj
            WRITE (unit_pool, '(ES25.16)') c_full(1:ni, :)
!$OMP END CRITICAL(pari_c3_dump)
         END IF
         DEALLOCATE (vdom, a2, c_full, vinv)
      END DO
!$OMP END PARALLEL
      ! fit-residual report (probe details above inside the solve loop);
      ! residuals of 1e-3..1e-1 are the expected Levenberg-Marquardt bias of
      ! the regularized solve (they track pari_eps_eig); O(1) residuals mean
      ! the home-only pair domains cannot span this system / interaction range
      WRITE (*, '(A,ES10.3)') " PARI| pair fit residual, worst over domains: ", res_max
      IF (res_max > 5.0E-1_dp) THEN
         CALL cp_warn(__LOCATION__, "PARI pair-domain fit residuals are O(1): the home-only "// &
                      "domains cannot span this system / interaction range")
      END IF
      IF (para_env%num_pe == 1) THEN
         CLOSE (unit_dump)
         CLOSE (unit_pool)
      ELSE
         ! MPI: reduce the rank-split pool; then the last rank walks the
         ! replicated pool and writes the regression dump (same bytes as
         ! the single-rank in-loop writer above)
         CALL para_env%sum(ri_data%pari_c3)
         CALL para_env%max(ri_data%pari_rcond_max)
         IF (para_env%mepos == para_env%num_pe - 1) THEN
            OPEN (NEWUNIT=unit_pool, FILE="PARI_POOL_DUMP")
            DO ipair = 1, npair
               iatom = pair_i(ipair); jatom = pair_j(ipair); b_img = pair_img(ipair)
               ni = ri_data%bsizes_RI(iatom)
               nAOi = ri_data%bsizes_AO(iatom); nAOj = ri_data%bsizes_AO(jatom)
               ncol = nAOi*nAOj
               offset = ri_data%pari_c3_offset(ipair) - 1
               ALLOCATE (c_full(ni, ncol))
               c_full = RESHAPE(ri_data%pari_c3(offset + 1:offset + ni*ncol), [ni, ncol])
               WRITE (unit_pool, *) iatom, jatom, b_img, ni, nAOi, nAOj
               WRITE (unit_pool, '(ES25.16)') c_full
               DEALLOCATE (c_full)
            END DO
            CLOSE (unit_pool)
         END IF
      END IF

      DEALLOCATE (img_to_mimg, img_cells)
      CALL timestop(handle)

   END SUBROUTINE pari_solve_pairs

! **************************************************************************************************
!> \brief ...
!> \param tensor ...
!> \param ind block indices (own atom, partner atom, aux atom)
!> \param dims expected block shape
!> \param blk block values
! **************************************************************************************************
   SUBROUTINE get_3c_or_zero(tensor, ind, dims, blk)
      TYPE(dbt_type), INTENT(INOUT)                      :: tensor
      INTEGER, DIMENSION(3), INTENT(IN)                  :: ind, dims
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: blk

      LOGICAL                                            :: found

      ALLOCATE (blk(dims(1), dims(2), dims(3)))
      blk = 0.0_dp
      CALL dbt_get_block(tensor, ind, blk, found)
      IF (.NOT. found) blk = 0.0_dp

   END SUBROUTINE get_3c_or_zero

END MODULE hfx_ri_kp_pari
