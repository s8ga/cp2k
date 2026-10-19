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
   USE cp_dbcsr_api,                    ONLY: dbcsr_copy, dbcsr_create, dbcsr_distribution_release, &
                                              dbcsr_distribution_type, dbcsr_get_block_p, &
                                              dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
                                              dbcsr_iterator_start, dbcsr_iterator_stop, &
                                              dbcsr_iterator_type, dbcsr_p_type, dbcsr_release, &
                                              dbcsr_replicate_all, dbcsr_type, dbcsr_type_no_symmetry
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_dist2d_to_dist
   USE dbt_api,                          ONLY: dbt_create, dbt_destroy, dbt_distribution_destroy, &
                                              dbt_distribution_new, dbt_distribution_type, &
                                              dbt_filter, dbt_get_block, dbt_mp_environ_pgrid, &
                                              dbt_pgrid_create, dbt_pgrid_destroy, dbt_pgrid_type, &
                                              dbt_type
   USE distribution_2d_types,           ONLY: distribution_2d_type
   USE hfx_ri_kp,                        ONLY: get_kp_and_ri_images
   USE hfx_types,                        ONLY: hfx_ri_type
   USE input_cp2k_hfx,                   ONLY: ri_pmat
   USE kinds,                            ONLY: dp
   USE kpoint_types,                     ONLY: get_kpoint_info, kpoint_type
   USE machine,                          ONLY: m_flush
   USE mathlib,                          ONLY: get_pseudo_inverse_diag
   USE message_passing,                  ONLY: mp_cart_type, mp_para_env_type
   USE particle_types,                   ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env, qs_environment_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_kind_types,                    ONLY: qs_kind_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type, release_neighbor_list_sets
   USE qs_tensors,                       ONLY: build_2c_integrals, build_2c_neighbor_lists, &
                                              build_3c_integrals, build_3c_neighbor_lists, &
                                              neighbor_list_3c_destroy
   USE qs_tensors_types,                 ONLY: create_3c_tensor, distribution_3d_create, &
                                              distribution_3d_type, neighbor_list_3c_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = 'hfx_ri_kp_pari'

   ! relative eigenvalue threshold of the per-pair truncated pseudo-inverse,
   ! oracle-calibrated (PHASE2_SLICE2B appendix 10: pair-domain condition numbers
   ! ~1e10 for Si union bases; the exchange matrix is insensitive for thresholds
   ! anywhere in [1e-8, no truncation])
   REAL(KIND=dp), PARAMETER, PRIVATE :: pari_eps_eig = 1.0E-10_dp

   PUBLIC :: hfx_ri_update_ks_kp_pari, pari_pair_pinv

CONTAINS

! **************************************************************************************************
!> \brief Truncated pseudo-inverse of a dense symmetric pair-domain RI metric.
!>        Always the eigendecomposition pseudo-inverse with a relative eigenvalue
!>        threshold (mathlib get_pseudo_inverse_diag), mirroring the reference
!>        semantics (symmetric eigendecomposition, keep w > eps_eig*wmax).
!>        A direct Cholesky inverse is deliberately not attempted: the PARI union
!>        basis is intentionally redundant (pair-domain condition numbers of order
!>        1e10 measured for Si), where a direct inverse deviates O(1) from the
!>        truncated solve in the near-null directions, while the contraction
!>        coefficients stay stable only under truncation.
!> \param vmat dense symmetric pair-domain metric (P|Q), untouched on exit
!> \param eps_eig relative eigenvalue threshold (oracle-calibrated default 1.0E-10)
!> \param vinv truncated pseudo-inverse, allocated by the caller (shape of vmat)
!> \param rcond reciprocal condition estimate 1/(||V||_1 ||V^+||_1), for logging
! **************************************************************************************************
   SUBROUTINE pari_pair_pinv(vmat, eps_eig, vinv, rcond)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: vmat
      REAL(KIND=dp), INTENT(IN)                       :: eps_eig
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: vinv
      REAL(KIND=dp), INTENT(OUT)                      :: rcond

      CHARACTER(LEN=*), PARAMETER :: routineN = 'pari_pair_pinv'

      INTEGER                                            :: handle, n
      REAL(KIND=dp)                                      :: ainv_norm, a_norm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)          :: work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)       :: vwork
      REAL(KIND=dp), EXTERNAL                            :: dlange

      CALL timeset(routineN, handle)

      n = SIZE(vmat, 1)
      CPASSERT(SIZE(vmat, 2) == n)
      CPASSERT(SIZE(vinv, 1) == n .AND. SIZE(vinv, 2) == n)
      CPASSERT(eps_eig > 0.0_dp)

      ! get_pseudo_inverse_diag destroys its input argument
      ALLOCATE (vwork(n, n), work(n))
      vwork(:, :) = vmat

      CALL get_pseudo_inverse_diag(vwork, vinv, eps_eig)

      ! condition estimate from the 1-norms of V and V^+ (the truncated spectrum
      ! makes dpocon-style bounds meaningless, so a norm product is used instead)
      a_norm = dlange('1', n, n, vmat, n, work)
      ainv_norm = dlange('1', n, n, vinv, n, work)
      rcond = 1.0_dp/(a_norm*ainv_norm)

      DEALLOCATE (vwork, work)

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

      CPASSERT(nspins > 0)
      CPASSERT(hf_fraction >= 0.0_dp)

      CALL timeset(routineN, handle)

      IF (geometry_did_change) THEN
         CALL cp_warn(__LOCATION__, "PARI k-points: SCF-side KS assembly not implemented yet, "// &
                      "exchange contributions are zero")
         CALL hfx_ri_pre_scf_kp_pari(ks_matrix(1, 1)%matrix, ri_data, qs_env)
      END IF

      !TODO(B3): four-term KS assembly from the cached C3 coefficients and the
      ! per-image 2c potential, mirroring the hfx_ri_kp.F contraction structure
      ehfx = 0.0_dp

      CALL timestop(handle)

   END SUBROUTINE hfx_ri_update_ks_kp_pari

! **************************************************************************************************
!> \brief PARI pre-SCF: build the per-image 2c metric, stream the 3c integrals
!>        through per-pair domain solves, and cache the pooled fit coefficients.
!>        Twin of hfx_ri_pre_scf_kp (hfx_ri_kp.F), but the RI direction is never
!>        extended over ncell_RI cells: all aux functions live in their main cell,
!>        and the per-atom global (bumped) metric inverse is replaced by the
!>        per-pair truncated pseudo-inverse.
!> \param dbcsr_template ...
!> \param ri_data ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE hfx_ri_pre_scf_kp_pari(dbcsr_template, ri_data, qs_env)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: dbcsr_template
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ri_pre_scf_kp_pari'

      INTEGER :: handle, i_img, iatom, jatom, natom, nkind, nimg, npair
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: img_to_RI_cell_pari, pair_i, pair_img, &
                                                            pair_j
      INTEGER, DIMENSION(:), POINTER                     :: col_bsize, row_bsize
      INTEGER, DIMENSION(3)                              :: pdims, pcoord
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk_2c
      TYPE(dbcsr_iterator_type)                          :: iter_2c
      TYPE(dbcsr_distribution_type)                      :: dbcsr_dist
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: t_2c_pari
      TYPE(dbt_type)                                     :: t_3c_dummy
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :)       :: t_3c_pari
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_2d_type), POINTER                :: dist_2d
      TYPE(distribution_3d_type)                         :: dist_3d
      TYPE(gto_basis_set_p_type), ALLOCATABLE, DIMENSION(:) :: basis_set_AO, basis_set_RI
      TYPE(mp_cart_type)                                 :: mp_comm_t3c
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: nl_2c
      TYPE(neighbor_list_3c_type)                         :: nl_3c
      TYPE(dbt_pgrid_type)                               :: pgrid
      TYPE(dbt_distribution_type)                        :: t_dist
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dist_1, dist_2, dist_3

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
      IF (para_env%num_pe > 1) THEN
         CPABORT("RI%KP_METHOD PARI is currently restricted to a single MPI rank")
      END IF

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

      CALL build_2c_integrals(t_2c_pari, ri_data%filter_eps_2c, qs_env, nl_2c, basis_set_RI, &
                              basis_set_RI, ri_data%hfx_pot, do_kpoints=.TRUE., do_hfx_kpoints=.TRUE.)
      CALL release_neighbor_list_sets(nl_2c)

      ! replicated so that every rank can harvest per-pair blocks (small natom^2 objects)
      DO i_img = 1, nimg
         CALL dbcsr_replicate_all(t_2c_pari(i_img))
      END DO

      ! ---------------------------------------------------------------- 3c: (mu^0 nu^b | P^0) per image
      ! plain (non-extended) RI direction: img_to_RI_cell maps every image to the
      ! main cell, so the builder folds all aux functions into main-cell blocks
      ALLOCATE (img_to_RI_cell_pari(nimg))
      img_to_RI_cell_pari(:) = 1

      pdims = 0
      CALL dbt_pgrid_create(para_env, pdims, pgrid, tensor_dims=[natom, natom, 1])

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
      CALL dbt_distribution_new(t_dist, pgrid, dist_1, dist_2, dist_3)
      CALL dbt_create(t_3c_pari(1, 1), "PARI_3c_ints", t_dist, [1, 2], [3], &
                      ri_data%bsizes_AO, ri_data%bsizes_AO, ri_data%bsizes_RI)
      DO i_img = 2, nimg
         CALL dbt_create(t_3c_pari(1, 1), t_3c_pari(i_img, 1))
      END DO
      CALL dbt_distribution_destroy(t_dist)
      DEALLOCATE (dist_1, dist_2, dist_3)

      CALL build_3c_integrals(t_3c_pari, ri_data%filter_eps, qs_env, nl_3c, basis_set_AO, &
                              basis_set_AO, basis_set_RI, ri_data%ri_metric, &
                              int_eps=ri_data%eps_schwarz, op_pos=2, desymmetrize=.FALSE., &
                              do_kpoints=.TRUE., do_hfx_kpoints=.TRUE., &
                              RI_range=ri_data%kp_RI_range, img_to_RI_cell=img_to_RI_cell_pari)
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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'pari_solve_pairs'

      INTEGER :: b_img, cell(3), handle, i, iatom, ipair, jatom, l, mb_img, nAOi, nAOj, ncol, ndom, &
                 ni, nj, npair, offset, unit_dump
      INTEGER, DIMENSION(3)                              :: ind
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: img_to_mimg
      LOGICAL                                            :: found, same_site
      REAL(KIND=dp)                                      :: rcond
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: a2, c_full, vdom, vinv
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: blk_3c
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: vij, vii, vjj
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
      DO b_img = 1, ri_data%nimg
         cell = 0
         DO i = LBOUND(cell_to_index, 1), UBOUND(cell_to_index, 1)
            DO jatom = LBOUND(cell_to_index, 2), UBOUND(cell_to_index, 2)
               DO l = LBOUND(cell_to_index, 3), UBOUND(cell_to_index, 3)
                  IF (cell_to_index(i, jatom, l) == b_img) cell = [i, jatom, l]
               END DO
            END DO
         END DO
         IF (ANY(cell /= 0)) THEN
            i = cell_to_index(-cell(1), -cell(2), -cell(3))
            IF (i >= 1 .AND. i <= ri_data%nimg) img_to_mimg(b_img) = i
         END IF
      END DO

      ALLOCATE (ri_data%pari_c3_offset(npair + 1))
      ri_data%pari_c3_offset(1) = 1
      DO ipair = 1, npair
         ni = ri_data%bsizes_RI(pair_i(ipair))
         ri_data%pari_c3_offset(ipair + 1) = ri_data%pari_c3_offset(ipair) + &
                                             ni*ri_data%bsizes_AO(pair_i(ipair))*ri_data%bsizes_AO(pair_j(ipair))
      END DO
      ALLOCATE (ri_data%pari_c3(ri_data%pari_c3_offset(npair + 1) - 1))
      ALLOCATE (ri_data%pari_nrow_i(npair))
      ri_data%pari_nrow_i(:) = ri_data%bsizes_RI(pair_i(:))

      !TODO(dev dump, remove before PR): C3 pool dump for the pyscf oracle gate
      IF (para_env%mepos == para_env%num_pe - 1) OPEN (NEWUNIT=unit_dump, FILE="PARI_C3_DUMP")
      DO ipair = 1, npair
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

         CALL pari_pair_pinv(vdom, pari_eps_eig, vinv, rcond)
         ri_data%pari_rcond_max = MAX(ri_data%pari_rcond_max, rcond)
         c_full = MATMUL(vinv, a2)
         IF (same_site) c_full = 0.5_dp*c_full
         ri_data%pari_c3(offset + 1:offset + ni*ncol) = RESHAPE(c_full(1:ni, :), [ni*ncol])
         !TODO(dev dump, remove before PR): main-cell pairs only (the pyscf oracle
         ! gate compares the molecular, i.e. same-cell, pair domains)
         IF (para_env%mepos == para_env%num_pe - 1 .AND. b_img == 1) THEN
            WRITE (unit_dump, *) iatom, jatom, b_img, ni, nAOi, nAOj, MERGE(1, 0, same_site)
            WRITE (unit_dump, '(ES25.16)') vdom
            WRITE (unit_dump, '(ES25.16)') c_full(1:ni, :)
         END IF
         DEALLOCATE (vdom, a2, c_full, vinv)
      END DO
      IF (para_env%mepos == para_env%num_pe - 1) CLOSE (unit_dump)

      DEALLOCATE (img_to_mimg)
      CALL timestop(handle)

   END SUBROUTINE pari_solve_pairs

! **************************************************************************************************
!> \brief fetch a 3c block; zero block if absent (below the integral filter)
!> \param tensor ...
!> \param ind block indices (own atom, partner atom, aux atom)
!> \param dims expected block shape
!> \param blk block values
! **************************************************************************************************
   SUBROUTINE get_3c_or_zero(tensor, ind, dims, blk)
      TYPE(dbt_type), INTENT(INOUT)                      :: tensor
      INTEGER, DIMENSION(3), INTENT(IN)                  :: ind, dims
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), INTENT(OUT) :: blk

      LOGICAL                                            :: found

      ALLOCATE (blk(dims(1), dims(2), dims(3)))
      blk = 0.0_dp
      CALL dbt_get_block(tensor, ind, blk, found)
      IF (.NOT. found) blk = 0.0_dp

   END SUBROUTINE get_3c_or_zero

END MODULE hfx_ri_kp_pari
