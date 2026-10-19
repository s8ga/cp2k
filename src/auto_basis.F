!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   Automatic generation of auxiliary basis sets of different kind
!> \author  JGH
!>
!> <b>Modification history:</b>
!> - 11.2017 creation [JGH]
! **************************************************************************************************
MODULE auto_basis
   USE ai_contraction_sphi,           ONLY: ab_contract
   USE aux_basis_set,                 ONLY: create_aux_basis
   USE basis_set_types,               ONLY: combine_basis_sets,&
                                              deallocate_gto_basis_set,&
                                              get_gto_basis_set,&
                                              gto_basis_set_type,&
                                              init_orb_basis_set,&
                                              sort_gto_basis_set
   USE bibliography,                  ONLY: Stoychev2016,&
                                              cite_reference
   USE gamma,                         ONLY: init_md_ftable
   USE kinds,                         ONLY: default_string_length,&
                                              dp
   USE libint_2c_3c,                  ONLY: eri_2center,&
                                              libint_potential_type
   USE libint_wrapper,                ONLY: cp_libint_cleanup_2eri,&
                                              cp_libint_init_2eri,&
                                              cp_libint_set_contrdepth,&
                                              cp_libint_t
   USE mathconstants,                 ONLY: dfac,&
                                              fac,&
                                              gamma1,&
                                              pi,&
                                              rootpi
   USE orbital_pointers,              ONLY: init_orbital_pointers,&
                                              nco,&
                                              ncoset,&
                                              nso
   USE periodic_table,                ONLY: get_ptable_info
   USE qs_kind_types,                 ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'auto_basis'

   PUBLIC :: create_ri_aux_basis_set, create_lri_aux_basis_set, &
             create_oce_basis, create_pari_aux_basis_set

CONTAINS

! **************************************************************************************************
!> \brief Create a RI_AUX basis set using some heuristics
!> \param ri_aux_basis_set ...
!> \param qs_kind ...
!> \param basis_cntrl ...
!> \param basis_type ...
!> \param basis_sort ...
!> \date    01.11.2017
!> \author  JGH
! **************************************************************************************************
   SUBROUTINE create_ri_aux_basis_set(ri_aux_basis_set, qs_kind, basis_cntrl, basis_type, basis_sort)
      TYPE(gto_basis_set_type), POINTER                  :: ri_aux_basis_set
      TYPE(qs_kind_type), INTENT(IN)                     :: qs_kind
      INTEGER, INTENT(IN)                                :: basis_cntrl
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: basis_type
      INTEGER, INTENT(IN), OPTIONAL                      :: basis_sort

      CHARACTER(LEN=2)                                   :: element_symbol
      CHARACTER(LEN=default_string_length)               :: bsname, kname
      INTEGER                                            :: i, j, jj, l, laux, linc, lmax, lval, lx, &
                                                            nsets, nx, z
      INTEGER, DIMENSION(0:18)                           :: nval
      INTEGER, DIMENSION(0:9, 1:20)                      :: nl
      INTEGER, DIMENSION(1:3)                            :: ls1, ls2, npgf
      INTEGER, DIMENSION(:), POINTER                     :: econf
      REAL(KIND=dp)                                      :: xv, zval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: zet
      REAL(KIND=dp), DIMENSION(0:18)                     :: bv, bval, fv, peff, pend, pmax, pmin
      REAL(KIND=dp), DIMENSION(0:9)                      :: zeff, zmax, zmin
      REAL(KIND=dp), DIMENSION(3)                        :: amax, amin, bmin
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set

      !
      CALL cite_reference(Stoychev2016)
      !
      bv(0:18) = [1.8_dp, 2.0_dp, 2.2_dp, 2.2_dp, 2.3_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, &
                  3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp]
      fv(0:18) = [20.0_dp, 4.0_dp, 4.0_dp, 3.5_dp, 2.5_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, &
                  2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp]
      !
      CPASSERT(.NOT. ASSOCIATED(ri_aux_basis_set))
      NULLIFY (orb_basis_set, econf)
      IF (.NOT. PRESENT(basis_type)) THEN
         CALL get_qs_kind(qs_kind, basis_set=orb_basis_set, basis_type="ORB")
      ELSE
         CALL get_qs_kind(qs_kind, basis_set=orb_basis_set, basis_type=basis_type)
      END IF
      IF (ASSOCIATED(orb_basis_set)) THEN
         ! BASIS_SET ORB NONE associates the pointer orb_basis_set, but does not contain
         ! any actual basis functions. Therefore, we catch it here to avoid spurious autogenerated
         ! RI_AUX basis sets.
         IF (SUM(orb_basis_set%nsgf_set) == 0) THEN
            CALL cp_abort(__LOCATION__, &
                          "Cannot autocreate RI_AUX basis set for at least one of the given "// &
                          "primary basis sets due to missing exponents. If you have invoked BASIS_SET NONE, "// &
                          "you should state BASIS_SET RI_AUX NONE explicitly in the input.")
         END IF
         CALL get_basis_keyfigures(orb_basis_set, lmax, zmin, zmax, zeff)
         !Note: RI basis coud require lmax up to 2*orb_lmax. This ensures that all orbital pointers
         !      are properly initialized before building the basis
         CALL init_orbital_pointers(2*lmax)
         CALL get_basis_products(lmax, zmin, zmax, zeff, pmin, pmax, peff)
         CALL get_qs_kind(qs_kind, zeff=zval, elec_conf=econf, element_symbol=element_symbol)
         IF (.NOT. ASSOCIATED(econf)) THEN
            CALL get_qs_kind(qs_kind, name=kname)
            CALL cp_abort(__LOCATION__, &
                          "AUTO_BASIS RI_AUX cannot process atom kind "// &
                          "<"//TRIM(ADJUSTL(kname))//"> due to missing "// &
                          "definition of potential or electron configuration; "// &
                          "consider setting keyword ELEC_CONF explicitly for "// &
                          "GHOST atom kind that has assigned a basis set")
         END IF
         CALL get_ptable_info(element_symbol, ielement=z)
         lval = 0
         DO l = 0, MAXVAL(UBOUND(econf))
            IF (econf(l) > 0) lval = l
         END DO
         IF (SUM(econf) /= NINT(zval)) THEN
            CPWARN("Valence charge and electron configuration not consistent")
         END IF
         pend = 0.0_dp
         linc = 1
         IF (z > 18) linc = 2
         SELECT CASE (basis_cntrl)
         CASE (0)
            laux = MAX(2*lval, lmax + linc)
         CASE (1)
            laux = MAX(2*lval, lmax + linc)
         CASE (2)
            laux = MAX(2*lval, lmax + linc + 1)
         CASE (3)
            laux = MAX(2*lmax, lmax + linc + 2)
         CASE DEFAULT
            CPABORT("Invalid value of control variable")
         END SELECT
         !
         DO l = 2*lmax + 1, laux
            xv = peff(2*lmax)
            pmin(l) = xv
            pmax(l) = xv
            peff(l) = xv
            pend(l) = xv
         END DO
         !
         DO l = 0, laux
            IF (l <= 2*lval) THEN
               pend(l) = MIN(fv(l)*peff(l), pmax(l))
               bval(l) = 1.8_dp
            ELSE
               pend(l) = peff(l)
               bval(l) = bv(l)
            END IF
            xv = LOG(pend(l)/pmin(l))/LOG(bval(l)) + 1.e-10_dp
            nval(l) = MAX(CEILING(xv), 0)
         END DO
         ! first set include valence only
         nsets = 1
         ls1(1) = 0
         ls2(1) = lval
         DO l = lval + 1, laux
            IF (nval(l) < nval(lval) - 1) EXIT
            ls2(1) = l
         END DO
         ! second set up to 2*lval
         IF (laux > ls2(1)) THEN
            IF (lval == 0 .OR. 2*lval <= ls2(1) + 1) THEN
               nsets = 2
               ls1(2) = ls2(1) + 1
               ls2(2) = laux
            ELSE
               nsets = 2
               ls1(2) = ls2(1) + 1
               ls2(2) = MIN(2*lval, laux)
               lx = ls2(2)
               DO l = lx + 1, laux
                  IF (nval(l) < nval(lx) - 1) EXIT
                  ls2(2) = l
               END DO
               IF (laux > ls2(2)) THEN
                  nsets = 3
                  ls1(3) = ls2(2) + 1
                  ls2(3) = laux
               END IF
            END IF
         END IF
         !
         amax = 0.0_dp
         amin = HUGE(0.0_dp)
         bmin = HUGE(0.0_dp)
         DO i = 1, nsets
            DO j = ls1(i), ls2(i)
               amax(i) = MAX(amax(i), pend(j))
               amin(i) = MIN(amin(i), pmin(j))
               bmin(i) = MIN(bmin(i), bval(j))
            END DO
            xv = LOG(amax(i)/amin(i))/LOG(bmin(i)) + 1.e-10_dp
            npgf(i) = MAX(CEILING(xv), 0)
         END DO
         nx = MAXVAL(npgf(1:nsets))
         ALLOCATE (zet(nx, nsets))
         zet = 0.0_dp
         nl = 0
         DO i = 1, nsets
            DO j = 1, npgf(i)
               jj = npgf(i) - j + 1
               zet(jj, i) = amin(i)*bmin(i)**(j - 1)
            END DO
            DO l = ls1(i), ls2(i)
               nl(l, i) = nval(l)
            END DO
         END DO
         bsname = TRIM(element_symbol)//"-RI-AUX-"//TRIM(orb_basis_set%name)
         !
         CALL create_aux_basis(ri_aux_basis_set, bsname, nsets, ls1, ls2, nl, npgf, zet)

         DEALLOCATE (zet)

         IF (PRESENT(basis_sort)) THEN
            CALL sort_gto_basis_set(ri_aux_basis_set, basis_sort)
         END IF

      END IF

   END SUBROUTINE create_ri_aux_basis_set
! **************************************************************************************************
!> \brief Create a RI basis set for pair atomic domain RI (PARI): the tempered RI_AUX
!>        basis combined with same-center Gaussian-product functions of the ORB basis,
!>        selected by principal component analysis (PCA) of the two-center metric. The
!>        rotationally invariant metric is block diagonal in the total angular momentum,
!>        so the PCA factorizes into small radial eigenproblems, one per angular
!>        momentum. Radial eigenvectors with eigenvalues above eps_pca times the largest
!>        eigenvalue are kept as contracted shells of pure angular momentum.
!> \param ri_aux_basis_set ...
!> \param qs_kind ...
!> \param basis_cntrl ...
!> \param basis_type ...
!> \param basis_sort ...
!> \param potential_parameter two-center metric operator (e.g. short range)
!> \param eps_pca relative PCA threshold for keeping product functions
!> \date    10.2026
! **************************************************************************************************
   SUBROUTINE create_pari_aux_basis_set(ri_aux_basis_set, qs_kind, basis_cntrl, basis_type, &
                                        basis_sort, potential_parameter, eps_pca)
      TYPE(gto_basis_set_type), POINTER                  :: ri_aux_basis_set
      TYPE(qs_kind_type), INTENT(IN)                     :: qs_kind
      INTEGER, INTENT(IN)                                :: basis_cntrl
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: basis_type
      INTEGER, INTENT(IN), OPTIONAL                      :: basis_sort
      TYPE(libint_potential_type), INTENT(IN)            :: potential_parameter
      REAL(KIND=dp), INTENT(IN)                          :: eps_pca

      CHARACTER(LEN=2)                                   :: element_symbol
      CHARACTER(LEN=512)                                 :: error_msg
      CHARACTER(LEN=default_string_length)               :: bsname
      CHARACTER(len=*), PARAMETER                        :: routineN = 'create_pari_aux_basis_set'

      INTEGER :: handle, i, iL, info, ii, iOrbSet, iOrbSh, isgf0, j, jOrbSet, &
                 jOrbSh, k, kL, la, lb, lmax_orb, Ltot, m, m0, maxL, maxne, &
                 msum, ne, nkeep_tot, nsetF, nsetL, nset_o, nshell_tot, nsoL, ntot
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kepF, kepL, keep_idx, Llist, ls1, &
                                                            ls1_f, ls2, ls2_f, neL, npgf_b, &
                                                            npgf_f, sh_l, sh_set
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: nl_b, nl_f
      REAL(KIND=dp)                                      :: eval_max, xtmp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: evals, rpgf0, sums, work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ev_all, exp_all, gpp, gsph, sij, &
                                                            sph, zet_b, zet_f
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: vec_all
      REAL(KIND=dp), DIMENSION(3)                        :: r0
      TYPE(cp_libint_t)                                  :: lib
      TYPE(gto_basis_set_type), POINTER                  :: boot_basis, orb_basis_set, &
                                                            pca_basis, tempered

! -------------------------------------------------------------------------

      CALL timeset(routineN, handle)
      CPASSERT(.NOT. ASSOCIATED(ri_aux_basis_set))
      NULLIFY (tempered, orb_basis_set, pca_basis, boot_basis)
      r0 = 0.0_dp

      ! tempered part: identical to the RI_HFX auto basis
      CALL create_ri_aux_basis_set(tempered, qs_kind, basis_cntrl, basis_type)

      ! product part: built from the ORB basis of this kind
      CALL get_qs_kind(qs_kind, basis_set=orb_basis_set, basis_type="ORB")
      IF (.NOT. ASSOCIATED(orb_basis_set) .OR. SUM(orb_basis_set%nsgf_set) == 0) THEN
         IF (.NOT. ASSOCIATED(tempered)) THEN
            CALL cp_abort(__LOCATION__, "PARI auto basis: neither ORB products nor a "// &
                          "tempered RI basis are available for this kind")
         END IF
         CPWARN("PARI auto basis: no ORB basis found, falling back to the tempered basis")
         ri_aux_basis_set => tempered
         CALL timestop(handle)
         RETURN
      END IF

      nset_o = orb_basis_set%nset
      lmax_orb = MAXVAL(orb_basis_set%lmax(1:nset_o))
      IF (2*lmax_orb > 9) THEN
         CALL cp_abort(__LOCATION__, "PARI auto basis: ORB basis l-max too large for "// &
                       "Gaussian-product expansion")
      END IF
      maxL = 2*lmax_orb
      CALL init_orbital_pointers(maxL)

      ! flatten the ORB shell list
      nshell_tot = SUM(orb_basis_set%nshell(1:nset_o))
      ALLOCATE (sh_set(nshell_tot), sh_l(nshell_tot))
      k = 0
      DO iOrbSet = 1, nset_o
         DO iOrbSh = 1, orb_basis_set%nshell(iOrbSet)
            k = k + 1
            sh_set(k) = iOrbSet
            sh_l(k) = orb_basis_set%l(iOrbSh, iOrbSet)
         END DO
      END DO

      ! upper bound on the number of product exponents
      ntot = 0
      DO i = 1, nshell_tot
         DO j = i, nshell_tot
            ntot = ntot + orb_basis_set%npgf(sh_set(i))*orb_basis_set%npgf(sh_set(j))
         END DO
      END DO
      ALLOCATE (sums(ntot), exp_all(ntot, 0:maxL))

      ! distinct product exponents per total angular momentum Ltot (sorted, unique)
      ALLOCATE (neL(0:maxL))
      DO Ltot = 0, maxL
         m = 0
         DO i = 1, nshell_tot
            la = sh_l(i)
            DO j = i, nshell_tot
               lb = sh_l(j)
               IF (Ltot < ABS(la - lb) .OR. Ltot > la + lb) CYCLE
               IF (MOD(la + lb - Ltot, 2) /= 0) CYCLE
               DO k = 1, orb_basis_set%npgf(sh_set(i))
                  DO ii = 1, orb_basis_set%npgf(sh_set(j))
                     m = m + 1
                     sums(m) = orb_basis_set%zet(k, sh_set(i)) + orb_basis_set%zet(ii, sh_set(j))
                  END DO
               END DO
            END DO
         END DO
         ! insertion sort
         DO i = 2, m
            xtmp = sums(i)
            j = i - 1
            DO WHILE (j >= 1)
               IF (sums(j) <= xtmp) EXIT
               sums(j + 1) = sums(j)
               j = j - 1
            END DO
            sums(j + 1) = xtmp
         END DO
         ! compact to unique exponents (relative tolerance)
         msum = 0
         DO i = 1, m
            IF (msum == 0 .OR. sums(i) - sums(msum) > 1.0E-6_dp*sums(msum)) THEN
               msum = msum + 1
               sums(msum) = sums(i)
            END IF
         END DO
         neL(Ltot) = msum
         exp_all(1:msum, Ltot) = sums(1:msum)
      END DO
      maxne = MAXVAL(neL)

      nsetL = COUNT(neL > 0)
      CPASSERT(nsetL > 0)
      ALLOCATE (Llist(nsetL), ls1(nsetL), ls2(nsetL), npgf_b(nsetL), kepL(nsetL))
      ALLOCATE (zet_b(maxne, nsetL))
      zet_b = 0.0_dp
      k = 0
      DO Ltot = 0, maxL
         IF (neL(Ltot) > 0) THEN
            k = k + 1
            Llist(k) = Ltot
            ls1(k) = Ltot
            ls2(k) = Ltot
            npgf_b(k) = neL(Ltot)
            zet_b(1:neL(Ltot), k) = exp_all(1:neL(Ltot), Ltot)
         END IF
      END DO
      DEALLOCATE (exp_all, sums)

      ! bootstrap basis: one uncontracted shell per product primitive (identity gcc)
      ALLOCATE (nl_b(0:maxL, nsetL))
      nl_b = 0
      DO iL = 1, nsetL
         nl_b(Llist(iL), iL) = npgf_b(iL)
      END DO
      bsname = "PARI-BOOT"
      CALL create_aux_basis(boot_basis, bsname, nsetL, ls1, ls2, nl_b, npgf_b, zet_b)
      boot_basis%gcc = 0.0_dp
      DO iL = 1, nsetL
         DO i = 1, npgf_b(iL)
            boot_basis%gcc(i, i, iL) = 1.0_dp
         END DO
      END DO
      CALL init_orb_basis_set(boot_basis)

      ! radial two-center metric per Ltot and its eigendecomposition
      ALLOCATE (ev_all(maxne, nsetL), vec_all(maxne, maxne, nsetL))
      ALLOCATE (work(64*maxne + 2))
      eval_max = 0.0_dp
      CALL init_md_ftable(2*maxL)
      CALL cp_libint_init_2eri(lib, maxL)
      CALL cp_libint_set_contrdepth(lib, 1)
      DO iL = 1, nsetL
         Ltot = Llist(iL)
         ne = npgf_b(iL)
         nsoL = nso(Ltot)
         ii = ne*ncoset(Ltot)
         ALLOCATE (sij(ii, ii), sph(ii, ne*nsoL), gsph(ne*nsoL, ne*nsoL), gpp(ne, ne), &
                   evals(ne), rpgf0(ne))
         sij = 0.0_dp
         rpgf0 = 0.0_dp
         CALL eri_2center(sij, Ltot, Ltot, ne, zet_b(1:ne, iL), rpgf0, r0, &
                          Ltot, Ltot, ne, zet_b(1:ne, iL), rpgf0, r0, 0.0_dp, lib, potential_parameter)
         isgf0 = boot_basis%first_sgf(1, iL)
         sph = boot_basis%sphi(1:ii, isgf0:isgf0 + ne*nsoL - 1)
         CALL ab_contract(gsph, sij, sph, sph, ii, ii, ne*nsoL, ne*nsoL)
         ! rotationally invariant metric: any fixed m column gives the radial block
         m0 = 1
         DO i = 1, ne
            DO j = 1, ne
               gpp(i, j) = 0.5_dp*(gsph((i - 1)*nsoL + m0, (j - 1)*nsoL + m0) + &
                                   gsph((j - 1)*nsoL + m0, (i - 1)*nsoL + m0))
            END DO
         END DO
         DEALLOCATE (sij, gsph)
         ! eigendecomposition of the product metric
         gpp = 0.5_dp*(gpp + TRANSPOSE(gpp))
         CALL dsyev('V', 'U', ne, gpp, ne, evals, work, SIZE(work), info)
         CPASSERT(info == 0)
         ev_all(1:ne, iL) = evals(:)
         vec_all(1:ne, 1:ne, iL) = gpp(:, :)
         eval_max = MAX(eval_max, evals(ne))
         DEALLOCATE (sph, gpp, evals, rpgf0)
      END DO
      DEALLOCATE (work)
      CALL cp_libint_cleanup_2eri(lib)
      CALL deallocate_gto_basis_set(boot_basis)

      ! PCA selection: keep radial eigenvectors above eps_pca*eval_max
      nkeep_tot = 0
      DO iL = 1, nsetL
         ne = npgf_b(iL)
         kL = COUNT(ev_all(1:ne, iL) > eps_pca*eval_max)
         kepL(iL) = kL
         nkeep_tot = nkeep_tot + kL*nso(Llist(iL))
      END DO
      IF (nkeep_tot == 0) THEN
         WRITE (error_msg, '(A,ES9.1,A)') "PARI auto basis: no product function survived the "// &
            "PCA threshold ", eps_pca, "; keeping the tempered basis only"
         CPWARN(TRIM(error_msg))
      END IF

      CALL get_qs_kind(qs_kind, element_symbol=element_symbol)
      bsname = TRIM(element_symbol)//"-RI-AUX-PARI-"//TRIM(orb_basis_set%name)

      IF (nkeep_tot > 0) THEN
         ! final PCA basis: kept eigenvectors as contracted shells;
         ! pack only angular momenta with kepL > 0 (create_aux_basis cannot
         ! handle sets without contracted shells)
         nsetF = COUNT(kepL > 0)
         ALLOCATE (keep_idx(nsetF), ls1_f(nsetF), ls2_f(nsetF), npgf_f(nsetF), kepF(nsetF))
         ALLOCATE (nl_f(0:maxL, nsetF), zet_f(maxne, nsetF))
         nl_f = 0
         zet_f = 0.0_dp
         k = 0
         DO iL = 1, nsetL
            IF (kepL(iL) > 0) THEN
               k = k + 1
               keep_idx(k) = iL
               ls1_f(k) = ls1(iL)
               ls2_f(k) = ls2(iL)
               npgf_f(k) = npgf_b(iL)
               kepF(k) = kepL(iL)
               zet_f(1:npgf_f(k), k) = zet_b(1:npgf_b(iL), iL)
               nl_f(Llist(iL), k) = kepF(k)
            END IF
         END DO
         CALL create_aux_basis(pca_basis, bsname, nsetF, ls1_f, ls2_f, nl_f, npgf_f, zet_f)
         DO k = 1, nsetF
            iL = keep_idx(k)
            ne = npgf_b(iL)
            kL = kepL(iL)
            DO i = 1, kL
               pca_basis%gcc(1:ne, i, k) = vec_all(1:ne, ne - kL + i, iL)
            END DO
         END DO
         DEALLOCATE (keep_idx, ls1_f, ls2_f, npgf_f, kepF, nl_f, zet_f)
         IF (ASSOCIATED(tempered)) THEN
            CALL combine_basis_sets(tempered, pca_basis)
            ri_aux_basis_set => tempered
            CALL deallocate_gto_basis_set(pca_basis)
         ELSE
            ri_aux_basis_set => pca_basis
         END IF
         ri_aux_basis_set%name = bsname
         ri_aux_basis_set%aliases = bsname
         WRITE (error_msg, '(A,I0,A)') "Automatic generation of PARI auxiliary basis "// &
            "(tempered RI_AUX union on-site product PCA, experimental): added ", nkeep_tot, &
            " product functions"
         CPWARN(TRIM(error_msg))
      ELSE
         ri_aux_basis_set => tempered
      END IF

      ! Coefficients are emitted in the basis-file convention (contractions of
      ! normalized primitives): the two-center metric used above for the
      ! Cholesky/eigendecomposition is exactly the overlap of solid-harmonic
      ! normalized primitives, and init_orb_basis_set folds the primitive
      ! normalization into gcc when the basis object is initialised, so the
      ! runtime basis equals the set constructed here. The basis-file printer
      ! divides the fold back out when writing the dump.

      IF (PRESENT(basis_sort)) THEN
         CALL sort_gto_basis_set(ri_aux_basis_set, basis_sort)
      END IF

      DEALLOCATE (ev_all, kepL, Llist, ls1, ls2, neL, nl_b, npgf_b, sh_l, sh_set, vec_all, zet_b)

      CALL timestop(handle)

   END SUBROUTINE create_pari_aux_basis_set
! **************************************************************************************************
!> \brief Create a LRI_AUX basis set using some heuristics
!> \param lri_aux_basis_set ...
!> \param qs_kind ...
!> \param basis_cntrl ...
!> \param exact_1c_terms ...
!> \param tda_kernel ...
!> \date    01.11.2017
!> \author  JGH
! **************************************************************************************************
   SUBROUTINE create_lri_aux_basis_set(lri_aux_basis_set, qs_kind, basis_cntrl, &
                                       exact_1c_terms, tda_kernel)
      TYPE(gto_basis_set_type), POINTER                  :: lri_aux_basis_set
      TYPE(qs_kind_type), INTENT(IN)                     :: qs_kind
      INTEGER, INTENT(IN)                                :: basis_cntrl
      LOGICAL, INTENT(IN), OPTIONAL                      :: exact_1c_terms, tda_kernel

      CHARACTER(LEN=2)                                   :: element_symbol
      CHARACTER(LEN=default_string_length)               :: bsname, kname
      INTEGER                                            :: i, j, l, laux, linc, lm, lmax, lval, n1, &
                                                            n2, nsets, z
      INTEGER, DIMENSION(0:18)                           :: nval
      INTEGER, DIMENSION(0:9, 1:50)                      :: nl
      INTEGER, DIMENSION(1:50)                           :: ls1, ls2, npgf
      INTEGER, DIMENSION(:), POINTER                     :: econf
      LOGICAL                                            :: e1terms, kernel_basis
      REAL(KIND=dp)                                      :: xv, zval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: zet
      REAL(KIND=dp), DIMENSION(0:18)                     :: bval, peff, pend, pmax, pmin
      REAL(KIND=dp), DIMENSION(0:9)                      :: zeff, zmax, zmin
      REAL(KIND=dp), DIMENSION(4)                        :: bv, bx
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set

      !
      IF (PRESENT(exact_1c_terms)) THEN
         e1terms = exact_1c_terms
      ELSE
         e1terms = .FALSE.
      END IF
      IF (PRESENT(tda_kernel)) THEN
         kernel_basis = tda_kernel
      ELSE
         kernel_basis = .FALSE.
      END IF
      IF (kernel_basis .AND. e1terms) THEN
         CALL cp_warn(__LOCATION__, "LRI Kernel basis generation will ignore exact 1C term option.")
      END IF
      !
      CPASSERT(.NOT. ASSOCIATED(lri_aux_basis_set))
      NULLIFY (orb_basis_set, econf)
      CALL get_qs_kind(qs_kind, basis_set=orb_basis_set, basis_type="ORB")
      IF (ASSOCIATED(orb_basis_set)) THEN
         CALL get_basis_keyfigures(orb_basis_set, lmax, zmin, zmax, zeff)
         CALL get_basis_products(lmax, zmin, zmax, zeff, pmin, pmax, peff)
         CALL get_qs_kind(qs_kind, zeff=zval, elec_conf=econf, element_symbol=element_symbol)
         IF (.NOT. ASSOCIATED(econf)) THEN
            CALL get_qs_kind(qs_kind, name=kname)
            CALL cp_abort(__LOCATION__, &
                          "AUTO_BASIS LRI_AUX cannot process atom kind "// &
                          "<"//TRIM(ADJUSTL(kname))//"> due to missing "// &
                          "definition of potential or electron configuration; "// &
                          "consider setting keyword ELEC_CONF explicitly for "// &
                          "GHOST atom kind that has assigned a basis set")
         END IF
         CALL get_ptable_info(element_symbol, ielement=z)
         lval = 0
         DO l = 0, MAXVAL(UBOUND(econf))
            IF (econf(l) > 0) lval = l
         END DO
         IF (SUM(econf) /= NINT(zval)) THEN
            CPWARN("Valence charge and electron configuration not consistent")
         END IF
         !
         linc = 1
         IF (z > 18) linc = 2
         pend = 0.0_dp
         IF (kernel_basis) THEN
            bv(1:4) = [3.20_dp, 2.80_dp, 2.40_dp, 2.00_dp]
            bx(1:4) = [4.00_dp, 3.50_dp, 3.00_dp, 2.50_dp]
            !
            SELECT CASE (basis_cntrl)
            CASE (0)
               laux = lval + 1
            CASE (1)
               laux = MAX(lval + 1, lmax)
            CASE (2)
               laux = MAX(lval + 2, lmax + 1)
            CASE (3)
               laux = MAX(lval + 3, lmax + 2)
               laux = MIN(laux, 2 + linc)
            CASE DEFAULT
               CPABORT("Invalid value of control variable")
            END SELECT
         ELSE
            bv(1:4) = [2.00_dp, 1.90_dp, 1.80_dp, 1.80_dp]
            bx(1:4) = [2.60_dp, 2.40_dp, 2.20_dp, 2.20_dp]
            !
            SELECT CASE (basis_cntrl)
            CASE (0)
               laux = MAX(2*lval, lmax + linc)
               laux = MIN(laux, 2 + linc)
            CASE (1)
               laux = MAX(2*lval, lmax + linc)
               laux = MIN(laux, 3 + linc)
            CASE (2)
               laux = MAX(2*lval, lmax + linc + 1)
               laux = MIN(laux, 4 + linc)
            CASE (3)
               laux = MAX(2*lval, lmax + linc + 1)
               laux = MIN(laux, 4 + linc)
            CASE DEFAULT
               CPABORT("Invalid value of control variable")
            END SELECT
         END IF
         !
         DO l = 2*lmax + 1, laux
            pmin(l) = pmin(2*lmax)
            pmax(l) = pmax(2*lmax)
            peff(l) = peff(2*lmax)
         END DO
         !
         nval = 0
         IF (exact_1c_terms) THEN
            DO l = 0, laux
               IF (l <= lval + 1) THEN
                  pend(l) = zmax(l) + 1.0_dp
                  bval(l) = bv(basis_cntrl + 1)
               ELSE
                  pend(l) = 2.0_dp*peff(l)
                  bval(l) = bx(basis_cntrl + 1)
               END IF
               pmin(l) = zmin(l)
               xv = LOG(pend(l)/pmin(l))/LOG(bval(l)) + 1.e-10_dp
               nval(l) = MAX(CEILING(xv), 0)
               bval(l) = (pend(l)/pmin(l))**(1._dp/nval(l))
            END DO
         ELSE
            DO l = 0, laux
               IF (l <= lval + 1) THEN
                  pend(l) = pmax(l)
                  bval(l) = bv(basis_cntrl + 1)
                  pmin(l) = zmin(l)
               ELSE
                  pend(l) = 4.0_dp*peff(l)
                  bval(l) = bx(basis_cntrl + 1)
               END IF
               xv = LOG(pend(l)/pmin(l))/LOG(bval(l)) + 1.e-10_dp
               nval(l) = MAX(CEILING(xv), 0)
               bval(l) = (pend(l)/pmin(l))**(1._dp/nval(l))
            END DO
         END IF
         !
         lm = MIN(2*lval, 3)
         n1 = MAXVAL(nval(0:lm))
         IF (laux < lm + 1) THEN
            n2 = 0
         ELSE
            n2 = MAXVAL(nval(lm + 1:laux))
         END IF
         !
         nsets = n1 + n2
         ALLOCATE (zet(1, nsets))
         zet = 0.0_dp
         nl = 0
         j = MAXVAL(MAXLOC(nval(0:lm)))
         DO i = 1, n1
            ls1(i) = 0
            ls2(i) = lm
            npgf(i) = 1
            zet(1, i) = pmin(j)*bval(j)**(i - 1)
            DO l = 0, lm
               nl(l, i) = 1
            END DO
         END DO
         j = lm + 1
         DO i = n1 + 1, nsets
            ls1(i) = lm + 1
            ls2(i) = laux
            npgf(i) = 1
            zet(1, i) = pmin(j)*bval(j)**(i - n1 - 1)
            DO l = lm + 1, laux
               nl(l, i) = 1
            END DO
         END DO
         !
         bsname = TRIM(element_symbol)//"-LRI-AUX-"//TRIM(orb_basis_set%name)
         !
         CALL create_aux_basis(lri_aux_basis_set, bsname, nsets, ls1, ls2, nl, npgf, zet)
         !
         DEALLOCATE (zet)
      END IF

   END SUBROUTINE create_lri_aux_basis_set

! **************************************************************************************************
!> \brief ...
!> \param oce_basis ...
!> \param orb_basis ...
!> \param lmax_oce ...
!> \param nbas_oce ...
! **************************************************************************************************
   SUBROUTINE create_oce_basis(oce_basis, orb_basis, lmax_oce, nbas_oce)
      TYPE(gto_basis_set_type), POINTER                  :: oce_basis, orb_basis
      INTEGER, INTENT(IN)                                :: lmax_oce, nbas_oce

      CHARACTER(LEN=default_string_length)               :: bsname
      INTEGER                                            :: i, l, lmax, lx, nset, nx
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: lmin, lset, npgf
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: nl
      INTEGER, DIMENSION(:), POINTER                     :: npgf_orb
      REAL(KIND=dp)                                      :: cval, x, z0, z1
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: zet
      REAL(KIND=dp), DIMENSION(0:9)                      :: zeff, zmax, zmin

      CALL get_basis_keyfigures(orb_basis, lmax, zmin, zmax, zeff)
      IF (nbas_oce < 1) THEN
         CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nset, npgf=npgf_orb)
         nx = SUM(npgf_orb(1:nset))
      ELSE
         nx = 0
      END IF
      nset = MAX(nbas_oce, nx)
      lx = MAX(lmax_oce, lmax)
      !
      bsname = "OCE-"//TRIM(orb_basis%name)
      ALLOCATE (lmin(nset), lset(nset), nl(0:9, nset), npgf(nset), zet(1, nset))
      lmin = 0
      lset = 0
      nl = 1
      npgf = 1
      zet = 0.0_dp
      !
      z0 = MINVAL(zmin(0:lmax))
      z1 = MAXVAL(zmax(0:lmax))
      x = 1.0_dp/REAL(nset - 1, KIND=dp)
      cval = (z1/z0)**x
      zet(1, nset) = z0
      DO i = nset - 1, 1, -1
         zet(1, i) = zet(1, i + 1)*cval
      END DO
      DO i = 1, nset
         x = zet(1, i)
         DO l = 1, lmax
            z1 = 1.05_dp*zmax(l)
            IF (x < z1) lset(i) = l
         END DO
         IF (lset(i) == lmax) lset(i) = lx
      END DO
      !
      CALL create_aux_basis(oce_basis, bsname, nset, lmin, lset, nl, npgf, zet)
      !
      DEALLOCATE (lmin, lset, nl, npgf, zet)

   END SUBROUTINE create_oce_basis
! **************************************************************************************************
!> \brief ...
!> \param basis_set ...
!> \param lmax ...
!> \param zmin ...
!> \param zmax ...
!> \param zeff ...
! **************************************************************************************************
   SUBROUTINE get_basis_keyfigures(basis_set, lmax, zmin, zmax, zeff)
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      INTEGER, INTENT(OUT)                               :: lmax
      REAL(KIND=dp), DIMENSION(0:9), INTENT(OUT)         :: zmin, zmax, zeff

      INTEGER                                            :: i, ipgf, iset, ishell, j, l, nset
      INTEGER, DIMENSION(:), POINTER                     :: lm, npgf, nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: lshell
      REAL(KIND=dp)                                      :: aeff, gcca, gccb, kval, rexp, rint, rno, &
                                                            zeta
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: zet
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: gcc

      CALL get_gto_basis_set(gto_basis_set=basis_set, &
                             nset=nset, &
                             nshell=nshell, &
                             npgf=npgf, &
                             l=lshell, &
                             lmax=lm, &
                             zet=zet, &
                             gcc=gcc)

      lmax = MAXVAL(lm)
      CPASSERT(lmax <= 9)

      zmax = 0.0_dp
      zmin = HUGE(0.0_dp)
      zeff = 0.0_dp

      DO iset = 1, nset
         ! zmin zmax
         DO ipgf = 1, npgf(iset)
            DO ishell = 1, nshell(iset)
               l = lshell(ishell, iset)
               zeta = zet(ipgf, iset)
               zmax(l) = MAX(zmax(l), zeta)
               zmin(l) = MIN(zmin(l), zeta)
            END DO
         END DO
         ! zeff
         DO ishell = 1, nshell(iset)
            l = lshell(ishell, iset)
            kval = fac(l + 1)**2*2._dp**(2*l + 1)/fac(2*l + 2)
            rexp = 0.0_dp
            rno = 0.0_dp
            DO i = 1, npgf(iset)
               gcca = gcc(i, ishell, iset)
               DO j = 1, npgf(iset)
                  zeta = zet(i, iset) + zet(j, iset)
                  gccb = gcc(j, ishell, iset)
                  rint = 0.5_dp*fac(l + 1)/zeta**(l + 2)
                  rexp = rexp + gcca*gccb*rint
                  rint = rootpi*0.5_dp**(l + 2)*dfac(2*l + 1)/zeta**(l + 1.5_dp)
                  rno = rno + gcca*gccb*rint
               END DO
            END DO
            rexp = rexp/rno
            aeff = (fac(l + 1)/dfac(2*l + 1))**2*2._dp**(2*l + 1)/(pi*rexp**2)
            zeff(l) = MAX(zeff(l), aeff)
         END DO
      END DO

   END SUBROUTINE get_basis_keyfigures

! **************************************************************************************************
!> \brief ...
!> \param lmax ...
!> \param zmin ...
!> \param zmax ...
!> \param zeff ...
!> \param pmin ...
!> \param pmax ...
!> \param peff ...
! **************************************************************************************************
   SUBROUTINE get_basis_products(lmax, zmin, zmax, zeff, pmin, pmax, peff)
      INTEGER, INTENT(IN)                                :: lmax
      REAL(KIND=dp), DIMENSION(0:9), INTENT(IN)          :: zmin, zmax, zeff
      REAL(KIND=dp), DIMENSION(0:18), INTENT(OUT)        :: pmin, pmax, peff

      INTEGER                                            :: l1, l2, la

      pmin = HUGE(0.0_dp)
      pmax = 0.0_dp
      peff = 0.0_dp

      DO l1 = 0, lmax
         DO l2 = l1, lmax
            DO la = l2 - l1, l2 + l1
               pmax(la) = MAX(pmax(la), zmax(l1) + zmax(l2))
               pmin(la) = MIN(pmin(la), zmin(l1) + zmin(l2))
               peff(la) = MAX(peff(la), zeff(l1) + zeff(l2))
            END DO
         END DO
      END DO

   END SUBROUTINE get_basis_products
! **************************************************************************************************
!> \brief ...
!> \param lm ...
!> \param npgf ...
!> \param nfun ...
!> \param zet ...
!> \param gcc ...
!> \param nfit ...
!> \param afit ...
!> \param amet ...
!> \param eval ...
! **************************************************************************************************
   SUBROUTINE overlap_maximum(lm, npgf, nfun, zet, gcc, nfit, afit, amet, eval)
      INTEGER, INTENT(IN)                                :: lm, npgf, nfun
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zet
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: gcc
      INTEGER, INTENT(IN)                                :: nfit
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: afit
      REAL(KIND=dp), INTENT(IN)                          :: amet
      REAL(KIND=dp), INTENT(OUT)                         :: eval

      INTEGER                                            :: i, ia, ib, info
      REAL(KIND=dp)                                      :: fij, fxij, intab, p, xij
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: fx, tx, x2, xx

      ! SUM_i(fi M fi)
      fij = 0.0_dp
      DO ia = 1, npgf
         DO ib = 1, npgf
            p = zet(ia) + zet(ib) + amet
            intab = 0.5_dp/p**(lm + 1.5_dp)*gamma1(lm + 1)
            DO i = 1, nfun
               fij = fij + gcc(ia, i)*gcc(ib, i)*intab
            END DO
         END DO
      END DO

      !Integrals (fi M xj)
      ALLOCATE (fx(nfit, nfun), tx(nfit, nfun))
      fx = 0.0_dp
      DO ia = 1, npgf
         DO ib = 1, nfit
            p = zet(ia) + afit(ib) + amet
            intab = 0.5_dp/p**(lm + 1.5_dp)*gamma1(lm + 1)
            DO i = 1, nfun
               fx(ib, i) = fx(ib, i) + gcc(ia, i)*intab
            END DO
         END DO
      END DO

      !Integrals (xi M xj)
      ALLOCATE (xx(nfit, nfit), x2(nfit, nfit))
      DO ia = 1, nfit
         DO ib = 1, nfit
            p = afit(ia) + afit(ib) + amet
            xx(ia, ib) = 0.5_dp/p**(lm + 1.5_dp)*gamma1(lm + 1)
         END DO
      END DO

      !Solve for tab
      tx(1:nfit, 1:nfun) = fx(1:nfit, 1:nfun)
      x2(1:nfit, 1:nfit) = xx(1:nfit, 1:nfit)
      CALL dposv("U", nfit, nfun, x2, nfit, tx, nfit, info)
      IF (info == 0) THEN
         ! value t*xx*t
         xij = 0.0_dp
         DO i = 1, nfun
            xij = xij + DOT_PRODUCT(tx(:, i), MATMUL(xx, tx(:, i)))
         END DO
         ! value t*fx
         fxij = 0.0_dp
         DO i = 1, nfun
            fxij = fxij + DOT_PRODUCT(tx(:, i), fx(:, i))
         END DO
         !
         eval = fij - 2.0_dp*fxij + xij
      ELSE
         ! error in solving for max overlap
         eval = 1.0e10_dp
      END IF

      DEALLOCATE (fx, xx, x2, tx)

   END SUBROUTINE overlap_maximum
! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \param n ...
!> \param eval ...
! **************************************************************************************************
   SUBROUTINE neb_potential(x, n, eval)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: x
      INTEGER, INTENT(IN)                                :: n
      REAL(KIND=dp), INTENT(INOUT)                       :: eval

      INTEGER                                            :: i

      DO i = 2, n
         IF (x(i) < 1.5_dp) THEN
            eval = eval + 10.0_dp*(1.5_dp - x(i))**2
         END IF
      END DO

   END SUBROUTINE neb_potential
! **************************************************************************************************
!> \brief ...
!> \param basis_set ...
!> \param lin ...
!> \param np ...
!> \param nf ...
!> \param zval ...
!> \param gcval ...
! **************************************************************************************************
   SUBROUTINE get_basis_functions(basis_set, lin, np, nf, zval, gcval)
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      INTEGER, INTENT(IN)                                :: lin
      INTEGER, INTENT(OUT)                               :: np, nf
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: zval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: gcval

      INTEGER                                            :: iset, ishell, j1, j2, jf, jp, l, nset
      INTEGER, DIMENSION(:), POINTER                     :: lm, npgf, nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: lshell
      LOGICAL                                            :: toadd
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: zet
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: gcc

      CALL get_gto_basis_set(gto_basis_set=basis_set, &
                             nset=nset, &
                             nshell=nshell, &
                             npgf=npgf, &
                             l=lshell, &
                             lmax=lm, &
                             zet=zet, &
                             gcc=gcc)

      np = 0
      nf = 0
      DO iset = 1, nset
         toadd = .TRUE.
         DO ishell = 1, nshell(iset)
            l = lshell(ishell, iset)
            IF (l == lin) THEN
               nf = nf + 1
               IF (toadd) THEN
                  np = np + npgf(iset)
                  toadd = .FALSE.
               END IF
            END IF
         END DO
      END DO
      ALLOCATE (zval(np), gcval(np, nf))
      zval = 0.0_dp
      gcval = 0.0_dp
      !
      jp = 0
      jf = 0
      DO iset = 1, nset
         toadd = .TRUE.
         DO ishell = 1, nshell(iset)
            l = lshell(ishell, iset)
            IF (l == lin) THEN
               jf = jf + 1
               IF (toadd) THEN
                  j1 = jp + 1
                  j2 = jp + npgf(iset)
                  zval(j1:j2) = zet(1:npgf(iset), iset)
                  jp = jp + npgf(iset)
                  toadd = .FALSE.
               END IF
               gcval(j1:j2, jf) = gcc(1:npgf(iset), ishell, iset)
            END IF
         END DO
      END DO

   END SUBROUTINE get_basis_functions

END MODULE auto_basis
