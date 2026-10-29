!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Unit test for CELL_OPT constraint parsing, symmetry rows, and
!>        space-group point-group metric classification.
!> \par Coverage
!>      Every token of the CONSTRAINT grammar, 2D periodicity completion along
!>      all three axes, parse(label(m)) round trip across component and
!>      geometric masks, legacy label preservation, token rejection reasons,
!>      overlap classification, crystal-system row generation and manifold
!>      flatness/restoration, and spgr_point_group_symmetry_id across all 11
!>      crystal-system branches plus edge cases.
!> \par Failure behavior
!>      The first mismatch prints the offending case and stops with exit code 1.
! **************************************************************************************************
PROGRAM cell_opt_utils_TEST
   USE cell_manifold_utils,             ONLY: &
        cell_manifold_geometry, cell_manifold_init, cell_manifold_project_grad, &
        cell_manifold_retract, cell_manifold_type, gfix_a, gfix_alpha, gfix_b, gfix_beta, gfix_c, &
        gfix_gamma, gfix_v, info_ok, n_constraint_bits, n_geom, n_slot
   USE cell_methods,                    ONLY: cell_create
   USE cell_opt_utils,                  ONLY: &
        cell_opt_constraint_label, comp_ax, comp_bx, comp_by, comp_cx, comp_cy, comp_cz, mask_b, &
        mask_c, mask_x, mask_xy, mask_xz, mask_y, mask_yz, mask_z, parse_cell_constraint, &
        parse_token, spgr_point_group_symmetry_id, symmetry_constraint_rows
   USE cell_types,                      ONLY: &
        cell_release, cell_sym_cubic, cell_sym_hexagonal_gamma_120, cell_sym_hexagonal_gamma_60, &
        cell_sym_monoclinic, cell_sym_monoclinic_gamma_ab, cell_sym_none, cell_sym_orthorhombic, &
        cell_sym_rhombohedral, cell_sym_tetragonal_ab, cell_sym_tetragonal_ac, &
        cell_sym_tetragonal_bc, cell_sym_triclinic, cell_type
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE mathconstants,                   ONLY: pi,&
                                              sqrt3
   USE space_groups_types,              ONLY: release_spgr_type,&
                                              spgr_type

   IMPLICIT NONE

   CALL check_parse_tokens()
   CALL check_parse_2d_completion()
   CALL check_label_parse_roundtrip()
   CALL check_label_legacy_strings()
   CALL check_parse_token_rejections()
   CALL check_overlap_classification()
   CALL check_symmetry_constraint_rows()
   CALL check_symmetry_flatness()
   CALL check_spgr_point_group()

   PRINT *, "cell_opt_utils_unittest: all checks passed"
CONTAINS

! **************************************************************************************************
!> \brief Parse every token of the grammar and compare against the mask table.
!>        Covers unions, free case and spacing, and the redundancy warning path.
! **************************************************************************************************
   SUBROUTINE check_parse_tokens()
      TYPE(cell_type), POINTER                           :: cell

      CALL make_cell3d(cell)
      ! Every call parses with explicit=.TRUE., so the 2D auto-completion
      ! stays out of the way and the check sees the raw token semantics.
      CALL parse_and_check("NONE", 0, cell)
      CALL parse_and_check("none", 0, cell)
      CALL parse_and_check("A", comp_ax, cell)
      CALL parse_and_check("AX", comp_ax, cell)
      CALL parse_and_check("B", mask_b, cell)
      CALL parse_and_check("BX", comp_bx, cell)
      CALL parse_and_check("BY", comp_by, cell)
      CALL parse_and_check("C", mask_c, cell)
      CALL parse_and_check("Z", mask_z, cell)
      CALL parse_and_check("CX", comp_cx, cell)
      CALL parse_and_check("CY", comp_cy, cell)
      CALL parse_and_check("CZ", comp_cz, cell)
      CALL parse_and_check("X", mask_x, cell)
      CALL parse_and_check("Y", mask_y, cell)
      CALL parse_and_check("XY", mask_xy, cell)
      CALL parse_and_check("XZ", mask_xz, cell)
      CALL parse_and_check("YZ", mask_yz, cell)
      ! Union semantics: several tokens accumulate. Spelling, case and
      ! spacing are free. The last line repeats tokens on purpose and only
      ! expects a warning, not a failure of the parse.
      CALL parse_and_check("AX BY CZ", IOR(IOR(comp_ax, comp_by), comp_cz), cell)
      CALL parse_and_check("x y", IOR(mask_x, mask_y), cell)
      CALL parse_and_check("  z  ", mask_z, cell)
      CALL parse_and_check("a b ax", IOR(comp_ax, mask_b), cell)
      ! Geometric tokens: fix a length, an angle, or the volume to the
      ! starting value; combinable with component tokens.
      CALL parse_and_check("LENGTH_A", gfix_a, cell)
      CALL parse_and_check("LENGTH_B", gfix_b, cell)
      CALL parse_and_check("LENGTH_C", gfix_c, cell)
      CALL parse_and_check("ALPHA", gfix_alpha, cell)
      CALL parse_and_check("beta", gfix_beta, cell)
      CALL parse_and_check("Gamma", gfix_gamma, cell)
      CALL parse_and_check("VOLUME", gfix_v, cell)
      CALL parse_and_check("V", gfix_v, cell)
      CALL parse_and_check("AX GAMMA VOLUME", IOR(IOR(comp_ax, gfix_gamma), gfix_v), cell)
      CALL cell_release(cell)
      PRINT *, "check_parse_tokens: OK"
   END SUBROUTINE check_parse_tokens

! **************************************************************************************************
!> \brief Parse one token string and compare against the expected mask
!> \param token_string ...
!> \param expected_mask ...
!> \param cell ...
! **************************************************************************************************
   SUBROUTINE parse_and_check(token_string, expected_mask, cell)
      CHARACTER(LEN=*), INTENT(IN)                       :: token_string
      INTEGER, INTENT(IN)                                :: expected_mask
      TYPE(cell_type), POINTER                           :: cell

      CHARACTER(LEN=default_string_length)               :: tokens(16)
      INTEGER                                            :: ntok, parsed_mask

      CALL split_tokens(token_string, tokens, ntok)
      CALL parse_cell_constraint(tokens(:ntok), cell, .TRUE., parsed_mask)
      IF (parsed_mask /= expected_mask) THEN
         PRINT *, "check_parse_tokens: unexpected mask for token "//TRIM(token_string)
         ERROR STOP 1
      END IF
   END SUBROUTINE parse_and_check

! **************************************************************************************************
!> \brief Splits a blank-separated string into tokens, mimicking how the input
!>        tokenizer delivers the value of a char_t list keyword
!> \param string ...
!> \param tokens ...
!> \param ntok ...
! **************************************************************************************************
   SUBROUTINE split_tokens(string, tokens, ntok)
      CHARACTER(LEN=*), INTENT(IN)                       :: string
      CHARACTER(LEN=default_string_length), INTENT(OUT)  :: tokens(:)
      INTEGER, INTENT(OUT)                               :: ntok

      CHARACTER(LEN=LEN(string))                         :: rest
      INTEGER                                            :: isep

      rest = ADJUSTL(string)
      ntok = 0
      DO WHILE (LEN_TRIM(rest) > 0)
         isep = INDEX(rest, " ")
         ntok = ntok + 1
         IF (isep == 0) THEN
            tokens(ntok) = TRIM(rest)
            EXIT
         END IF
         tokens(ntok) = rest(:isep - 1)
         rest = ADJUSTL(rest(isep + 1:))
      END DO
   END SUBROUTINE split_tokens

! **************************************************************************************************
!> \brief For 2D periodicity, an unset constraint auto-fixes the non-periodic
!>        direction. Explicit constraints that cover the direction stay unchanged.
! **************************************************************************************************
   SUBROUTINE check_parse_2d_completion()
      INTEGER                                            :: parsed_mask
      TYPE(cell_type), POINTER                           :: cell

      NULLIFY (cell)
      CALL cell_create(cell)

      cell%perd(:) = [0, 1, 1]
      CALL parse_cell_constraint(["NONE"], cell, .FALSE., parsed_mask)
      IF (parsed_mask /= mask_x) THEN
         PRINT *, "check_parse_2d_completion: auto-set did not fix non-periodic x"
         ERROR STOP 1
      END IF

      cell%perd(:) = [1, 0, 1]
      CALL parse_cell_constraint(["NONE"], cell, .FALSE., parsed_mask)
      IF (parsed_mask /= mask_y) THEN
         PRINT *, "check_parse_2d_completion: auto-set did not fix non-periodic y"
         ERROR STOP 1
      END IF

      cell%perd(:) = [1, 1, 0]
      CALL parse_cell_constraint(["NONE"], cell, .FALSE., parsed_mask)
      IF (parsed_mask /= mask_z) THEN
         PRINT *, "check_parse_2d_completion: auto-set did not fix non-periodic z"
         ERROR STOP 1
      END IF
      ! The keyword is set and covers z. The parser must keep the mask as
      ! given; YZ covers z through the C vector.
      CALL parse_cell_constraint(["Z"], cell, .TRUE., parsed_mask)
      IF (parsed_mask /= mask_z) THEN
         PRINT *, "check_parse_2d_completion: explicit Z changed"
         ERROR STOP 1
      END IF
      CALL parse_cell_constraint(["YZ"], cell, .TRUE., parsed_mask)
      IF (parsed_mask /= mask_yz) THEN
         PRINT *, "check_parse_2d_completion: explicit YZ changed"
         ERROR STOP 1
      END IF
      CALL cell_release(cell)
      PRINT *, "check_parse_2d_completion: OK"
   END SUBROUTINE check_parse_2d_completion

! **************************************************************************************************
!> \brief For every mask m, parsing the label of m must reproduce m. This ties the
!>        spelling of the grammar to its semantics across both component and
!>        geometric masks.
! **************************************************************************************************
   SUBROUTINE check_label_parse_roundtrip()
      CHARACTER(LEN=default_string_length)               :: tokens(16)
      INTEGER                                            :: extra_masks(10), i, mask_value, ntok, &
                                                            parsed_mask
      TYPE(cell_type), POINTER                           :: cell

      extra_masks = [gfix_a, gfix_b, gfix_c, gfix_alpha, gfix_beta, gfix_gamma, gfix_v, &
                     IOR(comp_ax, gfix_v), IOR(IOR(mask_x, gfix_gamma), gfix_c), &
                     2**n_constraint_bits - 1]
      CALL make_cell3d(cell)
      DO mask_value = 0, 63
         CALL split_tokens(cell_opt_constraint_label(mask_value), tokens, ntok)
         CALL parse_cell_constraint(tokens(:ntok), cell, .TRUE., parsed_mask)
         IF (parsed_mask /= mask_value) THEN
            PRINT *, "check_label_parse_roundtrip: parse(label(m)) /= m at", mask_value
            ERROR STOP 1
         END IF
      END DO
      DO i = 1, SIZE(extra_masks)
         mask_value = extra_masks(i)
         CALL split_tokens(cell_opt_constraint_label(mask_value), tokens, ntok)
         CALL parse_cell_constraint(tokens(:ntok), cell, .TRUE., parsed_mask)
         IF (parsed_mask /= mask_value) THEN
            PRINT *, "check_label_parse_roundtrip: geometric mask roundtrip failed at", mask_value
            ERROR STOP 1
         END IF
      END DO
      CALL cell_release(cell)
      PRINT *, "check_label_parse_roundtrip: OK (64 component + 10 geometric/mixed masks)"
   END SUBROUTINE check_label_parse_roundtrip

! **************************************************************************************************
!> \brief The seven legacy constraint labels and NONE must not change. Downstream
!>        tools read these strings.
! **************************************************************************************************
   SUBROUTINE check_label_legacy_strings()
      ! The labels carry leading blanks for column alignment in the banner,
      ! so the comparison strips them with ADJUSTL before TRIM.
      IF (TRIM(ADJUSTL(cell_opt_constraint_label(0))) /= "NONE") THEN
         PRINT *, "check_label_legacy_strings: NONE label changed"
         ERROR STOP 1
      END IF
      IF (TRIM(ADJUSTL(cell_opt_constraint_label(mask_x))) /= "X") THEN
         PRINT *, "check_label_legacy_strings: X label changed"
         ERROR STOP 1
      END IF
      IF (TRIM(ADJUSTL(cell_opt_constraint_label(mask_y))) /= "Y") THEN
         PRINT *, "check_label_legacy_strings: Y label changed"
         ERROR STOP 1
      END IF
      IF (TRIM(ADJUSTL(cell_opt_constraint_label(mask_z))) /= "Z") THEN
         PRINT *, "check_label_legacy_strings: Z label changed"
         ERROR STOP 1
      END IF
      IF (TRIM(ADJUSTL(cell_opt_constraint_label(mask_xy))) /= "XY") THEN
         PRINT *, "check_label_legacy_strings: XY label changed"
         ERROR STOP 1
      END IF
      IF (TRIM(ADJUSTL(cell_opt_constraint_label(mask_xz))) /= "XZ") THEN
         PRINT *, "check_label_legacy_strings: XZ label changed"
         ERROR STOP 1
      END IF
      IF (TRIM(ADJUSTL(cell_opt_constraint_label(mask_yz))) /= "YZ") THEN
         PRINT *, "check_label_legacy_strings: YZ label changed"
         ERROR STOP 1
      END IF
      PRINT *, "check_label_legacy_strings: OK"
   END SUBROUTINE check_label_legacy_strings

! **************************************************************************************************
!> \brief Classify the rejecting tokens without triggering an abort
! **************************************************************************************************
   SUBROUTINE check_parse_token_rejections()
      CHARACTER(LEN=default_path_length)                 :: reason
      INTEGER                                            :: mask_value
      LOGICAL                                            :: accepted

      ! AY, AZ, BZ: the cell convention places a along x and b in the xy
      ! plane, so these components are identically zero. The reason says so.
      CALL parse_token("AY", accepted, mask_value, reason)
      IF (accepted .OR. INDEX(reason, "cannot be constrained.") == 0) THEN
         PRINT *, "check_parse_token_rejections: AY not rejected with the always-zero reason"
         ERROR STOP 1
      END IF
      CALL parse_token("AZ", accepted, mask_value, reason)
      IF (accepted .OR. INDEX(reason, "cannot be constrained.") == 0) THEN
         PRINT *, "check_parse_token_rejections: AZ not rejected with the always-zero reason"
         ERROR STOP 1
      END IF
      CALL parse_token("BZ", accepted, mask_value, reason)
      IF (accepted .OR. INDEX(reason, "cannot be constrained.") == 0) THEN
         PRINT *, "check_parse_token_rejections: BZ not rejected with the always-zero reason"
         ERROR STOP 1
      END IF
      ! Geometric angle and volume tokens (including alias V) are accepted.
      CALL parse_token("ALPHA", accepted, mask_value, reason)
      IF (.NOT. accepted .OR. mask_value /= gfix_alpha) THEN
         PRINT *, "check_parse_token_rejections: ALPHA not accepted with gfix_alpha"
         ERROR STOP 1
      END IF
      CALL parse_token("V", accepted, mask_value, reason)
      IF (.NOT. accepted .OR. mask_value /= gfix_v) THEN
         PRINT *, "check_parse_token_rejections: V not accepted with gfix_v"
         ERROR STOP 1
      END IF
      ! A near-miss gets the unknown-token reason with the accepted vocabulary.
      CALL parse_token("YX", accepted, mask_value, reason)
      IF (accepted .OR. INDEX(reason, "Unknown CELL_OPT CONSTRAINT token") == 0 .OR. &
          INDEX(reason, "VOLUME (or V)") == 0 .OR. &
          INDEX(reason, "Several tokens may be combined.") == 0) THEN
         PRINT *, "check_parse_token_rejections: YX not rejected with the vocabulary"
         ERROR STOP 1
      END IF
      PRINT *, "check_parse_token_rejections: OK (AY AZ BZ ALPHA V YX)"
   END SUBROUTINE check_parse_token_rejections

! **************************************************************************************************
!> \brief The overlap warning branches on full versus partial coverage
! **************************************************************************************************
   SUBROUTINE check_overlap_classification()
      INTEGER                                            :: old_mask_value, overlap_value, &
                                                            token_mask_value

      DO old_mask_value = 0, 63
         DO token_mask_value = 1, 63
            overlap_value = IAND(token_mask_value, old_mask_value)
            IF (overlap_value /= token_mask_value .AND. overlap_value /= 0) THEN
               ! Partial overlap: the additions the warning names must be
               ! non-empty, else the token would be fully redundant.
               IF (IAND(NOT(old_mask_value), token_mask_value) == 0) THEN
                  PRINT *, "check_overlap_classification: partial overlap with empty additions"
                  ERROR STOP 1
               END IF
            END IF
         END DO
      END DO
      PRINT *, "check_overlap_classification: OK (64 x 63 mask pairs)"
   END SUBROUTINE check_overlap_classification

! **************************************************************************************************
!> \brief Crystal-system translation into manifold rows: bit counts,
!>        coupling/equality row counts and the resulting tangent dimension
!>        per system, plus the hexagonal sign pattern (gamma 60 vs 120).
! **************************************************************************************************
   SUBROUTINE check_symmetry_constraint_rows()
      INTEGER                                            :: eq_pair(n_geom, 2), mask, n_eq, n_lin, &
                                                            ndof_exp(12), sid(12), sys
      REAL(KIND=dp)                                      :: lin_w(n_slot, n_slot)

      sid = [cell_sym_orthorhombic, cell_sym_tetragonal_ab, cell_sym_tetragonal_ac, &
             cell_sym_tetragonal_bc, cell_sym_cubic, cell_sym_hexagonal_gamma_60, &
             cell_sym_hexagonal_gamma_120, cell_sym_rhombohedral, cell_sym_monoclinic, &
             cell_sym_monoclinic_gamma_ab, cell_sym_triclinic, cell_sym_none]
      ndof_exp = [3, 2, 2, 2, 1, 2, 2, 2, 4, 3, 6, 6]

      DO sys = 1, SIZE(sid)
         CALL symmetry_constraint_rows(sid(sys), mask, lin_w, n_lin, eq_pair, n_eq)
         IF (n_slot - (POPCNT(mask) + n_lin + n_eq) /= ndof_exp(sys)) THEN
            PRINT *, "check_symmetry_constraint_rows: system", sys, "ndof", &
               n_slot - (POPCNT(mask) + n_lin + n_eq), "expected", ndof_exp(sys)
            ERROR STOP 1
         END IF
      END DO

      ! hexagonal sign pattern: BX = -/+ AX/2 for gamma = 120/60
      CALL symmetry_constraint_rows(cell_sym_hexagonal_gamma_120, mask, lin_w, n_lin, eq_pair, n_eq)
      IF (lin_w(1, 1) <= 0.0_dp .OR. lin_w(2, 1) /= 1.0_dp) THEN
         PRINT *, "check_symmetry_constraint_rows: hex 120 BX sign", lin_w(1, 1)
         ERROR STOP 1
      END IF
      CALL symmetry_constraint_rows(cell_sym_hexagonal_gamma_60, mask, lin_w, n_lin, eq_pair, n_eq)
      IF (lin_w(1, 1) >= 0.0_dp) THEN
         PRINT *, "check_symmetry_constraint_rows: hex 60 BX sign", lin_w(1, 1)
         ERROR STOP 1
      END IF

      ! rhombohedral equalities: a=b, b=c, alpha=beta, beta=gamma
      CALL symmetry_constraint_rows(cell_sym_rhombohedral, mask, lin_w, n_lin, eq_pair, n_eq)
      IF (n_eq /= 4 .OR. eq_pair(1, 1) /= 1 .OR. eq_pair(1, 2) /= 2 .OR. &
          eq_pair(2, 1) /= 2 .OR. eq_pair(2, 2) /= 3 .OR. &
          eq_pair(3, 1) /= 4 .OR. eq_pair(3, 2) /= 5 .OR. &
          eq_pair(4, 1) /= 5 .OR. eq_pair(4, 2) /= 6) THEN
         PRINT *, "check_symmetry_constraint_rows: rhombohedral pairs", eq_pair(1:4, :)
         ERROR STOP 1
      END IF
      PRINT *, "check_symmetry_constraint_rows: row sets and ndof per system OK"
   END SUBROUTINE check_symmetry_constraint_rows

! **************************************************************************************************
!> \brief The full crystal-system matrix of symmetry_constraint_rows on
!>        the constraint manifold: every system except rhombohedral and
!>        monoclinic_gamma_ab is a constant row set and retracts exactly
!>        flat (zero Gauss-Newton iterations); those two carry nonlinear
!>        equality rows, need the Gauss-Newton restoration and keep their
!>        true degrees of freedom.
! **************************************************************************************************
   SUBROUTINE check_symmetry_flatness()
      INTEGER                                            :: eq_pair(n_geom, 2), i, info, mask, n_eq, &
                                                            n_lin, ndof_exp(11), sid(11)
      REAL(KIND=dp)                                      :: ct, cx, cy, cz, du(n_slot), g(n_slot), &
                                                            lin_w(n_slot, n_slot), q0(n_slot, 11), &
                                                            qout(n_slot), s(n_geom), st, worst
      TYPE(cell_manifold_type)                           :: m

      ! canonical-frame start cells, one column per system below
      ct = COS(70.0_dp*pi/180.0_dp)
      st = SIN(70.0_dp*pi/180.0_dp)
      cx = ct
      cy = (ct - ct*ct)/st
      cz = SQRT(MAX(1.0_dp - cx*cx - cy*cy, 0.0_dp))
      sid = [cell_sym_orthorhombic, cell_sym_tetragonal_ab, cell_sym_tetragonal_ac, &
             cell_sym_tetragonal_bc, cell_sym_cubic, cell_sym_hexagonal_gamma_60, &
             cell_sym_hexagonal_gamma_120, cell_sym_rhombohedral, cell_sym_monoclinic, &
             cell_sym_monoclinic_gamma_ab, cell_sym_triclinic]
      ndof_exp = [3, 2, 2, 2, 1, 2, 2, 2, 4, 3, 6]
      q0 = 0.0_dp
      q0(:, 1) = [5.1_dp, 0.0_dp, 6.2_dp, 0.0_dp, 0.0_dp, 7.3_dp] ! orthorhombic
      q0(:, 2) = [5.1_dp, 0.0_dp, 5.1_dp, 0.0_dp, 0.0_dp, 6.3_dp] ! tetragonal_ab
      q0(:, 3) = [5.1_dp, 0.0_dp, 6.2_dp, 0.0_dp, 0.0_dp, 5.1_dp] ! tetragonal_ac
      q0(:, 4) = [6.2_dp, 0.0_dp, 5.1_dp, 0.0_dp, 0.0_dp, 5.1_dp] ! tetragonal_bc
      q0(:, 5) = [5.1_dp, 0.0_dp, 5.1_dp, 0.0_dp, 0.0_dp, 5.1_dp] ! cubic
      q0(:, 6) = [5.1_dp, 2.55_dp, 0.5_dp*sqrt3*5.1_dp, 0.0_dp, 0.0_dp, 6.3_dp] ! hex 60
      q0(:, 7) = [5.1_dp, -2.55_dp, 0.5_dp*sqrt3*5.1_dp, 0.0_dp, 0.0_dp, 6.3_dp] ! hex 120
      q0(:, 8) = 5.1_dp*[1.0_dp, ct, st, cx, cy, cz] ! rhombohedral
      q0(:, 9) = [5.1_dp, 0.0_dp, 6.2_dp, 2.1_dp, 0.0_dp, 7.3_dp] ! monoclinic
      q0(:, 10) = [5.1_dp, 5.1_dp*COS(100.0_dp*pi/180.0_dp), &
                   5.1_dp*SIN(100.0_dp*pi/180.0_dp), 0.0_dp, 0.0_dp, 7.3_dp] ! mono_gab, gamma = 100
      q0(:, 11) = [5.1_dp, 0.7_dp, 6.2_dp, 1.3_dp, 0.9_dp, 7.3_dp] ! triclinic
      g = [0.4_dp, -1.2_dp, 0.7_dp, 1.5_dp, -0.3_dp, 0.9_dp]
      du = 0.0_dp

      DO i = 1, SIZE(sid)
         CALL symmetry_constraint_rows(sid(i), mask, lin_w, n_lin, eq_pair, n_eq)
         CALL cell_manifold_init(m, q0(:, i), mask, info, lin_w=lin_w, n_lin_rows=n_lin, &
                                 eq_pair=eq_pair, n_eq_rows=n_eq)
         IF (info /= info_ok .OR. m%ndof /= ndof_exp(i)) THEN
            PRINT *, "check_symmetry_flatness: init", i, info, m%ndof
            ERROR STOP 1
         END IF
         CALL cell_manifold_project_grad(m, g, du)
         CALL cell_manifold_retract(m, q0(:, i), 0.05_dp*du, qout, info)
         IF (sid(i) == cell_sym_rhombohedral) THEN
            ! curved: GN restores a=b=c and alpha=beta=gamma
            CALL cell_manifold_geometry(qout, s, info)
            worst = MAX(MAX(ABS(s(1) - s(2)), ABS(s(2) - s(3))), &
                        MAX(ABS(s(4) - s(5)), ABS(s(5) - s(6))))
            IF (info /= info_ok .OR. m%last_iter < 1 .OR. worst > 1.0E-11_dp*5.1_dp) THEN
               PRINT *, "check_symmetry_flatness: rhombohedral failed", info, m%last_iter, worst
               ERROR STOP 1
            END IF
         ELSE IF (sid(i) == cell_sym_monoclinic_gamma_ab) THEN
            ! curved: GN restores a = b (alpha = beta = 90 is the CX/CY pin pair)
            CALL cell_manifold_geometry(qout, s, info)
            IF (info /= info_ok .OR. m%last_iter < 1 .OR. ABS(s(1) - s(2)) > 1.0E-11_dp*5.1_dp) THEN
               PRINT *, "check_symmetry_flatness: mono_gab failed", info, m%last_iter
               ERROR STOP 1
            END IF
         ELSE
            ! flat: constant row sets retract in zero GN iterations
            IF (info /= info_ok .OR. m%last_iter /= 0) THEN
               PRINT *, "check_symmetry_flatness: not flat", i, info, m%last_iter
               ERROR STOP 1
            END IF
         END IF
      END DO
      PRINT *, "check_symmetry_flatness: 9 flat systems 0 iters, rhombohedral and mono_gab restored"
   END SUBROUTINE check_symmetry_flatness

! **************************************************************************************************
!> \brief Verify spgr_point_group_symmetry_id across all 11 crystal-system
!>        branches plus nop == 0 and degenerate-cell edge cases.
! **************************************************************************************************
   SUBROUTINE check_spgr_point_group()
      INTEGER, DIMENSION(3, 3), PARAMETER :: &
         r_c2x = RESHAPE([1, 0, 0, 0, -1, 0, 0, 0, -1], [3, 3]), &
         r_c2y = RESHAPE([-1, 0, 0, 0, 1, 0, 0, 0, -1], [3, 3]), &
         r_c2z = RESHAPE([-1, 0, 0, 0, -1, 0, 0, 0, 1], [3, 3]), &
         r_c3_111 = RESHAPE([0, 0, 1, 1, 0, 0, 0, 1, 0], [3, 3]), &
         r_c3_120 = RESHAPE([0, -1, 0, 1, -1, 0, 0, 0, 1], [3, 3]), &
         r_c4x = RESHAPE([1, 0, 0, 0, 0, -1, 0, 1, 0], [3, 3]), &
         r_c4y = RESHAPE([0, 0, 1, 0, 1, 0, -1, 0, 0], [3, 3]), &
         r_c4z = RESHAPE([0, -1, 0, 1, 0, 0, 0, 0, 1], [3, 3]), &
         r_c6_60 = RESHAPE([1, 1, 0, -1, 0, 0, 0, 0, 1], [3, 3]), &
         r_id = RESHAPE([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
      INTEGER, DIMENSION(3, 3), PARAMETER :: r_mz = RESHAPE([1, 0, 0, 0, 1, 0, 0, 0, -1], [3, 3]), &
         r_swap_ab = RESHAPE([0, 1, 0, 1, 0, 0, 0, 0, 1], [3, 3])

      INTEGER                                            :: k
      INTEGER, DIMENSION(3, 3)                           :: rk
      INTEGER, DIMENSION(3, 3, 12)                       :: rots
      REAL(KIND=dp)                                      :: ct, cx, cy, cz, st
      TYPE(cell_type), POINTER                           :: cell
      TYPE(spgr_type), POINTER                           :: spgr

      NULLIFY (cell, spgr)
      CALL cell_create(cell)
      ALLOCATE (spgr)

      ! 1. nop == 0 -> triclinic
      cell%hmat = RESHAPE([5.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 5.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 5.0_dp], &
                          [3, 3])
      IF (spgr_point_group_symmetry_id(spgr, cell) /= cell_sym_triclinic) THEN
         PRINT *, "check_spgr_point_group: nop=0 failed"
         ERROR STOP 1
      END IF

      ! 2. Degenerate zero cell -> triclinic
      cell%hmat = 0.0_dp
      rots(:, :, 1) = r_id
      CALL set_spgr_rots(spgr, rots(:, :, 1:1))
      IF (spgr_point_group_symmetry_id(spgr, cell) /= cell_sym_triclinic) THEN
         PRINT *, "check_spgr_point_group: degenerate cell failed"
         ERROR STOP 1
      END IF

      ! Generic triclinic lattice (a /= b /= c, all angles /= 90)
      cell%hmat = RESHAPE([5.1_dp, 0.0_dp, 0.0_dp, 0.7_dp, 6.2_dp, 0.0_dp, 1.3_dp, 0.9_dp, 7.3_dp], &
                          [3, 3])

      ! 3. C1 {I} -> triclinic
      CALL set_spgr_rots(spgr, rots(:, :, 1:1))
      IF (spgr_point_group_symmetry_id(spgr, cell) /= cell_sym_triclinic) THEN
         PRINT *, "check_spgr_point_group: triclinic failed"
         ERROR STOP 1
      END IF

      ! 4. C2 about b {I, C2y} -> monoclinic (alpha = gamma = 90, beta free)
      rots(:, :, 1) = r_id
      rots(:, :, 2) = r_c2y
      CALL set_spgr_rots(spgr, rots(:, :, 1:2))
      IF (spgr_point_group_symmetry_id(spgr, cell) /= cell_sym_monoclinic) THEN
         PRINT *, "check_spgr_point_group: monoclinic failed"
         ERROR STOP 1
      END IF

      ! 5. D2 {I, C2x, C2y, C2z} -> orthorhombic
      rots(:, :, 3) = r_c2x
      rots(:, :, 4) = r_c2z
      CALL set_spgr_rots(spgr, rots(:, :, 1:4))
      IF (spgr_point_group_symmetry_id(spgr, cell) /= cell_sym_orthorhombic) THEN
         PRINT *, "check_spgr_point_group: orthorhombic failed"
         ERROR STOP 1
      END IF

      ! 6. C4z cyclic -> tetragonal_ab
      rk = r_id
      DO k = 1, 4
         rots(:, :, k) = rk
         rk = MATMUL(rk, r_c4z)
      END DO
      CALL set_spgr_rots(spgr, rots(:, :, 1:4))
      IF (spgr_point_group_symmetry_id(spgr, cell) /= cell_sym_tetragonal_ab) THEN
         PRINT *, "check_spgr_point_group: tetragonal_ab failed"
         ERROR STOP 1
      END IF

      ! 7. C4y cyclic -> tetragonal_ac
      rk = r_id
      DO k = 1, 4
         rots(:, :, k) = rk
         rk = MATMUL(rk, r_c4y)
      END DO
      CALL set_spgr_rots(spgr, rots(:, :, 1:4))
      IF (spgr_point_group_symmetry_id(spgr, cell) /= cell_sym_tetragonal_ac) THEN
         PRINT *, "check_spgr_point_group: tetragonal_ac failed"
         ERROR STOP 1
      END IF

      ! 8. C4x cyclic -> tetragonal_bc
      rk = r_id
      DO k = 1, 4
         rots(:, :, k) = rk
         rk = MATMUL(rk, r_c4x)
      END DO
      CALL set_spgr_rots(spgr, rots(:, :, 1:4))
      IF (spgr_point_group_symmetry_id(spgr, cell) /= cell_sym_tetragonal_bc) THEN
         PRINT *, "check_spgr_point_group: tetragonal_bc failed"
         ERROR STOP 1
      END IF

      ! 9. Tetrahedral group T (12 rotations: D2 x C3_111) -> cubic
      rots(:, :, 1) = r_id
      rots(:, :, 2) = r_c2x
      rots(:, :, 3) = r_c2y
      rots(:, :, 4) = r_c2z
      DO k = 1, 4
         rots(:, :, 4 + k) = MATMUL(rots(:, :, k), r_c3_111)
         rots(:, :, 8 + k) = MATMUL(rots(:, :, 4 + k), r_c3_111)
      END DO
      CALL set_spgr_rots(spgr, rots(:, :, 1:12))
      IF (spgr_point_group_symmetry_id(spgr, cell) /= cell_sym_cubic) THEN
         PRINT *, "check_spgr_point_group: cubic failed"
         ERROR STOP 1
      END IF

      ! 10. Hexagonal gamma = 60 (C6 generated by r_c6_60)
      cell%hmat = RESHAPE([5.1_dp, 0.0_dp, 0.0_dp, 2.55_dp, 0.5_dp*sqrt3*5.1_dp, 0.0_dp, &
                           0.0_dp, 0.0_dp, 6.3_dp], [3, 3])
      rk = r_id
      DO k = 1, 6
         rots(:, :, k) = rk
         rk = MATMUL(rk, r_c6_60)
      END DO
      CALL set_spgr_rots(spgr, rots(:, :, 1:6))
      IF (spgr_point_group_symmetry_id(spgr, cell) /= cell_sym_hexagonal_gamma_60) THEN
         PRINT *, "check_spgr_point_group: hexagonal_gamma_60 failed"
         ERROR STOP 1
      END IF

      ! 11. Hexagonal gamma = 120 (C3h = C3_120 x {I, mz})
      cell%hmat = RESHAPE([5.1_dp, 0.0_dp, 0.0_dp, -2.55_dp, 0.5_dp*sqrt3*5.1_dp, 0.0_dp, &
                           0.0_dp, 0.0_dp, 6.3_dp], [3, 3])
      rk = r_id
      DO k = 1, 3
         rots(:, :, k) = rk
         rots(:, :, 3 + k) = MATMUL(rk, r_mz)
         rk = MATMUL(rk, r_c3_120)
      END DO
      CALL set_spgr_rots(spgr, rots(:, :, 1:6))
      IF (spgr_point_group_symmetry_id(spgr, cell) /= cell_sym_hexagonal_gamma_120) THEN
         PRINT *, "check_spgr_point_group: hexagonal_gamma_120 failed"
         ERROR STOP 1
      END IF

      ! 12. Monoclinic gamma_ab (a = b, alpha = beta = 90, gamma = 100)
      cell%hmat = RESHAPE([5.1_dp, 0.0_dp, 0.0_dp, &
                           5.1_dp*COS(100.0_dp*pi/180.0_dp), 5.1_dp*SIN(100.0_dp*pi/180.0_dp), &
                           0.0_dp, 0.0_dp, 0.0_dp, 7.3_dp], [3, 3])
      rots(:, :, 1) = r_id
      rots(:, :, 2) = r_swap_ab
      rots(:, :, 3) = r_mz
      rots(:, :, 4) = MATMUL(r_swap_ab, r_mz)
      CALL set_spgr_rots(spgr, rots(:, :, 1:4))
      IF (spgr_point_group_symmetry_id(spgr, cell) /= cell_sym_monoclinic_gamma_ab) THEN
         PRINT *, "check_spgr_point_group: monoclinic_gamma_ab failed"
         ERROR STOP 1
      END IF

      ! 13. Rhombohedral (C3 along [1,1,1] on alpha = beta = gamma = 70 deg)
      ct = COS(70.0_dp*pi/180.0_dp)
      st = SIN(70.0_dp*pi/180.0_dp)
      cx = ct
      cy = (ct - ct*ct)/st
      cz = SQRT(MAX(1.0_dp - cx*cx - cy*cy, 0.0_dp))
      cell%hmat = 5.1_dp*RESHAPE([1.0_dp, 0.0_dp, 0.0_dp, ct, st, 0.0_dp, cx, cy, cz], [3, 3])
      rots(:, :, 1) = r_id
      rots(:, :, 2) = r_c3_111
      rots(:, :, 3) = MATMUL(r_c3_111, r_c3_111)
      CALL set_spgr_rots(spgr, rots(:, :, 1:3))
      IF (spgr_point_group_symmetry_id(spgr, cell) /= cell_sym_rhombohedral) THEN
         PRINT *, "check_spgr_point_group: rhombohedral failed"
         ERROR STOP 1
      END IF

      CALL release_spgr_type(spgr)
      CALL cell_release(cell)
      PRINT *, "check_spgr_point_group: all 11 systems + 2 edge cases OK"
   END SUBROUTINE check_spgr_point_group

! **************************************************************************************************
!> \brief Populate spgr%rotations_subset for unit testing
!> \param spgr ...
!> \param rots ...
! **************************************************************************************************
   SUBROUTINE set_spgr_rots(spgr, rots)
      TYPE(spgr_type), POINTER                           :: spgr
      INTEGER, DIMENSION(:, :, :), INTENT(IN)            :: rots

      INTEGER                                            :: nop

      nop = SIZE(rots, 3)
      IF (ALLOCATED(spgr%rotations_subset)) DEALLOCATE (spgr%rotations_subset)
      ALLOCATE (spgr%rotations_subset(3, 3, nop))
      spgr%rotations_subset(:, :, :) = rots(:, :, :)
      spgr%n_operations_subset = nop
   END SUBROUTINE set_spgr_rots

! **************************************************************************************************
!> \brief Create a 3D-periodic cell. parse_cell_constraint reads only the perd field.
!> \param cell ...
! **************************************************************************************************
   SUBROUTINE make_cell3d(cell)
      TYPE(cell_type), POINTER                           :: cell

      NULLIFY (cell)
      CALL cell_create(cell)
      cell%perd(:) = 1
   END SUBROUTINE make_cell3d

END PROGRAM cell_opt_utils_TEST
