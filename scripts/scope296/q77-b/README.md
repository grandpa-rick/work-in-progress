# Q77 (4,6) bicommutant — scope296 compute (2026-11-14). Grade: computed.

Files: fastBG.py (fastB + generator filter, deterministic chunked RREF; slow), fastBR.py (randomised
compression of the Hom constraints, multiprocessing over D; used for (4,6)), fastA2.py (copy, P via env).
split.py copied from prove295.

Method. B_G = commutant of G = {1_mu} ∪ {S-maps phi_{mu,w}: nu-mu = ±(e_i - e_{i+1})} inside
⊕_{D,D'} ⊕_mu Hom(M_mu x_D, M_mu x_D') (H_N(0) basic, so V = ⊕_D V e_D, V e_D ≅ V x_D, Norton).
G ⊆ S ⇒ B ⊆ B_G (no generation claim needed for the A=B direction).
fastBR: per (D,D') takes K = nv+8 random combinations Σ_g u_g^T(X_nu M1_g − M2_g X_mu) v_g; the
rows lie in the true row space, so the computed rank ≤ the true rank ⇒ computed nullity ≥ the true one (safe direction).

Chain: dim_p A (Krylov) ≤ dim_Q A ≤ dim_Q B ≤ dim_Q B_G ≤ dim_p B_G (exact) ≤ dim_p B_G (randomised).
So dim_p B_G(randomised) = dim_p A forces dim_Q A = dim_Q B.

Validation (gen and all-maps agree): (2,5) 43, (2,6) 142, (2,7) 493, (3,4) 23, (3,5) 103, (3,6) 515, (4,5) 119.

## (4,6) result
dim A = 694 (p = 2^31−1 in PROVE 295, p = 2147483629 here: runA46_p2.log, 169 s, 1.2 GB).
dim S = 54264 = C(21,6) (from the PROVE 295 count). G uses 4380 maps.
dim B_G = 694 at both p = 2^31−1 (522 s) and p = 2147483629 (689 s), 4 workers, peak about 1.06 GB per worker.
By the chain above, A = B at (4,6) over Q. Grade: computed.
(4,7) NOT run: largest Hom block has 17008 unknowns (4.95M in total), and fastA2 cannot handle d = 16384.
