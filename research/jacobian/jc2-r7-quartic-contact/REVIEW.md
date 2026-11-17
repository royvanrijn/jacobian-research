# Author-side obligation audit

This is a self-audit of the new written proof, not independent review.
The source-finality route is preserved separately. No Theorem 3.8
conclusion is used in the present route.

| Obligation | Evidence and boundary |
|---|---|
| Entire ramification divisor really is one e2,f1 prime | Explicit hypothesis; not inferred from a terminal passport. Omitted e1 primes are retained. |
| Earlier three-escape theorem avoids hidden finality | PF2R7INF1 replaces the double-cusp singularity-at-infinity count by Borisov Theorem3.6 and the literal local pullback2E. It still imports the stated boundary and positivity theorems. |
| Degree-six primitive conic really has parameter degree six | Generic fibre maps birationally to its smooth conic, hence adding the pencil parameter generates the entire source field. Primitive Gauss irreducibility identifies its parameter degree with the field degree. |
| Cubic escape exclusions cover special coefficients | PF2R7V1 checks finite and infinite limiting base values, all local degrees1--3 and the allowed zero coefficients. PF2R7G1 only gives a multiple-root restriction; squarefreeness is checked separately for the two exact classes. |
| Remaining infinity points give conics rather than double lines | Two generic points of the index-one infinity curve have distinct target images because its target residue map is an isomorphism. A line cannot have both. |
| E infinity and infinity-curve point have the same base value | Both meet the connected exceptional fibre over A; all exceptional components are vertical by the actual H.F2 budget. No assertion that A is singular is imported from Theorem3.8. |
| Highest binary quadratic depends linearly on the base up to a factor | Its two generic roots describe the graph of the degree-two map on the index-one infinity curve; coprimality of its binary numerator and denominator leaves only the scalar polynomial factor h(z). |
| h order bound covers cancellations | The y^2 coefficient is a nonzero multiple of z*h(z); target-linear terms cannot cancel it. Actual parameter degree six forces deg(h)<=5. |
| Every local degree of the quartic map is covered | k1,2,3 contradict the two exact order equations; k4 implies unique zero at infinity and reciprocal of a quartic polynomial. |
| All quartic polynomial coefficients are included | Constant and leading coefficients are removed by affine base changes. Variables u,v,w range over all complex values; no real or numerical restriction. |
| Contact matrix is the correct necessary condition | In the generic smooth conic chart ramification makes the branch tangent to the conic. Separability converts four incidences plus tangencies to divisibility by the doubled quartic. Pivots use the actual monic degree-six p. |
| Symbolic matrix construction is not blindly trusted | A separate multivariate long division checks every reconstructed remainder. Selected determinants are recomputed from original p,q, then unit coefficients are replayed by dense Fraction multiplication. |
| Every determinant branch is excluded | The u^5 coefficient forces u0; the next product gives two branches per class; four explicit univariate Bezout identities eliminate all four over C. |
| Scope and assurance remain explicit | Only the two nonimmersed classes under the entire single-prime hypothesis are excluded. Immersed quartic/multiple-root cases and JC2 remain open. No formalization or external review. |
