# Author-side audit of the universal quartic exclusion

This is an obligation audit by the authoring agent, not independent review.

| Obligation | Evidence and limit |
|---|---|
| Every birational gap-21 pair is included | The retained gap normalization and parity identity force b!=0. Complex affine parameter/target scaling puts b1 and preserves all conic/quartic contact conditions. No squarefreeness assumption is made. |
| All quartic coefficients are retained | Leading and constant coefficients are normalized by affine changes of the pencil base; u,v,w initially range over all complex values. |
| Universal identities are actual polynomial identities | The doubled-quartic remainders are monic polynomial reductions. The grading bounds exactly which entry coefficients can contribute to a selected determinant coefficient. Filtering precedes every inhomogeneous specialization. |
| L-nonzero branch is completely excluded | Three successive square factors force w2/3, v4a/3-2d and c2ad-3d2; the next contact coefficient differs from a constant multiple of rho1 by a nonzero constant multiple of L. |
| The other-minor identity is not merely a normal form under a computed basis | After the explicit unit eliminations R=0 and rho1=0, its reconstructed coefficient is literally -40(3w-2)^4/81. All quotient multipliers in the exploratory ideal reduction were zero. |
| Final K/R0 branches include their overlap | The proof splits K0 versus K nonzero, with R0 forced in the latter. Both frozen unit identities also exclude K=R0=0; no division by a potentially zero polynomial is used in the certificate. |
| Each frozen generator represents a necessary equation | Every gap equation and contact coefficient is reconstructed before clearing only rational denominators and rational coefficient content. Those scalars are nonzero constants. E and D are eliminated only by coefficients1/3 and32400. |
| Units are checked independently of the discovery algorithms | Sparse Fraction multiplication replays both multivariate units directly, with mutations rejected. No Gröbner-basis generation, rational-root search, modular inference or proof-trace completion is needed in the checker. |
| Full matrix and weighted implementation agree | The checker verifies the column/remainder grading and compares exact symbolic specializations against the older untruncated matrix construction. These are implementation controls; universality follows from the grading argument and universal identities, not a sample extrapolation. |
| Geometry is not inferred from a green algebraic check | PF2R7QBASE1 is the written bridge forcing the quartic polynomial base. Its connected-fibre, incidence, field-degree and local-order inputs remain explicit. |
| The entire normal-r7 programme is not silently closed | The reciprocal cubic multiple-root case, deriving the entire single-ramified-prime hypothesis, and other ramification patterns remain open. |

The bounded failed cofactor and arithmetic-trace attempts are preserved
in the discovery record. Their printed bases, incomplete traces and
timeouts are not premises of the proof. The complete two-unit route
replaces them rather than claiming a successful replay of those attempts.
