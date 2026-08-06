# A boundary-selected classical-symbol census for `DC_2`

> **Status and scope.** This note starts the requested branch away from the
> normalized degree-five symbol. It constructs four exact classical
> incidence packets and attaches the same first bounded relative
> fiber/connection complex to each coefficient base. The resulting
> obstruction and strong-cocycle modules vanish on this canonical
> `4 -> 3 -> 2` truncation. This is a parameter-uniform local linear result,
> not a polynomial Weyl quantization, a boundary-descent theorem, or a result
> about `DC_2`.

The point of the census is to freeze inequivalent classical inputs before
increasing quantum correction support. A family enters the later search only
after its coefficient base, inverse incidence, boundary lattice, and
rank-two admission status have been recorded.

## 1. The four packets

| key | classical base and incidence | selected boundary | degree | rank-two status |
|---|---|---|---:|---|
| `W6` | normalized degree-six weighted seeds \(H=\sum_{j=3}^6h_j(W^j-W^2)\), \(\sum(j-2)h_j=-1\) | weighted \(\mathbb A^1\) normalization, controlled exponent one | 6 | polynomial exact-symplectic completion is proved |
| `R21` | cancellation type \((m,r)=(2,1)\), \(q^2-4q+6=0\) | reciprocal \(\mathbb G_m\) normalization, controlled exponent one | 4 | localized canonical rank-two candidate; polynomial four-dimensional descent is open |
| `C3` | \(H=W^2(1-W)\), selected from \(\mathbb Q[u^2,u^3]\subset\mathbb Q[u]\) | cusp conductor \(u^2\mathbb Q[u]\), quotient length one | 3 | the polynomial `c=-9` exact-symplectic completion is proved |
| `I3,6` | the tower \(F_6\circ F_3\), with the outer degree-six seed varying | imprimitive two-stage incidence, six blocks of size three | 18 | compose the two polynomial exact-symplectic completions |

These are not four specializations of the old quintic. Their displayed
fingerprints

\[
(\text{generic degree},\text{boundary type},\text{unit rank},
 \text{conductor length},\text{block system})
\]

are pairwise distinct. The degrees already separate all four selected
packets; the remaining entries retain the boundary information needed when
the census is enlarged to same-degree comparisons.

The reciprocal row is intentionally marked as a candidate. The
three-dimensional cancellation map is polynomial and Keller, but the
repository does not yet prove a polynomial rank-two completion for its
quartic reciprocal incidence. Treating its formal-local Darboux chart as a
global polynomial symbol would overstate the current result.

## 2. Exact defining data

For `W6`, eliminate \(h_3\):

\[
h_3=-1-2h_4-3h_5-4h_6.
\]

Then \(H(0)=H'(0)=H(1)=0\) and \(H'(1)=-1\) identically over

\[
B_W=\mathbb Q[h_4,h_5,h_6,(h_6(H''(1)+2))^{-1}].
\]

The rank-two descent theorem applies to this whole admissible family.

For `R21`, the unique cancellation jet over
\(B_R=\mathbb Q[q]/(q^2-4q+6)\) is

\[
h(A)=q+(4q-6)A.
\]

The inverse incidence is

\[
\Psi(T)=T-\frac{Q^2T^2}{2}
 \frac{2PQT^3}{3}-\frac{P^2T^4}{4}-R,
\]

and

\[
\partial_T\Psi=1-T(Q-PT)^2.
\]

Thus this is the exact degree-four reciprocal packet, not a heuristic
Laurent model.

For `C3`, translating the cubic discriminant by

\[
s=\frac13+S,\qquad t=\frac1{27}+\frac S3+V
\]

gives

\[
\operatorname{disc}_W(H(W)-sW+t)=-4S^3-27V^2.
\]

The normalization \(S=-3u^2,\ V=-2u^3\) has conductor
\(u^2\mathbb Q[u]\). This conductor selection recovers the known cubic
weighted symbol; it is not advertised as a new construction mechanism.

For `I3,6`, the strict function-field tower has degrees \(3\) and \(6\).
Its generic degree is \(18\), and the intermediate field gives six blocks
of size three. Composition of the two four-dimensional exact-symplectic
maps is again exact symplectic. This is the imprimitive packet; no expansion
of the large composite coordinates is needed to define it.

Canonical sources for these assertions are:

- [rank-two weighted descent](RANK_TWO_SYMPLECTIC_DESCENT.md);
- [the cancellation construction](../cancellation/CONSTRUCTION.md);
- [the conductor-first cusp realization](CONDUCTOR_FIRST_FOUNDATIONAL_CUSP_KELLER.md);
- [imprimitive Keller factorization](../verified/IMPRIMITIVE_KELLER_FACTORIZATION.md); and
- [the `c=-9` rank-two completion](QUADRATIC_LADDER_AND_POISSON_AUDIT.md).

## 3. Relative restricted complex before correction searches

Over each base \(B_\bullet\), use canonical fiber coordinates \(S,T\) and
central coordinate \(R\). At one homogeneous stage, a correction is
\((s,t,a)\), and

\[
d_1(s,t,a)=
\left(s_S+t_T,\ s_R-a_T,\ t_R+a_S\right).
\]

Hamiltonian gauge is

\[
d_0(h)=(-h_T,h_S,-h_R).
\]

The defects obey

\[
\partial_RF-\partial_SG-\partial_TH=0.
\]

Restrict gauges, corrections, and defects to total degrees \(4,3,2\).
Using the closed-defect module rather than three independent defect blocks
gives the integral matrix complex

\[
35\xrightarrow[\operatorname{rank}34]{d_0}
60\xrightarrow[\operatorname{rank}26]{d_1}26.
\]

It has

\[
(H^0,H^1,H^2)=(B_\bullet,0,0)
\]

after base change to every characteristic-zero census algebra. The proof is
stronger than a generic-rank calculation: a \(26\times26\) minor of \(d_1\)
is a nonzero rational unit. Hence \(d_1\) is split surjective over every
\(B_\bullet\), including the nonnormal cusp conductor algebra.

Let

\[
E^{(2)}_\bullet=\operatorname{coker}d_1.
\]

The coherent Fitting calculation is therefore

\[
E^{(2)}_\bullet=0,\qquad
\operatorname{Fitt}_0(E^{(2)}_\bullet)=(1).
\]

At this first bounded stage there are no earlier quantum lift parameters.
Consequently the strong dual-cocycle module is

\[
\mathcal P^{(2)}_\bullet
=\ker(d_1^\vee)=0.
\]

This says only that a closed defect of canonical degree at most two has no
linear obstruction in this truncation. Pullback to the Ore coordinates can
raise degrees and introduce poles, so this result does not imply that the
first Moyal defect of any displayed global symbol lies in the tested
module.

## 4. Branch order

The next computation is family-specific and should preserve this order:

1. write the exact Ore-to-incidence bridge for the packet;
2. intersect the canonical correction module with its polynomial boundary
   lattice;
3. compute \(E_m=\operatorname{coker}d_{1,m}\), its Fitting strata, and the
   coherent strong module before evaluating a large obstruction section;
4. only then enlarge correction support on a surviving horizontal
   component; and
5. after local lifting, prove support saturation and conductor descent.

For `R21`, step 1 also includes the missing polynomial rank-two admission
test. For `C3`, the conductor square must be kept even though the canonical
linear obstruction module is zero. For `I3,6`, corrections must respect the
intermediate block algebra; treating the degree-18 composite as a primitive
single incidence would discard its defining feature.

The first success criterion was reached on `W6` by the
[degree-six relative quantization theorem](DEGREE_SIX_RELATIVE_QUANTIZATION_OBSTRUCTION.md).
The subsequent
[classical-symbol family search](DC2_CLASSICAL_SYMBOL_FAMILY_SEARCH.md)
computes its exceptional strong-cocycle locus. A residue-degree-four
component is reconstructed over \(\mathbb Q\); its genuine order-five lift
scheme is \(\mathbb A^4\) over the quartic residue field, but its complete
inherited order-seven ideal is the unit ideal. Thus this component is closed
in the parity-preserving, root-weight-homogeneous filtration. Other weights,
odd corrections, wider filtrations, and genuinely different noninjective
classical families stay open.

## Reproduction

Run:

```bash
.venv/bin/python scripts/verify_dc2_classical_symbol_census.py \
  --output artifacts/generated-results/dc2_classical_symbol_census.json
```

The command verifies the four defining packets, their exact fingerprints,
the relative complex, one unit maximal minor, the unit Fitting ideal, the
zero strong-cocycle module, and ranks at three good primes.
