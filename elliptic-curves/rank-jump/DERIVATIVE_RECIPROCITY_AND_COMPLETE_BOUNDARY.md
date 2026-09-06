# A canonical derivative class closes the Selmer boundary

For each of the three completely covered fibres, the image of the **full**
2-Selmer group at the bad places is now exactly the image already generated
by the known rational points. The formerly unknown boundary bit is zero.

The certificate is a single point-blind class
\[
\beta=-\operatorname{disc}(f)\,f'(\theta)
       \quad\text{in }K^\times/K^{\times2},
\]
together with local nonmembership and global reciprocity. It gives
\[
\begin{array}{c|c}
\text{fibre}&\dim\mathrm{Sel}_2(E)\\ \hline
\text{MW16-05 }307/206&25+\epsilon\\
\text{R17 high }-2300/843&24+\epsilon\\
\text{R17 low }-1561/3133&17+\epsilon
\end{array}
\]
with a separate \(\epsilon\ge0\) for each curve. Here \(\epsilon\) is exactly
the excess of its \(S\)-class-group 2-rank above 10, 8 or 6, respectively.
It remains uncomputed. These are full-Selmer identities, not numerical
rank upper bounds.

This advances **incidence**: all further Selmer classes can now be represented
by strictly unramified characters after subtracting a known rational class.
Whether any such unknown characters exist, and whether they represent
rational points rather than Sha, remain unresolved.

## A global class defined by the equation

Use the same retained monic integral cubic \(f\), field
\(K=\mathbb Q(\theta)\), and set \(S\) as in the
[strict Selmer identification](STRICT_SELMER_AND_ARTIN_BLOCKS.md).
In all three cases, the complete factorization verifies that every prime
dividing the integral cubic discriminant lies in \(S\). Thus no
nonminimal-model prime has been omitted from the support argument below.

Put \(\delta=\operatorname{disc}(f)\ne0\). Since \(f\) is monic of degree
three,
\[
N_{K/\mathbb Q}(f'(\theta))=-\delta,\qquad
N_{K/\mathbb Q}(-\delta f'(\theta))=\delta^4.
\]
Hence \(\beta\) is a norm-square class and defines an element of
\(H^1(\mathbb Q,E[2])\).

The element \(f'(\theta)\) is integral with norm \(-\delta\).
Therefore \(\beta\) is a unit at every prime of \(K\) above
\(p\nmid\delta\). At odd such primes its quadratic extension is
unramified. Since \(2\in S\), \(\beta\) is unramified outside \(S\).

This construction uses no exceptional point. Under a rational
Weierstrass change \(x=u^2x'+r\), the corresponding cubic derivative
scales by \(u^4\) and its discriminant by \(u^{12}\); their product
scales by the square \(u^{16}\). Thus the class is unchanged under
these labelled model transports. Its large displayed coefficients
are not a new arithmetic source.

When \(\delta>0\), the three ordered real roots give derivative signs
\((+,-,+)\), so \(\beta\) has signs \((-,+,-)\).
The real point image has sign classes \((+,+,+)\) and \((+,-,-)\).
Thus \(\beta\) is outside the real Kummer image. For the R17 high fibre,
\(\delta<0\); it has one real embedding, and this particular real
argument does not apply. Its finite-place nonmembership is checked
separately.

## The reciprocity hyperplane

Let
\[
\Omega_S=\bigoplus_{v\in S}H^1(\mathbb Q_v,E[2]),\qquad
L_S=\bigoplus_{v\in S}\delta_v E(\mathbb Q_v).
\]
The sum of local Tate pairings is nondegenerate, and each local point
image is its own orthogonal complement. Thus
\[
L_S=L_S^\perp,\qquad \dim\Omega_S=2\dim L_S=2\ell.
\]
These are standard local-duality facts; see
[Morgan–Paterson, Section 3.1 and Lemma 4.5](https://arxiv.org/pdf/2011.04374).

Define the functional on the local point product
\[
\lambda_\beta(x)=\sum_{v\in S}\langle\beta_v,x_v\rangle_v.
\]
Global reciprocity forces it to vanish on the localization of every
Selmer class. Indeed, outside \(S\), both classes are unramified and
their local pairing vanishes. Equivalently, in the cubic Kummer
description this is the product formula for the quadratic Hilbert
symbols over all completions of \(K\).

If \(\beta_v\notin\delta_vE(\mathbb Q_v)\) at even one place, local
self-orthogonality implies that \(\lambda_\beta\) is nonzero on \(L_S\).
Let \(W\) be the known rational Kummer subspace. When
\(\dim\operatorname{loc}_S(W)=\ell-1\), containment and dimensions give
\[
\boxed{\operatorname{loc}_S(\mathrm{Sel}_2(E))
       =\operatorname{loc}_S(W)
       =\ker(\lambda_\beta|_{L_S}).}
\]

This is the elementary codimension-one certificate. The class \(\beta\)
is **not** a new Selmer direction: the local nonmembership proves that
it is not a Selmer class. It supplies a dual constraint on which tuples
of locally soluble classes can come from a global Selmer class.

## The entire relaxed boundary is also determined

Let \(R_S\) consist of all \(H^1(\mathbb Q,E[2])\) classes unramified
outside \(S\), with no local point condition imposed at \(S\).
Global reciprocity makes
\(\Lambda=\operatorname{loc}_S(R_S)\) isotropic in \(\Omega_S\).
Consequently \(\dim\Lambda\le\ell\).

Both \(W\) and \(\beta\) lie in \(R_S\). The image of \(W\) has dimension
\(\ell-1\), and \(\beta\) lies outside \(L_S\), hence outside that image.
Their images already span dimension \(\ell\). Therefore
\[
\boxed{\Lambda=
\operatorname{loc}_S(W)+\mathbf F_2\operatorname{loc}_S(\beta).}
\]
This supplies a complete global boundary without computing global units,
the full class group, or the full Selmer basis.
Intersecting with \(L_S\) recovers the previous boxed equality.
The additional derivative direction belongs only to the relaxed group.

## Exact local experiment

The [first protocol](DERIVATIVE_RECIPROCITY_PROTOCOL.json) attempted direct
Hilbert-symbol evaluation on the three completely covered curves.
It completed MW16-05, retaining 37 real/finite Hilbert terms and its
explicit functional. The first harness incorrectly imposed positive
discriminant on all selected rows; the negative-discriminant R17 high
row failed that gate. The R17 low calculation timed out at 40 seconds.
Both incomplete attempts remain UNKNOWN in the
[original inputs](../../artifacts/generated-results/elliptic-curves/rank_jump_derivative_reciprocity_inputs_v1.json)
and [original report](../../artifacts/generated-results/elliptic-curves/rank_jump_derivative_reciprocity_v1.json).

The [second protocol](DERIVATIVE_LOCAL_DUALITY_PROTOCOL.json) uses exact
local membership and Tate duality instead. It admits either discriminant
sign, has a 30-second per-curve cap, and makes no Hilbert-symbol calls.

At each place, the marked generic points already span the full local
point image, as certified in the preceding audit. The computation:

- evaluates the derivative squareclass in the same local cubic algebra;
- selects a local basis from those generic point images;
- checks membership by exact binary elimination;
- checks every possible generic-basis correction through PARI's separate
  local-square interface;
- retains a separating linear character for each outside place.

All three workers completed. The
[local certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_derivative_local_duality_inputs_v1.json)
contains all signatures, basis indices, corrections and separating
characters. The
[consequence report](../../artifacts/generated-results/elliptic-curves/rank_jump_derivative_local_duality_v1.json)
resolves the boundary bit on all three curves.

| Fibre | Discriminant sign | Places tested | Outside places | Simple outside witness | \(\ell\) | Full Selmer local-image dimension |
|---|---|---:|---:|---|---:|---:|
| MW16-05 high | positive | 14 | 10 | real place, also 13 | 16 | 15 |
| R17 high | negative | 13 | 10 | 3 | 17 | 16 |
| R17 low | positive | 10 | 8 | real place, also 5 | 12 | 11 |

There are 87 finite-place correction comparisons in total. Each checks
whether the corrected class is square in every completion above that
rational prime. The real comparisons use exact algebraic signs.
All replayed successfully. No point, parameter, factorization or
class-group search was performed.

## The remaining unknown is entirely strict

Write \(U=\mathrm{Sel}_2^S(E)\),
\(c_S=\dim\mathrm{Cl}(\mathcal O_{K,S_K})/2\), and
\(k=\dim(W\cap U)\). The strict identification gives \(\dim U=c_S\).
The complete image result gives a canonical isomorphism
\[
\mathrm{Sel}_2(E)/W\simeq U/(U\cap W),
\qquad
\dim\mathrm{Sel}_2(E)=\dim W+(c_S-k).
\]
This proves the formulas at the start, replacing the older
\(n+\epsilon+b\) identities by \(n+\epsilon\).

More strongly, let \(A=\delta E(\mathbb Q)\supseteq W\).
Since \(\operatorname{loc}_S(A)=\operatorname{loc}_S(W)\), every
\(\Sha(E)[2]\) class has a strictly unramified representative.
There is an exact sequence
\[
0\longrightarrow (U\cap A)/(U\cap W)
 \longrightarrow U/(U\cap W)
 \longrightarrow\Sha(E)[2]\longrightarrow0.
\]
The first term is isomorphic to \(A/W\). With no rational two-torsion,
\[
\boxed{\epsilon=
  \bigl(\operatorname{rank}E(\mathbb Q)-\dim W\bigr)+\dim\Sha(E)[2].}
\]
Neither summand is determined. In particular, if an independent
\(S\)-class-group upper bound proved \(c_S=k\), it would prove the
exact rank \(\dim W\) and \(\Sha(E)[2]=0\). That upper bound is still missing.

For the completely covered, scale-matched R17 pair this also gives the
exact comparison
\[
\dim\mathrm{Sel}_2(E_{\rm high})-
\dim\mathrm{Sel}_2(E_{\rm low})
=5+c_{S,\rm high}-c_{S,\rm low}.
\]
The cubic fields and sets \(S\) differ; this is a dimension identity,
not an identification of their class groups. Subtracting the marked
generic character contributions remains essential.

## Mechanism assessment

1. **Incidence, now exact:** a canonical equation-defined reciprocity
   constraint describes the complete Selmer localization on these fibres.
   The many separate local choices are coupled by one global condition.
2. **Incidence, still uncomputed:** every remaining Selmer direction is
   an excess \(S\)-class character beyond the retained strict rational
   block. The unknown local boundary dimension has been eliminated.
3. **Solubility, still missing:** the kernel of the map from those excess
   characters to \(\Sha(E)[2]\) determines any further rational directions.
   The derivative constraint does not make a class rational.
4. **Weak explanation:** interpreting an individual new local component,
   or the derivative class itself, as an additional Mordell–Weil point.
   The derivative is locally inadmissible and deliberately used as a dual
   certificate.

For Agent 1, this supplies a precise target for a future independent
incidence calculation: bound \(c_S\), rather than estimate an additional
local-image dimension or increase a chart budget. It does not authorize
a large class-group campaign or provide a current selector. The next
mathematical gate is a reproducible upper bound or independent source
of excess \(S\)-class characters on the R17 pair.

## Replay

From the repository root:
```sh
python3 elliptic-curves/rank-jump/derivative_local_duality.py check
sage -python elliptic-curves/rank-jump/derivative_local_duality.py verify --index 0
sage -python elliptic-curves/rank-jump/derivative_local_duality.py verify --index 4
sage -python elliptic-curves/rank-jump/derivative_local_duality.py verify --index 5
```
The original positive-discriminant Hilbert record is an optional regression:
`sage -python elliptic-curves/rank-jump/derivative_reciprocity.py verify --index 0`.
Its two unsuccessful rows are historical attempts, not replay prerequisites.

The duality and Hilbert implementations use the same labelled cubic.
[PARI documents the local Hilbert routine](https://pari.math.u-bordeaux.fr/dochtml/html-stable/General_number_fields.html#nfhilbert);
every call in the original attempt supplied an explicit prime ideal.
All construction and replay used Sage 10.9 and PARI 2.17.3.
No active-search script, output, protocol or mathematical-status entry
was changed.
