# Independent halving fields and arithmetic blocks, 6 September 2026

**All 91 main-panel observations have maximal joint halving extensions for
their specified independent witness subgroup.** The exceptional directions
do not share a smaller halving extension over the two-division field.
This is forced by their certified mod-2 independence and full \(S_3\) action,
not a newly discovered feature of high-rank curves.

More decisively, the same finite class fields persist in the fixed-cubic
control while eighteen rationally soluble anchor classes become a subspace
with CT rank sixteen. Shared first-level class fields therefore cannot
distinguish simultaneous rational solubility there.

The [protocol](HALVING_FIELD_PROTOCOL.json) and
[certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_halving_fields_v1.json)
include the finite-module proof checks, all panel degrees, and nine explicit
halving quartics from the three paired studies. Replay:

    python3 elliptic-curves/rank-jump/halving_fields.py check

The replay reevaluates point membership and finite Kummer signatures on 76
rows. Fifteen historic rows retain the earlier exact global-fingerprint
certificates and replay their binary algebra. No number-field compositum or
new point was computed.

## The exact first-level theorem

Let \(V=E[2]\), assume \(\operatorname{Gal}(\mathbf Q(V)/\mathbf Q)=S_3
=\mathrm{GL}_2(\mathbf F_2)\), and put \(L=\mathbf Q(V)\).
Let \(c_1,\ldots,c_n\) be independent in \(H^1(\mathbf Q,V)\).
Their restrictions to \(G_L\) are homomorphisms
\[
\kappa_i:G_L\longrightarrow V.
\]
Let \(F_W/L\) be the finite extension cut out by their joint kernel, where
\(W=\langle c_1,\ldots,c_n\rangle\). Then
\[
\operatorname{Gal}(F_W/L)=V^n,\qquad [F_W:L]=4^n.
\]
For \(n>0\), its Galois group over \(\mathbf Q\) is the full affine group
\[
V^n\rtimes S_3
\]
with the same natural linear action in each factor.

When \(c_i=\delta(P_i)\) for rational points, \(F_W\) is exactly
\(L(\tfrac12P_1,\ldots,\tfrac12P_n)\). These Kummer and affine representations
are described in [Lombardo–Tronto, §2.3](https://arxiv.org/pdf/1909.05376).
The following elementary argument gives the exact mod-2 conclusion needed
here, without using general asymptotic Kummer bounds.

**Proof.** The normal subgroup \(C_3\subset S_3\) has no invariants in \(V\),
and its positive-degree cohomology over \(\mathbf F_2\) vanishes by averaging.
Inflation-restriction gives \(H^1(S_3,V)=0\). Thus
\[
H^1(\mathbf Q,V)\longrightarrow H^1(L,V)
\]
is injective.

The joint image \(H\subset V^n\) is \(S_3\)-stable. The algebra spanned by
the natural \(S_3\) matrices is all of \(M_2(\mathbf F_2)\). Its matrix units
show that every stable subspace of \(V^n\) has the form \(V\otimes R\), for
some \(R\subset\mathbf F_2^n\): the projections onto the two coordinate rows
and the maps exchanging them identify a common row subspace.

If \(H\ne V^n\), a nonzero linear functional
\(\lambda=(\lambda_1,\ldots,\lambda_n)\) annihilates \(R\). Therefore
\(\sum_i\lambda_i\kappa_i=0\). Injectivity of restriction gives
\(\sum_i\lambda_ic_i=0\), contradicting independence. Hence \(H=V^n\).
The affine representation is contained in \(V^n\rtimes S_3\), contains the
full kernel, and surjects onto \(S_3\); it is the whole group. \(\square\)

The finite certificate independently enumerates all 1024 normalized
functions \(S_3\to V\): there are four 1-cocycles, exactly the four
1-coboundaries. It verifies the full matrix algebra and checks all 5, 67
and 2825 binary subspaces of \(V,V^2,V^3\); exactly 2, 5 and 16 are stable,
all of the predicted form. These checks support the general proof rather
than replace it by a bounded assertion.

If \(U,W\) are finite class subspaces, the same theorem and Galois
correspondence give
\[
F_UF_W=F_{U+W},\qquad F_U\cap F_W=F_{U\cap W},
\]
with \(F_0=L\). In particular, fields attached to independent summands are
linearly disjoint over \(L\).

## Applying it to the high-gain fibres

Let \(M\) be the marked generic subgroup and \(A\) the subgroup generated
by the independent indices explicitly selected in the panel. Write
\[
m=\dim\delta(M),\qquad r=\dim\delta(A),\qquad q=r-m.
\]
Then
\[
[F_{\delta(A)}:F_{\delta(M)}]=4^q.
\]
Every main-panel row has \(m=\operatorname{rank}M\) and
\(r=\operatorname{rank}A\). The field formula nevertheless uses the
**mod-2 dimensions**, not an unqualified free rank.

| Fibre | \(m\) | Certified \(q\) | Relative halving degree |
|---|---:|---:|---:|
| ICARM-398 | 16 | 14 | \(4^{14}=268435456\) |
| published-R17 control-r27 | 17 | 10 | \(4^{10}=1048576\) |
| 07ca9, \(-2507/3068\) | 17 | 9 | \(4^9=262144\) |
| MW16-04, \(-1647/91\) | 16 | 9 | \(4^9=262144\) |
| MW16-05, \(307/206\) | 16 | 9 | \(4^9=262144\) |
| 07ca9, \(505/794\) | 17 | 8 | \(4^8=65536\) |

These degrees are exact for the stated witness subgroups. They do not
describe the halving field of the unknown entire Mordell–Weil group.
Unselected points that share a finite signature with a selected point are
not silently added to \(A\).

An observed-zero control has relative degree 1 because this particular
witness subgroup supplies no new independent class. That remains a censored
search observation, not an upper bound on its curve's rank or class field.

The three older MW12 embeddings are excluded from this calculation.
Their free quotient ranks and Smith forms alone do not certify the ambient
global mod-2 images needed here. Likewise, no new family-relative statement
is made for the masked record-273/302 controls.

## Explicit quartics in the three paired studies

For \(E:y^2=x^3+Ax+B\) and \(P=(p,q)\), \(q\ne0\), the four x-coordinates
of halves of \(P\) are roots of
\[
h_P(Z)=Z^4-4pZ^3-2AZ^2-(4Ap+8B)Z+A^2-4Bp.
\]
The corresponding y-coordinate is recovered by
\[
y_Q=\frac{Z^3-3pZ^2-AZ-Ap-2B}{2q}.
\]
Thus this polynomial records the full four-point halving torsor. Its
discriminant and cubic resolvent satisfy
\[
\operatorname{disc}(h_P)=2^{12}q^4\operatorname{disc}(x^3+Ax+B),
\]
and a root of the standard cubic resolvent is
\[
\rho=2A+4p\theta+4\theta^2.
\]
The replay substitutes this generator into the resolvent exactly.
All halves therefore have the expected common cubic resolvent field; this
does not reduce the independent degree-four extensions above \(L\).

For each original high/low pair, the protocol selects the first generic
point on each curve and the first selected exceptional point on the high
curve. Every one of these nine quartics is independently certified \(S_4\)
using good reductions with cycle types \((4)\) and \((1,3)\).

| Fibre | Point role | Prime for \((4)\) | Prime for \((1,3)\) |
|---|---|---:|---:|
| MW16-05, \(307/206\) | generic | 37 | 29 |
| MW16-05, \(307/206\) | exceptional | 11 | 29 |
| MW16-05, \(-3158/1291\) | generic control | 11 | 19 |
| MW16-04, \(-1647/91\) | generic | 11 | 13 |
| MW16-04, \(-1647/91\) | exceptional | 11 | 13 |
| MW16-04, \(-2177/2397\) | generic control | 5 | 17 |
| published-R17, \(-2300/843\) | generic | 23 | 17 |
| published-R17, \(-2300/843\) | exceptional | 23 | 17 |
| published-R17, \(-1561/3133\) | generic control | 23 | 17 |

The coinciding cycle types do not identify coinciding fields. Joint
independence is established by the theorem, not inferred from a few primes.

## A shared first-level field can conceal opposite solubility

The theorem applies to **arbitrary independent classes in**
\(H^1(\mathbf Q,E[2])\), before rational solubility is known. Each class
defines a finite étale four-element \(E[2]\)-torsor. Its splitting field is
distinct in meaning from a rational point on its associated genus-one
two-cover. A rational point on the finite torsor would instead mean that
the class itself is zero.

In the fixed-cubic control, take the eighteen-dimensional class space
\(W_{-1}\) and compare \(E_0\) with \(E_{-1}\). The labelled \(E[2]\)-module
and all eighteen classes are identical. Their joint finite-torsor field
is **the same field**, of degree \(4^{18}\) over \(L\), at both parameters.
Yet all eighteen classes are represented by rational points on \(E_0\),
whereas the CT form on \(E_{-1}\) has rank sixteen and its restricted radical
has dimension two.

What varies is the embedding of the fixed module into the elliptic curve,
and hence the genus-one covers and higher descent. Neither the entire
first-level class field nor its abstract Galois group decides this
solubility change. This strengthens the earlier observation that the
two-division cubic alone is insufficient.

It also explains why field independence does not establish that directions
arise as independent accidents. A common auxiliary curve could construct
many rational points whose halving fields are maximally independent. Any
partition of an independent basis already gives independent first-level
fields; their existence supplies no preferred arithmetic block partition.

## The precise finite-prime correlation forced by shared two-torsion

There is one predictable correlation that should be removed from local
halving diagnostics. Conditional on the common linear Frobenius
\(g\in S_3\), the \(n\) translation coordinates are independent and uniform
in \(V\). A half of each point exists modulo a good unramified prime exactly
when each affine action \(x\mapsto gx+v_i\) has a fixed point.

| Linear Frobenius | Density among primes | Probability all \(n\) factors have a fixed point |
|---|---:|---:|
| identity | \(1/6\) | \(4^{-n}\) |
| transposition | \(1/2\) | \(2^{-n}\) |
| 3-cycle | \(1/3\) | \(1\) |

The finite-group count, and therefore the Chebotarev density, is
\[
\frac13+\frac{1}{2^{n+1}}+\frac{1}{6\cdot4^n}.
\]
For \(n=1,2,3\) it is \(5/8,15/32,51/128\), independently checked by
enumerating affine groups of orders 24, 96 and 384. At a prime with
irreducible two-division cubic, doubling on \(E(\mathbf F_p)\) is bijective;
all rational points reduce to divisible points together.

These are distributions over primes for a fixed curve and fixed independent
classes. They are **not probabilities of a rank jump across specializations**.
They also describe mod-prime halving, not global solubility of the
genus-one two-cover. Treating their shared \(1/3\) contribution as a new
multi-direction event would count the already fixed two-torsion action.

## Mechanism ranking and next test

1. **Secondary descent or a common auxiliary construction forcing several
   rationally soluble classes** remains the strongest candidate mechanism.
   The present first-level field calculation neither constructs nor excludes it.
2. **Correlated local incidence at bad primes** remains a concrete source
   of class-space changes. The fixed-field collision result explains one
   version, but the recent high-gain fibres still need their own bad-prime
   quotient-support comparison.
3. **Collapsed shared halving fields** are ruled out for the certified
   independent subgroups. Constant first-level class fields, their abstract
   groups, and ordinary prime-halving correlations cannot distinguish the
   soluble and obstructed fixed-field controls.

The next bounded retrospective experiment should use the same six paired
high/low fibres and the fixed odd primes \(3,5,7,11,13\). At primes where the
model has bad reduction, compare the exact local squareclass image of the
marked generic subgroup with the image after adding the certified quotient
directions. This fills a gap left by the earlier good-prime signature audit.
The endpoint is a new local quotient dimension beyond the generic image,
not a rank estimate. Full factorization, full Selmer groups and prospective
parameter sweeps are unnecessary. A null result would exclude only this
small support dictionary.

Agent 1 receives no new scoring rule from the field theorem. Its practical
use is to avoid spending selection effort on a field collapse or finite-prime
correlation already ruled out or forced by known independence.
