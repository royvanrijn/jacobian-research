# Strict soluble blocks force ordinary ideal-class 2-torsion

The known strict soluble blocks on the two completely covered high fibres
cannot consist mainly of unit squareclasses. Exact half-ideal certificates
force at least **7** and **5** ordinary ideal-class directions beyond the
image of the generic strict subspace. The matched R17 low control gives
no positive relative lower bound from this argument.

This is an **incidence** result derived retrospectively from known rational
points. It distinguishes concrete ideal-class 2-torsion from the unramified
characters and Artin-dual ideal classes constructed earlier. It does not
construct an unknown rational point, identify a particular independent
subset of these half ideals, or bound the full class group.

## Three different class-group objects

Keep the following objects separate:

- A strict Kummer class \(\beta\) defines an unramified quadratic
  character \(\chi_\beta\), trivial at all places over \(S\).
  Thus it is an element of
  \(\operatorname{Hom}(\operatorname{Cl}(\mathcal O_{K,S_K}),\mathbf F_2)\).
- The earlier [small-prime Artin dictionary](STRICT_SELMER_AND_ARTIN_BLOCKS.md)
  supplies ideal classes in \(\operatorname{Cl}(K)/2\) on which these
  characters can be evaluated. Those ideals need not be 2-torsion.
- Since every finite valuation of a strict \(\beta\) is even, its principal
  ideal has a unique fractional-ideal square root:
  \[
  (\beta)=\mathfrak J_\beta^2,\qquad
  \Phi(\beta)=[\mathfrak J_\beta]\in\operatorname{Cl}(K)[2].
  \]
  This last object is certified here.

Multiplying \(\beta\) by a square multiplies \(\mathfrak J_\beta\) by a
principal ideal. Therefore \(\Phi\) is a well-defined linear map on strict
squareclasses. It is neither the Artin pairing nor the Cassels–Tate pairing.

The earlier character argument already gave ordinary class-group 2-rank
lower bounds of 10, 8 and 6. The smaller bounds below are not improvements
to those total lower bounds. They locate substantial 2-torsion in the
specific half-ideal images of the retained rational classes.

## The unit-kernel bound

Let \(V=W\cap U\) be the known rational strict subspace and
\(G_0=G\cap U\subset V\) its marked generic strict subspace.
Write \(\dim V=k\), \(\dim G_0=g_0\).

If \(\Phi(\beta)=0\), then \(\mathfrak J_\beta=(a)\), so
\(\beta/a^2\) is an ordinary unit. Every strict class is totally positive,
and division by a square preserves that property. Hence
\[
\ker(\Phi|_V)
 \subseteq \mathcal O_K^{\times,+}/\mathcal O_K^{\times2}.
\]
For a cubic field of signature \((r_1,r_2)\), Dirichlet's unit theorem gives
\(\dim\mathcal O_K^\times/\mathcal O_K^{\times2}=r_1+r_2\).
The sign map has rank at least one because \(-1\) is negative at every real
place. Consequently
\[
\dim\ker(\Phi|_V)\le u:=r_1+r_2-1.
\]
The standard unit theorem is reviewed in
[Milne, *Algebraic Number Theory*, Chapter 5](https://www.jmilne.org/math/CourseNotes/ANT.pdf).
No fundamental units are computed here.

Linear algebra now gives
\[
\begin{aligned}
\dim\Phi(V)&\ge k-u,\\
\dim\Phi(G_0)&\ge \max(0,g_0-u),\\
\dim\bigl(\Phi(V)/\Phi(G_0)\bigr)
 &= k-g_0-\dim\ker(\Phi|_V)+\dim\ker(\Phi|_{G_0})\\
 &\ge \max(0,k-g_0-u).
\end{aligned}
\]
The subtraction is specifically by the **generic strict** image.
The whole marked generic subgroup need not have even valuations at every
bad prime, so this ordinary half-ideal map has not been defined on that
whole subgroup.

| Fibre | Certified known rank \(n\) | Marked generic rank \(m\) | Known quotient \(n-m\) | Signature | \(k\) | \(g_0\) | \(u\) | \(\dim\Phi(V)\ge\) | \(\dim\Phi(G_0)\ge\) | \(\dim\Phi(V)/\Phi(G_0)\ge\) |
|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|
| MW16-05 \(307/206\) | 25 | 16 | 9 | \((3,0)\) | 10 | 1 | 2 | 8 | 0 | 7 |
| R17 \(-2300/843\) | 24 | 17 | 7 | \((1,1)\) | 8 | 2 | 1 | 7 | 1 | 5 |
| R17 \(-1561/3133\) | 17 | 17 | 0 observed | \((3,0)\) | 6 | 6 | 2 | 4 | 4 | 0 |

All displayed curve ranks are certified lower bounds. Zero observed gain
is censored by the frozen experiment, not an exact zero jump.

On MW16-05 the entire nine-dimensional known exceptional quotient admits
generic corrections into the strict subspace. On R17 high, six of the
seven quotient directions do; one is visible in the joint bad-place image.
The ordinary ideal-class quotient bounds account for at least seven and
five of those strict relative dimensions, respectively. They do not
specify which displayed basis directions survive, nor imply that the
surviving points arise from a common low-degree rational construction.

## Constructing the half ideals without factoring point numerators

Use the retained integral short equation \(y^2=f(x)=x^3+Ax+B\) and
\(K=\mathbf Q(\theta)\). For a rational point write
\[
x=a/d^2,\qquad y=b/d^3,\qquad
\gamma=a-d^2\theta,\qquad N(\gamma)=b^2.
\]
The class of \(\gamma\) is the point's cubic Kummer class.
Form the integral ideal
\[
\mathfrak I=(b,\gamma).
\]

Outside the primes dividing \(2\operatorname{disc}(f)\), the cubic order
is maximal and étale. If a prime divides \(d\), \(\gamma\) is a unit there.
Otherwise at most one prime above a rational prime can divide \(\gamma\):
its residue root is \(a/d^2\), simple and of degree one. The norm identity
then makes its valuation \(2v_p(b)\). Thus
\(v_{\mathfrak p}(\mathfrak I)=v_{\mathfrak p}(\gamma)/2\) away from that
fixed support. This is why factoring the potentially enormous \(b\) is
unnecessary.

For a retained strict mask \(M\), multiply its point elements and ideals:
\[
\gamma_M=\prod_{i\in M}\gamma_i,\qquad
\mathfrak I_M=\prod_{i\in M}\mathfrak I_i.
\]
At each prime ideal \(\mathfrak p\) above the completely factored support,
put
\[
e_{\mathfrak p}
 =\tfrac12v_{\mathfrak p}(\gamma_M)-v_{\mathfrak p}(\mathfrak I_M),
\qquad
\mathfrak J_M
 =\mathfrak I_M\prod_{\mathfrak p}\mathfrak p^{e_{\mathfrak p}}.
\]
Strictness ensures even valuations. The computation checks exactly
\[
\mathfrak J_M^2=(\gamma_M),\qquad
N(\mathfrak J_M)=\prod_{i\in M}|b_i|.
\]
These final identities certify the result independently of the shortcut
used to find the ideal. The records retain the maximal-order basis,
point transports, each gcd ideal, all support corrections, the resulting
HNF matrix and the principal-product HNF.

The [protocol](STRICT_HALF_IDEAL_PROTOCOL.json) fixes three existing fibres,
24 existing strict masks and 40 seconds per fibre. There are no parameter
or point searches, coordinate factorizations, class-group calculations
or ideal enumerations. Every factor needed to certify the maximal order
and support was already retained by the preceding bounded audit.

The first run reached the final field metadata call, then all three
workers failed on an unsupported PARI accessor. Its
[UNKNOWN inputs](../../artifacts/generated-results/elliptic-curves/rank_jump_strict_half_ideal_inputs_v1.json)
and [source, protocol and failure logs](../../artifacts/generated-results/elliptic-curves/rank_jump_strict_half_ideal_failed_attempt_v1.json)
are preserved. Revision 2 fixes that accessor with the same mathematical
calculation and limits. All three cases completed.

## Independent replay

The [version 2 certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_strict_half_ideal_inputs_v2.json)
contains 10, 8 and 6 half ideals; the
[consequence report](../../artifacts/generated-results/elliptic-curves/rank_jump_strict_half_ideals_v2.json)
contains the bounds above.

The verifier does not repeat the gcd-ideal or correction algorithm.
It reads each ideal as a lattice in the integral basis and uses exact
number-field multiplication to check:

- closure under multiplication by all three integral-basis elements;
- containment of all six generating products of its square in
  \(\gamma_M\mathcal O_K\);
- equality of the indices through the exact determinant and norm;
- total positivity, exact point transports, and the integral-basis
  discriminant and maximal-order certificate.

Containment and equal indices prove equality of the ideals.
All 24 certificates passed. Two deliberately altered ideals per fibre
were rejected: one changes the determinant, while the other preserves it
and disturbs the lattice. This checks both the norm and containment gates.

```sh
python3 elliptic-curves/rank-jump/strict_half_ideals.py check
sage -python elliptic-curves/rank-jump/verify_strict_half_ideals.py --index 0 --negative-controls
sage -python elliptic-curves/rank-jump/verify_strict_half_ideals.py --index 4 --negative-controls
sage -python elliptic-curves/rank-jump/verify_strict_half_ideals.py --index 5 --negative-controls
```

Replay used Sage 10.9 and PARI 2.17.3. The checks are confined to these
retained curves and ideals.

## Mechanism assessment and next missing implication

1. **Incidence, strongest established structure:** most of the known
   exceptional quotient can be moved into one strict unramified block,
   and most of that block has a nontrivial ordinary half-ideal image.
   A unit-only explanation is impossible at these dimensions.
2. **Incidence, point-independent target:** the equation defines the cubic
   field and its \(S\)-class group before exceptional points are supplied.
   The [complete-boundary theorem](DERIVATIVE_RECIPROCITY_AND_COMPLETE_BOUNDARY.md)
   reduces all remaining Selmer dimensions to the excess \(S\)-class
   character space. A certified upper bound or independent character
   construction remains missing.
3. **Solubility, unresolved:** the present ideal classes are obtained from
   already known rational points. Generating suitable ideal classes
   independently would still require a square generator satisfying the
   strict conditions and a proof that its genus-one cover has a rational
   point. Simultaneous solubility does not follow from the ideal identity.
4. **Weak explanations:** small ordinary unit rank cannot account for the
   observed strict block; a common cubic field alone supplies no rational
   construction. Ordinary class-group size by itself also fails to
   distinguish generic directions, exceptional directions and Sha.
5. **Visibility:** half-lattice chart concentration remains a recovery
   diagnostic. None of the new ideal quantities is a prospective chart
   score or a validated rank predictor.

A critical limitation is localization:
\[
\operatorname{Cl}(K)\longrightarrow
\operatorname{Cl}(\mathcal O_{K,S_K})
\]
kills classes supported at \(S_K\). The certified half-ideal images may
lose dimensions under this map. Their lower bounds therefore cannot
replace an \(S\)-class-group upper bound or determine the unknown
\(\epsilon\) in the full-Selmer formulas.

For Agent 1, the useful information is a concrete arithmetic target with
a quantified unit ambiguity of at most two or one dimensions, not a
selector ready to deploy. The relative certificates use exceptional
points and must stay retrospective. Before any candidate-selection use,
one must construct the relevant classes from the equation alone and
separately test their rational solubility.

The missing chain is still
\[
\text{an independently detectable specialization condition}
\ \Longrightarrow\
\text{a large strict arithmetic block}
\ \Longrightarrow\
\text{many simultaneously rational covers}.
\]
These certificates establish a necessary structure in known soluble
examples. They do not establish either forward implication.
