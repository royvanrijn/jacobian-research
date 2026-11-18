# Good degree-twenty-six models exist locally at2,3 and over the reals

The high-genus accessory curve has infinitely many good points over
Q2 and Q3. Exact rescaled charts at infinity give simple residue roots
and uniform Hensel families. Its two simple real points at infinity
also give infinitely many good real points. Thus neither2,3 nor the
real place supplies a local obstruction to a rational point in this
affine-trace slice.

This is local solubility, not rational descent. It neither constructs
a good Q-rational model nor determines the curve's finite set of
rational points. Real coefficients do not imply that all critical
points are real.

## The normalization and its infinity chart

Let F(s,t) be the primitive integer accessory polynomial from the
[high-genus theorem](Q26_HIGH_GENUS_AFFINE_TRACE_CURVE.md). It is
geometrically integral, and

\[
 F=\sum_{i=0}^{12}a_i(t)s^i,\qquad \deg a_i=144-i.
\]

The frozen [coefficient input](q26_high_genus_affine_trace_curve_inputs.json)
binds it exactly to the normalized affine-trace construction.
Set

\[
 x=s/t,\quad v=1/t,\qquad
 G(x,v)=v^{144}F(x/v,1/v)
       =\sum_{i,j}a_{ij}x^iv^{144-i-j}\in\mathbb Z[x,v].
 \tag{1}
\]

For v!=0 this is an invertible rational coordinate change. At v0,
the polynomial is A(x)=sum_i LC(a_i)*x^i.

There is an important trap in the naive total-degree144 plane closure.
It always contains the exact rational point(S:T:U)=(1:0:0). In its
S1 chart every term has(T,U)-degree at least132, so this point is
singular. It gives a primitive homogeneous root modulo every prime
power regardless of whether a local branch of the normalization is
rational. Hence excluding all primitive projective roots cannot
provide an obstruction for this unresolved plane model. We use smooth
points of the explicit birational chart instead.

## A uniform toric Hensel criterion

**Lemma 1.** Let p be a prime and m,a,b be integers. Suppose

\[
 H(X,V)=p^mG(p^aX,p^bV)\in\mathbb Z_p[X,V]
\]

reduces to P(X), independent of V, and P has a simple nonzero root
r in Fp. Then for every V in Zp there is a unique X(V) congruent r
with H(X(V),V)=0. These points are smooth. For nonzero V they yield
points of F=0 with

\[
                t=p^{-b}/V,\qquad s=p^{a-b}X(V)/V.        \tag{2}
\]

Infinitely many of them lie in every nonempty geometric good open
of this integral accessory curve.

**Proof.** Fix V. The polynomial H(X,V) has integral coefficients,
and its root r modulo p has unit X-derivative. The simple-root
Hensel lemma, [Milne ANT3.00, Proposition7.31,p113](sources/milne-ant300.pdf),
gives the unique lift in r+pZp. The derivative remains a unit, so
the characteristic-zero curve point is smooth. For V!=0, (1)--(2)
give a smooth affine point and therefore a rational point on the
normalization. Varying V through units gives infinitely many distinct
t-values. The complement of a nonempty geometric open in an integral
curve is finite, so only finitely many of these points can be bad. ∎

No completeness of a residue search is needed for this positive
criterion. Its coefficient integrality and reduction are exact
identities, checked on every monomial of F.

## The two small-prime charts

**Theorem 2.** The following charts satisfy Lemma1:

| p | H(X,V) | Reduction modulo p | Simple root |
| --- | --- | --- | --- |
| 2 | G(X,2V) | X^11(X+1) | X1, derivative1 |
| 3 | 3^18 G(X/9,9V) | X^11(X+1) | X2, derivative2 |

In particular the normalization has infinitely many good Q2 and
Q3 points.

**Proof.** For a monomial c*x^i*v^k of G, the transformed coefficient
is c*p^(m+ai+bk). The exact checker processes all1807 nonzero
monomials. In the2-chart every coefficient is integral and the only
nonzero reduction terms are X^11 and X^12. In the3-chart it verifies

\[
       v_3(c)+18-2i+2k\ge0
                 \quad\text{for every monomial},
\]

and exact integer division when18-2i+2k<0. Again the only nonzero
residue terms are X^11 and X^12. Thus the claimed reductions are
global polynomial identities, not evaluations at selected points.
Both listed roots have unit derivative. Lemma1 now applies. ∎

For V a unit, the resulting coordinates have exact valuations

\[
 \begin{array}{c|cc}
       &v_p(t)&v_p(s)\\ \hline
 p=2&-1&-1\\
 p=3&-2&-4.
 \end{array}
\]

In particular the3-adic points escape integral affine searches.
At V0 the3-chart gives a genuine smooth Q3 branch above infinity,
whose image in the naive projective model reduces to its singular
point(1:0:0). This is exactly why that singular residue class cannot
be discarded.

The geometric good open is nonempty by the
[uniform theorem's exact Q29/number-field point](ODD_AFFINE_CRITICAL_TRACE_CURVES.md).
Geometric integrality follows from the high-genus theorem. Its bad
locus, including recurrence denominators and all critical-pair
nondegeneracy exclusions, is therefore finite. Hence the infinitely
many local smooth points above include infinitely many actual good
normalized paired pencils. Their Laurent monodromy is A26 by the
uniform theorem; this is distinct from the accessory cover's S12.

## The real place

**Proposition 3.** A has exactly two distinct real roots, both simple.
The normalization and its good open have infinitely many real points.

**Proof.** Divide the common integer coefficient content of A and
choose positive leading coefficient. The retained coefficient input
gives a degree-twelve polynomial. Exact rational Euclidean arithmetic
checks gcd(A,A')=1 and constructs the complete Sturm sequence.
The sign-variation difference at negative and positive infinity is
exactly two. Thus both real roots are simple.

At either root alpha, G(alpha,0)=0 and G_x(alpha,0)!=0. The real
implicit function theorem gives a smooth branch x(v) for v near0.
Its nonzero v-values yield infinitely many distinct finite points of
F=0. Removing the finite bad locus leaves infinitely many good real
points, just as in Lemma1. ∎

This argument makes no assertion about all critical points being real,
real-root-count spectra, or global real connected components.

## Evidence and limits

The [checker](verify_q26_small_prime_local_branches.py) binds all
chart coefficients to the frozen primitive F, verifies every integer
division and full residue polynomial, supplies unique Hensel digit
controls through precision8 at V1, and computes the full exact Sturm
sequence. Its [receipt](q26_small_prime_local_branches_verification.json)
and [source ledger](q26_small_prime_local_branches_source_ledger.json)
record the source pins and theorem imports. The written Hensel argument
supplies all p-adic parameters; finite digit controls do not establish
uniformity by sampling.

No complete local classification at every prime, rational Q point,
Hasse principle, effective rational-point determination or further
global existence theorem is established. The earlier high-genus,
operator, spectral and holonomy packages remain unchanged. Formal
verification, external human review, independent external replay and
priority claims remain absent.
