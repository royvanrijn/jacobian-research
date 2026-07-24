# Arithmetic Keller-fiber engineering

## 1. The programme

The organizing question is:

> Which finite etale schemes over a number field occur as complete fibers of
> polynomial Keller maps?

For the weighted threefold construction, the arithmetic object is not merely
a point count. On the chart \(C=1\), the entire fiber is

\[
 \operatorname{Spec}K[W]/(H(W)-sW+t).
\]

This identifies one subject linking the repository's rational-fiber,
real-sheet, Chebotarev, adelic, Hasse, and monodromy results. The basic
dictionary is:

| Fiber property | Polynomial property |
|---|---|
| finite etale of degree \(n\) | squarefree inverse polynomial of degree \(n\) |
| connected fiber | irreducible inverse polynomial |
| real signature | number of real roots |
| unramified local degrees | factor degrees modulo a good prime |
| rational point | linear factor over the ground field |
| everywhere local point | a root in every completion |
| Hasse failure | nontrivially intersective inverse polynomial |
| splitting-field group | Galois group of the inverse polynomial |

The subject has two levels which should remain separate:

1. **weighted occurrence**, where a sharp generator criterion is available;
2. **intrinsic Keller occurrence**, allowing every polynomial Keller map,
   where necessity beyond etaleness is presently unknown.

## 2. Exact transfer on the weighted chart

Let \(K\) be a characteristic-zero field and let \(P\in K[W]\) be squarefree
of degree \(n\ge3\). Define

\[
 H(W)=P(W)-P(0)-P'(0)W.
\]

The polynomial \(P\) is **tangent-admissible** if

\[
\begin{aligned}
P(1)-P(0)&=P'(0),\\
P'(1)&\ne P'(0),\\
P''(1)&\ne2\bigl(P'(1)-P'(0)\bigr).
\end{aligned} \tag{2.1}
\]

The first equation says \(H(1)=0\). The two inequalities say respectively
that the Keller determinant is nonzero and that the chosen weighted source
chart does not hit its exceptional value.

### Exact transfer theorem

Every tangent-admissible \(P\) occurs as a complete regular fiber of an
explicit polynomial Keller map of \(\mathbb A^3_K\). More precisely, put

\[
 c=P'(0)-P'(1),\qquad
 q=\left(\frac{P(0)}c,-P'(0),1\right).
\]

The weighted map attached to \(H\) has determinant \(c\), and its inverse
polynomial at \(q\) is exactly \(P\). Because \(C=1\) and \(P\) is
squarefree, reconstruction identifies the complete scheme fiber with

\[
 \operatorname{Spec}K[W]/(P).
\]

This is already the general transfer theorem needed for any intersective
polynomial which can be put on (2.1) by a \(K\)-affine change of generator.
It transfers all local roots, absence of a global root, real signatures,
and the splitting field unchanged.

## 3. What tangent normalization really restricts

Let \(A\) be a degree-\(n\) finite etale \(K\)-algebra. Since \(K\) is
infinite, \(A\) is monogenic. For a generator \(\theta\in A\), let

\[
 P_\theta(T)=\operatorname{Norm}_{A/K}(T-\theta).
\]

Then \(A\) occurs in the weighted chart exactly when its generator open
contains a \(\theta\) satisfying (2.1). Thus tangent normalization is:

- not a restriction on abstract etaleness;
- not visibly an isomorphism invariant of \(A\);
- one hypersurface equation, plus two open conditions, in the primitive
  element space of \(A\).

When \(\theta\) is a unit, the equation has the intrinsic norm-trace form

\[
 \operatorname{Norm}(1-\theta)
 =(-1)^n\operatorname{Norm}(\theta)
   \left(1-\operatorname{Tr}(\theta^{-1})\right). \tag{3.1}
\]

For an affine reparameterization \(T=a+bW\), \(b\ne0\), the equation becomes

\[
 P(a+b)-P(a)=bP'(a). \tag{3.2}
\]

Equivalently, the tangent line to the graph of \(P\) at \(a\) meets the graph
again at the \(K\)-rational point \(a+b\). After division by \(b^2\), (3.2)
is a degree-\((n-2)\) equation. This formulation is the correct starting
point for an exact classification.

The main open normalization question is therefore:

> Does every finite etale algebra of degree at least three have a primitive
> element satisfying (3.1), or what arithmetic obstruction prevents one?

A quadratic stabilization theorem now gives a uniform partial answer.
Every nontrivially intersective polynomial over a number field becomes
tangent-admissible after an affine change of variable and multiplication by
one irreducible quadratic. The tangent equation is linear in the three
quadratic coefficients. Its two-parameter solution plane has generically
nonsquare discriminant; Hilbert irreducibility supplies an irreducible
quadratic while avoiding common roots and both exceptional weighted loci.
Thus every Hasse-failing finite etale scheme transfers to a complete Keller
fiber after adjoining one quadratic etale component, with degree overhead
two. See
[INTERSECTIVE_POLYNOMIAL_TRANSFER.md](verified/INTERSECTIVE_POLYNOMIAL_TRANSFER.md).

## 4. Current frontier status

| Target | Status after this audit | Next mathematical step |
|---|---|---|
| general intersective transfer | **settled with quadratic stabilization and degree overhead two** | decide whether stabilization is ever intrinsically necessary |
| minimal Hasse-failing degree | **settled: \(d_{\mathrm{HP}}=5\)** | package the theorem and compare minimal presentations |
| one fixed map, infinitely many Hasse fibers | **settled by the quadratic-tilt family** | quantify and generalize the resulting target line |
| real signatures plus Frobenius types | already proved for every degree and every finite set of unramified types | extend from unramified cycle types to ramified local etale algebras |
| prescribed solvable Galois groups | arithmetic existence known for coverable solvable groups; weighted transfer needs normalization control | combine Sonn's realizations with tangent-admissible specialization or stabilization |
| tangent-normalization restriction | exact generator criterion (3.1)-(3.2) obtained; existence classification open | study rational points on the primitive-element tangent hypersurface |

## 5. The minimal Hasse degree is five

Berend and Bilu's polynomial

\[
 f(X)=(X^3-19)(X^2+X+1)
\]

is intersective, has no rational root, and has the minimum possible degree
among such polynomials. The affine change

\[
 X=10-27W
\]

puts it exactly on the tangent slice. After primitive scaling,

\[
\begin{aligned}
P(W)={}&-531441W^5+1003833W^4-758889W^3\\
      &\quad+286497W^2-53901W+4033,
\end{aligned}
\]

and

\[
P(1)-P(0)=P'(0),\qquad
\frac{H''(1)}{P'(0)-P'(1)}=-\frac{586}{79}\ne-2.
\]

The resulting map has determinant \(345546\), target

\[
\left(\frac{4033}{345546},53901,1\right),
\]

and complete fiber \(\operatorname{Spec}\mathbb Q[W]/(P)\). The detailed
proof and exact checker are in
[MINIMAL_HASSE_PRINCIPLE_KELLER_FIBER.md](verified/MINIMAL_HASSE_PRINCIPLE_KELLER_FIBER.md).

The arithmetic lower bound is independent of Keller geometry: a Hasse-failing
complete regular Keller fiber is a finite etale scheme with the same
local-global property, and no such scheme has degree below five.

## 6. A fixed map with infinitely many Hasse failures

This target is now settled.  The root-engineered quadratic gauge has inverse
pencil

\[
 G_P(S)-\frac{g_1}{2}(BS^2+C).
\]

Take

\[
G(S)=S^5-\frac32S^4+\frac32S^3-\frac54S^2+\frac9{16}S.
\]

The resulting single polynomial map has determinant \(-2\).  On the rational
target line

\[
(P,B,C)=\left(1,\frac{32a}{9},\frac{8a+1}{3}\right),
\]

its complete inverse polynomial is

\[
\left((S-\tfrac12)^3-a\right)\left(S^2+\frac34\right).
\]

After \(X=S-\tfrac12\), this is the classical intersective family

\[
(X^3-a)(X^2+X+1).
\]

For every rational prime \(a=\ell\equiv1\pmod {27}\), the polynomial is
squarefree and has no rational root, but has a root in every
\(\mathbb Q_p\) and in \(\mathbb R\).  Dirichlet's theorem supplies
infinitely many such primes.  Hence one fixed Keller map has infinitely many
everywhere locally soluble rational target fibers without a rational point.
For the standard projective height the targets have \(H(y_\ell)=32\ell\);
the constructed family therefore satisfies the quantitative asymptotic

\[
\#\{y_\ell:H(y_\ell)\le B\}
=\pi(B/32;27,1)
\sim\frac{B}{576\log B}.
\]

At every good finite-field prime, the same target line is almost entirely
contained in the image: the fixed quadratic splits when \(p\equiv1\bmod3\),
while cubing is bijective when \(p\equiv2\bmod3\).  This is an explicit
exceptional-locus correction to the ambient random-permutation law.
The complete proof and checker are in
[INFINITE_HASSE_KELLER_FIBERS.md](verified/INFINITE_HASSE_KELLER_FIBERS.md).

### The linear-pencil route remains arithmetically interesting

For a fixed seed \(H\), every target on \(C=1\) is the pencil

\[
 E_{s,t}(W)=H(W)-sW+t.
\]

In degree five, a Hasse-failing fiber must have the orbit pattern \(2+3\).
The unique minimal group-theoretic mechanism is the natural \(S_3\)-action:

- transpositions fix a point in the cubic orbit;
- \(3\)-cycles act trivially on the quadratic discriminant orbit.

This gives a concrete construction problem. Choose a monic quadratic

\[
 q(W)=W^2+aW+b.
\]

Division of \(H\) by \(q\) has a linear remainder. There is therefore a
unique pair \((s(a,b),t(a,b))\) for which \(q\mid E_{s,t}\), and the quotient
is cubic. Impose:

1. the quadratic field of \(q\) equals the discriminant field of the cubic;
2. both factors are irreducible;
3. every ramified prime has a local root in one factor;
4. the inverse polynomial is squarefree.

Condition 1 cuts out an explicit resolvent surface with a square parameter.
The first elliptic slice has positive rank, as described below.

For the new minimal seed, the first one-dimensional slice is already
explicit. Fix

\[
 a=-\frac79,
\]

the linear coefficient of the known quadratic factor. Equality of the
quadratic field with the cubic discriminant field becomes the genus-one
quartic

\[
\begin{aligned}
Y^2={}&3(57395628b^3-26749197b^2\\
     &\qquad+4181544b-219512)(49-324b). \tag{6.1}
\end{aligned}
\]

The known Hasse fiber gives the rational point

\[
 (b,Y)=\left(\frac{37}{243},19\right).
\]

An independent binary-quartic Jacobian calculation gives the minimal model

\[
 y^2+xy+y=x^3-x^2-59x-548301,
\]

whose Mordell--Weil rank is one and whose torsion subgroup is trivial
(PARI/GP `ellrank` and `elltors`).  Thus the slice has infinitely many
rational factorizations with matching discriminant squareclass. This does
**not** by itself prove that infinitely many are intersective: simultaneous
irreducibility is a thin arithmetic condition, and at primes ramified in
the varying quadratic discriminant field one must still exclude a full
\(S_3\) decomposition group. The quadratic-tilt construction above avoids
that moving-ramification problem by keeping the quadratic field
\(\mathbb Q(\sqrt{-3})\) fixed.

Two further rational points on this elliptic slice have now been checked at
every ramified prime and give two additional Hasse-failing fibers of the
same fixed weighted map. Their exact factors, targets, and local certificates
are in
[INTERSECTIVE_POLYNOMIAL_TRANSFER.md](verified/INTERSECTIVE_POLYNOMIAL_TRANSFER.md).

## 7. Signatures and Frobenius data

The repository's
[adelic theorem](verified/ADELIC_FIBER_ENGINEERING.md) already proves the
finite unramified version of the proposed target. For every degree \(N\ge3\),
every signature, and every finite collection of admissible Frobenius cycle
types at good primes, one rational target of the fixed map \(F_N\) realizes
all data simultaneously. One imposed \(N\)-cycle makes the fiber connected,
so it is a number field of the chosen signature.

The stratified version also allows the seed to lie on prescribed rational
geometric strata. What remains genuinely new is:

- ramified local algebras rather than squarefree residue types;
- decomposition and inertia groups, not only Frobenius conjugacy classes;
- Grunwald-Wang-type compatibility for global Galois fibers.

## 8. Prescribed solvable groups

The arithmetic input is mature. Sonn proves that every finite solvable group
which is a union of conjugates of suitable proper subgroups occurs as the
Galois group of an intersective polynomial:
[Polynomials with roots in \(\mathbb Q_p\) for all \(p\)](https://arxiv.org/abs/math/0612528).
Related work treats optimal numbers of factors and non-solvable families:
[Bubboloni-Sonn](https://arxiv.org/abs/1507.08593) and
[Koenig](https://arxiv.org/abs/1605.07802).

There are two distinct Keller goals:

1. realize a prescribed solvable group as the splitting-field group of an
   arbitrary complete fiber;
2. realize it while the fiber is Hasse-failing.

The second requires the group to be coverable by the stabilizers of the
fiber components and requires ramified-prime local conditions. The first has
no coverability requirement, but still needs tangent normalization without
unwanted enlargement of the splitting field. A stabilization by a quadratic
factor is therefore harmless for Hasse solubility but may change the Galois
group; this distinction must be explicit in any theorem.

## 9. Recommended paper sequence

### Paper A: minimal and transfer

- exact weighted transfer theorem;
- norm-trace and secant-tangent normalization criteria;
- the degree-five construction and \(d_{\mathrm{HP}}=5\);
- the quadratic-stabilized transfer theorem over arbitrary number fields.

This is the fastest paper because the central new example and its complete
Keller verification are already exact.

### Paper B: adelic Keller fibers

- consolidate real chambers, effective finite-field Chebotarev, constructive
  weak approximation, and stratified seed variation;
- add ramified local specifications;
- formulate the moduli problem of finite etale algebras in one Keller pencil.

### Paper C: infinite Hasse families and Galois engineering

- present the fixed quadratic-tilt family and quantify its target line;
- use the positive-rank linear-pencil resolvent curve to seek a second,
  geometrically distinct infinite Hasse family;
- combine coverable-group inverse Galois constructions with tangent transfer;
- classify the solvable groups obtainable in connected, disconnected, and
  Hasse-failing fibers.

## 10. Immediate computations

The highest-value computations are:

1. map generators of the rank-one Jacobian of (6.1) back to the quartic and
   test ramified local solubility along its Mordell--Weil sequence;
2. search the full quintic discriminant surface for rational genus-zero
   slices or positive-rank elliptic slices;
3. enumerate tangent-normalizable affine changes of standard intersective
   polynomials in degrees \(5\) through \(10\);
4. characterize finite etale algebras which meet the tangent hypersurface
   without quadratic stabilization;
5. build a local-condition checker using decomposition groups at the finite
   set of ramified primes.

The fixed-map existence claim is now proved.  The two remaining high-impact
claims are unconditional transfer of arbitrary finite etale schemes and a
quantitative theory of the rational target sets produced by the fixed map.
