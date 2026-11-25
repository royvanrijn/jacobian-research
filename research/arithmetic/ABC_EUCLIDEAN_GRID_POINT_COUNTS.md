# Bounded-height counts on fixed abelian varieties over varying prime-power extensions

Let \(A/K\) be a fixed abelian variety with a fixed symmetric ample
line bundle and its absolute logarithmic canonical height \(\widehat h\).
Fix a prime \(\ell\) and an integer \(d\ge1\). For every field \(F/K\)
contained in a Galois extension of \(K\) of \(\ell\)-power degree at most
\(d\), and every \(B\ge2\), one has
\[
 \#\{P\in A(F):\widehat h(P)\le B\}
       \le \exp\!\left(O_{A,K,\ell,d}\!\left(\frac{B}{\log(2+B)}\right)\right).
 \tag{EP1}
\]
The implied constant is uniform in the ramified primes and the field
discriminant. In a fixed model this gives a uniform subpower bound in
the exponential height. The geometric input to the proof is that an
algebraic point of small height cannot generate a field of large
discriminant. A small canonical-height ball in a large-discriminant
field therefore reduces, by translation, to points in proper subfields.

## A hereditary rank estimate

Enlarge the fixed base \(K\) to contain \(\mu_\ell\) and \(A[\ell]\),
and replace \(F\) by its compositum with this base. The Galois group
of the enlarged containing extension is a subgroup of the original
\(\ell\)-group, so the same degree cap works. From now on the base
contains these groups. Write \(D_E=|\operatorname{disc}(E)|\) for
any intermediate field \(K\subseteq E\subseteq F\). Every subgroup
of a finite \(\ell\)-group is subnormal, with a chain of steps of
index \(\ell\). Thus \(E/K\) is reached by a bounded number of cyclic
degree-\(\ell\) steps, even if \(E/K\) is nonnormal. Every rational
prime ramified in a step divides \(D_E\), by the discriminant tower
formula. Genus theory and the fixed base class group give
\[
 \log|\mathrm{Cl}(E)[\ell]|\ll_{K,\ell,d}
                          1+\omega(D_E).
 \tag{EP2}
\]
For each step one uses the bound for the \(\ell\)-rank of the upper
class group in terms of the lower \(\ell\)-rank and the number of
ramified primes. The proof by invariant classes also gives (EP2).
[Class group bound and cyclic-step proof](ABC_PRIMITIVE_NORM_RAY_PACKING.md#prime-power-towers-remove-the-class-number-cost)
· [Rosen, as stated and used in Dujella, Lemma 4.5](https://arxiv.org/html/2312.03655#S4)

The multiplication-by-\(\ell\) Kummer map injects \(A(E)/\ell A(E)\)
into \((E^\times/E^{\times\ell})^{2\dim A}\). Its image is unramified
away from \(\ell\) and the fixed bad reduction primes of \(A\):
multiplication by \(\ell\) is finite étale on the good abelian scheme
and rational points extend over the local valuation ring. The
unramified power classes have valuations divisible by \(\ell\)
outside these places. The group of these classes has size bounded
by the product of the corresponding unit-group quotient and the
\(\ell\)-torsion of the localized class group, which is bounded by
\(|\mathrm{Cl}(E)[\ell]|\). The number of places in the exceptional
set and the ordinary unit rank are bounded by the degree cap and the
fixed data. The resulting rank bound is
\[
 r(E)\ll_{A,K,\ell,d}1+\omega(D_E)
      \ll_{A,K,\ell,d}\frac{\log(3D_E)}{\log\log(10D_E)}.
 \tag{EP3}
\]
In particular this is a bound in terms of each subfield's own
discriminant.

## Point fields and covering by canonical-height balls

Fix a projective embedding of \(A\) over \(K\) coming from a very
ample symmetric line bundle that is a power of the chosen one.
The absolute logarithmic projective height satisfies
\(h(P)\ll_A 1+\widehat h(P)\).
For a point \(P\) in a field of bounded degree, the elementary
primitive-element discriminant estimate gives
\[
                  \log D_{K(P)}\ll_{A,K,d}1+\widehat h(P).
 \tag{EP4}
\]
Here \(K(P)\) denotes the field generated over \(K\) by coordinate
ratios of \(P\). To verify (EP4), use a chart with one nonzero
coordinate and adjoin its coordinate ratios to a fixed generator of
\(K/\mathbb Q\). Each ratio has absolute height at most \(h(P)\).
Among integer linear combinations of this bounded list with
coefficients in a fixed finite box, one separates every pair of
distinct embeddings of \(K(P)\). For each pair, equality defines
a proper linear hyperplane, occupying at most one fraction
\(1/(m+1)\) of a box with sides \(\{0,\ldots,m\}\). A box longer
than the number of pairs therefore supplies a generator \(\alpha\)
with \(h(\alpha)\ll 1+\widehat h(P)\). If \(a\) is the leading
coefficient of its primitive integer minimal polynomial, then
\(\beta=a\alpha\) is integral and \(h(\beta)\ll_d h(\alpha)\)
by the formula for the height using the leading coefficient
and conjugates. The absolute discriminant of \(K(P)\) divides
that of the order \(\mathbb Z[\beta]\). Bound the latter by
the product of squared differences of conjugates, using
\(|z-w|\le2\max(1,|z|)\max(1,|w|)\). This proves (EP4).
The discriminant-height inequality is classical.
[Silverman, *Lower bounds for height functions*](https://doi.org/10.1215/S0012-7094-84-05118-4)

Northcott finiteness at bounded degree, together with the fixed
height comparison, gives uniformly bounded torsion and a uniform
positive lower bound for the canonical heights of nontorsion
points. The quotient of \(A(E)\) by its torsion is a Euclidean
lattice of rank \(r(E)\) with squared norm \(\widehat h\).

Let \(J=\log D_E\) and suppose \(J\) exceeds a fixed constant.
Cover the radius-\(\sqrt B\) ball in the real span of this lattice
by balls of radius \(c\sqrt J\), for a fixed sufficiently small
\(c>0\). The number of balls is at most
\[
                         (1+C\sqrt{B/J})^{r(E)}.
 \tag{EP5}
\]
This follows by choosing a maximal separated set and comparing
volumes. If the whole ball fits in one of the smaller balls, use
one. In a nonempty ball choose an actual point \(P_0\) lying in it.
Every point \(P\) of the count in that ball satisfies
\(\widehat h(P-P_0)\le4c^2J\) and \(\widehat h(P-P_0)\le4B\).
By (EP4), with \(c\) small and the fixed threshold on \(J\)
large, \(K(P-P_0)\) is a proper subfield of \(E\) containing \(K\).
This deduction applies to every difference, including torsion.
Translation is injective. There are at most a bounded number of
intermediate fields, since the normal containing extension has
bounded degree. Thus each ball is bounded by the sum of the
radius-\(2\sqrt B\) counts over the proper subfields.

The logarithm of (EP5) is \(O(B/\log(2+B))\). If
\(J\le\sqrt B\), (EP3) gives a bound \(O(\sqrt B\log(2+B))\),
which is of the claimed size after increasing the constant. If
\(\sqrt B<J\le C_0 B\), then \(\log J\gg\log(2+B)\) and
\[
 J\log(1+C\sqrt{B/J})
   = B\, (J/B)\log(1+C(J/B)^{-1/2}) \ll_{C,C_0} B .
\]
Take \(C_0\) large enough that \(J>C_0 B\) uses just one ball.
For bounded \(J\), (EP3) and the uniform nontorsion spacing give a
polynomial bound in \(1+B\) by a volume comparison and the torsion
bound. The same last case includes the fixed base. Induct now on
\([E:K]\), with \(B\) replaced by at most \(4B\) at each step.
The number of steps and the intermediate-field counts are bounded,
so their exponent costs stay within \(O(B/\log(2+B))\).
This proves (EP1).

The covering and subfield argument applies whenever the same
hereditary rank bound (EP3) holds on a family of bounded-degree
extensions. A uniform hereditary bound \(r(E)=o(\log D_E)\)
would give the corresponding \(\exp(o(B))\) conclusion. To see
this, use the same split with \(\sqrt B\), the uniform
\(O(\log(3D_E))\) bound in the smaller-discriminant range, and the
uniform little-oh in the larger range. This version requires a
rank statement on all intermediate fields, including the fields
where a difference can be defined.

## Primitive quotients on two Euclidean signatures

For the signatures \((3,3,3)\) and \((2,4,4)\), in any orientation,
fix a squarefree integer \(S\). Count positive primitive triples
\(a+b=c\le X\) for which every prime outside \(S\) has valuation
on its own letter divisible by the corresponding order. Then,
uniformly for \(S\) and \(X\ge2\),
\[
 \#\{(a,b,c)\}\le
 \exp\!\left(O\!\left(1+\omega(S)+\frac{\log X}{\log\log(10X)}\right)\right).
 \tag{EP6}
\]
The same conclusions apply to multiples of the three corresponding
orders. A permutation uses one of the six fractional linear
transformations permuting \(0,1,\infty\).

For the displayed orientations let \(t=a/c\) and use the smooth
projective models of the covers
\[
                z^3=t(1-t),\qquad z^4=t^2(1-t).
\]
Their ramification orders at \(0,1,\infty\) are the indicated triples,
by the orders of the right-hand sides and the Kummer ramification
formula. The covers are geometrically connected because each
right-hand side has a simple zero, and Riemann–Hurwitz gives
genus one. Over a fixed number field containing an algebraic
basepoint each model is a fixed elliptic curve. Include the
relevant roots of unity in that field.

At a primitive triple, the support valuations of the rational
right-hand side outside \(S\) are all divisible by \(m=3\) or
\(m=4\), respectively. For the second cover the right-hand side
is \(a^2b/c^3\), explaining the order two on the first letter.
Write the right-hand side as \(u^m v\), where \(u\) is nonzero
rational and \(v\) is a signed integer with all exponents in
\(\{0,\ldots,m-1\}\) and support in \(S\). Signs allow the
permuted charts as well. There are at most \(2m^{\omega(S)}\)
choices of \(v\), and \(z=u\sqrt[m]{v}\) lies in a cyclic
extension of the fixed base of degree dividing \(m\). For the
lift, \(h(t)\ll 1+\log X\) and
\(h(z)=h(t(1-t))/3\) or \(h(t^2(1-t))/4\), respectively.
On the smooth models a point over the open set
\(t\notin\{0,1,\infty\}\) thus has canonical height
\(O(1+\log X)\). Indeed the change to any fixed projective
model is given by a finite set of rational maps on this smooth
open set, and the height bounds for their fixed coordinate polynomials
give the asserted bound. Each distinct primitive triple
determines a distinct \(t\). Apply (EP1) over each field,
and sum over the choices of \(v\), to prove (EP6).

In particular, if the support of valuations defective for one
of these signatures has radical at most \(X^{f\lambda}\), summing
(EP6) over the squarefree supports gives
\[
        \#\{(a,b,c)\}\ll_\varepsilon X^{f\lambda+\varepsilon}
                       \qquad(0\le f\lambda\le3).
 \tag{EP7}
\]
The constant here is uniform in the defect cap. With \(f<1/3\)
by a fixed margin this improves the small-budget coefficient.
For instance, if \(3\mid v_p\) on at least three quarters of the
logarithmic radical in the three letters, the count in the
\(abc\) budget \(\operatorname{rad}(abc)<c^\lambda\) is at most
\(X^{\lambda/4+\varepsilon}\). This is a discrete grid of fixed
spacing, compatible with diffuse limiting measures. On the
balanced minimally repaired quotient sector,
\(\operatorname{rad}(c)\ge c^{(1/3-\xi)\lambda}\) and
\(H_M(c)\le Y\) give
\(c\le Y^{1/(M(1/3-\xi)\lambda)}\). The same sector bound is
\[
             O_{\lambda,M,\xi,\varepsilon}
             (Y^{\,f/(M(1/3-\xi))+\varepsilon}).
\]
It concerns primitive quotients, so the argument introduces no
sum over common multipliers.

## Scope of the proof

The earlier hyperbolic prime-primary argument uses uniform
Mordell–Lang for a curve of genus at least two and obtains a
height-free estimate. For genus one the entire curve is an abelian
variety after choosing an origin. A positive-rank curve can have
infinitely many points, so the bounded-height qualifier is essential.
[Hyperbolic prime-power specialization](ABC_MISMATCHED_GRID_NUMBER_FIELDS.md#counting-the-specialization-fields-of-the-tower)

The two Euclidean counting consequences overlap the existing
literature. For \((3,3,3)\), Heath-Brown's coefficient-uniform
diagonal cubic estimate already implies the subpower bound per
weighted coefficient choice. Dujella's fixed-base-field theorem
also gives the stated quantitative bound on both signatures. For
a twist \(v u^3=t(1-t)\), the change
\(x=-4vu,\ y=4v(1-2t)\) gives \(y^2=x^3+16v^2\), with rational
3-torsion at \(x=0\). For a weighted \((2,4,4)\) relation divide by
a fourth power and rearrange signs to obtain \(w^2=\alpha u^4+\beta\).
Then \(x=\alpha u^2,\ y=\alpha u w\) gives
\(y^2=x^3+\alpha\beta x\) with rational 2-torsion. The map has
bounded fibres away from zero and the model has polynomial
coefficient height. One may restrict \(S\) to primes \(\le X\)
that occur in the actual defect record. More precisely each
such record for a triple has product at most \(X^3\), so the
finite twist labels or coefficient tuples used in these
comparisons have polynomial height. This supplies Dujella's
model-height threshold after enlarging the height parameter
by a fixed power. The contribution of (EP1) is a uniform
fixed-variety theorem over moving fields, with class torsion
allowed to grow, and its reduction to a bound on ranks over
all intermediate fields. The comparison with fixed-base
counting results is separate from that general theorem.
[Heath-Brown, Theorem 1.3](https://arxiv.org/html/2601.07817v1)
· [Dujella, Theorems 1.1 and 1.2](https://arxiv.org/html/2312.03655)
· [Related bounded-height point counts, Section 3](https://arxiv.org/html/2311.11266#S3)

The remaining Euclidean signature \((2,3,6)\) has the analogous
cover \(z^6=t^3(1-t)^2\), but its twists may require both a
moving quadratic and a moving cubic field. The present
rank argument does not bound the class torsion from a cyclic
step of unequal prime degree. A uniform hereditary
\(o(\log D_E)\) rank bound for this fixed elliptic curve over
the intermediate fields of the sextic extensions would suffice
for the same counting conclusion. Average bounds on new rank
in root extensions do not supply this pointwise hypothesis
on arbitrary support-restricted representatives.
[Shnidman and Weiss](https://arxiv.org/abs/2112.12864)

All the orders in (EP6) are fixed before the height limit. On the
formal balanced profiles supported on increasingly large prime
valuation degrees, these orders recover no positive limiting
radical mass. The result therefore leaves that part of the
unrestricted endpoint open.

For fixed coefficients, Achenjang, Arango-Piñeros and
Vemulapalli obtain a precise asymptotic for distinct primitive
images for every Euclidean signature, including \((2,3,6)\).
For a fixed equation \(F\) with points and positive Jacobian
rank \(r(F)\), the main term is a constant times
\((\log H)^{r(F)/2}\) at Weil height \(H\), and the primitive
points form finitely many cosets of a finite-index subgroup.
The rank, regulator, coset data and error constants depend on
\(F\). This does not supply a uniform rank or point bound
when the coefficients move with \(X\), as required above.
[Achenjang–Arango-Piñeros–Vemulapalli, Theorems 1.1 and 1.4](https://arxiv.org/html/2611.31699v1#S1)

## A parity count using a fourth determinant

A direct Legendre twist gives a complementary estimate using
Dujella's theorem over \(\mathbb Q\). For a triple \(a+b=c\le X\)
and a reduced signed rational \(r=a_0/c_0\) with \(c_0>0\), put
\(D_0=ac_0-a_0c\). Fix \(0<\lambda<1\) and nonnegative
\(f,g,\tau\le1\). Suppose
\[
 \operatorname{sf}_2(abc)\le X^{f\lambda},\qquad
 r\notin\{0,1,a/c\},\qquad H(r)\le X^{\tau\lambda},\qquad
 \operatorname{sf}_2(|D_0|)\le X^{g\lambda}.
\]
Here \(\operatorname{sf}_2\) selects the odd valuation primes.
Count primitive triples for which at least one such \(r\) exists.
Their number is
\[
                       O_\varepsilon(X^{(f+g+3\tau)\lambda+\varepsilon}).
 \tag{EP8}
\]
This permits imposing the original three-letter radical budget as
well. A uniformly bounded list of slopes costs at most
\(O_\varepsilon(X^{(f+g+\tau)\lambda+\varepsilon})\) if their
heights obey the same cap.

Indeed \(t(t-1)(t-r)=-abD_0/(c^3c_0)\) at \(t=a/c\). Its
square class has a signed squarefree representative \(D\) with
\(|D|\le X^{(f+g+\tau)\lambda}\), by multiplying the numerator
and denominator and comparing square classes. Choose rational
\(w\) such that \(D w^2=t(t-1)(t-r)\). For every fixed \(r,D\),
the point \(u=Dt,\ v=D^2w\) is on
\[
                         v^2=u(u-D)(u-Dr).
\]
This nonsingular elliptic curve has rational 2-torsion. The
curve coefficient heights and the point's \(u\)-height are
bounded by fixed powers of \(X\), uniformly over the stated
parameters. Dujella's Corollary 1.3 over the fixed field
\(\mathbb Q\), with the height cap raised to a fixed power,
gives \(O_\varepsilon(X^\varepsilon)\) points per \(r,D\).
Every \(u\) recovers \(t\) and its primitive triple. Enumerate
the \(O(X^{2\tau\lambda})\) reduced slopes and the
\(O(X^{(f+g+\tau)\lambda})\) square classes for each.
If \(0<|D_0|\le X^{\tau\lambda}\), one can take \(g=\tau\).
For fixed \(f+g+3\tau<1/3\) this yields a strict coefficient
saving. The estimate requires the actual fourth determinant
condition and does not provide a small-height slope satisfying
it at an unrestricted endpoint.
[Dujella, Corollary 1.3](https://arxiv.org/html/2312.03655)
