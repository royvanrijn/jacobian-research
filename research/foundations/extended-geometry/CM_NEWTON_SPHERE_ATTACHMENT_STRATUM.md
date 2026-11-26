# Coherent sphere detection on a fixed boundary stratum

Let \(z=x/t^2\), and consider
\[
 f=t^{18}(z^2+1)^4/z+t^{17}B(z)+t^{16}(-z^{-1}+c+ez),
     \qquad \deg B\le5.
\]
Write \(B=b+b_1z+\cdots+b_5z^5\). This note specifies a connected
Zariski-open coefficient stratum on which sphere variations contain both
core coherent period rows. It proves coherent detection throughout that stratum and exhibits
points where each initial finite determinant vanishes.

**Theorem (coherent detection on the fixed-boundary stratum).**
At every complex specialization in the stratum \(S\) below,
the kernel of all marked sphere responses on
\(t^{12}(z^2+1)^2(a_0+a_1z)+t^{11}C(z)\), \(\deg C\le5\),
has \(a_0=a_1=0\). The primitive order-18 infinity
eigenspaces are in the sphere-variation span throughout \(S\).
For each algebraic specialization in \(S\), a finite set of
shifts detects the coherent input at all sufficiently large
primes \(p\equiv7\bmod18\).

## Definition of the stratum

Use \(h=1/t\), \(P=(z^2+1)^4\), \(V=zB\), \(W=-1+cz+ez^2\), and
\(\Phi=P+hV+h^2W\). The level equation is \(\Phi=\lambda zh^{18}\).
Let
\[
 (a_0,a_1,a_2)=(18P,17V,16W),\qquad
 (g_0,g_1,g_2)=(zP'-P,zV'-V,zW'-W).
\]
The derivative numerators are \(a=a_0+a_1h+a_2h^2\) and
\(g=g_0+g_1h+g_2h^2\). Put
\[
 \alpha=a_2g_1-g_2a_1,\quad \beta=a_2g_0-g_2a_0,\quad
 \gamma=a_1g_0-g_1a_0,\qquad
 R=\frac{\beta^2-\alpha\gamma}{(z^2+1)^3}.
\]
The numerator is the resultant in \(h\). This formula shows that \(R\)
is a polynomial of degree at most 16, with leading coefficient
\(-1827e b_5^2\) and constant term 4.

Define \(S\) by the following conditions.

1. \(b_5eB(i)B(-i)\ne0\), and both \(W\) and \(V^2-4PW\) are
   square-free.
2. \(R\) is square-free, and \(z(z^2+1)W\alpha\beta\) is a unit in the
   quotient by \(R\). In that quotient put \(H=-\beta/\alpha\).
   The determinant \(a_z g_h-a_h g_z\), evaluated at \(h=H\), and
   the critical value \(\mathcal L=\Phi(z,H)/(zH^{18})\) are units.
3. The power basis \(1,\mathcal L,\ldots,\mathcal L^{15}\) has nonzero
   determinant in the quotient by \(R\).

These are Zariski-open conditions, given by the nonvanishing of
polynomials after clearing denominators that are already units. For
each root of \(R\), eliminating \(h^2\) gives the unique critical
coordinate \(h=H\). Indeed the linear relation
\(a_2g-g_2a=\alpha h+\beta\) determines \(h\) and the resultant
vanishes, with \(a_2\) a unit. At \(z=\pm i\) the value of \(R\) is
\(-17\alpha V(7z^2-1)\), so the same linear relation there forces
\(h=0\). Thus there are exactly sixteen critical points in the torus.
The Hessian condition makes them Morse, and the power-basis condition
makes their nonzero values pairwise distinct.

The stratum is nonempty. At \(B=1+z^5,c=e=1\) the coefficient vector
of \(R\) in increasing degree is
\[
(4,64,635,4992,27418,-25152,55283,-87180,25052,
 -39976,-14940,16420,-2159,29616,18390,-2320,-1827).
\]
It reduces to an irreducible polynomial at 103. The above unit tests
hold modulo 103 and the critical-value power-basis determinant in the
monomial basis is 33. The two derivative identities are checked exactly
in \(\mathbb Q[z]/(R)\). The sixteen values of the seed are therefore a
single Galois orbit. As a nonempty Zariski-open subset of complex
affine eight-space, \(S\) is path-connected.

## The uniform relative family

Compactify the levels over \(\mathbb C^*_\lambda\) in
\(\mathbb P^1_z\times\mathbb P^1_h\). The charts of the normalization
at the ends, uniformly locally in \(S\times\mathbb C^*_\lambda\), are as
follows.

* At \(z=0,h=\infty\), put \(T=1/h\). Writing the equation as
  \(T^{18}P+T^{17}V+T^{16}W=\lambda z\) gives a single chart
  solving for \(z\), with derivative \(-\lambda\). This is the
  small-\(t=0\) puncture, which is filled.
* At \(z=h=\infty\), put \(w=1/z\) and \(T=1/h\). The equation is
  \[
  T^{18}p(w)+w^2 T^{17}v(w)+w^6 T^{16}d(w)=\lambda w^7,
  \]
  where \(p(0)=1,\ v(0)=b_5,\ d(0)=e\). The substitution
  \(w=u^{18}, T=u^7 s(u)\) gives eighteen simple initial roots
  \(s(0)^{18}=\lambda\). The choices are related by reparametrizing
  \(u\) by eighteenth roots of unity, since 7 is coprime to 18.
  Thus they give one normalized chart with uniform parameters.
* At \(z=\pm i,h=0\), solve the level equation for \(h\) by its
  nonzero \(V\)-derivative. The chart is \(h\asymp(z\mp i)^4\).
* At \(z=0,h=\pm1\) the two regular points are the affine cap
  sections. The affine phase has \(y\)-derivative \(-2\) there.
  There are no other incidences with the removed torus divisors.

The normalized compact regular fibres accordingly form a proper smooth
family locally over the regular-value complement of the sixteen graphs.
One can obtain it directly by adjoining the displayed end charts to
the smooth torus fibres. Each local completion map is finite onto its
chart of the projective level closure (for the exceptional chart this
is given by \(u^{18}=w\)). The punctures and the two ordered cap
sections are smooth, disjoint after normalization, and submersive
over the parameter base. Relative homology of the partially filled
fibres with the cap sections as boundary therefore gives local
systems \(0\to V_{\rm hom}\to E_{\rm hom}\to\mathbb C\to0\) on
the regular-value complement. This also follows by applying
[Thom's first isotopy lemma, Mather Proposition 11.1](https://webhomes.maths.ed.ac.uk/~v1ranick/surgery/mather.pdf)
to the proper smooth family with strata the marked points and their
complement. The fibres are connected. For sufficiently large
\(|\lambda|\) the angular core has all eighteen sheets, and a loop
in \(z\) about zero permutes them transitively. Every component
meets the core, since a fixed \(z\) does not give an identically
vanishing polynomial in \(h\).

For a compact parameter path a common exterior \(\lambda_*\) can be
chosen larger in modulus than both the critical values and the
maximum of \(|f|\) on the real sphere throughout the path. The
distinguished relative sphere chain varies continuously at
\(\lambda_*\). In the angular residue formula this is the simple
Rouché root of \(x(f-\lambda_*)\) on each latitude disk, continued
through the filled puncture and between the cap sections. Its relative
class \(e\) is consequently flat along the exterior section.

The varying punctured value-spheres are locally topologically
trivial over the parameter path. Locally one can move the distinct
critical values by isotopies in disjoint disks avoiding \(0,\infty\)
and the exterior section. These isotopies identify the fundamental
groups of the complements, and the local system identifies the
relative monodromy actions along the path. The sphere-variation span
\(\mathcal V=\operatorname{span}\{T_\gamma e-e\}\) is carried
isomorphically along such transport, as are the monodromies at zero
and infinity up to the induced conjugacies. Thus an inclusion of a
selected generalized infinity eigenspace in \(\mathcal V\) at one
parameter persists at every parameter of \(S\).

## Attachment at the seed

For the rational seed, the unmarked denominator witness at
\(p\equiv25\bmod36\) gives a nonrational scalar sphere resolvent,
exactly as in the two preceding conductor-18 pencils. The sphere is
fixed by the outer loop, the resolvent has rational Laurent coefficients,
and the same regularity theorem applies. All finite nonzero singular
fibres have one ordinary node and the cap sections avoid the node.
Thus Picard--Lefschetz variation at each of the sixteen values has a
one-dimensional image, so the [arithmetic attachment lemma](CM_NEWTON_SPHERE_ESCAPING_VARIATIONS.md#the-arithmetic-attachment-lemma)
applies.

At zero there is no primitive order 18 eigenvalue in absolute
homology. The following checks give the A'Campo computation at the
seed, and in fact hold throughout \(S\). The strict poles are
\(h=0,z=0,z=\infty\) of orders \(18,1,7\), and the strict zeros
are \(h=\infty\) and the residual divisor of orders 16 and 1.
The residual curve is smooth in the torus by the square-free
quadratic discriminant, and is transverse to \(h=\infty\) away
from \(z=\infty\) by square-freeness of \(W\).
At \(h=\infty,z=0\), the exponents are \((16,-1)\);
at \(z=0,h=\pm1\) they are \((1,-1)\).
At \(h=\infty,z=\infty\) the two zero branches have \(T\)-orders
2 and 4 in \(w\), with nonzero leading terms since \(b_5e\ne0\).
The charts \(T=w^jU\), \(j=1,2,3,4\), produce multiplicities
\(11,29,46,63\), and the residual mixed crossing \((11,-7)\)
adds \(4,-3,1,-2,-1,0\).
At each \(h=0,z=\pm i\) the smooth zero branch has contact
order four with the pole, giving negative orders
\(-17,-34,-51,-68\) and a terminal mixed crossing with the
simple zero branch. A mixed crossing is resolved by repeatedly
adding opposite-sign orders; its positive orders cannot exceed
its starting positive order. The same-sign crossings of strict
poles are already normal crossings. Hence each positive
multiplicity in the resolution is among \(1,\ldots,16,29,46,63\)
and none is divisible by 18. [A'Campo's product, Section 2 and Theorem 3](https://www.e-periodica.ch/cntmng?pid=com-001%3A1975%3A50%3A%3A27)
excludes all
primitive eighteenth roots from the compact \(H_1\).
The remaining open punctures are the two fixed finite root
points and the single \(z=h=\infty\) point. Their boundary loop
space has only eigenvalue 1. The exclusion therefore holds for
partially filled absolute homology as well.

The arithmetic lemma now puts the primitive order-18 infinity
eigenspaces in \(\mathcal V\) at the seed. This inclusion transports
along all paths in \(S\). This uses one Galois orbit only at the seed.
A second way of producing attachment in a fixed pencil is to attach
one node using nonrationality and apply \( (T_d-1)\delta_c =
\pm\langle\delta_c,\delta_d\rangle\delta_d\) across nonzero
intersection edges. Each step stays in the stable space
\(\mathcal V\). A connected finite intersection graph is sufficient
in place of the Galois-orbit assumption. The propagation is the
Picard--Lefschetz step in [Movasati--Nakai, Section 2, Proposition 1](https://w3.impa.br/~hossein/myarticles/namo-2006-1.pdf);
their further spanning assertion assumes a strongly tame morphism.

## Coherent detection on the stratum

Take \(A=(z^2+1)^2(a_0+a_1z)\) and let \(C\) be any correction
polynomial of degree at most five. Over a compact angular core of the
central cover \(Y^{18}=Q=(z^2+1)^4/z\), the continued core has the
exact cyclic action at infinity
\[
 \zeta T_j(z,\zeta^{-1}\varepsilon)=T_{j+1}(z,\varepsilon),
 \qquad \zeta=e^{2\pi i/18},\quad \varepsilon=\lambda^{-1/18}.
\]
Consequently, a primitive projection of a core cycle is represented
on the same core and belongs to the primitive infinity eigenspace of
the actual fibre. At the coherent order \(\lambda^{-5/18}\),
the marker \(t^{12}A+t^{11}C\) has relative differential proportional to
\(Y^{-13}A\,dz/z\) on the core. The correction starts one order
later in \(\varepsilon\). The central differential extends
holomorphically over the compact cover. For the puncture coordinates
\(z=w^{18}\), \(z=w^{-18}\), and \(z\mp i=w^9\), the smallest
orders of this differential are respectively 12, 0, and 0. A compact
central cycle in the required primitive character lifts to the angular
core. Indeed the kernel of the core-to-compact map consists of
puncture loops with deck eigenvalues \(1,-1\). Continue this lift on
the actual fibre and project to the same character using the exact
sheet action.

The adjacent-sheet differences of the positive and negative rays
in the central cover give two independent coherent period rows
\((I_0,I_1)\) and \((I_0,-I_1)\) up to nonzero sheet factors, where
\[
 I_j=\int_0^\infty x^{j-5/18}(1+x^2)^{-8/9}\,dx,\qquad j=0,1.
\]
These cycles can be projected before being lifted because the
coherent differential pairs only with its character \(\zeta^{13}\).
After continuation the projected cycles belong to \(\mathcal V\)
by the transported attachment. Thus vanishing of all marked moments
forces both \(a_0=a_1=0\) for every specialization in \(S\).
For an algebraic point of \(S\), choose a finite row subset with
the same span as the full response. Its rank difference between
the full matrix and the correction block is two. It survives at
all sufficiently large primes, since a witnessing nonzero minor
remains nonzero and each larger correction minor is identically
zero at the point. Entrywise Frobenius preserves these ranks.
The difficult prime class has a nonzero coherent input and at
most three correction inputs for this support. The other prime
class already has a nonzero shift-zero residue. This gives a
finite witness window for each point. A computable uniform bound
follows from the [effective cutoff](CM_NEWTON_FULL_WIDTH_STRIP_MOMENTS.md#effective-bounds-for-the-reconstruction-window).
The period argument gives no numerical bound of sixteen shifts.

The stratum meets both determinant loci with \(b b_1b_2ce\ne0\). On the line \(B=b+z+z^2+z^5,\ c=e=1\), the
five roots of the [coherent quintic](CM_NEWTON_FULL_WIDTH_STRIP_MOMENTS.md#a-quintic-for-the-three-pole-response)
\(F\) and the twelve roots of the [full determinant](CM_NEWTON_SPHERE_ESCAPING_VARIATIONS.md#a-generic-eight-row-window)
\(\Delta_8\), counted
within their respective polynomials, are nonzero and all belong to
\(S\). To check this, reduce the two determinant polynomials modulo
1009; their leading coefficients are units, they are square-free
with nonzero constants, and their factor degrees are \(1,1,1,2\)
and \(2,3,7\), respectively. In the residue field of each factor,
each test defining \(S\) holds. The unit leading coefficients make
every algebraic root integral at places above 1009. Its reduction is
one of the tested points, so all the nonzero conditions lift to
characteristic zero. This certifies geometric detection on actual
examples where each of the small determinants vanishes. This stratum
argument leaves open specializations where the stated boundary or
critical-value conditions fail.


## A coherent direction missed by rows zero through seven

Set \(b_1=b_2=c=e=1\) and keep \(b\) variable. With
\(a_0=0,a_1=1\), solve rows 1 through 3 for the correction
coefficients and then rows 4 through 6 for \(b_3,b_4,b_5\).
These are the triangular pivots of the three-pole system.
The checker evaluates the literal rows over \(\mathbb Q(b)\)
and checks the first seven equations after substitution.
The seventh-shift residual has denominator \(b^6\)
and numerator of degree twelve. Its monic numerator has
\(1009\)-integral coefficients. Modulo \(1009\) they are, in
increasing order,
\[
 (109,847,249,565,381,59,859,657,853,49,620,0,1).
\]
The polynomial has a simple root at \(b=519\) in
\(\mathbb F_{1009}\). Its characteristic-zero Hensel lift is
algebraic, and none of the denominators in the solved
parameters or inputs vanish there. The reduced \(B\)-coefficients
are \((519,1,1,831,287,707)\); the reduced inputs on
\((a_1,c_0,c_1,c_2)\) are \((1,55,471,645)\).

All the \(S\) tests hold for these parameter residues, with
critical-value power-basis determinant \(869\). They also
hold at the lift. Thus this algebraic phase in \(S\) admits
an exact coherent blind direction for rows \(0,\ldots,7\).
Its row 8 response on that direction is \(532\bmod1009\).
Since the early pivots leave just the one coherent line,
rows \(0,\ldots,8\) detect the coherent component at this
example.

The same fixed base \(b_1=b_2=c=e=1\) gives a uniform exact
certificate. After the triangular solves, the numerators of
row 7 and row 8 have a verified rational Bézout identity.
Hence rows \(0,\ldots,8\) detect the coherent component for
all \(b,b_3,b_4,b_5\) on this coefficient slice, with \(b=0\)
handled by row 1. The lifted example shows that the shorter
initial window fails on this slice.

## Exact replay

With SageMath 10.9, run

    sage -python research/foundations/extended-geometry/verify_cm_newton_sphere_attachment_stratum.py

The checker verifies the generic seed and the two determinant
intersections and the eight-row failure. The boundary and transport arguments are written
proofs.
