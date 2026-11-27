# The ring class direction on a cubic tangent line

A fixed cubic tangent line puts a rational root in
every reflected cofactor cubic. The corresponding
cube root in the maximal quadratic order also
determines a projective residue at the moving
conductor. In the tangent arithmetic progression,
that residue has order three at every conductor
prime. The two consecutive classes can therefore
give independent order classes on a family where
their maximal-order classes coincide.

## Fixed-profile tangent squares and the cofactor root

Use \(p=(p_1,p_2,p_3)\) for distinct odd
prime vertices, and fix the normalized coefficients
\((x,y,z)\) in
[the signed reconstruction](ABC_COHERENT_NORM_PENCILS.md#exact-reconstruction-from-three-signed-coefficients).
Put \(s=(y,z,x)\), \(Q=|xyz|\), and
\(D=-d\delta^2=f^2D_K<0\) for the quadratic
order discriminant and its conductor \(f\).
Assume \(D_K<-4\). The Heron curve in prime
coordinates is the square image of
\(\sum_i s_i X_i^3=0\). Its cubic in the
dual coordinates of the prime line is
\[
             F(n)=\frac{n_1^3}{y^2}
                   +\frac{n_2^3}{z^2}+\frac{n_3^3}{x^2}=0 .
\]
The weights may be taken positive in this geometric
description, as the signs can be absorbed in \(X_i\).

**Proposition TL1.** Let a fixed rational prime line
with these weights carry genuine triples with
arbitrarily large distinct prime vertices. If the
restricted Heron sextic has a repeated factor, the
prime line has a fixed rational normal \(n\in F\)
with all three coordinates nonzero. The reflected
cofactor cubic of each genuine triple on the line
has the rational root \(n\). Thus the consecutive
virtual cube classes coincide. In particular a
fixed prime line carrying triangles with independent
virtual cube classes cannot lower the degree of
its Heron form by polynomial square extraction.

A line carrying these primes has all its
coefficients nonzero, since a zero coefficient
fixes a ratio between two primes. By
[HT1](ABC_WEIGHTED_HERON_TANGENT_SQUARES.md#rational-tangency-and-the-incident-prime-condition),
an axis intersection would fail the incident-prime
unit test at the corresponding large prime.
The remaining way to obtain a repeated factor is
a smooth tangency at the square image of a rational
point \(X\). The tangent normal is proportional
to \((s_i X_i)_i\), which lies on \(F\).
Equivalently this is the Gauss map of \(F\),
whose gradients have coordinates proportional
to \(n_i^2/s_i^2\).
For every prime triple on the line,
the section of \(F\) by \(\sum_i p_i n_i=0\)
has the fixed point \(n\). By
[RF2](ABC_COHERENT_NORM_CUBIC_REFLECTION.md#line-sections-and-a-global-restriction-on-the-labels)
this is the reflected cubic, after changing
its binary line coordinates. The Cardano
criterion of RF1 then gives equality of the
consecutive virtual cube classes.

## Which conductor primes see the root

Use the
[cofactor lattice and ideal](ABC_COHERENT_NORM_COFACTOR_CUBIC.md#exact-cubic-and-quadratic-forms).
Scale the fixed \(n\) of TL1 to be a primitive
vector of the cofactor lattice along its line.
This scale is fixed, because the cofactor lattice
consists of the plane \(p\cdot n=0\) with the
coordinate divisibilities \(y\mid n_1,\ z\mid n_2,\ x\mid n_3\).
Let
\[
 I=\overline{\mathfrak p}_1\mathfrak p_2^2
       \overline{\mathfrak p}_3,\qquad
 \Theta=\alpha_{12}\overline{\alpha}_{23},\qquad
 L(p)=\frac{p_1 n_2^2/z^2-p_2 n_1^2/y^2}{n_3}.
\]

**Theorem TL2.** The norm Hessian of the
cofactor cubic on its rational root satisfies
\(H_{\mathcal B}(n)/9=L(p)^2\).
Let \(\beta\in\mathcal O_K\) give
\(\beta^3=\Theta\), with \((\beta)=I\mathcal O_K\).
Up to the choice of sign of \(L(p)\), its
doubled irrational coefficient is
\(\delta n_2/(zL(p))\). At a prime
\(q\nmid6xyz n_1n_2n_3\) the following hold.

- If \(v_q(L(p))>0\), then \(v_q(f)=v_q(L(p))\).
  The projective residue of \(\beta\) at
  \(q^{v_q(f)}\) has order three, and
  \[
    \mathcal O_K/q^{v_q(f)}\mathcal O_K
       \simeq \mathbb Z[\zeta_3]/q^{v_q(f)}\mathbb Z[\zeta_3] .
  \]
- If \(q\mid f\) and \(v_q(L(p))=0\), the
  projective residue of \(\beta\) at \(q^{v_q(f)}\)
  is trivial.

Write \(a_2=n_1/y,\ a_3=n_2/z\) and use the exact
cofactor formula \(\gamma=(p_1p_2^2 a_2+\alpha_{12}a_3)/x\).
Thus \(H_{\mathcal B}(n)/9=N\gamma/(p_1p_2^2p_3)\).
Insert \(u_{12}=(p_1^3 y^2+p_2^3 z^2-p_3^3 x^2)/(yz)\)
in the norm expression, and group the terms using
\(p_1n_1+p_2n_2=-p_3n_3\). This gives
\[
 \frac{H_{\mathcal B}(n)}9
 =-\frac{n_1n_2n_3}{Q^2}\sum_i\frac{s_i^2p_i^2}{n_i}.
\]
Eliminate \(p_3\) by the same linear relation.
The equation \(F(n)=0\) then gives
\(H_{\mathcal B}(n)/9=L(p)^2\).
Since \(n\) is a root, \(\gamma^3/\Theta\) is
rational. Its valuations in \(K\) are multiples
of three, using \(I^3=(\Theta)\), and ramification
indices in a quadratic field are prime to three.
Hence \(\gamma^3/\Theta=(\pm L(p))^3\), where
the absolute value is given by the norm equality.
In particular \(L(p)\ne0\). The doubled irrational
coefficient of \(\gamma\) is \(\delta n_2/z\).
The cube root \(\beta=\gamma/(\pm L(p))\)
is integral in \(\mathcal O_K\) and generates
the extended ideal by ideal factorization.

For the stated primes \(v_q(\delta)=v_q(f)\).
Integrality gives \(v_q(f)\ge v_q(L(p))\).
The residue of \(\beta\) is a unit at a conductor
prime since the vertex product is prime to \(f\).
Its projective cube is trivial there because
\(\Theta\in\mathcal O_f\).
If \(q\mid L(p)\), then its irrational
coefficient is nonzero modulo \(q^{v_q(f)}\).
For \(q\ne3\) the kernel of reduction of
projective units from \(q^a\) to \(q\) is a
\(q\)-group, by the filtration of principal
units. The residue is therefore non-scalar
already modulo \(q\). Its irrational coefficient
is a unit and \(v_q(f)=v_q(L(p))\).
If \(t=\operatorname{Tr}(\beta)\) and \(b=N\beta\),
then \(\beta^3=(t^2-b)\beta-tb\).
The basis \(1,\beta\) modulo \(q^{v_q(f)}\)
now gives \(t^2=b\) in the scalar ring, with
\(t\) a unit. The element \(-\beta/t\)
provides the claimed cyclotomic presentation.
For \(q\mid f\) with \(v_q(L(p))=0\)
the coefficient is already zero modulo the
full prime power, proving triviality.

For example, use \((x,y,z)=(3,1,2)\),
\((p_1,p_2,p_3)=(127,37,67)\) and \(n=(1,2,-3)\).
Then \(L(p)=-30\), the signed field radicand
is \(-44823723=-3\cdot7\cdot2134463\),
and \(f=\delta=30\).
The Heron incidence gcd is one and the
three reconstructed traces are
\(-227936,143516,-1517546\) with doubled irrational
coefficients \(90,30,60\). Direct multiplication gives
the coherent edge relation and the indicated cube root.
The cofactor cube root is
\(\beta=-(1331+\sqrt{-44823723})/2\),
with norm \(127\cdot37^2\cdot67\).
Modulo \(5\) it itself satisfies
\(\beta^2+\beta+1=0\) and is non-scalar.
This example uses unequal tangent weights.

## The arithmetic progression at its full conductor

In the normalized
[CT1 family](ABC_CUBIC_TANGENT_NORM_RANK_AND_CONDUCTORS.md#the-tangent-identity-and-the-lost-cubic-direction)
take \(h\) divisible by \(6\), with distinct
primes \(p-h,p,p+h>3\) and squarefree signed
\(R=3p^4-6p^2h^2-h^4\).
Put \(Y^2=-R\), \(K=\mathbb Q(Y)\), and
\(\mathcal O_h=\mathbb Z+h\mathcal O_K\).
Here \(D_K=-R\equiv5\bmod8\).
Work with ordinary Picard groups.
For the fixed oriented vertices, let
\(t_{12},t_{23}\in\operatorname{Pic}(\mathcal O_h)[3]\)
be the consecutive prime-quotient classes.

**Theorem TL3.** The element
\(\beta=(Y-p^2-h^2)/2\) gives an isomorphism
\[
 \mathbb Z[\zeta_3]/h\mathbb Z[\zeta_3]
       \longrightarrow \mathcal O_K/h\mathcal O_K,\qquad
 \zeta_3\longmapsto \beta/p^2 .
\]
In the projective residue kernel of the ring
class exact sequence, \(t_{12}-t_{23}\) is
represented by this cyclotomic class. Its local
projective order is three at every prime power
in \(h\). For imaginary \(K\) the difference
is nonzero in \(\operatorname{Pic}(\mathcal O_h)\).
If the common maximal-order edge class is
nonzero, the two order classes are independent.

By CT1 the exact norm generators belong to
\(\mathcal O_h\), their principal ideals have
norms prime to \(h\), and their quotient gives
\(\beta^3=\Theta\) with
\((\beta)=I\mathcal O_K\). The norm and trace of
\(\beta\) are \(p^2(p^2-h^2)\) and \(-p^2-h^2\).
Modulo \(h\), its characteristic equation
thus gives the displayed cyclotomic relation.
In the basis \(1,(1+Y)/2\) the coefficient
of the second basis element in \(\beta\)
is one. Since \(p\) is a unit modulo \(h\),
this proves that the map is an isomorphism.
In particular \(\beta/p^2\) is non-scalar even
modulo each prime dividing \(h\), including
\(2\) and the ramified prime \(3\).

Under the ordinary ring class kernel map
from [the residue sequence](ABC_COHERENT_NORM_RING_CLASS_RESIDUES.md#ideal-class-lifts-and-exact-residue-patterns)
the principal maximal-order ideal generated
by \(\beta\) gives precisely \(t_{12}-t_{23}\).
On the imaginary side the global units are
\(\{\pm1\}\), so the nonzero projective residue
survives. The two edge lifts have order dividing
three and map to the same nonzero class.
A dependence upstairs would first require the
sum of its coefficients to be zero downstairs
and would then kill their nonzero difference.
This proves independence.

In this progression the vectors
\((1,-2,1)\) and \((p,2h,-p)\) form a
cofactor basis, since their cross product is
\(2(p-h,p,p+h)\) and the cofactor plane
requires the second coordinate to be even.
In its coordinates \(S,T\), the binary cubic
and its primitive quadratic Hessian are
\[
 \begin{aligned}
 \mathcal B(S,T)
  &=-T\{3hS^2+3(p^2-h^2)ST+h^3T^2\},\\
 H_{\mathcal B}(S,T)/9
  &=h^2S^2+h(p^2-h^2)ST+p^2(p^2-2h^2)T^2 .
 \end{aligned}
\]
Their discriminants are \(27h^2R\) and \(-h^2R\),
respectively. This displays the rational root and
the nonprincipal norm-form class in the quadratic
order together on every imaginary example.

## Real unit image and samples of the two outcomes

**Corollary TL4.** On the real side of TL3,
the ring class difference vanishes precisely
when the diagonal cyclotomic residue belongs
to the cyclic image of a fundamental unit.
This occurs at \((p,h)=(17,12)\) and
\((199,138)\). At \((p,h)=(37,30)\) the
difference survives.

The criterion is the
[real unit quotient in the residue sequence](ABC_COHERENT_NORM_RING_CLASS_RESIDUES.md#real-quadratic-fields-and-the-common-unit-image).
Write \(\omega=(1+Y)/2\) and
\(\epsilon=u+v\omega\) for a unit.
In the cyclotomic basis of TL3 the residue
coordinates are
\[
 \left(u+v\frac{p^2+h^2+1}{2},\ vp^2\right)\pmod h .
\]
At \((17,12)\) one can use the norm-one unit
\(3289+47\omega\) of discriminant \(19869\).
Its residue is a scalar multiple of \(\zeta_3\).
The real discriminants at \((199,138)\) and
\((37,30)\) are \(182924997\) and \(2580117\).
The fundamental-unit coordinates \((u,v)\)
modulo the respective \(h\) are \((96,37)\)
and \((29,0)\), verified by
[PARI/GP's exact quadratic unit calculation](https://pari.math.u-bordeaux.fr/dochtml/ref/Arithmetic_functions.html)
using `quadunit(D_K)` in the integral basis.
The first gives cyclotomic coordinates \((91,91)\),
a scalar multiple of \(\zeta_3^2\).
The second is scalar in the whole conductor,
so no unit can identify its nonzero cyclotomic
class with the trivial class.

## Ring class rank in moving quadratic orders

**Corollary TL5.** For each fixed integer \(r\ge4\),
there are \(\gg_r T^{1/3}/(\log T)^4\)
distinct imaginary quadratic orders of
absolute discriminant at most \(T\), each
from the tangent prime progression and
with three-torsion rank at least \(r\).
Their two displayed prime-quotient classes
are independent in the order. The underlying
quadratic fields are distinct.

By the Chinese remainder theorem, the
projective residue group modulo \(h\)
has independent three-torsion at every
prime dividing \(h\). Each such factor
contains the local class of TL3, and the
imaginary global units are scalar.
On the family of CT4 the edge class in the
maximal order is nonzero by the norm bound
\(p_1p_2<|R|/4\) for sufficiently large \(|R|\).
Thus the three-torsion rank in the quadratic
order is at least \(\omega(h)+1\), where \(\omega(h)\)
is the number of distinct prime divisors.
In particular at \(h=6k\) with prime \(k>3\)
it is at least four.

For the stated \(r\), fix a squarefree product
\(M\) of \(r-4\) distinct primes at least \(13\),
and use \(h=6M k\) with prime \(k>3\) outside \(M\).
The squarefree tangent sieve of CT4 applies
to the four nonproportional prime forms
\(p,p-6Mk,p+6Mk,k\).
At each prime of \(M\) take \(p,k\) units.
The quartic there reduces to \(3p^4\), so
this causes no fixed square divisor.
The conditions at \(2,3,5,7,11\) and the
remaining primes are the ones of CT4 after
rescaling \(k\) at the unit primes.
The binary quartic remains irreducible after
this fixed invertible rational rescaling, so
the squarefree and essential-representation
steps are unchanged. This gives
\(\gg_M X^{1/2}/(\log X)^4\) distinct
\(|R|\in[X,2X)\) with \(h\asymp_M X^{1/4}\).
Their order discriminants have absolute
value \(|R|h^2\asymp_M X^{3/2}\), which
gives the claimed count and the rank.

This rank comes from the moving conductor.
For the maximal-order torsion moments,
the tangent families still supply only
one class direction. TL1 identifies the
obstruction to using a fixed prime line
with polynomial square extraction for
two independent maximal-order directions.
