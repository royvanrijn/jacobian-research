# De Rham obstruction for leading Keller pairs

**Theorem KellerLeadingPairDeRhamObstruction.**
Let \(k\) have characteristic zero and let \(A,B\in k[x,y,z]\)
be algebraically independent. Set \(K=k(A,B)\), and let \(U\)
be a geometric component of the normalized affine generic fiber
of \((A,B)\), with smooth projective completion \(X\).
All differentials below are over the geometric generic base field.
The relative Jacobian differential \(\omega\) on this curve is
defined by
\[
          dA\wedge dB\wedge\omega=dx\wedge dy\wedge dz.
\]
For each point \(p\) of \(X\) let \(m_p\) be the pole order of
\(\omega\), or zero if \(\omega\) is regular there. Put
\[
           \Delta=\sum_p m_p p,\qquad
           D=\sum_p\max(m_p-1,0)p .
\]
Then \(\omega\) is exact on the function field of \(U\) if and
only if its class vanishes in \(H^1_{\mathrm{dR}}(U)\), or
equivalently if it belongs to the image of the finite-dimensional map
\[
        d:H^0(X,\mathcal O_X(D))\longrightarrow
                           H^0(X,\Omega_X^1(\Delta)).
\]
If the class is nonzero on even one geometric component,
\(A,B\) cannot be an ordinary or positive-weight leading pair
of a Keller map.

### The affine and complete curves

The generic fiber is smooth. Indeed on any irreducible
component of the rank-deficient locus of \(dA,dB\), the
restrictions of \(A,B\) are algebraically dependent in
characteristic zero, by the differential criterion for
transcendence degree. That component does not dominate the
target plane. Thus on the geometric generic fiber one can
use the local quotient by an invertible leading-pair minor
to define \(\omega\). It is regular and nowhere zero on
\(U\). The same formula in the source field characterizes it
by
\[
           df=\det J(A,B,f)\,\omega
\]
for any rational \(f\), with the derivative relative to
the chosen fiber.

In characteristic zero a pole of a rational function of
order \(n>0\) makes its differential have pole order \(n+1\).
Thus any rational primitive of \(\omega\) is regular on \(U\),
and its poles on \(X\) are bounded by \(D\).
On a smooth affine curve the degree-one de Rham group is
global regular differentials modulo differentials of global
regular functions, by affine vanishing for the de Rham
complex. This proves both equivalences. For a proposed
Keller completion, the
[leading-pair rational slice](KELLER_LEADING_PAIR_RATIONAL_SLICE.md)
gives \(df=\kappa\omega\) for a nonzero scalar \(\kappa\).
Restriction to each geometric component preserves this
identity, proving the obstruction.

On the complete curve the distinction is the boundary
residue. If any residue of \(\omega\) is nonzero, the
class on \(U\) is nonzero without further calculation.
When all residues vanish, \(\omega\) is a differential
of the second kind and has a class on the complete curve
\(X\), with the same exactness test. More precisely, for
the reduced boundary \(D_0=X\setminus U\) the sequence
\[
 0\longrightarrow H^1_{\mathrm{dR}}(X)
 \longrightarrow H^1_{\mathrm{dR}}(U)
 \xrightarrow{\mathrm{Res}}
 \{(r_p)_{p\in D_0}:\ \sum r_p=0\}
 \longrightarrow0
\]
separates the two contributions. For example, the quotient
of the meromorphic de Rham complex allowing poles on \(D_0\)
by the regular complex has only the local residue
cohomology in degree one, because nonlogarithmic
principal parts integrate locally. The meromorphic complex
computes the de Rham cohomology of \(U\): its sheaves have
no higher cohomology by allowing arbitrarily large poles
on the nonempty boundary and applying affine vanishing
or Riemann–Roch. The residue map lands in the zero-sum
space by the residue theorem and is onto by Riemann–Roch
for \(\Omega_X^1(D_0)\).

This gives a test using local principal parts and regular
forms on \(X\). If the residues vanish, integrate the polar
part of \(\omega\) formally at each pole, obtaining a
principal part \(f_p^{\mathrm{pp}}\). The parts are realized
by a single rational function \(f\) with poles bounded by
\(D\) exactly when
\[
                   \sum_p\operatorname{Res}_p
                         (f_p^{\mathrm{pp}}\eta)=0
                  \quad\text{for every }\eta\in H^0(X,\Omega_X^1).
\]
Indeed the connecting map for
\(\mathcal O_X\subset\mathcal O_X(D)\) is the
Mittag–Leffler obstruction, and this formula is its
Serre-duality pairing. If it vanishes, the remainder
\(\omega-df\) is a global regular form, independent of
\(f\) up to constants in the primitive. The differential
is exact precisely when this remainder is zero. In
particular a nonzero regular form on \(X\) is an
obstruction, but the test also treats forms with
higher-order poles and zero residues. With a curve
model and bases of the two Riemann–Roch spaces, the
same test is a single linear system.

### Transfer through finite covers and Hamiltonian derivations

Exactness can be tested after any finite extension of
the curve function field in characteristic zero.
Pullback preserves exactness. Conversely the trace on
differentials commutes with \(d\) and multiplies a pulled
back differential by the extension degree. One therefore
descends a primitive by taking the normalized trace.
This applies to the auxiliary cyclic covers on which a
relative differential has a simpler expression.

More generally let \(\delta\) be any nonzero derivation
of a curve function field over its constant base field.
There is a unique rational differential \(\omega_\delta\)
with \(\omega_\delta(\delta)=1\).
Then \(\delta f=h\) is equivalent to \(df=h\omega_\delta\).
For example on a Poisson surface take
\(\delta=\{a,\cdot\}\) on the generic level of \(a\).
The same de Rham test applies to a proposed Hamiltonian
slice \(\{a,f\}=1\), using an affine open on which the
dual differential is regular. For a scalar commutator pair in the first Weyl algebra,
the [weighted Hamiltonian slice theorem](../../dixmier/dixmier-formal-lifting/POLYNOMIAL_WEYL_LEADING_FACTOR_RESIDUE.md#weighted-poisson-symbols-of-a-scalar-commutator-pair)
produces the rational slice for every nonconstant
leading symbol at integer weights of positive sum.
It derives the equation from the full commutator,
using commuting lifts of the algebraic constants.
For a weighted homogeneous symbol, the same theorem
and its corollary give the cyclic curve and differential
on which to run the exactness test.

## A cyclic-cover calculation

**Theorem KellerCyclicCurveDifferentialRemainder.**
Over a field \(K\) of characteristic zero containing the
\(m\)-th roots of unity, let \(C\) be the smooth affine
curve \(u^m=f(t)\), where \(m\ge2\) and \(f\) is squarefree
of positive degree \(d\). For \(1\le r<m\) and \(Q\in K[t]\),
the differential \(Q(t)\,dt/u^r\) has a rational
primitive if and only if the remainder \(Q_{\mathrm{red}}\)
obtained by reducing \(Q\) by the image of
\[
             \mathcal D_r(L)=fL'+\frac{m-r}{m}f'L
                 \qquad(L\in K[t])
\]
is zero. Here the unique remainder has degree at most
\(d-2\), obtained by cancelling highest coefficients.
The corresponding character of \(H^1_{\mathrm{dR}}(C)\)
has basis \(t^i dt/u^r\), \(0\le i\le d-2\).

The simple roots of \(f\) make the cover irreducible and
the affine curve smooth. Under \(u\mapsto\zeta u\) the
form has character \(\zeta^{-r}\). A rational primitive
must be regular on the affine curve by its pole orders.
Projecting it to that character gives \(u^{m-r}L(t)\)
with \(L\in K[t]\). Differentiating yields
\[
                        d(u^{m-r}L)=\mathcal D_r(L)\,dt/u^r .
\]
For \(L\ne0\) of degree \(n\), the output degree is
\(n+d-1\), since the multiplier of the leading
coefficient, apart from that of \(f\), is
\(n+(m-r)d/m\), which is nonzero. This proves the
unique reduction and exactness test. Every regular
form of that character is represented by a polynomial
times \(dt/u^r\): the two character generators from
\(dt,du\) contribute the factors \(f,f'/m\), which
are coprime. This proves the basis assertion. The
polynomial reduction is the characteristic-zero
form of the differential reduction used in
[Kedlaya's hyperelliptic algorithm](https://arxiv.org/pdf/math/0105031)
and the
[Gaudry–Gürel superelliptic extension](https://www.iacr.org/archive/asiacrypt2001/22480482.pdf).

For a rational input form or a nonsquarefree model,
one can instead use the pole bound on the smooth
normalization and project the Riemann–Roch space to
the same character. For \(u^m=f(t)\), the character
\(\zeta^a\) reduces the primitive equation to
\[
                      L'+\frac{a}{m}\frac{f'}{f}L=Q
\]
for a differential \(u^a Q(t)\,dt\).
The normalized valuation at a point over \(t_0\)
has ramification index \(m/\gcd(m,\operatorname{ord}_{t_0}f)\)
when the cyclic cover has full degree. This converts
the allowable pole orders of \(u^a L\) into bounds on
\(L\) on the line. Thus one can test poles by valuations
on the cover without a separate Puiseux series for every
geometric point.

For example, \(u^2=t^3-1\) and \(t\,dt/u\) give a nonzero
class with zero residue, despite the pole at infinity.
The polynomial has degree one, already below the
reduction threshold. By contrast \(t^2\,dt/u=(2/3)\,du\)
is exact. This illustrates why a positive genus or a pole
alone cannot decide the slice obstruction.

For a Keller example, the [nonsquare axial exclusion](SEXTIC_KELLER_TANGENT_CONIC_QUADRATIC_OUTER.md#the-nonsquare-case) uses
\(y^4=S/[X(X-d)]\), \(d\ne0\), and the form \(y^3dX/X\).
With \(u=1/y\), the identity
\[
 d\bigl((X-d)/u^3\bigr)
             =\left(-\tfrac12+\frac{3d}{4X}\right)
                                      \frac{dX}{u^3}
\]
reduces its class to \( (2/(3d))\, dX/u^3\), the nonzero
degree-zero remainder for the squarefree quadratic
\(u^4=X(X-d)/S\). For the [mixed-frame fiber](SEXTIC_KELLER_TANGENT_CONIC_QUADRATIC_OUTER.md#constant-projective-ratio) with
\(\omega=dx/(xyz)\), the earlier normalization
calculation instead gives \(\Delta=0\), so the first
finite-space test has zero image and excludes it
directly.
