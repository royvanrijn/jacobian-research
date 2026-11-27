# Coherent three-prime paths in a product of elliptic curves

A triple of coherent cube norms with fixed quadratic field
and fixed edge cube classes lies on a surface in a product
of three \(j=0\) elliptic curves. Its positive-dimensional
abelian translates occur along the boundary, so uniform
Mordell–Lang bounds the number of coherent triples in a
fixed collection of twists without a height restriction.

## The multiplicative surface and its boundary

Fix \(a_i\ne0\) in characteristic zero, and write
\(E_i:y^2=x^3+a_i^2\). Over an algebraic closure put
\[
        t_i(P_i)=\frac{y(P_i)+a_i}{y(P_i)-a_i},
      \qquad \mathcal T:\quad t_1(P_1)t_2(P_2)=t_3(P_3).
 \tag{1}
\]
The equation means the inverse image of the closed
multiplication graph in \((\mathbf P^1)^3\).
Let \(\mathcal U\) be its open set where all points
are finite and all \(x(P_i)\ne0\).

**Theorem ML1.** The surface \(\mathcal T\) is geometrically
irreducible and has degree \(18\) for the product
polarization by the three origins. The open set
\(\mathcal U\) meets no positive-dimensional translate of
an abelian subvariety contained in \(\mathcal T\).
A smooth geometric model has ample canonical divisor,
\(K^2=15\), \(p_g=4\), and irregularity \(3\).
Its curves of geometric genus at most one consist of
sixteen smooth elliptic curves, all mapping to the
boundary of \(\mathcal T\).

Scaling the \(y_i\) by \(a_i\) and the \(x_i\) by a cube
root of \(a_i^2\) reduces the geometric statements to
the identical curves \(E:y^2=x^3+A\) with one common
square root \(a=\sqrt A\). Over the open set of the first
two factors where \(x_1x_2(y_1+y_2)\ne0\), the triple is
equivalent to
\[
 s^3=y_1+y_2,\qquad x_3=x_1x_2/s^2,\qquad
                     y_3=(y_1y_2+A)/s^3 .
 \tag{2}
\]
These identities follow by clearing the denominators
in (1) and using the curve equations. Conversely
\(s=s^3/(s^2)\) is recovered from (2), since
\(s^2=x_1x_2/x_3\). The cubic in \(s\) is irreducible
over the function field of \(E^2\), as \(y_1+y_2\)
has a simple zero divisor along each graph below.
Each \(t_i:E_i\to\mathbf P^1\) has degree three.
Their product is finite, so every component of the
inverse-image divisor in (1) dominates the
multiplication graph. Its generic open is the
irreducible cover (2), which proves geometric
irreducibility and generic reducedness. The divisor
is of class \(3L\) for the product polarization
\(L=\boxtimes_i\mathcal O_{E_i}(O)\), since the graph
has multidegree \((1,1,1)\). Thus \(L^2\cdot\mathcal T
=3L^3=18\).

For the smooth model, let \(B=E\times E\), let
\(H,V\) be the two origin fibres of its projections,
and let \(D=\sum_{j=0}^2\Gamma_j\) with
\(\Gamma_j=\{(P,-\zeta^j P):P\in E\}\), where \(\zeta\)
acts on \(E\) by \((x,y)\mapsto(\zeta x,y)\).
The divisor of \(y_1+y_2\) is \(D-3H-3V\).
The three graphs meet at the same three points
\((P,-P)\) with \((1-\zeta)P=0\), and their tangent
directions are distinct there. Each pair has
intersection number three and \(D^2=18\).
Blow up these three points on \(B\), writing
\(b:\widetilde B\to B\) and \(F_j\) for the exceptional
curves. The strict transforms of the graphs are disjoint
and give the branch divisor \(D'=b^*D-3\sum F_j\).
At a finite triple point the valuation of \(y_1+y_2\)
on the exceptional divisor is three; at \((O,O)\)
it is \(3-3-3=-3\). Hence there is no branch along the
exceptional divisors. The normalized cyclic cover
\(\pi:\widetilde{\mathcal T}\to\widetilde B\) is smooth,
since it is locally a cyclic cover branched simply onto a
single smooth branch component, or étale.

Its canonical divisor is the pullback of
\[
              (2/3)b^*D-\sum F_j .
 \tag{3}
\]
The square of (3) is \(5\). Its intersections with the
strict transforms of the \(\Gamma_j\) and with the \(F_j\)
are all one. For any other curve \(C\) on \(B\),
\(D\cdot C\) is at least three times the sum of the
multiplicities of \(C\) at the triple points, because
each of the three smooth branches contributes at least
the multiplicity there. Also \(D\cdot C>0\): if each
intersection were zero, the three maps with graph
fibres would all be constant on \(C\), and two of
those maps form an isogeny \(B\to E^2\). The
Nakai–Moishezon criterion gives ampleness of (3), and
\(K^2=3\cdot5=15\).

The cyclic root line bundle on the blown-up base is
\(M=b^*(H+V)-\sum F_j\). The trace splitting of
\(\pi_*\mathcal O\) is
\(\mathcal O\oplus\mathcal O(-M)\oplus\mathcal O(-2M)\).
Surface Riemann–Roch gives Euler characteristics
\(0,1,1\), since \(M^2=-1\) and \(M\cdot K_{\widetilde B}=3\).
For the canonical sections, \(K_{\widetilde B}+M=b^*(H+V)\)
has one section. The sections of
\(K_{\widetilde B}+2M=b^*(2H+2V)-\sum F_j\)
come from the span of \(1,x_1,x_2,x_1x_2\) on \(B\)
with the indicated poles. Vanishing at the two finite
triple points imposes the same condition, killing the
constant term, while vanishing at \((O,O)\) kills the
\(x_1x_2\) coefficient. There are two remaining
sections. Together with the unique base canonical
section, this gives \(p_g=4\), hence \(q=3\).

It remains to check the low-genus curves for the
arithmetic application. A curve with rational or elliptic
normalization maps to a point or an elliptic translate
in \(B\). The only curves mapping to a point lie over
the exceptional \(F_j\), each of which lifts to one
elliptic curve: its degree-three cover is branched at
the three distinct strict branch intersections.
If the image in \(B\) is one of the \(\Gamma_j\),
its preimage is the ramification elliptic curve.
For any other elliptic translate \(C\), every intersection
with a branch graph away from the three triple points
is transverse. Such a point makes the restricted cover
ramified and connected, so its genus is greater than one.
Thus a necessary condition for genus at most one is
that all intersections occur at the triple points.

Every elliptic subgroup of \(E^2\) has a parametrization
\(P\mapsto ([u]P,[v]P)\) for a coprime pair
\(u,v\in\mathbf Z[\zeta]\). Indeed it is isogenous to \(E\);
start with an isogeny from \(E\) and divide the two
endomorphisms by their common factor in the Eisenstein
integer ring, which is a principal ideal domain.
The resulting parametrization is injective by Bézout.
The intersection numbers with the three graphs are
\(|v+\zeta^j u|^2\), with sum \(3(|u|^2+|v|^2)\).
If \(C\) is parallel to a graph but distinct, it
contains no triple point and meets the other graphs,
so is excluded. Otherwise the intersection numbers
are each at most three if the condition holds.
When both \(u,v\) are nonzero, this forces both to
be Eisenstein units, and their ratio is either parallel
to a graph or gives one intersection number four.
Thus only the horizontal and vertical translates
through a triple point remain. For a nonorigin fixed
point \(P\), the restriction on the vertical fibre
has divisor \(3[-P]-3[O]\), which is not a cube divisor
of a function because \([-P]-[O]\) is not principal.
This produces one connected unramified cubic cover
of that fibre, an elliptic curve. For the origin fibre,
after removing the cube of a local parameter the generic
residue of \(y_1+y_2\) is a nonzero constant, so it
splits into three copies of the fibre. The same counts
hold for the horizontal fibres. This gives
\(3+3+2(2+3)=16\) elliptic curves and no rational curves.
The rational map from the smooth model to the ambient
abelian variety extends to a morphism. Resolve its
indeterminacies by point blowups and use the fact that
all the exceptional rational curves map to points.
The resulting morphism sends all sixteen elliptic curves
to the boundary of the open set in (1).

If an elliptic translate in \(\mathcal T\subseteq E^3\)
met \(\mathcal U\), its proper transform in the
birational model would have genus one, contradicting
the classification. A two-dimensional abelian translate
in \(\mathcal T\) would equal the entire surface and
would contain an elliptic translate through every
point, since every two-dimensional abelian subvariety
of \(E^3\) has a nontrivial elliptic subgroup. This
completes the exclusion for \(\mathcal U\).

## Coherent paths with any three fixed cube classes

**Corollary ML2.** Fix a negative quadratic radicand \(-d\)
with \(d>0\) squarefree and any three positive cube-free
integers \(w_{12},w_{23},w_{13}\). Put
\(E_{ij}:y^2=x^3-dw_{ij}^2/4\) and let
\(\rho=\sum_{ij}\operatorname{rank}E_{ij}(\mathbf Q)\).
The number of ordered coherent three-prime norm witnesses
with these cube-free edge parts is at most
\(C^{1+\rho}\) for an absolute \(C\). There is no bound
on the prime sizes in this statement.

Write \(v_{ij}=w_{ij}b_{ij}^3\) with signed integer
\(b_{ij}\ne0\), following
[the coherent norm equations](ABC_QUADRATIC_THREE_TORSION_FERMAT_TAIL.md#witness-selection).
Then
\[
 P_{ij}=\left(\frac{p_i p_j}{b_{ij}^2},
                          \frac{u_{ij}}{2b_{ij}^3}\right)
         \in E_{ij}(\mathbf Q).
\]
With \(a_{ij}=w_{ij}\sqrt{-d}/2\), the three \(t\)-coordinates
in (1) are the ratios of the edge generators to their
conjugates. Thus their relation is exactly the quotient of
the coherent product equation by its conjugate. All the
points are finite and have nonzero \(x\). The relation
is defined over \(\mathbf Q\), since conjugating
\(\sqrt{-d}\) inverts all three \(t_i\).
Apply [Gao–Ge–Kühne, Theorem 1.1](https://pmihes.centre-mersenne.org/item/10.5802/pmihes.26.pdf)
to the degree \(18\) irreducible surface in the
three-dimensional product abelian variety, with
\(\Gamma=\prod E_{ij}(\mathbf Q)\) of rank \(\rho\).
Every coset in the theorem meeting \(\mathcal U\) is
a point by Theorem ML1. For each such rational
point, the square classes of the three \(x\)-coordinates
give the three unordered prime pairs, hence determine
the ordered prime triple. The three absolute values
of \(b_{ij}\) then follow, so the remaining sign
choices have bounded multiplicity. Absorb this in \(C\).

An exact genuine triple with different cube-free parts
comes from the
[nonmaximal order example](ABC_COHERENT_NORM_MARKED_REPRESENTATIONS.md#a-nonmaximal-order-example).
Here \(d=10015\), the primes are \((7,17,23)\),
and the positive cube-free parts are \((22,2,5)\),
with \(b=(1,3,2)\). The three rational points are
\((119,688),(391/9,7244/27),(161/4,409/8)\).
The edge generators multiply as
\((688+11\sqrt{-10015})(7244+27\sqrt{-10015})
=17^3(409+20\sqrt{-10015})\).

## Uniformity across the curve parameters

**Proposition ML3.** For \\(d\\asymp X\\), \\(d\\ne3\\), and
\\(w_{ij}\\le X^B\\) with fixed \\(B\\), the multiplicity of a
fixed-colour fibre is bounded by (5) for every
\\(\\varepsilon>0\\), uniformly in those parameters.

After rescaling by \(x\mapsto4x,\ y\mapsto8y\) the
curves \(E_{ij}\) are \(y^2=x^3-d(4w_{ij})^2\).
The local conditions in
[the odd-prime Jacobian comparison](INERT_PRIME_TWISTS_AND_CYCLOTOMIC_SELMER_CHARACTERS.md#local-selmer-conditions)
give, at \(\ell=3\), a rank upper bound
\[
 \operatorname{rank} E_{ij}(\mathbf Q)
 \ \le\ r_3(\mathbf Q(\zeta_3,\sqrt{-d}))
                                    +O(1+\omega(w_{ij})) .
 \tag{4}
\]
For \(d\ne3\), this follows by mapping Selmer classes to the
finite valuation coordinates at \(2,3\) and the split
primes dividing \(w_{ij}\), with kernel controlled by the
class group and the fixed unit rank. Primes dividing \(d\)
away from \(2,3\) are nonsplit and contribute no coordinate.
The rational rank is bounded by the cyclotomic Selmer
dimension as in QN1, via semilinear descent and the free
cyclotomic Mordell–Weil module. The biquadratic character
decomposition at \(3\) gives
\(r_3(\mathbf Q(\zeta_3,\sqrt{-d}))=r_3(\mathbf Q(\sqrt{-d}))
+r_3(\mathbf Q(\sqrt{3d}))\) here: the cyclotomic quadratic
field has class number one. Thus for polynomially bounded
\(w_{ij}\) the fixed-colour multiplicity is at most
\[
       X^\varepsilon
       \big(h_3(-d)h_3(\operatorname{sf}(3d))\big)^\lambda
 \tag{5}
\]
for an absolute \(\lambda\) determined by the
Mordell–Lang constant, where \(d\asymp X\).
In particular the multiplicity is subpower uniformly
on fields where both \(3\)-torsion sizes are bounded by a
fixed power of \(\log X\). This bound covers all heights
in those fixed-colour fibres.

The original rank-one problem varies the field and all
three cube-free parts, and its high-torsion fields need
not have small Jacobian rank. To use this bridge there,
one needs a bound on the number of occupied colour
tuples weighted by \(C^{1+\rho}\); a rank bound
or a point bound in one fixed fibre does not give
that sum. Equations (1) and (2) apply without the
lattice rank-one restriction, and do not claim a new
exponent for the full prime-path tail.
