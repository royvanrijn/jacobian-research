# Torsion in Soulé's genus-six fibration

Soulé's [arXiv:2611.05217v1, Proposition 9 and §§6.2–6.3](https://arxiv.org/html/2611.05217v1#S6)
uses a genus-six family \(\phi:X\to Y\) that double covers a fixed
genus-three curve \(D\) on each fibre. If \(k=\deg(Y\to D)\), the
argument uses a lower bound \(-3k\) for the degree of
\(\phi_*(\omega_{X/Y}-2B_1-2B_2)\). The degree on each connected
geometric component is instead \(-13k/2\). In particular the
positive coefficient in (34)–(36) does not follow from this
fibration and auxiliary bundle.

## Square roots over the exceptional fibres

Work over an algebraic closure of the number field. Write
\(\chi:D\to C\) for the étale double cover of the genus-two base
curve and \(\nu\) for the free involution. The curve \(D\) is
hyperelliptic [by Proposition 4(i)](https://arxiv.org/html/2611.05217v1#S3).
Write \(\iota\) for its hyperelliptic involution and \(H\) for the
hyperelliptic pencil class on \(D\), so \(K_D=2H\). The involution
\(\nu\) descends to a nontrivial involution of the hyperelliptic
line. Choose its fixed points to be \(0,\infty\); neither is
a branch value since \(\nu\) is free. We can write
\[
 D:\ y^2=\prod_{i=1}^{4}(z^2-a_i),\qquad
 \nu(z,y)=(-z,-y),\qquad \sigma=\iota\nu:(z,y)\mapsto(-z,y),
\]
with distinct nonzero \(a_i\). The sign in \(\nu\) follows from
freeness over \(z=0\). Put
\[
 E:\ y^2=\prod_{i=1}^{4}(x-a_i),\qquad q:D\to E,\quad (z,y)\mapsto(z^2,y).
\]
The map \(q\) is a ramified double cover of the elliptic curve \(E\).
Its four ramification points are the two above \(z=0\) and the two
above \(z=\infty\). They are precisely the exceptional points \(p\)
with \(p+\nu p\sim H\), since this condition means \(\nu p=\iota p\).
All four lie outside the Weierstrass locus of \(D\).
The quotient \(C\) can be written \(v^2=x\prod_i(x-a_i)\), with
\(x=z^2\) and \(v=zy\). The two exceptional pairs have images
\(W_0,W_\infty\) on \(C\), and the eight Weierstrass points of \(D\)
map in pairs to the four other Weierstrass points on \(C\).

Over \(p\in D\), let \(M\in\operatorname{Pic}^1(D)\) satisfy
\(2M=p+\nu p\), using additive notation for line bundles. Fix
one pair \((p_0,M_0)\), put \(J=\operatorname{Pic}^0(D)\), and write
\(\xi=M-M_0-(p-p_0)\). Since \(\iota\) acts by \(-1\) on \(J\)
and \(q^*q_*=1+\sigma^*\) on \(J\), the equation on \(\xi\) is
\[
  2\xi=(\nu p-p)-(\nu p_0-p_0)
       =-q^*\bigl([q(p)-q(p_0)]\bigr).
\]
The map \(q^*:\operatorname{Pic}^0(E)\to J\) is injective.
A nontrivial kernel element would be of order two, by the norm map.
Its étale double cover of \(E\) would split after pullback to
\(D\), factoring the ramified degree-two map \(q\) through an
étale degree-two map. With \(P=q^*\operatorname{Pic}^0(E)\), the
inverse image of \(P\) under \([2]:J\to J\) consists of the
sixteen translates of \(P\) indexed by \(J[2]/P[2]\).
On each translate the above equation is the base change of
elliptic multiplication by two along \(q\), up to translation and
sign. This fibre product is connected: the ramified quadratic
function field extension for \(q\) is disjoint from the degree-four
Galois étale extension for multiplication by two. Thus the
square-root cover has sixteen connected components of degree four,
and on each fibre within one component the four roots differ by
\(P[2]\). This is the cover used in §5.1, via the change
\(M\mapsto M-M_0-(p-p_0)\). A choice of global root in §5.1
can shift the resulting root map by a constant element of \(J[2]\).
Indeed two such maps into \(\operatorname{Pic}^1(D)\) differ by a
map to the finite group \(J[2]\), hence by a constant on a connected
component. Further finite étale base change multiplies the point
counts and degrees below by its degree.

For exceptional \(p\), the effective square roots of \(H\) are
exactly the eight classes of Weierstrass points of \(D\).
Two distinct such classes differ by no element of \(P[2]\).
To see this, write \(Q_i=(a_i,0)\in E\) and denote its two
preimages by \(W_i^+,W_i^-\). Every nonzero element of \(E[2]\)
is \(Q_i-Q_j\) for some \(i\ne j\). Its pullback is represented
by a four-element subset of the eight Weierstrass points, while
the difference of two Weierstrass classes is represented by a
two-element subset. Here an even subset \(S\) represents
\(\sum_{W\in S}W-(|S|/2)H\); addition is symmetric difference
since \(2W\sim H\). The only trivial subsets are the empty and
full sets. In fact for a proper nonempty even subset one can
take its complement if needed, using \(\sum W\sim 4H\),
and assume \(|S|=2\) or \(4\). Every effective divisor in
\(|H|\) or \(|2H|\) is a sum of whole hyperelliptic fibres by
the splitting of the degree-two pushforward
\(\mathcal O\oplus\mathcal O(-4)\) on the line. A sum of
distinct branch points cannot be such a divisor. So at most
one effective root occurs above any specified exceptional \(p\)
on a component.

On a connected component the degree-zero class
\(\chi_*(M-p)\) is constant, since twice it is zero.
If \(M=W\) is effective over an exceptional point \(p\), this
class is \(\chi(W)-\chi(p)\). The four possible differences for
\(\chi(p)=W_0\) and the four for \(\chi(p)=W_\infty\) are disjoint.
Indeed on a genus-two hyperelliptic curve differences of distinct
Weierstrass points determine their unordered pair. This follows
by the same even-subset argument with six branch points, taking
complements so that the only nonempty proper even subsets to
check have size two. Hence each connected component has
effective roots over at most one exceptional pair, with at most
two effective points in total. The full cover has four exceptional
base points with eight effective roots each, giving thirty-two
effective points on sixteen components. The upper bound of two
is therefore an equality on every component. A degree-\(k\)
component, including after finite étale base change, has \(4k\)
exceptional points and exactly \(k/2\) points where \(M\) is
effective.

## Torsion lengths and the degree calculation

Write \(\rho:X\to D\times Y\) for the double cover,
\(f:D\times Y\to Y\) for the projection, and \(A_1,A_2\) for
the branch sections. Write \(\mathcal M\) for the global root
with \(2\mathcal M=A_1+A_2\) and \(B=2B_1+2B_2\) for the
pullback of \(A_1+A_2\). The two summands of \(R^1\phi_*(B)\)
are \(R^1 f_*(A_1+A_2)\) and \(R^1 f_*(\mathcal M)\).
The corresponding jumps in \(h^0\) both have size one.
For the first they occur at the \(4k\) exceptional points,
because the hyperelliptic class is the only degree-two
class with two sections. For the second they occur at the
\(k/2\) effective roots already counted. Outside the
exceptional points no effective \(M\) is possible:
\(M=W\) would imply equivalent effective divisors
\(2W\) and \(p+\nu p\) with distinct supports, forcing a pencil.

Each jump contributes length one. Here is the local check.
On a smooth base curve, represent cohomology by a two-term
complex of finite free modules over the local discrete
valuation ring. A jump of size one in \(h^0\) gives one
nonunit elementary divisor of the differential. Its exponent
is one if the first-order obstruction from the extra section
to lifting is nonzero on the special fibre. This follows by
reducing the diagonal form of the complex modulo the square
of the uniformizer. For a line bundle on the fixed curve \(D\),
the first-order obstruction is cup product with the tangent
to the corresponding map to the Picard variety, as follows by
differentiating the transition functions. Under Serre duality
it pairs with multiplication of sections on \(D\).

At exceptional \(p\) the tangent of \(p+\nu p\) in the Picard
variety is \(2\,d\operatorname{Abel}_p\). Indeed
\(p+\nu p\sim H+p-\sigma p\), and \(d\sigma=-1\) at its fixed
points. Under Serre duality this tangent acts on canonical
differentials by nonzero scalar times evaluation at \(p\).
For \(A_1+A_2\) at an exceptional point, the Petri pairing is
therefore the evaluation pairing on
\(H^0(H)\otimes H^0(H)\to H^0(2H)\).
It is a nonzero rank-one pairing because \(H\) is basepoint
free. Thus the unique jumping direction has nonzero
obstruction. For \(\mathcal M\) the tangent is
\(d\operatorname{Abel}_p\). At an effective \(M=W\), the products
of sections of \(W\) and \(K_D-W\) span \(H^0(K_D(-W))\) inside
\(H^0(K_D)\). These sections have common support only at \(W\),
as \(K_D=2H\) and \(W\) is a branch point of the hyperelliptic
map. Since \(p\ne W\), the obstruction is again nonzero. The
argument is unchanged by étale pullback.

Consequently the torsion has degree \(4k+k/2=9k/2\), in place
of the lower bound \(8k\) in (24). The ordinary pushforward
\(\phi_*B\) is trivial of rank one: the summand for
\(\mathcal M\) has generic \(h^0=0\), and the summand for
\(A_1+A_2\) has generic \(h^0=1\) with its section from the
effective divisor everywhere nonzero as a fibre section.
The intersection calculation of Proposition 8 and (18) uses
\(B^2=-16k\), \(\omega_{X/Y}\cdot B=8k\), and
\(\omega_{X/Y}^2=12k\), hence \(\deg R\phi_*B=-11k\).
Relative duality identifies
\(\phi_*(\omega_{X/Y}-B)\) with the dual of the torsion-free
quotient of \(R^1\phi_*B\), giving
\[
 \deg\phi_*(\omega_{X/Y}-B)=-11k+\frac{9k}{2}=-\frac{13k}{2}.
\]

For clarity, the generic generation used in the degree comparison
can be verified on a fixed root component. At a Weierstrass point
\(p\) of \(D\), also \(\nu p\) is Weierstrass. Proposition 4(iii)
gives \(h^0(M+p)=h^0(M+\nu p)=0\) for every root \(M\).
The residual divisor of \(K_D-p-\nu p\) is \(p+\nu p\)
on this fibre since both points are Weierstrass. In the basepoint
criterion of §3.5.3 the two residual points \(R_i\) are therefore
the original branch points, and the two vanishings suffice.
Generic generation follows by openness on the same component.
There is no need to replace the double cover of a fibre by one
defined from a different square root in §6.3.3.

## Consequence for the pointwise bound

With \(L=\omega_{X/Y}-B\), the intersection calculation gives
\(L^2=-20k\) and the fibre degree of \(L\) is \(6\).
Even granting the remaining bounds in §6.3, the coefficient
of the base height on the side from (31) is now
\(-6\deg\phi_*L=39k\), while its coefficient from \(-L^2\)
is \(20k\). Subtracting gives a negative coefficient, so this
comparison gives no upper bound for the base height in terms
of the discriminant. In particular v1 supplies no deduction
of its claimed finiteness for \(\lambda<1/24\) by this route.
The obstruction concerns the specified proof of the pointwise
bound, whose truth remains open here.

The freedom to change the root already exposes the erroneous
inference from Proposition 5(ii) to Proposition 9(iv) in the
[public discussion of v1](https://www.reddit.com/r/math/comments/1xnhcxe/comment/19s4d3xb1grr9/).
That discussion also shows the incompatibility of the proposed
numbers with a basepoint-free rank-two pencil by a top Chern class
argument. The component count and first-order calculation show why choosing
a different connected component cannot restore (24).
[Zeng's coefficient audit, §6 and Appendix I](https://hexagonmath.org/pdf/2611.00330v2)
treats a local model at an effective root and a separate norm
comparison, while the degree requires counting the effective
roots across the whole component.


## A geometric bound for the rank-two pencil comparison

Let \(f:S\to Y\) be a smooth projective surface fibration in
geometrically integral curves over a smooth projective curve.
Suppose \(E=f_*L\) is rank two and its evaluation generates \(L\)
on the generic fibre, where \(L\) has degree \(d>0\). Then
\[
                         L^2\ \ge\ d\deg E .
\]
In fact, \(E\) is locally free since it is torsion free.
At a closed fibre the natural map
\(E_y\to H^0(S_y,L_y)\) is injective, by the multiplication
sequence for a base uniformizer. Hence the common zero locus of
\(f^*E\to L\) contains no entire fibre. The generic generation
excludes horizontal components, so the base ideal has finite support.
On a resolution of its blowup \(\pi:\widetilde S\to S\), the
moving line bundle is \(\widetilde L=\pi^*L-D\) with \(D\) effective
and exceptional. It gives a map to \(\mathbb P_Y(E)\) in the
quotient convention of generic degree \(d\). There
\(\widetilde L\) is the pullback of \(\mathcal O(1)\), whose
square is \(\deg E\) by the projective bundle formula.
Thus \(d\deg E = \widetilde L^2=L^2+D^2\le L^2\), using
\(\pi^*L\cdot D=0\) and negativity of exceptional intersections.

The coefficient \(-L^2+d\deg E\) is therefore always
nonpositive for this use of the rank-two arithmetic degree
in the height comparison. There is also a cross-intersection
version. For two such pencils \((L_i,E_i)\), let \(d_i\)
be their fibre degrees and \(F\) the numerical class of
a fibre. The classes \(L_i-(\deg E_i)F/2\) have nonnegative
square and positive intersection with \(F\). They lie in the
same closed positive cone by the Hodge index theorem, hence
\[
       -L_1\cdot L_2+\frac{d_2\deg E_1+d_1\deg E_2}{2}\le 0.
\]
This uses only the generic basepoint-free pencil hypothesis and
the indicated degree comparison. It applies independently of
the explicit torsion calculation.

## The asymmetric pencil after removing its fixed divisor

For the neighboring choice \(L_{31}=\omega_{X/Y}-3B_1-B_2\),
Riemann–Hurwitz gives \(L_{31}=\rho^*(K_D-A_1)\), where \(K_D\)
denotes its pullback to \(D\times Y\). On a general fibre the
anti-invariant section space is \(H^0(K_D-p-M)\), which is dual
in dimension to \(H^0(p+M)\) by Riemann–Roch in degree two.
It vanishes generically by specialization to a Weierstrass
\(p\), where Proposition 4(iii) gives vanishing for every root
\(M\). The invariant space \(H^0(K_D-p)\) has dimension two
and moving part \(|H|\) with fixed divisor \(\iota p\).
Write \(A_\iota\) for the graph of \(\iota\) of the first branch
section over \(Y\), and \(\psi:Y\to D\) for that section's
coordinate on \(D\). The preimage of the diagonal of the
hyperelliptic line under the two hyperelliptic maps gives
\[
 A_1+A_\iota\sim \operatorname{pr}_D^*H+f^*\psi^*H .
\]
The moving bundle is therefore \(L_{\rm mov}=L_{31}-\rho^*A_\iota
=\rho^*\operatorname{pr}_D^*H+\phi^*(-\psi^*H)\).
Its two generic sections are the pullbacks from \(H^0(H)\),
with the same base twist. The extra summand of
\(\phi_*L_{\rm mov}\) is zero because its generic rank is zero
and the pushforward is torsion free. Thus
\(\deg\phi_*L_{\rm mov}=2\deg(-\psi^*H)=-4k\).
Its fibre degree is \(4\) and its self-intersection is
\(2\cdot 4\cdot(-2k)=-16k\). The coefficient after subtracting
the two sides of the same Minkowski comparison is exactly zero.
This choice gives no discriminant upper bound either.
