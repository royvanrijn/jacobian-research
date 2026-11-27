# Outer cusp restrictions for the order-\((12,42)\) Weyl transfer

The local sheet excess gives restrictions on the singular branches of the remaining outer curve. Rank-one germs have the largest local excess allowed by the common derivative. A branch of multiplicity two has an explicit analytic model, and its contribution depends on a divisor of the transverse ramification index. These two statements give separate restrictions on how the outer chart can realize its prescribed total excess.

Consider a finite holomorphic map germ from a smooth surface to a smooth surface. Assume its reduced ramification locus is a single smooth curve \(E\) with ramification index \(e>1\), and that its restriction to \(E\) normalizes its image branch \(B\). The branch in this statement is the image of the selected source germ, even when other source points over the same target value contribute other image branches.

## Rank one

**Proposition SmoothRamificationRankOneCuspDegree.**
Suppose the local map has differential of rank one at a point where the branch \(B\) has multiplicity \(a\ge2\). Then \(B\) has one Puiseux pair \((a,a(e-1)+1)\). The local degree is \(a(e-1)+1\), so the degree excess over \(e\) is \((a-1)(e-1)\).

Choose a target linear coordinate \(p\) whose differential pulls back nontrivially and whose order on the branch normalization is \(a\). Both conditions are open in the space of linear coordinates. Choose a parameter \(v\) on \(E\) so that \(p|_E=v^a\). Since \(a\ge2\), \(d p\) vanishes on the tangent to \(E\) at the origin. Thus any extension of \(v\) gives local source coordinates \((p,v)\). In these coordinates \(E=\{p=v^a\}\) and the other target coordinate satisfies
\[
        q(p,v)=q(p,0)+\int_0^v (p-z^a)^{e-1} L(p,z)\,dz,
                        \qquad L(0,0)\ne0.
\]
Indeed the Jacobian in these coordinates is \(q_v\), a unit times the indicated power. Subtract \(q(p,0)\) by a target coordinate change. At \(p=0\), \(q\) has order \(a(e-1)+1\) in \(v\), which is the local degree. At \(p=v^a\), the same order occurs, with nonzero coefficient \(L(0,0)\int_0^1(1-t^a)^{e-1}\,dt\). The two orders \(a,a(e-1)+1\) are coprime, giving the stated Puiseux pair.

## Monomial image branches and multiplicity two

**Theorem SmoothRamificationMonomialBranchDegree.**
If the image branch is analytically equivalent to \(q^a=p^b\), with coprime \(2\le a<b\), then
\[
                   \deg_{\rm loc}=\frac{ab e}{a+b-1}.
\]

Choose analytic target coordinates realizing this equation. Rotations \((p,q)\mapsto(\zeta^a p,\zeta^b q)\) preserve the branch and a small target bidisc. Lift their isotopy, starting at the identity, to the unramified complement of a finite representative of the germ. The lift is holomorphic in the source for each rotation. One turn acts as a deck transformation of finite order, so after a finite cover of the circle there is a circle action on the source complement. Each lift extends across the omitted analytic set. To see this, work in source coordinates on a bounded neighbourhood of the origin and shrink the target bidisc so its entire preimage is contained there. The lifted coordinate functions are bounded and extend by removal of analytic singularities. The extensions commute with the finite map and fix its sole point over the target origin. The inverses and action law extend by equality on the complement. Continuity of the extensions on smaller neighbourhoods follows by normal families, since their limits on the complement are unique.

Let \(L>0\) be the exponent of the cover of the target circle. On \(E\) the action has weight \(L\), since the normalization of \(q^a=p^b\) has parameter \(v\) with \(p=v^a,q=v^b\). Let \(K\) be the weight on the normal tangent line of \(E\). Averaging local coordinate functions against their tangent characters gives analytic coordinates \(r,v\) of weights \(K,L\), with \(E=\{r=0\}\). In particular one can average a defining function of \(E\) for \(r\). Its differential remains nonzero in the normal line. For \(v\) choose a function with nonzero tangent differential and average for weight \(L\).

The functions \(p,q\) in these coordinates have weights \(aL,bL\). Their determinant has weight \((a+b-1)L-K\) and is \(r^{e-1}\) times a unit. Evaluating the transformed unit at the origin gives
\[
                       eK=(a+b-1)L.
\]
Consequently \(K>0\), and \(p,q\) are weighted homogeneous polynomials. Finiteness makes them a homogeneous regular sequence. Their colength is the limit at \(t=1\) of \((1-t^{aL})(1-t^{bL})/((1-t^K)(1-t^L))\), namely \(abL/K\).

[Orevkov's classification, Theorem 2,](https://www.math.univ-toulouse.fr/~orevkov/k-en.pdf) has the additional hypothesis that every component above the image branch maps bijectively to it. The circle argument here uses birationality only on the ramification curve. [Kulikov's ADE rigidity theorem](https://arxiv.org/abs/1911.10848) concerns deformation equivalence. The degree calculation uses the circle action of the specified analytic equation directly.

**Corollary SmoothRamificationMultiplicityTwoNormalForm.**
For a branch of multiplicity two, write its Puiseux pair as \((2,2d-1)\), with \(d\ge2\). Necessarily \(d\mid e\). Put \(u=e/d\). In analytic coordinates its local map is
\[
       p=v^2+r^u,\qquad
       q=\frac{\displaystyle\int_0^v(p-y^2)^{d-1}\,dy}
               {\displaystyle\int_0^1(1-t^2)^{d-1}\,dt},
\]
where \(p\) is held fixed in the numerator. It has local degree \((2d-1)u\) and excess \(e-u\). Its differential has rank one exactly when \(u=1\).

For the analytic equation, choose a normalization parameter with \(p=v^2\). Subtract the even part of the other coordinate and divide its odd part by a unit in \(p\). This gives \(q=v^{2d-1}\) on the normalization, hence the monomial branch equation. The degree formula gives \(e(2d-1)/d\), so \(d\mid e\) and \(K/L=2/u\).

Normalize the averaged \(v\) so \(p|_E=v^2\) and \(q|_E=v^{2d-1}\). The weight \(2L\) allows only \(v^2,r^u\), and \(v r^{u/2}\) when \(u\) is even. Complete the square, preserving the coordinates on \(E\). The coefficient of \(r^u\) is nonzero since \(p=v^2\) would make the Jacobian divisible by \(v\). Scale \(r\) to make that coefficient one. The Jacobian in the weighted coordinates is a constant times \(r^{e-1}\), since the unit after division is invariant under positive weights and therefore constant.

At fixed \(p\), the derivative of the displayed integral with respect to \(v\) is \((p-v^2)^{d-1}\), up to its normalization constant. Thus the integral has Jacobian a nonzero constant times \(r^{e-1}\). Any two weighted \(q\)-polynomials with the same determinant differ by a weighted polynomial \(D\) that Poisson-commutes with \(p\). If \(D\ne0\), Euler's identity and the vanishing determinant imply \(D^2/p^{2d-1}\) is constant. The polynomial \(p=v^2+r^u\) is squarefree, making this impossible. Matching the determinant therefore fixes \(q\) up to a scalar multiple of the displayed integral, and its restriction \(q|_E=v^{2d-1}\) fixes the scalar. The formula also shows that the degree-one term of \(p\) is present precisely for \(u=1\), while \(q\) has no degree-one term.

## Consequences for the order-\((12,42)\) outer chart

**Corollary PlanarFortyTwoOuterRankZeroObstruction.**
Write \(H=2m\) for the number of high roots in the transferred profile. Except at \(m=1,4\), the outer chart must contain a point where the differential has rank zero. At \(m=1\) an everywhere rank-one chart would have singular branch Puiseux pairs in one of the following three lists:
\[
             (2,7)^3,\qquad (2,7),(3,10),\qquad (4,13).
\]

Here \(e=5m-1\) and the [total local degree excess](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_TRIANGULAR_FIBRES.md#local-sheet-counts-on-the-outer-curve) is \(9m\). Smooth image branches contribute zero, and the rank-one formula makes every singular contribution divisible by \(e-1=5m-2\). Thus \(5m-2\mid9m\), so \(5m-2\mid18\). The only positive \(m\) are \(1,4\). At \(m=1\), the sum of \(a-1\) over the singular branches must be three. Its three partitions give the lists above.

**Corollary PlanarFortyTwoMultiplicityTwoOuterProfiles.**
Suppose every singular analytic branch of the outer branch curve has multiplicity two. Then \(m\le335\), and the following are all possible \(m\) under the local degree and affine delta budgets:
\[
 1,2,3,4,5,8,10,11,17,20,23,29,31,38,47,65,80,101,119,173,227,335 .
\]
At \(m=1\), the list of singular branch pairs must be either three \((2,7)\) branches or three \((2,3)\) branches and one \((2,7)\) branch. At \(m=2\), it must consist of three \((2,5)\) branches. In particular, the degree-sixteen profile has a branch of multiplicity at least three or a branch with pair \((2,7)\). The budget allows the displayed combinations as numerical possibilities.

Indeed, for branch pairs \((2,2d_i-1)\) one has \(d_i\ge2\), \(d_i\mid(5m-1)\), and
\[
        \sum_i (5m-1)(1-1/d_i)=9m,\qquad
        \sum_i (d_i-1)\le99 .
\]
The second relation uses the [total affine delta invariant](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_TARGET_CYLINDER.md#the-plane-curve-of-outer-branch-values). The delta contribution of a multiplicity-two branch is \(d_i-1\); intersections of distinct branches only increase the total. There are at most four singular branches since each contributes at least \((5m-1)/2\). There must be at least two since the total \(9m\) exceeds \(5m-1\). The sum bound gives \(d_i\le100\). For each nondecreasing tuple of length two to four in this range, write \(s=\sum_i(1-1/d_i)\), solve \(m=s/(5s-9)\), and keep positive integral solutions meeting the divisibility tests. The [exact checker](verify_weyl_forty_two_outer_cusps.py) enumerates the 29 tuples. The two smallest cases also follow directly from \(2,4\mid4\) and \(3,9\mid9\) as the respective possible divisors greater than one.

For the smallest case the
[index-four jet calculation](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_INDEX_FOUR_CUSPS.md)
also treats rank-zero points of branch multiplicity greater than two.
It gives six necessary Puiseux-pair lists and the projection indices
of the ordered double locus at their cusps.
