# The index-four outer cusp profiles

The degree-sixteen order-\((12,42)\) transfer has transverse index four on its outer ramification curve. In this index the cubic Jacobian controls the first two map jets. A weighted comparison for a square quadratic jet then determines every singular branch type that fits the global sheet-excess budget.

Use the [outer local degree identity](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_TRIANGULAR_FIBRES.md#local-sheet-counts-on-the-outer-curve). In a finite smooth surface germ with one smooth ramification curve \(E\) of index four, write \(n\) for the local degree and \(\epsilon=n-4\) for its excess. Assume that the parametrization of its image branch by \(E\) is the normalization, as it is in the outer chart. If \(a\) is the multiplicity of a singular image branch, the common tangential derivative order is \(a-1\). In particular \(\epsilon\ge a-1\), and the local degree formula for a target coordinate \(p\) with order \(a\) on \(E\) is
\[
                        n=3a+1-\mu(p),
\]
where \(\mu(p)\) is its Milnor number at the source point. Indeed the polar identity \(I(p,J)=\mu(p)+n-1\) gives this formula since \(J\) is a unit times \(r^3\). On \(p_v=0\) the determinant gives \(\mu(p)+I(p_v,q_v)=3(a-1)\). Consequently \(\epsilon=I(p_v,q_v)\ge I(p_v,q_v,r)=a-1\). These local identities use just the smooth ramification germ.

## Rank-zero germs with small excess

**Proposition SmoothRamificationIndexFourRankZeroProfiles.**
For such a germ with rank zero at the source point and \(\epsilon\le9\), the image branch has one of the single Puiseux pairs in the table. Both rows with excess eight are excluded in an outer chart whose total excess is nine.

| Origin of the row | Orders on the branch | Local degree \(n\) | Excess \(\epsilon\) |
|---|---|---:|---:|
| A quadratic jet nonzero on \(E\) | \((2,3)\) | 6 | 2 |
| A normal square quadratic jet | \((3,4)\) | 8 | 4 |
| A normal square quadratic jet | \((4,5)\) | 10 | 6 |
| A normal square quadratic jet | \((5,6)\) | 12 | 8 |
| A mixed quadratic jet | \((4,9)\) | 12 | 8 |


### The quadratic and cubic jets

Choose source coordinates \(r,v\) with \(E=\{r=0\}\). The Jacobian is \(r^3\) times a unit. At rank zero both components start at degree at least two. At least one has order two since the Jacobian has order three. When both have a quadratic term, those terms have vanishing bracket and are proportional. Subtracting a scalar multiple makes the second component have order three, since order four or more would make the Jacobian have order at least four. Thus, writing the components as \(p,q\), their homogeneous terms satisfy
\[
                     J(p_2,q_3)=c_*r^3,\qquad c_*\ne0.
\]
The proportionality just used is the Euler identity for homogeneous forms of equal degree. We may choose \(p\) with the generic order on \(E\) by adding a multiple of \(q\), since this does not change \(p_2\).

If \(p_2|_E\ne0\), complete the square in \(v\) and scale \(p\) to get \(p_2=v^2+b r^2\). The case \(b=0\) would make \(J(p_2,q_3)\) divisible by \(v\). For \(b\ne0\), coefficient comparison gives
\[
                  q_3=A\left(v^3+\tfrac32 b r^2v\right),
                  \qquad A\ne0.
\]
Hence the image has pair \((2,3)\). The quadratic of \(p\) is nondegenerate, so \(\mu(p)=1\) and \(n=6\).

If \(p_2|_E=0\) and \(p_2\) has rank two as a quadratic form, a linear change preserving \(E\) gives \(p_2=rv\). If it has rank one, scaling gives \(p_2=r^2\). We treat the two cases separately.

### The mixed jet

Completing a nondegenerate quadratic direction analytically leaves a one-variable function with a nonzero quadratic term. Taking its unit square root and factoring the resulting sum of squares gives analytic coordinates \(s,t\) with \(p=st\). The two tangent lines to \(p=0\) are distinct, and \(E\) is tangent to one of them. Choose it to be \(s=0\), and write \(E\) as \(s=\theta(t)\), with \(\ell=\operatorname{ord}_t\theta\ge2\) finite since \(p|_E\) is nonconstant. The determinant identity becomes
\[
              (t\partial_t-s\partial_s)q=U(s,t)(s-\theta(t))^3,
                     \qquad U(0,0)\ne0 .
\]
No term \(s^i t^i\) occurs on the left. If \(\ell=2\), the coefficient of \(s^2t^2\) on the right is nonzero, coming only from the product of the unit's constant and \(-3s^2\theta(t)\). Thus \(\ell\ge3\).

Restricting to \(t=0\) gives \(\operatorname{ord}_s q(s,0)=3\); restricting to \(s=0\) gives \(\operatorname{ord}_t q(0,t)=3\ell\). Consequently \(n=3+3\ell\) and \(\epsilon=3\ell-1\). The budget \(\epsilon\le9\) leaves \(\ell=3\) only. To get the image pair, use weights \(3,1\) in \(s,t\), respectively. After subtracting polynomials in \(p=st\), the terms of \(q\) below weight nine vanish: the differential operator preserves weight and its kernel consists of powers of \(st\). If \(\theta(t)=d t^3+\cdots\), the weight-nine part of \(q\) is found by dividing the four coefficients of \(U(0,0)(s-dt^3)^3\) by the nonzero eigenvalues \(9-4i\), \(0\le i\le3\), for \(s^i t^{9-3i}\). Evaluating on \(s=d t^3\) gives a coefficient proportional to
\[
                -\frac19+\frac35-3-\frac13=-\frac{128}{45},
\]
so \(q|_E\) has order nine. The order of \(p|_E\) is four, giving \((4,9)\).

### The normal square jet

Write \(p=s^2+\gamma(v)\) by the one-variable analytic splitting lemma, leaving \(v\) as a parameter. Explicitly solve \(p_r=0\) for an analytic curve through the origin, subtract its critical value, and take the unit square root of the remaining quadratic factor. Since \(p_2=r^2\), the change satisfies \(s=r+O((r,v)^2)\) after initial scaling. The equation of \(E\) is \(s=\theta(v)\) with order \(\ell\ge2\), allowing \(\ell=\infty\) when \(\theta=0\). Set \(k=\operatorname{ord}_v\gamma\ge3\). This is finite: if \(\gamma=0\), then the critical curve \(s=0\) of \(p\) would have to be \(E\), making \(p|_E\) constant. In the new coordinates
\[
       2s q_v-\gamma'(v)q_s=U(s,v)(s-\theta(v))^3,\qquad U(0,0)=u_0\ne0 .
\]
The series \(q\) still has ordinary order at least three.

First suppose \(k>2\ell\), so that \(\ell\) is finite. For weights \((\ell,1)\) the leading \(p\) is \(s^2\) and the leading right side is \(u_0(s-dv^\ell)^3\), with \(d\ne0\). Any term of \(q\) of least weight below \(2\ell+1\) would commute with \(s^2\) and so be a power of \(s\) of ordinary degree at most two. Such a term cannot occur. At weight \(2\ell+1\) the leading bracket with \(s^2\) is divisible by \(s\), while the right side has a nonzero \(v^{3\ell}\) coefficient. This excludes \(k>2\ell\).

For \(k\le2\ell\) use weights \((k/2,1)\), and write \(p_*=s^2+c v^k\), \(c\ne0\), for the leading polynomial. It is squarefree. A homogeneous polynomial commuting with \(p_*\) has weight a multiple of \(k\) by the weighted Euler identity and squarefreeness. Below weight \(k+1\) the only possible nonconstant centralizer term is proportional to \(p_*\), which would give an ordinary quadratic term in \(q\). Thus the leading part of \(q\) is determined at weight \(k+1\) by the bracket equation.

If \(k<2\ell\), the leading right side is \(u_0s^3\), and the unique solution at this weight is
\[
                  q_*=\frac{u_0}{2}s^2v+\frac{u_0kc}{2(k+1)}v^{k+1}.
\]
Restricting to \(E\), where \(s\) has higher order than \(k/2\), gives the pair \((k,k+1)\).

If \(k=2\ell\), the right side is \(u_0(s-d v^\ell)^3\). Write
\[
                    q_*=\alpha s^2v+\beta s v^{\ell+1}+\gamma_* v^{2\ell+1}.
\]
Comparing the coefficients of \(s^3,s^2v^\ell,s v^{2\ell},v^{3\ell}\) yields
\[
   \alpha=\frac{u_0}{2},\qquad
   \beta=-\frac{3du_0}{2(\ell+1)},\qquad
   c=-\frac{(\ell+1)d^2}{3\ell},\qquad
   \gamma_*=\frac{(7-2\ell)d^2u_0}{6(2\ell+1)}.
\]
On \(E\) the leading \(p\)-coefficient is \(d^2+c=d^2(2\ell-1)/(3\ell)\), which is nonzero. The leading \(q\)-coefficient equals
\[
                \alpha d^2+\beta d+\gamma_*
                     =\frac{d^2u_0(2\ell-1)^2}{6(\ell+1)(2\ell+1)}.
\]
This gives the same pair \((k,k+1)\) in the equality case. In both cases \(a=k\) and the splitting presentation gives \(\mu(p)=k-1\), so \(n=2k+2\) and \(\epsilon=2k-2\). The budget gives \(k=3,4,5\), as in the table.

Finally every positive excess is at least two. For multiplicity at least three this is the bound \(\epsilon\ge a-1\). For multiplicity two it follows from the [multiplicity-two normal form](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_OUTER_CUSPS.md#monomial-image-branches-and-multiplicity-two) at index four, which allows excess two or three. A branch with excess eight would leave exactly one to be accounted for under the total budget nine. Hence neither excess-eight row can occur there.

## All singular-branch profiles at degree sixteen

**Theorem PlanarFortyTwoDegreeSixteenOuterProfiles.**
In a hypothetical degree-sixteen planar transfer of order \((12,42)\), the multiset of Puiseux pairs of the singular analytic branches of the affine outer branch curve must be one of the six lists below. The corresponding degree of the common derivative of its polynomial parametrization is shown on the right.

| Singular branch pairs | Local degrees | \(\deg\gcd(P_0',Q_0')\) |
|---|---|---:|
| \((4,13)\) | \(13\) | 3 |
| \((2,7),(3,10)\) | \(7,10\) | 3 |
| \((2,7),(4,5)\) | \(7,10\) | 4 |
| \((2,3),(2,7),(3,4)\) | \(6,7,8\) | 4 |
| \((2,7),(2,7),(2,7)\) | \(7,7,7\) | 3 |
| \((2,3),(2,3),(2,3),(2,7)\) | \(6,6,6,7\) | 4 |

Each profile includes a singular point of the parametrization at which the outer map has differential rank one. The sole \((4,13)\) branch in the first profile cannot share its target value with any other outer branch, even a smooth one. In each of the second and third profiles the two indicated singular branches have distinct target values.

For \(H=2\) the total excess is nine. At a rank-one singularity the [rank-one formula](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_OUTER_CUSPS.md#rank-one) gives pair \((a,3a+1)\) and excess \(3(a-1)\), so \(a=2,3,4\). At rank zero the proposition leaves pairs \((2,3),(3,4),(4,5)\) with excess \(2,4,6\), respectively. Smooth analytic branches contribute no excess. The partitions of nine using \(2,3,4,6,9\) are
\[
                      9,\quad 6+3,\quad 4+3+2,\quad
                      3+3+3,\quad 3+2+2+2 .
\]
The two possible pairs at excess six give the two lists with that part. The common derivative degree sums the \(a-1\) values across the singular branches. All odd excess comes from rank one. At any fixed target the local degrees of the distinct outer points over it sum to at most sixteen, since they contribute disjoint packets of sheets to the finite normalization. A smooth outer branch has local degree four. This proves the stated separations.

## The ordered double locus at these cusps

**Corollary PlanarFortyTwoDegreeSixteenCuspDoubleBranches.**
In the [ordered outer double algebra](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_POSITION_DOUBLE_LOCUS.md#specialization-in-the-outer-chart), the normalized branches arriving on the diagonal at a singular parameter have the following projection indices over the transverse parameter \(r=1/x\). The last column counts branches fixed by exchange of the two momenta at the origin of the normalization. A superscript denotes repetition.

| Puiseux pair | Projection indices on the ordered normalization | Fixed points |
|---|---|---:|
| \((2,7)\) | \(2^3\) | 3 |
| \((3,10)\) | \(3^6\) | 0 |
| \((4,13)\) | \(4^9\) | 3 |
| \((2,3)\) | \(1^2\) | 0 |
| \((3,4)\) | \(3^2\) | 0 |
| \((4,5)\) | \(2^6\) | 2 |

In the six profiles in the preceding theorem, in their displayed order, the total numbers of fixed points at the outer boundary are \(3,3,5,3,9,3\).

For a branch with coprime orders \((a,b)\), the local length of the special ordered double scheme at the diagonal is \((a-1)(b-1)\). Choose a parameter on \(E\) with \(p_0(v)=v^a\); by subtracting a function of \(p_0\) arrange \(q_0\) of order \(b\). The first divided difference factors into the \(a-1\) lines \(v_2=\zeta v_1\) with \(\zeta^a=1,\zeta\ne1\). The second has order \(b-1\) on each since \(a,b\) are coprime. Thus finding \((a-1)(b-1)\) distinct scaled ordered pairs with simple equations accounts for the entire local generic double scheme. Indeed the relative double algebra is finite free, and the length at the diagonal gives the rank of its local factor over the completed transverse ring.

At rank one, in coordinates \(p=v^a+s\) the new normal parameter \(s\) is the original \(r\) times a unit. The integral formula in the rank-one proof gives, for weights \((a,1)\) in \((r,v)\), leading \(p=v^a+\lambda r\) and leading \(q\) a nonzero scalar times
\[
                  \int_0^v(v^a+\lambda r-y^a)^3\,dy,
                          \qquad\lambda\ne0 .
\]
Set \(r=t^a\), \(v_i=tV_i\). For \(V_1\ne0\), the limiting first equation is \(V_2=\zeta V_1\), \(\zeta^a=1,\zeta\ne1\). The second reduces to \(I_a(V_1^a/\lambda)=0\), where
\[
                     I_a(z)=\int_0^1(1+z(1-y^a))^3\,dy .
\]
For \(a=2,3,4\) the discriminants of this cubic are \(-512/6125\), \(-19683/274400\), and \(-32768/570375\), respectively. This follows by expanding the integral; its coefficient of \(z^j\) is \(\binom 3j j!a^j/\prod_{i=1}^j(ia+1)\). The roots are simple and nonzero. At each of the \(3a(a-1)\) resulting pairs the first equation is simple in \(V_2\), and the derivative of the restricted second is nonzero. All lift analytically in \(t\). The substitutions \(t\mapsto\xi t\), \(\xi^a=1\), group them into orbits of \(a\); they fix the ratio \(\zeta\) and the chosen root of \(I_a\). Exchange takes \(\zeta\) to \(\zeta^{-1}\), so a branch is fixed exactly when \(\zeta=-1\). This proves the three rank-one rows.

For the rank-zero ordinary cusp \((2,3)\), the normal form has leading pair \(p=v^2+s^2,\ q=v^3+\tfrac32 s^2 v\), with \(s\) a unit times \(r\). At \(r=t\), \(v_i=tV_i\), the two limiting ordered pairs have \(V_2=-V_1\) and \(V_1^2=-3\lambda^2/2\), where \(\lambda\) is the nonzero leading coefficient of \(s/r\). They are simple and exchanged, with index one.

For rank-zero \((3,4)\), the strict-weight case in the square-jet proof applies with \(k=3\). Normalize \(s/r\) to leading coefficient one and write \(r=t^3,\ v_i=t^2 V_i\). The scaled polynomials are of the form \(1+cV^3\) and \(V(A+B V^3)\), with \(A,B,c\ne0\). There is one nonzero \(V^3\)-value giving equal second coordinates at the other cube-root ratios. The six simple ordered solutions give two branches of index three. Exchange interchanges the two ratios.

For \((4,5)\) in the strict-weight square case, \(r=t^2,\ v_i=tV_i\) similarly gives \(1+cV^4\) and \(V(A+B V^4)\). The twelve ordered solutions are simple. Quotienting the scaled solutions by \(t\mapsto-t\) gives six branches of index two, with the two branches of ratio \(-1\) fixed by exchange. In the equality case \(k=2\ell=4\), the leading relation is \(s=r+d v^2\) up to a nonzero rescaling of \(r\). The displayed square-jet coefficients, after rescaling \(v\) and the target, give the scaled polynomials
\[
                  1+2V^2+\tfrac12 V^4,\qquad
                  \tfrac12 V+\tfrac12 V^3+\tfrac1{10}V^5 .
\]
Equal values of the first imply either \(V_2=-V_1\) or \(V_2^2+V_1^2=-4\) off the diagonal. In the first case \(z=V_1^2\) satisfies \(z^2+5z+5=0\); these give four simple ordered solutions, none with \(z=-2\). In the second case put \(T=V_1V_2\). The divided difference of the second polynomial reduces to \((1+T-T^2)/10\). Each of its two roots gives two possible sums \(V_1+V_2\) and two orderings. None has \(V_1=\pm V_2\), and the Jacobian of \((V_1^2+V_2^2,V_1V_2)\) is nonzero there, so all eight solutions are simple. The same sign quotient gives six branches of index two with exactly two fixed by exchange.

In every row the number of scaled solutions equals the local length. Away from the diagonal the exchange has no fixed point, so the listed local counts sum to the asserted global counts at \(r=0\).
