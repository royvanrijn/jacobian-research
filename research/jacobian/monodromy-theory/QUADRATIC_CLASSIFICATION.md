# Complete quadratic primitive-pair reduction

Work over algebraically closed characteristic-zero constants.
As in KMFI48--53, write
\[
 P=A\begin{pmatrix}-bc&b^2\\-c^2&bc\end{pmatrix},\qquad
 \gcd(b,c)=1,\quad A\ne0,
\]
and impose
\[
 W(A,b,c)=(Ab^2)_{ww}+2(Abc)_{wz}+(Ac^2)_{zz}=0.
\]
The primitive pair has max(deg b,deg c)=2. Its scalar A has
**arbitrary finite degree**.

**KMFI54.** For every integer m>=0 the weighted operator
T_m(kappa) of KMFI52 is injective for every kappa.
Its top maximal minor has no real roots; its bottom maximal minor
has only simple real roots. Both assertions have all-weight proofs.
The earlier weight0..20 Bezout certificates remain independent
finite regression evidence, rather than the reason for this theorem.

**KMFI55.** Every primitive quadratic wave pair, with scalar of
any degree, is affinely conjugate to the completed k(w)[z] chart.
Its full geometric field action is natural alternating, natural
symmetric, or binary SL2(F3) on eight letters.
Both PSL2(11) actions and natural A4 are excluded throughout this
entire quadratic class, including independent stabilization and
birational recharts of the same field pair.

This is a classification of the specified sweep fields. It is not an
absolute exclusion for every polynomial Keller map, and it does not
provide an affine Keller realization of the binary field.

[Continuant helper](quadratic_continuant.py) ·
[exact checker](verify_quadratic_continuant.py) ·
[model](quadratic_continuant_model.json) ·
[receipt](quadratic_continuant_receipt.json).
The earlier [quadratic scalar note](QUADRATIC_SCALAR_OBSTRUCTION.md)
retains the finite certificates and scalar-linear coefficient audit.

## 1. The bottom tridiagonal matrix

Fix m>=0, n=floor(m/2)+1 and i_j=m-2j, for j=0,...,n-1.
KMFI52 gives T_m with n columns and n+1 rows:
column j has D_j on row j, E_j on row j+1, and F_j on row j+2,
where
\[
\begin{split}
 D_j&=(i_j+4)(i_j+3)+2\kappa(i_j+3)(j+1)
                     +\kappa^2(j+1)(j+2),\\
 E_j&=2(i_j+1)(i_j+2+\kappa(j+2)),\\
 F_j&=i_j(i_j-1).
\end{split}                                                    \tag{1}
\]
The final F_j is zero since i_(n-1) is0 or1.
Delete the first row. The square bottom matrix B_m is tridiagonal:
\[
 (B_m)_{j,j}=E_j,\quad (B_m)_{j,j+1}=D_{j+1},\quad
 (B_m)_{j+1,j}=F_j.                                             \tag{2}
\]
Let H_l be the determinant of its leading l-by-l submatrix.
Then H_0=1,H_1=E_0, and
\[
 H_l=E_{l-1}H_{l-1}-F_{l-2}D_{l-1}H_{l-2}\quad(2\le l\le n).
                                                                  \tag{3}
\]
In every adjacent edge i=i_j>=2 and
\[
\begin{split}
 D_{j+1}&=(i+2)(i+1)+2\kappa(i+1)(j+2)
                      +\kappa^2(j+2)(j+3),\\
 \operatorname{disc}_{\kappa}D_{j+1}
       &=-4(i+1)(j+2)(i+j+4)<0. 
\end{split}                                                       \tag{4}
\]
Its leading coefficient is positive, so D_(j+1)(x)>0 for every
real x. Also F_j>0. These strict inequalities hold at every
edge of every m, including the final even and odd boundary edges.

## 2. Every leading coefficient is positive

Equation(3) shows deg H_l<=l. The potential coefficient of kappa^l
is the leading principal determinant of the real tridiagonal matrix
L_m with entries
\[
 a_j=2(i_j+1)(j+2),\quad u_j=(j+2)(j+3),\quad
 v_j=i_j(i_j-1).                                                  \tag{5}
\]
This also follows immediately by applying recurrence(3) to the
highest-degree coefficients, since deg E_j=1,deg D_j=2.
All a_j,u_j,v_j are strictly positive where present.
A positive diagonal similarity symmetrizes L_m. After dividing its
diagonal to1 by positive diagonal congruence, the square of the
off-diagonal entry on edge j is
\[
 q_j={u_jv_j\over a_ja_{j+1}}
      ={i_j\over4(i_j+1)}<\frac14.                                \tag{6}
\]
Its normalized LDL pivots satisfy
\[
 p_0=1,\qquad p_{j+1}=1-q_j/p_j>\frac12.
\]
Induction proves the displayed strict bound: p_j>1/2 implies
q_j/p_j<1/2. Hence every leading principal determinant of L_m
is positive. It follows that every H_l has **exact degree l**
and positive leading coefficient.

This step is essential: positivity of the edge factors alone does
not ensure the extra two outer roots in the next step.

## 3. Elementary all-degree interlacing

We prove by induction that H_l has l simple real roots and that
they strictly interlace those of H_(l-1).
H_0=1 and H_1=2(m+1)(m+2+2kappa) give the base.

For l>=2, let x_1<...<x_(l-1) be the roots of H_(l-1).
Strict interlacing at the previous step makes H_(l-2)(x_a)
nonzero with alternating signs. By(3),(4),
\[
 H_l(x_a)=-F_{l-2}D_{l-1}(x_a)H_{l-2}(x_a).
\]
Thus H_l changes sign on each of the l-2 internal intervals.
At the rightmost x_(l-1), H_(l-2) is positive because it has
positive leading coefficient and all its roots lie to the left.
So H_l is negative there and positive at positive infinity:
there is one further root on the right.
At the leftmost x_1 the sign of H_(l-2) is (-1)^(l-2);
the sign of H_l there is (-1)^(l-1), whereas its sign at
negative infinity is (-1)^l. There is one further root on the left.
These l distinct real roots exhaust deg H_l=l. All are simple,
and their locations give strict interlacing.
This proves the induction without any bounded root computation.

The bottom maximal minor is H_n and therefore has only real roots.
The top maximal minor is the triangular product of the D_j.
For all its factors,
\[
 \operatorname{disc}_{\kappa}D_j
       =-4(i_j+3)(j+1)(i_j+j+5)<0.                               \tag{7}
\]
They have positive leading coefficient and no real roots.
Thus top and bottom maximal minors have no common complex root,
and they are coprime over Q[kappa].
At each complex kappa at least one maximal minor is nonzero,
so T_m(kappa) is injective.
The same conclusion holds over every characteristic-zero extension:
the two rational polynomials are coprime, hence have a rational
Bezout identity (even though no all-m explicit identity is needed).
This proves KMFI54.

## 4. Complete reduction and field actions

KMFI49 leaves two leading-direction cases for a primitive quadratic
pair. KMFI51 closes the constant-leading-direction case for all
scalar degrees. In the other case the possible obstruction is
\[
 b=w^2+\beta z+b_1w+b_0,\qquad
 c=\kappa wz+c_1w+c_2z+c_0.
\]
If beta!=0, weights(w,z)=(1,2) give, after scaling, the top
pair (w^2+z,kappa wz). A nonzero scalar's highest weighted component
would be in ker T_m(kappa). KMFI54 makes that kernel zero for every
m. Lower-weight pair or scalar terms cannot cancel the highest
wave. Therefore beta=0.
Both cases now have b depending only on w, so rho=c/b belongs
to k(w)[z]. KMFI46 applies. In this quadratic-pair subclass a
polynomial rho has degree at most2, so its higher even signed
V_n:A_n exceptions (which require polynomial degree at least4)
are absent. Its binary degree-eight exception remains.

For clarity, one can give the full unbounded action families.
All normalizations below allow the harmless affine terms already
removed in KMFI33,41,44. Write
\[
 C(r,s)=r(r-1)+2\kappa(r-1)(s+1)+\kappa^2(s+1)(s+2)=0.
\]
The following table imposes no scalar-degree cutoff:

| reduced primitive pair | parameters | full degree and geometric action |
|---|---|---|
| (1,h(w)), deg h=2 | heat top z-degree s>=0, epsilon=0 or1 | epsilon0: N=(s+1)(3s+5); epsilon1: N=(s+2)(3s+4); A_N, with the exception below |
| (w,w^2+kappa z) | r>=2,s>=0,C(r,s)=0 | N=(s+2)(r+2s+2) |
| (w^2,kappa wz+1) | r>=s+4,s>=0,C(r,s)=0 | N=(s+2)(r-s-1) |

In the last two rows G=S_N exactly when r is even and s is odd;
otherwise G=A_N.
In the first row the sole exception is s0,epsilon1:
the group is binary SL2(F3) on8 precisely when h is even about
the root of the affine lambda; otherwise it is A8.

Here is the case exhaustion, without extrapolating the old finite
heat-degree audit. Constant b gives exactly the quadratic polynomial
heat family, and KMFI44--45 give its degrees and groups for all s.
A nonconstant linear b has a single simple root. Primitivity and
KMFI33 make the z-dependent ratio a single simple pole, and the
nonaffine R is quadratic; this is the second row after affine
normalization. If b has degree2, the leading pair is diagonal.
The ratio's single simple pole and primitivity force b to be a
square of that pole's linear factor, with c=kappa wz+c_1w+c_0
and c_0!=0 after translation. This gives the third row after
removing affine R, and KMFI41's sharp pole bound is r-s>=4.
No z-independent rational ratio with nonconstant b remains:
KMFI33 makes such a ratio polynomial, contradicting primitivity.
Thus every wave pair lies in a displayed row.
Conversely each row has polynomial scalar witnesses from KMFI41
or44: lambda itself for constant b, lambda/w^2 for the second
row, and lambda/w^4 for the third row. The sharp bound in the
third row ensures polynomiality of the latter quotient.

Every displayed degree is composite except N=5 in the first row
s0,epsilon0. Thus the only prime-degree field action in this
entire primitive quadratic class is natural A5. In particular,
the degree-eleven PSL2(11) action is excluded already by degree;
the degree-twelve action is excluded by the full group classification.

The proofs concern the indicated field pairs. Independent
transcendental stabilization and birational recharts preserve
degree and monodromy action, so they cannot evade the exclusions.
A field witness still needs an actual polynomial affine Keller
completion to become an occurrence theorem.

## 5. Assurance and remaining goal

The universal proof consists of the explicit continuant recurrence,
strict positivity of its leading principal coefficients, and the
elementary interlacing induction. The checker verifies the generic
symbolic identities and exact finite controls, including weights
beyond the old bound. Finite controls do not replace the induction.
No formal verification or independent whole-proof review is claimed.

The remaining rational-pole frontier starts at primitive-pair degree
at least3, subject to KMFI49's necessary leading-tensor constraints.
Other direction fields, higher architectures and independent Shabat
completion remain possible routes. [KMFI60](BRANCH_CENTERED_COMPLETION.md)
now realizes the binary degree8 action. Absolute PSL2(11), natural A4
and the full action spectrum remain open.
