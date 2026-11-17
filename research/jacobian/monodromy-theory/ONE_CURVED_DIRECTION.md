# The one-curved-coordinate architecture in every rank and degree

The fixed-direction multiplier criterion closes a substantially larger
class than the quadratic binary shear. In every rank and every homogeneous
degree, arbitrary polynomial curvature in one last projective coordinate
is incompatible with full rank. The only surviving inputs reduce to the
scalar-times-linear class already classified by KMFI121--122.

This is a construction obstruction, not an absolute action obstruction.
It does not classify arbitrary primitive directions or decide exceptional11.

## 1. The full architecture and conclusion

Work over algebraically closed characteristic-zero constants. Let r>=2,
m=r-2, p=(u,x1,...,xm,t), and d>=1. For an arbitrary polynomial Q(x)
of total degree at most d, let Qh(u,x)=u^d Q(x/u) and set
\[
 V=(u^d,\ a_i u^{d-1}x_i\ (1\le i\le m),\
                    c u^{d-1}t+Qh(u,x)).                 \tag{1}
\]
The constants a_i,c are unrestricted initially; Q can have mixed terms,
all lower projective degrees, linear terms and a constant term.
We do not assume V primitive: this makes the surviving case precise.

**KMFI127.** Suppose a nonzero homogeneous polynomial F and polynomial
f,g give H=g-fp=FV, a normalized sweep detDS=gamma^r, and rankDX=r
for X=(f,g). Then Q is affine:
\[
                       Q=q_0+\sum_i q_i x_i.              \tag{2}
\]
Moreover c and the a_i are nonzero, the a_i are pairwise distinct
and different from1, and c is different from1 and every a_i.
Consequently
\[
 V=u^{d-1}Mp,\quad
 M=\begin{pmatrix}
 1&0&0\\
 0&\operatorname{diag}(a_i)&0\\
 q_0&(q_i)&c
 \end{pmatrix}.                                           \tag{3}
\]
Thus every full-rank member belongs exactly to the old factorized
architecture, with scalar G=F u^{d-1}. Its primitive coordinate degree
is1. In particular no primitive degree>=2 direction of form(1)
admits a full-rank multiplier, in any scalar degree.
Every degree11 field in this architecture is S11,A11 or regularC11
by KMFI122; neither exceptional PSL2(11) action occurs here.

The result includes all resonant cases a_i=d or c=d, repeated parameters
and arbitrary numbers of curved monomials. It does not assert that every
affine Q admits a multiplier or that no nonlinear Q admits a normalized
rank-deficient multiplier.

For a completely explicit description of the surviving multipliers,
put z_i=x_i and
\[
 z_0=u,\qquad z_t=t+\frac{q_0}{c-1}u+
                            \sum_i\frac{q_i}{c-a_i}x_i.
\]
Then M is diagonal with eigenvalues b=(1,a_1,...,a_m,c).
By KMFI121 every possible scalar is
\[
 G=\kappa u^{B_0}\prod_i x_i^{B_i}z_t^{B_t},
 \qquad B_j\in\mathbb Z_{>0},\qquad B_0\ge d-1,             \tag{4}
\]
\[
 F=G/u^{d-1},\quad f=-\lambda G,\quad
                         g=G(M-\lambda I)p,              \tag{5}
\]
where lambda is nonzero, differs from all b_j, and
\[
 \det(TI-[\operatorname{diag}(b)+b B^t])=(T-\lambda)^r.     \tag{6}
\]
Conversely (2)--(6) give a full-rank normalized polynomial input.
The positivity in(4) and the divisor inequality are essential.

## 2. An arrow matrix and identities independent of d

By KMFI125 (or its same proof without removing the coordinate gcd),
f=-Fh and detDV=d/(s+d)*h^r, where s=degF.
The derivative determinant of(1) is
\[
            \det DV=d c\prod_i a_i\ u^{r(d-1)}.
\]
Full rank forces c and all a_i nonzero, and
h=lambda*u^(d-1), lambda nonzero.
Set eta=dlogF. Then the characteristic identity is
\[
                  \det(TI-[DV+V\eta])=(T-h)^r.             \tag{7}
\]
Its Krylov matrix is invertible by KMFI124--125. These arguments
apply without primitivity: the determinant gives f/F polynomial,
and the Euler identity and full-rank cyclicity prove the same Krylov test.

On u=1 write A=DV|_(u=1), v=V|_(u=1), and eta_i for the transverse
components of eta there. The matrix A is an arrow:
\[
 A_{00}=d,\quad A_{i0}=(d-1)a_i x_i,\quad A_{ii}=a_i,
\]
\[
 A_{t0}=(d-1)ct+dQ-\sum_i x_iQ_i,\quad
                     A_{ti}=Q_i,\quad A_{tt}=c.           \tag{8}
\]
All unlisted entries vanish. Write
\[
 \chi=(T-d)(T-c)\prod_i(T-a_i),\quad
 K=c(c-1)\prod_i(c-a_i),
\]
\[
 R=c\left[Q\prod_i(c-a_i)+
       \sum_i(a_i-1)x_iQ_i\prod_{j\ne i}(c-a_j)\right].     \tag{9}
\]
Solving the arrow system gives the entire rational identity
\[
 [(TI-A)^{-1}v]_0=\frac1{T-d},
\]
\[
 [(TI-A)^{-1}v]_i
       =\frac{a_i x_i(T-1)}{(T-d)(T-a_i)},
\]
\[
 [(TI-A)^{-1}v]_t=
 \frac{ct(T-1)+TQ+
       \sum_i T(a_i-1)x_iQ_i/(T-a_i)}{(T-d)(T-c)}.          \tag{10}
\]
Multiplication by chi clears the displayed denominators. Hence these
are polynomial adjugate identities, valid also when a_i=d or c=d.
In particular
\[
 [\operatorname{adj}(cI-A)v]_j=0\quad(j\ne t),\qquad
                     [\operatorname{adj}(cI-A)v]_t=Kt+R.  \tag{11}
\]
For the Krylov matrix C=(v,Av,...,A^(r-1)v), the same identity gives
\[
 \det C=\left[\prod_i a_i(a_i-1)x_i\right]
       \left[\prod_{i<j}(a_j-a_i)\right]\,(Kt+R).            \tag{12}
\]
One proof for distinct d,a_i,c expands(10) into simple fractions:
C is the eigenvalue Vandermonde times the corresponding residue
matrix. Its factors d-a_i and c-d cancel. Both sides of(12) are
polynomials, extending the identity to all parameters.

Thus C invertible forces the a_i distinct, none equal1, and Kt+R
nonzero. No division by a_i-d or c-d is needed.

## 3. Polynomial logarithmic residues force an affine Q

Evaluate the determinant lemma in(7) at T=c on u=1. Equations(11) give
\[
                    \eta_t=-\frac{(\lambda-c)^r}{Kt+R}.  \tag{13}
\]
If lambda=c, the polynomial F has F_t=0. Then f=-lambda*F*u^(d-1)
and every coordinate of g=F(V-hp) is independent of t, so rankDX<r.
This is impossible.

If K=0, the right side of(13) is constant in t.
A nonzero polynomial's logarithmic derivative in t tends to0 at
infinity. Thus lambda=c, again impossible. This excludes c=1
and c=a_i, including cases where the Krylov matrix itself is invertible.

Now K!=0 and lambda!=c. Put S=R/K and q=t+S(x). Polynomiality gives
\[
 F(1,x,t)=J(x)q^\nu,\qquad
                 \nu=-(\lambda-c)^r/K\in\mathbb Z_{>0}.   \tag{14}
\]
Evaluation of(7) at T=a_i, using(10), gives
\[
 \eta_i+\frac{Q_i}{a_i-c}\eta_t
 =-\frac{(a_i-\lambda)^r}
 {a_i(a_i-1)x_i(a_i-c)\prod_{j\ne i}(a_i-a_j)}.            \tag{15}
\]
This identity also holds at a_i=d after clearing denominators.
Compare the poles at q=0 in(14)--(15):
\[
                              S_i=Q_i/(c-a_i).            \tag{16}
\]
Since the a_i are distinct, equality of mixed derivatives implies
Q_ij=0 for every i!=j. Hence
\[
                   Q=q_0+\sum_i P_i(x_i),\quad P_i(0)=0,
\]
\[
              S=q_0/(c-1)+\sum_i P_i(x_i)/(c-a_i).         \tag{17}
\]
The constant in S follows directly from(9).

Finally the trace of(7) says
\[
 V(F)/F=[r\lambda-d-c-\sum_i a_i]u^{d-1}.
\]
For a homogeneous F of degree s, its expression on u=1 has
projective derivative
\[
 {\cal W}=\sum_i(a_i-1)x_i\partial_i+
                                ((c-1)t+Q)\partial_t.
\]
Its logarithmic derivative under W is constant; the radial contribution
is s. Since J has no t dependence, the pole at q=0 in W(logF)
must vanish. The t coefficient determines the cofactor c-1, so
\[
                Q+\sum_i(a_i-1)x_i S_i-(c-1)S=0.          \tag{18}
\]
For the coefficient q_ik*x_i^k of P_i, (17)--(18) give exactly
\[
                       q_{ik}\frac{(k-1)(a_i-1)}{c-a_i}=0.
\]
For k>=2 it forces q_ik=0. Thus Q is affine.
For m=0 this conclusion is automatic; the same equations apply.
This proves(2), and(3)--(6) then follow from KMFI121.

## 4. What is closed, and what remains

This closes every homogeneous single-curved-coordinate input(1),
in every rank, primitive direction degree and scalar multiplier degree.
The quadratic binary shear KMFI126 is its degree2, rank4 special case.
Here Q may also contain terms involving u after homogenization.

The surviving affine-Q cases have genuine full-rank controls; they
must not be discarded merely because their displayed V has degree d.
Removing its coordinate gcd recovers degree1, not a counterexample
to the primitive obstruction.

The theorem does not cover multiple curved output coordinates,
nonlinear couplings in the earlier projective coordinates, arbitrary
nontriangular primitive directions, or nonhomogeneous matrix fields.
It does not supply a new abstract action, a PSL2(11) input, or an
absolute exceptional-action obstruction. The full faithful-transitive
Keller monodromy classification remains unresolved and active.

[Helper](one_curved_direction.py) · [checker](verify_one_curved_direction.py) ·
[model](one_curved_direction_model.json) ·
[receipt](one_curved_direction_receipt.json).

The unbounded proof is written above. Exact polynomial identities,
resonant controls, logarithmic consistency and positive controls supplement
it. Formal verification, independent whole-proof replay and external
review are not claimed.
