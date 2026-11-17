# Nilpotent normal derivatives and cubic multiplier obstructions

The zero-stratum obstruction now excludes every nonnilpotent normal
derivative, including singular ones. It also specifies the leading
multiplier when the normal derivative is a regular nilpotent block.
Together with a general polynomial translation test, this rejects
an explicit primitive cubic direction that passes both the
fourth-power Jacobian and cyclicity gates.

These results do not classify the remaining cubic directions or decide
the absolute exceptional11 action.

## 1. Complete linear feedback nilpotency test

Let B be any c by c matrix over a characteristic-zero field K, q
independent variables, and J in K[q] a nonzero polynomial. Put
\[
                       A=B+(Bq)d\log J.                  \tag{1}
\]

**KMFI158.** Matrix A is nilpotent if and only if B is nilpotent and
\[
                 dJ\,B^j q=0\quad(1\le j\le c-1).         \tag{2}
\]
No homogeneity assumption on J is needed.

The determinant lemma gives
\[
 d\log J\,(TI-B)^{-1}Bq
                   =1-\frac{T^c}{\det(TI-B)}              \tag{3}
\]
if A is nilpotent.
Use the Fitting decomposition K^c=E+K0, with B invertible on E
and nilpotent on K0. This decomposition is defined over K itself.
Suppose k=dimE>0, with l=c-k. Then det(TI-B)=T^l*chi_E(T),
chi_E(0)!=0. The right side of(3) is regular at T0 and tends to1.

The K0 part of the resolvent in(3) consists only of negative powers
of T, so those coefficients must vanish. The E part tends to -q_E.
Consequently
\[
                          dJ\,q_E=-J.                   \tag{4}
\]
This is an Euler eigenvalue -1 on the polynomial variables belonging
to E, impossible in characteristic zero. Hence B is nilpotent.
Now det(TI-B)=T^c and(3)'s right side is zero.
The finite nilpotent resolvent gives exactly(2).
Conversely(2) makes the determinant in(1) equal T^c,
hence A is nilpotent by Cayley--Hamilton.

If B is a regular nilpotent block, choose coordinates with
Bq=(q2,...,qc,0). The last nonzero condition in(2) first forces
J_q1=0; descending then forces J_q2=...=J_q(c-1)=0.
Thus J belongs to K[qc]. If it is homogeneous of degree m, then
\[
                              J=\kappa q_c^m.             \tag{5}
\]
The feedback for this J is again regular nilpotent: the last
superdiagonal factor is m+1, nonzero for every integer m>=0.

The polynomial condition is essential. B=diag(1,0) and J=1/q1
give zero feedback despite B being nonnilpotent.
Trace invariance alone is insufficient: for the regular nilpotent
3-block and J=2q1*q3-q2^2, dJ*Bq=0 but
dJ*B^2q=2q3^2, and(1) is not nilpotent.

## 2. All zero-stratum normal derivatives must be nilpotent

Keep the linear zero-stratum hypotheses of KMFI156:
p=(q,t), Z={q=0}, V|Z=h|Z=0.
Let B(t)=D_qV_q(0,t), now allowed to be singular.

**KMFI159.** If a nonzero polynomial F makes
N=DV+V*dlogF-hI nilpotent, then B is nilpotent over K=k(t).
For the least normal degree m and initial form J of F,
all the moment invariances(2) hold.
If B is regular nilpotent, J is a pure m-th power of its left-kernel
linear coordinate, up to a nonzero coefficient in k(t).

The exact normal rescaling from KMFI156 has limiting normal block
B+(Bq)d_qlogJ. It is nilpotent. Apply KMFI158.
The upper-right tangential correction remains present; both lower
blocks vanish. No argument from the full determinant alone can give
this refinement when B is singular.

There is an immediate all-degree homogeneous base-point test.
For any normalized full-rank homogeneous input H=FV with primitive
V of degree d>=2, every nonzero a with V(a)=0 satisfies
\[
                                DV(a)\text{ nilpotent}.   \tag{6}
\]
Indeed Euler gives DV(a)*a=dV(a)=0, and KMFI125's determinant identity
gives h(a)=0. The whole radial line k*a is a zero stratum of V,h.
In linear coordinates with that line tangential, the derivative has
block form [[B,0],[C,0]]. Its normal block B is nilpotent by the
just-proved result, so the whole derivative is nilpotent.
This applies at every base point, including singular points of the
projective base locus; the radial line itself is always linear.

The pure-power conclusion concerns only the initial normal form.
For example F=u+y^2 has initial form u along u=y0, but is not divisible
by u. No global divisibility of F is inferred from(5).

## 3. A general polynomial translation coefficient

Let a,b be independent constant vectors, and P a polynomial with
partial_bP=0 and partial_aP!=0. Suppose
\[
                        \partial_bV=P\,a,\qquad
                        \partial_bR=0.
\]

**KMFI160.** No nonzero polynomial F satisfies div(FV)=R*F.

Take coordinates a=partial_w,b=partial_z, so P is independent of z
and V=z*P*e_w+V0 with V0 independent of z.
For J the highest z coefficient of F, the highest coefficient of the
equation is partial_w(PJ)=0.
But P has positive w degree, so the product of P with a nonzero
polynomial J cannot have w degree zero. This is the contradiction.

This includes KMFI155 and all powers P=u^k, k>=1, when u(b)=0
and u(a)!=0. It also applies to nonlinear P with mixed coefficients.
It does not apply when partial_aP=0. Rational 1/P controls again
show why polynomiality is essential.

## 4. An exact primitive cubic direction is excluded in every scalar degree

Set
\[
 K=(u^3,u^2x,u^2y+ux^2,u^2z+xy^2),\qquad V=A K,\quad A\in GL_4.
                                                               \tag{7}
\]
Every such direction is primitive and detDV=3det(A)u^8.
A full-rank normalized homogeneous multiplier of degree s would have
h=lambda*u^2 and lambda^4=(s+3)detA.

**KMFI161.** Every such normalized multiplier requires
\[
                             A02=A03=A13=0.               \tag{8}
\]
These are necessary matrix conditions, not a classification of the
remaining matrices or scalar multipliers.

The plane Z1={u=y=0} is a zero stratum of V,h. Its normal derivative
in coordinates(u,y) is
\[
                    B_1=x^2\begin{pmatrix}A02&0\\A22&0\end{pmatrix}.
\]
Nilpotency by KMFI159 forces A02=0.
On Z2={u=x=0}, the normal derivative in coordinates(u,x) is
\[
                    B_2=y^2\begin{pmatrix}0&A03\\0&A13\end{pmatrix},
\]
so nilpotency forces A13=0.
Both normal blocks have identically zero determinant; the old
invertibility obstruction alone would not detect either eigenvalue.

Finally partial_zV=u^2*A*e_z and h is independent of z.
KMFI160 forces A03=0.
If A22!=0 after these restrictions, B1 is regular nilpotent.
The lowest (u,y)-homogeneous part of F must then be
kappa(x,z)*u^m over k(x,z), by KMFI159.

For a concrete all-scalar exclusion, take A swapping output rows0,2:
\[
                  V=(u^2y+ux^2,\ u^2x,\ u^3,\ u^2z+xy^2). \tag{9}
\]
Then detDV=-3u^8, and the Krylov determinant at p=(1,1,1,1) is1.
Thus both the perfect fourth-power Jacobian and generic cyclicity
gates pass. The last-coordinate translation has a=e_z and
partial_a(u^2)=0, so that particular translation test does not exclude it.
But at a=(0,1,0,0), V(a)=0 and
DV(a)=diag(1,0,0,0), violating(6).
No homogeneous polynomial scalar of any degree can make(9) a normalized
full-rank input; in fact the same zero-stratum obstruction excludes
any polynomial scalar with h=lambda*u^2.

Normalized lower-rank examples remain valid inside the retained matrix
conditions: A=diag(1,lambda,lambda,lambda), lambda=s+3,F=u^s,
for s>=0 give
\[
 f=-\lambda u^{s+2},\quad
 g=(-(s+2)u^{s+3},0,\lambda u^{s+1}x^2,\lambda u^sxy^2).
\]
The entire sweep Jacobian is gamma^4, while X is independent of z
and has rank3. These controls do not establish full-rank existence
for any remaining matrix in(8).

## 5. Evidence and the unresolved action objective

[Helper](nilpotent_normal.py) · [checker](verify_nilpotent_normal.py) ·
[model](nilpotent_normal_model.json) · [receipt](nilpotent_normal_receipt.json).

The exact identities supplement the all-field Fitting/Euler proof,
normal-limit proof, regular-block initial-form classification and
general translation proof. They also retain rational and
trace-only countercontrols, both entire cubic normal derivatives,
the exact cyclic excluded cubic, and actual normalized lower-rank
whole sweep Jacobians.
Formal verification, independent whole-proof replay and external
review are not claimed.

The previous goal turn was progress: KMFI155--157 closed every
primitive quadratic rank4 direction. This turn strengthens the
all-degree normal and translation tests and excludes an explicit
cyclic primitive cubic direction in every scalar degree.

The arbitrary cubic direction problem is not solved. Even family(7)
with(8) has unresolved full-rank scalar tests and generic groups.
For cubic V, all quadratic shapes of the fourth-power base h0 remain
in scope, including irreducible quadrics and reducible cases.
Higher primitive directions, genuinely nonhomogeneous fields,
actual exceptional11 or an absolute action obstruction, and every
remaining faithful-transitive Keller action still require resolution.
The user's complete goal remains active.
