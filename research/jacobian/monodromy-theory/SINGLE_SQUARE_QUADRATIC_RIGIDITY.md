# The entire single-square quadratic shear is excluded to full rank

The arbitrary-matrix single-square family is closed in every homogeneous
scalar degree. Its remaining unmatched case has curvature pointing
outside the critical hyperplane. That condition excludes a fixed
nonlinear hypersurface, then excludes both ways a nonlinear scalar
divisor could be ruled. Linear factors contradict a differentiated trace.

## 1. An elementary fixed-hypersurface lemma

Let n>=3, p in k^n over algebraically closed characteristic-zero constants,
M invertible, u and ell independent linear forms, and
\[
 V=uMp+a\,\ell^2,\qquad Mb=a,\quad u(b)=0,\quad u(a)\ne0.    \tag{1}
\]

**KMFI153.** The components of p wedge V have no common irreducible
factor of degree>=2. In particular the projective map [V] cannot fix
pointwise a nonlinear irreducible hypersurface. The stronger polynomial
statement also covers a hypersurface contained in the base locus.

Here is an elementary polynomial fixed-locus fact. If a homogeneous
quadratic map W fixes an irreducible degree-m hypersurface D pointwise,
then D divides every component of p wedge W.
If m>=3 these cubic components vanish identically: for m>3 use degree;
for m=3 their quotients are constant two-forms B, and
p wedge B=0 for all p forces B=0 in dimension>=3.
Thus W=L*p for a linear L.
If m=2, write p wedge W=D*B(p), with B a linear two-form.
The identity p wedge B=0 makes B=p wedge c for a constant c.
Indeed the coefficient of p_k^2 in the triple-index identity removes
every variable outside its two indices, and the remaining cross
coefficients identify B_ij=c_j*p_i-c_i*p_j.
Then p wedge(W-Dc)=0, hence
\[
                              W=L p+D c.                 \tag{2}
\]
This is polynomial coefficient algebra, with no normality assumption.

Apply this fact to(1). A radial W=L*p is impossible: on u0,
V=a*ell^2, with a outside that hyperplane, whereas L*p is radial in it.
In the quadratic case, restrict(2) to u0:
\[
                            a\ell^2=L_0p+D_0c.
\]
The hyperplane keru is not contained in span(a,c).
For n>=4 this follows from dimensions; for n3 equality of those two
planes would put a in keru. Modulo span(a,c), the last identity forces
L0=0. Since D is not u, D0 is nonzero; the identity makes c parallel a
and D0 a nonzero multiple of ell^2.
After rescaling D, write
\[
                 L=\kappa u,\quad D=\ell^2+uR,\quad c=a.
\]
Consequently M=kappa*I+a*r, where R(p)=r*p.
Apply u to Mb=a: u(b)=0 and u(a)!=0 give r*b=1.
Thus Mb=a also gives kappa*b=0, hence kappa=0.
This makes M rank<=1, a contradiction. The lemma follows.

The transverse-curvature hypothesis u(a)!=0 is essential:
V=(u^2,ux,uy,2uz+2x^2) fixes the irreducible quadric uz+2x^2
pointwise. It has M=diag(1,1,1,2), b=e_z and a=2e_z, with u(a)=0.
Thus all the other hypotheses hold, and a pointwise-fixed nonlinear
divisor cannot simply be omitted when transversality is dropped.

## 2. A single-square scalar cannot have nonlinear factors

Let
\[
 K_C=(u^2,ux,uy,uz+x^2),\quad V=A K_C,\quad A\in GL_4.
\]
Equivalently,
\[
 V=uMp+a x^2,\quad M=A,\quad b=e_z,\quad a=Mb,\quad
 B=\mathbb P(\ker u\cap\ker x).                            \tag{3}
\]
Here B is the entire projective base line.
Assume a homogeneous scalar F gives a normalized full-rank rank4 input.
KMFI150 gives u not dividing F. If A03=u(a)=0, then V_u=uL.
The logarithmic form dlogF is regular at generic u0; row u of
N=DV+V*dlogF-lambda*u*I there is (L|u0,0,0,0).
Nilpotency makes L|u0=0. Thus V_u=c*u^2 with c!=0,
the matched primitive case excluded by KMFI141.
It remains to treat u(a)!=0.

**KMFI154.** No homogeneous multiplier of any degree produces full rank
for this arbitrary-matrix single-square family.

First exclude any nonlinear irreducible factor D of F.
The global divisor homotopy KMFI144 puts the line joining p,V(p) in D.
KMFI153 excludes a pointwise-fixed D, so the joining line is generically
defined. The selected-line-family argument of KMFI140 gives a curve
of preserved projective lines with one selected line through a general
point of D. KMFI142 makes the restriction of [V] to each general line
degree1. Its quadratic sections have one common base point q in B.

Because u=x=0 at q, dividing the restricted sections by u gives at q
the projective value Mq: along the line x=t*u, so x^2/u=t^2*u vanishes.
Thus each selected line L contains q and Mq.

If q varies on B, then generically u(Mq)!=0, since u(Mb)=u(a)!=0.
Therefore L is the unique span(q,Mq).
If B and MB meet, all these lines lie in their projective plane span,
contrary to the nonplane D. Otherwise the two lines are skew, and
this corresponding-point family is a smooth ruled quadric.
There is a unique q0 in B with u(Mq0)=0.
The corresponding line L0=span(q0,Mq0) lies in u0.
Since the two directrices are skew, Mq0 is not in B,
so x(Mq0)!=0. At its nonbase points, V sends L0 to
a*x^2, outside u0. This contradicts preservation of L0.

Preservation holds also at this limiting line: the polynomial identities
q wedge Mq wedge V(alpha*q+beta*Mq)=0 hold for general q and all
alpha,beta, hence specialize to q0. No properness assumption on the
affine map or unproved extension of its inverse is used.
In coordinates, put Delta=M02*M13-M03*M12. Skewness is Delta!=0,
q0=M03*e_y-M02*e_z, and x(Mq0)=-Delta.
The offending u value is M03*Delta^2, nonzero.

If instead q is a fixed point w in B, the surface D is a cone with
vertex w. Necessarily Mw is parallel to w; otherwise all selected
lines would equal the single span(w,Mw).
Since M is invertible, it induces an invertible map on k^4/k*w.
The forms u,x descend independently, and so do a and b, with
Mbar*bbar=abar, u(bbar)=0 and u(abar)!=0.
The cone equation descends to a nonlinear irreducible curve in P2.
Because each selected line through w is preserved, the quotient
projective map fixes that curve pointwise. This contradicts KMFI153
in dimension3. The vector b cannot disappear in the quotient:
if b were parallel w, Mb=a would lie in u0.

These two cases exhaust the irreducible selected-line family.
Hence every factor of F is linear. This conclusion does not use the
individual factor-degree bound3; it applies in every factor degree.

## 3. Linear factors contradict the trace

For a linear factor D of F, weighted preservation gives
\[
 D(V)=c_D*u*D,\quad c_D\ne0.
\]
Compare x^2 terms in(3), then cancel u:
\[
                         D(a)=0,\qquad D M=c_D D.
\]
Since Mb=a, also D(b)=0. Therefore the whole scalar satisfies
\[
                              \partial_a F=\partial_b F=0. \tag{4}
\]
This includes a constant scalar if it had not already been excluded.

The complete trace equation is
\[
                        V(F)+F\,\operatorname{div}V
                                      =4\lambda uF.       \tag{5}
\]
Here u(b)=x(b)=0 gives
\[
 \partial_bV=u a,\qquad
 \operatorname{div}V=u\operatorname{tr}M+u(Mp)+2x\,x(a),
 \qquad \partial_b\operatorname{div}V=u(a).
\]
Differentiate(5) in the constant direction b.
By(4), its left side is u(a)*F and its right side is zero.
Since u(a)!=0 and F!=0, this is impossible.
The entire arbitrary-A family K_C is therefore excluded to full rank.

The full-rank qualification is indispensable, even with F coprime to u.
Let beta^2-beta+1=0, lambda=1+beta, F=y and
A=diag(1,lambda,beta,lambda).
Then f=-lambda*u*y and
\[
 g=((1-\lambda)u^2y,\ 0,\ -uy^2,\ \lambda x^2y)
\]
give a normalized input with whole sweep Jacobian gamma^4.
Its X is independent of z and has rank3:
the minor using rows f,g_y,g_z and columns u,y,x is
2lambda^2*u*x*y^3, nonzero.
The exact quadratic coefficient field and full Jacobian are retained.
No full-rank input or new action is claimed by this control.

## 4. Evidence and the remaining action problem

[Helper](single_square.py) · [checker](verify_single_square.py) ·
[model](single_square_model.json) · [receipt](single_square_receipt.json).

The elementary exterior identities, the transverse-curvature fixed-locus
obstruction, both ruling cases, the exact limiting-line coordinates,
all linear-factor and differentiated trace identities, and the actual
rank3 normalized control supplement the written unbounded proof.
No formal verification, independent whole-proof replay or external
review is claimed.

The remaining primitive quadratic rank4 direction families are K_B,
the mixed chain, and K_D, the product shear, with arbitrary output
matrices. Their scalars must be nonconstant and coprime to u, with
irreducible factor degrees1..3 and unbounded multiplicities.
Their exact generic groups, higher primitive directions,
nonhomogeneous fields, the absolute exceptional11 action and the full
faithful-transitive action classification remain open.
This is a complete constructor-family obstruction, not an absolute
action-level exclusion.

The previous goal turn was progress: KMFI151--152 closed the entire
arbitrary-matrix chained-square family. This turn closes the entire
arbitrary-matrix single-square family. The full goal remains active.
