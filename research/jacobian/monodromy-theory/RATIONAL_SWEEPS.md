# Nonfinite rational sweeps: full groups and polynomial affine lifts

Work in characteristic zero. Geometric groups use algebraically closed
constants. This note classifies the z-dependent rational-multiplier slice of KMFI33
in which lambda is independent of z. It does not classify every rational
multiplier or arbitrary Keller action.

**KMFI35.** Put kappa^2+kappa+1=0, lambda=w^2 and
rho=(kappa*z+R(w))/w, where R is any polynomial of degree D>=2.
The centered polynomial surface of KMFI34 and its sweep
S=X+gamma(1,w,z) have generic degree N=2D+4 and full geometric
group A_N in its natural action, although S is nonfinite.
Over any characteristic-zero coefficient field K containing kappa and
the coefficients of R, the arithmetic group is A_N or S_N according
as (-1)^D(2D+1) is or is not a square in K. The exact constant field
is K(sqrt((-1)^D(2D+1))).

Thus the entire genuine degree-at-least-two R class cannot realize
either specified PSL2(11) action. This conclusion is about these
field pairs. It does not exclude an action through a different field.
KMFI33 shows that, up to source translation/scaling, this is the entire
z-dependent rational-multiplier slice lambda=w^2.
For deg R<=1, including R=0, explicit affine equivalence gives the
regular C3 field, so KMFI1 excludes absolute Keller completion even
after independent stabilization. This settles the entire lambda=w^2 slice.

**KMFI36.** For every D>=3, an explicit choice R=alpha*w^D+w^2
below gives a polynomial Keller self-map of ordinary affine four-space
of degree 2D+4, geometric group A_(2D+4) and nonzero constant Jacobian.
Its source and target field pair is exactly S plus one independent
variable. These are actual affine maps with nonfinite input sweeps.
The alternating actions themselves were already realized in KMFI15;
this is an additional completed route, not a new abstract occurrence assertion.

**KMFI37.** More generally suppose rho depends on z, lambda is a
nonzero polynomial in w only, and the polynomial nilpotent matrix
satisfies the wave equation of KMFI33. Its single-pole normal form
has lambda=w^r,r>=2, rho=(kappa*z+R(w))/w,
kappa^2+(r-1)kappa+r(r-1)/2=0, with R Laurent polynomial and
every exponent at least 1-floor(r/2), after translations/scaling.
Remove its affine terms. If no term remains, its field is the
homogeneous KMFI34 cover of degree r+1 with natural A_(r+1)
for even r or S_(r+1) for odd r; its centered marked-origin stage
is impossible. Otherwise let j be its largest remaining exponent.
The full geometric group is natural A_(2(r+j)), even for negative j.
Over a coefficient field containing kappa and R, the arithmetic
sign class is (-1)^(r+j+1)(r+1)(r+2j-1), and its square root
gives exactly the constant field. Thus this entire rational slice,
including nonfinite sweeps and arbitrary allowed higher pole orders,
cannot realize either specified PSL2(11) action.

[Formula module](rational_sweep.py) ·
[exact checker](verify_rational_sweep.py) ·
[model](rational_sweep_model.json) ·
[receipt](rational_sweep_receipt.json).

The result determines this nonhomogeneous rational slice and completes
an infinite family of affine stages. Neither absolute PSL2(11)
occurrence nor the full action classification is settled.

## 1. Surface, eliminant and generic degree

Use zero-constant primitives I'=R, J'=wR, K0'=R^2 and put

\[
 X_1=-(\kappa+2)wz+wR+2\kappa I,\qquad
 X_2=\kappa^2w^2z+w^2R+2(\kappa-1)J,
\]
\[
 X_3=-wz^2+wRz-K_0,\qquad S=X+\gamma(1,w,z).
 \tag{1}
\]

Direct differentiation gives the polynomial normalization matrix
w^2*[[-rho,1],[-rho^2,rho]] and det JS=gamma^2.
Since X(0,z)=0, S is nonfinite. No finite-normalization inference
is drawn merely from that Jacobian.

For target coordinates Y1,Y2,Y3, set

\[
 Z=Y_2-wY_1+2\kappa wI-2(\kappa-1)J,\quad
 \Gamma=wY_1+(\kappa+2)Z-w^2R-2\kappa wI.
\]
\[
 P=w^3(Y_3+K_0)-w(Y_1-2\kappa I)Z-(\kappa+1)Z^2.
 \tag{2}
\]

On w!=0 the inverse is z=Z/w^2, gamma=Gamma/w.
Equation P=0 is exactly the remaining sweep equation.
If r_D is the leading coefficient of R, P has degree N=2D+4
and constant leading coefficient

\[
 L_D={D^2(D-1)^2r_D^2\over
             (2D+1)(D+1)^2(D+2)^2}\ne0. \tag{3}
\]

Its constant coefficient is kappa^2 Y2^2.
At fixed generic Y1,Y2, Y3 is a rational function of w with
poles of orders 3 and 2D+1 at zero and infinity. It has degree N.
Thus (2) gives the exact degree-N field extension, not a slice guess.

## 2. Birational ramification surface and distinct critical values

The surface has generic rank two: the coefficient of z in its
first-two-row derivative minor is
kappa^2(kappa+2)w^2, which is nonzero.
A normalized normal is

\[
 (\eta,\beta,1)=((\kappa-1)z+R,\;-(\kappa z+R)/w,\;1).
\]

At smooth image points eta,beta are rational functions of X.
Use w,beta as rational source coordinates. Then

\[
 \eta={R-(\kappa-1)\beta w\over\kappa},\qquad
 X_1=-(2\kappa+1)\beta w^2-2\kappa(wR-I).
 \tag{4}
\]

Over k(beta), these are polynomials in w of degrees D and D+1.
The degree of k(beta,w) over k(beta,eta,X1) divides both,
and is therefore one. Since k(eta,beta) is contained in k(X),
we obtain k(X)=k(w,z): the ramification surface is birational
onto its branch image.

Writing Y3=f(w) at fixed generic Y1,Y2, differentiation gives

\[
 f'=-\Gamma^2/w^4,\qquad
 wP'-3P=\Gamma^2. \tag{5}
\]

Gamma has degree D+2 with nonzero constant term for generic Y2.
Its Y1 and Y2 terms are respectively
-(kappa+1)Y1*w and (kappa+2)Y2, so its generic roots are simple.
They are the D+2 finite critical points, each of local degree three.
Their critical values are generically distinct by (4):
a general branch point has only one ramification-surface preimage.
The index sum is

2(D+2)+(3-1)+(2D+1-1)=2N-2.

All finite inertia generators of this rational curve are single
three-cycles. Their supports form a connected hypergraph because
the cover is transitive; the KMFI12 support lemma gives A_N.
This computes the group after algebraically closing k(Y1,Y2).
The discriminant identity below prevents any enlargement over k(Y).
This proves the geometric part of KMFI35 without assuming S finite.

## 3. The omitted sheets and the arithmetic constants

For completeness, the finite-normalization boundary can be audited
directly. P is monic up to the nonzero constant (3), and on w!=0
its hypersurface agrees with the sweep by the explicit inverse.
Possible omitted normalization sheets therefore lie over w=0,
which forces Y2=0.

The lowest homogeneous terms in w,Y2, at generic Y1!=0, are

\[
 -\kappa(Y_1w-Y_2)(Y_1w+\kappa Y_2).
 \tag{6}
\]

Put w=Y2*v and divide by Y2^2. At Y2=0 the two roots are
v=1/Y1 and v=-kappa/Y1. They are distinct and simple.
Hensel lifting makes both branches unramified over this generic
divisor. The other roots are nonzero and are controlled by the
sweep, whose branch image is not the plane Y2=0.
Thus no missing odd inertia was concealed by det JS=gamma^2.

Using deg P=N, deg Gamma=N/2 and (5), resultants give

\[
 \operatorname{Res}(P,P')
   =\left({\operatorname{Res}(P,\Gamma)\over\kappa Y_2}\right)^2,
\]
\[
 \operatorname{disc}P
   ={(-1)^D\over L_D}
     \left({\operatorname{Res}(P,\Gamma)\over\kappa Y_2}\right)^2.
 \tag{7}
\]

Indeed Res(P,w)=kappa^2 Y2^2; the equal degrees in (5) ensure
the leading-coefficient factors cancel in this resultant identity.
Its square class is exactly (-1)^D(2D+1).
Since the geometric group is A_N and N>=8, the arithmetic group
lies between A_N and S_N, with precisely that sign character.
The sign field is constant and there are no further constants.
This proves the arithmetic assertion in KMFI35.

For deg R<=1, write R=r0+r1*w and set
t0=-r0/kappa, t1=-r1/(kappa-1), z'=z-t0-t1*w.
The linear target change (Y1,Y2,Y3-t0*Y1-t1*Y2)
sends Delta to (1,w,z') and conjugates rho to kappa*z'/w.
The transformed X differs from the R=0 surface by a constant only:
if N JH=0 then both columns of JH are multiples of Delta,
and equality of mixed partials forces both multiples zero.
Evaluation at w=0 makes that constant zero.
KMFI34 gives its degree-three regular C3 action, which KMFI1
excludes as an absolute Keller field in every dimension.

For D=2, the complete centered class has
P(w;Y=0)=r_2^2*w^8/180, independently of the constant and linear
coefficients of R. Hence an origin preimage must have w=0 and
gamma=0. No member of this degree-eight centered class supplies
the unramified marked origin of a positive-weight Gao stage.
This does not contradict the actual A8 realization of KMFI15.

## 4. A marked origin for every D>=3

Let

\[
 c_D={D^2(D-1)^2\over(2D+1)(D+1)^2(D+2)^2},\qquad
 b_D={D(D-1)\over3(D+1)(D+2)(D+3)}.
\]
Choose either root alpha of

\[
 c_D\alpha^2+b_D\alpha+1/180=0,\qquad R=\alpha w^D+w^2.
 \tag{8}
\]

The discriminant is

\[
 -{D^2(D-2)^2(D-1)^2\over
 45(D+1)^2(D+2)^2(D+3)^2(2D+1)}<0.
\]

Thus alpha is nonreal and the coefficient field may be taken as
K_D=Q(sqrt(-3),sqrt(-5(2D+1))).
Put

\[
 z_0={2(\kappa+D+1)\alpha\over(D+1)(D+2)}
                        +{\kappa+3\over6},\qquad
 \gamma_0=-{D(D-1)\alpha\over(D+1)(D+2)}-1/6\ne0.
 \tag{9}
\]

At w=1, equation P(1;Y=0)=0 is exactly (8).
Equations (2) and (9) therefore give S(gamma0,1,z0)=0.
All inverse formulas are valid there and the Jacobian is gamma0^2.
For D=2, (8) degenerates to a repeated root that cancels R entirely;
it supplies no degree-eight map. No A8 occurrence is inferred
from this degenerate marked-point formula.

## 5. A complete polynomial affine stage

Put theta=(gamma,u,v), w=gamma*u,z=gamma^2*v, and

\[
 E_i(\theta)=S_i(\gamma,\gamma u,\gamma^2v)/
                    \gamma^{e_i},\qquad (e_1,e_2,e_3)=(1,2,3).
\]

These are polynomials: the X coordinates have lowest weights
3,4,5 for source weights (1,2), and gamma*Delta has weights 1,2,3.
Let theta0=(gamma0,1/gamma0,z0/gamma0^2),
J=DE(theta0), L=J^{-1}, c1 its first column, and

\[
 h=-\tfrac12(D_{c_1}^2E_3)(\theta_0),\qquad
 \theta=\theta_0+L(v_1,v_2,v_3+h v_1^2)^{\mathsf T}.
 \tag{10}
\]

The determinant chain at the marked point gives det J=1/gamma0,
so every coefficient is explicitly defined in K_D.
The stage is an affine isomorphism and a triangular shear.

On ordinary coordinates x,y,z_src,t_src put
v1=xy,v2=x^2*z_src,v3=x^3*t_src,C=x*gamma and define

\[
 F_D=(C,\ E_1(\theta)/x,\ E_2(\theta)/x^2,\
                              E_3(\theta)/x^3). \tag{11}
\]

To prove polynomiality, give v1,v2,v3 weights 1,2,3.
The linear terms of E(theta0+Lv) are v_i.
Below the three required output weights the complete lists are:
for E1 its constant; for E2 its constant and v1;
for E3 its constant, v1,v2,v1^2.
All vanish; (10) cancels the sole remaining quadratic coefficient.
Every higher Taylor term has sufficient weight.
Thus (11) is polynomial on every point of ordinary A4.

The full Jacobian chain is x^6*det L*gamma^6/C^6=gamma0.
The identity extends from the nonzero chart globally by polynomiality.
The source field inverse recovers x=C/gamma,u=w/gamma,v=z/gamma^2,
then applies J to theta-theta0, undoes the shear and divides by
x,x^2,x^3. The target inverse is Y_i=C^(e_i)*F_(i+1).
Hence the field pair is exactly S plus independent C.
KMFI35 supplies its full group and exact arithmetic constants.

For D=3, (8) is 81alpha^2+105alpha+35=0 and
gamma0=-3alpha/10-1/6. This is a compact actual degree-ten map.
For every D>=3 the arithmetic group and constants are
the formula in KMFI35 over K_D. No general specialization
or field-versal assertion is being made.

Universal geometric, field and affine-lift arguments are written here.
Exact controls are not formal verification or independent whole-proof review.
The full goal remains active.

## 6. All rational multipliers with lambda independent of z

KMFI33 gives lambda=C(w-a)^r, rho=kappa*z/(w-a)+q(w),
r>=2, after its single-pole conclusion; lambda has no lower
z coefficients because it is independent of z.
Translate a to zero and remove C by scaling.
Writing R=wq gives a finite Laurent polynomial with every
exponent at least 1-floor(r/2). Its residue equation is

\[
 \kappa^2+(r-1)\kappa+r(r-1)/2=0,\qquad
 (\kappa+r)^2=(r+1)(r/2+\kappa). \tag{12}
\]

The residue is nonreal: its rational discriminant is 1-r^2<0.
Use zero-constant polynomial primitives

I'=w^(r-2)R, J'=w^(r-1)R, K0'=w^(r-2)R^2.

All integrands have nonnegative exponents under the stated lower
bound. The complete centered surface is

\[
 X_1=-(r+\kappa)w^{r-1}z+w^{r-1}R+2\kappa I,
\]
\[
 X_2=-(r+\kappa-1)w^rz+w^rR+2(\kappa-1)J,
\]
\[
 X_3=-\tfrac r2w^{r-1}z^2+w^{r-1}Rz-K_0. \tag{13}
\]

It has the same normal (eta,beta,1) as in Section2.
All coordinates are polynomial and X(0,z)=0.
An affine term r0+r1*w in R is removed by
z'=z+r0/kappa+r1*w/(kappa-1) and the matching linear target
change of Section3. The normalization matrix is unchanged
apart from that conjugacy. The difference of centered surfaces
has N JH=0, so the uniqueness argument there makes it zero.
If R is thereby zero, KMFI34 gives the homogeneous conclusion
and its centered origin obstruction.

Suppose a non-affine term remains, and let j be the largest exponent.
It is either at least2 or at most-1.
Set F=w^(r-1)R+2kappa I, G=w^rR+2(kappa-1)J and

\[
 Z=Y_2-wY_1+wF-G,\qquad
 \Gamma=wY_1+(r+\kappa)Z-wF,\quad c=r/2+\kappa,
\]
\[
 P=w^{r+1}(Y_3+K_0)-w(Y_1-2\kappa I)Z-cZ^2.
 \tag{14}
\]

Again z=Z/w^r, gamma=Gamma/w and

\[
 wP'-(r+1)P=\Gamma^2.
\]

For a monomial R_l*w^l, its contribution to Gamma at Y=0 is
-u_l R_l*w^(r+l), where
u_l=l(l-1)/((r+l-1)(r+l)).
The l=0,1 terms vanish, as required by the affine equivalence.
It follows that P has degree N=2(r+j) and leading coefficient

\[
 L={u_j^2R_j^2\over r+2j-1}\ne0. \tag{15}
\]

The denominator is positive: the lower bound on j gives
r+2j-1>=1. Also N>r+1. For fixed generic Y1,Y2, the pole
orders in the rational w-curve are r+1 and r+2j-1.
Gamma has degree r+j, simple generic roots, and independent
parameter terms -(r+kappa-1)Y1*w+(r+kappa)Y2.
Thus there are r+j ordinary finite index-three critical points.

The ramification surface is birational onto its image.
Indeed k(X) contains beta. In coordinates beta,w, X1 and X2
are polynomials in w: the Laurent lower bound cancels every
possible pole. If j>=2 their degrees are r+j-1 and r+j;
their leading coefficients are respectively
j(r+2kappa)R_j/[kappa(r+j-1)] and
j(r-1+2kappa)R_j/[kappa(r+j)], both nonzero.
If j<=-1 their degrees are r and r+1 from the independent beta terms.
In either case the extension degree divides two consecutive
positive integers, hence is one. The rank-two hypothesis is
visible from the nonzero z coefficient
-(r+kappa)(r+kappa-1)w^(2r-2) of the first-two-row minor.

Therefore the finite critical values are distinct generically.
Their single three-cycles generate the full natural A_N.
The full Riemann--Hurwitz sum is
2(r+j)+r+(r+2j-2)=2N-2.

For the finite normalization, the w=0 quadratic is
(1-c)Y1^2*w^2+(2c-1)Y1*w*Y2-cY2^2.
Its two slopes are 1/Y1 and c/((c-1)Y1), distinct since c!=0,1.
Both missing boundary branches are unramified as before.

Finally Res(P,w)=-cY2^2, so the discriminant square class is
(-1)^(r+j+1)*c*(r+2j-1). Equation(12) replaces c by the
rational square class r+1, proving the stated arithmetic class
and exact constant field. The geometric group already fills A_N,
so no other arithmetic quotient or constant extension is possible.
This completes KMFI37. It is a whole constructor-class classification,
not an absolute PSL2(11) action obstruction.

