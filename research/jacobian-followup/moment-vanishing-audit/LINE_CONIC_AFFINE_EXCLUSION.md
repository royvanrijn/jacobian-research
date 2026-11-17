# Every line–smooth-conic leading cubic is GMC-safe

Current quantitative status,15 November2026: the [global pure-prefix proof](GLOBAL_TERNARY_CUBIC_PURE_PREFIX.md)
now gives a sufficient1855 pure-only gate and least-prefix interval[10,1855].
Earlier absence-of-global-upper or unquantified-node statements below are
historical. Their earlier class/certificate scopes remain unchanged; the
new proof remains internal and independently unreviewed.

The subsequent [weighted binary finite-prefix proof](SQUAREFREE_SINGULAR_FINITE_PREFIXES.md)
supplies the sufficient pure gate tangent935 (the transverse class remains unquantified). The qualitative
arguments below retain their scopes; earlier absence-of-prefix statements
are historical. No numerical arbitrary-cubic bound follows yet.

Current status,15 November2026: the [full ternary cubic theorem](TERNARY_CUBIC_GMC_CLOSURE.md)
closes the degree<=3 GMC decision as an internal written proof. Any
open cubic frontiers below are historical; the earlier results retain
their stated scopes. Independent review and numerical global pure-only
gates remain separate.


14 November 2026. Internal written formal/analytic proof with exact
auxiliary identities. No independent replay, external review or formal
verification. The checker does not certify the completeness arguments.

## Results

**G3LINECONIC20261114.** Every ternary cubic whose leading homogeneous
form is a line times a smooth conic satisfies GMC, for arbitrary lower
terms. All transverse configurations are pure-nonvanishing. Tangent
configurations are pure-nonvanishing if the tangency point is nonisotropic
or F2 is nonzero there; the remaining isotropic zero-F2 case uses the
retained W-linear GMC theorem.

**G3REDUCEDLEADING20261114.** Every ternary cubic with squarefree leading
form is GMC-safe. Moreover any pure-null member has a nonzero isotropic
leading singular direction v with F2(v)=0. Thus it admits a Gaussian
W-linear frame and satisfies the sufficient mixed cutoff m>3*deg Q.
Sharpness of3 within the squarefree-leading locus is not asserted.

The only remaining possible ternary cubic GMC failures have nonreduced
leading form L^2*M or L^3. No explicit global upper pure prefix follows;
coefficient schemes and lower-power GVC remain open.

Throughout, Q is the full affine Gaussian exponent, with invertible
complex symmetric Hessian. All means and cross coefficients survive
coordinate changes. We use the
[positive formal seed criterion](CUSPIDAL_AFFINE_EXCLUSION.md), the
[confluent numerator frame](AFFINE_TRIANGLE_EXCLUSION.md), and the
[uniform cut/current-limit argument](NODAL_CENTER_LINE_DECAY.md) only
after checking their respective completeness or coercivity hypotheses.

## A. A transverse line and conic

The two intersections are ordinary nodes. An isotropic zero-F2 node
gives a Gaussian W-linear frame with nondegenerate quadratic W coefficient
and is pure-excluded. If both nodes have F2!=0, the retained ordinary-node
coercive criterion applies. Otherwise choose a nonisotropic zero-F2 node,
put it on the T-axis, and normalize by linear changes and translations:

    P=T*(uv+h)+a*u^3+b*u^2+c*v^2+eu+fv+k, a!=0,
    Q=(A*u^2+B*v^2+m*T^2+2*Cuv+2*DuT+2*EvT)/2
        +l_u*u+l_v*v+l_T*T+q0,
    m!=0; if c=0 then B!=0.                       (1)

The other node is the v-axis. The leading form is u*(Tv+a*u^2),
whose quadratic cofactor is smooth precisely because a!=0. The
quadratic values of S_t=tP-Q at its two nodes are -m/2 and ct-B/2,
nonzero for large t. Each infinity point has length1, so the affine
critical and top twisted rank is8-2=6 by the retained tame comparison.

### A1. Finite critical sets

For h!=0, restrict to uv=-h:

    R(u)=a*u^3+b*u^2+eu+k-fh/u+c*h^2/u^2.

Its nonzero exponent support contains3 and possibly2,1,-1,-2.
The finite Milnor rank is3 minus the minimum support exponent, at most5.
Every point has exact local type alpha+s*w+y^(mu+1), by smooth conic
coordinates and a T shear. The T-node has one quadratic-action escape.
If c=0 the v-node supplies1 when f!=0,3 when f=0,e!=0,4 when f=e=0,b!=0,
and5 when f=e=b=0. These are the explicit nodal1/3/4/5 profiles, with
nonzero action and bounded physical Gaussian fluctuations. Their sum
with the finite rank is always6.

For h=0, the critical set is finite exactly when (c,f) is a nonzero
pair and at least one of e,f is nonzero. Its rank is

    F=(2 if e!=0, else1 if b!=0, else0)+indicator(c*f!=0).

The T-node supplies3 escapes when e*f!=0;4 when e!=0,f=0; and when
e=0,f!=0 it supplies4 if b!=0 and5 if b=0. The exposed v-node supplies
indicator(c=0). Again the total is6. Multiple finite roots retain
their full Milnor multiplicities and all finite zero-exponential seed
powers t^(-1-j/(mu+1)) are strictly positive. Every escaping factor
has a nonzero positive action. A generic linear form and the degree<6
confluent numerator basis give a complete formal frame: escape physical
fluctuations are bounded and all required separations grow. The seed
criterion excludes pure zeros, including equal-action resonances.

### A2. A critical hyperbola

The non-finite case with e!=0 is h=c=f=0. Then

    P=u*(Tv+a*u^2+b*u+e)+k.

Its critical hyperbola u=0,Tv=-e has value k. There are also the
nonzero roots of3a*u^2+2bu+e, of total Milnor rank2, on T=v=0.
Set s=Tv+e and scale u=rU,s=rS,r=t^(-1/2). Exactly,

    t(P-k)=US+bU^2+a*r*U^3,
    dT dv du=r^2*(dT/T) dS dU.

The leading Gaussian normal pair has determinant -1. The restricted
Gaussian on Tv=-e is a Laurent phase with nonzero quadratic poles
at T=0 and infinity, since B*m!=0. Its logarithmic twisted quotient
has rank4, by the four-step Euler recurrence proved in the triangle
note. Coefficientwise fixed-phase perturbation gives four formal
characters with seed exp(k*t)*t^-1*C[[r]]. Together with the finite
rank2 they exhaust6. Polynomial high-order cutoffs in u separate the
hyperbola u=0 from the nonzero finite roots;1,T,T^2,T^3 give its four
leading Laurent amplitudes, and finite Hermite jets give the other two.
The block leading matrix is invertible and every zero-action seed
has positive valuation.

### A3. Critical axes

The remaining case is h=e=f=0. P-critical values are only
k and, if b!=0, k+4b^3/(27a^2). There is a common complex linear
gradient bound for S_t outside one fixed ball, uniform for large |t|.
Near the T-node, the v equation first bounds u; the u equation then
bounds v. The resulting estimates
|u|<=C/|t|+C|v|/|T| and |v|<=C/|t|+C'|u|/|T|
absorb for large |T|, giving u,v=O(1/t). The T equation has leading
-mT, a contradiction. Near the v-node the T equation first gives
u=O(1/t), then the u equation T=O(1/t); the v equation has leading
(2ct-B)v, uniformly nonzero for large |t|. This proves the uniform
bound without discarding the cubic a*u^3.

Critical S_t branches are bounded, with actions alpha*t+O(1), alpha
in that at-most-two-element set. The rank6 carrier, tail and current-limit
argument therefore gives decay of every normalized scalar cut period.
Any two complex critical values fit in a closed half-plane through0;
choose its supporting t-ray. All six actual scalar periods decay,
contradicting the nonzero constant period supplied by all pure zeros.
No semisimple formal assertion is imposed on intersecting axes.

These alternatives exhaust32 support patterns:26 finite,2 hyperbola,
4 axes. Thus every transverse affine member is pure-nonvanishing.

## B. A tangent line and conic

Normalize F3=T*u^2+u*v^2. Its only leading singular direction is the
T-axis. If the T^2 coefficient alpha in P is nonzero, P is complex
linearly gradient-coercive. Indeed near that direction u,v=o(T);
small normalized v gradient gives uv=O(T), while small u gradient
gives v^2+2Tu=O(T). If u is unbounded, these imply u=O(|T|^(1/3));
if u is bounded the same conclusion u^2/T->0 is immediate. The T
gradient divided by T then tends to2alpha!=0, a contradiction.

The isotropic zero-F2 tangency is GMC-safe by the exact W-linear
criterion. In the remaining nonisotropic case retain m=Q_TT!=0.
Translations of u,T,v remove uT,u^2,uv respectively and give

    P=T*(u^2+g*v+h)+(u+c)*v^2+eu+fv+k.           (2)

The full Q is as in (1). For S_t, the T-node has nonzero quadratic
value -m/2. The preceding gradient argument proves coercivity for
each fixed nonzero t. Its projective local length is3: the T equation
solves the homogenizing w=O(u^2); the u equation gives u=-v^2/2+O(v^4);
the v equation then has leading -t*v^3. Hence the full twisted rank
is8-3=5. We now account for all five dimensions in every case.

### B1. g!=0

On u^2+gv+h=0 the restricted P has degree5 in u, so finite Milnor
rank4. Its local types are the one-variable powers plus a quadratic
pair, including all degeneracies. There is one escape

    u=tU,v=t^2 V,T=t^3 W,
    U=-5g^2/(8m), V=25g^3/(16m^2), W=125g^4/(64m^3).

Its limiting phase is W(U^2+gV)+UV^2-mW^2/2, with action
3125g^8*t^6/(8192m^5) and Hessian determinant625g^6/(128m^3).
Gaussian fluctuation scales are t^-2,t^-1,1 and prefactor t^-3.
All corrections have positive inverse powers. The numerator frame
combines finite4 and this escape into rank5, proving pure exclusion.

### B2. g=0,h!=0

Put Delta=c^2+h. The u=+sqrt(-h) and u=-sqrt(-h) critical branches
restrict to (u+c)*v^2+fv+eu+k. A nonzero u+c gives one finite Morse
point. A zero u+c gives no point if f!=0, or a critical parabola if f=0.

If c!=0,Delta!=0 there are two finite Morse dimensions and three
Morse escapes. One primary escape has T~(h/m)*t,u=O(1/t),v bounded,
and action h^2*t^2/(2m). The secondary pair has
u=-c+o(1),T~Delta*t/m,v~t^(1/2)*V,V^2=2c*Delta/m,
action Delta^2*t^2/(2m) and limiting Gaussian determinant8c*Delta.
Its fluctuation scales are t^-1,t^(-1/2),1; the primary has the same
scales and a nondegenerate limiting Gaussian. Generic projections
separate the pair at least at t^(1/2) and the primary at t. The complete
numerator frame has rank2+1+2=5, with only positive finite zero-action seeds.

If c=0, define L=f-Eh/m. For L!=0 there are two finite Morse dimensions
and three escapes near T=(h/m)*t. Write r=t^(-1/3),
u=rU,v=V/r,T=h/(m*r^3)+W/r. The leading critical system gives

    V^3=Lh/m, U=-mV^2/(2h), W=(U^2-EV)/m.

Its Jacobian is6mV^2!=0. After the common nonzero quadratic action,
the next action is (3LV/4)*t^(4/3). Fluctuations t^-1,t^(-1/3),1
give a nondegenerate limiting Gaussian; all separations needed for
the rank5 numerator frame grow, even if the W coefficients coincide.
If L=0, instead use u=U/t,T=(h/m)*t+W,v=V. After the same nonzero
quadratic action and a scalar linear action, the fixed phase is

    (h/m)U^2+U*(V^2+e-Dh/m)-mW^2/2-EVW
       -l_T*W-BV^2/2-l_v*V-q0.

Completing U,W gives a quartic V coefficient -m/(4h)!=0. This fixed
rank3 block lifts formally, with prefactor t^-1. Finite2 plus block3
is a full frame, by coarse T cutoffs. All three block actions retain
the nonzero quadratic leading term. Thus both L alternatives are closed.

If c!=0,Delta=0, one branch has one finite Morse point and the other
is near u=-c. There is still the primary quadratic-action escape.
For f!=0 the remaining three escapes have r=t^(-1/3),
u=-c+rX,v=V/r,T=W/r^2, with

    V^3=2c^2*f/m, X=-mV^2/(4c^2), W=V^2/(2c).

Their critical Jacobian is6mV^2; the action is (3fV/4)*t^(4/3), nonzero.
For Gaussian validity use the exact normal coordinates
x=u+c,s=(x-2c)T+v^2+e, in which P=k-ec+x*s+fv.
After x=rX,s=Y/r,v=V/r the leading phase is a quartic in V at order
r^-4, followed by a nondegenerate XY Gaussian at order r^-3.
Their fluctuation widths in these scaled variables are r^2 and
r^(3/2), with vanishing mixed corrections. Physical T fluctuations
remain bounded; treating them as independent diagonal u,v,T scales
would miss the normal-coordinate correlation. The volume is t^(-4/3).
The four escapes and finite1 form a full numerator frame.

For f=0 the same exact normal chart has P=k-ec+x*s. Scaling x,s
by t^(-1/2) gives a fixed quartic restricted Gaussian in v, of rank3,
with prefactor t^-1. Along the parabola T=(v^2+e)/(2c), its quartic
coefficient is -m/(8c^2)!=0. The finite point has u=c, the primary
escape u->0, and the parabola u=-c. Polynomial cutoffs at these three
distinct u values separate the fixed rank3 from the two Morse dimensions.
Thus1+1+3=5, and all zero-exponential seed powers are positive.

### B3. g=h=0,c!=0

Let v0=-f/(2c),delta=e+v0^2,alpha=k-f^2/(4c). The secondary pair
from B2 always exists, with Delta=c^2 and nonzero quadratic action.
If delta!=0 there are three additional escapes near v=v0:

    u=t^(-1/3)*U,T=t^(1/3)*W,v=v0+O(t^(-1/3)),
    U^3=-delta*m/2,W=U^2/m.

The v direction has a nondegenerate quadratic with coefficient c;
the U,W leading phase is WU^2+delta*U-mW^2/2. These are three Morse
characters with nonzero (3delta*U/4)*t^(2/3) action after alpha*t.
Their physical fluctuations are bounded. The numerator frame gives3+2=5.

If delta=0 the T-line is critical. With u=rU,v=v0+rV,r=t^(-1/2),
the leading phase after alpha*t is

    T U^2+2v0UV+cV^2-mT^2/2-(Ev0+l_T)T-constant.

Completing V,T gives a nonzero quartic U term1/(2m), hence a fixed
rank3 block with scalar seed t^-1. Its three polynomial amplitudes
1,t^(1/2)u,tu^2, times an eighth-power cutoff at u=-c, are separated
from the secondary pair by the reverse cutoff at u=0. Their worst
growth at the pair is t, while the cutoff is O(t^-4). The leading
matrix is invertible. Positive seed valuation closes the whole boundary.

### B4. g=h=c=0

If f!=0 all five dimensions escape. Set r=t^(-1/5),u=rU,v=V/r,T=W/r^3.
The leading phase r^-6*(WU^2+UV^2+fV-mW^2/2) has five Morse points

    U^5=-mf^2/8,V=-f/(2U),W=U^2/m.

Their actions are (5U^4/(2m))*t^(6/5), and Hessian determinants
-20U^3!=0. Gaussian scaled widths r^3 give physical fluctuations
r^4,r^2,1 and prefactor t^(-6/5). Their five W values are distinct,
so interpolation gives the full frame and no zero-exponential factor.

If f=0,e!=0 there are two critical T-lines at v_i^2=-e. Each gives
a fixed rank1 Gaussian after u=t^(-1/2)U,v=v_i+t^(-1/2)V: the apparent
cubic TU^2 is removed by V->V+TU/(2v_i), leaving a Gaussian pair and
the nondegenerate T Gaussian. Its scalar seed is t^-1. Three more
escapes use r=t^(-1/3),u=rU,v=rV,T=W/r, U^3=-em/2,W=U^2/m.
Their U,W block has nonzero action (3eU/4)*t^(2/3), while V has the
fixed nondegenerate quadratic U_i V^2-EW_i V. They give rank1 each,
prefactor t^-1. Cutoffs in v separate the two nonzero line centers
from the escaping center v->0; interpolation in the three U roots
gives the remaining block. Thus2+3=5 with positive zero-action seeds.

Finally f=e=0. Put u=r^2 U,v=rV,r=t^(-1/4), leaving T unchanged.
After k*t the fixed phase is

    F0=T U^2+UV^2-mT^2/2-l_T*T-q0.

Completing T gives U^4/(2m)+UV^2-l_T U^2/m plus one Gaussian.
The positive-weight leading part U^4/(2m)+UV^2 has isolated gradient
and Jacobian rank5, with basis1,U,U^2,V,V^2 (weights1/4,3/8).
Weighted leading-term reduction gives the same rank for the fixed
twisted quotient, including every value of l_T. All remaining terms
of S_t are positive powers of r. Fixed-phase perturbation supplies
the full five functionals with seed exp(k*t)*r^3*C[[r]]. For k=0
the valuation is positive; for k!=0 the exponential is nonzero.
This finishes all tangent cases.

## C. The squarefree-leading consequence

A squarefree projective plane cubic is either irreducible (smooth,
nodal or cuspidal), a line with a smooth conic (transverse or tangent),
or three distinct lines (independent or concurrent). The retained
theorems and Parts A/B therefore cover every case. All smooth, nodal,
triangle and transverse line-conic members are pure-excluded. For the
other three classes, the new proofs pure-exclude nonisotropic singular
directions and nonzero F2 values. Any surviving pure-null member must
therefore have an isotropic singular direction with F2=0. The exact
W-linear criterion gives the frame and its sufficient m>3*deg Q mixed
bound. This is a written case-cover corollary, not a new finite-prefix
certificate or a claim of a sharp mixed bound inside the reduced locus.

The [checker](verify_affine_line_conic_profiles.py),
[receipt](affine_line_conic_profiles_verification.json) and
[source audit](affine_line_conic_source_audit.json) distinguish the exact
model algebra from formal completeness and analytic cut arguments.
Only nonhomogeneous repeated-line leading forms L^2*M or L^3 remain
possible ternary cubic GMC failures. External review is still required.
