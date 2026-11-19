# Sharp minimum unsafe degree for every starting height

Date: 19 November 2026. Work over an arbitrary
characteristic-zero field K in the sphere algebra

    S=K[x,y,t]/(t²+2xy-1),
    p=t^a(C(x)+t²D(x))/x, a>=0 an integer, C,D inK[x].

Regularity is C(0)+D(0)=0. Purity means every positive angular moment
vanishes. Safety means eventual vanishing for every fixed sphere marker,
and for every fixed marker in the original two-pair contraction model.
The least sphere degree is the radial-primitive balanced lift degree.

**Theorem AllStartSharpMinimumUnsafeDegree.** Every pure phase of least
degree strictly less than2a+4 in this class is safe. In the unit-center
branch C(0)!=0, after normalizing C(0)=-D(0)=1, every such phase has
C=1,D=-1 and p=2t^a y. The minimum unsafe degree within this entire
class is exactly2a+4 for each integer a>=0. The upper bound is the already
retained F_k witness with k=a; no new witness is claimed here.

More strongly, every nontrivial pure unit-center phase has a geometric
root alpha of the numerator P of the rational selected square root
X=P/J, in reduced form, with

    v_alpha(H)>=a*v_alpha(P)+1,
    v_alpha(D)>=(a+2)*v_alpha(P)+1>=a+3,                (M1)

where C=HJ²,D=-HP² and gcd(P,J)=1. Thus the lower bound follows from
one LOCAL multiplicity, not only from a total-degree inequality.

## 1. Rational coordinate and exact/logarithmic separation

For odd a>=1 the retained OddStartRationalCenter theorem gives

    X=sqrt(-D/C) inK(x), X(0)=1,
    u=x X^a/C inK(x), u=x+O(x²),
    omega=W(x)dx, W=X^(a-1)/C inK(x).                 (M2)

Let v(u) be the inverse formal branch. The all-order endpoint criterion is

    Phi(u)=W(v(u))*v'(u), Phi(u)+Phi(-u)=2.             (M3)

For completeness, the same endpoint argument proves rationality when a is
even, including0, without any additional external input. The initial
selected X is its formal square root, and u=xX^a/C is already rational
because a is even. All positive m now have even am, so all endpoint
factors kappa_(a,m) are nonzero. With A=C/X^a,L=1/X, formal residue
substitution gives

    Lambda(p^m)=kappa_(a,m)[x^m]A^mL
                =kappa_(a,m)[u^m]Phi(u).

Endpoint order m+1 contact justifies replacing the formal integration
endpoints plus/minusX by plus/minus1, exactly as in the odd-start proof.
Purity therefore forces Phi=1, hence W=u'. Dividing u=xX^a/C by
x*u'=xX^(a-1)/C gives X=u/(x*u') inK(x), including a0. Thus(M2) and
(M3) hold for ALL integer a>=0; even a has the stronger Phi=1.
No finite moment-prefix test is used. Over k=Kbar, decompose the rational
differential into its rational primitive and logarithmic parts:

    omega=dR(x)+sum_i c_i dlog(x-r_i), R(0)=0.          (M4)

The polynomial part of W is included in R. All r_i are nonzero since
W is regular at0. The rational part R is unique after normalization,
and it descends to K: partial-fraction integration, or uniqueness under
coefficient-field automorphisms, gives this directly.

Equation(M3) says that the normalized formal primitive V of omega
satisfies V(v(u))-V(v(-u))=2u. In the common finite function field
k(u,v(u),v(-u)), the retained FunctionFieldExactLogSeparation lemma
therefore gives

    R(v(u))-R(v(-u))=2u.                              (M5)

This application needs no special shape of W. The exact/log lemma is
proved by trace/norm and the impossibility of a nonzero constant being
a sum of rational logarithmic derivatives at infinity.

Embed K(x) into K((u)) by x=v(u). The finite intermediate field
E=K(u,R(x)) is stable under the formal involution u->-u, because(M5)
sends R toR-2u. Let sigma be this genuine involution of E. It need NOT
extend to the original K(x). Since u is nonconstant, R is nonconstant.
The field tower and sigma give equality of rational MAP degrees:

    deg_x(R-2u)=[K(x):E]*[E:K(R-2u)]
               =[K(x):E]*[E:K(R)]=deg_x R.           (M6)

No Luroth theorem, rational return on the original line, or polynomial
quotient is used in this degree comparison.

**Lemma PoleDominanceObstruction.** Suppose u,R are nonconstant rational
functions satisfying(M6). There must be at least one pole of u where
the pole order of R is at least the pole order of u.

Indeed, otherwise at every u-pole R-2u has the full larger u-pole order;
at every other point it has exactly the R-pole order. Summing the pole
orders on the geometric projective line gives

    deg(R-2u)=degR+sum_(u-poles)(pole_u-pole_R)>degR,

contrary to(M6). A regular R has pole order0 in this formula. This
argument also shows why the field involution is used only for its exact
degree equality, without an unsupported involution extension.

## 2. All finite and infinite pole orders

Write the selected rational square root in reduced form

    X=P/J, P,J inK[x], gcd(P,J)=1, P(0)=J(0)=1.

Polynomiality of C,D and -D/C=P²/J² forces J²|C; hence

    C=HJ², D=-HP², H inK[x], H(0)=1.

Here H is generally NOT squarefree. Reducing P/J is essential for the
following local comparison. Equation(M2) becomes

    u=xP^a/(H J^(a+2)),
    omega=P^(a-1)/(H J^(a+1)) dx.                    (M7)

Let p=degP,j=degJ,h=degH. All finite roots under consideration are
nonzero. At a geometric point write p_alpha,j_alpha,h_alpha for the
respective multiplicities; at most one of p_alpha,j_alpha is positive.
If a rational differential has a pole of orderq>=2, its rational
primitive part has pole orderq-1. A simple differential pole contributes
only a logarithm, hence no pole of R. A regular differential also has
regular R there. These assertions are exact in characteristic zero.

At a J-root the u-pole and R-pole orders are respectively

    h_alpha+(a+2)j_alpha,
    h_alpha+(a+1)j_alpha-1.                           (M8)

The latter is strictly lower. At a root of H which is not a root of
P or J, their orders are h_alpha and max(h_alpha-1,0), again lower.
At a P-root there is a u-pole ONLY if h_alpha>a*p_alpha. Consequently,
unless(M1) holds at some P-root, R has strictly lower pole order at
EVERY finite u-pole.

It remains to check infinity, including polynomial parts of R. The
growth exponent of u there is

    q=1+a*p-h-(a+2)j.                                 (M9)

If q>0, u has a pole of orderq. The R-pole order is

    max((a-1)p-h-(a+1)j+1,0)=max(q+j-p,0).            (M10)

If j<p this is strictly less thanq. If j>=p, then

    q<=1-h-2p.

For q>0 this forces p=0,h=0; the original formula(M9) then forces
j=0 as well. This is exactly P=J=H=1, the exceptional normalized phase
C=1,D=-1. Thus for every NONTRIVIAL unit-center phase, R also has
strictly lower pole order at infinity whenever u has a pole there.
Finite or zero u(infinity) contributes no u-pole and needs no separate
assumption.

If(M1) failed at every P-root, all u-poles would have strictly smaller
R-pole order, contradicting PoleDominanceObstruction. Hence(M1) holds.
Since D=-HP², its local multiplicity there is at least(a+2)p_alpha+1,
so degD>=a+3. For ALL integer a>=0, write a=2q+epsilon with
epsilon0 or1. The unique normal form of t-degree at most1 is

    t^epsilon(1-2xy)^q*((C+D)/x-2yD).

Its two possible highest contributions have different highest y-powers,
so they cannot cancel when their total degrees agree. Reduction using
t²=1-2xy never increases total degree; uniqueness of this normal form
therefore proves the exact least-degree formula for even a as well as
odd a. It gives

    deg_sphere p=max(degC+a-1,degD+a+1)>=2a+4.         (M11)

This proves the unit-center lower bound over k. All conclusions about
original degrees hold over K without adjoining the root alpha.

## 3. Complete equality classification

**Theorem AllStartMinimumDegreeClassification.** Among nontrivial normalized
pure unit-center phases, the following are equivalent:

1. degD=a+3;
2. the least phase degree is2a+4;
3. for a unique lambda inK*, putting U=1+lambda*x, one has
   C=U^(a+1), D=-U^(a+3).

Condition2 implies1 by(M11). If1 holds, the compensating P-root in(M1)
accounts for ALL a+3 roots of D, counted with multiplicity. Necessarily
p_alpha=1,h_alpha=a+1, and there are no other roots of P or H. Thus
P=U,H=U^(a+1) in the reduced factorization. Both are defined over K,
so lambda belongs toK*. It remains to prove J=1.

For this remaining possible factor J, which is coprime toU,

    u=x/[U J^(a+2)], omega=dx/[U² J^(a+1)].            (M12)

At every J-root of multiplicityb>=1, the respective u-pole and R-pole
orders are(a+2)b and(a+1)b-1. Consequently the pole degree of R-2u
GAINS b+1>=2 there. At U=0 both R and u have a simple pole, so their
subtraction can LOSE at most1. There are no other poles: both u andR
are regular at infinity, and the differential has no other poles.
Hence any nonconstant J gives

    deg(R-2u)-degR>=sum_(J-roots)(b+1)-1>0,

contradicting(M6). This proves J=1. The argument in fact shows that
minimal degD already forces the unique family, with no additional
bound on degC or degJ.

Conversely the displayed family is exactly the retained sharp F_k
family with its allowed nonzero coordinate scale. It has
u=x/U and Phi=1. Equation(M3) proves all-order purity; its exact degree
is2a+4. This proves all implications, rather than only a necessary
parameter condition. For general lambda, the phase is lambda times the
lambda1 phase after the sphere automorphism x->lambda*x,y->y/lambda.
The original automorphism xi1->lambda*xi1,z1->z1/lambda fixes R and
sends Z->lambda*Z,W->W/lambda, preserving primitivity and nonpower
properties. The x-marked value gains precisely lambda^(m-1).

## 4. Safety below the bound and the retained sharp witnesses

For odd a the retained all-degree zero-center theorem already makes every
pure phase with C(0)=D(0)=0 safe. For even a this follows directly from
the first TWO pure moments. Write C=x(c+O(x)),D=x(d+O(x)). The torus
weight-zero part is t^a(c+d*t²), while every remaining term has strictly
positive weight. Consequently purity gives

    c/(a+1)+d/(a+3)=0,
    c²/(2a+1)+2cd/(2a+3)+d²/(2a+5)=0.

Substitute c=-(a+1)d/(a+3). The second coefficient of d² is

    4(2a²+6a+3)/[(a+3)²(2a+1)(2a+3)(2a+5)] !=0.

This nonzero rational forces c=d=0 in every characteristic-zero field.
Thus C,D are divisible byx² and the phase has strictly positive torus
weight, which gives safety for sphere markers and original E2 markers
by equivariance and the nonpositive weights of E2 outputs. This proves
the zero-center input for even a, including0, within this proof. Below2a+4 the unit-center proof leaves only
p=2t^a y. Its sphere safety follows from torus weight: x has weight1,
y weight-1 and t weight0, and any fixed marker has finitely many
weights, while p^m has weight-m.

For ORIGINAL markers, its balanced lift is2T^a W with W=2xi2*z1.
It has weight-1. Split a fixed original marker into bidegrees(r,s) and
weights. After multiplying by this balanced phase to the mth power,
the E2 output, if nonzero, has fixed z-degree s-r, independent ofm;
it is zero when r>s. Its possible output weights are the bounded
interval[-(s-r),0]. The input weight is the fixed marker weight minusm,
so equivariance forces eventual zero. This is the repaired bounded
output-weight argument used in the retained odd-start proof; E2 does
not kill every nonzero torus weight.

For sharpness reuse the existing F_k family in
TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md Theorem5.5/(5.11), with k=a. Its sphere phase is

    C=(1+x)^(a+1), D=-(1+x)^(a+3),
    u=x/(1+x), Phi=1, degree2a+4.

Every pure moment vanishes; the fixed marker x has nonzero moments for
every positive exponentm with am EVEN, since

    Lambda(x p^m)=kappa_(a,m),
    kappa_(a,m)=integral_0^1 t^(am)(1-t²)^m dt !=0.

The original marker Z likewise survives with the retained factorial
factor((2a+4)m+2)!. That existing proof also establishes radial
primitivity and failure of proper powers. Thus these are already known
unsafe witnesses at precisely the lower bound. This theorem supplies
the new universal lower bound and its local multiplicity mechanism.

In particular the exact minimum is4,6,8,10,12,... for a=0,1,2,3,4,...;
all a3 phases through degree9 are safe. There is no statement here
about arbitrary height profiles or arbitrary degree-d Image phases.

## 5. Classical scope and verification boundary

The odd rational-center and endpoint descent are imported from the retained
all-odd-start proof; the even rationality and zero-center steps are proved
above, including a0. The exact/log separation is the retained
self-contained Liouville trace/norm lemma, whose primary Rosenlicht
antecedent is already recorded there. Pakovich's rational-moment
inverse-branch work supplies a related method antecedent with different
hypotheses; no polynomial-moment composition classification or external
rational-map theorem is imported for the pole comparison.

The focused checker reconstructs exact partial-fraction primitives, all
local/infinite order formulas, the universal even zero-core determinant
numerator, map-degree obstruction controls and the already retained
witness moments using exact rational arithmetic. Such bounded controls are
separate from the written all-order proof(M3)--(M11).

The [source ledger](all_start_height_sharp_unsafe_degree_source_ledger.json) pins all three retained inputs and distinguishes the existing sharp witnesses from the new lower bound. The [checker](verify_all_start_height_sharp_unsafe_degree.py) and [receipt](all_start_height_sharp_unsafe_degree_verification.json) bind exact finite divisor and original-contraction controls to this universal proof.

```bash
python3 research/foundations/extended-geometry/verify_all_start_height_sharp_unsafe_degree.py \
  --output research/foundations/extended-geometry/all_start_height_sharp_unsafe_degree_verification.json
```
