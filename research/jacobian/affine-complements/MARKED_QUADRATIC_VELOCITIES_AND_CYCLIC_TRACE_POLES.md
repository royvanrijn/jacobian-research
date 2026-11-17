# Marked quadratic velocities and cyclic trace poles

PF2RULED143. 16 November 2026. ALL complex polynomial plane Keller degrees.
The full HC4/JC2 goal remains OPEN. The new reduction uses the ACTUAL
constant Jacobian and FULL finite normalization. It does not assume that
the marked quadratic traces below are regular.

## 1. A finite marked quadratic test for a generic first velocity

Let B=C[P,Q] subset R=C[x,y], Jac(P,Q)=1, F=Frac(B), L=Frac(R).
Let S be the FULL finite normalization of B in L and
Lambda={z in L:Tr(z*S) subset B}. Retain EVERY omitted prime D_i,
ramified or unramified, with target image C_i and index e_i.

Choose a constant source direction delta generically as in PF141.
Write q_i=-min(ord_Di x,ord_Di y)>0. Then choose a constant linear
target function f=alpha*P+beta*Q generically and put a=delta(f).
PF141's two-entry velocity calculation gives

    ord_Di(a)=-(q_i+e_i) for EVERY omitted D_i.       (1)

For each prime the lowest coefficients of delta(P),delta(Q) can cancel
for at most one constant target direction; hence only finitely many
directions need be excluded. The same a is polynomial on the WHOLE
source. No representative prime per image replaces the full list.

For each distinct image curve C, let H_C be its reduced target equation.
Define

    m_C=1+max_(i:C_i=C) ceil((2*q_i+1)/e_i).          (2)

The exact target ideal of marked quadratic denominators is

    {h in B:Tr(h*s*a^2) in B for EVERY s in S}
       = product_C H_C^m_C * B.                    (3)

Equivalently, for ANY finite B-module generating set s_1,...,s_N of S,

    max_k pole_C(Tr(s_k*a^2))=m_C.                 (4)

The maximum is zero off the nonproper image divisor. In particular m_C>=2;
an omitted unramified prime has q_i>=2 by PF141, and forces m_C>=6.

Here is the proof with the full codifferent threshold. By the DEFINITION
of Lambda, the condition on h is h*a^2 in Lambda. The retained canonical
module identity gives at EVERY height-one prime of S

    Lambda_D = u^(1-e_D)*S_D, up to a unit.          (5)

This includes all tame ramified primes and their residue extensions;
e_D-1 is not replaced by zero. At an omitted D_i, (1) and (5) give

    e_i*ord_C(h)-2*(q_i+e_i) >= 1-e_i,
    ord_C(h) >= 1+ceil((2*q_i+1)/e_i).              (6)

At primes meeting the source, a is regular and causes no condition.
Lambda is reflexive over the normal surface S, so these height-one tests
are globally sufficient, including singular finite points. Target B is a
UFD; taking the maximum over ALL primes above each C proves (3).
Testing module generators is equivalent to testing all S, proving (4).
The elementary trace-duality facts are LY0/LO1 inputs; the new step is
their application to PF141's exact first-velocity pole charge.

Thus, for these generic choices, the following are equivalent:

    every Tr(s_k*a^2) belongs to B;
    a belongs to S;
    the original plane Keller map is invertible.   (7)

Indeed a in S makes its square lie in S subset Lambda. Conversely the
quadratic condition and (6) leave no omitted prime. Normal Hartogs then
makes R=S; finiteness and plane etaleness force the inverse. More directly,
any omitted prime would contradict (1) and (5). A proper affine normal
enlargement cannot have only a codimension-two complement, by LO0.

This is a degree-TWO marked test, not a degree-independent number of
unmarked moments. The markers must generate the ACTUAL full S over B.
They cannot be replaced by 1, a chosen smaller order, a determinant, a
relative source-image ring or an arbitrarily rescaled lattice. Regularity
in (7) is UNSUPPLIED; (7) is not a proof of HC4/JC2.

## 2. Apply the actual coframe to MX0's surviving cyclic family

Retain the ENTIRE MX0 family, all g,k>=2, square-free polynomial b(u) of positive degree
with b(0)!=0. Put

    u=t^g, a0(t)=t^2*b(t^g)^k,
    sigma(t)=-b(u)+g*(1-k)*u*b'(u),
    W: a0*z=v^g-sigma^g.

The free action (t,v,z)->(epsilon*t,epsilon*v,epsilon^-2*z) gives
V=W/mu_g and the actual SURJECTIVE ETALE plane source with chart

    v=sigma+a0*x, z=(v^g-sigma^g)/a0.

MX0 supplies smooth affineness, exact unit volume, all whole charts and
C=O(V)=R intersect E, E=C(V), R=C[x,t]. No plane target is known.

Now impose a HYPOTHETICAL actual target pair P,Q in C with
Jac_(x,t)(P,Q)=1. This is precisely the missing target assumption;
it is used below, not inferred from the geometry of V.

For ANY F in C, regard F rationally as a function of t,v and define

    X_F=a0*F_v=F_x,
    A_F=F_t|v+(a0'/a0)*v*F_v,
    B_F=(sigma'-(a0'/a0)*sigma)*F_v.

Then F_t=A_F+B_F. Under the cyclic action X_F,A_F,B_F have respective
characters +1,-1,-2. In particular B_F=K(t)*X_F for the SAME K for all F.
At the retained source line t=0,

    K(t)=2*b(0)^(1-k)*t^-3 + higher terms.          (8)

If F_x(x,0) is not identically zero, A_F and B_F have EXACT pole order
three there and opposite leading coefficients, although F_t is a
polynomial. Their separate rationality must not be discarded.

For g>=5, put delta_lambda=lambda*partial_x+partial_t, lambda!=0, and
w_F=delta_lambda(F). Character projection gives EXACTLY

    Tr_(L/E)(w_F)=0,
    Tr_(L/E)(w_F*w_G)=g*lambda*(X_F*A_G+A_F*X_G). (9)

All other quadratic character sums are nonzero modulo g. For the actual
pair P,Q, B_F=K*X_F and Jac(P,Q)=1 imply

    X_P*A_Q-A_P*X_Q=1.                             (10)

Consequently the TWO-BY-TWO relative moment matrix in (9) has

    det(Tr_(L/E)(w_i*w_j))=-g^2*lambda^2.           (11)

Its determinant is a nonzero CONSTANT, while at least one diagonal entry
has a pole of order THREE at t=0. Indeed the constant Jacobian makes
(P_x(x,0),Q_x(x,0)) generate the unit ideal in C[x]. For a generic constant
linear combination F of P,Q, X_F(x,0) is not identically zero; (8) gives

    Tr_(L/E)(w_F^2)=2*g*lambda*X_F*A_F
                  with exact pole order three.     (12)

This pole is on a RETAINED divisor of V: the plane source is etale there,
so its t-valuation equals the intrinsic divisor valuation. It does not
mean that all source sheets over that curve are missing.

The small cyclic orders are also explicit. For w=X*lambda+A+B,

    g=2: Tr(w)=2*B,                                  pole 3;
    g=3: Tr(w^2)=6*(lambda*X*A+A*B),                  pole 6;
    g=4: Tr(w^2)=4*(2*lambda*X*A+B^2),                pole 6. (13)

For g=2 the quadratic trace is 2*((lambda*X+A)^2+B^2), with pole SIX:
A and B have opposite leading terms, whose squares ADD. Thus for EVERY
g>=2 a hypothetical actual pair has a generic first velocity with a
relative QUADRATIC moment that is not regular on V. It follows that this velocity is not integral
over C, and therefore not integral over B. This is an explicit obstruction
a claimed integrality proof MUST overcome; it is not an actual plane
counterexample and does not exclude the target pair by itself.

In the exceptional transverse direction lambda=0, all relative powers
through j<ceil(g/2) vanish. The first possible power j0=ceil(g/2) is

    g even: Tr(F_t^j0)=g*B_F^j0;
    g odd:  Tr(F_t^j0)=g*j0*A_F*B_F^(j0-1),          (14)

with exact pole order 3*j0 when F_x(x,0)!=0. A short unmarked prefix may
therefore hide the obstruction in that direction; the generic mixed
direction (9) detects it quadratically.

## 3. Full finite markers survive passage from the relative field

Let S_E be the full finite normalization of the original B in E.
It is a finite normal B-algebra contained in C; it is NOT identified with C.
The open etale V embeds in Spec(S_E) by Zariski's Main Theorem. If r in E
has a pole on a divisor of V, then r is not in the codifferent of S_E/B:
that divisor is unramified, so its threshold is zero in (5).

Therefore there is an integral marker s in S_E with Tr_(E/F)(s*r) not in B.
This is simply the defining trace-dual condition, without an informal
no-cancellation assertion. Apply it to r=Tr_(L/E)(w_F^j) from (12)--(13) and the g=2 quadratic formula, always taking j=2.
Trace transitivity gives

    Tr_(L/F)(s*w_F^2) not in B.               (15)

The marker lies in the ORIGINAL S as well, since S_E subset S. Neither
this argument nor (3) deletes the original ramified boundary primes.
They may lie outside W and outside the relative source image.

Thus an actual all-degree theorem making the FULL S-marked quadratic
first-velocity traces regular would close the generic integrality gate
and exclude this entire cyclic target family. Equations (8)--(15) show
why unmarked relative traces or unit moment determinants cannot supply it.
This is a specific finite-lattice/marked-energy obligation, rather than a
renewed mixed 2+3 search. PF142's minimal-degree descent remains current.

## 4. An explicit constant-Jacobian open control and its full normalization

There is an unconditional local-global control for (11). It is on a
Laurent OPEN source, not the whole plane. Use the same a0,sigma,b and set

    p=(x+sigma/a0)/(g*t^(g-1)), q=t^g.              (16)

On R_open=C[x,t,1/(t*b(t^g))], (p,q) has Jacobian ONE.
Its function field extension has degree g. The FULL finite normalization
of C[p,q] is C[p,t], q=t^g. The omitted prime t=0 has index g; primes
over every root of b(q) have index ONE and are also retained in the ledger.

For EVERY g>=5 and k,b as above, (9) gives exactly

    (Tr(w_i*w_j)) =
       [ -2*lambda*(g-1)*p/q,  g*lambda ]
       [  g*lambda,            0       ].          (17)

The determinant is -g^2*lambda^2, but the first entry is polar at q=0.
The integral marker t^4 in the FULL C[p,t] detects a pole of order TWO:

    lim_(q->0) q^2*Tr(t^4*w_p^2)
       =4/(g*b(0)^(2*k-2)) !=0.                    (18)

There can be an additional lower simple pole when g=6, which does not
alter (18). Formula (2) agrees: q_source=2, e=g, m=1+ceil(5/g)=2.
Every other omitted b-root prime has q_source=k, e=1 and marked charge
2*k+2. Nonconstant source units and the denominators in (16) prevent
this example from being a plane Keller counterexample.

## 5. Evidence and remaining work

The unbounded statements follow from the written valuation, duality and
character arguments. The new exact checker supports local thresholds,
finite cyclic projection, the constant-determinant identity, explicit
global cyclic functions, and the Laurent control with its integral marker.
PF141, PF142, MS0, LO1, LY0 and MX0 remain sealed; no predecessor checker
or expensive search is replayed. Their source/classification assurance
limits are retained, with no new external theorem import.

The outstanding step is a proof of the regularity in (7) from the ACTUAL
whole-plane target coframe, or another exclusion of the cyclic target.
Exact volume, affine primitive descent, Picard classes and unit partial
Gram determinants do not establish it. The original full normalization,
all ramified primes, non-equivariant cases and full trace-dual multiplication
remain in scope. Formal verification, independent whole-proof replay,
external review and full-goal completion are not claimed.
