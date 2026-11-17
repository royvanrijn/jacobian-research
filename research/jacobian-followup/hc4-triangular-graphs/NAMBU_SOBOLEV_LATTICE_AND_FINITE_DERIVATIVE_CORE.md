# One target-additive derivative detects the normalization

2026-11-16. HC4MYLX0/1. Complex polynomial PLANE Keller maps, ALL DEGREES.

FULL THEOREM TARGET: prove HC4/JC2 or construct an actual counterexample.
CURRENT LAST UNSOLVED IMPLICATION: ACTUAL polynomial plane source forces
Tr(Lambda_L*Lambda_L) subset B, or full marked-energy integrality; OPEN.
The previous goal turn made PROGRESS with LW0/1. External review/outreach
is DEFERRED. The parity source theorem remains unadopted.

This note tests the ACTUAL commuting Keller differential coframe against the
Bergman obstruction. One first derivative recovers the FULL normalization
inside its finite trace dual, and its defect is exactly the finite trace
discriminant module. All successive derivatives leave only target polynomials.
Arbitrarily long finite derivative prefixes still survive under putative
failure. Thus increasing derivative depth or invoking a common dense domain
does not supply the missing global estimate.

## 1. Setup, target additive fields, and the component-level source input

Keep B=C[P,Q] subset R=C[x,y], full finite normalization S, finite embedded
trace dual Lambda, and Omega=dP wedge dQ. Use LW0's bounded target ball V
meeting EVERY omitted image, U=K^-1(V), and H=A2(U). Then

    R intersect H=Lambda, omega_S=Lambda*Omega,
    div_Z(Omega)=sum_i(e_i-1)E_i.                  (1)

Let eta be ANY NONZERO locally nilpotent derivation of the TARGET polynomial
ring B. The already retained plane Rentschler input, used in KX0, gives a
polynomial target coordinate frame (P,Q) with

    eta=alpha(Q)*partial_P, alpha in C[Q] nonzero.

Its unique lift delta to the finite source field is a polynomial SOURCE
derivation: it is alpha(Q) times the Nambu derivative D_P. It preserves the
source volume, commutes with D_Q when alpha is constant, and has a polynomial
Hamiltonian H0(Q) with H0'=alpha. No local finiteness, completeness or algebraic
flow on the source is assumed. Every nonzero constant Nambu direction in ANY
polynomial target frame is included. Changing target coordinates leaves the
embedded B,S,Lambda unchanged; boundedness on the fixed ball is unaffected.

We use the PREVIOUSLY ADOPTED Chau no-line statement from JY0, not a new
external theorem. The Remark after Theorem 1 of *Two remarks on non-zero
constant Jacobian polynomial maps of C2*, Ann. Polon. Math.82(2003),43,
excludes an irreducible exceptional component isomorphic to A1. For a Keller
map the exceptional set is its nonproper set. Therefore EVERY omitted image
D_i has Q|D_i nonconstant: a component contained in Q=c would be an entire
coordinate line. In characteristic zero its derivative is generically nonzero.
Hence eta is GENERICALLY TRANSVERSE to EVERY omitted image. The finitely many
zeros of alpha(Q) do not change this generic statement.

Scope audit: Chau's Theorem 2 concerns the WHOLE exceptional curve and alone
does not exclude one line component. The component statement is the separate
Remark, following Lemmas 5 and 6(ii). The cached primary HTML was read through
those lemmas and their proofs; the published PDF pp.42–43, including the
Remark, was visually inspected. Its content agrees with the retained JY0
import. The other theorem's scope was not strengthened. No correction to a
registered claim or new primary cache is needed.

## 2. LX0: the EXACT first-derivative lattice is S

For EVERY eta above,

    {r in R : r in H and delta(r) in H}
         =Lambda intersect delta^-1(Lambda)=S.     (2)

If s belongs to S, then ds and dH0(Q) are global regular Kahler forms.
The Hamiltonian identity gives

    delta(s)*Omega=ds wedge dH0(Q).

The right side defines a regular reflexive canonical form, even at singular
finite points. Thus delta(S) subset Lambda, and S is in the left side of (2).

Conversely let r in Lambda have a pole of order q>0 at E_i. Such a pole can
only occur for e_i>=2, and q<=e_i-1. At a generic point use t=0 for D_i,
v=Q as a parameter on it, and a tame root coordinate u with t=u^e. Since
delta(t) is a unit generically and delta(v)=0,

    delta=(unit)/(e*u^(e-1))*partial_u,
    v_E(delta(r))=-q-e < 1-e.                     (3)

The coefficient -q is nonzero in characteristic zero. Thus delta(r) is not
in Lambda. Normality/reflexivity now give (2) GLOBALLY: no boundary prime or
finite singular point is dropped. Unramified omitted primes cause no gap,
since Lambda is already regular there.

There is consequently a canonical B-LINEAR INJECTION

    Lambda/S -> R/Lambda, [lambda] -> delta(lambda)+Lambda,

whose image gives the exact GLOBAL identification

    Lambda/S isomorphic (Lambda+delta(Lambda))/Lambda. (4)

It is well defined because delta(S) subset Lambda; Leibniz and delta(B)
subset B give B-linearity. Injectivity is exactly (2). Thus EVERY target
additive direction exhibits the SAME finite discriminant quotient as its
first canonical-lattice derivative defect. This is an equality of embedded
modules at all finite points, not only a generic length calculation.

The module Lambda/S is the cokernel of the finite trace-pairing map
S->S^vee. Locally its presentation is the trace Gram matrix, and its zeroth
Fitting ideal is the finite discriminant. Its vanishing is the retained
inverse criterion. Equation (4) supplies no automatic vanishing.

In particular, for ANY finite B-generators lambda_j of Lambda,

    all delta(lambda_j) in L2(U)
       iff Lambda=S iff the full Keller inverse exists. (5)

Indeed (2) puts every generator in S. Conversely an inverse makes the
preimage of the closed target ball compact, so all polynomial derivatives
have finite norm. This replaces LW's fourth-power requirement by an exact
FIRST-DERIVATIVE requirement; it does not prove that requirement automatic.

## 3. Every lifted target additive field has a NON-DENSE natural domain under failure

Let Dom(delta)={f in H:delta(f) in H}, with the maximal natural differential
operator. If the map is noninvertible, KN0 supplies a positive canonical
boundary prime E with e>=2. Choose a generic smooth point p over V and a
local polydisk. Every f in H has an L2-removable canonical form

    f*Omega=F_f(u,v) du wedge dv,

whose coefficient extends holomorphically to p. Evaluation ell_p(f)=F_f(p)
is a NONZERO bounded functional by LW1's mean-value argument and a global
canonical section. Write delta=A(u^e,v)u^(1-e)partial_u with A(0,v)!=0.
It preserves Omega, hence

    delta(f)*Omega=L_delta(f*Omega).

The coefficient on the right is

    A u^(1-e) F_u + partial_u(A u^(1-e))*F.

If F_f(p)!=0, its leading term is

    (1-e)A(0,v)F_f(0,v)u^-e,

which cannot be the coefficient of an L2 holomorphic top form. Therefore

    Dom(delta) subset ker(ell_p) proper CLOSED subset H. (6)

Its natural domain is NOT DENSE, for EVERY nonzero target additive field,
including EVERY nonzero constant Nambu direction in every target frame.
Conversely, under inversion target polynomials pull back to a dense set
in H contained in the domain. Thus density for even ONE such operator would
force the complete inverse. Density is UNSUPPLIED for actual Keller maps.

Formal commutation, polynomial coefficients, polynomial target slices and
divergence zero do not produce a dense common domain. No self-adjointness,
unitary representation, integration-by-parts boundary cancellation or
Nelson analytic-vector theorem is imported. Those are additional hypotheses,
not consequences of the algebraic coframe. Equation (6) identifies the
specific missing-end functional that any proposed cutoff/core argument
would have to overcome.

## 4. LX1: the largest target-additive-stable finite lattice is the target ring

More generally, if a FINITE B-submodule W of R is stable under ONE delta
as above, then

    W subset B.                                   (7)

No canonical-module hypothesis or full generic rank is required. This is
finite generation over B, not over C. The assertion concerns target additive
fields; arbitrary polynomial source fields do not satisfy it.

Proof. Set k=C(Q). The source Q is a polynomial submersion, so HX0 gives
k relatively algebraically closed in L=Frac R. Over k, W_k is a finite
k[P]-module stable under partial_P, since alpha(Q) is a nonzero scalar in k.
At ANY finite point of the normalized generic curve over A1_P, a finite
list of module generators gives a UNIFORM lower valuation bound for every
partial_P^n(w). Coefficients in k[P] have nonnegative valuation there.

After algebraically closing k, put t=P-c=u^e at such a point. A Laurent
Puiseux expansion w=sum a_j t^(j/e) cannot contain a negative integer power
or ANY noninteger power: its successive derivatives have nonzero falling
factorials and arbitrarily negative valuations. Nonnegative integral powers
differentiate regularly, and distinct powers cannot cancel. Consequently all
branches of w are regular and invariant under inertia at every finite P.
The intermediate field generated by w is unramified over the finite P-line.

There is no nontrivial connected finite cover of A1 in characteristic zero.
For completeness its projective completion has degree m, ramification only
over infinity and s points above infinity. Riemann--Hurwitz gives

    2g-2=-2m+(m-s)=-m-s,

so m+s<=2 and hence m=s=1. Thus w belongs to kbar(P). Relative algebraic
closedness from HX0 descends this to k(P): its regular constant extension is
linearly disjoint, and Galois descent fixes exactly k(P). Finally the retained
actual-source intersection R intersect Frac(B)=B puts w in B, proving (7).

The curve proof uses only bounded pole orders of all derivatives, not that
the derivatives lie in a particular codifferent. It is plane-only: HX0 for
one complementary component supplies the needed constant-field regularity.
No corresponding n>=3 statement or direct application to an HC4 gradient
is claimed; universal HC4 use remains via retained Tang.

In particular the polynomial all-derivative Bergman core for ONE delta is

    C_delta={r in R:delta^n(r) in H for EVERY n>=0}=B. (8)

To see this without an infinite-generation assumption, C_delta is a B-module
by Leibniz, is contained in finite Lambda, and is therefore finite over the
Noetherian ring B. It is delta-stable. Apply (7); B is contained in it because
eta is locally nilpotent on B. The largest delta-stable B-submodule of Lambda
is thus exactly B. Also

    delta(S) subset S iff inverse,
    delta(Lambda) subset Lambda iff inverse.        (9)

For the first implication use (7) with S; for the second use either (2) or
(7). All statements include EVERY supplied nonzero target additive direction.

Since target polynomials are dense in the Bergman space of a ball, (8) has
the exact Hilbert closure

    closure_H(C_delta)=C(A2(V)),

whose orthogonal projection is LW0's automatic Pi=CT/d. Under failure the
orthogonal complement is NONZERO: finite Lambda has a nonzero trace-zero
polynomial section. Thus the automatic Hilbert trace projects precisely onto
the closed base space containing all algebraic vectors with every derivative
integrable, not onto a dense common polynomial source core. No classification
of all nonpolynomial smooth vectors is asserted.

## 5. Arbitrarily long finite derivative prefixes in the SAME actual source

Under putative noninvertibility d>1. Choose 0!=s in S with Tr(s)=0 and let
H_D in B be the reduced equation of the ENTIRE nonproper curve. For EVERY
integer N>=0 put

    r_N=H_D^N*s.

Then r_N is a nonzero source polynomial outside B, with trace zero, and

    delta^j(r_N) in S subset Lambda subset H,
    Tr(delta^j(r_N))=0, for 0<=j<=N.               (10)

Indeed v_E(H_D)=e_E, and the local derivative bound is
v_E(delta(a))>=v_E(a)-e_E. Hence these valuations are at least
(N-j)e_E+v_E(s)>=0. All source divisors are regular because delta is a
polynomial source derivation. Normality gives the GLOBAL S membership,
including singular finite points. Trace commutes with the lifted derivation.
But (8) makes some later derivative fail to belong to H.

Thus no fixed finite derivative order, WITHOUT a degree bound or a genuine
stable finite-module assertion, can detect the entire obstruction in all
degrees. This is a conditional theorem inside ANY actual counterexample,
not a constructed Keller counterexample. Degrees of r_N grow. The module
condition in (5) is different: it controls ALL finite Lambda generators
and therefore their entire B-module, not one growing-degree element's prefix.

## 6. Sharp actual Laurent controls for directions and all derivative orders

Reuse LW1's explicitly constructed connected Laurent maps, without replaying
their registered checker:

    P=u^e, Q=v/u^(e-1), e>=2,
    R_e=C[u,u^-1,v], S_e=C[u,Q], Lambda_e=u^(1-e)S_e.

The target branch image is the LINE P=0; the source is Gm times A1.
In coordinates (u,Q) its two commuting target derivatives are

    D_P=(1/e)u^(1-e)partial_u, D_Q=partial_Q.

On the bounded target bidisk preimage, D_P has the same nondense domain
and first-derivative lattice S_e. In contrast D_Q PRESERVES Lambda_e;
its domain IS DENSE. The weighted Laurent monomials u^j Q^m,
j>=1-e,m>=0, are an orthogonal Bergman basis and their finite span is in
the D_Q domain. Its all-derivative polynomial core is Lambda_e, not B_e.
This explicitly shows why the whole-plane no-line and parameter regularity
hypotheses must not be dropped from a one-direction assertion.

There are NEW closed-form all-order prefixes. Set r_N=P^N*u. Then

    D_P^k(r_N)=(N+1/e)_(falling k)*u^(e*(N-k)+1).   (11)

It belongs to S_e for k<=N, to Lambda_e for k<=N+1, and fails L2 at the
EXACT next order k=N+2. Its trace and all those derivative traces are zero.
At integrable orders its squared norm is

    pi^2*|(N+1/e)_(falling k)|^2 /
                 (e*(N-k+1)+1).                   (12)

All e>=2 and N>=0 are covered. The next falling factorial is nonzero,
so the divergence is genuine; this is not cancellation of a trace.
These are actual Laurent controls with a constant source Jacobian and
finite codifferent. They are NOT whole-plane Keller or HC4 counterexamples.

## 7. Source audit, verification, and the unresolved global estimate

LW's canonical L2 and bounded evaluation arguments, KU's finite duality,
HX0's relative algebraic closedness, the retained Rentschler theorem, and
JY0's Chau no-line import are the exact structural inputs. The component
source remark was checked separately from the whole-curve theorem. No new
external theorem is adopted or primary file cached. Riemann--Hurwitz and
the Puiseux/falling-factorial proof are written; no analytic operator theorem
is substituted for its dense-domain hypotheses.

New exact controls check the new derivative/canonical exponents, exact finite
defect kernel in the tame model, all-order induction and norm thresholds,
and Laurent tangent/transverse distinction. They do not replay ANY registered
mathematical checker or certify the universal geometric/analytic proofs.
Independent replay, formal verification and external review remain false.

The remaining universal task is to force the finite discriminant quotient
Lambda/S to vanish, or equivalently obtain the actual first-derivative/L4
estimate. The actual polynomial commuting coframe already has the exact
defect (4); its formal algebra does NOT supply a dense common core. Finite
derivative prefixes can miss the obstruction to arbitrary depth by (10).
No automatic Sobolev estimate, ordinary trace-pairing regularity or full
marked-energy integrality is proved. Full HC4/JC2 remains OPEN and the goal
remains active.
