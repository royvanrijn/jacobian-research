# Branch clearing and Haar automorphisms of the exceptional trace models

15 November 2026. Internal written proofs, with exact algebraic controls.
No novelty, independent specialist replay, external review or formal
verification is claimed. The [checker](verify_branch_clearing_rigidity.py),
[receipt](branch_clearing_verification.json) and
[source audit](branch_clearing_source_audit.json) delimit the assurance.

The [quadratic trace theorem](QUADRATIC_TRACE_MODELS_AND_HAAR_RESIDUES.md)
provides the full algebra and normalized functional used here. Two new
all-degree conclusions constrain the next construction:

* The sharp universal branch-clearing power always produces an unstable,
  pure-null, eventually mixed-safe phase. It cannot produce an MZ failure.
* Every automorphism of the quadratic ambient algebra preserving its Haar
  trace preserves the exact admissible subalgebra. It cannot turn any
  nonadmissible expression into an admissible one.

These statements cover arbitrary expressions, rather than selected old
circuits. They do not classify all admissible pure-null phases. Both
exceptional MZ decisions remain **OPEN**.

## 1. Exact setup and the inherited torus

For Gamma=2O or2I, set

    (n,m,r)=(8,12,18) or(12,20,30), r=n+m-2,
    A_*=C[X,Y,F,G]/(XG-YF-1),
    K_O=F^3 G/36-3G^3, K_I=12F^5+G^3/144,
    R=A_*[eta]/(eta^2-K), delta=F partial_X+G partial_Y,
    B=C[F,G,eta]/(eta^2-K), U=Spec B minus{0}.

The map pi:Z=Spec R -> U is an additive torsor, with slices X/F and
Y/G. Write S for the image of C[SL_2/Gamma] in R. Precisely,

    S={a in R:delta^j a belongs to eta^j R for all j>=1}.       (1)

The normalized functional is

    L(a_0+eta a_1)=h_*(a_0),                                 (2)

where h_* is ordinary normalized SL_2 Haar. Its restriction to S is
exactly quotient Haar. Every phase and marker in the original question
belongs to S, not just R.

The original left diagonal torus on SL_2/Gamma induces the following
integer weights on R:

    wt(X)=m, wt(Y)=n, wt(F)=-n, wt(G)=-m, wt(eta)=-r.          (3)

Both defining relations are homogeneous. This follows also from the
covariant map: eta X=a_f has weight2-n and eta Y=a_g has weight2-m.
The derivation delta has weight-(n+m); delta/eta on S is the original
left unipotent derivation, of weight-2. Thus weights in (3), when
restricted to S, are actual original left-group weights. They are
not an artificial grading inferred from the trace formula.

For nonzero a in R, write dep(a)=max{j:delta^j a!=0}, and let W(a)
be its maximum weight in (3). Put M=max(n,m) and l=min(n,m).

**Filtered weight bound.** For every a in R with dep(a)<=d,

    W(a)<=M d.                                               (4)

Proof. On A_* put E=X partial_F+Y partial_G. The standard SL_2
decomposition spans A_* by E^k q, where q is homogeneous of ordinary
degree N in F,G and0<=k<=N. The depth filtration restricts this to
k<=d. This is a direct weight decomposition for the left SL_2 ladder;
there is no cancellation of different ladder depths. Decompose q
further into monomials. Its weight is at most-l N, and E has weight
n+m. Therefore

    wt(E^k q)<=-l N+(n+m)k<=M k<=M d.

Every a=a_0+eta a_1 has both coefficient depths<=d. The extra eta
lowers weight by r, so the same bound holds on R. Equality is attained
by X^d in the two cases here, where m=M.

## 2. Sharp universal clearing is always safe

**Theorem.** Let0!=a in R have depth d>=1. Then eta^d a belongs to S,
and every one of its original left weights is at most

    -(l-2)d = -6d on2O, or-10d on2I.                         (5)

Consequently its inverse image P in C[SL_2/Gamma] is in the original
left-group nullcone, h(P^j)=0 for every j>=1, and for every regular
marker Q,

    h(Q P^j)=0 whenever j>W(Psi^*Q)/((l-2)d).                 (6)

A negative maximum marker weight already forces vanishing for all
nonnegative j. Formula (6) is a sufficient weight cutoff; no optimal
marked cutoff is asserted.

Indeed delta^j(eta^d a)=eta^d delta^j a is divisible by eta^j for
j<=d, and is zero for j>d. This proves exact admissibility in (1).
The weight bound is

    W(eta^d a)<=Md-rd=-(l-2)d<0.

A one-parameter subgroup, with parameter tending to infinity, sends
P to zero. All nonzero weights in QP^j are negative under the condition
in (6); Haar kills them by left torus invariance. This proves all
pure and marked assertions without a moment prefix.

The exponent d is **sharp as a uniform clearing exponent on the whole
depth<=d filtration**. For a=X^d,

    delta^d(eta^N X^d)=d! eta^N F^d.

F is nonzero on eta=0: F=eta=0 would force G=0, contradicting the
determinant relation. At every generic branch component F is a unit
and eta has order1. The last expression is divisible by eta^d only
if N>=d. Also eta^d X^d=(a_f)^d has weight exactly-(n-2)d, so the
uniform weight bound in (5) is sharp. Individual expressions can
have better divisibility and require smaller clearing powers. This
theorem does not classify those individual minimal powers or their
moment behavior. It treats every universal clearing eta^N a with
N>=dep(a), all of which have the same strictly one-sided conclusion.

For d=0 the expression lies in B. Under (3), each element of B minus
its constant term has strictly negative weights; this gives the
same pure-null and safety conclusion. A nonzero constant is not pure-null.

For example, the old target Hopf and normalizer phases have depths6
and20. Clearing by eta^6 or eta^20 makes each admissible, but forces
strict weight gaps at least36/120 on2O and60/200 on2I. This is a safe
regularization of each expression, not a transfer of its old bad tail.

There is also no denominator that clears **every** ambient expression:

    {a in R:a R subset S}=0.                                (7)

To prove this, work at the generic point of any branch component of
U, where its local ring O is a discrete valuation ring and eta is a
uniformizer. Since F is a unit there, the two algebras on the additive
chart are

    R_O=O[s], S_O=O[eta s], s=X/F.

For any nonzero a, let its s-degree be e and let its leading
coefficient have finite valuation v. The leading coefficient of
a X^N has the same valuation v, while membership in O[eta s] would
require v>=e+N. Large N violates this. Thus no nonzero a can satisfy
aR subset S. This is the zero conductor of this birational inclusion;
it does not contradict finite clearing on each depth filtration.

## 3. The period paired with a nowhere vanishing volume

On Z, the form

    Omega_Z=p^*Omega_*/eta                                  (8)

extends across the branch and is nowhere zero. In local coordinates
this is the residue form for eta^2-K=0. Smoothness of U ensures dK
is nonzero on its branch curve, so dividing the branched pullback
volume by eta cancels its simple zero. The additive torsor charts
then prove the global assertion.

Every unit of R is constant: restrict a unit to each additive fibre,
where a polynomial unit is constant. It descends to a unit on U;
normality across its missing codimension-two point identifies these
with units of the positively graded ring B, which are C^*.

Let T be the period functional for Omega_Z on H_3(Z;C), normalized by
T(eta)=1. This is possible since [eta Omega_Z]=[p^*Omega_*]!=0.
The cohomologically trivial deck involution negates Omega_Z, so

    T(a_0+eta a_1)=h_*(a_1), L(a)=T(eta a).                  (9)

In particular T(1)=0; T itself is not normalized Haar. The bilinear
pairing T(ab) on R is nondegenerate. If a_0!=0, test against eta b
and use nondegeneracy of h_*(a_0 b). If a_0=0 and a_1!=0, test against
b instead. Nondegeneracy of ordinary Haar follows, for example, by
testing a nonzero polynomial against the regular polynomial extending
its complex conjugate on SU_2: the integral of its squared absolute
value is positive. This argument applies to arbitrary complex
coefficients and does not impose a real-form condition on the phase.

## 4. Haar-preserving automorphisms preserve the branch

Let alpha be any C-algebra automorphism of R. Since Omega_Z generates
its canonical module and R^*=C^*, its pullback is u Omega_Z for some
u in C^*. The one-dimensional H^3(Z;C) is acted on by a nonzero scalar.
Consequently there is c in C^* with

    T(alpha a)=c T(a) for every a in R.                       (10)

Suppose first that L alpha=L. From (9) and (10), for every a,

    0=T((c alpha^(-1)(eta)-eta) a).

Nondegeneracy gives alpha(eta)=c eta. Conversely, if alpha(eta)=t eta,
then (9)-(10) give L alpha=(c/t)L. Evaluate at1 to obtain c/t=1.
We have proved the exact intrinsic characterization

    Aut(R,L)={alpha in Aut(R):alpha(eta) in C^* eta}.          (11)

Equivalently these are the automorphisms preserving the reduced branch
divisor eta=0. This conclusion concerns all algebraic automorphisms,
including nonlinear ones. It is not limited to the visible additive
shears or a chosen finite list of generators.

## 5. Positive-genus fibres force exact admissibility preservation

For a nonzero generic value b of eta, the curve under its additive
torsor is

    C_b={K(F,G)=b^2}.

Its smooth projective completion has genus3 on2O and4 on2I.
Here is an elementary calculation. On2O,

    F^3=108G^2+36b^2/G.

The cyclic degree-three map to the G-line has three simple zeros,
a pole of order1 at0 and a pole of order2 at infinity. All five
places are fully ramified. Riemann--Hurwitz gives
2g-2=-6+5*2=4. On2I,

    G^3=144b^2-1728F^5.

There are five simple zeros and a pole of order5 at infinity, all
fully ramified. Thus2g-2=-6+6*2=6. The assertion also holds over the
algebraic closure of the generic eta field. The checker verifies the
equations, squarefreeness and ramification arithmetic; the universal
curve and automorphism arguments remain written proofs.

An A^1 cannot map nonconstantly to either affine curve: extend such
a map to P^1 and the smooth projective completion, then apply
Riemann--Hurwitz in characteristic zero. Therefore every isomorphism
between these generic additive-torsor surfaces maps each additive
fibre into an additive fibre. This does not require a classification
of the automorphism group of the surface or a Makar--Limanov claim.

Now alpha in Aut(R,L) scales eta by (11). Its induced isomorphism of
generic eta fibres must preserve their curve fibration. Hence
delta(alpha b)=0 for every b in B, first generically and then
identically. Applying the same argument to alpha^(-1) proves

    alpha(B)=B.                                             (12)

The conjugate alpha delta alpha^(-1) is then a vertical vector field
for Z->U. Its relative tangent bundle is freely generated by delta,
because this is an additive torsor. The conjugate field has no zeros,
as delta has none. Its coefficient is therefore a unit of R, so

    alpha delta alpha^(-1)=v delta, v in C^*.                 (13)

Equations (11),(13) preserve every derivative divisibility in (1).
Thus

    alpha(S)=S for every alpha in Aut(R,L).                  (14)

In particular each such automorphism descends to an automorphism of
SL_2/Gamma and preserves its normalized Haar. No Haar-preserving
ambient automorphism can make a nonadmissible phase **or marker**
admissible. This closes that proposed repair mechanism for every
ambient expression, not just the old Hopf and normalizer circuits.
It does not assert that every automorphism downstairs lifts to R.

Combining (5) and (14), any expression first cleared by a universal
branch power and then transformed by a Haar-preserving ambient
automorphism is still pure-null and eventually mixed-safe for every
admissible marker. Its image lies in the automorphism saturation of
the original nullcone, a previously established safe class. Such
transformations cannot turn this universal clearing procedure into
an exceptional MZ failure.

## 6. What remains

The full exceptional question is still the constrained trace problem
(1)-(2). Individual admissible expressions can have lower branch
valuation in some derivative coefficients than universal clearing
supplies. They need not lie in the original left nullcone or in a
proved exhaustive safe normal form. No exhaustion by automorphisms
is established here. Regular multivariable constructions must meet
the exact derivative tests through genuine coefficient cancellation,
or a positive proof must classify all such admissible pure-null pairs.
No new degree, support or moment-prefix sweep is proposed.
