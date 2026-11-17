# Algebraic flow descent and integral source coordinates

HC4MYKX0/1. Full HC4/JC2 remains OPEN. This note proves an all-degree
obstruction to EVERY locally finite source derivation descending to a
noninvertible plane Keller target field. It also computes the exact global
translation multipliers available on the finite normalization when an
integral source coordinate is supplied. Neither that coordinate nor descent
of the resulting actions is asserted automatic.

## 1. No nonzero locally finite source field can descend: KX0

Let B=C[P,Q] subset R=C[x,y] be an ACTUAL complex plane Keller inclusion,
F=Frac(B), L=Frac(R). For every NONZERO locally finite C-derivation delta
of R,

    delta(F) subset F ==> B=R.                     (1)

Local finiteness means every polynomial has a finite-dimensional span of
its iterated derivatives. Thus, under noninvertibility, the entire Lie
algebra of source derivations preserving F has ZERO intersection with
the nonzero locally finite derivations. In particular this excludes:

* EVERY nonzero locally nilpotent source derivation;
* EVERY polynomial replica p(f)*partial_g, p!=0, in ANY polynomial
  source coordinates (f,g), including multipliers not in the target field;
* every nontrivial connected linear algebraic group action on the source
  which preserves F.

The existing PF2RULED52 theorem excludes unmultipled source translations
and integer diagonal Euler fields in arbitrary source coordinates. The
present argument closes the multiplier loophole and passes from those
normal forms to ALL locally finite fields. It does not prove that any
one of them preserves F. The resulting inverse criterion is not itself
reported as a new exclusion of otherwise unqualified Keller packets.

First, the retained quasi-finite intersection lemma gives

    F intersect R=B.

Hence delta(F) subset F implies delta(B) subset B. Conversely preservation
of B implies preservation of F. If delta kills B, it kills F and therefore
L, by finite separability, so delta=0. A nonzero descending derivation
thus induces a nonzero target derivation.

### Every translation replica is excluded, without target linearization

Suppose delta=p(f)*partial_g in actual polynomial source coordinates, with
p!=0, and delta(B) subset B. It is locally nilpotent on R and on B, so
its finite exponentials give compatible algebraic Ga actions and the
Keller map is equivariant. For any c with p(c)!=0 the ENTIRE source line
f=c is one transitive Ga orbit. The source vector field is nonzero there;
the Keller differential makes the target vector field nonzero at its
image. A nontrivial Ga orbit in characteristic zero has trivial stabilizer:
the only proper algebraic subgroup of the additive line is {0}.
Consequently the restriction of the Keller map to this source line is
INJECTIVE. The already retained Gwozdziewicz one-line theorem forces the
inverse, including after the polynomial source-coordinate change.

This proof DOES NOT divide by p(f) in F, does not assume p(f) is a target
function, and does not linearize an action on the target. Merely knowing
the unmultipled translation cannot descend would not justify (1) for
these replicas. The orbit argument supplies the missing step.

### Every locally nilpotent field reduces to a replica

We explicitly retain the classical Rentschler input already used in GV0,
as stated in Furter Lemma4.2: a polynomial Hamiltonian whose derivation
is locally nilpotent is a polynomial in a source coordinate. This is an
imported classification theorem, not a consequence of finite tests here.

For clarity, it suffices in exactly that form. If delta is any LND, its
polynomial Ga automorphisms have Jacobian1: the Jacobian is a unit in
C[t,x,y] and equals1 at t=0. Differentiating gives div(delta)=0.
Polynomial Poincare integration writes delta as a Hamiltonian derivation
D_a. The retained criterion gives a=H(f) for a polynomial coordinate f.
In complementary coordinates (f,g), delta is therefore p(f)*partial_g.
The preceding argument proves (1) for every LND, without importing
any new target-action classification.

### The locally finite case includes all semisimple weights

Here are the Jordan steps, so that no extra unproved decomposition
assumption is hidden. On each finite delta-invariant space the ordinary
linear Jordan decomposition is compatible with inclusions. It defines
commuting operators delta_s,delta_n on R. If a,b are generalized
eigenvectors of eigenvalues alpha,beta, the identity

    (delta-alpha-beta)^N(ab)
       =sum_i binomial(N,i)*(delta-alpha)^i(a)
                            *(delta-beta)^(N-i)(b)

puts their product in the generalized eigenspace alpha+beta. Thus delta_s
acts there by alpha+beta and is a derivation; delta_n=delta-delta_s is a
locally nilpotent derivation. For any b in B its finite orbit space lies
in B, and both Jordan parts are polynomials in delta on that space.
Therefore BOTH parts preserve B. These standard steps are also proved
in Flenner--Zaidenberg Lemma2.2; the above argument is written directly.

If delta_n!=0 the LND result proves the inverse. Otherwise delta is a
nonzero semisimple derivation. Choose finitely many eigencomponents of
the source generators; they generate R. Their finitely many eigenvalues
generate a free abelian subgroup Gamma of (C,+). Grading by Gamma defines
a finite-dimensional algebraic torus action on R. Spectral interpolation
puts every eigencomponent of every b in B back in B, so B is also Gamma
graded. Choose a homomorphism Gamma->Z nonzero on an occurring nonzero
weight. It gives a NONTRIVIAL algebraic Gm action preserving B.

Use the Gutwirth plane linearization theorem in the already retained
Flenner--Zaidenberg Corollary4.6 (PF2GRADE1). After a polynomial source
change its generator is a nonzero integer Euler field. PF2RULED52 now
gives the inverse. This handles arbitrary complex semisimple eigenvalues,
not just commensurable integer weights. No new classification theorem
is imported in this note; the Rentschler and Gutwirth imports remain
explicit and are not independently reproved.

Finally a rational action of a connected linear algebraic group is locally
finite on its coordinate ring. Every Lie algebra element acts by a locally
finite derivation on R; if the group preserves F, it preserves B by the
intersection lemma and so do these derivations. A nontrivial connected
action in characteristic zero has a nonzero infinitesimal action. Equation
(1) proves the group assertion, without a classification of the group.

## 2. Exactly which translations extend to the full normalization: KX1

Let S be the FULL finite normalization of B in L. Suppose f is an ACTUAL
POLYNOMIAL SOURCE COORDINATE and is integral over B, equivalently f in S.
Write R=C[f,g]. Integrality of one coordinate is a hypothesis, NOT proved
to imply inversion here. PF2RULED40 already retains it in a common-ruling
sector, with genuine affine non-Keller completion controls; the one-
coordinate hypothesis must not be mistaken for the source-finite generic
linear-coordinate property of PF2OMIT3.

Every omitted prime E_i of S is VERTICAL over f. To see this, localize S
at C[f] minus0. The resulting normal affine curve over C(f) contains
A1_(C(f),g) openly. Its smooth projective completion is P1_(C(f)); the
only possible added point is its single rational infinity. Adding that
point would give the complete P1, which is not affine. Thus no omitted
prime dominates f. Write

    f|E_i=c_i in C, r_i=v_(E_i)(f-c_i)>0,
    e_i=the ramification index over the target.

There are finitely many such primes and boundary values c. Define

    a_c=max_(c_i=c) max(0,ceil(e_i/r_i)-1),
    p_min(f)=product_c(f-c)^a_c.                    (2)

Then, for every p in C[f],

    p(f)*partial_g preserves S
          iff p_min divides p in C[f].             (3)

Every nonzero such restriction is a GLOBAL LND on S with kernel EXACTLY
C[f]. This is the complete multiplier ideal, not merely a sufficient
large power, a generic boundary check, or an assertion that codimension-
two singular points can be ignored without normality.

Proof. The actual source volume Omega=df wedge dg has ramification order
e_i-1 along E_i by KN0. At its generic smooth point use parameters u,v:

    f-c_i=u^r_i*a(u,v), Omega=u^(e_i-1)*b(u,v)du wedge dv,
    a,b units.

The Hamiltonian identity i_(partial_g)Omega=-df makes the tangential
coefficient of partial_g have valuation r_i-e_i and its normal coefficient
have at least r_i-e_i+1. Characteristic zero and a being a unit make the
first order EXACT. Multiplication by p adds
r_i*ord_(f=c_i)(p), so regularity of the vector field along E_i is exactly

    r_i*(ord_c p+1)>=e_i.                         (4)

This local threshold is already retained in PF2RULED99 (and e=2 in
PF2RULED40); it is credited rather than claimed new. On source divisors
the vector field is polynomial automatically. Since S is NORMAL, an
element belongs to S iff it has no pole at any height-one prime. Applying
this to the derivative of each s in S proves the SUFFICIENCY GLOBALLY,
including all singular finite points. Conversely preservation of S
preserves its height-one localizations; the exact tangential coefficient
proves necessity. Conditions(4) at all primes are precisely (2)--(3).
Local nilpotence restricts from R, and its kernel is S intersect C[f]=C[f].

The DIVISOR of the extended vector field is also explicit. On the
retained source line f=c its order is ord_c p. On E_i it is

    r_i*(ord_(c_i)p+1)-e_i.                       (5)

Equation(5) concerns codimension-one zero orders, not a claim to classify
all isolated fixed points on singular S. No smoothness of the entire
normalization or prescribed residue degree is assumed.

Under noninvertibility, NONE of these nonzero extensions can preserve
the target FIELD F, by KX0. Thus no polynomial multiplier repairs the
descent, even though such multipliers ALWAYS produce actual algebraic
actions on S under the integral-coordinate hypothesis. Finiteness of
S/B does not make its source action equivariant over an action on B.
This closes a whole family of descent shortcuts; it does not prove
that f exists or that any constructed action descends.

## 3. The dimension boundary is witnessed by a known whole Keller map

The KNOWN three-dimensional map from KO2/KW1 has the exact weights

    wt(x)=1, wt(y)=-1, wt(z)=-2,
    wt(A)=-2, wt(B)=-1, wt(C)=1.

Hence E=x*partial_x-y*partial_y-2z*partial_z is a NONZERO locally finite
source field preserving C[A,B,C] and its fraction field, while the map
is noninvertible. It is a genuine equivariant whole-Keller countercontrol,
not an unconstructed graded candidate. This map and its noninvertibility
are PRIOR credited results; only the weight identities are used here.

Therefore (1) is PLANE-specific. It cannot be applied directly to the
four-dimensional gradient in HC4. Indeed the retained cotangent form
under hypothetical HC4 failure has two liftable constant target directions
already classified by IT1/2. The result applies to its genuine plane
Keller base through Tang, with the existing assurance limits intact.

## 4. Full target and what still has to be supplied

FULL THEOREM TARGET: prove full HC4/JC2 or construct an ACTUAL counterexample.
CURRENT LAST UNSOLVED IMPLICATION remains

    ACTUAL polynomial plane source
       ==> Tr(Lambda_L*Lambda_L) subset B.          OPEN

The alternative full marked canonical-energy integrality gate remains
open. A source symmetry route must now supply ACTUAL field descent of
some nonzero locally finite field; existence of the source action, a
normalization action, its multiplier ideal or a formal lift is insufficient.
The theorem does not assert such descent automatic and is not a proof
of HC4. An integral source coordinate alone is NOT proved sufficient.
External review and outreach remain DEFERRED; no assurance flag is upgraded.
