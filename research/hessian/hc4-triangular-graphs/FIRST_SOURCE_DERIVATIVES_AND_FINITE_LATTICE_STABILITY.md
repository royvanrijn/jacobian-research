# First source derivatives detect the full additive trace-lattice stabilizer

HC4MYMS0. 16 November 2026. All actual complex polynomial plane Keller
maps, all degrees, every polynomial source frame and every source LND.
Full HC4/JC2 remains OPEN. Stability/integrality is not proved automatic.
No novelty, independent replay, formal verification or review is asserted.

The preceding MR0 theorem recovers the full field from a bare translation's
finite derivatives. Here the global finite-module condition is reduced to
integrality of just TWO source polynomials. We also identify the normalized
directional saturation as an actual open of the FULL finite normalization.
This supplies a global algebra argument under an explicit premise; it does
not derive that premise from field reconstruction.

## 1. The first-derivative criterion

Let B=C[P,Q] subset R=C[x,y] be ACTUAL plane Keller, with Jacobian one
after a constant target scaling. Let F=Frac B, L=Frac R, d=[L:F], S the
FULL finite normalization of B in L, and Lambda its FINITE trace-dual.
All ramified and unramified omitted primes are retained. Fix any nonzero
source locally nilpotent derivation. In its retained Rentschler frame write

    delta=p(y)*partial_x, p!=0,
    a=delta(P), b=delta(Q), H'(y)=p(y).

The following conditions are EQUIVALENT:

1. delta(Lambda) subset Lambda.
2. delta(S) subset S.
3. BOTH a and b are integral over B, equivalently a,b in S.
4. The finite source-orbit algebra

       J_delta=C[delta^j P,delta^j Q:0<=j<=M],
       M=max(deg_x P,deg_x Q),

   is a FINITE B-MODULE.
5. For EVERY 1<=j<=d, BOTH Tr_(L/F)(a^j) and Tr_(L/F)(b^j)
   belong to B.

The zero derivation satisfies the same equivalences trivially. For nonzero
delta, existence of the frame is the already retained classification, not
a new theorem inferred from the checker.

### Proof of the new implication 3=>2

Assume a,b in S. The actual volume Omega=dP wedge dQ=dx wedge dy gives

    dH = i_delta Omega = a*dQ-b*dP.                         (1)

The right side is a GLOBAL regular Kahler one-form on Spec S. At a generic
omitted prime a pole of H of order q>0 would give a du coefficient of dH
of exact order -q-1 in characteristic zero, impossible. H is already a
source polynomial. Normal Krull intersection therefore puts H in S.
The polynomial equation H(T)-H(y)=0, scaled to be monic, makes y integral
over S. Normality gives y in S.

KX1 now applies at its full stated scope: EVERY omitted prime E_i is
vertical over this ACTUAL source coordinate. Write

    y|E_i=c_i, r_i=v_Ei(y-c_i)>0,
    e_i=target ramification index, ell_i=ord_(c_i) p,
    k_i=r_i*(ell_i+1).

At a generic smooth point of E_i use coordinates u,v. The KN0 volume and
KX1 calculation give

    Omega=u^(e_i-1)*unit*du wedge dv,
    delta's tangential coefficient has EXACT order k_i-e_i,
    delta's normal coefficient has order >=k_i-e_i+1.       (2)

The finite map E_i->its target curve is nonconstant and separable. Hence
at least one of P|E_i,Q|E_i has nonzero derivative in v. For that function,
the tangential term in delta(P) or delta(Q) has exact order k_i-e_i;
the normal contribution has strictly higher order. If k_i<e_i this is
a pole, contradicting a,b in S. Thus k_i>=e_i for EVERY omitted prime.
The retained KX1 necessity-and-sufficiency theorem now gives delta(S)
subset S, including singular finite points by normality. This is a global
argument, not an assertion that arbitrary derivations extend integrally.

Conversely 2=>3 is immediate. LO1 proves 1<=>2 using
omega_S=Lambda*Omega and (Lambda:Lambda)=S; every source LND preserves
Omega. Under 2, all source iterates lie in finite S, and the B-submodule
J_delta is finite because B is Noetherian. This proves 2=>4; 4=>3 is
integrality in a finite algebra. No field-generation premise is needed
for these four equivalences, so proper replica-stable factors are included.

For 3<=>5, use the characteristic polynomial of multiplication by a on
the d-dimensional F-space L. Newton identities in characteristic zero
recover its coefficients from its first d power traces. Regular traces
make this polynomial monic in B[T], proving integrality. Conversely integral
a has integral characteristic coefficients and traces, which lie in the
normal ring B. Repeat for b. Trace ONE, trace TWO, or a bounded prefix
independent of d is not substituted for the stated first d tests.

Equivalently, 3 forces an integral source coordinate and the exact full
boundary inequalities r_i*(ell_i+1)>=e_i. These are not automatically true.
Under noninvertibility all nonzero passing source LNDs have LO1's ONE
invariant coordinate pencil. Passing tests for two distinct invariant
source fields force the inverse. In particular integrality of all FOUR
entries P_x,Q_x,P_y,Q_y is equivalent to inversion: (1) gives both x,y
in S, so R is finite over B. No claim that those four entries are integral
is made. A single passing source LND still need not imply inversion.

## 2. The normalized finite orbit is the maximal additive open in the full normalization

Suppose additionally Frac J_delta=L. MR0 proves this for EVERY bare source
translation in EVERY polynomial source frame, and supplies other sufficient
cases via gcd(d,deg p+1)=1 or trivial decks and a monomial replica.
Let T be the normalization of J_delta in L and Z=Spec T.

Then Z is an AFFINE OPEN SUBSET of N=Spec S containing the WHOLE specified
source X=A2. It is the LARGEST open of N containing X on which the given
source Ga action extends as an algebraic Ga action. Its omitted primes in N
are detected already by the first two derivatives a,b.

Proof. Normalization is finite for the finitely generated complex algebra
J_delta. Normality of R gives T subset R. The map X->Z is birational and
quasi-finite, because its fibres are subsets of Keller fibres. The retained
normal-target Zariski-Main argument makes it an OPEN IMMERSION.

The Ga coaction on J_delta extends to T. Indeed its field coaction sends
an element of T to an element of L(t) integral over J_delta[t]. The normal
ring T[t], with fraction field L(t), contains that image. Group-law identities
extend from the dense field/source. Thus this is a polynomial action, not
merely a derivation or a formal exponential. It preserves the given X and
every boundary prime of Z, since the group is connected.

Also S subset T: an element integral over B subset T is integral over
normal T. Hence Z->N is a birational morphism. We show it is quasi-finite.
Any curve contracted by Z->Y must lie in Z minus X, because K is quasi-finite
on X. It must therefore be a whole boundary prime E. The Ga action preserves
E, so if P,Q are constant on E, all their positive delta iterates restrict
to zero on E. EVERY generator of J_delta is then constant on E. This would
contract E under the FINITE normalization map Z->Spec J_delta, impossible.
Thus no positive-dimensional fibre exists. Z->N is quasi-finite and birational,
and normal-target Zariski Main again makes it an open immersion.

LO0 applied to this normal affine whole-plane enlargement puts y in T.
Every retained boundary prime is therefore vertical over y. At a prime of
N where y has a pole, (1) makes at least one of a,b polar: H(y) has a pole
and its differential cannot be regular there. At a prime where y is regular,
y is constant on that prime. One can see this by localizing S[y], then
normalizing: its normal affine spectrum contains the whole X and the regular
coordinate y, so KX1's elementary affine generic-fibre argument applies.
Equivalently, a horizontal finite y and a pole of x would give a pole of
dx wedge dy, contradicting regularity of the actual target volume.

At such a vertical prime E, (2) shows that a,b are both regular EXACTLY when

    r_E*(ord_c p+1)>=e_E.                                  (3)

If (3) holds, delta preserves the DVR O_(N,E), so ALL iterates are regular
there, and T subset O_(N,E). Since Z is an affine open of N, its separated
valuative centre then contains the generic point of E. If (3) fails,
one first derivative has a pole and E cannot be retained.

The complement N minus Z has no isolated component outside these curves.
Otherwise normal Hartogs extends every function of affine Z across such
a point, giving a morphism back to Z over N and contradicting its omission.
Therefore the complete description, including finite points, is

    Z=N minus the UNION OF WHOLE primes where a or b is polar. (4)

To see maximality, take any Ga-invariant open W of N containing X. All
delta iterates of P,Q are regular on W. Integral elements of J_delta are
regular on normal W, so T subset Gamma(W,O). This gives a morphism W->Z
over N, and hence W subset Z. W need not be affine for this argument.

In particular T=S iff condition3 in Section1 holds. Field reconstruction
alone gives (4), not T=S. The normalizations of a raw first-jet order and
of the full finite orbit are not silently identified; no assertion that the
raw J_delta itself is normal or injective is needed.

## 3. Why the Mathieu resolvent criterion does not automatically prove this premise

The user's coordinator working note is credited:
[two explicit bridges](../../jacobian/cross-programme-connections/LINEAR_MARKER_SYMMETRY_AND_A5_QUOTIENT_COUNTS_2026_11_16.tex).
Its A5 count is a conditional knot test requiring twenty distinct quotient
kernels for ONE double-cover group. It supplies no Keller map, finite lattice
or knot certificate here. Its centred Gaussian corollary retains its supplied
weighting hypothesis. Neither deduction is promoted to a new HC4 theorem.

The retained MVALGORDSAFE20261116 theorem is a genuine potential global
period bridge: under purity and its EXACT time-zero source hypothesis,
algebraicity of ALL ORDINARY original marked resolvents implies all-marker
safety by a Weyl-linear retraction and indecomposability. It is not a statement
that algebraic source functions make their ordinary period resolvents algebraic.
The ordinary Cauchy variable and target-prime valuations must stay distinct.

Here is a small finite-order control. Put B0=C[t], S0=C[u], t=u^2,
F0=C(t), L0=C(u), and z=1/u=u/t. Then

    chi_z(T)=T^2-1/t,
    Tr(z)=0, Tr(z^2)=2/t,
    (1/2)*Tr(1/(s-z))=s*t/(t*s^2-1)
                            =1/s+1/(t*s^3)+1/(t^2*s^5)+... . (5)

The resolvent is RATIONAL over F0(s), yet z is not integral over B0 and
its moment coefficients have target-prime poles. This is a non-Keller DVR
control, not an actual plane counterexample. Algebraicity of the finite
source field, or even rationality of its field-trace resolvent, therefore
does not supply the B-regularity demanded in Section1. A field-trace moment
functional is not substituted for the exact Gaussian/Haar period source;
purity of a phase for finite field trace is itself very restrictive by
Newton identities, and is not provided by the HC4 field reconstruction.

Any proposed Mathieu-to-HC4 bridge must exhibit the ACTUAL marked period
source, its purity and time-zero hypothesis, prove ordinary marked-resolvent
algebraicity (or the appropriate zero-value logarithmic condition), and then
identify the resulting conclusion with the ACTUAL B-regular finite lattice.
None of those transfers follows merely from L/F being finite. Exponential
series, ordinary moment resolvents and target-prime regularity are different
objects. The full exceptional Mathieu theorem is not assumed here.

## 4. Remaining gap and assurance

The NEW global reduction is: TWO first source-derivative integrality tests
are necessary and sufficient for the ENTIRE additive finite-module/trace-lattice
stability condition. Under full-field reconstruction, their poles specify
exactly which whole primes the normalized directional saturation deletes
from the FULL finite normalization. No more field generators or increasing
jet-order census is needed to formulate this gate.

What is still missing is automatic passing of these tests for enough source
directions, or another proof of Tr(Lambda*Lambda) subset B/full marked-energy
integrality. Neither the finite trace-moment criterion nor the maximal-open
construction establishes that premise. Ramified primes remain present in
all statements; the etale completion is not substituted for S.

The retained MR0/KX1/LO1/KN0/ordinary-resolvent proofs are used at their exact
scope, with Zariski Main, normalization finiteness and normal polynomial rings
as already used foundations. No new external theorem is adopted. The new
checker only checks local differential/valuation, Newton/field-trace controls
and provenance; the unbounded arguments above establish the deductions.
No registered mathematical checker is rerun. Full HC4/JC2 stays OPEN,
published quartic PLANE exclusion is retained, and external review/outreach
remains DEFERRED.
