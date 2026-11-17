# The Long purity derivation cannot be transported by a regular volume-preserving lift

16 November 2026. Internal written proof. Independent replay, formal
verification and external review are not claimed; external review remains
deferred. Both FULL exceptional Mathieu--Zhao statements remain **OPEN**.
The family-transfer theorem is already recorded and is unchanged.

This note tests a direct all-order counterexample mechanism. A regular
volume-preserving derivation with eigenfunction f and nonzero constant
eigenvalue proves every positive pure moment zero. Long's known phase has
such a derivation. However, on BOTH exceptional quotients, that derivation
has NO regular volume-preserving lift through ANY of the retained regular
Bezout maps, in ANY degree or parameter. This remains true after ANY
algebraic automorphism of the source. The obstruction is a calculation at
the unequal-exponent cusp component of the Klein branch divisor.

This excludes that specific transport of a purity certificate. It does
NOT prove that every Bezout pullback is nonpure, exclude other derivations
with the same eigenfunction, split the intrinsic boundary extension, or
construct an exceptional pure counterexample.

## 1. A direct sufficient mechanism for all positive pure moments

Let X be a smooth affine complex threefold with a regular nowhere-zero
volume form Omega, and let h be its normalized period on a compact
three-cycle. For a regular algebraic derivation E define divergence by

    L_E Omega = (div_Omega E) Omega.

Cartan's formula and Stokes give, for every regular A,

    h(E A) + h(A div_Omega E) = 0.                       (1)

Consequently, if E is regular, div_Omega E=0, and E f=lambda f for a
NONZERO CONSTANT lambda, then

    0=h(E(f^m))=m lambda h(f^m), m>=1.                 (2)

Thus f is all-positive pure. Completeness of E is NOT required.
Completeness must not be assumed: this is a possible mechanism beyond
contracting algebraic torus actions, whose safety is already known.
A variable eigenvalue a gives only h(a f^m)=0 and does not give (2).

For the original exceptional Bezout phase g_b below, finding ANY such
original derivation E with E g_b=g_b would produce an ACTUAL pure failure:
the fixed original marker C_b is already supplied by the retained
[uniform-marker theorem](COHOMOLOGICAL_PULLBACK_AND_A_UNIFORM_ORIGINAL_MARKER.md).
Equation (2) would settle purity at every order, not by testing a prefix.
No such original derivation is constructed here.

## 2. Long's target derivation

Use target coordinates [[X,Y],[F,G]], XG-YF=1, put t=1+Y, and set

    g=t+F t^2=(1+Y)(1+F+YF).

The following derivation is globally regular on the target SL2:

    delta X=0, delta Y=-Y(1+Y),
    delta F=1+F+2YF, delta G=YG.                        (3)

It kills XG-YF-1. On the chart X!=0 take

    Omega_*=dX/X wedge dY wedge dF.

Its divergence is

    partial_Y[-Y(1+Y)]+partial_F[1+F+2YF]=0.           (4)

Since (3) and the volume are regular globally, (4) holds globally.
A direct identity gives delta g=g. Equations (1)-(2) therefore certify
all positive pure moments of g. Another first integral is

    B=(1+Y)G, delta B=0.                               (5)

More precisely delta is the Hamiltonian derivation defined by

    dX wedge dB wedge dR = (delta R) Omega_*.

This last formula can be checked on the same dense chart. It explains
why simply pulling back its two Hamiltonians produces a Jacobian factor.

## 3. The entire original Bezout family and its branch component

Use the full original algebra A_Gamma=C[SL2]^Gamma, Gamma=2O or2I,
and the [retained Bezout construction](REGULAR_BEZOUT_MAPS_AND_WEIGHTED_HAAR_TAILS.md):

    Phi_b=[[X_b,Y_b],[F,G]],
    X_b=U+F b, Y_b=V+G b, UG-VF=1,
    b in A_Gamma arbitrary regular.

Let Omega_Gamma be the original nowhere-zero volume. For a nonzero
constant c depending only on the volume normalization, its Jacobian is

    Phi_b^* Omega_* = J_b Omega_Gamma,
    J_b=c W_b, W_b=C(H+D b), D=z partial_x+w partial_y. (6)

The retained weighted identity h(W_b)=1 proves J_b!=0 for EVERY b.
In particular Phi_b is dominant and its function-field extension is
finite separable. The Klein relation is

    C^2=K(F,G),
    K_O=G(F^3-108G^2)/36,
    K_I=12F^5+G^3/144.                                (7)

The original quotient maps surjectively to the punctured Klein surface
as an affine-line torsor. Select the following cusp component of C=0:

    Gamma   G^q=lambda F^p       (p,q)      r=q/p
    2O      G^2=F^3/108         (3,2)       2/3
    2I      G^3=-1728F^5        (5,3)       3/5.         (8)

Each punctured cusp is an irreducible curve. Its inverse image Delta in
the original quotient is an irreducible divisor: an affine-line torsor
over that curve. On Delta, F and G are nonzero and F is NONCONSTANT.
The other 2O component G=0 is not used or silently discarded.
Equation (6) implies ord_Delta(J_b)>0, independently of how much extra
vanishing H+D b has. No minimum ramification-order assumption is needed.

## 4. A volume-preserving lift would preserve this divisor

Suppose a regular original derivation E lifts (3), meaning

    E(Phi_b^*R)=Phi_b^*(delta R) for EVERY target R,     (9)

and suppose div_Omega_Gamma E=0. Lie differentiation of (6), using
(4) and (9), yields

    E(J_b)=0.                                        (10)

An elementary valuation fact now matters. At the generic point of any
prime component of the zero divisor of J_b, write J_b=pi^e u in its
DVR, with e>=1 and u a unit. Since E is regular it extends to this
localization. Equation (10) gives

    0=e u E(pi)+pi E(u), hence E(pi) in (pi).

Thus E preserves the prime ideal of Delta and induces a derivation of
its function field. This uses BOTH regularity AND zero divergence.
A rational lift, or a regular lift with nonzero divergence, does not
satisfy the same conclusion.

## 5. The contradiction on the cusp, for both groups simultaneously

Restrict (9) to Delta. Write X,Y for the restrictions of X_b,Y_b.
The cusp equation in (8) gives

    q E(G)/G=p E(F)/F.

Using (3) and writing k=2-r, this becomes

    r Y=(1+F)/F+2Y,
    Y=-(1+F)/(k F),
    E(F)=-r(1+F)/k, E(G)=-(1+F)G/(k F).              (11)

The determinant then forces

    X=(1+YF)/G=(k-1-F)/(k G).                         (12)

But differentiation of (12) by (11) gives the exact identity

    E(X)=(r-1)(1+F)(F-1)/(k^2 F G).                  (13)

Both r=2/3 and r=3/5 differ from 1 and 2. Since F is nonconstant on
Delta, the right side is NONZERO in its function field. The coefficients (r-1)/k^2 are -3/16 on2O and -10/49 on2I.
This contradicts E(X)=Phi_b^*(delta X)=0. It proves the promised universal obstruction:

**For EVERY regular b on either exceptional quotient, the Long
Hamiltonian derivation (3) has NO regular divergence-free lift through
Phi_b.** No phase-degree or weight-support bound occurs in this proof.

The calculation also gives an abstract unequal-exponent-cusp obstruction:
whenever a positive-Jacobian-order divisor maps onto G^q=lambda F^p
with p!=q, p!=q/2 and nonconstant F, the same hypotheses and equations
force (13). We use only the two explicit components (8) here; no claim
about arbitrary finite-group quotients is needed.

## 6. Why division by the Jacobian does not transport purity

Pull back the Hamiltonians X,B from Section2, and define the ORIGINAL
Hamiltonian derivation V_b by

    dX_b wedge d(Phi_b^*B) wedge dA=(V_b A)Omega_Gamma.

It is regular and divergence-free. Equations (3),(6) imply

    V_b(Phi_b^*R)=J_b Phi_b^*(delta R),
    V_b(g_b)=J_b g_b.                                 (14)

Hence Stokes supplies h(J_b g_b^m)=0, exactly a WEIGHTED identity.
The rational derivation E_b=V_b/J_b is the unique function-field lift
of delta: separability gives uniqueness of derivation extension.
Its divergence is

    div E_b=-E_b(J_b)/J_b.                            (15)

The theorem says that E_b cannot be BOTH regular and divergence-free.
Treating E_b as an original volume-preserving derivation would therefore
lose a genuine obstruction at the branch divisor or in the divergence.
One cannot infer the unweighted pure identities by this division.
This proof allows arbitrary cancellations in b and arbitrary extra
Jacobian vanishing; neither changes (13).

## 7. Automorphisms, the fixed marker and the full obligation

The obstruction survives precomposition Phi_b o phi by ANY algebraic
source automorphism. Such automorphisms multiply the original volume by
a nonzero constant, by the retained
[intrinsic Haar theorem](INTRINSIC_REYNOLDS_AND_FINITE_QUOTIENT_DESCENT.md).
Pushing a hypothetical regular divergence-free lift through phi would
give the forbidden lift through Phi_b. Nonzero constant multiples of
the target derivation are excluded as well, by dividing that constant.

The retained functions and identities remain in force:

    g_b=(1+Y_b)(1+F+Y_b F), C_b=W_b(1+Y_b),
    h(C_b R(g_b))=integral_0^1 R(-u)du                 (16)

for EVERY polynomial R, with ONE fixed original regular marker C_b.
Thus any eventual-pure H(g_b) would still be a genuine failure, with
that SAME marker. The present no-lift theorem neither constructs such
an H nor proves that none exists. A DIFFERENT original derivation E
with E(g_b)=g_b need not intertwine all four target coordinates (9),
and is not excluded. The sufficient eigenfunction mechanism (2) is
not claimed to describe every pure phase.

The whole remaining target is unchanged: in ALL degrees, every pure
original phase must lie in the safe locus, equivalently its intrinsic
boundary extension must split, or one actual pure failure with a fixed
original bad marker must be constructed. The sharp capacity inputs
(failure support at least six weights on2O/eight on2I after every allowed
automorphism and torus), primitive/Euler restrictions and all earlier
results are retained. No finite-degree campaign, additional decidability
theorem, arbitrary-ring promotion or external-review gate is introduced.
