# Uniform integration and the exact relative quotient

Date: 2026-11-16. Claims HC4MYLG0 and HC4MYLG1.

FULL THEOREM TARGET: prove HC4/JC2 or construct an actual counterexample.
CURRENT LAST UNSOLVED IMPLICATION: ACTUAL polynomial plane source forces
Tr(Lambda_L*Lambda_L) subset B, equivalently regular inverse Gram entries;
the alternative full marked-energy integrality implication is also OPEN.
No such automatic regularity is established here. External review/outreach
remains DEFERRED. The previous goal turn was PROGRESS: LF0/1 ruled out a
general trace-to-complete-flow inference in every odd degree.

This note changes direction from averaging flows to integration on the
whole polynomial source. The uniform integration criterion below is proved
by the simplicity of polynomial differential ideals. The surjectivity
criterion for plane Jacobian derivations is CLASSICAL, not a new result.
The exact relative quotient retains the special-fibre torsion that is lost
on passing to a rational Gauss-Manin connection.

## 1. Setup and the classical zero-quotient criterion

Let R=C[x,y], B=C[P,Q] subset R, with j=Jac(P,Q) a nonzero constant.
Define the uniquely lifted target derivatives

    E=D_P=(Q_y partial_x-Q_x partial_y)/j,
    D=D_Q=(-P_y partial_x+P_x partial_y)/j.

They are polynomial derivations, commute, and satisfy

    E(P)=D(Q)=1, E(Q)=D(P)=0,
    d r=E(r)dP+D(r)dQ.

Their coefficient matrix is invertible over R; explicitly

    partial_x=P_x E+Q_x D,
    partial_y=P_y E+Q_y D.                             (1)

Put J=D(R), an additive C-subspace, and M=R/J. Multiplication by P
and application of E descend to M, and [E,P]=1. Multiplication by an
arbitrary source polynomial does NOT automatically descend.

We will use the classical equivalence

    M=0  <==>  (P,Q) is a polynomial automorphism.       (2)

Here is a proof adequate for the present setting. The forward implication
gives surjectivity of D on R/(P-c) for every c. Choose a smooth connected
generic fibre, using retained LE0 relative algebraic closedness. On that
affine curve dQ is a nowhere-zero generator of its regular one-forms.
Surjectivity says every such form r*dQ is exact. Algebraic de Rham theory
then gives H1=0. If its smooth projective completion has genus g and s>0
punctures, dim H1=2g+s-1, so g=0,s=1 and the fibre is A1. Its regular
function Q is a polynomial with nowhere-zero derivative, hence linear.
Thus the plane map has generic degree one. Retained PF52 gives
Frac(R)=Frac(B) ==> R=B, proving invertibility. Conversely on R=B,
D=partial_Q is surjective by polynomial integration.

Friedland (1999), Theorem5.1, recalls this criterion and credits Krasinski
and Spodzieja, with an earlier analytic version by Stein. Gurjar, Masuda
and Miyanishi (1211.0744v1), Corollary2.4 and Remark2.5, give the related
classical surjective-derivation formulation. We do NOT claim priority for
(2) or use a general classification of surjective derivations.

## 2. LG0: a uniform integration multiplier forces inversion

For a C-subspace V subset R define its ideal core and multiplier ring by

    I(V)={a in R : a*R subset V},
    A(V)={h in R : h*V subset V}.

The first is an ideal, the second a subring. If V is stable under BOTH
commuting Jacobian derivatives E,D, its ideal core is stable under them:
for a in I(V), r in R and T=E or D,

    T(a)*r=T(a*r)-a*T(r) is in V.                     (3)

Equation(1) makes I(V) stable under partial_x and partial_y. The polynomial
ring over a characteristic-zero field has no proper nonzero ideal stable
under all coordinate derivatives: differentiating a nonzero element of
least total degree produces a nonzero constant. Consequently

    I(V)=0 or R;
    if I(V)!=0, then V=R.                             (4)

This argument works in EVERY DIMENSION for a whole polynomial source
with an invertible polynomial frame of lifted target derivatives. It
does not assert that one derivative being surjective forces inversion in
higher dimensions. The curve argument in (2) is specifically PLANE.

Apply (4) to J=D(R). Commutation gives E(J) subset J and D(J) subset J.
Combining with (2), in EVERY generic and polynomial degree:

    exists 0!=a in R with a*R subset D(R)
        <==> D(R)=R <==> the plane map is invertible.   (5)

This is a UNIFORM statement over EVERY source polynomial r. The primitive
of a*r may depend on r; no linear choice of primitives is assumed. It
applies to a anywhere in R, not just in B or C[P]. Under plane failure
the exact uniform integration ideal is ZERO, not merely a proper ideal.

The multiplier ring is also exact. If h*J subset J, then for all r,

    D(h)*r=D(h*r)-h*D(r) is in J.                      (6)

Thus D(h) belongs to I(J). Under noninvertibility (5) forces D(h)=0.
Conversely a D-constant preserves J by the product rule. Retained LE0
gives ker(D|R)=C[P]. Therefore under plane failure

    A(D(R))=C[P].                                    (7)

Equivalently, a single multiplier h not in C[P] preserving the ENTIRE
derivative image forces inversion. This does not mean checking h on a
finite list of derivatives suffices. In the inverse case J=R and A(J)=R.
The same conclusions hold in every polynomial source-coordinate frame
and every polynomial target-coordinate frame. No degree cutoff is used.

For example, Q*D(R) subset D(R) alone would imply D(R)=R directly from
D(Qr)=r+Q*D(r). It is precisely this image stability that is NOT automatic.
The identity 1=D(Q), exactness of the full de Rham complex and separate
local primitives do not establish it.

## 3. LG1: the unlocalized quotient has a bijective transverse derivative

On the ACTUAL unlocalized M=R/D(R), E is BIJECTIVE, whether or not M=0.
Surjectivity follows from H2_dR(A2)=0: for any r choose a polynomial
one-form alpha=a*dP+b*dQ with

    d alpha=r*dP wedge dQ,
    r=E(b)-D(a).

Then [r]=E[b]. For injectivity, E[r]=0 means E(r)=D(s) for some s.
The form s*dP+r*dQ is closed. H1_dR(A2)=0 makes it d h for a polynomial h;
comparison of coefficients gives r=D(h), hence [r]=0. In particular

    ker(E|M)=coker(E|M)=0.                            (8)

This is stronger than merely saying J+E(R)=R. It does NOT say M=0.

Let T be the C[P]-torsion submodule of M, and V=M/T. These are E-stable.
For c in C put z=P-c and let r_c be the number of irreducible components
of the fibre P=c. Smoothness makes them reduced and disjoint. On each
component D is nonzero and its regular constants are C. Applying the
kernel/cokernel sequence to multiplication by z on (R,D) gives

    ker(z|M) = C^(r_c)/C,
    dim ker(z|M)=r_c-1.                              (9)

Here ker(D|R)=C[P] maps diagonally to the constants on the fibre.
No claim that EVERY fibre is irreducible has been introduced.

Write delta_c for the point Weyl module with basis u_l=E^l u_0 (l>=0),
z u_0=0 and z u_l=-l*u_(l-1). The locally z-nilpotent part of any Weyl
module is a sum of copies of delta_c. To see this directly, the finite sum

    Pi_c(w)=sum_(l>=0) E^l*z^l*w/l!

lies in ker z, and the inverse expansion is

    w=sum_(l>=0) (-1)^l*E^l*Pi_c(z^l w)/l!.

Commutation [E,z]=1 verifies both identities and independence. Thus

    T = direct_sum_c delta_c^(r_c-1).                 (10)

Only finitely many r_c exceed one: generic geometric irreducibility and
the bounded polynomial degrees of possible fibre factors make the
reducible-value set constructible and nongeneric in A1, hence finite.
This is the classical total-reducibility setting, not an assertion of
irreducibility for all special fibres.

On delta_c, E is injective with one-dimensional cokernel. The exact
sequence 0->T->M->V->0 and (8) now give

    E:V->V is SURJECTIVE,
    ker(E|V) = coker(E|T) = C^delta,
    delta=sum_c(r_c-1).                              (11)

The connecting isomorphism is canonical; the displayed coordinate-space
identification depends on bases. If delta>0 the extension cannot split
as a Weyl module, since a splitting would make the nonzero coker(E|T)
a summand of coker(E|M)=0. M has no nonzero point-supported quotient,
for the same surjectivity/cokernel reason.

Localizing in C[P] kills T and yields the rational relative de Rham
connection. V is the image of M in that localization. Friedland's
unlocalized period image on p245 is this saturated quotient, rather than
the raw M. Its transverse derivative is surjective; rational horizontal
classes can survive. Our explicit raw/saturated calculation explains why
discarding special-fibre torsion changes injectivity. The generic
connection, the saturated original-form image, and raw M are distinct.

## 4. Exact nonzero controls and failed shortcuts

All properties in (8) allow NONZERO modules. For example take

    K_alpha=C[z,z^-1]*e, alpha rational and noninteger,
    P(z^m e)=z^(m+1)e,
    E(z^m e)=(m+alpha)*z^(m-1)e.

Then [E,P]=1 and E is bijective with an explicit inverse on every basis
element, but K_alpha is nonzero. This is a regular-singular algebraic
Kummer connection; acyclicity alone does not give vanishing.

There is also a control WITH the point torsion in (10). Set

    M0=C[t,t^-1], E(w)=t*w, P(w)=-d w/dt.

Again [E,P]=1 and E is bijective. Its P-torsion is exactly T0=C[t],
one copy of delta_0. The quotient V0=M0/T0 has E surjective and
ker E=C*[t^-1]. Its extension is nonsplit by the cokernel argument above.
These are abstract module controls, NOT polynomial plane Keller sources.
Neither their nonzero existence nor a formal completion proves anything
against the actual whole-source HC4 implication.

Dropping the full commuting polynomial frame also breaks (4): for
D=x*partial_x on C[x,y], D(R)=xR has the nonzero proper ideal core xR.
The slice, Hamiltonian, constant-Jacobian and whole-source hypotheses must
not be silently omitted from the inversion criterion.

## 5. Literature, assurance and the remaining gap

Primary source inspected: Shmuel Friedland, *Monodromy, differential
equations and the Jacobian conjecture*, Ann. Polon. Math.72(1999),219-249,
DOI10.4064/ap-72-3-219-249. Actual complete pages230 and242-245 were
rendered and visually read, together with extracted context225-232.
Theorem5.1, Proposition5.3, Corollary5.6 and Proposition5.7 with proofs
are the classical baseline. Proposition2.9 explicitly warns that affine
fibre monodromy need NOT be semisimple, in contrast to the compact part.
We do not import semisimplicity, Kaliman's change-of-target reduction, or
an unrestricted fibre-invariance theorem from these remarks. The present
identities, ideal-core argument and point-module decomposition are written
out above. The classical surjectivity criterion is credited and reproved,
not registered as a newly discovered theorem. No new major external
classification theorem is imported.

The checker gives exact transcription controls for polynomial frames,
de Rham homotopies, point-module projectors and the two nonzero controls.
Finite samples are NOT a proof of the universal identities or the
absence of an actual Keller counterexample. No registered mathematical
checker is rerun. Independent replay, formal verification, external review
and original goal completion remain false.

The useful new sufficient test is a SINGLE NONZERO UNIFORM source
integration multiplier, or a single nonconstant-along-D multiplier of
the WHOLE image D(R). Neither is supplied automatically by polynomiality.
The original actual-source trace-pairing/full marked-energy implication
remains OPEN. This is progress on the full goal, not a reduced definition
of success or an HC4/JC2 proof.
