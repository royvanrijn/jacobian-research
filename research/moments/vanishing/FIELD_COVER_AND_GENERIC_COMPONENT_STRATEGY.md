# Field-valued safety: radical covers and generic component certificates

Preferred explicit global gate: the
[nilpotent contraction-zero theorem](NILPOTENT_CONTRACTION_ZERO_AND_THREE_GLOBAL_TARGETS.md)
classifies every pure phase with p=r=0 and Delta^2(H2*H3)=0, with
unrestricted linear harmonic, as one-sided or Rodrigues. This reduces a
sufficient complete-chart safe-cover certificate to three radical tests
of degrees(2,3,2) and positive weights(6,9,5). Nine fixed rational
rotations recover the entire contraction vector from one normalized
coefficient. All three moment memberships remain open. Nonzero B needs
only twelve moments under the contraction condition; B0 uses the
all-order odd theorem, so odd prefix false positives need not be removed
by this gate. The earlier31-generator ideal remains a stronger optional
Rodrigues/nullcone prefix-exhaustion target.

New global reduction: the [complete positively graded highest-cubic chart](HIGHEST_CUBIC_GLOBALIZATION_AND_GRADED_GATE.md)
now covers all H3-nonzero phases by seven fixed rotations and two explicit
normalization parameters, including every D0 component. Prefix14 makes
the normalized chart finite over five base variables. Both the first14
zero locus and the full pure locus have dimension7. Every cubic pure
component contains z^3, and the normalized reduced moment algebra embeds
in the completed germ there. Marked-tail identities in that complete
reduced germ would therefore give a global finish. This remains a
component/marked-safety obligation; neither the dimension theorem nor
the finite central fibre establishes the required safe cover.

New complete input option: the [five quadratic normal forms](QUADRATIC_ORBIT_NORMAL_FORMS.md)
cover the whole coefficient space over algebraic field extensions, with
all linear/cubic coefficients retained and exact D factors. They allow
a global five-family campaign without first localizing at D. These are
coefficient containers, not actual pure components, so the purity and
generic-safety obligations below remain in force. The odd family is
closed; the other four full mixed-linear families remain open. The
[axial exclusion](AXIAL_QUADRATIC_ZERO_LINEAR_EXCLUSION.md) and nilpotent
zero-linear theorems supply further certified subcases.

17 November 2026. The user proposed allowing larger moment prefixes,
arbitrary radical exponents and safe components beyond Rodrigues, and
pointed to the already retained generic-to-special theorem. These are
valid alternatives to the current fixed-prefix length computation.
This note records the precise obligations and updates the work priority;
it does not prove a new global exhaustion or safety statement.

The existing [finite-family decision theorem](ALGEBRAIC_SPECIALIZATION_AND_FINITE_FAMILY_DECISION.md),
Section5 step4, already proves that generic all-original-multiplier safety
descends to every point of an irreducible pure component. Its general
component/holonomic algorithm is written, not implemented or executed
for the dense(3,3) family. The simpler polynomial-identity argument below
is all that is needed to justify this particular descent.

## 1. The primary field-valued certificate

Let \(R\) be the relevant affine chart algebra, including all localization
relations, and put \(I_\infty=(\mu_1,\mu_2,\ldots)\). The full pure point
locus is \(X_\infty=V(I_\infty)\). Suppose \(Y\) is a closed set of known
safe pure points and \(J_{\rm safe}\) defines it. It suffices to prove

\[
 J_{\rm safe}\subseteq\sqrt{I_N}
\]

for any finite workable \(N\). Indeed
\(X_\infty\subseteq V(I_N)\subseteq V(J_{\rm safe})\).
There is no requirement that this prefix generate the full ordinary
pure ideal, or even have the same radical as that ideal.

More generally \(Y\) can be a closed safe container: it is enough that
every **pure** point in \(Y\) is already known safe. This distinction is
needed for a slice theorem that proves purity-implies-safety on a slice;
the entire coefficient slice need not itself be pure.

For a finite union \(Y=\bigcup_jY_j\), use
\(J_{\rm safe}=\bigcap_j I(Y_j)\). The sum of component ideals would
instead define their intersection. Equations for a chosen safe cover
and complete coverage are still required; naming some safe families
does not supply them automatically.

For each generator \(f\) of \(J_{\rm safe}\), the exact unit test is

\[
 1\in I_N+(1-vf)\quad\hbox{in }R[v].                 \tag{1}
\]

This is equivalent to \(f\in\sqrt{I_N}\): the quotient is the localization
\((R/I_N)[1/f]\), which is zero exactly when a power of \(f\) is zero in
\(R/I_N\). See [Stacks localization](https://stacks.math.columbia.edu/tag/00CM).
Given \(f^e=\sum_j A_jg_j\) with \(g_j\in I_N\), a direct unit identity is

\[
 1=\sum_jv^eA_jg_j+(1-vf)\sum_{k=0}^{e-1}(vf)^k.
\]

Conversely substitute \(v=1/f\) into a unit identity and clear the
finitely many powers of \(f\). The exponent need not be prescribed.
The tests for different generators are separate: imposing every
\(1-v_if_i\) simultaneously would test a different open set and is
insufficient for this cover. If separate identities use prefixes
\(N_i\), their finite maximum works for all of them.

For an unnormalized \(D\ne0\) chart the inverse relation must be retained
explicitly; on the current normalized chart \(D=1/8\) it is already built
into the presentation. No such normalization covers the \(D=0\) boundary.

## 2. A concrete radical fork for the existing seven targets

The [selected-locus support theorem](SEVEN_LOW_DEGREE_GLOBAL_GATE.md) gives
\(P=I_9+(K_1,\ldots,K_7)\), where \(P\) defines the three reduced
normalized Rodrigues points. Hence for every \(N\ge9\),

\[
 P\subseteq\sqrt{I_N}
 \quad\Longleftrightarrow\quad
 K_i\in\sqrt{I_N}\ (i=1,\ldots,7).
\]

The radical fork can therefore use the \(K_i\) themselves, of degrees
\((1,1,1,1,1,1,2)\), in seven separate tests (1). There is no need to
insert their local nilpotence powers or the cubic replacement first.
The initial available moment input is through14; later exact moments
may be added as needed. No global stopping bound is known.

The [Rabinowitsch input exporter](export_cyclic_field_radical_gate.py)
retains the full source ideal, including \(\mu_{13},\mu_{14}\), and
exports [seven graded unit-test inputs](cyclic_field_radical_gate_program.json).
For \(K_i\) of character \(w_i\), assign \(v\) character \(-w_i\pmod{12}\).
Then \(1-vK_i\) has character0, so the exact rational grading survives.
The multiplier on a source column of character \(c\) can be restricted
to monomials of character \(-c\), with \(v\)'s assigned character included.
Finite witness degrees and actual finite matrix column counts remain unknown.

There is a useful limitation to the word "weaker" here. At any **fixed**
\(N\ge14\), the local certificate already gives length17 at each of the
three Rodrigues points. Thus

\[
 \begin{aligned}
 &K_i\in\sqrt{I_N}\text{ for every }i\\
 &\quad\Longleftrightarrow\quad V(I_N)=\{\text{three Rodrigues points}\}\\
 &\quad\Longleftrightarrow\quad\dim_{\mathbb C}A_N=51.
 \end{aligned}                                                    \tag{2}
\]

Finite support forces an Artinian algebra, and the known local lengths
then determine its total length. Arbitrary radical exponents may still
make the computation easier, but at the same prefix this is an
alternative certificate for the same Rodrigues-only conclusion. The
main logical flexibility comes from allowing a later prefix, or
allowing additional safe components in the cover. Neither Rodrigues-only
exhaustion nor length51 is necessary for the full Image theorem.

## 3. Generic safety already covers every specialization

For this descent, use the original phase coefficient family, or a frame
presentation carrying it and all original marked outputs. Fixed original
markers need not be scalar invariant functions on the unmarked quotient.
The identities to descend are their regular polynomial output coefficients
in the original family.

Let \(C\) be an irreducible pure component with prime ideal
\(\mathfrak p\), domain \(A=R/\mathfrak p\), and fraction field \(K\).
For an original monomial multiplier \(Q\), write each output coefficient
of \(\mathcal E(P^mQ)\) as \(c_{m,Q,\gamma}\in A\).

Generic all-marker safety says that for this fixed \(Q\) there is a
threshold \(M_Q\) with

\[
 c_{m,Q,\gamma}=0\quad\hbox{in }K
 \quad(m\ge M_Q,\ \text{all output indices }\gamma).
\]

The map \(A\hookrightarrow K\) is injective, so these are zero elements
of \(A\), not just identities on a generic open subset. They remain
zero at every point of \(C\). Denominators in a generic annihilator,
Euler field or other certificate do not create exceptional points for
these polynomial tail identities. No direct specialization of such a
rational certificate is required.

Every original polynomial multiplier is a finite linear combination of
original monomials; take the maximum of their finitely many thresholds.
This proves all-marker safety at every point of \(C\), with no need for
one threshold valid for all multiplier degrees. The output variables
remain independent; scalar evaluation of the output at a selected point
must not replace the full polynomial-valued contraction.

The premise matters. A component of a finite-prefix zero set need not
be pure. For a generic-component classification route, use the actual
pure locus, or certify that the chosen finite prefix has the correct
**radical**. Ordinary scheme equality is unnecessary for this field-valued
step. If a finite-prefix component is generically nonpure, a later
moment cuts it to a proper closed subset; its special pure points must
still be retained, rather than discarding that whole component.

One can organize a finite cover certificate incrementally: split a
finite-prefix locus into components, certify any actual pure components
generically safe, and cut the remaining generically nonpure pieces by
explicit later moments. Each witnessed proper cut lowers dimension.
This is a valid proof format, not a computed component decomposition
or an implemented decision campaign.

## 4. Boundary and nilpotent obligations

Generic safety of a full pure component meeting \(D\ne0\) covers its
\(D=0\) specializations automatically by the same polynomial identities.
The remaining boundary task includes every pure component contained in
\(D=0\) and any boundary points not covered by the already certified
closures. Keep the retained cyclic-rank/Gram-degeneracy strata and every
cubic coupling. Finishing the normalized cyclic chart alone does not
complete the full(3,3) problem.

Field generic points detect the reduced locus. A nilpotent coefficient
can vanish in every field while remaining nonzero in the ordinary
coefficient algebra. Thus this reorientation proves field-valued and
reduced-family conclusions; arbitrary coefficient-ring extensions still
require the appropriate ordinary identities or nilpotent analysis.

The known product of three length17 local algebras is a quotient of
\(A_{14}\), giving a lower bound. Checking moments in that product alone
does not give an upper bound, global support coverage or exhaustion.
Its exact local information is valuable once coverage is proved, as in
(2), but it is not itself a coverage certificate.

## Current work status

- Primary finish line: field-valued all-original-multiplier safety, by
  a finite radical safe cover and/or generic safety on every actual pure
  component. Additional pure components are acceptable when safe.
- Cyclic Rod-only option: seven separate Rabinowitsch units, allowing
  later prefixes and unrestricted radical exponents.
- Stronger benchmark: the existing seven fixed-power/cubic memberships
  and exact \(A_{14}\) length51.
- Independent residual obligation: \(D=0\) components and uncovered
  boundary strata.

No radical unit identity, additional generic safety certificate, full
component decomposition, global finiteness result or full Image theorem
is established by this note. The input exporter records zero verified
unit certificates and no witness-degree bound.
