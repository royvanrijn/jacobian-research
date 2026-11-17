# Ramification and admission of finite quantum closure carriers

15 November 2026. **QGC1.** Long larger S3 closure has a quadratic
ramification divisor inside its original affine source. Its natural
division action cannot pass through a full finite carrier with the
correct classical graph at ordinary quantum parameter, or through an
equivariant carrier of the correct torsion-free classical rank at any
ramified quantum parameter. Saturation need not be assumed in either
case: a finite valuation level or finite character shifts recover the
residue multiplier obstruction.

The core theorem says that a full finite residue lattice preserved by
the induced scalar-potential Hamiltonian connections has a finite
etale multiplier ring over the physical target polynomial ring.
Ramification forbids it. The arguments cover arbitrary fractional
classical lattices and quantum correction degrees without a free frame
or quantum centralizer finite generation.

Saturation itself is not automatic, even for full Weyl actions with
the correct identity graph. The positive constant-cover controls show
its failure and the finite equivariant repair. Arbitrary matrix repairs,
nonequivariant ramified actions lacking the quantum-order graph
congruence, and special modules with torsion or excess rank retain
separate admission obligations. The ordinary rank-one Long graph is
not excluded.

## The finite carrier and its classical lattice

Let K have characteristic zero, R=K[[t]], and let
hbar=a(t)=t^e*u(t), with e>=1 and u a unit in R.
Let A_R be the finite-PBW Weyl algebra on n pairs, with
[p_i,x_j]=hbar*delta_ij. Its generic division ring D has
the t Gauss valuation and residue field
k=K(x_1,...,x_n,p_1,...,p_n)=Frac S, where S is the polynomial ring.

Let Dtilde be a valued division extension containing D, with the
same central uniformizer t and commutative residue field ell,
finite separable over k. Let Otilde be its valuation ring.
Suppose a left A_R-submodule M subset Otilde satisfies
\[
 \begin{gathered}
 M\text{ is finitely generated over }A_R,\qquad
                    M\cap tO_{\widetilde D}=tM,\\
 J=\operatorname{res}(M)\subset\ell,\qquad kJ=\ell .
 \end{gathered}                                                \tag{1}
\]
Thus M/tM embeds as the full finite S-lattice J in the residue
field. No choice of basis, ideal presentation or projective
structure is imposed. Saturation at all t powers follows
inductively from the displayed saturation at t.

Let F_1,...,F_(2n) lie in O_D and preserve M on the right.
Their residues psi_i in S are a full polynomial symplectic frame:
\[
            \{\psi_i,\psi_j\}=J_{ij}^{\rm symp},\qquad MF_i\subset M.
                                                               \tag{2}
\]
Their Hamiltonian fields span Der_K(S) over S. The necessary
condition below does not require quantum commutator relations
between the F_i; it therefore applies in particular to a full
quantum Weyl tuple with these classical limits.

## The multiplier ring is finite with the correct fraction field

Define the commutative multiplier ring
\[
                  B=(J:J)=\{b\in\ell:bJ\subset J\}.             \tag{3}
\]
It contains S and injects into End_S(J). The latter is a finite
S-module, and S is noetherian, so B is finite over S.

Moreover Frac B=ell. For any theta in ell, express theta times
each of finitely many S-generators of J in a k-basis selected
from J. Clearing finitely many S-denominators gives nonzero
s in S with s*theta*J subset J. Hence s*theta belongs to B.
Finally
\[
                         B\cap k=S,                           \tag{4}
\]
because B is integral over the normal polynomial ring S.
These facts hold for every fractional lattice J, without
assuming that J is a ring or contains 1.

## Ramifying the quantum parameter creates no escape

Commutators of integral elements of Dtilde are divisible by t^e.
Here is the reason, rather than an assumption about a chosen
presentation of the closure. Commutativity of ell first gives
divisibility by t. If divisibility by t^j is known for all
integral elements and j<e, the leading t^-j commutator gives
a biderivation on ell. It restricts to zero on k, because
the physical Weyl division ring has commutators divisible by t^e.
Every derivation of a finite separable extension is uniquely
determined by its restriction to the base. First fix one
argument in k and then let it vary in ell; this makes the
entire biderivation zero. Induction gives the assertion.

Consequently the leading normalized commutator
\[
            \{\bar b,\bar c\}_{\ell}
                    =\operatorname{res}([b,c]/a(t))            \tag{5}
\]
is a Poisson bracket extending the standard bracket on k.
It is the unique such extension to ell.

Choose a finite polynomial PBW lift P_i in A_R of psi_i.
Initially F_i-P_i is divisible by t. If, for some j<e,
v(F_i-P_i)>=j, set c_ij=res(t^-j*(F_i-P_i)) in k.
For m in M with residue b, both P_i*m and m*F_i belong to M,
and their difference has valuation at least j. Saturation
puts that difference in t^j M. Its normalized residue is
\[
                    -c_{ij}b\in J\qquad(j<e),                 \tag{6}
\]
since its commutator part has valuation at least e.
Thus c_ij belongs to B intersect k=S. Add t^j times a polynomial
PBW lift of c_ij to P_i. This increases the difference valuation.
After finitely many steps, P_i still belongs to A_R, has
residue psi_i, and v(F_i-P_i)>=e.

At quantum order e, division by a(t) is division by t^e and
a scalar unit, so saturation now defines an operator on J:
\[
 \begin{split}
 \nabla_i(b)
   &=\operatorname{res}\bigl((P_i m-mF_i)/a(t)\bigr)\\
   &=\{\psi_i,b\}_{\ell}-c_i b\in J,\qquad
       c_i=\operatorname{res}((F_i-P_i)/a(t))\in k .
 \end{split}                                                  \tag{7}
\]
Changing m by t times an element of M changes the numerator
by valuation at least e+1, so this is well-defined.
The c_i may be arbitrary rational functions; their degrees
and pole orders have not been bounded.

## Scalar connection terms cancel on the multiplier ring

For b in B and z in J, use (7) and multiplication in ell:
\[
 \nabla_i(bz)-b\nabla_i(z)
                          =\{\psi_i,b\}_{\ell}z\in J.          \tag{8}
\]
The scalar c_i cancels. When e=1, no sub-quantum descent was
needed: the F_i may even lie in Otilde rather than O_D.
Then c_i lies in ell and still cancels in exactly the same way.
The arbitrary-parameter-order statement retains F_i in D.
Hence every Hamiltonian X_i={psi_i,-} preserves B. Since the X_i span the coordinate
derivations over S, all coordinate partial derivatives extend
as derivations preserving B.

For completeness, a finite module over a characteristic-zero
polynomial ring with connections for all coordinate partials
has derivative-stable Fitting ideals. Over dual numbers,
1+epsilon*nabla is an invertible semilinear map for the
base automorphism f -> f+epsilon*partial(f).
Fitting base change proves derivative stability.
A nonzero derivative-stable ideal of S is the unit ideal,
by repeated differentiation of a polynomial of smallest
degree. The generic-rank Fitting criteria therefore make B
locally free over S. This is the same finite-connection
argument retained in
[QFC1](FINITE_CENTRALIZERS_AND_COMPLETION.md);
see [Stacks, Lemma15.8.4](https://stacks.math.columbia.edu/tag/07ZA)
for the Fitting facts.

The finite S-module Omega_(B/S) has induced connections
nabla_i(b*dc)=partial_i(b)*dc+b*d(partial_i(c)).
It is generically zero, because ell/k is separable.
The same Fitting argument makes it zero. Therefore
\[
                     B/S\text{ is finite etale}.              \tag{9}
\]
Finite local freeness and vanishing relative differentials
give etaleness by
[Stacks, Lemma29.37.15](https://stacks.math.columbia.edu/tag/02GU).

A finite etale algebra over the regular ring S is normal, so B
is the normalization of S in ell. In particular ell/k cannot
be ramified along a divisor of S.
Using the standard characteristic-zero affine-space covering
theorem, B=S tensor_K K' for a finite separable constant
field extension K'/K. This last global conclusion uses
[Milne, Remark3.3(c), Theorem3.4 and the fundamental group sequence](https://www.jmilne.org/math/CourseNotes/LECc.pdf#page=27),
as in QFC1. Over algebraically closed K it gives ell=k.
The local ramification exclusion does not require that
external global covering theorem.

## The actual larger closure of Long is ramified inside its source

Use QDI1's full birational coordinates r,v,a,d and QPA1's
classical cubic:
\[
 \begin{gathered}
 s=2v+4a-3ra^2,\qquad t_0=(ra^3-a^2-va)/2,\\
 f(X)=rX^3-2X^2+sX+4t_0,\qquad f(a)=0,\\
 Q=(3ra-2)^2-8rv,\qquad
             \Delta(s,t_0)=4v^2Q .
 \end{gathered}                                                \tag{10}
\]
Here t_0 is a classical output, separate from the quantum
parameter t. The cubic source has degree three over the
output field; its closure has group S3 by QGD1.

The other two roots and the closure field are
\[
 b={2-ra+w\over2r},\qquad c={2-ra-w\over2r},\qquad
        w^2=Q,\qquad
 \ell=K(r,a,v,d)(w)=K(r,a,w,d).                                \tag{11}
\]
The last equality uses v=((3ra-2)^2-w^2)/(8r).
Thus the larger closure is rational as a field, with degree
two over the original source field. Rationality alone does
not provide a regular Poisson order or a finite carrier.

For an explicit subgroup diagram put A=3ra-2. The reflection
fixing the chosen root is (A,w)->(A,-w). A cycle is
\[
                    (A,w)\longmapsto
                  ((-A+3w)/2,\;-(A+w)/2).                    \tag{12}
\]
They generate S3 and fix r,s,t_0,d. The reflection-fixed
field is the original source field.

The divisor Q=0 is a genuine divisor in Long's affine source,
not only an omitted boundary of the cubic chart.
Indeed, after localizing the original source at x,
QDI1's inverse formulas give
\[
 S[x^{-1}]=K[r,v,v^{-1},a,d],\qquad v=x^{-1}.                   \tag{13}
\]
For example q=(2v-rv^2)/3, Y=a-v and
beta=5v^2-rv^3-3av, and the inverse system for the original
momenta has determinant minus one. All inverse entries are
Laurent polynomials in v, with no further position divisor.
After also localizing at r, Q is linear in v with derivative
-8r. Its zero divisor meets v!=0 densely, so it defines a
height-one place of the original affine source.

The extension (11) is ramified of index two there. Equivalently,
the unique field derivation extending partial_v satisfies
\[
                             \partial_v w=-4r/w,              \tag{14}
\]
which has a pole at the branch. On the generic branch r and v
are units. Since this cover is ramified, (9) contradicts the
existence of any M satisfying (1)--(2) for a valued quantum
closure with residue field (11).

Thus the natural S3-closure route cannot pass through a finite
left target-Weyl carrier with a faithful full classical lattice,
no matter which fractional lattice is chosen. All images may
have arbitrary rational quantum corrections, including a
nonpositional first image; the proof used only their full
classical Hamiltonian frame. Arbitrary positive parameter
contact order e was included in (5)--(7).

The same theorem applies after identity stabilization: the
branch divisor remains ramified after adjoining spectator
variables. Regular polynomial Weyl source/target changes
with regular inverses transport the module and its classical
embedding, preserving the exclusion. The following admission theorems
remove the saturation assumption for classical graph carriers at
ordinary quantum parameter and for equivariant abelian closures of
the correct torsion-free classical rank at every parameter order.

## Controls and essential boundaries

A local model for the obstruction is x=z^m. On the basis
1,z,...,z^(m-1), multiplication Z and the unique derivative
connection C satisfy
\[
       Z^m=xI,\qquad
       C=\operatorname{diag}(0,1/(mx),\ldots,(m-1)/(mx)),
       \qquad Z'+CZ-ZC=Z/(mx).                                \tag{15}
\]
The fractional differences of connection residues cannot be
removed by one scalar potential. The proof using (J:J)
excludes every full fractional lattice, rather than only
this displayed basis. Finite controls check m=2,3,5.

Two positive or exception controls prevent overextension.
A constant extension theta^2=2 over Q is finite etale;
its derivative connection is regular and
A_R tensor_Q Q(sqrt(2)) is a finite saturated natural carrier.
The ramified quantum-parameter identity shear
x -> x,p -> p+t*x^2 with hbar=t^3 is also an exact finite
positive control. Its sub-quantum coefficient t*x^2 is
removed by the finite P_i correction above.

A nonsaturated commutative lattice can conceal its cover:
over R[x] inside K(x,sqrt(x))((t)), take
M_0=R[x]*1+R[x]*(x+t*sqrt(x)).
It is finite and t-torsion-free, but the nonzero class of
t*sqrt(x) in M_0/tM_0 has zero ambient residue.
It therefore fails (1). This control is not a Weyl right
action and does not prove that such a concealed closure
can complete Long.

Finally, a free matrix module with a componentwise derivative
connection is regular even when its underlying module is
identified with S[z]/(z^2-x). That connection sends the basis
vector z to zero and is not a derivation of the field
multiplication: differentiating z^2=x would give 0=1.
It does not have the natural scalar-potential form (7).
The theorem therefore excludes no arbitrary matrix repair
that changes the division-induced action or the faithful
classical field embedding.

The ordinary Long carrier has classical left lattice S itself;
its multiplier ring S is etale over S and passes this necessary
test. No contradiction for that original rank-one carrier
is obtained. The admission results below distinguish recoverable
nonsaturation from changes of the natural division action.

## Classical graph admission supplies a finite valuation level

There is a second route to the multiplier theorem that does not
assume saturation. Let M be finite on the left inside Dtilde,
and suppose its left D-span is all Dtilde. After one central
parameter shift, arrange M subset Otilde. Assume there are
polynomial left lifts P_i of psi_i such that
\[
                         (P_i m-mF_i)\in a(t)M
                         \quad\text{for every }m\in M.         \tag{16}
\]
At quantum order e=1, this is exactly the identified classical
graph condition: right source generators act on M/tM by the
left scalar polynomials psi_i. It is stronger than merely
set-theoretic support on that graph. At e>1, (16) is the
corresponding scalar graph congruence modulo the quantum
parameter, and is a separate hypothesis.

Condition (16) can be checked on finite left generators of M:
the extra term on a left multiple is [P_i,b]m, which is
divisible by a(t). It also forces v(F_i-P_i)>=e. Choose m of
valuation zero; its commutator with P_i has valuation at least e,
and the left/right difference lies in a(t)M subset a(t)Otilde.
Subtracting the two expressions gives the claimed valuation
for F_i-P_i.

Now nabla_i=(L_(P_i)-R_(F_i))/a(t) preserves M and its ambient
valuation filtration. Put
\[
 N_d=M\cap t^dO_{\widetilde D},\qquad
             J_d=\operatorname{res}(t^{-d}N_d)\subset\ell
                                      \quad(d\ge0).           \tag{17}
\]
Each N_d is a finite left A_R-module, because A_R is noetherian.
The quotient N_d/N_(d+1)=J_d is killed by t and is therefore
finite over S. Moreover J_d subset J_(d+1). The divided
operators induce {psi_i,-}-c_i on every J_d, with
c_i=res((F_i-P_i)/a(t)). The F_i may lie in Dtilde here:
the ell-valued scalar potentials still cancel on multipliers.

Some finite J_d spans all ell over k. Choose integral lifts
of a residue-field basis of ell/k. Fullness D*M=Dtilde and
left Ore denominator clearing put a nonzero polynomial
multiple of each lift in M after a central parameter power.
Their normalized residues are nonzero k multiples of the
chosen basis. All lie in one J_d for sufficiently large d.

Apply the multiplier argument (3), (8)--(9) to that finite
full J_d. Its multiplier ring is finite etale over S.
Thus a ramified ell/k admits no full finite natural carrier
satisfying (16). In particular the ordinary quantum-parameter
Long closure is excluded by its correct classical graph,
without saturation, equivariance, a free frame or a bound
on the quantum correction degrees.

For a ramified quantum parameter, classical graph agreement
only modulo t need not imply (16). A surviving ramified
construction would need a genuinely nonscalar action modulo
a(t), or a different admitted action. This is a necessary
condition, not a construction of such a model.

## A specialization rank bound for finite carriers

We next remove saturation in a different way, valid for every
parameter order and requiring only the classical rank rather
than (16). For any finite left A_R-module N,
\[
 \dim_D(D\otimes_{A_R}N)
              \le\operatorname{rank}_S(N/tN).                 \tag{18}
\]
Choose a finite presentation A_R^b -> A_R^c -> N.
If the special relation matrix has rank q over k, a square
q block has invertible reduction over k. Its lift is invertible
over the completed division ring by a geometric-series inverse,
hence has rank q already over D. The generic relation rank
cannot be smaller than the special relation rank. This proves
(18). It uses no projectivity theorem. It is the general-rank
version of the presentation argument in QMC1.

## Equivariance repairs the lattice by finitely many character shifts

Suppose D=Dtilde^H for a finite abelian group H fixing the
central coefficient field. Assume K contains all its character
values. Let ell/k have degree delta. Let M subset Dtilde be
full, finite as a left A_R-module, stable under H and under
right multiplication by the F_i in D. Assume only that M/tM
is torsion-free over S and has rank delta. The ambient residue
map need not be faithful, and (16) is not required.

QGD1's finite-quotient theorem gives
[Dtilde:D]_l=|H| and delta=|H_v| for its valuation stabilizer.
Fullness identifies D tensor M with Dtilde. Together with
(18) and residue-basis independence this gives
\[
                        |H|=\delta,\qquad H_v=H.              \tag{19}
\]
Thus H preserves the valuation; a parameter-pole action cannot
hide extra branches inside this minimal classical carrier rank.
This rank step does not use abelianness. For a degree-two
residue cover it already forces H to have order two, so the
character splitting below is automatic in Long application.

For each character chi, use the central-coefficient projector
\[
        e_\chi={1\over|H|}\sum_{\sigma\in H}
                     \chi(\sigma)^{-1}\sigma,\qquad
                             M_\chi=e_\chi M.                 \tag{20}
\]
They are left A_R-linear and commute with the right F_i.
The corresponding eigenspace of ell has k-dimension one.
Indeed the character projectors are nonzero by independence
of field automorphisms, and their |H| nonzero dimensions sum
to delta=|H|. Lifting and averaging gives a nonzero eigenvector
in Dtilde for each character. Since [Dtilde:D]=|H|, every
generic character space has left D-dimension one as well.

Apply (18) separately to the finite M_chi. Their special
ranks are at least one and sum to delta, so each is exactly one.
Each special module is torsion-free, being a direct summand of
M/tM. Choose the minimum ambient valuation m_chi of M_chi;
it is attained on finite generators. The module
\[
                  M^{\rm new}=
                       \bigoplus_\chi t^{-m_\chi}M_\chi       \tag{21}
\]
is finite on the left and stable under every right F_i.
It lies in Otilde. Each character summand has a nonzero
specialization map into its one-dimensional residue eigenspace.
A nonzero map from a torsion-free rank-one S-module into a field
is injective. Distinct character residue spaces are independent.
Consequently Mnew/tMnew embeds faithfully in ell and spans it.

This proves the saturation condition (1) for Mnew. The strict
finite-carrier theorem therefore makes ell/k unramified over S.
No assertion that the original M was saturated is needed.
No right action is changed, and only finite character
projections and central parameter powers are used.

For Long's quadratic source closure, all character values are
already +1,-1. Even if a proposed larger quantum closure starts
with more group elements, its degree-two special carrier rank
and (18)--(19) force the relative fixing group to have order two
and preserve the valuation. Hence no equivariant finite natural
closure carrier with a torsion-free degree-two special module
exists across Q=0, for any ramification order e and any
quantum correction degrees.

This closes the nonsaturated equivariant version of the
natural larger-closure route with the correct classical rank.
It does not exclude every matrix repair: the right action
must still be multiplication in the division closure, and
the residue extension must be the stated one. Nonabelian
relative groups of higher degree, torsion or excess rank in
the special module, and nonequivariant ramified models
without (16) require separate admission arguments.

## A genuine Weyl action shows why saturation itself is not automatic

There are full finite nonsaturated positive controls, beyond the
earlier commutative example. Let theta be a central constant
quadratic element, theta^2=2 over Q, and work in D(theta).
For every e>=1, with [p,x]=a(t), take
\[
                M=A_R\,1+A_R(x+t\theta).                      \tag{22}
\]
It is stable under the involution theta -> -theta and under
natural right x,p. In particular
(x+t*theta)*p=p*(x+t*theta)-a(t).
Its special module is free of rank two, with the correct
identity classical graph. Yet t*theta has a nonzero class
in M/tM and zero ambient residue.

The character summands are A_R*1 and A_R*t*theta.
Their valuation minima are 0,1; (21) replaces M by
A_R*1+A_R*theta. This repaired carrier is saturated.
Its multiplier ring is the allowed constant etale extension.

Over Q(omega), omega a primitive cube root, the same test with
theta^3=2 and basis 1,x+t*theta,x^2+t^2*theta^2 has character
valuation minima 0,1,2. The finite repair yields the constant
cubic carrier. The exact checks include the full Weyl right
relations, group equivariance, projector identities and
classical graph, so neither parameter flatness nor the full
right action may be cited as automatic saturation.

Torsion in the special module is a separate issue. The actual
bimodule I=A_R*x+A_R*t is finite and stable under both natural
Weyl actions. Its special module has a nonzero class of t killed
by x, so the torsion-free hypothesis of the character repair fails.
At a(t)=t, [p,x]=t is not in tI; this module also lacks the
identified classical graph, despite its full right Weyl action.
Indeed 1 is not in I after reduction modulo t, so t is not in tI.

The ordinary rank-one Long graph remains outside the
ramified-cover contradiction: its residue field is k and
its multiplier ring is S. These stronger admission tests
construct no such Long carrier.

## Verification boundary

The checker verifies the complete cubic splitting formulas,
S3 action identities, discriminant pullback, the interior
branch derivative, companion connections for m=2,3,5,
a constant etale extension, the ramified-parameter finite
shear, and both the nonsaturation and matrix-connection
boundaries. It also checks genuine quadratic/cubic constant-cover
Weyl carriers, equivariance, character projectors and finite repairs.
It binds the retained QGD1, QFC1, QDI1 and QMC1 evidence.

The arbitrary-field, arbitrary-lattice, parameter-ramified
finite-carrier theorem, multiplier finiteness, potential
cancellation, Fitting/etale argument, interior-divisor
identification, finite valuation-level graph admission, specialization
rank bound, arbitrary split-abelian character repair and every-degree
Long closure exclusion are written proofs. No universal automatic
saturation or descent to a scalar frame is inferred. Finite controls do not exhaust
lattices, orders, gauges or matrix actions.

No actual finite Long carrier, universal saturation theorem,
arbitrary matrix exclusion, all-representative nonlifting,
scalar/matrix counterexample, DC2/JC2 resolution,
independent theorem replay, formal verification, external
review or novelty priority is claimed.

Run the declared exact controls from the repository root:

```sh
.venv/bin/python research/jacobian-followup/dixmier-formal-lifting/check_ramified_closure_carriers.py
```
