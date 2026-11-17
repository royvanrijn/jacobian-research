# Finite quantum symmetries and Galois residue fields

15 November 2026. **QGD1.** No actual full quantum lift of Long can
realize its image division ring as the fixed ring of a finite group
of ambient division automorphisms. This includes parameter-pole
automorphisms, nonpositional quantum corrections, nonfree rank-one
projective carriers and identity stabilizations. The obstruction
uses the full classical function-field extension, not a chosen
position centralizer or an assumed invariant finite order.

More generally, finite group quotients of discretely valued division
rings with a central uniformizer and commutative characteristic-zero
residue have Galois residue extensions. A direct approximation argument
accounts for group elements moving the valuation. The generic indices
equal the group order; the selected completed index is the order of
the valuation stabilizer. Their difference counts the other valuation
branches. No finite-dimensional-over-the-center hypothesis is used.

For an actual bounded Long lift, all automorphisms fixing its full
image form a finite group. If its order is g and the generic indices
are N_l,N_r, then g divides both indices and N_l,N_r>=3g. The selected
completion of the fixed division ring of this group is the entire
ambient completion. In particular either generic index below six
forces this automorphism group to be trivial. The completed Long
extension always has trivial automorphism group, even without a
spatial bound.

These statements exclude finite symmetry-quotient constructions.
They do not construct a finite Long carrier or exclude every
nonsymmetric division realization.

## A finite quotient and its valuation branches

Let D be a division ring with a discrete valuation v, value group Z,
central uniformizer h, and commutative residue field k of
characteristic zero. Let G be a finite group of automorphisms fixing h.
Set
\[
 W=D^G,\qquad e=\overline W\subset k,\qquad
 H=\{g\in G:v\circ g=v\},\qquad r=|G|/|H|.                    \tag{1}
\]
There is no assumption that G preserves v. The distinct valuations
v_i=v composed with g_i form an orbit of size r. Their restrictions
to W are identical, with value group Z and uniformizer h.

We first construct branch selectors. For v_i!=v_j choose a in D^*
with v_i(a)>v_j(a), replacing a by its inverse if necessary. Put
\[
 b=a^2h^{-v_i(a)-v_j(a)},\qquad
                       t_{ij}=(1+\lambda b)^{-1}.             \tag{2}
\]
Then v_i(b)>0 and v_j(b)<0. Choose nonzero rational lambda so
1+lambda*b has nonzero residue at every v_l where v_l(b)=0.
Only finitely many rational values are forbidden, because the
residue fields have characteristic zero. At all other branches
there can be no cancellation. Hence every v_l(t_ij)>=0,
res_i(t_ij)=1 and v_j(t_ij)>0.

Take c_i to be an ordered product of t_ij over j!=i. It satisfies
\[
 v_i(c_i)=0,\quad\operatorname{res}_i(c_i)=1,\qquad
                            v_j(c_i)>0\quad(j\ne i).           \tag{3}
\]
The products need not commute. Multiplicativity of valuations and
of each residue map gives exactly the stated properties.

## Averaging after selecting one branch

The stabilizer H acts on k, and residues of W are fixed by that
action. Conversely choose alpha in k fixed by H and a lift a with
v(a)=0, unless alpha=0, which is immediate. Select c for the branch v.
For sufficiently large n, z=c^n*a has positive valuation at every
other branch and still has residue alpha at v. Therefore
\[
 b={1\over|H|}\sum_{g\in G}g(z)\ \in W,\qquad
                      v(b)\ge0,\quad\overline b=\alpha .       \tag{4}
\]
Only the |H| terms from the stabilizer contribute to the residue;
all other terms have positive valuation. This proves
\[
                              e=k^{\overline H}.              \tag{5}
\]
Unweighted averaging before selecting a branch would not prove this:
terms from other valuations can have parameter poles.

The action of H on k is faithful. Suppose sigma in H has finite
order d and identity residue action. Choose a of valuation zero
with sigma(a)!=a and write sigma(a)-a=h^m*b with m>=1 and v(b)=0.
Since sigma fixes h, preserves v and is identity on residues,
every sigma^j(b) has the same nonzero leading residue as b. The sum
sigma^d(a)-a=sum_(0<=j<d) sigma^j(sigma(a)-a) then has leading
coefficient d*res(b), nonzero in characteristic zero. This
contradicts sigma^d=1. Thus sigma=1. By (5) and the ordinary
Artin fixed-field theorem,
\[
              k/e\text{ is Galois},\qquad [k:e]=|H|=\delta .   \tag{6}
\]

## Generic indices and the precise completion loss

The elementary Artin dimension bound also works for division rings:
\[
                            [D:W]_l,\ [D:W]_r\le |G|.         \tag{7}
\]
Here is the right proof, so no noncommutative Galois theorem is
imported. If m>|G| elements b_j were right W-independent, the
matrix (g(b_j)) with |G| rows would have a nonzero right null
vector over D. Choose a null vector of smallest support and
normalize one nonzero coefficient to 1. Applying any element of G
permutes the equations and gives another null vector with that
coefficient 1. Their difference has smaller support, so is zero.
All coefficients are therefore G-fixed. The identity row is now a
nonzero right W-relation among the b_j, a contradiction. Apply this
argument to the opposite ring for the left bound.

For the lower bound, choose at each branch lifts b_ia of an e-basis
of its residue field, of size delta. Multiply these lifts by a
sufficient power of c_i to obtain z_ia of valuation zero at branch i
with those basis residues, and of positive valuation at every
other branch. A right W-relation among the r*delta elements z_ia
can be normalized by a central h power so all coefficients have
nonnegative valuation and at least one has valuation zero.
Reduction successively at each v_i forces all coefficients to have
positive valuation, a contradiction. The same left argument gives
\[
                         [D:W]_l,\ [D:W]_r\ge r\delta .        \tag{8}
\]
Combining (6)--(8) yields, without any general index-symmetry theorem,
\[
 \begin{split}
 [D:W]_l=[D:W]_r&=|G|=r\delta,\\
 [D:W]-\delta&=\delta(r-1).
 \end{split}                                                  \tag{9}
\]
Thus the selected residue extension must be Galois even when G
moves the valuation. Such movement can increase the generic index;
it cannot change a non-Galois residue cover into a possible quotient.

Write Dhat and What for the selected completions. The stabilizer H
extends to Dhat and fixes What. The fixed subring Dhat^H is complete
and has residue e: averaging a lift over H proves residue
surjectivity. The shared-central-uniformizer basis lemma of
[QDI1](QUANTUM_DIVISION_AND_COMPLETION.md) gives
[Dhat:Dhat^H]=delta=[Dhat:What], so
\[
                               \widehat W=\widehat D^H.       \tag{10}
\]
This is a statement about division rings and their h-valuations.
It does not identify finite polynomial orders or supply a finite
left Weyl-module lattice.

## Rigidity beyond finite order

Suppose now only that W subset D has finite separable residue
extension k/e, still with the shared central uniformizer.
Every valuation-preserving automorphism fixing W is determined by
its residue action:
\[
 \operatorname{Aut}_{W,v}(D)\ \lhook\joinrel\longrightarrow\
                                      \operatorname{Aut}_e(k). \tag{11}
\]
For the kernel, if sigma is identity on residues, sigma-id maps the
valuation ring into h times itself. If sigma!=1, choose the smallest
integer m>=1 attained by v(sigma(a)-a) on the valuation ring.
The map res(a) -> res(h^-m*(sigma(a)-a)) is a well-defined nonzero
derivation of k over e. Changing a by h*b changes this expression
by positive valuation; multiplicativity gives the derivation rule.
A derivation of a separable algebraic field extension over its base
is zero, a contradiction.

If either generic index [D:W]_l or [D:W]_r equals delta, every
W-automorphism is automatically isometric. Lifts of a residue basis
are a generic basis with the exact valuation formula (QDI1).
Writing an automorphism in this finite basis gives
v(sigma(a))>=v(a)+C for one constant C. Apply this to a^n and divide
by n to get v(sigma(a))>=v(a). Its inverse gives equality.
In particular this holds for every automorphism of the completed
extension, whose dimension is delta. No continuity assumption is
needed in that completed situation.

## The full quantum frame has scalar centralizer

For a full formal polynomial symplectic lift over L=K((h)),
the centralizer of its completed image in the ambient completion
is exactly L. A commuting element, normalized to valuation zero,
has residue a satisfying {psi_i,a}=0 for all classical images.
Their Hamiltonian vector fields form a full frame, since the
symplectic Jacobian is invertible. Hence a is in K.
Subtract that scalar residue, divide by h and repeat.
Completeness identifies the element with a scalar Laurent series.
Consequently an actual full lift inside the Weyl division ring has
\[
                           \operatorname{Cent}_D(W)=L.        \tag{12}
\]
The same conclusion holds in a faithful formal rank-one frame of a
projective carrier, as in [QMC1](FINITE_GRAPH_CARRIERS_AND_MATRIX_REALIZATION.md).
No joint-position-centralizer assumption is required.

If the generic indices N_l,N_r are finite, (12) implies that
Aut_W(D) is finite of order at most min(N_l,N_r). To see this,
distinct W-automorphisms cannot differ by an inner automorphism:
an implementing element would centralize W and lie in L.

Automorphisms inequivalent modulo inner automorphisms are
left D-linearly independent as functions on D. For a shortest
relation sum a_i*sigma_i(x)=0, substitute y*x and subtract
a_1*sigma_1(y)*a_1^-1 times the relation. Minimality forces
a_i*sigma_i(y)=a_1*sigma_1(y)*a_1^-1*a_i for every y, making
sigma_i an inner conjugate of sigma_1, a contradiction.
The right W-linear maps D->D have left D-dimension N_r,
and the opposite-ring argument gives N_l. This proves
\[
                         |\operatorname{Aut}_W(D)|
                                      \le\min(N_l,N_r).       \tag{13}
\]

## Long has no finite quantum quotient

QDI1 identifies Long's full classical field extension with
K(r,s,t,d) subset K(r,s,t,d)(a), where
\[
 \begin{split}
 f(a)&=ra^3-2a^2+sa+4t=0,\\
 \Delta&=4s^2-4rs^3+128t-144rst-432r^2t^2 .
 \end{split}                                                  \tag{14}
\]
The independent variable d is Long's corrected fourth image.
The extension has degree three, since t is a cubic polynomial
in a over K(r,s,d), as proved by the retained birational chart.
Moreover
\[
                  \operatorname{disc}_t(\Delta)
                                  =256(4-3rs)^3\ne0.          \tag{15}
\]
Thus Delta is squarefree of degree two in t and is not a square
in K(r,s,t,d), even after algebraically extending K. The irreducible
cubic has Galois closure group S3 and trivial deck group.
Adding independent spectator variables preserves these properties.

Every full formal Long lift has this residue extension by QDI1.
Equation (11) and the completed dimension formula imply that
Aut_What(Dhat) is trivial. For an actual finite lift, or a finite-left
full graph carrier in QMC1's rank-one division realization, suppose
W=D^G for a finite group of L-linear division automorphisms.
Equation (6) would make the classical extension Galois, contradicting
(14)--(15). This rules out every such finite quotient, regardless
of quantum degrees, parameter poles or whether G preserves the
chosen valuation.

There is also a restriction on symmetries when W is strictly smaller
than their fixed division ring. By (12)--(13), take the entire finite
group A=Aut_W(D), of order g, and put F=D^A. Its valuation stabilizer
is trivial by (11), because the classical deck group is trivial.
Apply (5)--(9) to D/F: the residue field of F is all of k,
[D:F]_l=[D:F]_r=g, and F completes to Dhat. Residue independence
for F/W still gives [F:W]_l,[F:W]_r>=3. Therefore
\[
       g\mid N_l,\quad g\mid N_r,\qquad
                  N_l,N_r\ge3g,\qquad\operatorname{cl}_h(F)=\widehat D.
                                                               \tag{16}
\]
Here cl_h denotes closure in the selected completion. Any
nontrivial quantum deck symmetry must move
the valuation and occur entirely in the generic dimensions lost
at completion. If either N_l or N_r is below six, g=1.
No equality of arbitrary generic left and right indices is used.

All statements survive identity stabilization. Source division
automorphisms leave the full image division ring unchanged; ambient
division automorphisms transport it and conjugate any proposed group.
Thus the finite-quotient exclusion also transports through arbitrary
L-linear rational source/target division changes, even ones with
parameter poles. No extension of those changes to the selected
completion is assumed. This concerns the full image field; it does
not require fixing a position field as QCI1's rational-polarization
statement does.

## Controls with and without valuation movement

Over constants containing a primitive mth root zeta, the rational
pair F_m=x^(1-m)*p/m,T_m=x^m of
[QCI1](CENTRALIZER_INDEX_AND_RATIONAL_POLARIZATION.md) is fixed by
x -> zeta*x,p -> zeta^-1*p. Its whole division index is m by
QCI1, so the Artin bound and the division tower identify its image
with the fixed division ring. Its completed residue degree is m,
the group preserves v, and its classical cover is cyclic Galois.
This is allowed by (6). For m=3 it demonstrates that degree three
alone is not an obstruction. QCI1 separately excludes a finite
Weyl lattice for this displayed rational action.

A second rational pair has a finite symmetry with parameter poles:
\[
 \begin{split}
 U&=x+hx^2,\qquad V=(1+2hx)^{-1}p,\qquad [V,U]=h,\\
 \sigma(x)&=-h^{-1}-x,\qquad \sigma(p)=-p .
 \end{split}                                                  \tag{17}
\]
The involution fixes U,V and their image division ring has index two,
so that image is D^<sigma>. Its classical map is the identity,
of completed degree one. The two valuation branches have trivial
stabilizer, giving N=2,delta=1,r=2 as in (9).

For a direct check of (2)--(4), c=(1+hx^2)^-1 has residue one at v
and positive valuation at v composed with sigma. With z=c^2*x,
the invariant z+sigma(z) has residue x. This verifies why valuation
movement cannot simply be discarded or repaired by averaging an
arbitrary lift before controlling its other branches. The rational
operator V has Bernstein degree zero and a Weyl mate, so QCI1 also
excludes a finite diagonal carrier for this second control in every
positive rank. Neither rational control is a polynomial Long realization.

## Relation to fixed ring rigidity and remaining work

[Tikaradze, Theorem1](https://arxiv.org/html/1708.07923)
classifies domains having a differential-operator ring as a finite
group fixed ring. In particular a complex Weyl algebra cannot be
the fixed ring of a nontrivial finite group acting faithfully on
a domain. That theorem concerns actual rings. The present
division argument allows cyclic fixed Weyl division rings, as the
first control shows, and uses their Galois residue property to
exclude the Long reduction. It does not assume that the group
preserves a finite Weyl order or that its invariant order equals
the polynomial image. Tikaradze's theorem is comparison context,
not an input to (1)--(17).

A useful remaining construction must realize the non-Galois
residual extension with a finite carrier and the full Long action.
A cyclic kernel on a fixed elliptic energy curve is not, by itself,
a finite quotient of the full quantum division realization.
A different nonsymmetric embedding or a deformation moving the
energy remains possible. A larger Galois closure with D as an
intermediate subgroup-fixed division ring is also outside the
direct-quotient exclusion. Classically, Long has an S3 closure
with source point stabilizer S2; realizing that diagram quantum
mechanically would still require a compatible finite carrier
and descent across its branch boundary. No finite-group normal
form for all division embeddings is asserted.

## Verification boundary

The exact checker binds the retained QDI1/QCI1 inputs, verifies
the cubic discriminants and squarefree condition, the cyclic
rational positive control, the pole-bearing involution with
its full Weyl relation, its branch selector and its averaged
residue, and representative finite monodromy actions.
The branch approximation, characteristic-zero kernel rigidity,
Artin dimension bounds, generic/completed Galois theorem,
scalar centralizer, automorphism finiteness, all-degree Long
exclusion, all identity stabilizations and rational transport
are written proofs. Finite controls do not exhaust possible
groups, operators or degrees.

No actual finite Long carrier, scalar/matrix counterexample,
all-representative nonlifting, general higher-rank matrix
exclusion, DC2/JC2 resolution, independent theorem replay,
formal verification, external review or novelty priority is claimed.

Run the declared exact controls from the repository root:

```sh
.venv/bin/python research/dixmier/dixmier-formal-lifting/check_quantum_galois_residues.py
```
