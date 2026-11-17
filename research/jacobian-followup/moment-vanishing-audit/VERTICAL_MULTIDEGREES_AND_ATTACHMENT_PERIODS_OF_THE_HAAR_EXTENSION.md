# Vertical multidegrees and attachment periods of the Haar extension

16 November 2026. Written local geometric calculation in every degree.
Both FULL exceptional Mathieu--Zhao cases remain **OPEN**. External
review remains deferred. No degree campaign or holonomic run is used.

The corrected intrinsic coupling can be calculated from vertical
intersection degrees and the components met by the node sections.
This removes the smooth-normalization restriction of the preceding
[good-node lemma](BOUNDED_GRAPH_DESCENT_AND_GOOD_NODE_EXTENSION.md).
It also explains an important model dependence: after special-fibre
blowups ANY corrected toric line bundle can extend invertibly. Its
nonzero component then survives in vertical multidegrees, which must
not be dropped or summed away.

## 1. Keep the actual original corrected point

Start with the canonical ORIGINAL boundary section for a phase on
SL2/2O or SL2/2I. As in the
[toric correction theorem](ACTUAL_TORUS_CORRECTIONS_AND_RAMIFIED_COMPONENT_VALUATIONS.md),
on a finite splitting and labeling cover s=u^e choose

    t=[n]sigma-b in the ACTUAL generic torus,
    b in G^0(C[[u]]), pi(b)=[n]pi(sigma),
    v=val(t)/(e*n).                                  (1)

No volume or Jacobian is substituted. The generic torus and its character
comparison are those of the ORIGINAL uncut fibre source. The preceding
clutching-symbol cut graph is not being identified with it. A supplied
isogeny comparison can be cleared by the same integer/cover choices;
all invariants and zero tests below are rational.

Put R=C[[u]], K=C((u)). Let Gamma be the generic boundary graph realizing
the relevant actual Picard torus, and let L_K be the generic toric line
bundle representing t. It is trivial on each generic normalization.

For every generic vertex v, suppose a proper regular STRICT SEMISTABLE
normalization model Y_v/R has been supplied. Its special fibre is a
connected reduced union of smooth proper curves E_(v,j), with transverse
nodes. The marked half-node sections are disjoint and meet smooth points
of this special fibre. Write j(a,v) for the component met by the
half-node of edge a at v. Gluing the marked sections gives a proper
relative curve C/R. Extra vertical nodes INSIDE Y_v are retained.

These are explicit model hypotheses. We do not infer them from volume
regularity, an ambient tree or a picture of the generic graph.

## 2. The exact intersection formula

Suppose L_K extends to an invertible line bundle L on C. Pull it back
to Y_v. Since its generic restriction is trivial, a generic trivialization
has a vertical divisor

    D_v=sum_j d_(v,j) E_(v,j),
    L|Y_v ~= O(D_v).                                 (2)

The coefficients are integers. A different generic trivialization multiplies
by a scalar in K*, adding the same integer to every d_(v,j).

Let M_v be the intersection matrix (E_(v,j) dot E_(v,l)), and let

    delta_(v,j)=degree(L|E_(v,j)).

Then

    delta_v=M_v d_v.                                 (3)

For an oriented generic edge a:tail->head, use the convention that its
gluing scalar carries the tail generic basis to g_a times the head
generic basis. In local R-bases the gluing of the INVERTIBLE extension
is a unit. The generic bases have orders d at their marked sections.
Consequently

    ord_u g_a=d_(tail,j(a,tail))-d_(head,j(a,head)).    (4)

Thus for EVERY generic graph cycle z the ORIGINAL corrected coupling is

    v(z)=(1/(e*n)) sum_a z_a
           [d_(tail,j(a,tail))-d_(head,j(a,head))].    (5)

No parameter-dependent unit contributes an extra valuation once the
actual invertible extension has been retained. This does not drop such
units from a period identity; it only evaluates their orders in (1).

Formula (5) is independent of all generic trivializations: the change
is a vertex coboundary, whose sum on z is zero. It is also independent
of the chosen invertible extension and further models, because (4)
recovers the character valuations of the SAME generic point t. Further
ramified covers scale these orders and e together, exactly as in (1).

## 3. Multidegrees and attachment positions suffice

The whole special fibre of Y_v is div(u), so M_v times the all-ones
vector is zero. Distinct components meet in their number of nodes.
It follows directly that

    d^t M_v d=-sum_(j<l) #(E_j intersect E_l)(d_j-d_l)^2.  (6)

Connectedness gives ker(M_v)=Q*1. Hence delta_v determines d_v up to
the constant already harmless in (5). This gives a numerical interface
to the existing geometric extension calculation: supply the ACTUAL
vertical multidegrees and the components met by the marked sections.
Solve (3), then evaluate (5). No phase-degree cutoff is involved.

Equivalently define J_v(z) on the special components by adding z_a at
the tail attachment and subtracting z_a at the head attachment. For a
cycle, sum_j J_(v,j)=0. Solve M_v q_v=J_v over Q; a solution exists
by (6), unique modulo constants. Symmetry of M_v gives the same coupling:

    v(z)=(1/(e*n)) sum_v delta_v dot q_v(z).           (7)

Both (5) and (7) keep each vertical component and each attachment
separately. Total degree zero, or cancellation after summing an orbit
or taking a scalar norm, is insufficient.

A useful sufficient geometric case is immediate. If, at each generic
vertex, all half-nodes belonging to the cycles under consideration meet
ONE special component, their contribution to J_v is zero by cycle
balance. Then every such pairing vanishes, independently of unused
degrees elsewhere. This needs a verified invertible extension on that
model; a modification which moves one attachment to a new component
cannot be silently contracted in this test.

## 4. Invertibility alone is not invariant under resolution

Given the supplied normalization models and ANY generic toric point t,
one can make an invertible extension after modifications at special
marked points. Represent t by g_a=u^(m_a) epsilon_a, epsilon_a in R*.
For each edge choose its tail half-node. Blow up its marked special
point on Y_tail once. The marked section now meets the exceptional
P1 E_a at a smooth point, away from the strict special-fibre intersection.
Different half-nodes use disjoint points. The blowup is regular, proper,
unchanged generically, and its special fibre is still reduced semistable.

On each modified Y_v put

    L_v=O(sum_(a with tail v) m_a E_a).               (8)

Its distinguished generic trivialization has order m_a at the chosen
tail section and order 0 at the other half-node of that edge. Thus
the desired identification g_a is, in the local R-bases, exactly a
unit (including the invertible local factors of these trivializations).
Glue by those units. The resulting L is INVERTIBLE on the modified
relative curve, with generic fibre the original L_K.

This construction also works for negative m_a using negative divisors,
and for a self-edge using its two distinct marked half-nodes. One
blowup at the chosen half-node suffices for arbitrary m_a; the exponent
is retained as a line-bundle coefficient, not as a discarded component.

Therefore the existence of an invertible extension on SOME expanded
model does NOT imply zero intrinsic component. The obstruction can
move from a nonfree node to nonzero vertical degrees. The preceding
good-node lemma remains correct on its stated smooth-normalization,
persistent-node model; dropping that hypothesis without tracking
(3)--(5) would be incorrect.

## 5. The model-stable zero test and contraction gate

Under the supplied normalization-model hypotheses, the following are
equivalent for a good generic subgraph and its actual corrected point:

    (a) its good-cycle coupling is zero;
    (b) its generic line bundle has an invertible extension
        whose degree is zero on EVERY vertical normalization component.
                                                               (9)

For (b)->(a), (3),(6) make each d_v constant. Equation (5) is then a
generic vertex coboundary. Conversely (a) means the good torus point
belongs to T_good(R), by the preceding bounded-descent lemma. Choose
unit gluing representatives and glue the TRIVIAL line bundles on all
Y_v. This is an invertible extension with all multidegrees zero.

For a given extension, zero multidegree is sufficient but not necessary:
nonzero degrees can have zero attachment periods in (5). What is
equivalent in (9) is EXISTENCE of an extension with all these degrees
zero. The exact test on any supplied extension is always (5) or (7).

Once these actual good-cycle periods vanish, the already proved bounded
descent gives v=p^*v_contracted. Contracting a good piece is legitimate
only AFTER that test, with every remaining attachment retained. A
forest contracted FIBRE graph then implies zero full component and,
under eventual purity, safety of every original marker. If cycles
remain, their corrected periods must still be evaluated.

This is a local component calculation, not an assertion of a surface
blowdown or a global unit-retaining quotient of variations.

## 6. Actual Long data: the lost component has degrees +1 and -1

The preceding [original-volume control](BOUNDED_GRAPH_DESCENT_AND_GOOD_NODE_EXTENSION.md)
uses the known nonexceptional pure phase p=(1+y)(1+z+yz). Its ACTUAL
generic fibre has the nearby model A2_(X,W) minus {XW=-s}, with boundary
triangle H_s,D_X,D_W in P1 times P1. The actual corrected torus section,
in the retained Kummer normalization, is

    t=s/(s+1).                                       (10)

The original volume in these coordinates was checked as
(s+1)dX wedge dW wedge ds/(XW+s)^2. Its two infinity components form
a volume-good tree, but contracting that tree leaves a genuine cycle.

At s=u=0 the NORMALIZATION FAMILY of H_s is

    Y_H={XW+u=0} in P1_X times P1_W times Spec R.

It is the blowup of P1_X times Spec R at X=u=0: it resolves the rational
map W=-u/X. It is a regular semistable surface with two special curves,
the strict C and exceptional E. Its intersection matrix and the divisor
D_H=E give

    M_H=[[-1,1],[1,-1]],
    d_H=(0,1), delta_H=(1,-1).                      (11)

The H-D_X attachment X=infinity lies on C. The H-D_W attachment X=0
lies on E. The other two normalization line bundles may be trivial.
Choose the generic gluing t on the oriented edge H->D_W and unit
gluings on the other two edges. The line bundle O(E) on Y_H absorbs
the factor u in t=u/(1+u), so all node identifications are units and
the generic toric line bundle extends INVERTIBLY on this model.

Nevertheless (5) on the oriented triangle gives

    ord_0 t=1,                                      (12)

because its two H attachments meet components with potentials 0 and 1.
Summing the degrees in (11) gives zero and loses the obstruction.
Discarding E does the same. This computes the surviving cycle coupling
from explicit actual vertical intersection data, rather than inferring
it from the shape of an ambient resolution tree.

The example remains on ordinary SL2 and is NOT an exceptional pure
counterexample. The original special v=0 phase component omitted by the
generic v!=0 chart remains part of the completed actual source which
supplied (10). This boundary-model calculation does not replace that
source or assert that the generic chart is a complete central model.
The existing exceptional Bezout pullbacks retain the nonzero intrinsic
component and their fixed original marker, but unweighted purity and
volume regularity are not transferred by dropping their Jacobian.

## 7. The all-degree exceptional question is now a numerical geometric one

For the ACTUAL corrected original section, obtain its Picard line-bundle
extension on a supplied normalization model and retain its degrees on
ALL vertical components. The formula (5), or equivalently (7), is the
precise bridge from that model to the intrinsic cycle coupling. It
applies equally to degenerating positive-genus normalizations after the
actual abelian correction; no independent genus obstruction is added.

The remaining full implication is that unweighted purity forces these
actual attachment periods to vanish for every relevant good cycle and
every residual contracted cycle, in BOTH exceptional cases. Original
volume regularity and the unique dominant leading root have not yet
been proved to force this numerical condition. Capacity, automorphism
and primitive constraints remain inputs to such a proof, not substitutes.

No actual exceptional pure phase with nonzero period and an explicit
fixed original bad marker is constructed. No ordinary marked-resolvent
algebraicity or generating-function transfer is inferred. The completed
family-transfer result and all prior goals remain intact. Both FULL
cases remain OPEN; external review remains deferred.
