# Prefix excess closes the genuine iterated quartic class

15 November 2026. PF2RULED96 excludes the remaining complementary
A4 earlier-reset families, at ALL companion depths and ALL finite
numbers of nonsolid cyclic3 loops and simple tails. Together with
PF82--95, the ENTIRE genuine positive iterated-torus quartic end
class is excluded in BOTH source pencils, at every finite length.

This statement retains the whole smooth normalization, entire source
boundary charts, actual preferred peripherals, genuine positive
stages and simple-meridian full-final-S4 hypotheses. It is not a
classification of degenerate stages or higher-degree packets, a
finite Keller-map theorem, or completion of the full research goal.

## 1. Local prefix identities at every characteristic depth

Use the full unramified D chart and contact diagonal n_D of
PF77/79/89. Every support index is at least sigma+1, and its
terminal free order kappa is at least2. Truncating after a
characteristic node and dividing sigma,p and retained support
indices by the remaining gcd gives a primitive LOCAL chart.
Such a prefix need not separately extend to an affine surface.

Conversely, append a genuine last winding a>=2 to a prefix chart
with sigma_0,p_0 and free order kappa_0. Write kappa for the new
free order. The new last support index and its predecessor are

    j_new=a*sigma_0-1+kappa,
    j_old=a*(sigma_0-1+kappa_0).

If the old chart has a single node, its relevant last index is p_0;
the same formula uses kappa_0=p_0+1-sigma_0. Lower noncharacteristic
terms may remain, but cannot weaken j_new<j_old. Thus

    D=j_old-j_new=a*(kappa_0-1)-kappa+1 >=1,
    kappa<=a*(kappa_0-1).                             (1)

D is exactly the positive primitive fibre intersection at this edge,
by PF95's characteristic-difference law. A valid appended chart also
has gcd(a,kappa-1)=1. We keep that primitivity condition when using it.

For a FIXED positive INTEGER h and integer q, set
delta=h*n_D-q^2*p_D^2. PF89's exact
strip identity gives, for every such append,

    delta_new=a^2*delta_old-h*(a-1)*(a+kappa).          (2)

This is an integer identity on the local charts, at arbitrary depth.
No global affineness, group or monodromy is assigned to a prefix.
In particular if the FULL chart has delta>0, EVERY primitive prefix
has positive delta: solve(2) backwards, adding a positive term.

Define the terminal excess E=delta-h*kappa. Equations(1)--(2) give
the exact positive-edge recurrence

    E_new=a^2*E_old+h*a*D.                            (3)

Therefore E_old>=0 forces E_new>0, and this property persists through
every later valid extension. If the full chart ends with E<0, then
EVERY one of its prefixes has E<0. Index order and positive D are
essential; the determinant recurrence alone does not imply this.

A second general consequence is LOW-DETERMINANT TERMINALITY. If
1<=delta_old<h, no further genuine node with kappa>=2 can keep
positive determinant. Indeed

    delta_new<=a^2*(h-1)-h*(a-1)*(a+2)
             =-a^2-h*(a-2)<=-4.                     (4)

Once the determinant is negative, every later append in(2) keeps it
negative. These two prefix lemmas are independent of field degree,
support size and characteristic depth. They require the stated
positive h, genuine node and whole-chart index conditions; no
arbitrary higher-degree packet is assumed to satisfy them.

## 2. A balanced odd block has determinant exactly5

Consider a full source chart with h=3, delta_final=2. Since its
terminal free order is at least2, its final excess is negative.
By Section1 every prefix has positive determinant and negative excess.

Suppose its root ratio is

    Q=3b+2, Q,b odd, gcd(Q,b)=1,                      (5)

and an initial block of its rooted D path has only odd windings
prime to3 and EVEN fibre intersections between adjacent block nodes.
The block can have any positive finite length, including just its root.
The winding orders are the characteristic gcd ratios, as in PF83.

Cut immediately after this block and divide by the remaining gcd.
PF95's edge formula is unchanged under this operation: both the
characteristic difference and its gcd denominator are divided by
the same integer. The root ratio and all normalized outer cone
residues are likewise retained. This is only a local characteristic
cut, not a new affine source or a replacement finite cover.

In the cut chart p and sigma are odd: p is the product of the odd
windings, and sigma=b times their inner product. Starting from p,
each next characteristic index differs by EVEN D times an odd gcd.
All its relevant indices are therefore odd. Every primitive prefix
of the block consequently has ODD kappa as well. This is exactly
where the even-intersection hypothesis is used.

At its single-node root, p=Q,sigma=b, so

    n_D=(b+1)*Q-b+1,
    kappa_0=Q+1-b=(2Q+5)/3,
    delta_0=3*n_D-Q^2=5.                             (6)

Now assume a block prefix has determinant5 and append its next
odd winding a, prime to3. It has a>=5 and odd free order kappa.
Equation(2) becomes

    delta_new=5a^2-3(a-1)(a+kappa)
             =3*kappa+a*m,
    m=2a+3-3*kappa EVEN.

The negative-excess condition forces m<0, hence m<=-2. If m<=-4,

    delta_new=2a+3+(a-1)*m<=7-2a<0,

contrary to prefix positivity. The ONLY possibility is m=-2, giving

    delta_new=5, kappa=(2a+5)/3.                     (7)

Induction proves that EVERY prefix of the whole odd block has
determinant EXACTLY5, however long the block. Each of its non-root
windings is in fact5 modulo6. No finite parameter census is involved.

## 3. Appending a winding2 arrival leaves no room for a reset

After a block in Section2, append a genuine winding2 node. Its
primitive free order is EVEN because gcd(2,kappa-1)=1. Formula(2)
with delta_old=5 gives

    delta_arrival=14-3*kappa.

Prefix positivity restricts the even integer kappa>=2 to2 or4.
At kappa=2 the determinant is8 and the excess is8-6=2>=0,
forbidden by(3) and the negative final excess. Hence necessarily

    kappa_arrival=4, delta_arrival=2.                 (8)

But delta_arrival=2<h=3. Formula(4) makes ANY further genuine
characteristic node impossible in a full positive-determinant chart.
Thus no such full source can contain an odd/even-bridge initial
block, then a winding2 node, and then ANY additional genuine node.

This is a source obstruction with no bound on the odd block, the
later windings, the final free order or field degree, under its
explicit local chart and determinant hypotheses.

## 4. Every remaining A4 earlier-reset trajectory has that shape

PF95 leaves only complementary incoming A4 histories with an EARLIER
nonsolid reset. Their actual closed end is

    H1=Z/2, [R]=0.

Every passing later simple tail has odd winding and odd parameter;
even-parameter tails are already excluded by the canonical cone
budget. Their normal/fibre distance2 forces the rho=2 source collar,
so h=3,q=1,delta_final=2. The zero-R full source lattice makes
3 not divide p_D, hence every D winding is prime to3 (PF83).

The source root is the FINAL simple tail, with order Q and first
ratio(Q,b). PF89's whole contact bound excludes Q=3b-2;
therefore Q=3b+2. Since Q is odd, b is odd. The outer block of
the canonical target path consists precisely of ALL simple tails,
read from the last one back toward the first. It has at least one
node. The actual primitive fibres of two consecutive odd-parameter
simple tails have intersection

    2*(r*B_previous-k)>0,                             (9)

by PF83, an EVEN integer. All these edges are canonical and survive.
Full peripheral recognition therefore matches this entire outer
block to the source block in Section2, with every required parity.

Next comes the first simple full-S4 arrival from the A4 state. PF89
forces its winding to be EXACTLY2; its drilled root is an annulus
with one genuine cone2. It is the winding2 node of Section3.
At least one MORE canonical genuine node lies inside it: the earlier
nonsolid reset terminal disk. If nonsolid cyclic3 loops occur after
that reset, they supply still more nodes between the reset and this
arrival. The reset disk has negative Euler base and incompressible
boundary; all actual fibre determinants are nonzero. No product or
fibre merger deletes this part of the path (PF83/89).

This is precisely the forbidden source shape of Section3.
Contradiction. It excludes BOTH earlier-reset terminal assignments
(2,8) and(4,10), ALL four actual(r,J) assignments, arbitrary all-plus
pre-reset solid histories, ZERO or arbitrarily many post-reset cyclic3
loops, and EVERY positive finite later-tail length, all at once.
PF95's finite first-post-reset winding lists are superseded by this
closure. No source packet is assumed at any target prefix.

## 5. Complete genuine iterated-quartic exhaustion

Use PF82's complete ten-state reduction and its exact absorbing
simple full-S4 state3. A passing trajectory either starts there or
first enters it from states2,4,5,6,7,9. The exclusions are now:

| Initial or incoming class | All-length exclusions in BOTH pencils |
| --- | --- |
| Initial simple full-S4 state3 | PF88 |
| Incoming S4 states4/5 | PF86 |
| Incoming cyclic states6/9 | PF85 |
| Incoming S3 state2 | PF90/91/93/95 |
| Incoming A4 state7, still solid at arrival | PF87/89/95 |
| Incoming A4 state7, earlier nonsolid reset | PF87/89/94 and Section4 here |

All later passing stages are actual odd-winding cyclic2 tails.
The complete state graph and these all-length proofs exhaust
EVERY genuine positive iterated-torus quartic end match in the
retained scope, in COMMON and COMPLEMENTARY source pencils.
There is no remaining genuine iterated A4 reset family in that scope.

The conclusion is still conditional on the mathematical hypotheses
already retained throughout: the whole smooth normalization and
entire-line charts (exactly two omitted lines in the complementary
comparison), exact fold volume, actual preferred/peripheral maps,
genuine positive initial torus and cable stages, and simple-meridian
full-final-S4 action. These are not asserted for all degenerate
quartic boundary data or arbitrary higher-degree normalizations.

## 6. Controls and essential limits

The [checker](verify_prefix_excess_closure.py),
[receipt](prefix_excess_closure_receipt.json), and
[manifest](sources/prefix-excess-closure/manifest.json) retain the
unbounded algebraic proofs, exact local prefix chains, predecessor
seals and actual reset/arrival/tail covers. Primary characteristic,
contact, JSJ and full-peripheral inputs are reused with their
existing audits, not independently reproved or newly assumed.

Countercontrols make the scope concrete:

* The FORMAL update(8,2) --a2,kappa8--> (2,8) would evade excess
  absorption, but violates index order. It is not a valid prefix
  extension; determinant recurrence alone is insufficient.
* The actual local chart sigma28,p140,support{140,48,46,35} has
  delta2 and terminal(2,8). Its virtual odd block has determinant11,
  free order6 and an ODD bridge23, so it fails Section2's required
  parity. The chart is retained; it is not declared nonexistent.
* Actual local charts with a balanced odd block and a LAST
  winding2/free4 node, determinant2, exist. They have no node
  inside the winding2 arrival, so Section3 does not exclude them.
  PF95 uses separate target slope data to close still-solid arrivals.
* A primitive valid extension of a determinant2 prefix can exist
  as a LOCAL chart with NEGATIVE determinant. Such a chart cannot
  be a prefix of the required full positive-determinant source.

Finite controls do not prove universality. No full higher-degree
packet theorem, degenerate-stage classification, full quartic/Keller
exclusion, formal verification, external review or novelty claim is
made. KIG28's full generic A4 theorem is not used to identify local
or intermediate A4 actions. The full research goal remains ACTIVE:
next establish which remaining actual quartic boundary data can
escape the genuine-stage hypotheses, and retain the full compatible
finite-map/completion and higher-degree obligations.
