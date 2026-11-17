# An intrinsic gate for one-boundary normalizations

14 November 2026. PF2RULED36 separates the general valuation statement
from the extra geometry needed to apply the marked-end obstruction.
No identification with the Wright ring is assumed in Sections 1--3.

The algebra statement holds in every finite field degree. The geometric
statements concern a finite normalization of a plane Keller map of field
degree four whose complement of the source plane is **one omitted prime
of generic index two**. They do not cover an arbitrary multi-prime packet,
or silently dispose of a higher-index omitted prime.

For that quartic setting the preferred longitude has identity monodromy.
If its end is a lens space with a ramification core, its homology order
cannot be odd and divisible by three. In particular the normalization
cannot admit a source-coordinate ruling with multiplicity-one omitted
fiber. Multiple fibers and other end pairs remain genuine open gates.

The topology below extends the written PF2RULED35 argument; neither that
proof nor this extension has independent external review or formal
verification. Exact controls check identities, finite groups and source
seals, not the entire three-manifold argument.

## 1. The associated-prime condition follows from one omitted prime

Let A be a normal affine domain, with an open immersion
U=Spec B -> N=Spec A, where B=C[x,y], Frac A=Frac B and N minus U
is precisely one prime divisor R. Let p be its prime ideal and v its
normalized divisorial valuation. The Krull intersection description gives

    A = B intersect O_(N,eta_R) = {b in B : v(b)>=0}.       (1)

Indeed every other codimension-one point lies in U, so regularity on U
already checks every other valuation in the normal-domain criterion.

For b in B minus A put n=-v(b)>0. Then

    Ann_A(b+A) = (A:b) = {a in A : v(a)>=n}
                         = p^(n),                       (2)

where the last expression is a symbolic power, not an assertion about
ordinary powers. This ideal is p-primary: all valuations of elements of
A are nonnegative, and multiplication by an element outside p leaves v
unchanged. Therefore every nonzero cyclic submodule of M=B/A has sole
associated prime p. Associated primes of a module are detected in cyclic
submodules; p-primary ideals in a Noetherian ring have associated prime
p. Consequently

    Ass_A(B/A)={p},   Supp_A(B/A)=V(p).                    (3)

Here B/A is nonzero: otherwise the given open immersion would be an
isomorphism. Multiplication by a outside p is injective on M;
multiplication by a in p is locally nilpotent on M. The latter means
that its exponent may depend on the element of M.

In fact Ann_A(M)=0. Choose b of negative valuation; powers b^j have
unbounded negative valuation. No nonzero fixed a in A can put all ab^j
back in A. Thus a uniform annihilating power must not be substituted
for local nilpotence. No Keller equation or degree cutoff was used.

This proves directly the associated-prime formulation discussed in
Rodriguez Diaz Lemma2.3. It supplies no ruling multiplicity.

## 2. The quartic longitude is intrinsic

Now assume the actual finite normalization pi:N->A2 of a polynomial
Keller map F:U=A2->A2 has degree four, and N minus U=R is one prime
with generic ramification index two.

[PF2OMIT7](../jc2-omitted-divisors/SIMPLE_RAMIFICATION_AND_RADIAL.md)
makes N smooth along the whole R and gives residue degree one.
[PF2BATT1](../jc2-boundary-attachment/PROOF.md) then gives R=A1.
N is smooth everywhere, since U is smooth. Write t for its boundary
parameter and C=pi(R). This is the entire branch/nonproper curve;
gamma:R->C is its finite birational normalization.

The [Euler ledger](../jc2-omitted-divisors/GLOBAL_LEDGER.md) gives
4-1=2+E_R, hence E_R=1. The local index-two normal form in PF2OMIT7
therefore supplies exactly one ordinary (2,3) cusp. Every other boundary
parameter is immersive. A collision uses two local blocks of length
two, so has exactly two immersed branches; the cusp cannot collide
with another parameter. Write (a_i,b_i) for each pair and m_i>=1 for
its intersection multiplicity. No transversality is assumed.

Consider the reduced residual divisor T in

    pi^*C=2R+T.

Its normalized finite cover of the normalization of C has degree two.
At a collision, the fold chart on one ramified sheet pulls the other
branch back as z^2=f(t), with ord f=m_i. Normalization branches this
quadratic cover at a_i exactly when m_i is odd, and likewise at b_i.
At the ordinary cusp the local cubic block contributes one residual
branch of degree one and the fourth, etale sheet contributes the other.
Thus both points above its normalization parameter are unramified.
There is no other finite branching.

A degree-two etale algebra over C(t) is determined by its odd
valuations. Since C is algebraically closed, there is no additional
unramified quadratic factor over A1. The complete cover is therefore

    w^2 = S(t) = product_(m_i odd) (t-a_i)(t-b_i),          (4)

with the split algebra when there are no odd m_i. Its squarefree
degree is even. Its two infinity points are unramified, also in the
split case. This part of the residual passport uses no Wright chart,
normal coordinate or conormal coefficient identity.

The global monodromy is S4: it is transitive and generated by the
simple meridians of C, all transpositions. The infinity group maps
onto the global group by the primary input retained in
[PF2RULED34](TORUS_INFINITY_OBSTRUCTION.md).
For a meridian mu and preferred longitude lambda of the knot at
infinity, rho(mu)=tau is a transposition. The residual infinity
calculation says lambda fixes the other two sheets. It commutes with
tau, hence is 1 or tau; being a knot longitude it maps into A4.
Therefore

    rho(lambda)=1.                                      (5)

A change of push-off framing multiplies by a power of tau, fixing
both residual sheets. Thus the residual conclusion applies to the
preferred longitude without identifying it with a normal framing.

## 3. The cabling obstruction extends to odd orders divisible by three

Here is a purely topological extension. Let n>1 be odd and divisible
by three. No iterated torus knot K admits a connected simple
four-sheeted cover with S4 monodromy, identity preferred longitude,
source a lens space with H1=Z/n, and ramification knot a core
(its exterior is a solid torus). A core represents a primitive
generator of Z/n, a stronger condition than merely nonzero.

Use the complete cabling case analysis of
[PF2RULED35 Sections 2--5](QUARTIC_CABLING_EXCLUSION.md), with its
actual boundary slopes, lattice identifications and compression
argument. That analysis depends on the order three only in the
places audited here. All incompressible-torus, graph-cycle,
nonfiber filling and noncyclic free-product contradictions continue
to apply in every lens space with an essential marked knot.

* Auxiliary A4 type, torus base: its H1 order is 2uv, even, whereas n
  is odd. The unknot group is cyclic.
* Auxiliary A4 type, nontrivial fiber: after the forced solid-torus
  filling the closed Seifert base has coprime cone orders p,3.
  The regular fiber generates H1, but equals three times the marked
  class. Multiplication by three cannot generate Z/n when 3 divides n.
  This replaces the stronger n=3 statement that this class is zero.
  The other A4 cable cases use only incompressibility or at least
  three genuine cone orders.
* Symmetric type, four-letter p=3 case: its marked regular fiber
  is three times the solid-torus core, so cannot be primitive in
  Z/n. This replaces the n=3 null-class argument. The p=2 and
  annular cases are unchanged.
* Three-letter transposition pattern: the p=2 boundary lattice
  m0=3mu_J, l0=lambda_J still identifies the actual meridian of the
  solid torus B. Replacing it by the usual triple branched filling
  preserves the lens/core condition and gives the excluded A4 type.
  This uses the precise gluing slope, not equal homology alone.
  For p>=4 the marked-exterior arguments are unchanged.
* Three-letter three-cycle pattern: the companion dual-graph cycle
  or incompressible marked-exterior torus is independent of n.
* Two-letter pattern: p is odd. Even q is excluded as before.
  For odd q, the marked branch has index two. Fiber filling would
  kill twice its primitive class; that class is nonzero since n is
  odd. The remaining filling has distance one. The coprime cones
  p,2 give torsion-free H1(B)=Z and the exact boundary map
  m0=2mu_J ->0, l0=lambda_J ->generator. Replacement by the standard
  simple branch solid torus preserves H1=Z/n. The new source is a
  lens space with the companion ramification core, and retains S4
  and identity longitude. The class identity [K_R]=+/-p[l0] forces
  gcd(p,n)=1 since [K_R] is primitive. Thus gcd(p,2n)=1; in particular
  the removed winding remains coprime to six.

Each surviving simple step strictly lowers finite cabling depth.
The A4 alternative has already been excluded. The terminal simple
torus case is excluded by PF2RULED34's finite-H1 torus lemma: its
ramification class is zero. The unknot cannot have S4 monodromy.
This proves the extension for every n=3,9,15,21,..., with no bound
on winding, signed second cable parameter, tangency or cabling depth.

Combining this with Section 2 excludes those lens/core ends for the
actual one-simple-prime quartic normalization. It does not prove
that an arbitrary such normalization has a lens/core end.

## 4. A multiplicity-one ruling would already suffice

The relevant primary geometric input is Borisov--Gabber--Vasiu,
Theorem11.2(1)--(2): a normal affine surface containing A2 with an
irreducible complement admits an A1 ruling to P1; its omitted fiber
may be qR, and q=1 gives local affine-line triviality.
Lemma11.1 extends coordinate projections to P1.
[Author manuscript](https://www.ihes.fr/~/gabber/BGV160.pdf).
These are explicit imports, not independently reproved here.

Suppose the actual normalization in Section 2 has a source-coordinate
ruling with fiber infinity equal to R scheme-theoretically. It is an
affine-line torsor over P1 with linear bundle O(-m). In coordinates
x,y on the source plane its second trivialization has the form

    u=1/x,   v=x^m*y+g(x),

with g Laurent on the overlap; affine terms do not affect the
Jacobian. Hence

    dx wedge dy = -u^(m-2) du wedge dv.

The Keller equation and sole index-two ramification give the divisor
of this rational two-form as R. Thus m-2=1 and m=3. A smooth section
identifies the torsor with O(-3), preserving the whole fiber R.
Its end is a lens space of order three with fiber a core. The stable
end-pair group and lens/core recognition argument in PF2RULED35
applies verbatim. Sections 2--3 now give a contradiction.

Consequently the omitted fiber of every BGV ruling has q>=2 in this
quartic setting. More specifically, no source polynomial coordinate
can have v_R(x)=-1: its extended projection would have precisely the
multiplicity-one fiber R at infinity. Coordinates with nonnegative
valuation are not excluded by this assertion.

The canonical calculation m=3 also holds in every field degree with
one index-two omitted prime and a multiplicity-one ruling. In that
setting the marked end and PF2RULED34 exclude torus knots at infinity
in every degree. The stronger finite cabling descent here is quartic.

No explicit isomorphism with the Wright ring is needed for this
conditional exclusion. Establishing q=1 for an actual normalization,
or another excluded end condition, is sufficient; neither follows
just from R being a smooth reduced affine line.

## 5. Audit of the multiplicity gap, and a retained control

Rodriguez Diaz, arXiv2403.02219v2, Theorem3.5 passes from the reduced
complement to the scheme-theoretic infinity fiber in its bundle
argument. The cited fibration result permits multiplicity. We do not
adopt that passage, or the whole claimed completion theorem.
[Primary manuscript](https://arxiv.org/html/2403.02219v2).

The issue can be seen within the actual Wright surface itself, with
smooth N, R=A1 and K_N represented by R. On its source plane put

    a=x^2+y,   d=x+a^2.
    x=d-a^2,  y=a-(d-a^2)^2,   J_(x,y)(a,d)=-1.

In the entire omitted chart x=1/u, y=-u^2+u^3 v,
the two coordinate functions a,d have pole orders two and four.
Their numerators after multiplication by u^2,u^4 restrict to one
on R. Both projections therefore extend to P1 with infinity fibers
2R and 4R, despite having the same reduced complement R.
This is a control against the unsupported geometric inference, not
a counterexample satisfying finite Keller-extension hypotheses.
The polynomial coordinate pair itself does not extend to regular
target functions on N.

Section 1 means the associated-prime formulation alone adds no
missing multiplicity constraint in the one-prime case. An argument
using the full finite Keller algebra could still supply one; it
has not been given here.

## Scope and next gate

Proved here: the all-degree valuation/primary-module lemma, the
intrinsic quartic residual-infinity longitude, the odd-three-divisible
lens/core obstruction, and exclusion of multiplicity-one rulings under
the stated one-simple-prime quartic hypotheses.

Still open: multiple omitted fibers, general non-lens or non-core
ends, higher-index or multi-prime quartic packets, and higher finite
field degrees on Wright. The full research goal remains active.
No new coefficient search, full packet exclusion, general comparison,
formal verification, independent full reproof or novelty priority is
claimed. Sources are pinned in
[sources/intrinsic-normalization/manifest.json](sources/intrinsic-normalization/manifest.json).
