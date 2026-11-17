# Torus-cusp excess from a solid-torus source complement

PF2RULED114, 16 November 2026. This uses the knot-theory route at a FINITE
cusp, with its actual meridian fillings. It gives a proportional, all-degree
ramification bound; no iteration-depth or permutation census is involved.
All global statements retain the ORIGINAL finite normal Keller completion,
the WHOLE marked U=A2, and ALL omitted primes. Singular N is allowed.

**Local theorem.** Let a connected finite normal surface germ h have just one
reduced critical branch D. Assume the source pair is a smooth plane with a
smooth D, or a small cyclic quotient with D the image of a smooth invariant
coordinate axis. Suppose its image is a torus cusp C of coprime type (a,b),
2<=a<b. Write e for transverse ramification, q for local degree D->C^nu,
ell for local surface degree, w=e*q, and eta=ell-w. Then

    w <= 5*eta-2*delta <= 5*eta-2*q,               (1)

where delta>=1 is the order of regular-fibre monodromy on the ell sheets.
In particular eta=1 forces w<=3. No scalar-Morse or homogeneous leading-form
hypothesis is imposed on h, and q need not be one. The TARGET cusp is a
one-Puiseux-pair cusp; an arbitrary plane branch is not covered.

**Global theorem.** For an actual omitted prime of index e with 2e>d,

    e+ceil((e+2)/5)+r <= d-1,                    (2)

where r counts ALL omitted primes. Consequently EVERY omitted index obeys

    e <= floor((5*d-17)/6) for d>=9.              (3)

The entire former e=d-4/d>=9 two-ramified-prime profile of PF113 is
EXCLUDED, including singular sources. The new frontier is (2), the
shared-image bound below, other image arrangements, and actual finite-target
algebra. It is not replaced by a proposed census of the next few indices.
The retained d=7,8 bounds remain e<=3,4 respectively; no inverse theorem
in those degrees or in general is claimed.

## 1. The local source complement is a solid torus

Take the link pair (W,K_D) of the source germ. In the smooth case K_D is
the unknot in S3. In the cyclic case, straighten the smooth invariant lift
of D and diagonalize the small cyclic action, as in PF106 Section1. Its
link pair is the diagonal quotient of (S3, one coordinate-axis circle).
Deleting a tubular neighborhood of that circle gives a solid torus; its
free cyclic diagonal quotient is again a solid torus. Thus

    W minus int(nu K_D) is a solid torus.         (4)

This is a LOCAL link statement. It does not identify the end of the whole
affine N from the type of one of its singularities. The cyclic-axis input
is the retained actual-boundary result of Borisov/PF105/PF106; it is NOT
being assumed for arbitrary curves on arbitrary normal surface germs.

Use the standard local-link construction for a finite complex analytic
map: a sufficiently small transverse target link, pulled back, is a
branched cover representing the source link pair. This follows from local
conical structure and link independence for analytic pairs. A smooth
invariant axis in a cyclic quotient has the pair just described regardless
of the particular finite defining map. Removing the complete inverse image
of C gives a connected ell-sheeted unbranched cover E of its knot exterior.
Connectedness is irreducibility of the local normal germ off a divisor.

The meridian permutation has q cycles of length e and eta fixed sheets.
The q ramified cycles form ONE JOINT meridian-longitude orbit: the link of
the single branch D maps to the cusp knot with longitudinal degree q.
The other sheets represent ordinary unramified preimage curves, NOT
additional omitted primes. Filling just their boundary tori by the actual
lifted meridians recovers the solid torus (4), leaving K_D drilled.

Other global branch components may meet C at the cusp image. They have
no critical preimage in this local factor, because actual omitted primes
are disjoint. The local factor's branch is C alone. No other global prime
is removed from the global ledger in using this local statement.

## 2. Keep regular-fibre monodromy and the unramified filling slopes

The exterior of T(a,b) has orientable base orbifold a disk with cone orders
a,b. In positively oriented meridian-longitude notation its regular fibre
is h=mu^(a*b)*lambda. Pull its Seifert fibration back to E. Brin,
Lemma2.4.21 and the complete discussion on PDF145-147, supplies this
finite-cover fact; it is also the retained PF69 fibre-orbit construction.

The image H of h is central in the TRANSITIVE local sheet action. Every
cycle of H has the same length delta: a power fixing one sheet fixes all
sheets by centrality and transitivity. Therefore

    q divides delta, delta divides w, delta divides eta,
    base-cover degree = ell/delta.               (5)

The q meridian cycles of D are permuted as a q-cycle by lambda, hence
also by H=mu^(a*b)*lambda; this proves q divides delta. The divisibility
delta|w uses the single ramified peripheral orbit, which is H-invariant
and has w sheets. On the eta meridian-fixed sheets,
lambda=H, so there are exactly eta/delta unramified boundary tori, each
with lattice basis

    m=lift(mu), L=lift(lambda^delta).

On such a torus the primitive lifted regular fibre is

    lift(h^delta)=a*b*delta*m+L.                  (6)

Its intersection with the actual filling meridian m is ONE, even when
delta>1 and longitude joins ordinary sheets. Hence each unramified
filling extends the Seifert fibration with an ORDINARY core, adding no
exceptional fibre. This filling assertion follows directly from the
primitive slope (6), or Brin Lemma1.3.7/Theorem1.3.8, PDF26-27.
It is not a consequence of unbranched-cover theory alone.

The filled manifold (4) is a solid torus. Brin Lemma2.1.5 and its full
proof, PDF86-87, imply its inherited base is a disk with at most ONE
exceptional fibre. Alternatively its orbifold group is a quotient of Z;
an orientable base with boundary has orbifold group a free product of
handle generators and finite cone groups, so it can have genus zero and
at most one nontrivial cone group. The base of E has the SAME cone points
and genus, with eta/delta extra boundary circles. Let alpha be the order
of its sole exceptional fibre, setting alpha=1 when there is none.
Every inherited cone order divides a or b (look at the disk cover over
one target cone value). In particular

    alpha in {1} union divisors(a) union divisors(b), alpha<=b.

The argument retains nonidentity regular fibre H, the actual ordinary
meridian fillings, and rational/quotient source topology. It makes no
inference from end H1 or an uncolored knot alone.

## 3. Orbifold Euler characteristic pays for the cusp excess

The base of E has 1+eta/delta boundary circles and at most one cone of
order alpha. Multiplicativity of orbifold Euler characteristic gives

    (ell/delta)*(1/a+1/b-1) = -eta/delta+1/alpha,

or the exact necessary identity

    (w+eta)*(1-1/a-1/b)=eta-delta/alpha.          (7)

For completeness, remove the two cone disks in the target. The ordinary
base cover multiplies Euler characteristic by ell/delta. Restoring a disk
for each permutation cycle over a cone of order a contributes that
cycle's local cone weight (cycle length)/a, and similarly for b. This
proves the weighted formula without assuming an unramified cover of
the underlying disks. The distinction between fibre degree delta and
base degree ell/delta is essential.

Solving (7) gives

    (a+b)*eta=(a*b-a-b)*w+a*b*delta/alpha.        (8)

Here w>=delta. Define

    A=5*a*b-6*a-6*b
     =(5*a-6)*(b-3)+9*(a-2)>=0,
    C=5*a*b-3*a-8*b
     =(5*a-8)*(b-3)+12*(a-2)>=0.

Multiplying 5*eta-w-2*delta by a+b and using (8) gives EXACTLY

    A*(w-delta)+C*delta+5*a*(b/alpha-1)*delta.    (9)

Every term is nonnegative, since a>=2,b>=3,alpha<=b. This proves (1)
at ALL indices, local degrees, quotient orders and residue degrees.
For the type (2,3) with alpha=3 the sharp arithmetic relation is
w=5*eta-2*delta. Arithmetic equality alone does not construct a finite
analytic germ, and no optimality/realizability claim is made.

## 4. Apply the local theorem to high-index images

PF112 shows that 2e>d makes the image C_R injectively normalized by A1,
hence topologically contractible, and gives f_R=1. Chau's component
Remark excludes a smooth A1 image. The already retained Palka TheoremA
therefore supplies ONE global polynomial target coordinate system in
which C_R is x^b=y^a, coprime 2<=a<b. It has a single torus cusp. This
classification is IMPORTED, not independently reproved here.

At the unique point p of R over that cusp, the retained actual source
pair is smooth-axis or cyclic-axis. Its local residue q is one, so (1) gives

    eta_p >= ceil((e+2)/5).                      (10)

PF113 also supplies another RAMIFIED prime with a DIFFERENT image.
In the original all-prime Euler ledger, the other r-1 primes cost at
least r (one is ramified); keep all other nonnegative excess terms.
Then d-1>=e+eta_p+r, proving (2). Equivalently

    6*e <= 5*d-7-5*r <= 5*d-17.                 (11)

If 2e<=d, simply use e<=floor(d/2). For d>=9 this is at most
floor((5*d-17)/6), proving (3) for EVERY prime. Odd/even integer
endpoints are checked separately; strict half is not assumed at equality.

PF113's former e=d-4/d>=9 profile forced eta_R=1 and e>=5. Equation(1)
instead requires e<=3, so that ENTIRE class is now excluded. The sealed
PF113 necessary-profile proof remains valid as a historical predecessor.
Its former OPEN label must not be kept as a current queue.

This is a proportional bound, approaching 5d/6, rather than a new
constant offset obtained by expanding one more near-maximal strip.
No claim that every image component is contractible is made: that step
uses the strict high-index/high-mass condition.

## 5. Shared-image mass, including nontrivial residues

Suppose C carries t actual omitted primes D_i with (e_i,f_i), mass
w_C=sum e_i*f_i, and 2*w_C>d. PF112/PF113 make C contractible, exclude
unramified primes over it, and require a different RAMIFIED image.
Chau and Palka again give its unique torus cusp. For each i, the finite
normalized map D_i^nu->C^nu has, over the cusp parameter, points p_ij
of local degrees q_ij, with sum_j q_ij=f_i. Apply (1) in EACH local
factor, then add. Each delta_ij>=1, and every D_i has at least one such
point. Thus the excess charged to this image satisfies

    eta_C >= sum_i ceil(f_i*(e_i+2)/5)
          >= max(t,ceil((w_C+2*F_C)/5)),          (12)

where F_C=sum_i f_i. Here delta_ij>=q_ij by (5), so adding local costs
charges the full residue sum, even when several cusp preimages occur.

All excesses away from these cusp points remain nonnegative. This does
not assume a nontrivial residue map has several cusp preimages or that
q_ij is its global degree; it only uses the sum of local degrees.
Writing h_C for the first sum in (12), the complete ledger requires

    w_C+h_C+r-t <= d-2,
    w_C+h_C <= d-3,
    6*w_C <= 5*d-15-2*F_C.                       (13)

The second/third use r>=t+1. Known same-image sublists still supply
necessary filters, retaining rather than deleting unknown primes.
PF113's irreducible FULL branch half-mass bound persists. The new
estimate does not force all branch unions to be weighted or all
residue degrees to be one.

## 6. Genuine finite controls and verification boundary

The local theorem is not a prohibition on finite torus-cusp covers.

* For every e>=2,m>=2, the actual finite smooth map
  P=u, Q=integral_0^v (u-z^m)^(e-1) dz has degree
  ell=m*(e-1)+1, sole smooth critical branch u=v^m and generic index e.
  Its cusp is (m,ell), q=1, delta=1, alpha=m. It satisfies (7);
  eta=(m-1)*(e-1). Its etale open is C* times C, not a whole Keller plane.
  The cubic cusp e=m=2 has eta=1 and remains ALLOWED locally.
* PF105's genuine finite A1 quotient family, s>=1, has
  e=2s+1,q=s,ell=4s^2,eta=2s^2-s, cusp(2,4s-1),
  delta=s,alpha=4s-1. The regular-fibre lift has order s, so this is
  an essential longitude-joining control. At s=1 its e3/eta1 germ
  remains allowed. At s>=2 its residue>1 remains allowed. Its actual
  etale open is C* times C. The construction and finiteness are the
  retained PF105 proof, not inferred from these numerical parameters.
* PF107's actual degree12/e9 finite A1 quotient has cusp(2,3), q1,
  eta3, delta1,alpha1. It also satisfies (7) and is not excluded as
  a local finite germ. It lacks the required whole-plane Keller open.

The checker pins these predecessors and primary sources; checks (7)-(9),
the all-parameter finite-family identities, degree/parity budget algebra,
and selected exact peripheral permutation controls. Those controls check
implementation and signs; the universal theorem rests on Sections1-5
and the declared topological/classification inputs. A necessary arithmetic
gate is not a cover realization certificate. No formal verification,
external review, historical novelty or full JC2 solution is claimed.

The global Laurent finite-algebra implication, noncontractible image
arrangements and compatibility of all local covers remain OPEN. Published
Domrina2000 quartic retirement and the shared trace priorities persist.
The complete user goal remains ACTIVE; no external outreach is authorized.
