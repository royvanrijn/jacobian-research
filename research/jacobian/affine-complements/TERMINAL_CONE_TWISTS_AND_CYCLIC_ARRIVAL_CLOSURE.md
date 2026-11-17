# Terminal cone twists and complete cyclic first-arrival exclusion

15 November 2026. PF2RULED85 closes PF84's remaining FINAL cyclic4
arrival of winding3. Consequently ALL first arrivals at the simple
full-S4 state from CYCLIC order3/4 companion images are excluded in BOTH
source pencils, at arbitrary companion depth and arbitrary later
tail length, under the actual quartic packet hypotheses. In PF82
state labels these incoming states are6(C3) and9(C4).

The new mechanism is a degree-free comparison of NORMALIZED cone
coefficients, not only their orders. The last source winding gives
a coefficient inverse to kappa-1; an unbranched fibre cover of a
positive torus-knot exterior gives an inverse to f*kappa.
Combined with the all-depth terminal modulus, the two allowed signs
force a winding divisor. At f=4,n_R=13 the divisors are8 or16,
impossible for an odd winding greater than1. The repeated-order25
case separately fails4=+/-1 mod25. No deeper coefficient search is used.

Other noncyclic incoming states2/4/5/7, and initial state3 trajectories
outside the previously closed seed classes, remain open. This is not
a full quartic or higher-degree exclusion.

## 1. Terminal residue pair from the source extraction

Retain PF77-84's full smooth exact complementary source with TWO
entire omitted lines and their actual whole charts. Let its terminal
D disk have winding order a>=2, free-tail order kappa>=2 and at
least two characteristic nodes. Then a divides sigma,p_D and all
support indices before the last gcd drop. The final index is

    j_*=sigma-1+kappa, gcd(a,j_*)=1,
    j_*=kappa-1 mod a.                                 (1)

Let a negative plumbing leg have order alpha, and b be the positive
cofactor obtained by deleting its vertex adjacent to the Seifert node.
Its continued fraction is alpha/b. In an oriented cone presentation
the normalized coefficient is b or -b modulo alpha, according to
the boundary convention. We work with the UNORIENTED necessary
class {b,-b}. Allowing each leg its own sign only WEAKENS the true
matching requirement, so a contradiction here excludes either
orientation convention and any actual oriented identification.

The two terminal legs have the following residues:

    order a:      b_a = (kappa-1)^(-1) mod a, up to sign;
    order kappa:  b_kappa = kappa-1, hence +/-1 mod kappa. (2)

Here is the all-depth local derivation. Resolve the preceding
characteristic terms first. Their exponents are multiples of a.
At the final characteristic pair the smooth local projection has
order a and the transverse characteristic order m differs from
+/-j_* by a multiple of a: previous monomial divisions and subtraction
of truncated terms change it only by preceding, divisible orders.
On the boundary of the order-a fibred solid torus its primitive
regular fibre is m*mu+a*ell. A primitive section c*mu+d*ell satisfies
a*c-m*d=+/-1; after filling mu its cone relation gives normalized
coefficient +/-d. Thus m*d=+/-1 mod a, proving the first part of(2).

The contact between conjugate branches at the last drop is2p_D-j_*.
The last free creation of D has contact2p_D+1-sigma; their difference
is kappa. The remaining smooth free blowups leave the terminal
chain of kappa-1 vertices of weight-2. Its determinant is kappa and
its adjacent cofactor kappa-1. This proves the second part, at any
earlier depth. Extra noncharacteristic terms below j_* cannot change
the last gcd or these residue classes.

This can equivalently be read from the determinant-one Euclidean
matrices [[c,-1],[1,0]]: the head and tail cofactors multiply to1
modulo the leg determinant. The parameter supplied at the opposite
end is j_*=kappa-1 mod a. The local argument uses no bounded extraction.

The normalization is checked against complete Hatcher PDF25-27,
especially Proposition2.1: with boundary, changing a section changes
a cone coefficient by an INTEGER, not by an arbitrary unit. Reversing
orientation changes its sign. The negative-leg sign is checked against
Neumann--Raymond Theorem5.1, complete PDF11-12. We allow both signs
rather than infer an orientation from homology.

## 2. Cone coefficients of an actual fibre cover

The positive torus-knot exterior T(a,kappa), gcd(a,kappa)=1, has
normalized coefficient at its order-a cone

    +/-kappa^(-1) mod a.                               (3)

Use the weighted circle action: on a small exceptional solid torus
the regular fibre has slope kappa*mu+a*ell. A determinant-one
section gives(3) by the same calculation as above.

An unbranched cover over that cone with base degree1 and generic
fibre degree f sends h_new to h_old^f and a base section to the
old section times an integral fibre power. Thus

    f*b_new=b_old mod a.                               (4)

For a base-degree1 fibre cover, when gcd(f,a*kappa)=1, its two cone
coefficients are consequently

    order a:      +/-(f*kappa)^(-1) mod a,
    order kappa:  +/-(f*a)^(-1) mod kappa.               (5)

This is a LOCAL primitive-torus calculation; it is valid at every
cover degree. Neumann--Raymond Lemma1.3 on complete PDF4 independently
checks the multiplication in the fibre quotient; solving it in
the opposite direction gives(4). The cover direction must not be
reversed. Neither cone order alone nor an Alexander determinant
determines these residues.

All target terminal pieces below have negative orbifold Euler
characteristic. Their canonical fibration is therefore the one
recognized by PF79, with its product/Klein exceptions already
excluded. Full peripheral recognition, not a homology isomorphism,
compares these actual terminal Seifert pieces.

## 3. A degree-free terminal divisor obstruction

Suppose the full complementary source has determinant delta=1 and
ramified diagonal n_R, and its actual terminal piece matches the
base-degree1 fibre-f cover in(5). PF84's all-depth terminal formula
gives

    n_D=a^2*Z+a-(a-1)*kappa,
    1=delta=n_R*kappa mod a,                           (6)

since a divides p_D and hence the full cross term q_src*p_D.
Matching the order-a coefficient in(2) with(5), even allowing either
sign, forces

    f*kappa=+(kappa-1) OR -(kappa-1) mod a.

Combine with(6) and gcd(a,kappa)=1:

    a divides n_R+f-1 OR n_R-f-1.                      (7)

All source orders, characteristic depths and fibre-cover degrees
are unrestricted under the displayed hypotheses. If either integer
in(7) vanishes, that alternative gives no obstruction. This is a
conditional source/terminal-cover lemma, not a classification of all
higher-degree normalization packets.

## 4. Close the both-odd cyclic4 residual at EVERY depth

PF84 reduced every remaining cyclic4 arrival to a FINAL winding3.
Its companion has an arbitrary nonempty chain of odd cyclic4 cable
levels over T(p,q). In the both-odd terminal case the actual terminal
cover has base degree1, generic fibre degree4, cones(p,q), and no
branched homology. Source recognition forces rho=4,q_src=3,n_R=13;
the closed target has H1=0, hence delta=1.

Let a be whichever of p,q is the source terminal winding, and
kappa the other. The source D side has at least two nodes because
the companion has a genuine cable, so(1)-(2) and(6) apply.
Use f=4,n_R=13 in(7):

    a divides16 OR8.

But a is ODD and at least3. Impossible. This closes BOTH possible
terminal assignments, and allows both orientation signs. It does
not need PF84's smaller-winding or finite-divisor pruning.

## 5. Close the repeated T(2,25) residual

PF84 reduced the even-first-pair terminal to T(2,25), whose actual
connected cyclic4 exterior cover has base degree2 and fibre degree2.
There are TWO order25 cones. At each, the target torus-knot residue
is +/-2^(-1), and the local base cover has degree1 at each preimage.
Equation(4) gives target residues

    +/-4^(-1) mod25.

The source terminal disk has orders(a,kappa)=(25,25), so its free-tail
cofactor has residue +/-1 mod25 by(2). It cannot match either target
cone: the necessary condition would be

    4=+1 OR -1 mod25,

which is false. The two identical target orders cause no ambiguity;
both carry the same residue class up to sign. This also excludes
orientation reversal without choosing a global section.

Thus every PF84 final cyclic4 winding3 residual is excluded in
the complementary source. PF84 already excludes all common sources,
all cyclic3 arrivals, cyclic4 arrivals with later tails and larger
arrival windings. The torus-companion two-stage case is PF81.
Together these close the ENTIRE incoming CYCLIC first-arrival class
at all actual companion depths and tail lengths, with arbitrary
lower-weight key terms and no prefix eta restriction.

## 6. Scope, evidence and remaining problem

The target statement is about actual iterated torus-knot quartic
packets of PF46/50, full simple S4 final periphery, all true fixed-sheet
and ramified meridian fills, and either the common source or the
explicit smooth/exact/whole-chart two-line complementary source.
The all-degree lemmas concern actual terminal pieces under their
stated last-gcd, homology and local-cover hypotheses. Smoothness
and these hypotheses are not automatic for arbitrary higher degree.

The [checker](verify_terminal_cone_twists.py),
[receipt](terminal_cone_twists_receipt.json) and
[manifest](sources/terminal-cone-twists/manifest.json) retain full
Euclidean-cluster leg matrices including a non-self-inverse order5
cofactor, determinant-one section/cover identities, both sign branches,
and the sealed all-length PF84 reduction. The written local-slope,
free-tail and congruence arguments prove the unbounded statements;
finite controls are not an extrapolation in degree or depth.

Fresh complete Hatcher PDF25/26/27 and Neumann--Raymond PDF3/4/11/12
were rendered and visually inspected. Their classification/cover/
plumbing results are imported with the explicit conventions above.
PF77/79's complete fibre-uniqueness, JSJ and full-peripheral Waldhausen
audits remain retained, not independently reproved.

PRIORITY is now first arrival through NONCYCLIC incoming states2/4/5/7,
plus trajectories already at state3 on the initial torus outside the
closed intrinsic seed ranges. Keep the two-noncyclic-level bound and
actual cyclic annulus geometry. Do not reopen cyclic incoming states6/9,
or delete their earlier geometric levels. Full longer-stage/quartic/
higher-degree exclusion and the full goal remain OPEN. No finite
Keller map, formal verification, external review or novelty/priority
claim is asserted.

