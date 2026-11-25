# Two-cable augmentation ranks and the odd-cable compatibility obstruction

**Theorem TwoCableGapBarrier.** Let K be a member of either the four-bridge
family K_n (n>=0) or the uniform family K_(b,k) (b>=5,k>=0), with the
explicit braids and colors defined below. For every odd integer Q, let
C_(2,Q)(K) be its Seifert-framed winding-two cable. Then

```
ar_C(C_(2,Q)(K))>=2b-1,
mr(C_(2,Q)(K))<=bridge(C_(2,Q)(K))<=2b.
```

Here b=4 for K_n. Thus every such cable has ACTUAL complex augmentation
gap at most one. Moreover ar_C(C_(2,Q)(K))>2 ar_C(K).
This strict inequality is an explicit family consequence; strictness of
augmentation products was already known, and no priority claim is made.

**Corollary OddCableSameRepresentationObstruction.** Let rho_K be the
specific retained A_(2b) quotient, followed by its natural complex
sum-zero representation of dimension 2b-1. For odd p>=3 and gcd(p,Q)=1,
there is no rank-one meridional representation of the (p,Q)-cable whose
restriction to the companion group is conjugate to rho_K. This concerns
that SAME companion representation; it is not an augmentation-rank upper
bound for odd-p cables using other companion representations.

The common braid prefix is

```
U=[1,2,-3,-3,2,-3,-3,2,1,1,2,2,2,-3].
```

## 1. Primary statements inspected

Cornwell--Hemminger, arXiv:1408.4110v1, Theorem1.2, Proposition2.1,
Proposition5.2 and Section5.2 were inspected in the actual PDF.
[primary PDF](sources/cornwell-hemminger-1408.4110v1.pdf)
Their multiplication theorem assumes both companion and pattern attain
their displayed braid indices. Their general conjecture is a LOWER bound;
Proposition5.2 also supplies a lower bound. Neither supplies the upper
bound ar(cable)<=p ar(companion) needed for simple deficit multiplication.
Strictness of augmentation products is already explicit there: printed
p22 gives alpha=sigma1^5,gamma=sigma1 with satellite augmentation rank4,
while the companion/pattern-closure product is2*1. This review makes no
novelty claim for strict supermultiplicativity itself.
Proposition2.1 identifies blackboard and Seifert framings: the two-parallel
braid followed by sigma1^q is the (2,Q)-cable with Q=q+2w(beta).

Gao--Gu--Li--Wang, arXiv:2611.15673v1, Proposition2.16,
Proposition2.22, Lemma2.10 and Theorem1.3 were inspected.
[primary PDF](sources/gao-gu-li-wang-2611.15673v1.pdf)
Their Whitehead construction gives augmentation LOWER bounds under
explicit peripheral conditions; their unconditional meridional bound
points UPWARD. The square-knot example has only a one-unit proven gap.
A finite orthogonal rank-three point is insufficient for the doubling
hypotheses: its distinguished longitude eigenvalue is +/-1, failing the
first condition lambda^2!=1. An appropriate actual peripheral curve is a
separate obligation for these companions.

Boileau--Dutra--Jang--Weidmann, arXiv:1608.01570, Definition1 and
Theorem1 were inspected. [primary PDF](sources/boileau-dutra-jang-weidmann-1608.01570v1.pdf)
Their braid-satellite result starts from meridionally tame companions.
The released knot's equality mr=bridge alone does not verify that premise.
No unrestricted cable-mr multiplication formula is imported here.

## 2. A two-parallel crossing and ordered packets

Let iota2 replace every strand by two blackboard parallels. With the
left-to-right Hurwitz convention, substitute for sigma_i the literal word

 (2i,2i-1,2i+1,2i),

and substitute its inverse for sigma_i^-1. On four group colors
(a1,a2,b1,b2), this word sends them to

 (b1,b2,B^-1 a1 B,B^-1 a2 B),       B=b1 b2.             (1)

This follows by four substitutions and holds in a FREE group. Thus we
can compute cabling independently by treating each ordered pair as a
packet, with its product as the conjugating element.

Our packets are disjoint transpositions, so each packet commutes with its
own product. Standard Hurwitz band and Coxeter formulas therefore hold
for these ordered packets, with simultaneous conjugation of both entries.
This can also be proved directly from (1). In particular a band between
mobile packets restores each intervening protected packet, in its
individual order, and the product R of consecutive mobile Hurwitz moves
satisfies R^d=conjugation by p^-1 on a d-packet tuple of product p.
For detail, one R moves the first packet past all others and conjugates
it by their product S. Since its own product commutes with both of its
entries, S^-1 P1 S=p^-1 P1 p. The total product p is preserved.
After d such moves every packet has been conjugated once by p^-1.
One additional move gives precisely the R^(d+1) formula used below.
The band formula follows by moving the first packet rightward through
the intervening packets and moving the last packet back leftward by
inverse moves; each intervening packet is restored. This is a direct
induction using (1), not merely a group-product calculation.

## 3. Four-bridge family: all n and all odd slopes

For [the four-bridge family](PRIME_AUGMENTATION_GAP_FAMILY.tex),
beta_n=U sigma2^(3+6n), n>=0, use the eight TRANSPOSITION colors

 t=((01),(23),(23),(14),(05),(16),(12),(67)).              (2)

These are the factors of the retained four A8 colors, reversing only the
factor order of the second color. Direct use of (1) gives

 t*iota2(U)=((23),(01),(05),(14),(23),(16),(12),(67)),
 t*iota2(U sigma2^3)=
             ((23),(01),(23),(14),(05),(16),(12),(67))=:e.

Directly e*iota2(sigma2^6)=e. These are identities between explicit
finite permutation arrays, not a finite inference about n. Induction
therefore gives t*iota2(beta_n)=e for ALL n>=0. The first two entries
commute, so any odd power sigma1^q changes e back to t.

The released meridional-rank equality and displayed four-strand braid
also verify minimal braid index4, as used by the cited framing proposition.
Consequently (2) colors iota2(beta_n) sigma1^q for all n>=0 and all
odd integers q, positive or negative. Since Q=q+2w(beta_n), these are
ALL knotted Seifert-framed two-cables, not just selected slopes.

## 4. Uniform families: all bridge numbers and twist parameters

Use [the uniform construction](UNIFORM_PROTECTED_AUGMENTATION_GAPS.md), b>=5, d=b-2,

 delta_b=sigma2 sigma3 sigma4 sigma3^-1 sigma5...sigma_(b-1),
 r=(d+1)(1+d^2 k),       beta_(b,k)=U delta_b^r, k>=0.

Order the factors of its top colors by

 P0=((01),(23)), P1=((23),(14)),
 P2=((05),(16)), P3=((12),(67)),
 P4=((08),(79)),
 Pj=((2j-2,2j),(2j-1,2j+1)), j>=5.                     (3)

The prefix U has the same ordered-packet calculation as Section3. Its
first packet reverses, its fourth packet and all later ones remain
unchanged, and its mobile packets start as

 A=((05),(14)), H=((23),(16)).

Frame the mobile tuple exactly as in the released uniform proof. Put
X=P1, Y=P2, and let Z_i be the later framed packets. Their product is

 p=(0 8 10 ... 2b-2 5)(1 6 9 11 ... 2b-1 4)(2 3)(7).

Its order is d(d+1). The displayed cycles verify AS ORDERED PACKETS

 p X p^-1=H,    Z_i=p^i Y p^-i,    p^d Y p^-d=A.         (4)

For a d-packet tuple of product p, the exact Coxeter formula gives

 R^(d+1)(A,H,Z1,...,Z_(d-2))
  =(p^-1 H p,p^-1 Z1 p,...,p^-1 Z_(d-2) p,p^-2 A p^2)
  =(X,Y,Z1,...,Z_(d-2)),

by (4). Thus delta_b^(d+1) returns all mobile packets in their original
individual order. Moreover R^(d^2(d+1)) conjugates by
p^-d(d+1)=1, so the additional k-dependent power fixes every packet.
All protected packets are restored individually by the band formula.

It follows for EVERY b>=5,k>=0 that iota2(beta_(b,k)) reverses only
P0's internal order. Appending any odd sigma1^q returns the entire tuple
(3). This is an unrestricted proof in b,k, using the explicit cycle
formula, rather than an extrapolation from checked bridge numbers.

## 5. Actual rank and strict nonmultiplicativity consequences

The input transpositions (3) generate S_(2b). Their distinct edges form
a spanning tree: the first four packets supply the connected old seven-edge
tree on eight symbols, and every additional packet attaches two fresh
symbols by two edges. The natural sum-zero complex representation of
S_(2b) has dimension 2b-1, is irreducible, and a transposition has nonzero
rank-one displacement. For irreducibility, subtract a transposed vector
from any nonconstant vector in an invariant subspace; this produces a
coordinate difference, whose conjugates span the sum-zero space.
Cornwell's representation-to-augmentation correspondence
(arXiv:1310.7526v3, Proposition2.12 and Lemma3.14, with Definition1.8,
[primary PDF](sources/cornwell-kch-1310.7526v3.pdf)) therefore gives

 ar_C(C_(2,Q)(K)) >= 2b-1.

The explicit cable braid has 2b strands, so

 mr(C_(2,Q)(K)) <= bridge(C_(2,Q)(K)) <= 2b.

Hence its ACTUAL augmentation gap is at most one. This conclusion needs
no cable meridional-rank equality. Since the companion has ar_C(K)<=b-1,
we also have the strict inequality

 ar_C(C_(2,Q)(K)) > 2 ar_C(K).                          (5)

This proves TwoCableGapBarrier. Thus exact augmentation-rank multiplication
fails here outside the
maximal-rank hypotheses; (5) is consistent with Cornwell--Hemminger's
lower-bound conjecture. It also explains concretely why multiplying the
companion's one-unit deficit is invalid.

### Exact bridge and braid indices for the two-cables

Each companion in the theorem is nontrivial and has bridge number \(b\),
by the four-bridge and uniform companion theorems used above. Its
preferred \((2,Q)\)-cable with \(Q\) odd is a knot in the companion solid
torus with algebraic winding number two. Every meridian disk thus
meets the cable at least twice, and a meridian disk of the standard
two-strand pattern meets it in exactly two points. The geometric index
is therefore two. Schubert's satellite bridge inequality, in
[Schultens's proof](sources/schultens-bridge-additivity.pdf), Definition 2
and Theorem 1, gives \(\operatorname{bridge}(C_{2,Q}(K))\ge 2b\).
The displayed cable braid has \(2b\) strands. Since the bridge number is
at most the braid index, both indices equal \(2b\).

The complex augmentation rank is at most the meridional rank in the
stated \(U=1\), \(\mu\ne1\) convention. Indeed, the converse direction of
Cornwell's correspondence gives an irreducible rank-one meridional
representation of dimension equal to the augmentation rank. If
\(m_1,\ldots,m_r\) generate the cable group by meridians, the sum of the
lines \(\operatorname{im}(\rho(m_i)-I)\) is nonzero and stable under each
\(\rho(m_i)\) and its inverse. It is therefore the full representation,
of dimension at most \(r\). Hence for every stated companion and slope
\[
  2b-1\le \operatorname{ar}_{\mathbb C}(C_{2,Q}(K))\le
  \operatorname{mr}(C_{2,Q}(K))\le
  \operatorname{bridge}(C_{2,Q}(K))=2b.
\]
On this family the meridional rank conjecture is equivalent, cable by
cable, to excluding generation by \(2b-1\) meridians. An augmentation
of rank \(2b\) suffices. If a cable has meridional rank below its bridge
number, both its meridional rank and its complex augmentation rank are
exactly \(2b-1\).

## 6. Independent satellite-group check and odd-p limitation

A p-strand closed braid pattern has group F_p semidirect Z. If its fiber
meridians are x_i, their product is the outer meridian M, and its outer
longitude L satisfies L^-1 x_i L=phi(x_i). For the (p,Q) torus braid,
phi=R^Q and R^p is conjugation by M^-1. The companion group is glued by
identifying this M,L with its meridian and preferred longitude.

For b4,n0 an independent ORIGINAL four-strand color transport computes
a blackboard longitudinal conjugator

 J=(3,2,0,1,5,7,4,6),       M=(01)(23),       w=7.

At each undercrossing the conjugating path acquires the over-meridian
with the crossing sign. Concatenating these paths along the one strand
cycle gives J, the blackboard parallel. Its exponent sum is the writhe7;
removing M^7 therefore gives the zero-linking preferred longitude.
The preferred longitude is L=J M^-7 in this convention. It commutes
with M and swaps the two transpositions x=(01),y=(23) by conjugation.
Also xy=M. For every odd Q, the two-strand pattern action swaps x,y
because they commute. Thus the mapping-torus amalgam relations hold,
independently of the expanded two-parallel word. In particular this
checks the blackboard-q1/Seifert-Q15 example by a different presentation.
Inverting the longitude convention or multiplying by a meridian power
does not change its factor-swapping property or the all-odd-slope scope.

The same factor-swapping conclusion holds for EVERY family parameter,
not only this independently checked example. Here is the required transport
argument. Let the knot braid carry top strand i to bottom strand pi(i),
and let J_i be its ordered-packet conjugator; thus the output in position
pi(i) is J_i^-1 P_i J_i. Our proved all-parameter packet identity says
that output packet j equals P_j except at j=0, where its two entries
are reversed. The permutation pi is one b-cycle. Follow that cycle from
strand0 and compose the conjugators J_i in traversal order. Internal
packet reversal commutes with simultaneous conjugation, so the resulting
blackboard longitudinal conjugator J satisfies

                    J^-1 P0 J=reverse(P0).

Exactly one reversal occurs around the cycle. The preferred longitude is
obtained by multiplying J (or its inverse if the longitudinal orientation
is reversed) by the appropriate power of the meridian M. Since M is the
product of the two DISJOINT transpositions in P0, every power of M fixes
each of them individually under conjugation. Thus the preferred longitude
still swaps the two factor lines for every n and every b,k.
The writhes are w(beta_n)=7+6n and
w(beta_(b,k))=4+(b-3)(d+1)(1+d^2 k), respectively. In particular the
blackboard-to-preferred correction Q=q+2w is an even shift in every case.

Odd p>=3 does NOT inherit this reflection lift. There is even a stronger
necessary condition for extending the SAME companion representation.
In the torus-pattern group the element Z=M^Q L^-p is central. In any
representation whose pattern meridians are conjugate rank-one matrices,
Z acts by the SAME scalar on every meridian root line, hence on their
span. As im(M-I) lies in that span, Z must act as a scalar on im(M-I).

For the natural S_(2b) companion representation restricted to A_(2b),
im(M-I) is the two-space spanned by e0-e1,e2-e3. Here M=-I, while L
swaps those two lines; therefore L has two opposite eigenvalues. If p
is odd, M^Q L^-p has two distinct eigenvalues on that two-space,
contradicting the necessary scalar condition. Thus this particular
companion representation has NO rank-one extension to an odd-p cable,
regardless of the satellite meridian eigenvalue. Merely changing it
from a transposition to a complex reflection cannot avoid the obstruction.
This proves OddCableSameRepresentationObstruction. It does not bound
augmentation ranks of odd-p cables obtained from other companion representations.

## 7. Remaining route and hyperbolicity

The primary formulas alone do not establish a gap-two odd-p cable or
Whitehead double of these companions: the needed augmentation UPPER
bound and a sufficiently strong meridional LOWER bound are missing.
A genuine satellite with essential companion torus stays nonhyperbolic
even if its pattern exterior is hyperbolic. A hyperbolic deformation
would have to destroy that torus and independently preserve the two
rank inequalities; neither is a continuity consequence of the inspected
satellite theorems. No such preservation theorem or new hyperbolic example is established here.


## 8. Exact verification and assurance

The [standard-library checker](verify_two_cable_augmentation_rank_barrier.py)
checks the two-parallel crossing formula in a free group, the exact
four-bridge induction identities, and 64 parameter controls using both
literal expanded crossings and ordered-packet transport. It separately
checks the original four-strand longitude and the satellite mapping-torus
relations at blackboard q=1/preferred Q=15. Parity and changed-crossing
controls are rejected. Optimization with `python -O` is rejected because
the exact assertions are required.

```
python3 research/meridional-rank/representations/verify_two_cable_augmentation_rank_barrier.py \
  --write-receipt research/meridional-rank/representations/two_cable_augmentation_rank_barrier_verification.json
```

The [receipt](two_cable_augmentation_rank_barrier_verification.json) and
[source ledger](two_cable_augmentation_rank_barrier_source_ledger.json)
pin the exact primary PDFs and existing companion proofs. No duplicate
companion proof edition is introduced. The written ordered-Coxeter and
peripheral-transport arguments establish the unrestricted b,k,n,Q scope;
the bounded controls do not replace them. This is internal proof and
computational evidence, without external review or formal verification.
