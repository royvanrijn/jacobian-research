# Actual Seifert sections force common winding three

15 November 2026. PF2RULED75 determines the ACTUAL section twists in
PF74's three-cone quotient, without inferring them from homology.
It excludes every non-trefoil S3/S3 common source with r>=9, uniformly
in both first exponents, winding and all source orders and poles.
PF74 already excluded every trefoil common source. With the earlier
representation closures, every surviving COMMON source in the actual
primitive-binomial TWO-STAGE QUARTIC class therefore has

    r=3; e=2b,o=3t coprime, t odd; (b,t)!=(1,1).            (1)

This is a necessary reduction, not closure of r=3. The full problem,
passing complementary sources, other keys, longer sequences and higher
degree stay open. The source obstruction used here is PF74's ALL-DEGREE
no-left-orderable-image theorem under its explicit common-coordinate,
smooth, entire disjoint-line hypotheses. No extra higher-degree source
classification is asserted.

## 1. Retain the actual cone and boundary sections

Use PF73/74 notation

    e=2b,o=3t,r=3s, A=bt, c=6bt-2b-3t,
    k=3sc+1, D=k/2-As, V=c-2A=4bt-2b-3t,
    B=k/2-rA=(3sV+1)/2.

s,t are positive odd integers, gcd(2b,3t)=1 and k is even.
If (b,t)!=(1,1), V is POSITIVE: for b=1 it is t-2>=1;
for b>=2 it is at least 2b-3>=1. Thus B>=2. Also

    k+B=3D, D>0.                                         (2)

The main three-sheet cable AFTER its actual branch-meridian fillings
has base annulus, one cone of order s, and regular fibre h. Let z be
PF63's original root with z^r=h, beta the original boundary word, and
mu=z*beta the filled branch meridian. The main action has z a 3-cycle,
beta a transposition. With oriented arcs and conjugate lifts chosen,
write

    X=lift(z^3), Y=lift(beta on its fixed sheet),
    Z=lift(beta^2 on its two-sheet orbit).

These are the actual cone section and the two boundary sections.
In particular X^s=h, not merely an abelian relation. The fibre h is
central. The base-annulus one-cone presentation has

    XYZ=h^j                                                (3)

for an integer j; all three exponents are +1, by the oriented base
cover of z*beta*mu^(-1)=1. The remaining mu boundaries are filled
with their actual meridians. We now determine j integrally.

PF63's words M=h*beta^r and L=h^(-c)*beta^(-k) give on the fixed
and two-sheet boundary tori respectively

    n=h*Y^r, L_1=h^(-c)*Y^(-k),
    m=h^2*Z^r, L=h^(-c)*Z^(-k/2).                         (4)

Here m=M^2 on the central torus, n=M on the fixed torus. PF73's
lifted meridian disk gives m+n=0 in integral H1 of this two-boundary
piece. Its annulus one-cone presentation, with X^s=h, has H1 FREE
of rank2 and h nonzero: eliminating h and one boundary generator
in (3) leaves two unrestricted generators. Abelianizing (3)-(4),

    0=m+n=3h+r(Y+Z)=r*j*h.

Therefore j=0. This determines the actual relative section Euler
integer; no exterior pushforward is descended through a meridian
filling, and no separate individual-meridian vanishing is assumed.

## 2. Exact group presentation of the quotient

PF74 simultaneously abelianizes both actual companion leaves.
The degree1 companion imposes L_1=1. The reduced main companion
imposes eta=m^A*L=1, its null slope, and kills its independent torsion.
Using (4), these relations are

    Y^k=h^(-c), Z^B=h^(-V).

Consequently the actual target end group SURJECTS onto the group

    G=<X,Y,Z,h | h central,
         X^s=h, Y^k=h^(-c), Z^B=h^(-V), XYZ=1>.             (5)

For non-trefoil parameters B>0. Put f=h^(-1), X'=X*f. Then (5) becomes

    (X')^s=f^(s-1), Y^k=f^c, Z^B=f^V, X'YZ=f, f central.

This is EXACTLY the standard negative-star plumbing group with
normalized data

    e0=-1; (s,s-1),(k,c),(B,V).                            (6)

For s=1 omit the order1 pair (1,0). When s>1 every fraction is strictly
between0 and1. Each pair is primitive: gcd(s,s-1)=gcd(k,c)=gcd(B,V)=1,
since k=3sc+1 and 2B=3sV+1. PAIRWISE coprimality of cone orders
is not assumed or needed. In particular k and B may share primes.
The group identification follows from its generators and relators,
not from the determinant or a guessed oriented homeomorphism.

The rational Euler number of (6) is

    -1+(s-1)/s+c/k+V/B
      =-(k+B)/(3skB)=-D/(skB)<0.                          (7)

Thus (6) gives a connected negative-definite plumbing graph. Its
abelianization is Z/D, consistently with PF74's independent gluing
calculation. Deep rational-link/orderability results are imported
with the precise audited hypotheses retained from PF74.

## 3. One-step non-rationality at every s>1

s is odd, so s>1 means s>=3. The three cone orders s,k,B all exceed1;
their normalized fractions in (6) yield three NONEMPTY negative
continued-fraction legs. The central vertex has weight-1 and valency3.
For Laufer's initial reduced cycle l0=sum E_v, its intersection with
that central vertex is

    (l0,E_centre)=-1+3=2.                                 (8)

Laufer's criterion in Nemethi Section2.1 requires every positive
intersection in a computation sequence to equal1 for a rational graph.
Choosing the centre at the first step violates that criterion. Hence
the graph is NONRATIONAL, uniformly in every exponent and winding.
Equivalently Nemethi equation(1.13) has the SINGLE witness l=1:
the three floors are all -1, so their sum is -3=l*e0-2.
Neither argument searches an interval of increasing denominators.

Nemethi Theorem1.6 implies that this group's nontrivial fundamental
group is LEFT-ORDERABLE. Nontriviality also follows from its quotient
onto the nontrivial base triangle group (all orders>=2 and at least
s>=3,k>=10); no assumption D>1 is required for this step. Alternatively,
1/s+1/k+1/B<1: k>=3s*3+1 and B>=(3s+1)/2 give a strict bound
already at s=3. Thus a hyperbolic triangle quotient is available.
The actual end group has this nontrivial left-orderable image via (5),
contradicting PF74's all-degree source theorem. This closes all
non-trefoil S3/S3 common sources for s>=3, or r>=9.

PF74 handles the exceptional negative-V case b=t=1 at EVERY winding.
PF68/71/72 had closed every non-S3/S3 common representation sector.
Together they prove (1) for the ENTIRE stated actual-binomial common
class. They do not close the complementary source pencil.

## 4. Why the remaining r=3 case needs the actual two pieces

At s=1, X'=1 in (6), leaving two cone factors. The quotient group is
ABELIAN: YZ=f central makes Y and Z commute. Its abelian presentation
on Y,f has rows

    [k,-c], [B,-(B-V)].                                   (9)

The determinant has magnitude D, and gcd(k,c)=1 makes it cyclic Z/D.
Thus this quotient genuinely loses the orderability obstruction at
r=3. A finite cyclic quotient is not evidence of an actual source;
the full target group still contains the companion information.

There is a useful exact replacement of the residual problem. The
ACTUAL main cable at s=1 is an oriented circle bundle over an annulus
with no exceptional fibre, hence T2 times I. This is a homeomorphism
of the actual cable, unlike the virtual fillings in Sections1-2.
Eliminate X=h in (3), giving Z=h^(-1)*Y^(-1). All generators of this
annulus group commute. For r=3,k=3c+1, (4) gives

    h=n^k*L_1^3, Y=n^(-c)*L_1^(-1),
    m=n^(-1), L=n^(k/2)*L_1.                              (10)

The pair (n,L_1) is a genuine integral basis: its exponent determinant
in (h,Y) is -k+3c=-1. The central pair (m,L) is also a basis.
Therefore the actual r=3 target is exactly the gluing of

    E': base D, cones (b,2b,t),
    E_1: the actual T(2b,3t) knot exterior,
    m -> -n, L -> (k/2)*n+L_1 in oriented H1.               (11)

The ramification knot is represented by the primitive torus slope
h=k*n+3*L_1; gcd(k,3)=1. Retain that mark. No exterior is replaced
by a solid torus in (11), and no JSJ classification is inferred solely
from its matrix. The unrestricted r=3 family, at all first exponents,
is now a concrete two-Seifert-piece marked gluing problem.

The next useful gate is orderability/L-space gluing or rational
plumbing recognition for these ACTUAL pieces, or actual affine finite
algebra. A primary literature lead is Rasmussen's
[L-space intervals for graph manifolds and cables](https://arxiv.org/abs/1511.04413);
only its primary record/abstract was consulted here. No theorem from it
is imported without a complete-page audit and checks of its hypotheses.
The bounded arithmetic diagnostic in research/.scratch/wright75 is explicitly
exploratory; a passing binary lattice there is not full marked linking,
an affine source, or a finite map. It does not establish (1).

## Evidence and exact limits

[The checker](verify_seifert_sections_winding_three.py) seals PF74/73
and their predecessors, checks all64 S3/S3 incidence types, the exact
relative-section and group generator changes, primitive normalized
coefficients even with shared cone primes, the unbounded Euler and
one-step witnesses, ordinary satellite controls, and the genuine r=3
annulus gluing basis and mark. [The receipt](seifert_sections_winding_three_receipt.json)
and [manifest](sources/seifert-sections-winding-three/manifest.json)
retain the already completed primary-page audits. Universal claims
come from the written proofs, not the finite controls.

All-degree input is the explicit common-pencil source theorem, not an
unconditional higher-degree classification. New target closure is the
actual-binomial TWO-STAGE QUARTIC COMMON sector with r>=9, and the
resulting necessary r=3 reduction. Residual r=3 non-trefoil sources,
COMPLEMENTARY passing cases, other keys, longer sequences, higher degree
and full goal remain OPEN. PF63 polynomial immersions still exist.
No finite Keller map, independent full reproof, formal verification,
external review or novelty/priority claim is asserted.
