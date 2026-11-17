# All primitive binomial two-stage keys: exact periodic gate and polynomial existence

15 November 2026. PF2RULED63. Programme archived; full goal active.

The exact reduction of PF2RULED62 extends to every coprime first pair
(p,q), with p,q>=2, and every last winding r>=2. All three parameters
are unbounded. The entire eta0 peripheral S4 problem is periodic modulo12
in (p,q,r), and determinant-square9 adds only modulo9 data. The complete
certificate retains 244 of 1152 admissible ordered residue cells modulo12;
the other 908 cells are excluded for the actual-binomial quartic class.
The whole orbit-tree test removes 112 of 596 representation classes.

There is also a necessary correction to the research direction.
Zannier's original existence theorem supplies an immersed polynomial A1
curve of eta0 for EVERY p,q,r in this class, without constructing its
coefficients. For coprime odd p,q and every r=9 modulo36 its knot passes
the peripheral S4, orbit-tree and determinant9 gates. This is an unbounded
family of actual polynomial immersions with passing infinity data, not
a family of finite Keller maps. Its finite multiple-value types are not
classified here. The missing affine-cover and marked-source constraints
remain essential.

## 1. Uniform parameters and primary peripheral conventions

Keep the immersive primitive ACTUAL binomial key from PF2RULED60:

    deg X=p*r, deg Y=q*r, gcd(p,q)=1, p,q,r>=2,
    F=X^q-Y^p, deg F=k<p*q*r, gcd(k,r)=1.

PF60 and PF59 force eta0 in an actual surviving quartic. Put

    c=p*q-p-q,  k=r*c+1,  K=C_(r,k)(T(p,q)).            (1)

For these coprime p,q, c is positive and ODD. In particular r*k is even,
and gcd(r,k)=1 automatically. General keys with essential lower-weight
terms are not replaced by this binomial, and longer key sequences are
not assumed to have two stages.

The companion group and its positive meridian/longitude are

    <a,b | a^p=b^q>,  M=a^u*b^v, q*u+p*v=1,
    L=a^p*M^(-p*q).                                   (2)

The abelianization sends a to q and b to p. A primitive transverse
meridional curve on the torus gives the displayed Bezout word. Changing
u,v by p*t,-q*t leaves the word unchanged because a^p=b^q is central.

For the cable set

    h=M^k*L^r, z^r=h, beta=M^(-c)*L^(-1),
    mu=z*beta, lambda=h*mu^(-r*k).                     (3)

Clay--Watson, complete arXiv-v1 PDF6, gives the root amalgam and both
peripheral words for general positive cabling. Their winding is r and
their other parameter is k. In their positive Bezout convention choose
v_C=r-1 and u_C=c*(r-1)+1; then r*u_C-k*v_C=1. Their meridian becomes
beta*z after using z^r=h. Conjugating by beta^-1 gives (3), while fixing
M,L,h. Their longitude gives precisely the displayed lambda. Thus the
generalization retains both framing and the common companion basis.
No left-orderability or Dehn-surgery theorem is imported.

Equivalently, in the cable-space group,

    M=beta^r*h, L=beta^(-k)*h^(-c),
    [beta,h]=[mu,h]=1, (mu^-1*beta)^r*h=1.              (4)

The positive unknot specialization is z^r=M^k, as required. These
identities generalize the c1 formulas in PF61/62 without changing them.

## 2. The complete representation problem is finite-periodic

The actual quartic infinity action is onto S4, with mu a transposition
tau and lambda1. Since r*k is even, (3) forces h1. Therefore the exact
necessary and sufficient equations for this PERIPHERAL representation are

    a^p=b^q, M=a^u*b^v, L=a^p*M^(-p*q),
    M=beta^r, L=beta^(-k), (tau*beta)^r=1,
    tau a transposition, <a,b,beta,tau>=S4.             (5)

For sufficiency, z=tau*beta^-1 is conjugate to (tau*beta)^-1,
so its r-th power is1, exactly the root relation with h1.
The definition beta=M^-c*L^-1 follows back from these equations, since
k-r*c=1. It must not be chosen independently from a different companion
peripheral basis. The three subgroup images used for the orbit tree are

    A=<a,b>, B=<beta,tau>, P=<M,L>.                    (6)

S4 has exponent12. For p0,q0 modulo12, actual coprimality implies
gcd(p0,q0,12)=1, giving 96 ordered residue pairs. Choose ANY u0,v0
satisfying q0*u0+p0*v0=1 modulo12. This modular choice is sufficient:
for any actual coprime lifts p,q, write q*u0+p*v0=1+12*T, and choose
m,n with q*m+p*n=-T. Then u=u0+12*m,v=v0+12*n is an exact Bezout pair.
All such choices give the SAME permutation M by (2). The checker verifies
this independently across every modular Bezout pair and every torus image.

Thus (5), the entire simultaneous-conjugacy classes, and the complete
orbit graphs depend only on (p,q,r) modulo12. Residues0 and1 remain
valid residues of actual degrees >=2; they are not interpreted as literal
zero-degree or unknot parameters. This exact exponent-and-Bezout argument
is the unbounded proof. The following finite enumeration is exhaustive
after that proof, rather than a test over bounded polynomial degrees.

| Certificate quantity | Exact value |
|---|---:|
| Admissible ordered cells (p0,q0,r0) modulo12 | 1152 |
| Cells with a peripheral S4 representation | 244 |
| Cells with a tree representation | 244 |
| Peripheral representations, counting all cells | 14304 |
| Simultaneous-conjugacy classes | 596 |
| Tree classes retained | 484 |
| Cycle classes rejected | 112 |

Every conjugacy orbit has24 representations. The cycle classes have
graph rank1; 64 have companion/cable image orders24/6 and 48 have12/24.
The tree test therefore improves the representation classification but
does not remove an additional entire residue cell in this family.
All orbit data are in the machine-readable receipt; a nonempty cell is
not treated as an affine-cover certificate.

The small companion images C3 and C4 MUST be retained. In particular,
a C3 companion and one transposition can generate S4. The old exploratory
script research/.scratch/test_wright60_peripheral.py discarded companion groups of
order below4 and did not certify its general peripheral convention.
Its residue count is not adopted. The new proof and exhaustive certificate
derive the convention from (3) and include ALL companion images.

## 3. The determinant gate and the exact modulo36 predicate

The first torus determinant is

    D0=1 if p,q odd;
       q if p even;
       p if q even.                                  (7)

At most one of p,q is even. Since c is odd, r even gives k odd and
r odd gives k even. The satellite formula consequently becomes

    D_K=k if r even,  D_K=r*D0 if r odd.               (8)

PF2RULED55 requires 9|D_K for the actual quartic. Combining (5)-(8)
depends only on (p,q,r) modulo36. Among the 31104 cells satisfying the
NECESSARY modular coprimality gcd(p,q,36)=1, exactly2772 pass both the
representation/tree-existence and determinant-square tests. The actual
integers p,q must still be coprime. These are residue counts, not counts
or densities of curves, maps or realizations.

The reusable gate is: reduce actual p,q,r modulo12 and retain the full
tree representation list from the receipt; then require9|(r*c+1) for
even r, or9|(r*D0) for odd r. No lift through an affine complement follows.
The p2/q3 slice reproduces PF62 exactly: r3mod6 or r8mod18, with its
additional representation-level cycle exclusions.

## 4. Zannier supplies eta0 polynomial immersions for every parameter

Zannier, Acta Arithmetica71(1995), Theorem1 on printed114/PDF8, together
with equations(3)-(6) on printed113-114/PDF7-8, proves attainment of the
polynomial-difference degree bound for prescribed root multiplicities.
We import that existence theorem, not an independent reproof or a claim
of rational-coefficient realization.

Apply it with common degree N=p*q*r, the first partition q repeated p*r
times, and the second partition p repeated q*r times. Their common gcd
is1. The attained bound is

    N-(p+q)*r+1=r*c+1=k>0.                            (9)

Over C, absorb constants into roots to write the two attained polynomials
as X^q and Y^p. Their degrees are p*r and q*r. Equality in the preliminary
bound on PDF7 forces them coprime and their difference squarefree.
The identity from PF60 gives

    Y*F'-p*Y'*F=X^(q-1)*(q*X'*Y-p*X*Y'),
    degree(q*X'*Y-p*X*Y')=0.

The Wronskian is a nonzero constant, so X',Y' never vanish together.
The normalization is immersive. It is finite because a nonconstant
polynomial Y makes t integral over C[Y]. Its field degree divides
degX,degY,degF, whose gcd is gcd(r,k)=1, so it is birational. It has
normalization A1 and one place at infinity. The first nondivisible
binomial defect k drops the pole-order gcd directly from r to1, giving
the two-stage knot in (1) under the retained infinity convention.

Hence EVERY coprime p,q>=2 and r>=2 has such a complex polynomial eta0
immersion. Its finite self-intersections are NOT classified by the
existence theorem. In particular no nodality, double-value bound, affine
S4 representation or marked plane finite normalization is inferred.
PF60's separately certified six-node curve retains those stronger local
properties at its specific parameters.

## 5. A closed all-odd classification and an infinite passing family

Suppose p,q are coprime odd integers >=3. Every S4 image of (2) is cyclic.
Indeed if the common central power a^p=b^q is a 4-cycle or a 3-cycle,
both a,b lie in its cyclic centralizer. If it is an involution, odd
exponents force a=b to be that involution. If it is1, both a,b have odd
orders1 or3, and coprimality prevents both exponents p,q being divisible
by3, so at least one generator is1. The image is again cyclic.

For a full S4 cable image, this cyclic companion must be C3 or C4;
a trivial or order2 companion with one more involution cannot generate
S4. The meridian M generates the cyclic companion and L=1.

For C3, beta=M^-c has order3. Full S4 requires tau to join its three-letter
orbit to the fourth letter, so tau*beta has order4. Equations(5) are
equivalent to4|r and3|k. The orbit graph is a tree: the companion and
edge have two orbits and the cable has one. The determinant-square gate
upgrades3|k to9|k.

For C4, beta is a 4-cycle since c is odd. Full S4 requires an adjacent
transposition tau, making tau*beta have order3. Equations(5) require
3|r and4|k. Since odd p,q give c=3 modulo4, this means r=9 modulo12.
The determinant gate is9|r, leaving r=9 modulo36. All three subgroup
actions in (6) are transitive, so the graph is again a tree.

Thus for coprime odd p,q the combined abstract peripheral/tree/determinant
gate passes EXACTLY in the following two alternatives:

    r=9 modulo36;  OR  4 divides r and 9 divides r*c+1.        (10)

Sufficiency has a uniform explicit representation. Choose M a 4-cycle
in the first alternative, or a 3-cycle in the second, and put

    a=M^q, b=M^p, beta=M^-c, z=tau*beta^-1.

Choose tau adjacent to the 4-cycle, or joining the 3-cycle to its fixed
letter. Then z has order3 or4, respectively; z^r=h=1, the cable meridian
is tau, and its preferred longitude is1. The image is S4 and the graph
is a tree. All formulas are independent of polynomial coefficients.

Combining (9) with the FIRST alternative of (10) gives an infinite family
of realized polynomial eta0 immersions, for ALL coprime odd p,q and ALL
r9mod36, whose infinity knots simultaneously pass the peripheral S4,
orbit-tree and determinant-square9 tests. These passing knots do not
certify sign Picard torsion or its linking form, a finite quartic cover,
the quartic double-value condition or an affine-plane source.

## Remaining obligation

Use the complete finite target algebra, sign Cartier Picard/linking data,
source lattice and canonical marking to distinguish the passing curves.
The binomial class now has an exact bounded monodromy test and general
polynomial existence; increasing coefficient depth or sampling larger
parameter boxes is not the missing proof. Non-binomial keys, longer
key sequences, and higher-degree cover exclusions remain in the goal.

[Checker](verify_binomial_periodic.py), [receipt](binomial_periodic_receipt.json)
and [primary manifest](sources/binomial-periodic/manifest.json) retain the
complete finite residue census, independent explicit passing controls,
and primary existence/peripheral imports. No full quartic or higher-degree
exclusion, finite Keller pair, affine-cover extension, formal verification,
independent full reproof, external review, novelty priority or goal
completion is claimed.
