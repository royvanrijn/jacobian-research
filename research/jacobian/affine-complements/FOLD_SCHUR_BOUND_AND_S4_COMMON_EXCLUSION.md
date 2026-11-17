# A strict fold Schur bound and the S4/S4 common-source exclusion

15 November 2026. PF2RULED71 strengthens PF67's source-wide Gram
argument for a sole fold. At ALL finite field degrees and with any
number of unramified omitted lines, under PF67's explicit smooth,
whole disjoint-line and common-coordinate hypotheses, the Schur
complement satisfies S>2. Hence

    0 < lambda_end(K_R,K_R) = 1/S < 1/2.                    (1)

For every remaining S4/S4 representation in the actual primitive-
binomial two-stage quartic class, the even first exponent is2 by
PF70. We prove that its ramification knot has order4 and self-linking
EXACTLY1/2, uniformly in both odd parameters and all representations.
Thus its COMMON source pencil is excluded. Together with PF69/70,
the ENTIRE S4/S4 common-source sector is closed, at every first pair,
winding and source order/pole. Complementary sources remain open.

## 1. Positive log orders cannot have zero-order descendants

Use PF67's blowup model X->P1_x times P1_f, resolving only outside
the actual smooth N. All retained omitted-line closures are FREE
leaves of the full boundary Gamma. For a boundary vertex E set

    b(E)=1+ord_E(dx wedge df).

Both initial infinity sections have b=-1. A free blowup creates
b_new=b_parent+1; a crossing blowup creates b_new=b_1+b_2.
Induction shows that no boundary edge has one strictly positive
endpoint and one strictly negative endpoint: the assertion holds
initially, free creation cannot jump across0, and crossing creation
between two nonnegative or two nonpositive endpoints preserves it.

Consequently a new descendant of a vertex with b>=1 has b>=1:
at a free centre its value increases, and at a crossing its other
parent cannot have negative b. Here descendant is the curvette-path
ancestry, including both parents at a crossing. A b=0 vertex cannot
descend a positive vertex. Moreover, a new FREE b=1 component cannot
descend any earlier positive vertex: its b=0 creation parent would
already have to be such a descendant. Thus no retained unramified
leaf D_i (a free b=1 component) descends any other positive vertex.

This extends the particular ancestry observations used in PF67.
It uses the initial boundary and the full blowup rules, not a
bounded valuation census or a normalized skewness convention.

## 2. A b=1 parent stays at squared distance greater than1

Keep the integral curvette vectors m_E in PF67's orthogonal total-
exceptional basis. If e_E is the new unit coordinate of a creation,

    free:     m_E=m_parent+e_E;
    crossing: m_E=m_parent1+m_parent2+e_E.                   (2)

These are the same proximity recurrences that give the whole Gram
matrix there. Initial-section exceptional vectors are0. Let W be
the real span of the m_D_i. Their terminal coordinates are private
among the retained D vectors, which are pairwise incomparable.

Claim: for ANY b=1 exceptional component P distinct from the retained
D_i, with no retained D_i created later on P, one has

    dist(m_P,W)^2>1.                                       (3)

The extra wording about later creation is automatic by Section1;
it emphasizes that no D terminal coordinate can lie below P. The
coordinate e_P is therefore orthogonal to W. A b=1 vertex was
created either freely from b=0, or at a crossing with endpoint
values (0,1). Other possibilities are ruled out by Section1.

If no D_i is an ancestor of P, put w=m_P-e_P. It is nonzero:
obtaining b=1 from the initial b=-1 sections needs at least two
exceptional creations. It is zero at every private D coordinate.
Such a nonzero vector cannot belong to W, since those coordinates
would force every coefficient in a D-vector combination to be0.
Splitting off e_P gives (3).

If some D_i is an ancestor of P, follow the unique b=1 parent
backwards. P cannot be free from a b=0 parent descending D_i.
It must be a crossing of a b=1 descendant and a b=0 vertex. Repeat
along the b=1 parent until reaching D_i. This chain cannot merge
two retained D ancestors: there is only one positive parent at
each b=1 crossing. Along it (2) gives

    m_P = m_D_i + e_P + w,

where w is nonnegative and nonzero. Even at its first crossing,
the b=0 parent contributes a NONZERO exceptional vector; it cannot
be an initial section, since those have b=-1. Later crossings also
contribute their nonzero vectors and intermediate unit coordinates.
Every private D coordinate of w is0: the b=0 parents cannot descend
a D, and the coefficient at D_i's terminal stays exactly1 along
the unique b=1 chain. Hence w is not in W. Again e_P is orthogonal
to both w and W, proving (3). This proves the claim at every depth.

## 3. The all-degree sole-fold source inequality

Now impose PF67's source hypotheses with SOLE positive volume order1:
N smooth affine contains U=A2_(x,f); its complement is the ENTIRE
DISJOINT affine lines R,D_1,...,D_n; f is regular and constant at
finite values on each line; x:N->P1 is regular with a pole on each;
and div_N(dx wedge df)=R, with zero order on each D_i. The values
of f need not coincide. These include the common source of the
surviving quartic packet, but are not automatic in higher degree.

The R leaf has b=2 and was created freely, so its creation parent
P has b=1 and m_R=m_P+e_R. P cannot be a retained D_i, since
that free blowup would increase its boundary degree and destroy
the retained leaf property. Section1 also makes e_R private from
all D vectors. Apply (3):

    S=dist(m_R,W)^2=1+dist(m_P,W)^2>2.                       (4)

By PF61/PF67 the WHOLE source lattice is Q=-Gram(m_R,m_D_i),
and (1) follows by block inversion. No first-pair or target key
restriction was used in this source theorem. The statement here
requires volume order EXACTLY1; no S>a+1 theorem for general higher
ramification order a is asserted.

The constant2 is an infimum already in admissible boundary-cluster
models: a b=0 vertex with vector u of norm1, n free b=1 siblings D_i,
another free b=1 vertex P, and a free R on P give

    G_RR=3, G_R,D_i=1, G_D=I+J,
    S=2+1/(n+1).

Explicit local curvettes are R: x=z^-1+T,f=z^2 and
D_i: x=(i+1)w^-1+W,f=w. Their retained components are free leaves
with volume orders1,0. This is a sharp cluster control, not a new
claim of affine-source or finite-cover realization for those models.
It shows why the strict uniform bound is not replaced by S>=3
when the number of D leaves is unrestricted.

## 4. Integral peripheral homology of the remaining S4/S4 cover

PF69 excludes the genus-one S4/S4 generator types, and PF70 excludes
every even first exponent e>2. Write the remaining odd first exponent
and winding as

    q=3t, r=3s, s,t positive ODD integers.

The actual peripheral equations impose r(q+2)=1 mod4. The companion
meridian M is a 4-cycle, its torus generator a is a transposition
and b is a 3-cycle with a fixed sheet. Set

    H=a^2=b^q, L=H M^(-2q), m=M^4,
    j=(1-r(q+2))/4, h=M^k L^r=H^r m^j,
    k=r(q-2)+1.                                            (5)

H,m are a BASIS of the connected lifted companion boundary torus:
the peripheral action is cyclic4 with H identity, and their exponent
determinant in (M,L) coordinates is4. No preferred longitude is
silently substituted for this actual basis. The ramification mark
after the meridian fillings is [K_R]=h.

We use the following elementary Seifert fact. For an orientable
base-D piece with cone orders alpha_i and primitive twisting
coefficients beta_i, its H1 has rank1. The torsion order is
product(alpha_i)/lcm(alpha_i): take the gcd of maximal minors of
the matrix with rows alpha_i*x_i+beta_i*h=0. At a prime dividing
exactly two alpha_i, the mod-prime dimension is2 (one free and
one torsion generator). These facts follow from the explicit
presentation audited in PF70, and determine the groups used below.

### The companion piece E

Its base is D with cone orders 2,2,t,3t; orders1 are harmless.
The preceding calculation gives

    H1(E)=Z*g direct-sum (Z/(2t))*delta.                    (6)

The map to the target knot exterior on free homology is surjective:
a closed lift of a on a fixed sheet maps to q times the meridian,
and a fixed-sheet lift of b maps to2 times it; gcd(q,2)=1. Choose
g mapping to one meridian. H maps to2q=6t times it. Each cone
relation makes H divisible by its cone order, since beta_i is
coprime to alpha_i. In particular H is divisible by2t, so has no
torsion component in (6). Therefore

    H=6t*g, m=4g+epsilon*delta.                             (7)

PF70's exact peripheral quotient for this piece is Z/(2t).
Reducing (6)-(7) mod2 shows epsilon is ODD: otherwise killing
H,m would leave TWO mod2 generators. It follows from the 2x2
Smith presentation that

    H1(E)/<m> = Z/(8t).                                    (8)

Indeed its matrix is [[4,epsilon],[0,2t]], whose entries have gcd1
because epsilon is odd, and whose determinant is8t.

### The filled cable piece P

The cable cover has a single connected four-sheet piece. Its root
permutation has cycle lengths3,1, so its lifted cone orders are
s,3s. The three meridian-cycle fillings leave a base D with exactly
these cone orders. Each filling meridian is a primitive section
curve on its lifted base boundary and intersects the regular fibre
once; it introduces no additional exceptional order. This argument
keeps the true meridian fillings, including both fixed sheets.
Thus

    H1(P)=Z direct-sum Z/s.                                (9)

Let the target filled cable piece be its original companion solid
torus V. Its homology is generated by L, with M=0. The branched
map P->V sends h to r*L and H to L. Its free homology map is onto
because of the latter boundary loop. The order-3s cone relation
makes h divisible by3s; the torsion in (9) has order s. Consequently

    h=r*H in H1(P).                                        (10)

The boundary relation in (5) also gives h=rH+jm. Hence jm=0.
The pushforward of m is0, so m is torsion of order dividing s.
But gcd(j,s)=1, since 4j=1-r(q+2) and r=3s. Therefore

    m=0 in H1(P),
    H1(P)=Z*H direct-sum Z/s, h=rH.                         (11)

This is the required INTEGRAL peripheral map. It is proved here;
no replacement of an arbitrary noncyclic piece by a homology torus
is assumed. Its torsion summand has zero peripheral coordinates.

## 5. Full end homology, the order4 mark and its half self-linking

Glue (6) and (11) in the common ACTUAL boundary basis H,m. The
gluing graph has one edge. Mayer--Vietoris, retaining the free H
generator and the untouched Z/s summand, gives

    H1(M)=Z/(8t) direct-sum Z/s,
    [K_R]=rH=6rt*g in its first summand.                    (12)

This is a full group formula; the linking form of the two displayed
summands is not asserted to split in those chosen coordinates.
From m=0 and (7),

    4H=24t*g=-6t*epsilon*delta=0.

The 2-primary part of (12) is cyclic of order8. In that part g
is a generator (epsilon is odd), and r,t are odd, so 6rt*g has
order4. Thus K_R has EXACT order4 integrally, with no odd-primary
part. Its quotient group is consequently

    H1(M)/<K_R>=Z/(2t) direct-sum Z/s.                      (13)

Nonsingularity of the linking form on the cyclic 2-primary group
Z/8 means a generator has self-linking u/8 for some ODD u. Every
order4 element is twice an odd multiple of a generator, and hence
has self-linking 4*(odd)^2*u/8=1/2 modZ. Therefore

    lambda_M(K_R,K_R)=1/2.                                 (14)

Odd-primary components are orthogonal to the 2-primary component,
and K_R has no such component. No orientation ambiguity affects
(14), because negating1/2 leaves1/2 modZ. Group orders alone would
not justify it: the exact mark and cyclic 2-primary group are used.

The source strict bound (1) contradicts (14). This closes every
common source in the entire S4/S4 actual-binomial two-stage quartic
sector: e>2 was already excluded by PF69/70, and all e=2 cases
are now excluded in the common pencil. For the complementary source,
PF70's remaining condition gcd(s,t)=1 is still necessary from (13).
It is not shown sufficient. No complementary S4/S4 exclusion is
claimed by (1), whose common-coordinate hypotheses are essential.

## Evidence and exact scope

[verify_fold_schur_and_s4.py](verify_fold_schur_and_s4.py) checks
the blowup/curvette recurrences, inherited whole boundary controls,
the sharp all-count Schur family, Seifert maximal-minor groups,
integral peripheral and marked presentations, and full ordinary
satellite controls through PF69. It audits all 64 planar S4/S4
periodic cycle classes. The unrestricted arguments above, rather
than the finite controls, prove the universal claims.
[fold_schur_s4_receipt.json](fold_schur_s4_receipt.json) and
[the source manifest](sources/fold-schur-s4/manifest.json) retain
the audited PF67 surface factorization and PF69/70 Seifert models.
No new external classification or uninspected PDF is imported.

The source theorem is all-degree with ONE positive volume order1,
any unramified count and the explicit global-coordinate/whole-line
hypotheses. It is not extended here to arbitrary higher volume orders
or complementary pencils. The target formula is the stated e=2
S4/S4 representation sector only, not arbitrary noncyclic pieces.
Passing complementary S4/S4, A4/S3, S3/S3, S4/C2, remaining cyclic4,
general/non-binomial keys, longer sequences and higher degrees remain
open. PF63 polynomial immersions still exist. No full Keller
exclusion, finite Keller map, independent full reproof, formal
verification, external review or novelty/priority claim is made.
