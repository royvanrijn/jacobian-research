# Universal regular-coordinate affine classification and surgery end

14 November 2026. PF2RULED43 closes the affine-converse gap left explicitly
open in PF2RULED40 and PF2RULED42. It treats every Laurent support and every
characteristic gcd chain. Its geometric statements do not depend on a finite
Keller field degree. It does not exclude the remaining multi-characteristic
quartic models.

## 1. Complete marked affine classification

Work over C. Mark U=A2_(x,y) inside a smooth affine N with sole complement
R=A1 and div(dx wedge dy)=R. Suppose a source coordinate is regular on N.
It is constant on R: otherwise its restriction is a nonconstant polynomial
A1 -> A1, and a transverse component of a suitable coordinate fiber through
R would extend the source coordinate line, although that line is already
closed in N. Shift the constant to zero. Write v_R(y)=r and v_R(x)=-q.
The r=1 case is PF2RULED38-39. Here r>=2.

PF2RULED40 gives the necessary root-chart form

    y=z^r, x=z^(2-r)*t_i+h_i(z),
    h_i(z)=h(zeta^i*z), h(z)=sum_(j=r-1)^q c_j*z^(-j),
    c_q!=0, q>=r-1, gcd(r,{j:c_j!=0})=1.                (1)

There are r boundary charts A2_(z,t_i) and the source chart A2_(z,x),
all identified over z!=0. The distinct boundary components are cycled
by the root deck action. Conversely, EVERY datum (1) gives a smooth
affine N with precisely these properties. No extra condition on mixed
Laurent coefficients is required for affineness.

Proof of the converse. Regard the root gluing as a line-bundle torsor over
the affine line with r+1 origins. Its weights in Dubouloz's notation are
mu_source=0 and mu_i=r-2. Section2.8 of his *Danielewski-Fieseler surfaces*
states that the glued torsor is affine iff both directional translation
functions have a pole at zero for every pair of origins. The source/boundary
translations are h_i and -z^(r-2)h_i: their pole orders are q and q-r+2,
both positive. For two boundary charts the translation is
z^(r-2)(h_j-h_i). Support gcd one says h_i!=h_j for every i!=j.
Its largest differing support index is at least r-1, so this translation
and its negative have a pole. Thus the criterion applies to all pairs,
including the source origin. This is an affine theorem, not merely a
separatedness test. Smoothness follows from the plane charts.

The action z -> zeta*z, x fixed, sends chart i to chart i+1 and
rescales t by zeta^(r-2), with the indexing direction adjusted if desired.
It acts freely near the r boundary lines and permutes them transitively.
On the source chart it is the reflection (x,z)->(x,zeta*z). The finite
affine quotient is therefore smooth, with source A2_(x,y=z^r) and one
complement R=A1. In a boundary chart,

    dx wedge dy = r*z*dt_i wedge dz.                    (2)

Hence its divisor downstairs is exactly R. This proves sufficiency and,
with PF2RULED40, the full marked classification in this sector.

The standard volume is globally exact iff c_r=0. This is the common
residue theorem of PF2RULED40, now applied to every actual affine datum,
not to hypothetical gluing schemes. Consequently q=r is impossible for
an exact-volume model. Source translation extends as the regular locally
nilpotent field delta=partial_x; its boundary expression is
z^(r-2)*partial_(t_i). Its zero order on R is r-2.

Marked equivalence means an isomorphism preserving the specified source
open and omitted line, without fixing its source coordinates pointwise.
For r>=2 the classification of the data is also explicit: r must agree,
and, for some a,b in C*,

    c'_j=b*a^j*c_j for all j.                           (3)

It is induced by y'=a^r*y, x'=b*x+p(y), where p is any polynomial.
Indeed the boundary substitution is

    z'=a*z,
    t'=b*a^(r-2)*t+a^(r-2)*z^(r-2)*p(z^r),             (4)

which is polynomial, as is its inverse. Conversely PF2RULED38's
regular-coordinate uniqueness forces every regular source coordinate
to be affine in y; its boundary value fixes the translation. A plane
automorphism over that coordinate has x'=b*x+p(y). Comparing the unique
negative Laurent parts gives (3). The r=1 restriction on the constant
term of p is different and remains exactly as in PF2RULED38. This is a
marked classification; no assertion about arbitrary unmarked isomorphisms
or preservation of a finite target algebra is included.

## 2. The support determines all characteristic exponents

At the original compactification corner use w=1/x and y. Let

    Gamma: w=1/h(z), y=z^r.                             (5)

This is a convergent primitive plane branch: the cyclic extension
C((z))/C((z^r)) has no nonidentity element fixing h, by the support gcd.
Gamma is a SOURCE reference branch. It is distinct from the nonproper
TARGET image curve of any hypothetical finite Keller extension.

For a=1,...,r-1 put

    j(a)=max{j:c_j!=0 and r does not divide a*j},
    b(a)=2q-j(a).                                      (6)

The exact identity for the difference of two reciprocals gives

    ord_z(1/h(z)-1/h(zeta^a*z))=b(a).                  (7)

Its numerator has a unique highest nonzero differing Laurent power;
its denominator has pole order 2q. Thus there is no genericity
assumption on the nonzero coefficients.

Read the support in descending order. Start g_0=r and retain j_i
exactly when g_i=gcd(g_(i-1),j_i) strictly decreases. Stop at g_s=1.
Then, after subtracting integral terms in y where needed, the
characteristic exponents and their contact multiplicities are

    beta_i=2q-j_i, multiplicity g_(i-1)-g_i.
    S=sum_i (g_(i-1)-g_i)*beta_i=sum_a b(a).
    c=2*delta(Gamma)=S-r+1.                            (8)

For the exponent assertion, at each gcd stage all higher support indices
and q are divisible by the current gcd. Inverting the Laurent series
therefore has only terms divisible by that gcd until the first differing
index j_i, whose coefficient is -c_(j_i)/c_q^2. It cannot cancel.
For (8), factor the reduced Weierstrass equation in w into its r conjugate
roots. Its w derivative on Gamma has order S. Teissier's polar identity
for the smooth coordinate y gives S=mu+r-1, and the irreducible-branch
conductor formula gives mu=2*delta in characteristic zero. Precise primary
imports are pinned in the source manifest; no positive-characteristic
version is used.

For explicit knot conventions set p_i=g_(i-1)/g_i,
alpha_1=beta_1/g_1 and, for i>1,

    alpha_i=p_i*p_(i-1)*alpha_(i-1)
                    +(beta_i-beta_(i-1))/g_i.           (9)

The reference knot is the positive iterated torus knot with successive
pairs (p_i,alpha_i). In the smooth borderline case a first T(r,1)
represents the unknot. For s>=2 the strict algebraic cabling inequality
alpha_i>p_i*p_(i-1)*alpha_(i-1) is automatic. Its last cabling fiber
framing is

    F=p_s*alpha_s=S+beta_s.                            (10)

One can prove (10) inductively from (9): at stage i its framing is
(S_i+g_i*beta_i)/g_i^2, with S_i the partial contact sum. Equivalently,
two curvettes differing first in their last characteristic coefficient
have intersection S+beta_s: r-1 conjugate contacts contribute S and
the remaining matching contact contributes beta_s. This identifies the
framing in the actual embedded-resolution graph, rather than just its
absolute homological determinant.

## 3. Universal free tail and surgery coefficient

Let j_*=j_s and beta_*=beta_s. After the characteristic resolution the
whole omitted chart has exactly

    kappa=j_*-r+2=2q-r+2-beta_* >=1                    (11)

additional free blowups. Indeed the actual root chart gives

    w=1/(h(z)+z^(2-r)*t)
     =1/h(z)-c_q^(-2)*z^(2q-r+2)*t+higher terms.        (12)

The last characteristic curvette parameter first varies at order beta_*.
Its difference from (12) is precisely kappa free smooth contacts.
The coefficient of t in the final free coordinate is nonzero. Hence
that coordinate covers the entire omitted A1, not just a dense open
or a formal valuation. Finitely many reference coefficients determine
all blowup centers. As a discrepancy check, at the last characteristic
vertex v(dw wedge dy)=r+beta_*-1; each free blowup adds one, while
w^-2 subtracts 2q. Requiring v(dx wedge dy)=1 gives exactly (11).

For completeness retain the original two deleted boundary curves.
The original w=0 curve detaches from Gamma at its initial integral or
first nonintegral contact, before the last characteristic vertex. It
meets one earlier exceptional divisor transversely. The other original
deleted curve is its zero-weight leaf. Neumann's retained R6 in its
punctured form removes this pair with one remaining connected branch,
one attachment and no extra summands, also when the attachment is at
an internal vertex. Neither the free tail nor its omitted puncture
is removed. This is the same move as PF2RULED42, now with an arbitrary
characteristic-resolution tree beyond that attachment.

The exceptional tree C_0 of the characteristic blowups has boundary S3,
is unimodular, and its last curvette is Gamma. Its last marked inverse
entry is F, the linking of two curvettes computed in (10). Performing
kappa free blowups, omitting the last exceptional R, decreases that
vertex weight by one and adds kappa-1 deleted vertices of weight -2,
with the puncture at the free endpoint. Write mu,lambda for the
meridian and zero-linking longitude of Gamma and h=F*mu+lambda for its
last cabling fiber. Continued elimination of this free chain gives

    mu_R=kappa*(mu+h)-(kappa-1)*h
         =(F+kappa)*mu+lambda.                         (13)

Thus the compactification PLUMBING END PAIR is

    n integer surgery on Gamma, omitted knot its surgery dual,
    n=F+kappa=2*delta(Gamma)+2q+1.                      (14)

A simultaneous orientation reversal is harmless; the cabling signs and
all slopes must otherwise stay consistent. Schur complementation of
C_0 and the -2 tail gives determinant n and marked cofactor n-1.
Therefore H1=Z/n, n odd >=3, and the omitted class is primitive. A
primitive surgery-dual class also follows directly from numerator n
and denominator one in (13).

The stable peripheral end comparison already used in PF2RULED36-42
transfers this H1 and marked class to an actual target-radial section
when a finite map exists. This suffices for Section5. We do NOT claim
here a homeomorphism classification of every multicharacteristic
radial section from its abstract group, or import its model JSJ tori
into a radial section without a recognition argument. PF2RULED42
already proves the needed recognition in its single-pair sector.

## 4. A universal gap in the terminal filling, and the distance-one case

One entire terminal possibility disappears:

    kappa != 2.                                       (15)

Otherwise j_*=r, but an index divisible by r cannot decrease any
current gcd, all of which divide r. Consequently the terminal filling
is either kappa=1 or kappa>=3, for every admissible support, regardless
of field degree or coefficient sizes. Moreover

    kappa=1 iff c_(r-1)!=0 and
       gcd(r,{j>r-1:c_j!=0})>1.                        (16)

The empty higher-support set has gcd r. Thus distance one means an
essential lowest term z^(-(r-1)) breaks the remaining higher-support
symmetry, rather than an arbitrary fine tuning of coefficients.

For a genuine multicharacteristic reference branch, let
p=p_s=g_(s-1)>1 and F_prev=p_(s-1)*alpha_(s-1). When kappa=1,
the preceding drop index is divisible by p, exceeds r-1, and cannot
equal r. Hence j_(s-1)>=r+p. Equations (9)-(10) give

    alpha_s-p*F_prev=beta_s-beta_(s-1)>=p+1,
    n=p*alpha_s+1,
    n-p^2*F_prev >= p*(p+1)+1 >=7.                    (17)

Ito's cable-surgery formula therefore collapses the outer cable space
to companion filling n/p^2, in lowest terms. Its distance from the
companion cabling fiber is the positive numerator n-p^2*F_prev,
at least seven. There is no consecutive distance-one collapse.

When kappa>=3, the surgery model instead retains the terminal Seifert
piece over a disk with two exceptional orders p_s and kappa, attached
to the companion exterior. Its dual knot is the kappa-order filling
core. These are genuine model JSJ data by the primary cable-surgery
formula; they are not atoroidal small-Seifert data. The existing
quartic atoroidality theorem cannot be applied to this sector on the
basis of (14) or cyclic H1 alone.

Examples, with h the sum of the indicated Laurent monomials:

| r | q | support | characteristic beta | 2 delta | F | kappa | n |
|---|---|---|---|---|---|---|---|
|4|6|6,5|6,7|16|26|3|29|
|4|6|6,3|6,9|18|30|1|31|
|6|10|10,8,5|10,15|50|70|1|71|
|6|12|12,10,5|14,19|70|94|1|95|
|8|12|12,10,9|12,14,15|84|106|3|109|
|8|12|12,10,7|12,14,17|86|110|1|111|

Every example is an actual smooth affine exact-volume model by Section1
(c_r=0), with no claimed finite Keller target. In the second example
Gamma is the (2,15)-cable of the trefoil and the surgery collapses to
31/4 trefoil filling, at distance 7 from its regular fiber. In the first,
Gamma is the (2,13)-cable of the trefoil and the order-3 terminal filling
piece remains. These small controls illustrate the two universal cases;
no growing coefficient search is being proposed.

## 5. All-degree consequence and the remaining finite-map gate

The whole regular-coordinate geometric sector, including all
multicharacteristic supports, now has finite odd cyclic end homology
with nonzero primitive omitted class. For r=1 this is PF2RULED39's
order 2q+1; for r>=2 it is (14). Apply PF2RULED34's general topological
lemma: a connected simple S_d branched cover of a torus knot or unknot,
d>=3, with finite source H1 has zero ramification class. Thus in EVERY
finite field degree allowed by the one-index-two-prime Keller setup,

    target boundary infinity is neither torus nor unknot;
    target boundary-normalization degrees cannot be coprime;
    each nonconstant target-coordinate boundary degree is composite >=4.

The last assertion uses the same elementary degree-lowering shears as
PF2RULED34. These are target boundary degrees, not the source r,q or the
characteristic pairs of Gamma. The marked source class, rather than an
assumed torus type of the source, is the hypothesis needed for this
all-degree obstruction. The d<=2 cases are already excluded by the
retained finite-plane degree/branch inputs; no new degree-two assertion
is inferred from the S_d lemma.

PF2RULED42's quartic frontier stays composite r>=4, q>r, gcd(q,r)>1,
q>=r+ell(r), and at least two source characteristic pairs. The affine
existence and end-presentation gates are now closed for every support.
What remains is a global comparison with the actual finite target
algebra and its residual/peripheral data. In a quartic JSJ approach,
prove that the ramification fiber slopes in the S4 cabling census
match the source surgery model, handle product pieces and the
kappa=1 collapse, and justify radial recognition before using source
JSJ incompressibility. Equations (15)-(17) are available necessary
structure, not a proof that every remaining quartic model is impossible.
The norm/conductor conditions, exactness c_r=0 and all omitted primes
remain essential; deleting one generally loses finiteness. The
no-regular-coordinate case and higher finite degrees remain separate.

## 6. Sources, exact controls and assurance

The new primary imports are Dubouloz Section2.8 (both directional
translation poles), Garcia Barroso--Ploski Proposition2.1(iii)-(iv) and
Corollary3.4 (conductor and characteristic-zero polar formula), and
Lokteva Proposition5 and the proof of Proposition6/Figures3-4
(characteristic cabling and the surgery tail). The retained primary
Neumann and Ito statements are used with their marked/peripheral data.
Lokteva Proposition6 is stated for n>=F+2; the kappa=1 calculation
here is our marked-chain/cable-filling argument, not an application
of that proposition outside its range. Their exact bytes, page audits
and restricted uses are in
[sources/universal-regular-end/manifest.json](sources/universal-regular-end/manifest.json).
No novelty priority is asserted.

Run [verify_universal_regular_end.py](verify_universal_regular_end.py).
Its [receipt](universal_regular_end_receipt.json) checks exact cyclic
contacts, all transition-pole directions, characteristic recurrences,
independent numerical-semigroup Apery conductors, actual reference
resultants, determinant/cofactor identities and the distance-one gap.
Finite controls support the written universal proofs; they do not
implement the full affine theorem or three-manifold argument.
Formal verification, independent full reproof and external review
remain false. No remaining whole quartic-sector exclusion, finite
Keller construction or full-goal completion is claimed.
