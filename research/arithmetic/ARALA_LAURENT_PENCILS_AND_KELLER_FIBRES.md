# Alternating Laurent pencils in admissible degrees and fixed Keller fibres

Rational square-discriminant Laurent pencils with full alternating monodromy
already exist in every admissible degree. Arala's explicit construction gives
them after an affine change of source and a scaling of target. Their finite
ramification is concentrated in an odd long cycle and a three-cycle. It does
not supply paired simple critical points.

Combining that construction with the fixed-map compiler gives one rational
polynomial Keller map in three variables, for each admissible degree, whose
complete fibres contain infinitely many alternating field specializations.
The ambient map has generic symmetric group. An elementary index argument
also explains why splitting the branch values of Arala's map cannot produce
the paired-simple passport while keeping the degree and pole orders.

## The known construction in Laurent coordinates

Let n>=6 be even and set

\[
 c=\frac{(n-2)^2}{n},\qquad P(X)=X^n-nX^{n-1},
 \qquad H(z)=P(z+c)-P(c),\qquad a=P(c).
 \tag{1}
\]

Then a is nonzero and

\[
 F_s(z)=H(z)-sz+a=P(z+c)-sz,\qquad
 f(z)=\frac{P(z+c)}z.
 \tag{2}
\]

The rational function f has degree n, with poles of orders n-1 and one
at infinity and zero respectively. Its two finite branch values are

\[
                     0,\qquad -n(n-2)^{n-2},                 \tag{3}
\]

with inertia types (n-1,1) and (3,1^(n-3)). Its geometric group is A_n.
Furthermore

\[
 \boxed{\operatorname{Disc}_z F_s
 =(-1)^{n/2+1}(n-1)^{n-1}s^{n-2}
                       \bigl(s+n(n-2)^{n-2}\bigr)^2.}
 \tag{4}
\]

This is Arala's family
P(X)-TQ(X), Q(X)=-nX+(n-2)^2,
under X=z+c and T=-s/n. See
[N. Arala, *Lines of polynomials with alternating Galois group*,
arXiv:2205.01502v1, Section 4, Theorems 1.2--1.3 and equation (4.4)](https://arxiv.org/pdf/2205.01502).
The [retained primary PDF](sources/arala-lines-2205.01502v1.pdf)
is the version inspected here. The Laurent coordinate formula (1), its
low Taylor coefficients and the compiler consequence below are made explicit
for comparison with the paired-critical problem; the underlying alternating
construction is prior work.

Here is a direct check of the assertions needed in this application. With
X=z+c,

\[
 z^2f'(z)=(X-c)P'(X)-P(X)
           =(n-1)X^{n-2}(X-(n-2))^2.                       \tag{5}
\]

The two displayed critical points avoid z=0, and their values give (3).
Their local degrees are n-1 and three, so all finite ramification is
accounted for. At infinity f has order n-1, and a!=0 gives its other
simple pole.

The equation F_s is irreducible over Qbar(s), since it is linear in s
and its two coefficients -z and P(z+c) are coprime. The geometric group
is transitive. Its (n-1)-cycle forces primitivity: the block containing
the fixed letter is preserved by that cycle and is either a singleton
or the entire set. The group contains a three-cycle with at least three
fixed letters. The prime-cycle theorem in
[G. A. Jones, *Primitive permutation groups containing a cycle*,
Theorem 1.1](sources/jones-primitive-cycle-1209.5169v1.pdf)
therefore makes it contain A_n. All the stated inertia permutations
are even, so it equals A_n. This also follows from Arala's Section 4.

For an independent resultant derivation of (4), let
R=(X-c)P'-P and set B=n(n-2)^(n-2). Since n is even,

\[
 \operatorname{Res}(F_s,R)
 =(n-1)^n (sc)^{n-2}
                \left(-\frac{2(n-2)}n(s+B)\right)^2.
\]

On the other hand R=(X-c)F_s'-F_s, hence
Res(F_s,R)=a Res(F_s,F_s'). Using

\[
             a=-\frac{4(n-1)}n c^{n-1}
\]

and c=(n-2)^2/n gives (4), including its sign.

If n=m^2+1 for an odd integer m>=3, (4) is the literal rational square

\[
 \boxed{\operatorname{Disc}_z F_s
 =\left[m^{n-1}s^{(n-2)/2}
                     \bigl(s+n(n-2)^{n-2}\bigr)\right]^2.}
 \tag{6}
\]

The arithmetic group is then contained in A_n and contains the geometric
A_n. Both groups equal A_n, and the splitting field is regular over Q.
This applies in particular to n=10,26,50,82,...

There are only two finite branch values. Their critical multiplicities
are n-2 and two, not the simple multiplicities required for n/2
double-transposition values. Also H need not be even; no point of the
paired even-numerator chart is obtained by (1).

## Corollary 1. Uniform alternating fibres of fixed three-dimensional maps

For every odd m>=3 and n=m^2+1 there is one polynomial map

\[
                    \mathcal F_n:\mathbb A^3_{\mathbb Q}
                                    \longrightarrow\mathbb A^3_{\mathbb Q}
\]

with Jacobian determinant one, generic degree n, and generic arithmetic
and geometric groups S_n. It has a rational target curve along which
infinitely many complete fibres are degree-n fields with splitting-field
group A_n.

More precisely, use H,a from (1) and put

\[
 A=-(3n-4)c^{n-2},\quad
 c_2=-(n-1)(n-2)c^{n-3},\quad
 c_3=-\frac{(n-1)(n-2)(n-4)}6c^{n-4}.                       \tag{7}
\]

The fixed-map theorem applies with tau=0, and its target curve is

\[
             \Gamma_n(p)=
       \left(p,\frac{c_2}{p^2},-\frac{2a}{p^8}\right),
                    \qquad p\ne0.                         \tag{8}
\]

Every separable complete fibre above Gamma_n(p) is

\[
 \boxed{
 \mathcal F_n^{-1}(\Gamma_n(p))
 \simeq\operatorname{Spec}\mathbb Q[z]/
              \bigl(H(z)-(A-p^5)z+a\bigr).}                \tag{9}
\]

Its coordinate degrees are at most 20n-54, giving 466 in degree 26.

**Proof.** The relevant import is precisely Theorem 1 of the
[fixed Keller map theorem](../jacobian/monodromy-theory/FIXED_KELLER_MAP_ALTERNATING_FAMILY.md):
it permits any degree-n polynomial H in characteristic zero once a source
translation has nonzero third Taylor coefficient. Differentiating P at c
gives (7), with c_3!=0 for every even n>=6, so tau=0 is permitted.
The constant denoted kappa there is a. Thus that theorem gives a fixed
polynomial map, its full scheme-theoretic fibre isomorphism (9), generic
S_n and the degree bound. No assumption of even H or simple critical
points occurs in the compiler theorem.

For clarity the fixed profiles are g1=1,g2=0,g3=c3 and
g_j(V)=c_j V^(2j-8) for j>=4, where c_j is the coefficient of z^j in H.
Its inverse polynomial E satisfies the full identity

\[
 p^8 E_{\Gamma_n(p)}(S)
       =H(p^3S)-(A-p^5)p^3S+a.                            \tag{10}
\]

This identifies the actual complete fibre, not a selected subalgebra.

The extension in (6) is regular A_n. The normal closure of the
reparametrization s=A-p^5 is solvable. Its intersection with the
A_n splitting field has a solvable quotient of the perfect group A_n,
and is trivial. The same argument after extending constants proves that
the pulled-back extension is regular A_n over Q(p). Proposition 2 of
the fixed-map theorem gives this precise argument. Hilbert irreducibility
as imported there now gives infinitely many rational p with full
A_n specialization after deleting the finite branch set and p=0.
Transitivity makes their root algebras fields. ∎

The generic group of the single ambient map remains S_n. This gives a
uniform supply of alternating complete fibres in dimension three;
it does not realize the Laurent extension as the generic extension of
a polynomial map, nor solve the paired-critical degree-26 problem.

## Theorem 2. A coalescence obstruction for the paired passport

Let n=2b>=8. Start with a connected degree-n genus-zero cover having
b finite branch cycles that are double transpositions, and infinity
type (n-1,1). Combine the finite branch cycles into disjoint nonempty
blocks, replacing each block by its ordered product. Suppose the
resulting tuple is again the complete branch-cycle tuple of a
degree-n genus-zero cover with the same infinity type.

Then none of its finite branch permutations is a pure three-cycle.
In particular Arala's cover is not obtainable from this paired passport
merely by coalescing finite branch values while retaining a degree-n
smooth genus-zero source and the same pole orders.

**Proof.** Define ind(sigma)=n minus the number of cycles of sigma,
including fixed letters. This is the minimal number of transpositions
whose product is sigma: each transposition changes the cycle count by
one, and a cycle of length d is a product of d-1 transpositions.
Consequently

\[
             \operatorname{ind}(\sigma\tau)
             \le \operatorname{ind}(\sigma)
                         +\operatorname{ind}(\tau).        \tag{11}
\]

Each original finite permutation has index two, giving total finite
index 2b=n. A genus-zero cover with the stated infinity permutation
must still have finite index n by Riemann--Hurwitz. Hence (11)
is equality separately in every block: a block of k double
transpositions has product index exactly 2k.

A pure three-cycle has index two. Its block would therefore have k=1.
A single double transposition is not a three-cycle, a contradiction. ∎

This is a restriction on branch-value coalescence with the stated
degree and smooth-source conditions. It does not exclude degenerations
with additional source components, changes of degree or pole passport,
or a different construction of rational paired-critical covers.


## Theorem 3. Uniform Hasse counterexamples from the concentrated passport

Let n=m^2+1 with m>=3 odd, and b=n/2. Write

\[
 H_k=(S_k\times S_{n-k})\cap A_n\quad(1\le k<b,\ k\ {\rm odd}),\qquad
 W_b=(S_b\wr S_2)\cap A_n.
 \tag{12}
\]

Here W_b stabilizes an unordered partition into two b-element subsets.
There are infinitely many pairwise linearly disjoint A_n extensions
K_j/Q for which

\[
 \mathcal A_j=
 \left(\prod_{\substack{1\le k<b\\ k\ {\rm odd}}}K_j^{H_k}\right)
                   \times K_j^{W_b}                       \tag{13}
\]

has rank 2^(n-2), has a point over every completion of Q, and has no
Q-point. The extensions can be obtained by specializing Arala's one
fixed cover along a rational arithmetic progression.

In degree ten this gives three field factors of degrees 10,120,126
and rank 256. This total rank is minimal among everywhere locally
soluble, globally pointless finite etale algebras whose splitting
field has group A10. The minimality imports the classified subgroup
calculation already established in Proposition 4.2 of the
[alternating subgroup-cover theorem](ALTERNATING_SUBGROUP_COVERS_AND_HASSE.md).
It is not a new classification proof.

**Proof of the uniform cover.** Since m is odd, b is odd.
Every even permutation fixes a coset of a subgroup in (12).
If it has an odd cycle of length d<b, use H_d. If d>b, its complement
has odd size n-d<b. If d=b, its support gives the partition fixed by
W_b. These cases exhaust permutations having an odd cycle. If all
cycles have even length, color the letters alternately on each
cycle. The two color classes have size b, and the permutation
interchanges them, so it belongs to a conjugate of W_b.
Thus every cyclic subgroup is contained in a conjugate of a listed
subgroup. Their coset degrees sum to

\[
 \sum_{\substack{k<b\\k\ {\rm odd}}}\binom nk+\frac12\binom nb
                    =2^{n-2},                             \tag{14}
\]

by pairing complementary odd subsets.

The same cover contains every actual ramified decomposition group
at all sufficiently good primes of an Arala specialization.
Indeed, suppose a nontrivial inertia group I is a subgroup of a
cyclic group generated by an (n-1)-cycle. Every nonidentity power
of that cycle still moves all its n-1 letters. Consequently the
common fixed set of I is its unique remaining letter. Since
I is normal in the decomposition group D, that letter is D-invariant,
and D lies in H1. If instead I is nontrivial three-cycle inertia,
its support has size three and is D-invariant, giving H3.
Both subgroup classes occur in (12). If I is trivial, D is cyclic
and the preceding element-cover argument applies. The same argument
handles the (n-1)-cycle inertia at infinity.

Apply [Conrad, *Inertia groups and fibers*, Theorem 0.2,
printed page 2](sources/conrad-inertia-groups-fibers-1999.pdf)
to the regular A_n Galois closure. The source of that Galois cover
is a smooth projective curve in characteristic zero, and the target
is the regular curve P1. The morphism is finite flat and tame
along its branch divisor, so the theorem applies. It supplies a
fixed finite set S of rational primes such that every specialization
inertia outside S is a subgroup of a geometric cyclic inertia group.
Enlarge S to contain any desired presentation and denominator
exceptions. The preceding fixed-set argument proves local
solubility outside S, including every ramified prime.

All primes in S can be controlled simultaneously without imposing
a density assumption on an arbitrary local open set.
Use the original parameter T and rational source point X=1:

\[
 T_0=-\frac1{n-4},\qquad
 (P-T_0Q)(1)=0,\qquad
 (P'-T_0Q')(1)=-\frac{n(n-3)^2}{n-4}\ne0.                  \tag{15}
\]

For p in S put e_p=v_p(n(n-3)^2) and choose

\[
       M_0=\prod_{p\in S}p^{2e_p+1},\qquad
                 T=T_0+M_0t,\quad t\in\mathbb Z.           \tag{16}
\]

The integer polynomial

\[
 G_t(X)=(n-4)(P(X)-TQ(X))
\]

satisfies

\[
 G_t(1)=-M_0t(n-1)(n-4)^2,\qquad
 G_t'(1)=-n(n-3)^2+n(n-4)M_0t.                            \tag{17}
\]

Thus v_p(G_t(1))>2e_p and v_p(G_t'(1))=e_p, except that
G_0(1)=0 already gives the required root. The strict Hensel
inequality gives a Q_p-root of the degree-n polynomial for every
integer t and every p in S. Its natural permutation action
therefore fixes a sheet, so D_p lies in H1. This calculation
allows denominators in T0 and does not require the polynomial
P-TQ itself to be p-integral.

The affine reparametrization (16) preserves the regular A_n
extension. Integral Hilbert irreducibility in the precise form of
[Zywina, *Hilbert's irreducibility theorem and the larger sieve*,
Theorem 1.2, printed page 2](sources/zywina-quantitative-hit-2010.pdf)
gives infinitely many integer t with full group after deleting
the finite branch exceptions. To obtain linearly disjoint fields
inductively, take the compositum E of all previously chosen
Galois fields. Regularity makes the generic compositum have
group A_n x Gal(E/Q). Apply the same arithmetic Hilbert theorem
over Q to a primitive element of that compositum, with t still
an integer. Its full product-group specializations give the next
A_n field linearly disjoint from E. No application of Hilbert
irreducibility over E with an uncontrolled parameter field is used.

At the real place the decomposition group is cyclic of order
at most two. We have now covered every local decomposition group.
The standard fixed-coset criterion identifies a Q_v-point of
Spec(K_j^H) with a fixed coset for D_v on A_n/H. Hence (13)
has a point over every completion. All factors have degree
greater than one, so none admits a Q-algebra homomorphism to Q;
the product spectrum has no rational point. Each subgroup is
proper, and simplicity of A_n makes its core trivial, so each
factor has normal closure K_j. This also makes the algebras
pairwise nonisomorphic. Formula (14) proves the rank. ∎

For the degree-ten minimum, any finite etale algebra split by an
A10 field is a product of fixed fields of subgroups. Absence of a
Q-point makes each subgroup proper. By
[Chebotarev's theorem, in Lenstra's notes, Theorem 3.1, page 3](sources/lenstra-chebotarev-density-notes.pdf),
every conjugacy class occurs as an unramified Frobenius.
Local solubility at all such primes forces those subgroup classes
to cover every cyclic subgroup. The minimum 256 in the cited
Proposition 4.2 therefore applies to every such algebra. The
construction attains it.

The degree-ten rank-256 conclusion uses the concentrated
(n-1)/(3) passport. It does not change the decomposition-group
obstruction for rank 256 attached to the earlier double-transposition
pencil. Those are different covers.

## Corollary 4. One fixed map containing the Hasse family

For each admissible n as in Theorem 3, there is one polynomial
Jacobian-one map of A4_Q, of generic degree R=2^(n-2) and generic
arithmetic and geometric group S_R, whose complete fibres along a
rational target curve include infinitely many of the algebras (13).
For n=10 the generic degree is 256.

**Proof.** Over Q(T) form the product of the fixed subfields in (13)
inside the regular Arala splitting field. This is a finite etale
Q(T)-algebra of rank R. Corollary 2 of the
[fixed finite-etale-family theorem](../jacobian/monodromy-theory/FIXED_KELLER_MAP_FINITE_ETALE_FAMILIES.md)
makes it monogenic and realizes its separable specializations
over a nonempty parameter open set as whole fibres of one fixed
map in four variables. Its hypotheses require only rank at least
three and generic separability, both satisfied here.

After deleting the finite presentation exceptions, the full-group
specializations in Theorem 3 identify the specialized factors with
K_j^H for the corresponding H. They survive this finite deletion.
The whole-fibre isomorphism preserves local solubility and absence
of rational points. The generic symmetric group follows from the
same compiler theorem. ∎

Corollary 4 uses monogenicity for the fixed four-dimensional compiler.
The following argument constructs a particular monogenic presentation
that also gives a fixed three-dimensional realization of the product.
The degree-\(n\) fibres in Corollary 1 are the individual root fields.

## Corollary 5. One fixed threefold for the concentrated Hasse family

For every admissible \(n=m^2+1\) with \(m\geq3\) odd, the product family
in Theorem 3 occurs along a rational target curve as complete
scheme-theoretic fibres of a single polynomial Jacobian-one map
\(\mathbb A^3_{\mathbb Q}\to\mathbb A^3_{\mathbb Q}\). Its generic
degree is \(2^{n-2}\) and its generic arithmetic and geometric groups
are \(S_{2^{n-2}}\). The curve uses the original Laurent parameter \(s\)
of (2), on a nonempty open subset. In particular the degree-ten
case gives infinitely many everywhere locally soluble but globally
pointless rank-256 fibres, with the \(A_{10}\) splitting fields
of Theorem 3.

Let \(X_1,\ldots,X_n\) be the roots over \(\mathbb Q(s)\) in the original
\(X=z+c\) coordinate, so \(\sum_i X_i=n\). Take the following explicit
resolvent roots:
\[
 \frac{1}{n}\sum_{i\in I}X_i
       \quad (|I|=k<b,\ k\ \text{odd}),\qquad
 \left(\sum_{i\in J}X_i-\frac n2\right)^2-\frac{n^2}{4}
       \quad (\{J,J^c\}\ \text{a middle partition}).
\]
The latter expression is unchanged by replacing \(J\) by \(J^c\).
The product \(\mathcal P(s,T)\) of \(T\) minus these expressions
is a monic polynomial over \(\mathbb Q[s]\): all coefficients are
symmetric polynomials in the \(X_i\), and their elementary symmetric
polynomials belong to \(\mathbb Q[s]\).

To check all orbit sizes, no linear combination \(\sum_i u_i X_i\)
belongs to the base field for a nonconstant \(u\in\mathbb Q^n\).
Subtract the mean of the coefficients and average over \(A_n\).
This reduces to a zero-sum coefficient vector giving zero.
For an index \(i\) with \(u_i\ne0\), averaging this relation over
the stabilizer of \(i\) gives
\(u_i(X_i-(n-X_i)/(n-1))=0\), since it is transitive on the other
labels. This would force \(X_i=1\), contrary to the nonconstant
irreducible polynomial for \(X_i\). Thus different indicated subsets
give different sums, also when their sizes differ. An equality
of two middle-partition expressions forces the centered sums to
agree up to sign, and the same argument yields identical or
complementary subsets. These expressions therefore have respective
stabilizers \(H_k\) and \(W_b\). Stabilizers from the two types cannot
coincide, since \(W_b\) is transitive on the original labels and
\(H_k\) is intransitive. For the transitivity of \(W_b\), any
permutation of the labels preserving the unordered partition can
have its parity changed by a transposition within a block while
fixing a prescribed image label. This is possible since \(b\geq5\).
Consequently \(\mathbb Q(s)[T]/(\mathcal P)\) is the full product
of fixed fields in Theorem 3, and \(\mathcal P\) is generically separable.

At \(s=0\) the original root multiset is \(0^{n-1},n\).
All middle-partition expressions specialize to zero; each odd-subset
expression specializes to one if and only if its subset contains
the unique label of \(n\). Therefore
\[
 \mathcal P(0,T)=T^A(T-1)^B,\qquad
 A=\sum_{\substack{k<b\\k\ {\rm odd}}}\binom{n-1}{k}
          +\tfrac12\binom nb,\qquad
 B=\sum_{\substack{k<b\\k\ {\rm odd}}}\binom{n-1}{k-1}.
\]
In degree ten, \((A,B)=(219,37)\). In every admissible degree,
\(A>B>0\) since \(k<b=n/2\). Set \(\tau_0=A/(A+B)\).
The logarithmic derivative of \(T^A(T-1)^B\) at \(\tau_0\) is
zero. Its derivative there is
\(-A/\tau_0^2-B/(\tau_0-1)^2\ne0\), and its second derivative
is \(2A/\tau_0^3+2B/(\tau_0-1)^3\ne0\) by \(A\ne B\).
As \(\mathcal P(0,\tau_0)\ne0\), these give a simple critical point
with nonzero third derivative.

Apply [Theorem 4 of the fixed-family compiler](../jacobian/monodromy-theory/FIXED_KELLER_MAP_FINITE_ETALE_FAMILIES.md#theorem-4-a-critical-section-gives-a-fixed-compiler-in-dimension-three)
to \(\mathcal P\). Its finite exclusions discard only finitely many
members of the original progression of Theorem 3, using \(s=-nT\)
to translate the progression parameter \(T\). For each retained
full-group specialization, the resolvent factors specialize to
those same fixed fields, as in Corollary 4. Their complete local
and global point properties transfer by the scheme-theoretic
fibre identity. This construction leaves the paired-critical
double-transposition product problem separate.

## Evidence and scope

The [checker](verify_arala_laurent_keller.py) recomputes full exact
discriminants for n=6,8,10,26, checks every coefficient of (10) for those
degrees and n=50, verifies the low Taylor coefficients, and checks all
ordered pairs of double transpositions in S6 as finite controls for the
index-loss mechanism. It also checks all even cycle types in degrees ten
and twenty-six, the rank identity, every nontrivial long-cycle power,
and the exact integer Hensel identities. The uniform proofs are written above.
The [receipt](arala_laurent_keller_verification.json) and
[source ledger](arala_laurent_keller_source_ledger.json) give exact inputs
and replay scope. The finite controls are not proofs of the quantified
claims. No paired-simple rational degree-26 polynomial, formal
verification or external human review is asserted.
