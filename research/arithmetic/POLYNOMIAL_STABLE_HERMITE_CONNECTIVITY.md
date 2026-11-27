# A simple-branch threshold for arbitrary polynomial passports

Consider nonidentity finite factors with fixed product
\(P_N=(0\,1\,\ldots\,N-1)\) and total permutation index \(N-1\).
Separate all literal transpositions into \(s\ge1\) factors \(T\).
Call the remaining \(k\ge1\) factors \(g_1,\ldots,g_k\). Here
each \(g_i\) has index at least two; it may have any disconnected
cycle type. Write its nontrivial cycle lengths as
\[
 q_{i1},\ldots,q_{ip_i}\ge2,\qquad
 p=\sum_i p_i,\qquad I=\sum_{i,j}(q_{ij}-1),\qquad N=I+s+1.
                                                               \tag{1}
\]
The prescribed partition of factor \(i\) includes
\(N-\sum_j q_{ij}\) fixed points, which must be nonnegative.
Full braids allow every placement of the prescribed cycle types;
absolute equivalence is simultaneous conjugacy.

The main result is one full Hurwitz component whenever
\[
                             s\ge p-k.                         \tag{2}
\]
With at least two non-transposition factors, its actual monodromy
is \(S_N\). A single such factor has precisely the cyclic power
obstruction described in K3. The degree-six passport
\(\mathrm{DT},\mathrm{DT},T\), with \(p-k=2,s=1\), has two components.
Thus the bound cannot be lowered uniformly. This is a sufficient
range, not a classification of the passports below it.

## Primary comparison and the existence input

Mullane, [published Theorem 7.11, p1413](sources/mullane-second-kind-hurwitz-2022.pdf),
proves genus-zero Hurwitz irreducibility with at least \(N-3\)
transpositions and at most one branch partition having more than
one ramified point. Its purity hypothesis includes the branch at
infinity. It does not cover arbitrary multiple disconnected finite
factors. The proof below instead uses their total number \(p\) of
nontrivial cycles minus the number \(k\) of their branch values,
and polynomial Hermite interpolation. No novelty
claim is made.

Existence and the total number of literal tuples use the
Goulden--Jackson cactus formula, [Section 3, pp361--363](sources/goulden-jackson-cacti-1992.pdf).
For partitions \(\lambda_1,\ldots,\lambda_r\) of \(N\), put
\(c_i=\ell(\lambda_i)\) and let \(m_{i,a}\) be their part multiplicities.
Under \(\sum_i(N-c_i)=N-1\), the fixed ordered-pattern count is
\[
 A(\lambda_1,\ldots,\lambda_r)
     =N^{r-1}\prod_{i=1}^r
                  \frac{(c_i-1)!}{\prod_{a\ge1}m_{i,a}!}>0.     \tag{3}
\]
There is a printed slip in Theorems 3.1 and 3.2 of the retained
1992 paper: the cycle-count hypothesis is written \(N+1\).
For \(r\) factors it must be \((r-1)N+1\), as stated in the last
line of its proof on p363. The displayed formula and proof use
this corrected condition. In particular, the proof's final
coefficient calculation is
\[
 \prod_i\frac{(c_i-1)!}{\prod_a m_{i,a}!}
 \left\{N^r-N^{r-1}\sum_i\sum_a(a-1)m_{i,a}\right\}.
\]
The index sum \(N-1\) gives (3). This is also the formula stated
for valid passports in Jones--Zvonkin
[Section 1.3, retained p133](sources/jones-zvonkin-braid-cacti-2002.pdf).
The actual original statements, the corrected proof, and this
comparison were inspected. Positivity, together with the polynomial
Riemann-existence correspondence, supplies nonemptiness. This
count is not a component theorem.

Explicitly, append the inverse full cycle as the monodromy at
infinity. The factors generate a transitive group because their
product is a full cycle. Riemann existence gives a connected
cover, and Riemann--Hurwitz gives genus zero: the total index
including infinity is \(2N-2\). The cover has a unique point over
infinity, so source and target coordinates make it a polynomial.

## Theorem K1. Exact coefficient spaces, including the boundary

Give the \(p\) nontrivial critical points separate labels \(u_{ij}\).
Require them to be pairwise distinct. Let \(v_i\) be the value of
the \(i\)-th non-transposition branch, and put
\[
 Q(z)=\prod_{i,j}(z-u_{ij})^{q_{ij}},\quad
 Z(z)=\prod_{i,j}(z-u_{ij})^{q_{ij}-1},\quad
 M=\deg Q=I+p.
                                                               \tag{4}
\]
Let \(\Phi\) be the unique polynomial of degree less than \(M\)
with
\[
           \Phi\equiv v_i\pmod{(z-u_{ij})^{q_{ij}}}.             \tag{5}
\]
Its coefficients are rational functions of the \(u_{ij},v_i\),
with denominators nonzero on the distinct-point configuration.
All monic degree-\(N\) polynomials with these congruences are as follows.

If \(s\ge p-1\), then \(d=N-M=s+1-p\ge0\), and
\[
                  f=\Phi+QR,\qquad R\ \hbox{monic},\quad
                                     \deg R=d.                 \tag{6}
\]
Degree zero means \(R=1\).

If \(s=p-2\), then \(N=M-1\), and
\[
                       f=\Phi,\qquad L_u(v)=1,                 \tag{7}
\]
where \(L_u(v)=[z^{M-1}]\Phi\) is a linear form in the \(k\)
values. Only configurations with \(L_u\ne0\) occur.

More generally, if \(N<M\), set
\[
 h=M-N-1=p-s-2,\qquad
 F_{bi}(u)=\sum_j\operatorname{res}_{z=u_{ij}}
                                      \frac{z^b\,dz}{Q(z)}
 \quad(0\le b\le h),\qquad e=(0,\ldots,0,1)^t.
                                                               \tag{M1}
\]
Then the complete coefficient space is
\[
                         f=\Phi,\qquad F(u)v=e.               \tag{M2}
\]
In particular (7) is the case \(h=0\). All rows of \(F\)
annihilate the constant value vector \((1,\ldots,1)\).
No assumption that \(F\) has full rank is imposed in (M2).

In either case put
\[
 A(z)=\frac{f'(z)}{NZ(z)},\qquad
 C(Y)=\operatorname{Norm}_{\mathbb C[z]/(A)}(Y-f(z)).             \tag{8}
\]
Then \(A\) is monic of degree \(s\). Necessary and sufficient
additional conditions for the exact passport are
\[
 \begin{gathered}
 v_i\ne v_h\quad(i\ne h),\qquad
 \gcd(A,A')=1,\qquad
 \gcd\!\left(A,\prod_{i,j}(z-u_{ij})\right)=1,\\
 C(v_i)\ne0\quad(1\le i\le k),\qquad \operatorname{disc}(C)\ne0.
 \end{gathered}                                                 \tag{9}
\]
Thus the exact coefficient locus is a nonempty algebraic open
subset of the parameter space (6) or (M2), and includes every
polynomial in this passport after choosing critical-point labels
and making it monic.

**Proof.** The pairwise coprime powers in (5) give the unique
Chinese-remainder interpolant. Division by \(Q\) gives (6).
When \(N=M-1\), no nonzero multiple of \(Q\) can occur and monicity
is precisely (7). The congruences imply \(Z\mid f'\); the leading
term and (1) give the degree assertion for \(A\).

For (M2), the sum of residues of \(z^b\Phi/Q\) at the finite
poles equals the coefficient of \(z^{-1}\) in its expansion
at infinity. At a pole of color \(i\), replacing \(\Phi\) by
\(v_i\) does not change that residue. Consequently
\[
                 F_bv=[z^{-1}]\,z^b\Phi/Q.                   \tag{M3}
\]
These equations are triangular in the coefficients of
\(\Phi\): the first \(h\) set its coefficients of degrees
\(M-1,\ldots,N+1\) to zero, and the last sets its coefficient
of degree \(N\) to one. This proves (M2). Since \(b\le h<M-1\),
the same coefficient for \(z^b/Q\) is zero, proving the
constant-vector assertion.

The second gcd in (9) makes every prescribed ramification index
exact. The first makes the remaining critical points simple.
The roots of \(C\), with multiplicity, are their critical values.
Its discriminant and the evaluations \(C(v_i)\) give precisely
the remaining value exclusions. There are \(I+s=N-1\) derivative
zeros with multiplicity, so none have been omitted. An additional
multiple preimage of any \(v_i\) would be one of these forbidden
critical points; consequently all other preimages in that fibre
are simple. These facts prove both directions.

All exclusions are algebraic open conditions. Formula (3) and
Riemann existence give an actual polynomial of the passport,
which lifts to a point in this space. In the boundary case this
also proves that the open set \(L_u\ne0\) is nonempty. No finite
list of chart witnesses is being used to prove this assertion.
\(\square\)

## Theorem K2. One component and unirationality in the stable range

Every admissible passport (1) satisfying \(s\ge1\) and (2) has
one absolute and one literal full-braid orbit. Its polynomial
Hurwitz space, modulo source and target affine changes, is
unirational over \(\mathbb Q\), of complex dimension \(k+s-2\).

**Proof.** In (6), the ambient space is the product of the ordered
distinct-point configuration, the \(k\) value coordinates, and
the \(d\) free coefficients of the monic polynomial \(R\).
It is irreducible and rational, of dimension
\[
                   p+k+d=k+s+1.                              \tag{10}
\]
In (7), the base is the nonempty open set of configurations on
which \(L_u\ne0\). Above it \(L_u(v)=1\) is an affine hyperplane
bundle of rank \(k-1\). This base is irreducible, and a locally
trivial affine bundle over it is irreducible. Solving for any
value with a nonzero coefficient in \(L_u\) gives a rational
dense chart. Its dimension is \(p+k-1=k+s+1\).
The good locus (9) is a nonempty open subset, hence irreducible
and connected in either of these cases. For \(h\ge1\), condition
(2) is \(k\ge h+2\); K5--K6 below prove that the locus (M2)
is smooth and irreducible on (9), with a rational dense chart
over the full-rank configuration open. This also gives (10)
and the same conclusion in every remaining case.

The finite quotients forgetting equal-index point labels and
equal-profile branch labels have the entire polynomial passport
as image. Source and target affine quotients are also images
of this connected irreducible locus. Their dimension drops by
three from the monic coefficient space: in the four-dimensional
source/target affine group, preservation of monicity imposes
one relation. Stabilizers are finite for a polynomial with at
least two distinct finite branch values. The quotient dimension
is therefore \(k+s-2\). A dominant rational parametrization
survives all these quotients, proving unirationality; rationality
of these quotient fields is not asserted.

Over the configuration space of distinct branch values the
polynomial Hurwitz map is a covering (with the usual finite
automorphism convention). A path between two covers projects
to a braid, so connectedness is exactly transitivity on the
absolute Nielsen fibre, with every factor-type placement allowed.
For a fixed product, the remaining relabelings are powers of
\(P_N\). The full twist acts by conjugation by that product.
Thus each absolute orbit lifts to a single literal orbit as well.
\(\square\)

For \(k=1\), admissibility already requires \(M\le N\):
its sole branch has \(M\) prescribed preimages counted with
multiplicity. Every case using the moment equations therefore
involves multiple non-transposition branch values.

## Theorem K3. Actual groups, all blocks, deck groups and sizes

In the range of K2, if \(k\ge2\), every tuple generates \(S_N\).
Its only block-count divisors are \(1,N\), and its deck group
is trivial.

If \(k=1\), let \(D\) be the gcd of all parts of its full
partition, including fixed-point parts, and set \(n=N/D\).
Then \(n\ge2\) and the actual embedded group is
\[
                              S_n\wr C_D.                     \tag{11}
\]
Its block-count divisors are exactly
\(\{d:d\mid D\}\cup\{N\}\). The deck order \(e\) is two
when \(n=2\), and one otherwise.

For either case let \(\nu_\lambda\) count the occurrences of each
complete factor partition among all \(r=k+s\) factors and put
\[
 W=\frac{r!}{\prod_\lambda\nu_\lambda!},\qquad
 F=W\,N^{r-1}\prod_{i=1}^r
                         \frac{(c_i-1)!}{\prod_a m_{i,a}!}.     \tag{12}
\]
The unique full component has \(F\) literal tuples and \(Fe/N\)
absolute tuples, where \(e=1\) for \(k\ge2\). These are full
type-placement counts, not fixed ordered-pattern counts.

**Proof.** Use L6 of the
[return-transport note](POLYNOMIAL_PASSPORT_LIFT_TRANSPORT.md#lemma-l6-monodromy-containing-a-transposition).
The normal transposition subgroup has blocks of size \(b\ge2\),
say \(N=bm\), and the group is the full symmetric base
\(S_b^m\) over its quotient. If \(m=1\), it is \(S_N\).
Suppose \(m>1\). Every literal transposition fixes the quotient
and contributes one to the total return-index budget \(b-1\).

Let \(c\ge1\) count all nontrivial quotient cycles of the
non-transposition factors. A quotient cycle of length \(\ell>1\)
with return \(h\in S_b\) lifts to \(b-\operatorname{ind}(h)\)
nontrivial original cycles: all their lengths are \(\ell\) times
return-cycle lengths.

Let \(a\le c\) be the number of factors moving quotient blocks;
the other \(z=k-a\) factors have total index \(I_o\) and
nontrivial-cycle count \(p_o\ge z\). Their entire indices
come from returns. If \(R\) is the sum of return indices
over nontrivial quotient cycles, H3 gives
\[
 R\le b-1-s-I_o,\qquad
 p-k\ge cb-R-a+p_o-z
       \ge(c-1)(b-1)+s+I_o.                                 \tag{13}
\]
If \(c\ge2\), this is at least \(s+1\), contrary to \(p-k\le s\).
Thus \(c=1\). The quotient index is \(m-1\), so its unique
nontrivial cycle is an \(m\)-cycle, in one original factor.
Write \(t\) for that slot's return index. This original factor
has exactly \(b-t\) nontrivial cycles and no fixed points.
All other non-transposition factors fix quotient blocks; write
\(I_o,p_o\) as above.
Then
\[
 b-1=s+t+I_o,\qquad
 p-k=s+I_o+p_o-(k-1).                                       \tag{14}
\]
If \(k\ge2\), an additional factor has index at least two and
at least one nontrivial cycle, so \(p-k\ge s+2\), another
contradiction. This proves the full symmetric assertion.

For \(k=1\), any proper quotient has just the one full quotient
cycle, so its block count divides every part of the full
partition, and therefore divides \(D\). Conversely, translating
its branch value to zero writes
\[
             f(z)=\prod_x(z-x)^{\lambda_x}=H(z)^D.             \tag{15}
\]
This gives a preserved cyclic \(D\)-block quotient, including
the trivial quotient if \(D=1\). The normal transposition
blocks refine every block system with nonsingleton blocks, because a transposition
cannot move whole blocks of size greater than one. Their block
count must therefore be at least \(D\); as it divides \(D\),
it is exactly \(D\). L6 gives (11).

There are at least two sheets in each transposition block.
The standard fixed-product labels identify all blocks with
residue classes. Every nonsingleton block system coarsens
the normal transposition partition, and cyclic \(C_D\) has
exactly the partitions indexed by divisors of \(D\). L6 gives
the deck assertions, also for \(D=1\).
Finally (3) counts each ordered pattern and \(W\) counts its
distinct placements. K2 puts all \(F\) tuples in one orbit.
Simultaneous conjugation preserving \(P_N\) has stabilizer the
deck group and orbit size \(N/e\), proving the absolute count.
\(\square\)

## Proposition K4. The first excluded level already splits

On six sheets the passport
\[
                    (2^2 1^2),\quad(2^2 1^2),\quad T          \tag{16}
\]
has \(p=4,k=2,s=1=p-k-1\). Its complete literal domain has 243 tuples.
It has two components: 216 with group \(S_6\), deck one and
36 absolute tuples; and 27 with group \(S_2\wr D_3\), deck
two and nine absolute tuples. Here \(D_3\) has order six, so
the second group has order 48.

**Proof.** This is a case of B1--B5 in the
[two-involution classification](POLYNOMIAL_INVOLUTION_COMPONENT_CLASSIFICATION.md).
The checker below independently takes the two signed braid
closures, checks their actual groups, blocks and deck orders,
and checks that their disjoint union exhausts the ordinary
hook-character count 243. This finite verification alone is
also a proof of this particular counterexample. \(\square\)

## Lemma K5. Moment derivatives and the polynomial product bound

Assume \(h\ge1\) and use the labeled critical points of K1;
write them as \(u_1,\ldots,u_p\), with local indices \(q_j\).
Put \(P(z)=\prod_j(z-u_j)\). For any value vector \(w\) with
\(F(u)w\in\mathbb C e\), let \(g=\Phi_w\), so that
\(\deg g\le N\), and put
\[
                           B_g=g'/Z\in\mathbb C[z]_{\le s}.
\]
For a polynomial \(\beta\) of degree at most \(h\), interpreted
as a row vector of moment coefficients, one has
\[
 \beta\, (dF)_u(\dot u)\,w
             =-\sum_{j=1}^p
                \frac{\beta(u_j)B_g(u_j)}{P'(u_j)}\,\dot u_j.  \tag{J}
\]
In particular the Jacobian of \(F(u)v=e\) has row rank \(h+1\)
at every admissible point of (M2).

If \(U,V\) are nonzero finite-dimensional subspaces of
\(\mathbb C[z]\), then
\[
                  \dim\operatorname{span}(UV)
                                \ge\dim U+\dim V-1.           \tag{P}
\]

**Proof.** Vary one point \(u_j\), keeping \(w\) fixed. The
variation \(\dot g\) has zero jets of orders below \(q_j-1\)
at every prescribed point, and its jet of order \(q_j-1\)
at \(u_j\) is \(-\dot u_j g^{(q_j)}(u_j)\). Writing
\(Q_j=Q/(z-u_j)^{q_j}\), the unique degree-less-than-\(M\)
polynomial with these jets is
\[
 \dot g=-\dot u_j
       \frac{g^{(q_j)}(u_j)}{(q_j-1)!}
       \frac{(z-u_j)^{q_j-1}Q_j(z)}{Q_j(u_j)}.               \tag{D}
\]
Since \(g'=ZB_g\), its scalar coefficient after dividing
by \(Q\) becomes \(-\dot u_j B_g(u_j)/P'(u_j)\).
Also
\(\dot Q=-q_j\dot u_j Q/(z-u_j)\).
In differentiating the coefficient
\([z^{-1}]\beta g/Q\), the term from \(\dot Q\) contributes
zero: \(\deg\beta\le h,\deg g\le M-h-1\), so that
\(\beta g/[Q(z-u_j)]\) is \(O(z^{-2})\) at infinity.
The term from (D) contributes \(\beta(u_j)\) times its
scalar, proving (J).

For an admissible monic \(f\), \(B_f=NA\) and \(A(u_j)\ne0\).
The \(u\)-Jacobian in (J) for
\(\beta=1,z,\ldots,z^h\) is a Vandermonde matrix with nonzero
column multipliers. Since \(p\ge h+1\), it has rank \(h+1\).

For (P), choose bases with distinct increasing leading
degrees \(a_1<\cdots<a_d\) and \(b_1<\cdots<b_e\).
The products with degrees
\(a_1+b_1,\ldots,a_d+b_1,a_d+b_2,\ldots,a_d+b_e\)
have distinct leading degrees, hence are linearly independent.
\(\square\)

## Theorem K6. Rank loss cannot create another component

Suppose \(h\ge1\) and \(k\ge h+2\). The admissible locus (M2)
is irreducible, smooth of dimension \(p+k-h-1=k+s+1\),
and rational over \(\mathbb Q\). Its full-rank configuration
open is nonempty. Every feasible incidence locus of rank
\(r<h+1\), with deficiency \(\delta=h+1-r\), has dimension
at most \(p-\delta+1\), strictly below the dimension of
the admissible locus.

**Proof.** Let \(\Sigma_r\) be the locally closed set of
configurations with \(\operatorname{rank}F=r\) and
\(e\in\operatorname{im}F\). At a general smooth point of
any of its irreducible components the kernels below are
vector bundles and can be differentiated within that component.
Set
\[
 U=\ker(F^t)\subset\mathbb C[z]_{\le h-1},\qquad
 S=F^{-1}(\mathbb C e),\qquad
 V=\{\,\Phi_w'/Z:w\in S\,\}\subset\mathbb C[z]_{\le s}.
                                                               \tag{S}
\]
Here \(\dim U=\delta\). The smaller degree bound in (S)
follows from \(\beta e=0\) for \(\beta\in\ker(F^t)\).
Also \(\dim S=k-r+1\). The derivative map in (S) has
one-dimensional kernel, precisely the constant polynomials:
the CRT values of a constant are all equal. Thus
\(\dim V=k-r\).

For tangent vectors to \(\Sigma_r\), differentiating
\(\beta F=0\) and \(Fw=\alpha e\), with \(\beta e=0\),
gives \(\beta(dF)w=0\). Formula (J) therefore puts the
evaluation vectors of every product \(\beta B_g\) in the
conormal space. These products have degree at most
\[
                         (h-1)+s=p-3.
\]
Evaluation at the \(p\) distinct points, followed by division
by the nonzero \(P'(u_j)\), is injective on this polynomial
space. By (P) the conormal rank is at least
\[
                           \delta+k-r-1.                     \tag{C}
\]
This proves
\(\dim\Sigma_r\le p-(\delta+k-r-1)\).
The affine solution fibre \(Fv=e\) has dimension \(k-r\);
the feasible rank-\(r\) incidence locus consequently has
dimension at most \(p-\delta+1\).

K5 shows directly that the admissible incidence locus is
smooth of pure dimension
\[
                         D=p+k-h-1.
\]
Since \(k\ge h+2\), one has \(D\ge p+1>p-\delta+1\).
Thus no irreducible component can be supported on a deficient
rank stratum. Nonemptiness from K1 now implies the existence
of a full-rank admissible point.

Over the nonempty configuration open where \(F\) has rank
\(h+1\), the solutions \(Fv=e\) form an affine bundle of rank
\(k-h-1\) over an irreducible rational base. Its good open
is irreducible. Every component of the smooth admissible
locus must meet this one open, so there is only one component.
Solving any nonzero maximal minor gives a rational dense
chart, proving rationality of this root-marked coefficient
space. The unmarked affine quotient is only asserted
unirational, as in K2. \(\square\)

## Exact support, audit and remaining scope

Run from the canonical checkout:

    python3 research/arithmetic/verify_polynomial_stable_hermite.py

The checker constructs exact rational CRT and moment charts,
including \(h=1,2,3,4\), verifies critical-value norms by two exact methods,
rejects the relevant collisions and invalid profiles, and
compares selected complete signed braid orbits with an independent
ordinary hook-character total. It checks normal transposition
supports, quotient orders, all block divisors and deck orders.
It independently differentiates moving-contour residues to
test (J), checks the polynomial product/conormal dimensions
at a feasible rank-deficient configuration, and checks the
one-dimensional conormal in the excluded degree-six case.
Three additional complete controls use pure braid generators
and same-type exchanges on a fixed ordered pattern, quotienting
by product rotations. These generate the colored braid group:
the pure generators generate its permutation kernel, and
same-type exchanges generate its finite permutation image.
The three orbit sizes exhaust the independent fixed-pattern
absolute totals, giving full domains of 15,360, 700,000 and 1,771,470
literal tuples after type placements.

| Degree | Finite passport | \(p-k=s\) | Fixed-pattern absolute tuples | Full literal tuples |
|---|---|---:|---:|---:|
| 8 | \(\mathrm{DT},3,3,T\) | 1 | 160 | 15,360 |
| 10 | \(\mathrm{DT},3,3,3,T\) | 1 | 3,500 | 700,000 |
| 9 | \(\mathrm{DT},\mathrm{DT},3,T,T\) | 2 | 6,561 | 1,771,470 |

Each row is one \(S_N\) component. Here \(3\) denotes one
three-cycle and all other sheets fixed. These three examples
lie below the former \(s\ge p-2\) sufficient range.
The degree-nine row has two disconnected non-transposition
factors, so it also lies outside the one-arbitrary-factor
and single-cycle-tail family.

The bounded domains and all inputs are retained in
[the receipt](polynomial_stable_hermite_verification.json);
[the source ledger](polynomial_stable_hermite_source_ledger.json)
separates the written proof and imported existence/count input
from these finite controls. The all-parameter statements use
K1--K6, not an extrapolation of the finite orbits.

The general component goal remains open below \(s=p-k\),
and at \(s=0\), outside the subfamilies already classified.
Power phases, intrinsic branch classes and genuine reduced
lifting splits already occur outside this range. Neither
unrestricted connectivity nor a general period-cell gluing
algorithm follows from this theorem. Independent mathematical
review and replay remain pending; no external review or formal
verification is claimed.
