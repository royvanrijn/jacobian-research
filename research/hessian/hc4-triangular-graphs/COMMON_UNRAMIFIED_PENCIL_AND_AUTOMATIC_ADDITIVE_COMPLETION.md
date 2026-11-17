# The common unramified pencil makes the etale completion additive

HC4MYMP0. 16 November 2026. ALL actual complex polynomial plane Keller
maps, in every generic degree, with all omitted primes retained.
This combines the previously proved PF2RULED46/50 with HC4MYLM0/1.
It closes their previously uncombined rigid-completion/noncyclic-deck
alternative. Full HC4/JC2 and the full trace/marked-energy gate remain
OPEN. No novelty, independent replay, formal verification or external
review is asserted; all predecessor proofs and checkers stay sealed.

## 1. Precise retained inputs

Let K:X=A2->Y=A2 be an actual complex polynomial Keller map. Write
R=C[X], B=C[K], L=Frac(R), F=Frac(B), and N for the FULL finite
normalization of Y in L. Let U be its ENTIRE etale locus. KU/LA,
as retained in HC4MYLA0 and HC4MYLM0/1, give

    U=Spec(A), A=B[Lambda] subset R,
    U minus X = D_1 disjoint-union ... disjoint-union D_b0,
    D_i=A1, U smooth affine, pi1(U)=1, chi(U)=1+b0.

The D_i are precisely the omitted UNRAMIFIED primes. Ramified omitted
primes lie outside U; they are neither erased from N nor counted in b0.
The source volume equals a nonzero constant times dP wedge dQ, so it
extends to an exact nowhere-zero form on U.

The crucial retained input is **PF2RULED50**, Sections2--4:

    ALL unramified omitted primes of the ACTUAL N have the SAME
    intrinsic regular SOURCE-coordinate pencil.                  (1)

This has no field-degree or component-count bound. Its proof uses
PF46's full component charts, the arbitrary-word source-coordinate
lemma, Xie's adopted intersection-matrix theorem and the exact-volume
bound n_i<=p_i^2. For distinct pencils their cross pairing m*p_i*p_j
would violate strict negative definiteness. We reuse that written
all-degree theorem and its stated classical-input limits; finite
checker samples are not its proof. No new Xie/bundle/automorphism
classification is adopted here.

For each i, M_i=X union D_i is a smooth AFFINE component open of N,
by the component-open theorem already used in PF46/50. PF2RULED46
Sections6--7 gives its WHOLE normalized root charts. Align the pencil
in (1) by choosing source coordinates (x,y), so that on M_i a
normalized coordinate pair is

    y_i=b_i*(y-lambda_i), x_i=a_i*x+p_i(y), a_i*b_i!=0.

The complementary-coordinate change has this form because it is an
automorphism of the polynomial plane over the same coordinate y.
On the degree-r_i root chart of the component open,

    y_i=z^r_i, x_i=h_i(z)+z^(1-r_i)*t_i, r_i>=1,       (2)

with h_i the retained finite negative Laurent principal part. The
root cover is etale near D_i and covers the ENTIRE boundary A1 by
plane charts. It is branched along the retained source line y_i=0;
we check the derivation separately on X there. Exactness removes the
coefficient of z^(-r_i), as in PF46/50, but imposes no new condition
in the extension argument below. The lambda_i may be DIFFERENT.

## 2. One source translation extends across ALL unramified primes

Assume b0>0 and choose the coordinates (x,y) above. Consider the
SOURCE locally nilpotent derivation

    delta=partial_x on R=C[x,y].

In (2), delta(z)=0: differentiating z^r_i=y_i in the function field
gives this in characteristic zero. Since delta(x_i)=a_i, we obtain
the EXACT boundary expression

    delta=a_i*z^(r_i-1)*partial_(t_i).                         (3)

The exponent r_i-1 is NONNEGATIVE for every boundary order. Thus
(3) is a regular algebraic vector field on each whole root plane,
not just a formal/generic boundary germ. It is invariant under the
root deck action because it is the same field derivation partial_x;
etale descent near D_i proves regularity at every point of D_i.
On X it was regular from the start. This also avoids any invalid
descent of vector fields through the branched source line.

The component opens M_i cover U: every point is in X or one of the
D_i. Their vector fields agree on the common function field, hence
glue to a regular vector field on U. Therefore

    delta(A) subset A.                                           (4)

Local nilpotence follows from the ACTUAL inclusion A subset R:
for each a in A, its polynomial x-degree is finite, so some power
of partial_x annihilates it. No integration of an arbitrary complete
field, infinite Lie bracket algebra or merely formal flow is used.
The finite exponential on every element defines an algebraic Ga
action on U, and its restriction to the WHOLE source X is exactly

    (x,y) -> (x+s,y).

This action PRESERVES X and every D_i. It is nontrivial because
Frac(A)=C(x,y), so delta cannot kill all A. Also y is regular on
every M_i, hence belongs to A, and

    ker(delta|A)=A intersect C[y]=C[y].                          (5)

For b0=0 we have U=X and source translations already suffice.
Consequently EVERY actual canonical etale completion U admits a
nonzero additive action. The proof uses the actual common-pencil
theorem (1), not exact volume alone or a conjectural automatic ruling.

## 3. Proper completions: all additive actions preserve the source

For b0>0, HC4MYLM0 proves that a proper actual U with an additive
action has a UNIQUE LND kernel; two distinct kernels would force
U=A2 by its retained classification/period argument. Combining with
(5) yields

    ML(A)=C[y], not A.                                           (6)

The former rigid alternative ML(A)=A is therefore EMPTY for actual
canonical completions under the retained inputs. It remains a valid
abstract alternative for surfaces not satisfying those inputs.

There is a stronger explicit description. Let theta be ANY nonzero
LND on A. Its kernel is C[y], by LM0. Put

    H(y)=product_(distinct lambda_i)(y-lambda_i).

All omitted primes lie over its zero set, so U_H=X_H and
A[H^-1]=C[x,y,H^-1]. Localizing further at C[y] minus zero gives

    A tensor_(C[y]) C(y)=C(y)[x].

Since theta kills y, localization preserves local nilpotence. A
nonzero LND of the one-variable polynomial ring C(y)[x] is
q(y)*partial_x. Here is the elementary degree reason: if theta(x)
has x-degree k>=1, successive theta^n(x) have nonzero leading term
of degree 1+n*(k-1); in characteristic zero the leading coefficient
recurrence multiplies by positive integers, so no power vanishes.
Thus k=0.

Restricting the regular vector field theta from U to its open X
gives regular coefficients on X, even though an arbitrary restricted
flow need not preserve an open. In the present expression its x
coefficient is q(y), so

    q(y) in C(y) intersect C[x,y]=C[y].

Conversely every polynomial q(y)*delta preserves A by (4),(5) and
is locally nilpotent. Including the zero derivation, we have exactly

    LND(A)= {q(y)*partial_x : q in C[y]}   when b0>0.              (7)

Each such field has the polynomial SOURCE flow
(x,y)->(x+s*q(y),y), so EVERY additive action on a proper canonical
completion preserves the marked source X. This is stronger than
assuming that an arbitrary action restricts to a complete field:
we proved its explicit form first. It does not say that every
automorphism or generic deck transformation preserves X.

## 4. All generic rational deck groups are cyclic, without an index bound

Let D=Aut_F(L), the FULL generic rational deck group. HC4MYLA0 makes
its action on U regular and FREE; D need not act polynomially on X.
HC4MYLM1 proves that the existence of a nonzero additive action on
this actual U forces D CYCLIC. Section2 supplies that premise in
EVERY degree. Hence

    Aut_F(L) is CYCLIC for EVERY actual complex plane Keller map. (8)

For clarity, the retained LM1 argument is short. In the proper case
D preserves the intrinsic algebra (6), so acts affinely on A1_y.
An element in the kernel acts on the whole generic affine-line fibre
by x->a(y)*x+b(y). A nontrivial finite-order affine transformation
in characteristic zero has a fixed point on that generic A1.
The algebraic fixed locus then has a closed point on U, contradicting
freeness. The action on the base is faithful. A finite subgroup of
Aut(A1) fixes a point by averaging and embeds into C*, hence is cyclic.
If b0=0, LA0 already gives |D| dividing chi(U)=1, so D is trivial.

Thus EVERY noncyclic subgroup of the permutation centralizer of
actual transitive plane monodromy is excluded, in ANY action and
ANY degree. In particular any relative GALOIS extension L/E at the
top of a genuine intermediate tower has CYCLIC group, because
Aut_E(L) is a subgroup of Aut_F(L). No relative-etaleness or Euler
bound is needed for this statement.

This supersedes the noncyclic-deck survivor regions in LA1, LM1,
PF132--134 and MN0: there is no remaining noncyclic deck case at
index five, at binary index eight, or at any larger index. Their
sealed conditional inequalities and exclusions stay true; the new
combination closes their antecedent globally. PF134's dihedral and
platonic DECK alternatives are eliminated. Its noncyclic fundamental
groups of arbitrary genuine intermediate etale surfaces are NOT
thereby eliminated: their universal covers need not be the actual U,
and their automorphisms need not lift to L. General non-Galois relative
monodromy and the full original monodromy remain distinct from D.

The cyclic order still obeys |D| dividing gcd([L:F],1+b0), and may
be nontrivial. Known exact-volume W_k controls have free cyclic
symmetry that moves their chosen plane chart, consistent with (8).
They are reused, not reconstructed or replayed, and still lack the
actual etale target pair W_k->A2.

## 5. What remains, and the source audit

This is a new combination of previously retained ALL-DEGREE inputs.
The proof reads PF46 Sections6--7, PF50's complete common-pencil and
matrix proof, and the complete LM0/LM1/LA0 proofs. It adopts no new
external theorem and does not rerun any registered mathematical checker.
In particular the classical whole-chart, Xie and ML-classification
foundations keep their original adopted/not-independently-reproved
assurance boundaries. A NEW preregistered checker only supports the
coordinate/gluing arithmetic, LND degree argument, finite affine
fixed-point transcription, provenance and scope.

The direct search for new coframe/divisibility classifications supplied
no new usable theorem. The primary 2026 Matysiak SSRN record7229458
asserts a plane proof, but only its abstract was accessible; the paper
and its hypotheses/proofs were NOT read or adopted. The recent
Gurjar--Thandar finite-fundamental-group-at-infinity record was screened
at its primary introduction only; finite pi1 at infinity was NOT
established for U. Neither screen is evidence of HC4 closure. The
new result instead uses the complete canonical proofs already retained.

Even U=X (b0=0) can leave omitted RAMIFIED primes in N, so (4)--(8)
do not prove the Keller inverse. A nontrivial cyclic D also remains
unexcluded. The action constructed here need NOT preserve B or F,
need NOT extend over ramified primes, and need NOT descend to a target
Ga action. No second independent LND, target-coordinate ruling,
full trace-dual multiplication, marked-energy integrality or automatic
solvable-monodromy inversion follows.

The remaining global condition concerns those actual target/finite
normalization conditions. Reopening noncyclic deck families or a rigid
canonical completion would repeat a now-closed branch. The published
quartic PLANE exclusion is retained; no quartic intermediate is assumed
plane. Full HC4/JC2 remains OPEN; independent external review is unclaimed.
