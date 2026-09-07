# Branch divisibility bounds whole native twist blocks

Every one of the **37 retained native quadratic branch supports** has
arithmetic generic twist rank at most **three**. This holds for every
nonzero rational scalar multiple of its branch polynomial, including
scalars not used in the original construction. It bounds all rational
sections of each twist, not just the displayed native section.

For any product of at least two distinct supports in this set, the twist
rank is at most **two**. Consequently a fibre product of two of these
quadratic covers can add at most **eight** generic directions in total.
A proposed generic +9, +10 or +14 block carried by this full native
fibre-product construction therefore needs at least three supports;
their fibre-product normalization has genus at least five.

These are **incidence capacity bounds**. They do not exclude an individual
fibre from gaining further rank, and do not show that the upper bounds
are attained. The new obstruction is divisibility at the *branch fibres
of the proposed parameter cover*, not a point-search statistic at the
successful specialization.

## A stronger arithmetic gate than geometric room

Let F=Q(t) and let E/F be the published R17 family. For a nonsquare
d in F*, identify E^(d)[2] with E[2]. Define S_d to be the subspace of
H1(F,E[2]) satisfying the local Kummer condition for E^(d) at every closed
place of the parameter line. This definition uses completions k(b)((u))
with number-field residue fields; it is not the ordinary number-field
Selmer group of a specialized curve. Rational twist points inject into
S_d modulo two.

Let L be the original global pool from
[the cubic-root-curve theorem](ROOT_CURVE_TORSION_AND_REAL_CAPACITY.md),
and let G be the image of the 17 marked generic sections. Suppose all
branch places b of d are smooth fibres and

\[
 E_b(k(b))[2]=0.
\tag{1}
\]

Equivalently the branch-fibre cubic is irreducible over k(b). Then

\[
 S_d\subset L,\qquad
 S_d\subset\ker\left(L\longrightarrow
       \prod_{b\mid\operatorname{branch}(d)}H^1(k(b),E_b[2])\right).
\tag{2}
\]

If a finite character calculation on the specialized generic sections
has rank r, the specialization map in (2) has rank at least r. Hence

\[
 \boxed{\operatorname{rank}E^{(d)}(F)
 \le\dim S_d\le\dim L-r.}
\tag{3}
\]

No completeness claim for the branch-fibre Mordell-Weil group or its
Selmer group is required. Even global independence of the supplied
generic classes is not needed for (3): a rank-r image already supplies
the required r-dimensional subspace. The independent generic calibration
below is used to describe the branch kernels in the marked 17-coordinate
basis.

## Why the local condition vanishes at an irreducible branch fibre

Here is a direct proof of the new local step. Complete at b, writing
d=c*u with c a unit; square factors have been removed. The original
smooth equation has coefficients A,B in k(b)[[u]], with irreducible
residue cubic f_b and therefore nonzero residue B. The twisted equation is

\[
 y^2=x^3+d^2Ax+d^3B.
\]

If v(x)=1, cancellation of the order-three leading term would require
a root of f_b in k(b); without that cancellation the valuation is odd,
which is impossible for y^2. If v(x)>1, the term d^3B alone has lowest
valuation three, also impossible. The case x=0 is excluded the same way.

Thus every finite rational point has v(x)=2a<=0. Put x=u^(2a)X and
y=u^(3a)Y. Reduction gives Ybar^2=Xbar^3, so Xbar is a nonzero square.
In the unramified cubic algebra defined by f, the Kummer representative is

\[
 x-d\theta=u^{2a}\bigl(X-c u^{1-2a}\theta\bigr).
\]

Its unit factor has square residue and is a square by Hensel lifting.
The point at infinity has trivial class too. The whole local twist
Kummer image is therefore zero. This is also the I0* component-group
description: its nonidentity components are indexed by the three
two-torsion roots, none of which is rational over k(b).

At an unbranched place, d is a square after passing to the geometric
completion, so the twist has the original geometric local Kummer
condition. At the branch places the stronger local vanishing just proved
holds. This proves S_d subset L and its zero specialization in (2).
There is no additional branch at infinity for the even-degree products
used here.

The norm interpretation explains the divisibility language. If a class
of P in G is also represented by a rational twist point Q, restriction
to F(sqrt(d)) gives P-iota(Q)=2R. There is no rational two-torsion over
this quadratic extension, since the generic cubic remains irreducible.
Thus R+sigma(R)=P. At a branch place sigma acts trivially on the residue
field, and smooth proper specialization gives P(b)=2R(b). In particular,
shared rational classes must lie in the branch-divisibility kernel.
This is consistent with the standard Kummer-intersection/norm identity;
see Mazur--Rubin, [Lemma 2.9](https://arxiv.org/pdf/0904.3709). The argument
above proves the needed function-field case directly.

## The bounded calculation on the 37 retained supports

The [protocol](BRANCH_DIVISIBILITY_CAPACITY_PROTOCOL.json) uses precisely
the 37 quadratics in the existing five-fibre native census. Each is
matched coefficient by coefficient to the earlier complete equation-only
atlas. Historical selection by rational lifts is retrospective; it is
not a prospective selector. The eleven oracle-fitted quartics are
excluded. Arithmetic receives only A,B, the generic section polynomials,
branch polynomials and labels. No exceptional coordinates or quotient
classes are supplied, and no point or section is searched for.

For primes through 2003, every simple rational reduction of a quadratic
branch root is tested. A reduction where the fibre cubic is irreducible
certifies (1) over the quadratic branch field: a characteristic-zero
root would reduce to a root at this good degree-one prime. Split cubic
reductions provide quadratic characters of the 17 generic sections.
Blocks where a section ordinate vanishes are skipped. All remaining
blocks in the fixed dictionary are retained; collection does not stop
when a desired rank is reached.

The [completed arithmetic](../../artifacts/generated-results/elliptic-curves/rank_jump_branch_divisibility_capacity_v2.json)
and [independent verification](../../artifacts/generated-results/elliptic-curves/rank_jump_branch_divisibility_capacity_verification_v1.json)
give:

| Quantity | Result |
|---|---:|
| Native quadratic supports | 37 |
| Branch cubics certified irreducible over their quadratic fields | 37 |
| Branch character rank, on every support | 16 |
| Generic branch-divisibility kernel dimension | at most 1 |
| Distinct one-dimensional finite fingerprint kernels | 37 |
| Pairwise disjoint branch-support pairs | 666 |
| Real components of the original cubic root curve | 10 |
| Arithmetic global-pool capacity | 19 |
| Generic rank of every scalar twist on a single support | at most 3 |

The genus-ten root curve has 18 real simple branch values. Rational
Sturm isolation and ordered-root gluing give ten real components, hence
dim J[2](R)=10+10-1=19. This computes the real bound directly on the
published model; it does not assume that a different fibration chart is
the same root curve.

The independent verifier checks all **1604** finite blocks, comprising
**81804** quadratic-character bits. It uses exhaustive finite cubic root
tests and residue-square sets instead of the worker's polynomial
factorization and Euler tests, and Sage F2 ranks instead of its bit
elimination. It separately checks the Sturm intervals and counts graph
components. Forty-seven blocks certify generic dimension 17; the other
1557 are branch blocks.

The v1 UNKNOWN artifact is preserved: the finite fingerprint at the
initial t=0 calibration did not certify generic dimension 17. The
[completion protocol](BRANCH_DIVISIBILITY_COMPLETION_PROTOCOL.json)
uses the already retained parameter -2/377 for that generic-only
calibration. Branch supports, branch arithmetic and prime cap are
unchanged. Failure of the first fingerprint is not an exact rank claim
about its fibre.

## Products require no enumeration of exponentially many covers

The 37 finite kernel lines are pairwise distinct. For any two supports,
their joint generic character map has rank 17. Every product involving
at least two distinct supports contains such a pair in its branch
divisor, so (3) gives rank at most 19-17=2. Pairwise coprimality ensures
that no branch place cancels. This proves the statement for all nonempty
subsets and all nonzero rational scalar multipliers.

For a chosen set of k supports, the generic Mordell-Weil group over the
full multiquadratic cover decomposes over Q into its quadratic character
spaces. Integral gluing can change indices, but cannot change these
rational ranks. There are k singleton characters and 2^k-1-k larger
characters. The [capacity certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_branch_character_product_capacity_v1.json)
therefore gives

\[
 \boxed{\operatorname{rank}E(\mathbb Q(C_k))
 \le17+3k+2(2^k-1-k)=15+k+2^{k+1}.}
\]

| Supports k | Cover genus | New generic rank, at most | Full generic rank, at most |
|---:|---:|---:|---:|
| 1 | 0 | 3 | 20 |
| 2 | 1 | 8 | 25 |
| 3 | 5 | 17 | 34 |

The genus formula is 1+2^(k-1)(k-2), since the k quadratics have disjoint
branch divisors. The lower bound genus five for a proposed generic gain
of at least nine applies to this **full native fibre-product
construction**. It is not a bound on every possible auxiliary curve or
on an unrelated quotient of such a curve.

Thus the retained published-R17 +9, +10 and +11 fibres cannot have their
whole gain explained by the complete generic section group on one or
two of these covers, even if the needed rational lift exists. A fibre
with gain J would still require at least J-3 further directions after
one cover, or max(0,J-8) after two. At the observed +4 control, even a
single cover leaves at least one further direction. None of these
statements bounds a fibre's actual rank or assumes injectivity of the
entire specialized generic group.

## What class-creation structure remains possible

There is a concrete arithmetic distinction at a potential branch place:
irreducible cubic, no new ramification; a cubic with rational two-torsion,
possible new ramification. In general the arithmetic residue of a
global E[2] class at a smooth closed place b lies in E_b(k(b))[2]. Its
dimension is 0, 1 or 2 according as the cubic is irreducible, splits as
linear times irreducible quadratic, or splits completely. This is a
capacity for new classes, not an existence theorem. It concerns closed
places of the parameter line, not prime divisors of the specialized
number-field discriminant.

For the tested supports that residue capacity is zero. A large rational
twist block would instead have to fit within L and lose many generic
directions in the branch-specialization map. Specifically, a block of
n directions requires branch character rank at most 19-n. A proposed
ten-direction block would require rank at most nine, whereas all 37
measured ranks are sixteen. This is a falsifiable construction gate,
not a correlation between a score and fibre rank.

The strongest remaining mechanisms are therefore: branch fibres with
new rational two-torsion over their residue fields; substantially larger
simultaneous generic divisibility at branch fibres; or carriers outside
this quadratic-character construction. Each still needs independent
class and rational-solubility certificates. A large twist Hodge bound,
a square branch value at t0, or a single visible native section does
not supply those certificates.

The condition on t0 that creates the unramified class block in the
fresh +10/+11 and historic +12/+14 panel is still missing. This is not a
replacement measurement of that panel's additional classes or CT form.
The result sharply limits a candidate carrier mechanism before such
points are supplied. It introduces no visibility feature or candidate
score, and changes no active search work.

```sh
timeout 60 sage -python elliptic-curves/rank-jump/verify_branch_divisibility_capacity.py check
python3 elliptic-curves/rank-jump/branch_character_product_capacity.py check
```
