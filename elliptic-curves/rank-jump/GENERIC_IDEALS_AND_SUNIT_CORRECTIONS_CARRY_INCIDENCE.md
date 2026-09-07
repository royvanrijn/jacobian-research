# Generic half ideals and S-unit corrections carry a strict incidence block

The full generic half-ideal image detects **all** known strict characters
on the three earlier controls: dimensions **10, 8 and 6**. On MW16 high
this improves the previously certified elementary S-class factor from
nine to **ten**, using half ideals of generic points as its generators.

The S-unit-only exclusion does not extend to arbitrary generic
corrections. More positively, there is an exact sufficient incidence
criterion:

\[
 \boxed{\dim\Phi_S(G)\ge d
 \ \Longrightarrow\
 \dim\bigl((G+\mathcal E_S)\cap U\bigr)\ge d.}
\tag{1}
\]

Here G is the marked generic Kummer subgroup on a rational fibre, U is
the strict Selmer space, Phi_S takes localized half ideals in the cubic
field, and E_S denotes **norm-square S-unit squareclasses**. The latter
classes need not be locally soluble or strict. Equation (1) says that
enough generic half ideals remaining independent after localization
force enough combinations of generic classes and S-units to become
strict together.

On the retained controls this forces at least **9, 6 and 0** strict
directions beyond the generic strict subgroup in this particular
constructor. These are incidence conclusions. The new combinations
have not been written explicitly, and their rationality versus Sha is
not decided. The known rational strict block need not equal this
constructed-in-principle block.

The independence certificates below use **old exceptional-point-derived
characters retrospectively**. The generic ideals themselves require only
the equation and generic points. There is no point-independent
independence certificate or new selector for the fresh panel yet.

## What the larger arithmetic test measures

For K=Q(theta), S containing 2, infinity and the retained bad support,
write C=Cl(O_K,S). Every rational point has an even-valuation Kummer
class outside S, even if its local class at S is nonzero. Thus the map

    Phi_S : Sel2(E) -> C[2]

is defined on the **whole** generic subgroup, not just G intersect U.
For x=a/d², y=b/d³ on y²=x³+Ax+B, take

    gamma=a-d² theta,       I=(b,gamma).

Away from S, I²=(gamma). Consequently I represents Phi_S([gamma])
in C[2]. Removing S factors changes no class in C. The new calculation
checks that the quotient I²/(gamma) is supported exactly in S.

For each existing strict character chi_i and every point-column j,
including all non-strict generic points, compute

    B_ij = chi_i(Phi_S(gamma_j)).

The [frozen protocol](GENERIC_SUNIT_CARRIER_PROTOCOL.json) uses only the
three earlier controls with complete ideal and Artin data. The result is
entirely retrospective. There are no new points, parameters, integer
factorizations or class-group computations.

| Earlier fibre | Generic rank m | Known strict dimension k | Generic strict dimension g0 | Artin rank on strict half ideals | Artin rank on all generic half ideals |
|---|---:|---:|---:|---:|---:|
| MW16-05 high, 307/206 | 16 | 10 | 1 | 9 | **10** |
| Published R17 high, -2300/843 | 17 | 8 | 2 | 8 | **8** |
| Published R17 low, -1561/3133 | 17 | 6 | 6 | 6 | **6** |

These are exact ranks of the displayed evaluation matrices, hence lower
bounds on the full ideal-class images. The curve ranks remain lower
bounds, and the low fibre remains a censored zero-gain control.

The quotient of the known strict half-ideal image by the detected full
generic image has rank **zero in all three cases**. This is not proof
that every strict half ideal equals a generic half ideal in C: their
difference could lie in the common annihilator of the retained characters.
In particular it does not supply an S-unit correcting any specified
exceptional class.

## An explicit elementary factor from generic ideals

All generic half-ideal classes have order dividing two in C. Choose k
independent columns from the full generic matrix and invert that square
minor. This produces explicit generic ideal words D_j satisfying

    chi_i(D_j)=delta_ij,       2D_j=0.

The [independent certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_generic_sunit_carrier_verification_v1.json)
records each word as a mask in the marked generic basis. Therefore

\[
 C=H\oplus C',\qquad H=\langle D_j\rangle\simeq(\mathbf Z/2)^k,
\tag{2}
\]

with retraction given by the k retained strict characters. This proves
both independence and nondivisibility by two. No full class-group
computation is needed for (2).

For MW16 high, one possible independent set uses zero-based generic
positions 0,1,2,3,4,5,6,8,9,11. Its dual ideal masks are

    547, 2406, 887, 2915, 845, 54, 74, 784, 332, 515.

It gives a ten-dimensional elementary factor. The old right-kernel mask
637 among the *strict* half ideals is still unresolved as a full ideal
class. Its zero Artin evaluations do not become nonzero merely because
new generic columns were supplied. What improved is the detected factor
from the enlarged column population.

This is a useful distinction for class creation: new strict characters
need not come with half ideals outside the span of those attached to
generic points. Class characters, half-ideal classes, and elliptic Kummer
classes are different objects and different maps.

## Proof of the generic-plus-S-unit incidence criterion

Let F_S be the norm-square classes in K*/K*² with even valuation outside
S_K, and let E_S be its subgroup represented by S-units. There is an
exact sequence

\[
 0\longrightarrow\mathcal E_S\longrightarrow F_S
   \stackrel{\Phi_S}{\longrightarrow}C[2]\longrightarrow0.
\tag{3}
\]

For surjectivity, a two-torsion S-ideal class has a representative I with
I²=(alpha) after localization. The norm of alpha has even valuation
outside the rational S primes, so N(alpha)=c*r² for a rational S-unit c.
Replacing alpha by c*alpha makes its norm c⁴r² square, without changing
its localized half-ideal class. The kernel statement is precisely the
S-unit criterion proved in
[the preceding note](STRICT_BLOCKS_NEED_IDEALS_OUTSIDE_BAD_SUPPORT.md).

Put c_S=dim C/2=dim C[2] and e_S=dim E_S. If K has signature (r1,r2),
Dirichlet's S-unit theorem and the split odd-degree norm map give

    e_S = r1+r2-1 + |S_K,finite| - |S_Q,finite|.

Indeed, S-unit squareclasses have dimension r1+r2+|S_K,finite|,
and the rational S-unit group has dimension 1+|S_Q,finite|. Norm is
surjective on these squareclasses because a rational unit has norm
its cube. For the three rows above, e_S is 16,17,12 respectively.

Let lambda be localization at every place in S, including infinity.
Its kernel in F_S is U. The
[strict class-field identification](STRICT_SELMER_AND_ARTIN_BLOCKS.md)
gives dim U=c_S: these are exactly the unramified quadratic characters
of K split at S. Their norm squareclass over Q is automatically trivial;
it is locally square at S and unramified elsewhere. Therefore (3) implies

    dim F_S=e_S+c_S,       rank lambda(F_S)=e_S.

For T=G intersect E_S and b=dim Phi_S(G)=m-dim T, put H0=G+E_S.
Then

    dim H0=m+e_S-dim T=e_S+b,
    rank lambda(H0)<=e_S.

Subtracting proves (1). In particular, with g0=dim(G intersect U),

\[
 \boxed{\dim\frac{(G+\mathcal E_S)\cap U}{G\cap U}
       \ge\max(0,b-g_0).}
\tag{4}
\]

This theorem is equation-defined and contains no exceptional-point
assumption. The retrospective matrices supply b>=10,8,6 in the three
examples, yielding the bounds 9,6,0. They prove those lower bounds using
oracle characters; they do not supply a blind method to certify b.

An explicit S-unit basis would turn the existence statement into a
finite linear construction. For each generic word g, solve

    lambda(u)=lambda(g),       u in E_S.

Then g+u is a strict class. All simultaneous corrections are the kernel
of one joint localization matrix; it is not necessary to solve unrelated
point-search problems. Completeness or a sufficient certified span of
the S-unit squareclasses is the missing arithmetic input. The earlier
null test of individual P² principalizations does not exhaust the mixed
ideal relations needed for that input.

## Factorization-free regularization makes the matrix reproducible

After removing S factors, every retained gcd ideal has cyclic quotient
O_K/I=Z/N. The certificate checks its ring multiplication, root residue r,
f(r)=0 modulo N, and gcd(f'(r),N)=1. Each point element has residue
a-d²r. When this residue is a unit, its Jacobi symbol gives its Artin
evaluation on I.

For a nonunit residue, repeated gcd splits N=UV with gcd(U,V)=1,
the residue a unit modulo U, and every prime of V dividing the residue.
There is no factorization of N. The exact cubic identity is

    (a-d² theta) * (a²+a*d²*theta+d⁴*(theta²+A)) = b².

At primes of V, d is a unit and the second factor reduces to
3a²+A*d⁴, which is a unit because the cubic is étale there. Thus the
quadratic character of the first factor equals that of the second,
even if the first has a positive even valuation. The full symbol is

\[
 \left(\frac{a-d^2r}{U}\right)
 \left(\frac{3a^2+Ad^4}{V}\right).
\tag{5}
\]

The new worker computes all 1490 point symbols, including **70** such
nonunit regularizations. Multiplying them by the retained strict masks
reproduces all **200** entries of the old strict Artin matrices. The
independent verifier checks the polynomial identity, cyclic rings,
ideal transports, and partitions, and uses an integer Jacobi algorithm
instead of PARI's symbol implementation. Both calculations pass.

## Mechanism ranking and the remaining implication

1. **Incidence, now a concrete sufficient condition:** independent
   localized half ideals of the *generic* points force a strict block
   after S-unit correction, by (4). This replaces the stronger and
   unsupported requirement that new strict characters must have new
   half ideals outside the full generic image.
2. **Incidence, not yet prospective:** the elementary part of Phi_S(G)
   is a mathematically interpretable feature. Certifying it from ideals
   requires nonprincipality and, for a direct factor, nondivisibility
   by two. The present certificates use exceptional characters for
   those tests. No score or selector is supplied to Agent 1.
3. **Solubility, still missing:** the constructed-in-principle strict
   classes are Selmer classes. Equation (4) does not put them in the
   rational Mordell–Weil image. CT and the remaining Sha obstruction
   stay separate; no new CT values were computed.
4. **Weak explanation:** S-units alone remain too small in the earlier
   observed blocks, but adding non-strict generic corrections is a
   materially larger constructor. The previous exclusion must not be
   applied to it. The low control also has a substantial generic ideal
   factor; subtracting g0 is essential, and its zero lower bound is
   not an exact zero-gain result.

The next construction target is a sufficiently large certified S-unit
squareclass span, together with its joint local matrix, on a fixed
matched pair. A direct ideal-independence certificate would make (4)
point-independent. Neither a fresh parameter sweep nor another isolated
prime-square generator box addresses that missing input.

The [arithmetic artifact](../../artifacts/generated-results/elliptic-curves/rank_jump_generic_sunit_carrier_v1.json)
and verifier bind their retained sources by hashes. Only rank-jump files
are changed; active searches and mathematical-status entries are untouched.

```sh
timeout 60 sage -python elliptic-curves/rank-jump/verify_generic_sunit_carrier.py check
```
