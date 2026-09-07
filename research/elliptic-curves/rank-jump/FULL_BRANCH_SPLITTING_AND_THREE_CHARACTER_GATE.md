# Full branch splitting and a three-character Picard gate

This is a bounded, coefficient-only test on the retained rank-20 anchor.
It completes the generic-rank calculation left open by
[the norm-square analysis](NORM_SQUARE_EVENT_AND_EXTREMAL_K3.md).
It does not recompute the ranks of production fibres or change their status.

**Result.** Adjoining all three square roots `sqrt(1-u*theta_i)` gives
geometric generic Mordell–Weil rank **1**, entirely inherited from the original
pencil, and arithmetic generic rank **0**. The three previously unresolved
character K3 surfaces all have geometric Picard rank **19**, hence
Mordell–Weil rank **0**. Two small reductions and independent point counts
certify this.

The cover has infinitely many rational parameters. Consequently the absence
of new generic directions is not explained by an empty rational parameter
locus. It is a statement about function fields; ranks of its individual
specialized elliptic curves remain unbounded by this calculation. Possible
finite-index changes of the inherited lattice have not been computed.

## 1. Objects and rank accounting

Use the pinned coefficients in
[the input projection](../../artifacts/generated-results/elliptic-curves/rank_jump_local_collision_inputs_v1.json):
\[
\begin{aligned}
A&=-5750886029903523759416717668139307,\\
B&=167347710468055045100164888198438918505621536951206,\\
f(t)&=t^3+At+B,\quad \delta=-4A^3-27B^2,\\
\gamma_i&=1-u\theta_i,\quad D=\prod_i\gamma_i=1+Au^2+Bu^3.
\end{aligned}
\]
Here `f` is irreducible with Galois group `S3`, `B` is nonsquare, and
`B*delta != 0`. The deformation is
\[
E_u:\ y^2=x^3+2Au x^2+(A+3Bu+A^2u^2)x+B+ABu^2-B^2u^3.
\]
Its roots are \(\theta_i+u\theta_i^2\) and its discriminant is
\(16\delta D^2\). Its arithmetic generic rank is zero, although the
special fibre at `u=0` has the retained independent rank-20 subgroup.
Those twenty directions are not generic sections of this pencil.

Over \(\overline{\mathbb Q}(u)\), let
\(L=\overline{\mathbb Q}(u)(\sqrt{\gamma_1},\sqrt{\gamma_2},\sqrt{\gamma_3})\).
The three private simple zeroes make their squareclasses independent.
The rank over `L` is the sum of the ranks of the eight quadratic characters
\(d_I=\prod_{i\in I}\gamma_i\).

| Character weight | Number | Singular fibres | Trivial lattice rank | Geometric MW rank |
|---|---:|---|---:|---:|
| 0 | 1 | `3 I2 + I0*` | 9 | 1 |
| 1 | 3 | `I2* + 2 I2` | 10 | 0 |
| 2 | 3 | `2 I2* + I2 + I0*` | 19 | 0, certified below |
| 3 | 1 | `3 I2*` | 20 | 0 |

For weights zero and one the minimal surface has \(\chi(\mathcal O)=1\)
and is rational, so its geometric Picard rank is 10. For weights two and
three it is K3. The weight-three trivial lattice already has rank 20.
The symbolic verifier checks discriminants, vanishing orders, infinity
degrees, and this bookkeeping. The fibre/root-lattice and Shioda–Tate
inputs are standard; see
[Schütt–Shioda, ``5–6](https://arxiv.org/pdf/0907.0298).

The inherited rank-one space is generated over the algebraic closure by
\((Bu^2,\sqrt B D)\). Its Galois action has the nontrivial sign character of
\(\sqrt B\). The rational cover below is geometrically connected, hence
does not add constants. Once all other character spaces vanish, its
arithmetic generic rank is therefore zero.

## 2. The two-place Picard certificate

Fix \(K=\mathbb Q(\theta)\) and the weight-two twist
\(d=D/(1-u\theta)\). All three such surfaces are Galois conjugates.
The [frozen protocol](BRANCH_SPLIT_PICARD_PROTOCOL.json) chose the first two
good completely split primes above 3, at most 59: 23 and 59. It used every
root at each prime and extensions of degrees one and two. No exceptional
points entered this selection.

These are good reductions of the resolved K3 surfaces: the residue
characteristics exceed 3, \(B\delta\) is a unit, and the three branch roots
remain distinct and nonzero. The same minimal Tate configurations resolve
over these local rings. All nineteen trivial classes are defined over the
residue field. In particular, full rational 2-torsion injects into the
order-four component groups of the additive fibres (the identity component
has no prime-to-characteristic 2-torsion). It fixes their outer components
and hence the entire `D6` or `D4` diagram. The `A1` component is also fixed.

For \(q=p,p^2\), write
\[
s(u)=\sum_{x\in\mathbb F_q}\chi(F_u(x)),\qquad
T_q=\sum_{u\in\mathbb F_q}\chi(d(u))s(u).
\]
Here `F_u` is the right side of the original cubic equation.
Resolving the two `I2*` fibres adds `12q` points, the `I2` adds `q`,
and infinity adds `4q`. Thus
\[
\#X(\mathbb F_q)=q^2+19q+1+T_q.
\]
The three-dimensional orthogonal complement of the trivial lattice in
\(H^2\) consequently has trace \(T_q\).

Its eigenvalues have the form \(\epsilon p,\alpha,\beta\), with
\(\epsilon=\pm1\) and \(\alpha\beta=p^2\). If \(t=\alpha+\beta\), then
\(T_p=\epsilon p+t\) and \(T_{p^2}=t^2-p^2\). Both traces determine
the unique sign in each retained row. The two decisive places of the
**same surface over K** are:

| Place | \(T_p\) | \(T_{p^2}\) | \(\epsilon p\) | Residual quadratic | NS discriminant squareclass |
|---|---:|---:|---:|---|---:|
| \((23,\theta-1)\) | −7 | 371 | 23 | \(X^2+30X+529\) | −19 |
| \((59,\theta-32)\) | −85 | −2805 | −59 | \(X^2+26X+3481\) | −23 |

In both cases the residual quadratic is noncyclotomic after normalization
by `p`: its integer trace is outside
\(\{-2p,-p,0,p,2p\}\). Hence the geometric Picard rank of each reduction is
20 by Tate for elliptic K3 surfaces. Over \(\mathbb F_{p^2}\) all twenty
algebraic classes are fixed. The Artin–Tate formula gives the **signed**
Néron–Severi discriminant squareclass
\[
-\frac{4-t^2/p^2}{p^2}\equiv t^2-4p^2
\pmod{\mathbb Q^{*2}}.
\]
The Brauer-group order is a square. For the reduction/squareclass method,
see [van Luijk, `2](https://arxiv.org/pdf/math/0506416) and
[Schütt, `7.2](https://arxiv.org/pdf/1202.1066).

If the characteristic-zero Picard rank were 20, specialization at both
places would embed lattices of equal rank, so their discriminants would
agree modulo squares. They do not. The known trivial lattice supplies
rank 19, proving \(\rho=19\). Galois conjugacy gives the result for all
three character surfaces.

All six root/prime rows are retained. Two rows have cyclotomic residual
quadratics and are labelled `UNKNOWN` for this rank-19 endpoint; they
are not used as witnesses. The other root at each prime independently
gives the same discriminant squareclass as its selected witness.

## 3. The full splitting locus is a positive-rank genus-one curve

A rational square root of \(\gamma=1-u\theta\) in `K` can be written
\(z=(a+b\theta+c\theta^2)/h\). Its coefficient equations are
\[
C_\gamma:\quad a^2-2Bbc=h^2,\qquad b^2+2ac-Ac^2=0,
\]
with
\[
u=-\frac{2ab-2Abc-Bc^2}{h^2}.
\]
The smooth projective curve is connected of degree 8 over the `u`-line.
Its four branch inertia vectors are \(100,010,001,111\), so
Riemann–Hurwitz gives genus one. The point \((1:0:0:1)\) is rational
and smooth. The affine equations and tangent rank are checked exactly.

The norm map to \(C_D:v^2=D(u)\) has degree 4 and is unramified:
both covers ramify with index 2 at precisely those four branch points.
After translating the image of the rational origin to zero it is an
isogeny. Its geometric kernel is the even-sign `V4`, hence the entire
2-torsion. Factoring through multiplication by 2 gives a rational
isomorphism \(C_\gamma\simeq C_D\). An explicit formula for that
isomorphism is not required or supplied here.

On \(E_D:Y^2=X^3+AX^2+B^2\), put \(Q=(0,B)\).
The change \(X=Bu,\ Y=Bv\) identifies it with `C_D`. Exact addition gives
\[
2Q=(-A,-B),\qquad
3Q=(4B^2/A^2,-B-8B^3/A^3).
\]
The reductions of `Q` at the good primes 23 and 59 have orders 12 and
18. If `Q` were torsion, injectivity of prime-to-residue-characteristic
torsion at both primes would force its order to divide 6. But `3Q` is
finite and has nonzero ordinate, so `6Q != O`. Thus both parameter curves
have positive rational rank, and their maps supply infinitely many
rational full-splitting parameters.

A nonzero coefficient-only example is
\[
u^\dagger=4B/A^2,\quad z=1+2\theta^2/A,\quad
z^2=1-u^\dagger\theta,\quad N(z)=1+8B^2/A^3\ne0.
\]
This certifies full branch splitting without supplying exceptional points
on \(E_{u^\dagger}\). Its Mordell–Weil rank is `UNKNOWN`. The norm point
has the negative of the ordinate of `3Q`; both have the same `u`.

## 4. What mechanism survives?

There is a precise **geometric incidence** gate for any characteristic-zero
irreducible `S3` anchor with the same nondegenerate configuration:
\[
\operatorname{rank}E(L)=1+3(\rho-19),\qquad \rho\in\{19,20\},
\]
where \(\rho\) is the common geometric Picard rank of the three weight-two
K3s. A Picard increase to 20 would therefore create a block of three
geometrically independent directions, one in each character space.
The retained anchor fails this gate: \(\rho=19\).

This is a condition on the anchor/surface coefficients, not a condition
varying with `u` inside this one pencil. It is not yet an explanation for
the large jumps in the R17/MW17 and A1/MW16 production fibres.

The mechanism ranking after this falsifiable experiment is:

1. **Incidence candidate:** a three-character Picard increase gives a
   proved geometric block implication. Missing: an explicit coefficient
   criterion, a relevant anchor satisfying it, rational Galois invariants
   of the block, and independence after specialization.
2. **Solubility candidate:** the strict cup-product/CT obstruction from
   [the previous paired analysis](CUP_IDEAL_AND_STRICT_LIFTING_OBSTRUCTION.md).
   It distinguishes rational and Sha classes with the same ordinary
   class-field data. Missing: independent production norm witnesses and
   the remaining local corrections, followed by a rational-point theorem.
3. **Weak explanation:** a smaller torsion field, square norm, or even
   full branch splitting alone. None creates generic directions here.
   No statement about their statistical association with specialized
   high rank has been proved.
4. **Visibility only:** chart exposure and half-lattice recovery. This
   calculation provides no new visibility feature or search-coordinate
   score.

Agent 1 could eventually use a certified Picard/character construction
to select **families**, after rational descent is solved. It cannot use
this rank-zero result to veto individual fibres. No selector, scoring
policy, candidate population, or worker setting is changed.

## 5. Reproduction and evidence

The [raw finite-field arrays](../../artifacts/generated-results/elliptic-curves/rank_jump_branch_split_picard_inputs_v1.json)
and [Picard certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_branch_split_picard_v1.json)
bind the input, protocol and counting source by SHA-256.
The primary count uses direct NumPy character sums. The
[independent replay](../../artifacts/generated-results/elliptic-curves/rank_jump_branch_split_picard_verification_v1.json)
uses Sage 10.9 / PARI 2.17.3 elliptic cardinalities for every smooth fibre,
up to explicitly paired degree-two Frobenius orbits, and exact nodal
counts at singular fibres: 2,116 smooth cardinality calls and 12 singular
checks, covering all 4,092 base parameters across the four fields.

The [geometry certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_branch_split_geometry_v1.json)
records all eight character configurations, the full-split example, and
the auxiliary non-torsion proof inputs. A root-order preflight failure
and [two coefficient-domain verifier failures](../../artifacts/generated-results/elliptic-curves/rank_jump_branch_split_geometry_attempts_v1.json)
retain their source and explanation; no failed attempt supplies evidence
for the theorem.

Replay from the repository root (no point searches or new parameters):

```sh
python3 elliptic-curves/rank-jump/branch_split_picard.py check
python3 -m unittest discover -s elliptic-curves/rank-jump -p test_branch_split_picard.py
sage -python elliptic-curves/rank-jump/verify_branch_split_geometry.py
sage -python elliptic-curves/rank-jump/verify_branch_split_picard.py --prime 23
sage -python elliptic-curves/rank-jump/verify_branch_split_picard.py --prime 59
```

Construction commands refuse to overwrite certificates. The frozen
calculation used two primes, degrees at most two, a 30-second cap per
prime, and checkpoints after each field. No class group, global descent,
prospective parameter sweep, or rational point search was run.
