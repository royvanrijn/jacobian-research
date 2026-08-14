# Exact rank-lower-bound reproduction of the Fermigier--Mestre family

## Source construction

Fermigier fixes

\[
(a_1,\ldots,a_6)=(0,55,314,378,1007,1036)
\]

and, with \(p_6(x)=\prod_i(x-a_i)\), writes

\[
q(x)=p_6(x-s)p_6(x+s)=g(x)^2-r(x),\qquad \deg_x r=4.
\]

Consequently \(y^2=r(x)\) contains twelve visible affine points with abscissas
\(a_i\pm s\).  The source's two-parameter construction specializes at
\((u,v)=(3,5)\).  Substituting those values in its displayed \(A+B t\) formula
and recovering the fixed-root coordinate gives the material thirteenth point

\[
x_{13}=\frac{1256}{5}-\frac{17}{35}s,
\]

with exact rational ordinate

\[
y_{13}=\frac{50616}{1225}s
 \left(936s^3-254422s^2-283436139s+34925066050\right).
\]

The repository checks the square identity exactly.  The machine-readable
roots, equations, and normalization caveats are in the
[family record](../families/fermigier_mestre_rank12.json).

## Canonical adapter and discriminant

Put \(u=s/2\).  Exact binary-quartic invariants connect the quartic at
\(s=2u\) to the canonical Weierstrass model

\[
C_u=[1,a_2(u),1,a_4(u),a_6(u)]
\]

stored in the family record.  If \(I,J\) are the classical invariants of the
quartic, coefficientwise exact calculation gives

\[
I=101232^4u^4c_4(C_u),\qquad
J=2\,101232^6u^6c_6(C_u).
\]

Thus the raw quartic Jacobian differs by the weighted scale
\((202464u)^{12}\).  The factor \(u^{12}\) and the twelfth-power constant are
nonminimal coordinate artifacts, not usable bad-fiber geometry.  The
canonical discriminant is instead the primitive irreducible even degree-20
polynomial \(\Phi(u)\) whose exact coefficients are stored and checked in the
family record and `ecsearch/fermigier.py`.

The classical binary-quartic covariants give an explicit rational map from
every non-ramification quartic point to the canonical model.  Taking one of
the thirteen images as origin produces twelve section differences.  At the
E22 reconstruction specialization, an exact finite-reduction certificate
proves these twelve differences independent.  A relation between the generic
sections would survive this defined specialization, so this also replays the
generic rank-at-least-twelve claim.

## Exact arithmetic generic rank

The arithmetic generic Mordell--Weil rank of the canonical adapter is
unconditionally **exactly twelve** over \(\mathbf Q(u)\).  This is distinct
from the geometric rank over \(\overline{\mathbf Q}(u)\), which the present
calculation bounds to the interval \([12,13]\).

Indeed, the canonical coefficient degrees are \((0,4,0,8,12)\), and its
degree-twenty discriminant is squarefree and coprime to \(c_4\).  Thus there
are twenty geometric \(I_1\) fibers at finite \(u\).  With \(v=1/u\) and the
K3 scaling \(x=v^{-4}X, y=v^{-6}Y\), the fiber at infinity is

\[
Y^2=X^3-8X^2-64X+512=(X-8)^2(X+8).
\]

The scaled discriminant has order four, the scaled \(c_4\) is \(4096\), and
the node has rational tangent slopes \(\pm4\).  Hence this is a split
\(I_4\) fiber.  The trivial lattice therefore has rank
\(2+\operatorname{rank}(A_3)=5\).  The twelve independently certified
sections give seventeen independent divisor classes defined over
\(\mathbf Q\).

At the good prime \(p=41\), exact character sums on the generalized
Weierstrass equation give

\[
\#S(\mathbf F_{41})=2244,
\qquad
\#S(\mathbf F_{41^2})=2856000.
\]

The corresponding \(H^2\) traces, using
\(\#S(\mathbf F_q)-1-q^2\), are \(562\) and \(30238\).  After removing the
seventeen known \(+41\) eigenvalues, the residual traces are \(-135\) and
\(1661\).  Reciprocity and the exact Weil interval uniquely reconstruct the
residual factor as

\[
(X+41)
\left(X^4+94X^3+4428X^2+158014X+41^4\right).
\]

Its value at \(X=41\) is \(2136275316\ne0\), so there is no additional
\(+41\) eigenvalue.  After normalizing the quartic by \(X=41Z\), its gcd with
every cyclotomic polynomial of degree at most four is one.  Thus the only
residual eigenvalue of the form \(41\zeta\) is \(-41\).

Smooth proper specialization injects characteristic-zero divisor classes,
and a divisor defined over \(\mathbf Q\) contributes a \(+41\) eigenvector.
Consequently an additional \(\mathbf Q(u)\)-section is impossible, and
Shioda--Tate gives arithmetic generic rank exactly twelve.  A geometric
divisor need only contribute \(41\) times a root of unity, so the single
possible \(-41\) class leaves geometric generic rank at most thirteen.  This
argument uses the cycle-class eigenvalue condition, not the Tate conjecture;
it does not decide whether the geometric rank is twelve or thirteen.

The complete replay is
[`verify_fermigier_generic_rank_exact.py`](../cas/verify_fermigier_generic_rank_exact.py),
with its pinned output in
[`elliptic_fermigier_generic_rank_exact.json`](../../artifacts/generated-results/elliptic_fermigier_generic_rank_exact.json).

## Exceptional quotients and two-anchor transport

The exact generic-rank theorem turns the two record-neighborhood fibers into
well-defined exceptional quotients.  At E22, the reconstructed accidental
representatives \(P_{13},\ldots,P_{22}\) give ten independent directions
modulo the generic subgroup.  The remaining representative \(P_6\) satisfies
the exact group-law relation

\[
3P_6-(G_1+G_2+G_3+G_4+G_7+G_8+G_{11})
+2(G_5+G_6+G_9+G_{10})-3G_{12}=O.
\]

At the low-conductor rank-20 anchor, the twenty-point certificate decomposes
as the twelve generic directions plus eight independent exceptional
directions.  This makes the comparison intrinsic rather than dependent on a
raw parameter string or a numerically selected point count.

Every one of the \(11\cdot8=88\) affine abscissa interpolants between the
published E22 accidental representatives and the eight rank-20 directions
gives an irreducible squarefree sextic, hence a genus-two square condition.
For the complete one-parameter quadratic interpolation through each pair,
the branch discriminant has signature \(k^{16}\) times an irreducible
degree-32 factor; its only rational collision is \(k=0\), which returns the
same genus-two affine member.

The complete projective Möbius pencil

\[
x(T)=\frac{aT+b}{cT+d}
\]

through each pair was also classified exactly.  In the finite chart \(d=1\),
the degree-72 branch discriminant has factor signature

\[
(1,12),(1,12),(1,16),(32,1).
\]

The two multiplicity-12 rational roots put a pole at one of the prescribed
anchors and are invalid interpolants; the multiplicity-16 root is \(c=0\),
the affine genus-two limit.  The degree-32 factor is irreducible over
\(\mathbf Q\), and the missing \(d=0\) chart is an irreducible degree-ten
genus-four member.  Thus none of these 88 complete Möbius pencils contains a
genus-zero or genus-one transport.

Finally, the calculation constructs the actual covers

\[
y_1^2=f_1(T),\qquad y_2^2=f_2(T),
\]

for all \(\binom{80}{2}=3160\) unordered pairs of independent affine
transports, including pairs sharing one anchor endpoint.  Every pair has
disjoint branch loci; its third character quotient has genus five and the
connected fiber product has genus nine.  This avoids the weaker
product-is-a-square surrogate, which can accept points for which neither
factor is a square.

These are exact finite classifications of the stated affine, quadratic,
Möbius, and pair-product ansatzes.  They do not exclude higher-degree
multisections, higher-support representatives beyond the finite ball below,
or rational points on the resulting higher-genus curves.  The replay and its
machine-readable manifest are
[`analyze_fermigier_exceptional_transport.py`](../cas/analyze_fermigier_exceptional_transport.py)
and
[`elliptic_fermigier_exceptional_transport.json`](../../artifacts/generated-results/elliptic_fermigier_exceptional_transport.json).

The representative choice has also been enlarged exactly.  In the E22
exceptional quotient, all 200 signed vectors of support at most two in the ten
element basis were enumerated; the analogous eight-element rank-20 quotient
gives 128 vectors.  Every one of the 328 transported representatives
round-trips through the pointed quartic, agrees with the canonical group law,
and has a mod-5 certificate column independent of the twelve generic columns.
All 25,600 cross-anchor affine interpolants are irreducible squarefree
sextics, hence genus two.  Thus the affine negative result is independent of
the original representative choice throughout this complete signed
support-at-most-two ball.  The exact replay is
[`classify_fermigier_exceptional_quotient_ball.py`](../cas/classify_fermigier_exceptional_quotient_ball.py),
with artifact
[`elliptic_fermigier_exceptional_quotient_ball.json`](../../artifacts/generated-results/elliptic_fermigier_exceptional_quotient_ball.json).

A first genuinely nonlinear interpolation was classified for the independent
pair \(P_{13}\times R20E1\).  In the finite denominator chart the complete
bidegree-\((2,1)\) pencil is

\[
x(T)=\frac{a(c)T+b(c)+k(T-T_{22})(T-T_{20})}{cT+1},
\]

where \(a(c),b(c)\) enforce both endpoint values.  Its generic squareclass
kernel has degree ten and genus four.  The exact branch discriminant has one
irreducible total-degree-32 factor and five rational lines.  Two lines are
invalid anchor poles; the cancellation line returns the known affine
genus-two sextic; and the two valid degree-drop lines have squarefree
degree-eight kernels of genus three.  The nonlinear component generically has
kernel degree eight and genus three.  Hence this complete rational-component
classification contains no genus-at-most-one member.  Rational special points
on the degree-32 component were outside this pilot; the other 79 independent
pairs are handled by the complete replay below.
The bounded exact pilot is
[`analyze_fermigier_bidegree21_pilot.py`](../cas/analyze_fermigier_bidegree21_pilot.py),
with artifact
[`elliptic_fermigier_bidegree21_p13_r20e1_pilot.json`](../../artifacts/generated-results/elliptic_fermigier_bidegree21_p13_r20e1_pilot.json).

The same classification has now been completed for all
\(10\cdot8=80\) independent cross-anchor pairs.  Across their 400 rational
discriminant components, every valid kernel has genus two or three.  The 80
residual total-degree-32 factors are all irreducible over \(\mathbf Q\): 79
have an irreducibility witness modulo 101, while the one collision at that
prime is witnessed modulo 103.  There are no unresolved pairs or
characteristic-zero fallback assumptions.  This exhausts the rational
components in the declared finite bidegree-\((2,1)\) charts; it does not
classify rational points on the irreducible degree-32 curves or their
intersections.  The complete replay is
[`analyze_fermigier_bidegree21_all80.py`](../cas/analyze_fermigier_bidegree21_all80.py),
with artifact
[`elliptic_fermigier_bidegree21_all80.json`](../../artifacts/generated-results/elliptic_fermigier_bidegree21_all80.json).

For the pilot pair's irreducible degree-32 component, a separate exact
projective sieve covers every primitive \((C:K:D)\) with \(D\ge0\) and
\(\max(|C|,|K|,D)\le1024\).  Five fixed residue primes reduce
2,098,176 \((C,D)\) pairs to 21,819 primitive exact evaluations, with no
affine point.  The top binary form is irreducible of degree 32, so there is
no rational point at infinity; exact restrictions to the five known lines
also have no rational intersection.  This is bounded negative coverage of
one nonlinear component, not a rational-point theorem for that curve or the
other 79 components.  The replay is
[`search_fermigier_bidegree21_nonlinear_points.py`](../cas/search_fermigier_bidegree21_nonlinear_points.py),
with artifact
[`elliptic_fermigier_bidegree21_p13_r20e1_nonlinear_points_h1024.json`](../../artifacts/generated-results/elliptic_fermigier_bidegree21_p13_r20e1_nonlinear_points_h1024.json).

The 3,160 genuine pair covers were also searched directly, testing the two
square conditions separately rather than using their product.  In the exact
projective box \(T=a/b\), \(b>0\),
\(\max(|a|,b)\le200000\), every simultaneous intersection contains only the
two prescribed anchors.  Two isolated parameters occurred on one cover each,
but failed the second square test.  This is a complete bounded negative result
for that box, not a global rational-point theorem on the genus-nine covers.
The replay and manifest are
[`search_fermigier_exceptional_pair_simultaneous_h200000.py`](../cas/search_fermigier_exceptional_pair_simultaneous_h200000.py)
and
[`elliptic_fermigier_exceptional_pair_simultaneous_h200000.json`](../../artifacts/generated-results/elliptic_fermigier_exceptional_pair_simultaneous_h200000.json).

## Exact independence certificates

For points \(P_1,\ldots,P_n\), the checker uses cyclic good reductions
\(E(\mathbb F_p)=\langle G_p\rangle\).  If the matrix of discrete logarithms
\(P_j=k_{p,j}G_p\) has rank \(n\) modulo a prime \(\ell\), every integral
relation among the \(P_j\) has coefficients divisible by \(\ell\).  A separate
good-reduction prime where \(\ell\nmid\#E(\mathbb F_p)\) rules out rational
\(\ell\)-torsion.  Infinite descent then proves independence.

The pinned replay uses \(\ell=5\) for the twelve section differences and
\(\ell=2\) for all twenty-two points printed for E22.  For the latter,
\(\#E(\mathbb F_{31})=41\) excludes rational 2-torsion, and 22 cyclic
reductions give a full-rank parity matrix.  The standard-library checker
replays point counts, generator orders, scalar multiples, and matrix ranks;
PARI/GP independently replays every stored finite-field equality.  Thus E22's
rank-at-least-22 lower bound is local exact evidence, not merely a numerical
height determinant.

## Factor-two source discrepancy

The paper prints both the product with shift \(s=t\) and the E22 parameter
\(t=19754/39\).  Literal exact substitution gives the minimal model

```text
[1,0,1,
 -1223348097402005168062873899944,
 -213263015130965060475376699543914227367884158]
```

and conductor

```text
3336936695055698757544757721801363002721636124100955091377369358240007970.
```

It does not give the displayed E22 curve.  Substitution at the doubled shift
\(s=39508/39\), equivalently adapter coordinate \(u=19754/39\), gives exactly
Fermigier's displayed minimal model and conductor.  The two literal
specializations have different exact \(j\)-invariants.  No intervening change
of parameter or published erratum was found, so this remains an explicit
reproduction discrepancy rather than a resolved normalization.

## Search and evidence boundary

`verify_benchmarks.py` checks the invariant bridge, both literal
specializations, the exact E22 model and conductor, and the literal
`log(N)<182.72` cutoff.  `verify_family_data.py` cross-checks the stored family
equation, discriminant, and thirteenth-point metadata.
`verify_fermigier_rank_certificates.py` replays both independence certificates.

For a new parameter, `evaluate_fermigier_specialization.py` reconstructs the
twelve baseline differences and can call PARI's `hyperellratpoints` for a
bounded quartic search.  It can select and exactly certify a modularly
independent subset of the resulting point cloud.  The height bound is not a
complete point search, and failure of a modular test does not prove
dependence.  The certificates do not prove saturation or an upper rank bound;
E22's exact-rank statement remains conditional as recorded in the benchmark
metadata.  The factor-two source normalization discrepancy also remains open.

## Low-conductor rank-20 near miss

The same evaluator finds, at adapter parameter (u=28917/20), twenty points
whose finite-reduction matrix has full rank modulo 5.  PARI/GP gives the global
minimal model

```text
[1,1,1,
 -4437412060110743641525245114305,
 3586842216822165612930264910099076801587288127]
```

and exact conductor

```text
2876153493562761211278364526603564191699143885403233935132057708367930,
```

whose natural logarithm is
`159.9348252255254533984... < 182.72`.  The
[pinned certificate](../../artifacts/generated-results/elliptic-curves/fermigier_rank20_near_miss_v1.json)
records the bounded `ratpoints` box, all retained abscissas, the exact minimal
model and conductor, and the twenty finite-reduction rows.  It proves rank at
least 20, not rank 21; it supplies no upper bound or saturation theorem.

The specialization now also has a canonical candidate record at
[`elliptic_curve_candidate_fermigier_mestre_v1_u28917_20.json`](../../artifacts/generated-results/elliptic_curve_candidate_fermigier_mestre_v1_u28917_20.json).
Its stable identity is `fermigier-mestre-v1:u=28917/20`; the literal shift
`T=28917/10` and their negatives are aliases, not separate candidates.  The
record joins the generalized, quartic-Jacobian, short, and global-minimal
models and their exact transformations; the complete discovered point pool;
the selected basis and certificate hashes; bounded saturation evidence; and a
promotion/rejection ledger covering every recorded search region.  The
bounded saturation result is labeled as such and is not a global saturation
theorem.  Final rank claims retain both the imported cyclic-log verifier and
an independent quotient-enumeration implementation.

A separate sinc-squared explicit-formula diagnostic uses
\(\Delta=11/5\) and every prime through
\(\lfloor\exp(22\pi/5)\rfloor=1007525\).  Its conservative value, including
the declared `0.001` numerical allowance, is
`21.0335328229846198389...<22`.  Since the exact root number is `+1`, GRH
would force analytic rank at most 20; BSD together with GRH would therefore
make the algebraic rank exactly 20.  This is a conditional fixed-fiber
closure, not an unconditional rank upper bound.  It redirects the search to
nearby parameters rather than promoting more bounded point searches on this
one curve.
