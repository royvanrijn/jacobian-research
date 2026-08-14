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
