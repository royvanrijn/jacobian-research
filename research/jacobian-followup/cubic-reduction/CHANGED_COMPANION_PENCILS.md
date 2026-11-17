# Two changed rank-six companion pencils (KCUB21)

Two explicit cubic deformations of KCUB3 give genuinely changed19-variable
homogeneous companions. For every nonzero parameter, each padded tensor has
**derivative-space rank37**, so it needs at least37 shared cubes. There is a
verified98-cube identity at one rational point in each pencil, plus uniform
100- and101-cube identities respectively. These are conditional tensor
bounds and search benchmarks. They do **not** improve the existing
[89-variable cubic-linear bound](POWER89.md), prove a global lower bound,
or prove optimality of either finite search pool.

No dominated98-variable cubic-linear matrix is presented as a size improvement.
The full goal, including smaller dimensions and fewer terms, remains open.

## Exact cubic source changes

Let K be the integral12-variable cubic of [KCUB3](PROOF12.md), with variables
(a,b,c,R,x5,x6,S,x8,x9,T,U,x12)=(x1,...,x12). Define two source automorphisms:

```
E_s: x5 -> x5+s*b*c/2,  U -> U+s*b*R;
E_w: x8 -> x8+w*c*R,    x12 -> x12+w*R^2.
```

Each pair consists of commuting triangular quadratic shears with determinant1.
The full compositions K_s=K E_s and K_w=K E_w are cubic; quartic terms from
the two shears cancel exactly. Their nonidentity changes are

```
K_s:
 row3 += 4s*a*b*c-4s*b*R*S,
 row5 += s*b*c/2,
 row9 += s*a*b*c,
 row11 += s*b*R;

K_w:
 row2 += 12w*a*c*R-12w*R^2*S,
 row8 += w*c*R,
 row10 += w*a*c*R,
 row12 += w*R^2.
```

All other rows agree with K. These are the s- and w-axes of the completed
[KCUB20 family](FULL_CUBIC_CARRIER_CLASSIFICATION.md). The independent checker
reconstructs the source compositions directly from K, without trusting the
listed differences. Both maps are Keller and inherit distinct rational
collision points through E_s^-1 or E_w^-1. They are not new ordinary size
records: the existing12-variable36-term map remains smaller.

For nonzero parameter the cubic output rank is exactly6. In the s-pencil
the active cubic rows are1,...,5,9; in the w-pencil they are1,...,5,10.
Six cubic coefficient pivots have determinant-288 times the parameter,
which is verified exactly. The earlier all-degree one-axis obstruction
KCUB19 does not exclude these coordinated two-axis changes.

## The19-variable companions and collisions

Write K_theta=X+Q_theta+C_theta, and let C_* be its six active cubic rows.
Let iota include those six coordinates in the original twelve. In variables
Z=(X,W1,...,W6,tau), define

```
H_theta(Z)=(tau*Q_theta(X)-tau^2*iota(W), C_*(X), 0),
J_theta(Z)=H_theta(Z)+e19*tau^3.
```

J is the padded tensor used for shared-cube decompositions; the Keller map
uses I+H. Both are homogeneous cubic tensors. In each case J=J0+theta*J1,
with all coefficients retained in the input JSON. For every scalar lambda,
block elimination gives

```
det D(I+lambda H_theta)
 = det(I+lambda*tau*DQ_theta+lambda^2*tau^2*DC_theta)
 = det DK_theta(lambda*tau*X)=1.
```

The last tau coordinate is fixed, so its row first removes that block;
the W block is the identity, and its Schur complement is the displayed
12-by12 matrix. This is the standard companion proof from KCUB5, applied
to the exactly checked new source maps. A19-variable determinant is not
expanded, and no formal proof assistant claim is made.

For any base collision p_theta,q_theta of K_theta,

```
p_hat=(p_theta,-C_*(p_theta),1),
q_hat=(q_theta,-C_*(q_theta),1)
```

collide under I+H_theta. The independent checker evaluates both full
19-coordinate points exactly. At s=-6 and w=-3/2 the complete points and
common images are retained in the validation receipt. These19-variable
examples do not improve the pre-existing18-variable homogeneous construction.

## Exact derivative rank37, uniformly in the parameter

Form all361 partial derivatives of the19 padded tensor outputs, as
quadratic polynomials in Z. Their coefficient matrix has entries affine in
theta. The certificate performs Gaussian row operations over
Q[theta,theta^-1]. Every pivot is a single nonzero rational multiple of a
power of theta; exactly37 pivots leave a zero residual matrix. Therefore
specialization at **every** theta!=0 has rank exactly37, with no hidden
exceptional roots. At theta=0 the independently recomputed rank is35.

- [s-pencil rank certificate](companion-pencil19/derivative-rank-certificate.json)
- [w-pencil rank certificate](companion-w19/derivative-rank-certificate.json)

The independent checker reconstructs every derivative, quadratic column,
row label, pivot unit and elimination step using separate sparse Fraction
arithmetic. It never trusts a primary numerical rank or a generic-field rank
as proof at exceptional parameter values.

If J=sum_i b_i*l_i^3 uses m scalar linear forms, every derivative of every
output is in span{l_i^2}. Its dimension is at most m. Thus m>=37 for both
pencils at every nonzero parameter, over any characteristic-zero field.
The derivative rank is invariant under invertible linear source/output
changes; these nonzero tensors are not merely linear presentations of
the zero-parameter rank35 tensor. This is a lower bound for these tensors,
not for all Keller counterexamples or all cubic-linear maps.

## Verified identities and bounded searches

The existing89-cube tensor identity is first transported exactly to KCUB3's
coordinates and embedded in19 variables. New cubes polarize the changed
monomials and the additional tau^2*W6 coupling. The s-pool has277 forms and
cube rank165; its enlarged mixed pool has427 forms and rank171. The w-pool
has271 forms and rank163.

The basis search treats the parameter as common to every output. In any
independent cube basis, each coefficient row is affine in the parameter.
A row can disappear only at one exact rational cancellation value, unless
it is identically zero. The search scores the generic support minus the
largest common nonzero cancellation group, while allowing bounded temporary
support increases and exact zero-output basis pivots.

| Pencil/pool | Initial | Best | Parameter | Exact exchanges | Budget |
| --- | ---: | ---: | --- | ---: | ---: |
|s /277|104|98|-6|5313|40 seconds|
|s /427 mixed|98|98|-6|2040|40 seconds|
|w /271|101|98|-3/2|5821|40 seconds|

These bounded searches prove no finite-pool minimum. The separately retained
linear/cardinality queries in the277-form s-pool, with mandatory pure-tau
padding, returned **UNKNOWN** after20 seconds at both bounds88 and97.
Neither is UNSAT. The exact executed SMT source is retained separately after
a later fix to its unexecuted SAT numeric-value validation branch. No
unchanged UNKNOWN query was rerun.

The independent checker expands every cube of the two98-form identities
and verifies every output coefficient, the pure-tau padding row and matrix
rank19. It also verifies both coefficients of the entire affine tensor
pencil in the saved uniform identities: at most100 cubes for all s and
at most101 for all w. Thus the least cube count over nonzero parameters
in either pencil lies in[37,98], while the saved uniform bounds are100/101.
This interval is conditional and does not replace the global89 upper bound.

## Replay and next construction gate

```
.venv/bin/python research/jacobian-followup/cubic-reduction/verify_companion_pencils19.py
```

[Independent validation](companion-pencils19-validation.json) checks both
source maps, cubic rank6, full companion formulas, all37 Laurent-unit
pivots per pencil, zero-parameter ranks, all cube coefficients, uniform
pencil identities and full rational base/companion collisions. The pinned
KCUB3 determinant and the written companion Schur identity form the
explicit determinant trust boundary. No CAS/SMT is imported by the replay.

This tests a changed-companion route, rather than another basis exchange
inside the fixed89-form tensor. It has not improved the actual bound.
The next specific gate is the **coupled w,z rank-six component**, not another
unchanged search in these three pools. Its source shears are

```
x8 += w*c*R, x12 += w*R^2,
x9 += z*c*R, U += z*R*(3b-S).
```

The extra cubic rows are proportional to a*c*R, so the rank remains6 when
(w,z)!=(0,0). A normalized companion can use extra cubic output a*c*R and
embedding w*e10+z*e11, keeping its tensor affine in both parameters. In
particular z=-4/3 makes

```
3b^2+3z*b*S-z*S^2=(3b-2S)^2/3,
```

suggesting a grouped polarization and cube forms that follow this repeated
factor. This is a proposed changed model, not a result of the present
search. Other base constructions, dimensions and normalized sparsity remain
open. Current positive bounds stay ordinary cubic12/36, cubic-linear89,
raw3D13 and normalized3D[10,48].

The proposed fixed-z=-4/3 gate is now completed in
[KCUB22](COUPLED_RANKONE_COMPANION.md): the478-form pool stays at101,
with exact derivative ranks39 offzero and37 atzero. This is not a closure
of arbitrary z or of the full coupled component. Consult that note for
the updated next gates rather than repeating this proposed experiment.
