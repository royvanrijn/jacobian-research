# H3 q=8 module-intersection continuation — 2026-08-22

This note continues
[`H3_Q8_MODULE_INTERSECTION_2026-08-22.md`](H3_Q8_MODULE_INTERSECTION_2026-08-22.md)
and **supersedes its current-frontier statements where they conflict with this
file**.

Corrected historical q8 state:
[`H3_Q8_CURRENT_FRONTIER.md`](H3_Q8_CURRENT_FRONTIER.md).

The target remains

```text
H3 E7+E8/MW2 --q6--> E8+E6/MW3 --q8--> D13/MW4.
```

The literal q8 source-nef vertical difference is still

```text
(-11,0,2,3,4,6,5,5,6,-4,-5,-7,-10,-8,-6,-4,-2,0,0)
```

relative to `9*O+9*(-P1)`, with complete generic degree-18 frame

```text
B =
1,m,...,m^9,
x,x*m,...,x*m^7,

m=(y-y(P1))/(x-x(P1)).
```

The q8 `-11F` coefficient is global and must be applied exactly once.

## 1. E8 is now closed exactly

This corrects the previous ledger, which still treated resolved E8 saturation
as open.

Use integral II* coordinates

```text
u=1/t,
X=u^4*x,
Y=u^6*y,
Q=u^2*m=(Y-Y(P1))/(X-X(P1)).
```

The horizontal integral frame is

```text
C =
1,Q,...,Q^9,
X,XQ,...,XQ^7.
```

The complete exceptional ideal is

```text
I=(u^2,X,Y),
R/I=<1,u>.
```

With

```text
s=Y(P1)/X(P1),
```

its full preimage in `C` is

```text
u^2,
Q^b-s^b      (1<=b<=9),
X*Q^b        (0<=b<=7).
```

The missing affine saturation argument was completed.  On `u=0`,

```text
Y^2=X^3,
X=r^2,
Y=r^3.
```

Writing

```text
s0=Y(P1)(0)/X(P1)(0),
```

gives

```text
Q0=(r^3-s0^3)/(r^2-s0^2)
  =(r^2+s0*r+s0^2)/(r+s0).
```

The 18 specializations

```text
Q0^0,...,Q0^9,
r^2*Q0^0,...,r^2*Q0^7
```

have rank 18 over `QQ`, so `C` is already saturated on the affine component.

Exact output:

```text
Q8E8RESOLVED|
affine_rank=18|
exceptional_colength=2|
C_over_B_det=178|
reduced_E8_det=180|
saturated=1|
status=PASS_EXACT_Q8_E8_RESOLVED_REDUCED_LATTICE
```

Therefore

```text
reduced resolved E8 determinant valuation relative B = 180.
```

The helper-normalized value `342` is superseded; its extra `162` came from the
q6^9 helper `u^9` factor.

## 2. E7 coefficient provenance correction

One early corrected script still contained a stale assertion because the q6
quotient artifact was misread.

Important:

```text
translated_coordinates["c3"]
```

in the q6 E7 quotient artifact is the `c3` coefficient of the **formal y=0
branch**, not the actual marked section `P1`.

`c2` happens to coincide and is valid.

The actual marked-section coefficients

```text
c3_P1,
c4_P1,
d3_P1,
d4_P1
```

must be expanded directly from the exact rational functions of `P1`.

All valid step-4-and-later calculations use the actual P1 series.

Permanent regression rule:

> never reuse the formal y=0 branch `c3` as `c3_P1`.

## 3. Reduced affine/generic E7 saturation sequence

After stripping the q6 helper fibre normalization, the E7 affine component
satisfies

```text
y^2=x^3,
m=y/x,
x=m^2.
```

Set

```text
z=x-m^2.
```

The exact small branch follows from the chord quadratic

```text
z^2+(m^2+xP)z+2m^2*xP+xP^2+a-2yP*m=0.
```

Writing

```text
z=k2*t^2+k3*t^3+k4*t^4+k5*t^5+...
```

the adapted coordinates used by the successful reverse-saturation sequence
are

```text
z1 = z-k2*t^2,

z2 = m^2*(z1-k3*t^3),

z3 = z2-m^2*k4*t^4,

z4 = m^2*(z3-m^2*k5*t^5).
```

Here

```text
k2=-2*c2
```

and, using actual P1 jets,

```text
z1 = z+2*c2*t^2,

z2 =
m^2*z1
+t^3*(A1-2*d3_P1*m+2*c3_P1*m^2),

z3 =
z2
+t^4*(2*c4_P1*m^2-2*d4_P1*m+A+3*c2^2).
```

### 3.1 Two-prime kernel sequence

The successive one-`t` obstruction kernels at both `43` and `59` were

```text
step 1: 8
step 2: 8
step 3: 8
step 4: 6
step 5: 4
step 6: 2
step 7: 0.
```

At step 5 the affine condition alone would have allowed six directions, but
the resolved exceptional rows become genuinely active and cut to four.

At step 6 only two directions remain.

A convenient terminal rank-18 basis for this calculation is

```text
m^0,...,m^9,
t^-3*z1, t^-3*z1*m,
t^-4*z2, t^-4*z2*m,
t^-5*z3*m^2, t^-5*z3*m^3,
t^-6*z4*m^2, t^-6*z4*m^3.
```

Step 7 has rank 18 / kernel 0 at both primes.

### 3.2 Critical retraction

The previous conversation temporarily called this “reduced E7 saturated
exactly”.  That statement is **retracted**.

The step-1..7 drivers test an affine + generic-component approximation.  They
do not impose the complete finite ninth-power non-Cartier module on the six
actual resolved charts, the marked smooth branch, and overlap gluing.

Therefore these are also retracted as statements about the true q8 line
bundle:

```text
true E7 saturation index = 36
true E7 poles = (3,3,4,4,5,5,6,6).
```

The adapted-coordinate tower and two-prime kernel sequence remain valid
diagnostics for the generic/affine layer.

## 4. Correct marked-E7 fibre normalization

The exact repository marked-chart q8 frame is

```text
m^b/t^6       (0<=b<=9),
x*m^b/t^8     (0<=b<=7)
```

inside the q6^9 helper normalization.

The q6 marked module contains one `t` from its `-F` representative, so its
ninth power contains `t^9`.

The **reduced q8** marked frame, with the helper fibre factor removed, is
therefore

```text
m^b/t^15,
x*m^b/t^17.
```

This is the correction that fixed the global envelope.

The existing generic E7 compiler evaluates

```text
w + (4*k-i-9)*ord(t) + a*ord(x) + b*ord(m).
```

For an actual global coefficient

```text
u^d/h^18
```

use

```text
i=d-9,
k=18,
```

inside that helper compiler.  Then

```text
4*k-i-9=72-d,
```

which is the correct reduced order.

## 5. Superseded 1422-column global assembly

A first global assembly incorrectly used the diagnostic affine/generic E7
pole cap `6`:

```text
h^-18 * sum P_i(u) B_i,
deg(P_i)<=78,
ambient=1422.
```

At both `43` and `59` the full matrix had rank `1422` and kernel zero.

This is **not** a q8 obstruction.

The p43 block diagnostic was

```text
H        rank=1068, kernel=354
E8       rank=378,  kernel=1044
E7diag   rank=72,   kernel=1350

E8+H     rank=1374, kernel=48
E7+H     rank=1140, kernel=282
E8+E7    rank=450,  kernel=972

ALL      rank=1422, kernel=0.
```

The useful localization was

```text
E8+H kernel = 48,
```

showing that the E8 and smooth gates were coherent and that the failure was at
E7.

The same false sheaf had degree

```text
1404-(378+72+1068)=-114
```

instead of the K3 target `-16`, independently diagnosing the E7
normalization error.

The 1422-column envelope and its kernel-zero result are superseded.

## 6. True 1600-column marked global envelope

The worst `h`-denominator in the generic `B` frame is `h^18`.

Since

```text
deg(h^18)=72
```

and the corrected reduced marked E7 poles are `15` and `17`,

```text
m^b families:
  degree <= 72+15 = 87

x*m^b families:
  degree <= 72+17 = 89.
```

Thus the current global ambient has

```text
10*(87+1)+8*(89+1)=1600
```

columns.

The literal q8 `-11F` twist is represented once at E8.

## 7. True p43 global generic-E7 result

The repository's actual seven-component E7 residue machinery was run on this
1600-column ambient with the helper shift `i=d-9`.

Compiler summary:

```text
ambient=1600
negative_groups=230
singleton_groups=15
singleton_coordinate_rank=6.
```

Exact function-field residue rows:

```text
E7_1..E7_3: 1435
E7_5..E7_6:  791
E7_4..E7_7:  255

total non-singleton rows = 2481
+ singleton rows          = 6
total E7 rows             = 2487.
```

p43 ranks:

```text
H:
  rows=1944
  rank=1068
  kernel=532

E8:
  rows=378
  rank=378
  kernel=1222

E7 generic:
  rows=2487
  rank=218
  kernel=1382

E8+H:
  rows=2322
  rank=1446
  kernel=154

E8+H+E7 generic:
  rows=4809
  rank=1582
  kernel=18.
```

The strongest current true global result is therefore

```text
1600
  -- E8 + smooth h -->
154
  -- all 7 generic E7 components -->
18.
```

The full true 1600-column pipeline has not yet been repeated at `59`; defer
that until the finite E7 gate is assembled so the final system can be checked
twice.

## 8. All six leading E7 node screens are vacuous

The exact two-parameter leading bidegree frames were applied to the current 18
global survivors at all six actual nodes:

```text
E7_1--E7_4
E7_4--E7_3
E7_3--E7_7
E7_7--E7_2
E7_3--E7_6
E7_2--E7_5.
```

Only singleton Pareto-minimal negative bidegrees were allowed to cut the
space; no guessed cancellation was used.

Result:

```text
Q8TRUEALLNODES|
prime=43|
generic_survivor=18|
constraints=0|
node_survivor=18|
sweeps=1|
by_node=
E7_1--E7_4:0,
E7_4--E7_3:0,
E7_3--E7_7:0,
E7_7--E7_2:0,
E7_3--E7_6:0,
E7_2--E7_5:0.
```

Therefore the missing q8 E7 information is not visible in unique leading
node monomials.

This does **not** mean the exact node conditions are vacuous.  Equal leading
terms, the surface relation, local units, the marked branch, and overlap
gluing must be handled in the finite local quotient.

## 9. Degree checksum predicts exactly 16 remaining conditions

This is a checksum, not a derivation.

For the isotropic q8 class on a K3,

```text
rank(pi_*O(D))=18,
deg(pi_*O(D))=-16.
```

Known degree contributions relative to `B` are

```text
reduced resolved E8:  -180
smooth h-collisions:  +228
global -11F:          -198.
```

So the complete E7 contribution must be

```text
+134
```

to reach `-16`.

The reduced marked frame permits

```text
10*15 + 8*17 = 286
```

units before finite resolved-chart conditions.  Hence the complete finite E7
correction should have codimension

```text
286-134=152.
```

The true global generic E7 gate already contributes

```text
154 -> 18
```

or `136` independent conditions on the `E8+h` survivor.

The checksum therefore predicts

```text
152-136=16
```

remaining finite E7 conditions and

```text
18 -> 2.
```

This striking match is evidence that the normalization is now aligned, but it
must still be proved by the exact local quotient/gluing calculation.

## 10. Exact next gate already exists in the repository

Do not do more leading-order node screens.

Use

```text
scripts/derive_h92_q8_e7_node_principal_clearings.sage
```

on the true 1600-column ambient.

At each resolved node the q6 marked module has local frame `t*R`, so the q8
condition is

```text
g*f/t^9 in R.
```

The clearing script derives the actual chord in each blow-up chart and clears
only genuine local units.

For the true ambient its common parameters should be

```text
K=18,
T=17.
```

Then use

```text
scripts/probe_h92_q8_e7_node_principal_local_normal_form_modp.sage
```

for each of the six nodes.

That checker uses a Singular local `ds` standard basis and emits the sparse
normal-form image of every ambient column in the genuine local quotient

```text
R/(t^T)
```

for the ideal

```text
(surface,t^T).
```

This is the correct next layer after the vacuous leading screens.

### Immediate algorithm

1. build the six principal-clearings for the true 1600-column ambient;
2. compute each node local-normal-form image modulo `43`;
3. reconstruct its sparse matrix;
4. multiply/restrict it immediately to the current 18-dimensional global
   survivor;
5. combine the incremental node ranks;
6. impose the marked smooth-branch condition and sibling/overlap gluing;
7. target 16 new independent conditions and a two-dimensional kernel;
8. repeat the complete system at `59`;
9. lift the two sections to `QQ`.

Sibling transition data already exists in

```text
scripts/derive_h92_q8_e7_sibling_chart_transitions.sage.
```

Do not replace the six resolved-chart conditions by one complete ideal in the
singular downstairs E7 germ: the integral twist is non-anti-nef.

## 11. Local exploratory driver hashes

These continuation drivers were generated outside the repository.  Hashes are
recorded for reproducibility.

| driver | SHA-256 | result |
|---|---|---|
| `certify_q8_e8_resolved_reduced_lattice.sage` | `95845081d8ccfd06cce62f6fcc30e1c536c8d97f8cc6b6d9ebf636be2e1e728f` | E8 exact, determinant 180 |
| `q8_e7_reduced_normalization_step1.py` | `386d5ec726273bd43806d200f14f00386bbf9d49b9b79994d3a0faa022ded3ee` | kernel 8 |
| `q8_e7_reduced_normalization_step2.py` | `fd7a2aae371947d25cbf0f599f9bd7bc9196f80b450280f9da40e6a57ee50385` | kernel 8 |
| `q8_e7_reduced_normalization_step3.py` | `8f9437f733943be6f0b3396b2e09e4e8faa52e2dc00babe340ec49c47ae5f617` | `z1`, kernel 8 |
| `q8_e7_reduced_normalization_step4_corrected.py` | `01c049bffbabc94b0025649c74399ff6806229a183c62c754a5979dae2c6c504` | actual-P1 `c3`, kernel 6 |
| `q8_e7_reduced_normalization_step5_corrected_v2.py` | `c436ea28f8bd881c889e77e35322f319e335da1c0aa119f5b64184bb3fabfbc5` | kernel 4 |
| `q8_e7_reduced_normalization_step6_v2.py` | `486faf13954cead0a5f6f97c050bb8143fa9494e1c4f50274ae77c271fd66594` | kernel 2 |
| `q8_e7_reduced_normalization_step7_v2.py` | `0807a2a033173917787964561995d3a7a4e84e1cb9bb9648ff17979042f86203` | diagnostic kernel 0 at 43/59 |
| `q8_global_lattice_intersection_modp_v2.sage` | `234589d3b4fa73d360aa749356ced6916e8179fd71f01bd6c9e9db1fb8c09272` | superseded 1422-column kernel 0 |
| `q8_global_lattice_intersection_diagnostic.sage` | `a2221bda24c96e9c318d8c0f1cd42703aac076324f8201d1b4e447aed51b2738` | localized failure to E7 |
| `q8_global_true_e7_generic_gate_v2.sage` | `66e8ab8ad11deeb08bfd7a9b1a270addcb509cec558db7d133fd2b9fcd691c76` | true `1600 -> 154 -> 18` |
| `q8_global_true_e7_all_nodes_gate.sage` | `ec2125b6f35e82ed691d26d1e40501278a2a95763f51372bfa304d12c980b2d9` | six leading node screens all vacuous |

These are exploratory drivers, not repository certificates.

## 12. Retractions / do-not-repeat list

Explicitly superseded:

- promoting local `q=(m-p)/h` to a global coordinate;
- ignoring the degree-10 `d0` transition;
- applying q6 `-F` nine times and q8 `-11F` again;
- using the q6 formal `y=0` `c3` as `c3_P1`;
- treating the affine/generic E7 step-7 kernel zero as complete resolved E7
  saturation;
- using E7 index `36` or poles `(3,3,4,4,5,5,6,6)` as true q8 E7
  elementary data;
- the 1422-column / degree-78 global envelope;
- interpreting its two-prime kernel zero as a q8 obstruction;
- helper-normalized E8 determinant `342`;
- claiming a leading node screen is sufficient after all six were vacuous on
  the true 18-dimensional survivor.

## 13. Claim boundary

### Exact / structurally certified

- q8 source class and generic degree-18 frame;
- smooth h-adic q/X lattice and degree-10 `d0` transition guard;
- q8 E7 and E8 resolved target cycles;
- q8 E7 ninth-power non-Cartier module/gluing specification;
- actual P1 jet provenance correction;
- exact reduced resolved E8 lattice and determinant `180`;
- repository generic E7 actual component/function-field residue machinery;
- corrected reduced marked E7 bounds `15/17`;
- true 1600-column marked global envelope;
- p43 reduction `1600 -> 154 -> 18`;
- vacuity of all six exact singleton Pareto-leading node screens on those 18.

### Checksum / prediction only

- total true E7 degree contribution `+134`;
- finite E7 codimension `152`;
- 16 remaining finite E7 conditions;
- final kernel dimension `2`.

### Open

- exact local-normal-form quotient matrices on the 18 survivors;
- marked smooth-branch and overlap/sibling gluing restrictions;
- final p43 and p59 two-dimensional kernel;
- characteristic-zero lifting;
- q8 pencil and D13/MW4 child;
- downstream rootless H3 equation and bisection-collision work.
