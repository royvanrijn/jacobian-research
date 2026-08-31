# First-jet elimination for the five-fibre `wgxli` target

Status: **exact finite-field necessary-condition elimination in one complete
distinct-parameter chart; literal rootless-K3 interpolation rejected modulo
17 unless the unknown model has bad or colliding reduction**.

This is not a characteristic-zero nonexistence theorem.  It assumes that the
first seventeen public points on curves 351, 356, 376, 377, and 385 are
literally corresponding sections, with their displayed signs, in the
rootless-K3 degree bounds.

## Outcome

All polynomial coefficients of the proposed surface and sections can be
removed before a nonlinear solve.  After the gauges

```text
t_351=0,  t_356=1,  t_385=-1,  u_351=1,
```

the only intrinsic continuous unknowns are

```text
t_376, t_377, u_356, u_376, u_377, u_385.
```

The exact first-jet equations introduce ten linear auxiliary values
`A'(t_k),B'(t_k)` but no unknown section or surface coefficients.  Exhausting
the complete distinct-node chart over `GF(17)` gives

```text
ordered (t_376,t_377) pairs tested: 182
pairs with a geometric solution:       0
timeouts:                              0
```

Each fixed-pair ideal was solved over the algebraic closure, not merely by
enumerating `GF(17)`-rational points.  Thus a literal characteristic-zero
model with integral good reduction of the parameters and scalings at 17
would force a collision among the five reduced base parameters.  The result
does not exclude such a collision or another bad denominator at 17.

## Removing the 52 unknowns

Write the five canonical short fibres as

```text
Y^2 = X^3 + A_k X + B_k
```

and put `w_k=u_k^-1`.  Family coordinates at the proposed parameter `t_k`
are

```text
x_i(t_k)=w_k^2 X_{i,k},       y_i(t_k)=w_k^3 Y_{i,k},
A(t_k)=w_k^4 A_k,             B(t_k)=w_k^6 B_k.
```

Let

```text
L(t)=product_k (t-t_k).
```

Five values determine the quartic `x_i(t)` uniquely.  Every sextic ordinate
through its five values has the unique form

```text
y_i(t) = ybar_i(t) + L(t)*(r_i+s_i*t),
```

where `ybar_i` is the degree-at-most-four Lagrange interpolant.  This already
isolates the 34 free ordinate coefficients as the pairs `(r_i,s_i)`.

Differentiate the section identity

```text
y_i^2=x_i^3+A*x_i+B
```

at each fibre.  With

```text
alpha_k=A'(t_k),   beta_k=B'(t_k),
c_ik=2*y_i(t_k)*L'(t_k),
d_ik=2*y_i(t_k)*ybar_i'(t_k)
     -3*x_i(t_k)^2*x_i'(t_k)-A(t_k)*x_i'(t_k),
R_ik=x_i(t_k)*alpha_k+beta_k-d_ik,
```

the differentiated identity is simply

```text
c_ik*(r_i+s_i*t_k)=R_ik.
```

Because the right side divided by `c_ik` must be affine in `t_k`, the two
ordinate unknowns disappear.  Using the fibres at `0` and `1`, each of the
other three fibres gives

```text
c_i0*c_i1*R_ik
 - c_ik*((1-t_k)*c_i1*R_i0+t_k*c_i0*R_i1) = 0.
```

There are three such equations for each of seventeen sections, hence 51
first-jet equations.  The five values and five `alpha_k` have a unique
degree-at-most-nine Hermite interpolant; `deg(A)<=8` removes its leading
coefficient.  One saturation equation makes the four remaining `w_k`
nonzero.  At 17 the Hermite leading equation vanishes identically after
reduction, leaving 52 nonzero equations in each chart.

Relative to the earlier 52-variable interpolation count, this removes

```text
34 free section-ordinate coefficients,
12 free coefficients of A and B,
```

and enumerates the two base parameters.  Each modular chart consequently has
only four nonlinear scaling variables, ten linear derivative auxiliaries,
and one saturation inverse.

## Positive control and interpretation

The same differentiated identities replay exactly on all seventeen sections
of the certified published R17 surface at the five control nodes

```text
0, 1, 2, 3, -1  (mod 17).
```

The eliminator is therefore detecting a failure of the literal five-fibre
input, not an algebraic error visible on the known family control.

The mod-17 result leaves three live explanations:

1. the unknown rational parameters or scalings have bad/colliding reduction
   at 17;
2. the public points preserve a nearby rank-17 lattice basis but not literal
   corresponding sections;
3. the common lineage is not a rootless K3 with section degrees `(4,6)`.

## Reproduction

From the repository root:

```bash
sage -python \
  elliptic-curves/cas/eliminate_icarm_wgxli_rank17_first_jet.sage \
  --prime 17 --threads 4 --pair-timeout 10 \
  --output \
  artifacts/generated-results/elliptic-curves/icarm_wgxli_rank17_first_jet_mod17_v1.json

sage -python \
  elliptic-curves/cas/eliminate_icarm_wgxli_rank17_first_jet.sage \
  --prime 17 --threads 4 --pair-timeout 10 --check
```

The replay requires Sage and `msolve`.  The generated artifact has SHA-256

```text
de1ebb881e881f29a4e2850e6d75bb234f6a5acd1e55ffb8e87aed8b786f0971
```

and records the hash of every one of the 182 generated solver inputs.
