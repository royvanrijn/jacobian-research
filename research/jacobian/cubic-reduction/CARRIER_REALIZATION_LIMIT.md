# Exact limits for changing the twelve-variable cubic carriers (KCUB19)

Fix the integral36-term twelve-variable cubic K of [KCUB3](PROOF12.md).
Two construction routes now have exact exclusions in this displayed frame.
Neither gives a global term or dimension minimum.

1. For any index i and any polynomial h independent of x_i with degree d>=2,
   the one-axis source shear E_i=x_i+h makes an unchanged output of K E have
   degree>3. No subsequent change confined to target coordinate i can fix it.
   This includes all one-axis nonlinear conjugations, in every degree.
2. In the47-coefficient triangular carrier realization specified below,
   every cubic map has exactly the displayed four-parameter form.
   Its minimum nonlinear occurrence count is36, attained only on the
   two-parameter diagonal conjugacy family of the existing K.

The proof uses43 independently multiplied polynomial memberships with123
multiplier monomials, maximum multiplier degree3, and written elementary
elimination and degree arguments. No new counterexample or size bound emerges.

## One-axis nonlinear changes

If i=1, the unchanged sixth output after substitution is
x6+(x1+h)^2, of degree2d. If i=2, the unchanged fifth output is
x5+(x2+h)^3, of degree3d.
For i=3,...,12 choose the unchanged row r listed here. It is affine in x_i,
and its x_i coefficient D_i is independent of x_i and has degree exactly2:

| Source i | Unchanged r | D_i |
| ---: | ---: | --- |
|3|4|x1^2|
|4|2|12*(x10+x1*x8)|
|5|3|8*x1+8*x1*x4|
|6|1|x7+x1*x3|
|7|1|x1^2+x6|
|8|2|12*x1+12*x1*x4|
|9|3|4*x1*(3*x2-x7)+4*x4|
|10|2|12*x1*(3*x2-x7)+12*x4|
|11|3|-4*(x7+x1*x3)|
|12|2|-12*(x7+x1*x3)|

Thus (K E)_r=K_r+h D_i has degree d+2>3: the product of the nonzero
leading forms cannot vanish in a polynomial domain, and deg K_r<=3.
A change of only target coordinate i leaves this row unchanged.
This proof allows arbitrary h, including lower-degree terms; it does not
allow a prior change of frame or coordinated changes to other target rows.

The initial exact screen of all792 monomial quadratic conjugacy lines is
retained in `quadratic-conjugacy12/`. It motivated this stronger written
argument; the all-degree claim does not depend on exhaustive searches in h.
The checker verifies every polynomial coefficient D_i and the two special
base-output formulas directly from K.

## The47-coefficient carrier realization

Write X=(x1,...,x12), with weights

```
(1,-1,-2,0,-3,2,-1,-2,-2,-1,-1,0).
```

K is weighted for these source/output weights. Let C=(K4,...,K12) and fix
x1,x2,x3. These nine carrier outputs give a triangular polynomial
change of the other variables, with unit diagonal.

Construct new carrier outputs C_z, with unchanged linear terms x_i, by
allowing the independent coefficients of the union of:

- every original nonlinear carrier monomial (12 atoms);
- every weight-compatible monomial in x1,x2,x3 of total degree2 or3;
- every weight-compatible quadratic monomial in x1,...,x_(i-1) for output i.

There are exactly47 distinct atoms. Original atoms are multiplied by z_j
relative to their old coefficient; added atoms have coefficient z_j.
`carrier-reparam-wide/input.json` fixes the complete labeled list.
This is a defined triangular support class, not all weighted cubic carrier
polynomials: cubic monomials involving carriers and changes of order/frame
are outside it.

Solve recursively for E_z fixing the first three variables, with

```
C(E_z(X))=C_z(X),       G_z=K(E_z(X)).
```

Every E_z is a triangular automorphism with unit diagonal, E_z(0)=0,
JE_z(0)=I. Therefore G_z is Keller with identity linear part and inherits
the collision through E_z^-1, for every parameter value. No determinant or
collision variables need to be solved again. Being cubic is equivalent to
the vanishing of every coefficient of source degree>3 in G_z.
The independent replay reconstructs E_z, the whole composition and these
93 distinct high-degree coefficient equations from the pinned map.

## Necessary identities and complete elimination

Set

```
A={1,2,3,4,5,6,7,8,9,10,11,12,13,17,18,22,38,42}.
```

The exact membership certificate expresses all29 z_j with j not in A,
and the following14 polynomials, in the ideal of the93 high equations:

```
z3-1, z6-1, z1-z2, z7-z9, z8-z9, z10-z12,
z2*z9-1, z5*z12-z1, z4-z9*z12,
z11-z4-3*z22, z5*z38-z17, z2*z17-2*z13,
z2*z18-z5*z22, z5*z42+z22.
```

No generic parameter-field solve or discarded denominator is used.
From z2*z9=1, z9 and z2 are nonzero; z1=z2 is nonzero, and
z5*z12=z1 makes both z5,z12 nonzero. Put
p=z9, q=z12, s=z17, t=z22. Then pq!=0 and all coefficients are forced to be

```
z1=z2=1/p, z3=z6=1, z4=p*q, z5=1/(p*q),
z7=z8=z9=p, z10=z12=q, z11=p*q+3*t,
z13=s/(2*p), z17=s, z18=t/q, z22=t,
z38=p*q*s, z42=-p*q*t,
all other z_j=0.
```

Conversely the replay substitutes this Laurent parameterization into all93
high equations and verifies zero exactly. Therefore it is the complete
cubic family over every characteristic-zero field. No extra denominator
stratum exists: the two divisions are justified by the exact product
identities above. The written division/elimination proof is not formalized.

The smaller12-atom and22-atom base-only models are retained as discovery
steps. Both give only the p,q diagonal family. Their successful primary
bases are not separate proof premises; they are subsets of the47-atom class.

## Exact support minimum and its scope

After substitution the generic map has45 nonlinear occurrences. Removing
nonzero rational and Laurent monomial factors in p,q from their coefficients
gives the exact factor inventory

| Factor | Occurrences |
| --- | ---: |
|1|33|
|s|4|
|t|5|
|3p+2s|1|
|pq+3t|2|

The independent replay reconstructs all these coefficients and their factors.
Consequently, with brackets denoting indicator functions,

```
N=36+4*[s!=0]+5*[t!=0]-[s=-3*p/2]-2*[t=-p*q/3].
```

If s!=0 its net contribution is at least3, and if t!=0 its contribution
is at least3. Hence N>=36, with equality exactly when s=t=0. There are
no37- or38-term maps in this family; the next possible count is39.
At s=t=0 the entire map is exactly D^-1 K(DX), where

```
D=diag(1,1,1,p,1,1/(p*q),p*q,1,1/p,1/p,1/(p*q),1/(p*q)).
```

This equality is independently checked coefficient by coefficient.
It establishes the minimum only in the displayed47-atom realization,
before arbitrary further changes of source/output coordinates. It does not
bound the affine orbit, arbitrary nonlinear sequences, other triangular
orders, different carrier supports, other constructions or dimensions.

## Verification and changed-frame experiments

Run from the repository root:

```
.venv/bin/python research/jacobian/cubic-reduction/verify_carrier_reparam12.py
```

[Independent validation](carrier-reparam-wide/validation.json) uses sparse
standard-library Fraction arithmetic. It reconstructs all47 atoms and
triangular source changes, all high and low coefficient equations,43
memberships (123 multiplier monomials,maxdegree3), all four-parameter
identities, the exact factor inventory, the minimum diagonal conjugacy
and all one-axis degree witnesses. The base Keller/collision theorem is
imported from KCUB3 with a pinned base validation; it is not re-proved here.
SymPy and Singular are used only by the primary producers/exporter.
The [seven-module helper-domain audit](substitution-domain-audit.json)
rejects Laurent input to the shared polynomial substitution helper and
intercepts checker writes. It caught and led to correction of one earlier
KCUB18 source-normalization audit; the new checker and all six other listed
modules now pass that domain audit. This does not claim repository-wide
helper coverage or formal proof completeness.

Changed raw graph restrictions of all nine carriers in the p=q=1 family
have primary unit ideals for the remaining high-degree coefficients.
These are retained in `carrier-family-graphs12/`, without claiming
independent all-degree target-repair exclusions for the new family.
Six fixed (s,t) frames followed by at most six greedy rational linear
conjugacy passes are retained in `carrier-family-linear12/`. They start
at39,39,42,40,41,45 and end at39,36,39,40,36,40 respectively. No frame
improves36. These primary bounded experiments do not prove the theorem
about the47-coefficient class, nor extend it to arbitrary subsequent frames.

The subsequent [KCUB20 theorem](FULL_CUBIC_CARRIER_CLASSIFICATION.md)
classifies all135 degree2/3 carrier monomials in the original order,
including the cubic atoms outside this earlier47-atom model. Next change
the base/frame/order, use simultaneous nonlinear transformations involving
more than one target row, or use another stable construction. A higher degree cap on
single-axis h does not help. The ordinary cubic bounds12 variables/36
terms, cubic-linear89 variables, raw3D13 terms and normalized3D[10,48]
remain separate. All global simplicity optima remain open.
