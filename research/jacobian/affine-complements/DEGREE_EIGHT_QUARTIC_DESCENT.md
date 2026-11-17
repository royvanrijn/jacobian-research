# Sign descent of the rooted degree-eight core to quartic fields

14 November 2026. PF2RULED30 reduces the coefficient field needed for
an exact degree-eight certificate from degree eight to degree four.
It transports the complete quadratic presentation and the complete
four-parameter core by invertible constant changes of coordinates.
The existing degree-three unit rows also descend exactly.

Both signs of each marked root are covered. Thus the sixteen octic
embeddings give eight quartic embedding problems. No nonzero-Jacobian
candidate is excluded by this descent. The characteristic-zero unit
certificate for the full core remains missing.

## The coefficient subfields

The two marked polynomials are even. With t=y^2 their quartic polynomials
are

    sector (2,7): t^4 - (32/25)t^3 + (32/7)t^2
                        - (384/35)t + 8192/1225,
    sector (4,5): t^4 + 112t^3 - (896/5)t^2
                        + (384/5)t - 256/25.

In each case let K=Q[y]/P8 and F=Q[t]/P4. Both defining polynomials
are irreducible, their constant terms are nonzero, and the embedding
F -> K sends t to y^2. Consequently [K:F]=2. The nontrivial automorphism
sigma sends y to -y and fixes exactly F. In the power basis of K,
a coefficient belongs to F precisely when its odd powers of y vanish.
The checker tests this and the exact inverse embedding for every
transported coefficient.

This is a coefficient-field statement. It supplies no degree-four
finite map of surfaces or degree-four target subfield.

## The sign symmetry preserves the actual source problem

On the Wright ring

    A=C[Y,XY,X^2Y,X^3Y+X] inside C[X,Y],

use the source involution iota(X,Y)=(-X,Y). It fixes Y and X^2Y and
negates XY and X^3Y+X, so it preserves A. In the omitted chart
X=1/u, Y=u^3*v-u^2 it is (u,v) -> (-u,-v). Its action on the four
actual chart generators is checked directly.

The complete Cox basis has chart expression

    f_(m,j,i)=(-v)^(m-j)*(1-v^2-u*v)^j*(v+u)^(m+j-i).

Each factor gives the exact identity

    f_(m,j,i)(-u,-v)=(-1)^i*f_(m,j,i)(u,v).

Thus each full pole space V_m is preserved. This identity is valid for
all allowed m,j,i; the checker additionally verifies every whole basis
polynomial at m=3,6,9. It retains all lower-pole terms.

The marking changes x,y,s,delta to their negatives. In particular
H_(-y)(w)=-H_y(-w), with H monic of degree nine. On the original
source, precomposing P,Q by iota negates their constant Jacobian and
preserves its nonvanishing. It preserves finiteness of A over C[P,Q]
and the field-extension degree because it is an algebra automorphism.
The clean cusp, actual-root condition and all later Jacobian equations
are carried to the corresponding conjugate marking. This transports
the obligations; it does not prove they have solutions.

## Exact transport of all thirty-nine quadratics

Keep the [PF2RULED29 presentation](DEGREE_EIGHT_QUADRATIC_PRESENTATION.md),
including its complete root family a0,...,a3, perturbation coefficients
c2,...,c7 and auxiliary polynomials B,M,N. Set v=y*w and use the
following substitutions from new variables to old ones:

| Old variable | Substitution |
| --- | --- |
| b_j, j=0,...,8 | y^(6-j) times new b_j |
| m_j, j=0,...,5 | y^(17-j) times new m_j |
| n_j, j=0,...,10 | y^(22-j) times new n_j |
| c_j, j=2,...,7 | y^(8-j) times new c_j |
| a0,a1,a2,a3 | y times new a0, new a1, y times new a2, new a3 |
| kappa, invkappa | y^29 times new kappa, y^(-29) times new invkappa |

For the coefficient of v^j in the first polynomial identity, divide
by y^(17-j). For the second identity divide by y^(29-j). Leave the
inverse equation unchanged. Every resulting coefficient is fixed by
sigma and lies in F. The checker transports every sparse coefficient
and then applies the inverse substitution to recover each original
polynomial exactly.

The result still has 38 variables and 39 equations of degree at most
two over its coefficient field. The equation kappa*invkappa-1 is
unchanged. All powers of y used here are nonzero field constants;
no parameter polynomial is inverted and no rank-drop, coefficient-zero
or critical-root collision case is lost. The prior monic-elimination
proof therefore remains valid after this change and field extension.

The [retained certificate](degree8_quartic_descent_certificates.json.gz)
contains all quadratic polynomials, both quartic defining polynomials,
every variable and equation scaling, and the compact descended core.

## The four-parameter core and unit rows descend as well

Starting independently from the primary PF2RULED27 data, set

    q4=y^4*q4_new, q5=y^3*q5_new,
    q6=y^2*q6_new, q7=y*q7_new, Z=y^12*z,
    s_new=s/y,
    W_new(w)=y^(-12)*W(y*w),
    M_new(w)=y^(-17)*M(y*w),
    L_new(w)=y^(-6)*L(y*w).

All coefficients lie in F. The parameters remain free; the new s even
generates the entire quartic subfield, so a rational-coefficient claim
is not made. With h_new=w-s_new,

    D_new=(M_new(s_new)-M_new)/h_new = y^(-16)*D(y*w),
    E_new=rem_w(M_new*W_new/h_new,L_new) = y^(-28)*E(y*w).

The equalities follow from exact division and the nonzero leading
scalar 8 of L_new. W_new vanishes doubly at s_new and the degrees in
w are still 12, at most 5, and 6 for W_new,M_new,L_new.
The entire necessary system is therefore

    E_new[5]=0,
    D_new[i]*z-E_new[i]=0, i=0,...,4,
    M_new(s_new)*z != 0.

The normalized kappa is y^(-29) times the old one. Hence its nonzero
condition is equivalent over every complex embedding.

The original exact multipliers C_i from
[PF2RULED28](DEGREE_EIGHT_GLOBAL_CONSTANT_RECOVERY.md) transport by

    C_new[i]=y^(16-i)*C_i(y^4*q4_new,...,y*q7_new).

Every transported multiplier coefficient is already fixed by sigma;
no averaging or new search is required. The checker verifies the full
characteristic-zero polynomial identity

    sum_(i=0)^4 C_new[i]*D_new[i]=1

in F[q4_new,q5_new,q6_new,q7_new]. Their degrees remain at most three.
Thus the all-zero-D stratum stays excluded, and the global recovery
z=sum C_new[i]*E_new[i] gives the same necessary ideal after eliminating
the root constant and preserving saturation. Keep E and the recovered
constant as exact circuits; no dense degree-fifteen expansion is needed.

## Why all sixteen embeddings are covered

Each of the four complex embeddings of F extends to the two choices
of y with y^2=t. For either choice, the displayed diagonal changes
are invertible. They identify the scalar extension of the descended
ideal with the original ideal at that octic embedding. Consequently
emptiness, nonemptiness and any exact unit certificate over F apply
to both sign choices. Across the two sectors it suffices to handle
eight quartic embedding problems instead of sixteen octic ones.
No embedding or parameter locus is discarded.

For a coefficient search with N field unknowns, restriction of scalars
now requires 4N rational unknowns instead of 8N. This is an exact
size reduction, not a measured runtime or a claim that a full-core
certificate has been found.

## Remaining certificate gate and replay

Seek an exact F-coefficient unit identity in either the descended
quadratic ideal with its inverse equation, or the descended saturated
four-parameter core. An appropriate power-of-kappa certificate in
the unsaturated ideal also suffices. Each sector must be covered.
The prior finite-field unit observations remain separately labelled
as finite-field evidence.

A retained [verbose native diagnostic](sources/degree8-quartic-descent/native_unsaturated_diagnostic.log)
replayed the existing first-sector unsaturated residue input. It reached
[1] in about 24.6 seconds, with a largest matrix of 10,639 by 10,493.
This describes the modular algorithm only; it supplies no exact unit
multipliers or characteristic-zero inference. Its process is terminal.
The sign descent and descended unit rows use exact characteristic-zero
arithmetic and do not depend on that diagnostic.

Run from the repository root:

```sh
sage -python research/jacobian/affine-complements/verify_degree8_quartic_descent.py
.venv/bin/python research/jacobian/affine-complements/audit_pole_integration.py
```

The [checker](verify_degree8_quartic_descent.py) pins its predecessor,
checks every coefficient transport and inverse, verifies the exact
unit rows and source parity identities, and writes the
[receipt](degree8_quartic_descent_receipt.json). It performs no full-core
Groebner calculation or expanded E calculation.

The actual-global-root hypothesis, unsolved nonzero-kappa core, general
root-free determinant and rank-drop locus, later normal coefficients
through u^9, contact cases, other height-nine pairs, higher heights,
all-height pole-four, field degree four, global surface finiteness and
actual compatible embedding remain in the original goal. No source
Keller pair is constructed, no whole marked candidate is excluded,
and the goal remains active. There is no formal verification,
independent whole-proof implementation, external review, new imported
theorem or novelty-priority claim.
