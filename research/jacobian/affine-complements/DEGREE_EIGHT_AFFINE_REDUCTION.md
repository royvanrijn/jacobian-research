# A global affine model for the rooted degree-eight gate

14 November 2026. PF2RULED27 replaces the nonlinear high-remainder
presentation in [PF2RULED26](SQUARE_BOUNDARY_DEGREE_GAP.md) by a global
polynomial parametrization. In each marked coefficient field the entire
degree-eight high locus is isomorphic to affine five-space. No parameter
discriminant or rank minor is inverted in this parametrization.

Four parameters determine the nonconstant part of the root and the
boundary perturbation. The fifth is Z=S0(s). The full remaining necessary
divisibility test is six equations linear in Z, with the condition
M(s)*Z!=0 retained. It admits five explicit pivot charts on four core
parameters, together with a separately retained zero-coefficient stratum.
This is an exact reduction of the saved necessary equations, not a
solution or exclusion of their nonzero-kappa locus.

## The original necessary problem

Keep PF2RULED26's actual global-root hypothesis, all normal freedoms and
its exact marked coefficient fields K=Q[y]/(P8). Normalize the leading
coefficient of the boundary difference R to1, as in that proof. Write

    R=h^8+sum_(j=2)^7 c_j*h^j, h=v-s,
    C_R=J1(S,T_R)=S0'*U+M.

The complete four-parameter S family uses coordinates a0,...,a3. The
high necessary equations say deg M<=5. Their saved ring has ten variables

    (a0,a1,a2,a3,c2,c3,c4,c5,c6,c7).

The remaining necessary equations are

    kappa_hat=2*M(s)*S0(s),
    kappa_hat-2*M*S0 divisible by h*L, L=R'/h,           (1)
    kappa_hat!=0.

Here kappa_hat is the original constant divided by the nonzero leading
coefficient of R. No source or target leading symbol is changed by that
algebraic normalization. A solution of (1) still requires the actual
section-lifting and later Jacobian/finiteness obligations identified in
PF2RULED26.

## Five variables can be eliminated without a parameter denominator

Let g0,...,g5 be the saved high Groebner generators, in their retained
order. The two equations g1,g2 are affine linear in a1,a2. Their2-by2
coefficient determinant is a nonzero element of K, independent of all
ten parameters. It is not a generic rational-function determinant.
Consequently these two variables are polynomial expressions in c4,...,c7.

The linear equations g3,g4,g5 then solve a0,c2,c3 polynomially. The
remaining g0 vanishes identically after those substitutions. Thus the
complete high ideal is the polynomial graph

    (a0,a1,a2,c2,c3)=Phi(c4,c5,c6,c7),

with a3 free. Every component of Phi has degree at most3. The inverse
map simply reads the five free coordinates. In particular special values
of c6,c7, colliding critical roots and vanishing later pivot coefficients
are included; no exceptional parameter set has been discarded.

The checker verifies explicit identities expressing each graph generator
as a combination of g0,...,g5, and verifies that every high generator
vanishes identically on the graph. It also replays the saved identities
in both directions between the original high coefficients and g0,...,g5.
Together these establish equality of ideals and a polynomial isomorphism,
not just a birational parameter count.

For reference, the constant elimination determinants reduce modulo1009
to681 and604 in the two sectors, respectively. The proof itself checks
the exact nonzero elements and their inverses in K. It uses no modular
rank inference. Exact identities in both irreducible fields cover all
sixteen complex marked points.

## Replace the free constant by the root value at the cusp

The free coordinate a3 changes S0 by a fixed nonzero scalar mu. It
changes none of S0', M or L. Write the resulting root as

    S0=Z+W(v;q), q=(c4,c5,c6,c7), Z=S0(s),
    W(s)=W'(s)=0.

This is another globally invertible polynomial coordinate change:
a3=(Z-beta(q))/mu. Both mu and its inverse belong to K. The functions
W,M,L have coefficients polynomial in the four q parameters alone.
The polynomial L has degree6 and fixed leading coefficient8.

The saved constant becomes exactly kappa_hat=2*M(s)*Z. Therefore its
nonzero condition is M(s)*Z!=0, including every special q value where
either factor might vanish. A nonzero-kappa solution cannot lie at
M(s)=0 or Z=0; no additional genericity assumption is made.

## Six linear equations replace the expanded remainder

Define, by monic polynomial operations in the actual v coordinate,

    D(v)=(M(s)-M(v))/h, deg D<=4,
    E(v)=remainder_v(M(v)*W(v)/h, L(v)), deg E<=5.       (2)

The divisions by h are exact. Division by L uses its fixed leading
scalar8, so introduces no parameter-dependent denominator. All
coefficients D_i,E_i are polynomials in q over K.

For F=M*W/h=L*T+E, the full numerator in (1) is

    2*M(s)*Z-2*M*(Z+W)=2h*(D*Z-E)-2h*L*T.

Because deg(D*Z-E)<=5, its remainder modulo h*L is exactly

    2h*(D*Z-E).                                       (3)

Thus the full saved divisibility ideal, after the global high-locus
parametrization, is equivalent to the six equations

    E5=0, D_i*Z-E_i=0 for i=0,...,4,                   (4)
    M(s)*Z!=0.

Multiplication by the monic h does not change the coefficient ideal:
its coefficients can be recovered recursively without a division by
any parameter. The checker verifies (3) as an exact polynomial
identity and exact Euclidean remainder in each coefficient field.
The old expanded remainder is not merely sampled at chosen q values.

A single polynomial presentation of the nonzero-kappa necessary locus
is now

    (E5, D0*Z-E0,...,D4*Z-E4, t*M(s)*Z-1)
    in K[c4,c5,c6,c7,Z,t].                             (5)

This keeps the nonzero condition algebraically and has only six variables,
instead of the earlier ten high variables plus a saturation variable.
No emptiness, solution, independence or dimension is claimed for (5).

## Complete pivot charts and the retained zero-coefficient case

If D_i!=0 for some i in0,...,4, solve

    Z=E_i/D_i.

On this chart (4) is equivalent to

    E5=0,
    D_i*E_j-E_i*D_j=0 for j!=i, j=0,...,4,
    D_i!=0, M(s)*E_i!=0.                              (6)

These are five basic compatibility polynomials on four core parameters,
with the stated inequalities. They are not asserted independent or
inconsistent. Five possibly overlapping charts cover the whole D!=0
part. If a disjoint stratification is desired, additionally set
D0,...,D_(i-1)=0 on the i-th chart. The certificate records this option;
it does not discard those zero coefficients.

If every D_i is zero, then M is constant. The complete remaining case is

    D0=...=D4=0, E0=...=E5=0, M(s)!=0,
    Z arbitrary and nonzero.                          (7)

Equivalently, since M=M(s) is a unit on this case, L must divide W/h.
This stratum has not been excluded by the present proof. In particular
one may not divide by a D coefficient and silently lose (7).

An exact bounded exploration in sector(2,7) reported a unit Groebner
basis for the constant-M equations, but extraction of an explicit
identity did not finish within its time cap. The second sector was
not certified by that attempt. The [input](sources/degree8_constant_remainder_unit_attempt.py.txt)
is preserved. That observation is not used to remove (7), and neither
certificate extraction nor a full elimination is run by this checker.

## Compact reproducible input and remaining scope

The [compressed certificate](degree8_affine_reduction_certificates.json.gz)
retains the exact global graph, its ideal identities, the scalar constant
change, W,M,L,D and the arithmetic circuit (2) defining E. It records
hashes and sizes of the expanded E coefficients for replay, without
maintaining a second multi-megabyte expanded remainder. These circuit
inputs are actual retained polynomials and are sufficient to reconstruct
the six equations exactly.

The [checker](verify_degree8_affine_reduction.py) and
[receipt](degree8_affine_reduction_receipt.json) replay the field and
polynomial identities, old-input seals, the complete parametrization,
(3), and every chart identity. The ideal-isomorphism and chart-coverage
arguments are written proofs using those exact identities. Shared exact
utilities and the same number-field conventions as the predecessor
are used. No formal verification, independent whole-proof implementation,
external review, new imported theorem or novelty priority is claimed.

Run from the repository root:

```sh
sage -python research/jacobian/affine-complements/verify_degree8_affine_reduction.py
.venv/bin/python research/jacobian/affine-complements/audit_pole_integration.py
```

The nonzero-kappa solutions of (5), all chart cases (6) and (7), actual
section lifting, the general determinant locus without the actual-root
hypothesis, later normal coefficients u^3 throughu^9, contact cases,
other height-nine pairs, higher heights, all-height pole-four, field
degree4, surface finiteness and actual compatible embedding remain open.
The first open pairs remain (5,9),(6,9),(7,9),(8,9). This reduction
advances the original full goal and leaves it active.
