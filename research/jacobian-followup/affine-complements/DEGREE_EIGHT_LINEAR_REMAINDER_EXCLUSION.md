# The rooted degree-eight critical remainder cannot be linear

14 November 2026. PF2RULED31 excludes the entire stratum deg M<=1 of
the complete rooted degree-eight necessary problem, at all sixteen
marked points. The zero polynomial is included. Every surviving point
must satisfy

    deg_w M >= 2, equivalently deg_w D >= 1.

This is a characteristic-zero exclusion. In both quartic coefficient
fields an explicit degree-three polynomial U satisfies

    U*E5 = 1 modulo (M2,M3,M4,M5),                    (1)

where Mj is the coefficient of w^j in the actual critical remainder
and E5 is the actual highest compatibility coefficient. Since E5=0
is necessary, (1) rules out every constant or linear M. The proof does
not require nonzero kappa, a generic pivot, or a sampled parameter point.

The remaining degrees 2,3,4,5 and the full nonzero-kappa core are open.
No whole marked candidate, pole pair or quartic packet is excluded.

## Exact scope and primary equations

Use the [PF2RULED30 quartic descent](DEGREE_EIGHT_QUARTIC_DESCENT.md)
of the complete necessary system. The coefficient field F=Q[t]/P4 is
irreducible of degree four, separately in sectors (2,7) and (4,5).
The four parameters q4,q5,q6,q7 remain unrestricted. Put

    h=w-s,
    D=(M(s)-M)/h,
    E=rem_w(M*W/h,L).

The complete root family and all normal freedoms of
[PF2RULED26](SQUARE_BOUNDARY_DEGREE_GAP.md) are still retained.
W is divisible by h, L has degree six with leading coefficient 8,
and the necessary equations include

    E5=0, D_i*Z-E_i=0 for i=0,...,4.

The actual-global-root hypothesis, monic nonzero boundary-difference
normalization and both marked fields remain explicit. This theorem
does not classify determinant points without that root hypothesis.

For a linear or constant M, the four coefficients M2,...,M5 vanish.
Let I=(M2,M3,M4,M5) in F[q4,q5,q6,q7]. We prove I+(E5)=(1)
without imposing the other compatibility equations or nonzero kappa.

## Thirty-five exact polynomial reduction rules

In each field the retained data contain 35 monic polynomials G_a,
one for every monomial q^a of total degree four:

    G_a=q^a-lower_a, deg(lower_a)<=3.

Every rule comes with four explicit polynomial coefficients R_(a,j)
satisfying the full exact identity

    G_a=sum_(j=2)^5 R_(a,j)*Mj.                      (2)

The checker multiplies and verifies all seventy identities over the
two characteristic-zero fields. It also checks the primary W,M,L,s
coefficientwise against the PF2RULED30 source, so these are rules for
the actual coefficient equations.

Replacing a multiple of q^a by the corresponding multiple of lower_a
strictly decreases total degree and preserves congruence modulo I by
(2). Every polynomial therefore has a computable congruent expression
in the 35 monomials of degree at most three. This congruence statement
is sufficient for the proof. No confluence or unverified quotient
dimension assumption is needed.

Discovery used a Gröbner basis of the four coefficient equations and
a 35-by-35 multiplication matrix. Replay uses only the displayed
polynomial memberships, terminating reductions and exact identities.

## Reducing the actual E5 circuit

Reduce the coefficients of W,M,L by those rules. On the stratum, M
reduces to degree at most one in w. Perform the two univariate divisions
in the original E circuit, reducing coefficient products as needed:

    divide W by h,
    divide M*(W/h) by L and retain its remainder.

The checker verifies that the first division has zero remainder and
that L retains its fixed leading scalar 8. The divisions use only
inversion of nonzero field constants. Coefficient reductions preserve
congruence modulo I, and division by a polynomial with fixed unit
leading coefficient preserves that congruence. Thus the resulting
coefficient e5 satisfies

    e5 congruent to the full primary E5 modulo I.    (3)

The retained e5 has parameter degree at most three. Its equality with
the replayed circuit is checked exactly, with all four parameters
present. The dense degree-fifteen E polynomial need not be expanded.

## A full exact inverse certificate

The [retained certificates](degree8_linear_remainder_certificates.json.gz)
contain U of parameter degree at most three in each field. The checker
multiplies every term of U and e5, reduces the complete product using
(2), and obtains exactly 1. Hence

    U*e5-1 belongs to I.

Together with (3), this proves (1) and I+(E5)=(1). At any purported
linear-M point, I and E5 vanish, so (1) evaluates to 0=1. This excludes
the entire stratum over an algebraic closure of F. The field embeddings
and invertible changes of PF2RULED30 cover all sixteen original complex
markings. No particular M1 or M0 value is assumed.

Since D=(M(s)-M)/h, a nonconstant M of degree d gives D of degree d-1.
The new bound is therefore deg D>=1. In particular, the full constant-D
stratum is now excluded on the necessary compatibility locus, extending
the earlier exclusion of the all-zero-D stratum. This is not an
unconditional unit identity for (D1,...,D4) alone.

## Construction, replay and assurance

The producer forms the exact 35-by-35 multiplication matrix of e5 in
the reduction monomials. Restricting scalars on 1,t,t^2,t^3 gives a
140-by-140 rational system for an inverse. FLINT Dixon arithmetic solves
it and verifies the full rational matrix identity. The final checker
then verifies the complete polynomial inverse identity directly; its
proof does not depend on matrix rank flags or modular unit outputs.

Run from the repository root:

```sh
.venv/bin/python research/jacobian-followup/affine-complements/verify_degree8_linear_remainder.py
sage -python research/jacobian-followup/affine-complements/compute_degree8_linear_remainder.py
.venv/bin/python research/jacobian-followup/affine-complements/audit_pole_integration.py
```

The [checker](verify_degree8_linear_remainder.py) uses rational sparse
polynomial arithmetic modulo P4, implemented with python-flint. It
replays the primary circuit and every reduction membership and inverse
identity. The [producer](compute_degree8_linear_remainder.py) invokes
its [rational Dixon stage](compute_degree8_linear_remainder_dixon.py)
through the repository virtual environment. Temporary dense matrices
are removed after use. The [receipt](degree8_linear_remainder_receipt.json)
pins all proof inputs and records the exact assurance scope.

Earlier weighted Singular reduction, experimental Macaulay2 change
matrices, direct number-field inversion and direct Singular product
reduction were unsuccessful within their caps or failed their small
compatibility checks. The [retained attempts](sources/degree8-linear-remainder/attempts.json)
and [inputs](sources/degree8-linear-remainder/degree8_Mlinear_matrix.py.txt)
record the change to finite reduction, exact rational linear solving
and direct FLINT polynomial verification. No such unsuccessful attempt
is a live computation or proof of an exclusion.

The original root-free determinant and rank-drop problem, higher normal
equations through u^9, contact cases, other height-nine pairs, higher
heights, all-height pole-four, field degree four, global finiteness and
actual compatible embedding remain open. Degrees M=2 through 5 remain
in the rooted degree-eight core. No global source Keller pair or finite
map is constructed. Written ideal and congruence arguments are supported
by exact characteristic-zero identities, without formal verification,
independent whole-proof implementation, external review, new imported
theorem or novelty priority. The full research goal remains active.
