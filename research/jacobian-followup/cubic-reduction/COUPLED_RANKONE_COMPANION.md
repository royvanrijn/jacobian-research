# Coupled rank-six companion at the square-factor parameter

KCUB22 completes the specific next construction gate from KCUB21. Fix
z=-4/3 and allow every w, including zero. This family has cubic-output
rank6, derivative-space rank39 for w!=0 and rank37 for w=0. Its saved
101-cube identity at w=0 and uniform104-cube identity do not improve the
89-variable cubic-linear construction. No global optimality follows.

## Model and determinant trust boundary

Use KCUB3 coordinates (a,b,c,R,x5,x6,S,x8,x9,T,U,x12). Compose its
integral12-variable cubic map K with the four commuting source shears

```
x8 += w*c*R,   x12 += w*R^2,
x9 += z*c*R,   U += z*R*(3b-S),     z=-4/3.
```

All right sides use unchanged coordinates, so the inverse subtracts the
same expressions and the Jacobian is1. Independent guarded sparse
Fraction substitution verifies all12 composed outputs and exact
cancellation of every term above cubic degree. Write K E_w=X+Q_w+C_w.
Let Cstar consist of the first five cubic outputs and a*c*R. The embedding
iota has its first five columns e1,...,e5 and final column w*e10+z*e11.
The complete coefficient identity C_w=iota*Cstar is checked, including
all zero rows. Since z!=0, iota is injective for every w. On columns
(a^2*b,a*R*x8,a*R*x5,a^2*c,b^3,a*c*R), the Cstar coefficient minor is
the constant -288. Therefore the cubic-output rank is6 everywhere.

In19 coordinates (X12,W6,tau), form

```
H_top = tau*Q_w(X) - tau^2*iota*W,
H_bottom = Cstar(X),       H_last=0,
J = H + e19*tau^3.
```

Both J=J0+w*J1 and every coefficient of the formula are reconstructed.
The determinant argument is the standard companion Schur factorization
written in KCUB21: det D(I+lambda*H)=det D(K E_w)(lambda*tau*X)=1.
It imports the pinned independently checked KCUB3 determinant, the
unit triangular source maps and this written block argument. The
19-variable determinant is not expanded and no formal proof is claimed.
The padding is removed by the same target shear used for KCUB5; a
101-variable matrix is not constructed or presented as an improvement.

For a base collision p,q, apply E_w inverse and lift each point u to
(u,-Cstar(u),1). The independent checker verifies the full12-coordinate
base collision and full19-coordinate companion collision at the saved
w=0 identity; complete rational coordinates are in validation.json.

## Uniform derivative ranks

Differentiate all19 outputs of J with respect to all19 variables.
The certificate applies Gaussian row operations over Q[w,w^-1] to the
complete quadratic coefficient matrix. Every pivot is a nonzero rational
monomial in w;39 pivots leave a zero residual. Independent arithmetic
reconstructs every derivative, row label, coefficient, pivot and elimination.
Hence the rank is exactly39 at every nonzero complex w, without omitted
exceptional roots. At w=0 a separate exact rational rank is37.

If J=sum b_i*l_i^3 uses m shared scalar linear forms, every partial
derivative of every output lies in span{l_i^2}. Thus m is at least39
offzero and at least37 at zero. These are bounds for these particular
padded tensors, not for arbitrary Keller maps or cubic-linear dimensions.

## The repeated factor and bounded experiment

The anticipated simplification does occur:
3b^2+3z*b*S-z*S^2=(3b-2S)^2/3 at z=-4/3. The producer groups this
square before polarization, and the independent checker reconstructs
all recorded products of three linear factors and matches the full tensor.
It also adds two involutive rotations of the old89-form list and the
three preceding finite pools. The final pool has478 forms, scalar cube
rank191 and227 relevant monomials. Forty seeded termwise restarts give
101 cubes at w=0. The separate40-second exact basis search makes1551
exchanges and stays at101. This is bounded failure to improve, not an
optimality certificate for the pool or the pencil.

The independent checker expands every cube in the saved101-form point
identity, verifies all19 output rows, coefficient-matrix rank19 and the
isolated pure-tau padding row. It also verifies both parameter coefficients
of the saved uniform104-cube identity. Therefore the least cube count over
the whole fixed-z pencil is in[37,101]; every nonzero member separately
has bounds[39,104]. Neither interval changes the global upper bound89.

## Replay and next gate

```
.venv/bin/python research/jacobian-followup/cubic-reduction/verify_coupled_rankone19.py
```

- [Independent validation](coupled-rankone19/validation.json)
- [Exact derivative certificate](coupled-rankone19/derivative-rank-certificate.json)
- [Saved101-cube identity](coupled-rankone19/best-identity.json)
- [Bounded search checkpoint](coupled-rankone19/search-checkpoint.json)

Do not repeat this unchanged478-form pool. Varying z, changing the
remaining carrier parameters or a genuinely different construction remains
open; these checks do not cover them. A separate higher-priority baseline
gate is the [classical cubic theorem source review](CLASSICAL_CUBIC_SOURCE_GATE.md),
followed by auditing the new quartic theorem. This could close entire dense
degree profiles rather than only another carrier slice. The full goal is
still open: ordinary cubic12/36, cubic-linear89, raw3D13 and normalized
3D[10,48] are separate current positive bounds.
