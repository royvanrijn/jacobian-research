# Intrinsic Nambu directions and the inverse-Hessian quadratic cone

15 November 2026. For a four-variable unit-Hessian potential, three
independent locally finite Nambu directions force a polynomial gradient
inverse. A single locally finite direction also suffices if it has
nonzero inverse-Hessian square. Under hypothetical nonautomorphy the
entire locally finite direction space is exactly the annihilator of
Tang's unique affine plane. Its directions all have inverse-Hessian
square zero. This gives intrinsic criteria for arbitrary actual HC4
potentials, without choosing ambient observables, principal generic
components, graph coordinates or a degree bound.

The quadratic cone of constant covectors `c` satisfying
`c^T (Hess h)^(-1)c=0` must, under nonautomorphy, be exactly that
linear two-plane. Every such direction produces an explicit polynomial
quasi-translation. A cone of any other shape forces the inverse.
A linear two-plane can also occur for an automorphic gradient, so the
necessary condition is not a converse or a solution of HC4/JC2.

Claims HC4MYIT0/1/2 are recorded in
[MATH_STATUS.json](../../MATH_STATUS.json). The complete admitted
overlap closure and its residual-free JC2 boundary remain unchanged.
Neither plane endpoint nor unrestricted HC4/JC2 is resolved.

## 1. A locally finite Keller direction is already locally nilpotent (HC4MYIT0)

Let `F=(F_1,...,F_n)` be a complex Keller map, `R=C[x]`,
`A=C[F]`, `k=Frac(A)` and `L=Frac(R)`. The extension `L/k`
is finite separable. Define polynomial derivations

```
delta_i(F_j)=delta_ij,
delta_i=sum_j (JF^(-1))_(j,i) partial_(x_j),
D_c=sum_i c_i delta_i, c in C^n.                    (1)
```

They commute: their commutators vanish on `k`, and a derivation of
a finite separable extension which vanishes on its base field is zero.
This also proves uniqueness of their field extensions without any
assumption that `F` is invertible.

Call a derivation locally finite if every polynomial has a
finite-dimensional complex span of its iterated derivatives. If
`D_c` is locally finite, it is locally nilpotent.

For `c=0` this is immediate. Otherwise suppose a nonzero polynomial
`a` is an eigenvector, `D_c(a)=lambda a`, with `lambda in C*`.
Its field norm `r=N_(L/k)(a)` is a nonzero rational target function.
Uniqueness of separable extension of derivations makes them commute
with all conjugate embeddings, so the product defining the norm gives

```
D_c(r)=[L:k] lambda r.                              (2)
```

Choose constant linear target coordinates in which `D_c=partial_t`.
For a nonzero rational function of `t` over the other target-coordinate
field, its logarithmic derivative is `O(1/t)` at infinity. It cannot
be a nonzero constant. This contradicts (2). Thus `D_c` has no
nonzero eigenvalue on any finite-dimensional invariant subspace.
Over `C` its restriction to each such subspace is nilpotent. Local
finiteness therefore implies local nilpotence on all of `R`.

The locally finite directions form a complex linear subspace: the
`D_c` commute, and a linear combination of commuting locally nilpotent
derivations is locally nilpotent. For a particular polynomial killed
by `D^m` and `E^n`, the binomial expansion shows it is killed by
`(D+E)^(m+n-1)`. No common bound over the whole ring is claimed.

### A single nonzero direction suffices in the plane

For a plane Keller pair, if `D_c` is locally finite for some
`c!=0`, the pair is a polynomial automorphism. Make a constant
linear target change to `(U,V)` with `D_c U=1,D_c V=0`.
The preceding argument makes `D_c` locally nilpotent.

Its field of constants is exactly `C(V)`. Indeed any other constant
transcendental over `C(V)` would make `L` algebraic over a field
killed by `D_c`, forcing `D_c=0` and contradicting `D_c U=1`.
Thus every constant is algebraic over `C(V)`. HC4MYHX0 applies
because the polynomial Keller component `V` has no critical point,
and says `C(V)` is relatively algebraically closed in `L`.
Intersecting with `C[x1,x2]` gives precisely `C[V]`: if a reduced
fraction `p(V)/q(V)` is polynomial, a Bezout identity for `p,q`
makes `q(V)` a polynomial unit, hence `q` constant.

The elementary locally nilpotent slice argument now gives

```
C[x1,x2]=(ker D_c)[U]=C[V,U].                       (3)
```

Here is the needed slice argument rather than an imported
classification. If `D(s)=1`, then for each `a` the finite sum
`pi(a)=sum_n (-s)^n D^n(a)/n!` is killed by `D`, by cancellation
of adjacent terms. Equivalently, induction on the largest nonzero
derivative subtracts `(D^m(a)/m!)s^m` and expresses `a` as a
polynomial in `s` with constant coefficients for `D`. Applying the
highest derivative proves uniqueness. This proves (3) and the inverse,
without Rentschler, cancellation or a plane degree classification.
The one-direction inverse conclusion here is plane-only; it is not
asserted for arbitrary `n`.

## 2. One actual Legendre energy tests a whole direction (HC4MYIT1)

For a Keller potential `h`, put `F=grad h`, `H=Hess h` and

```
G=x dot grad h-h.
```

Since `grad G=H x`, direct multiplication by `H^(-1)` gives

```
delta_i(G)=x_i,
D_c^2(G)=c^T H^(-1)c=:q_h(c).                       (4)
```

Primary recent context is [Adamus--Hajto,
arXiv:2610.36736v1, Section 2 and Lemma 4.1](https://arxiv.org/html/2610.36736v1).
The identities and the assertions used below are proved here. We do not
import their differential-Galois criterion, inverse-degree bound,
Pascal-finiteness claims or an assertion that the open cases are settled.

For each constant `c`, the following conditions are equivalent:

```
span_C{D_c^n G:n>=0} is finite-dimensional;
D_c is locally finite on R;
D_c is locally nilpotent on R;
D_c^m G=0 for some m.                               (5)
```

Only the first implication needs explanation beyond Section 1.
A polynomial `p` annihilates the finite orbit of `G`. Since the
Nambu derivations commute,

```
p(D_c)x_i=delta_i(p(D_c)G)=0.
```

Thus every source generator has a finite-dimensional invariant span.
Leibniz makes the span of products of these spaces of any fixed
bounded product degree finite-dimensional and invariant, proving local
finiteness on every polynomial. The other implications follow from
Section 1 or directly. No finite derivative window is accepted unless
an exact recurrence or termination is proved.

Write `Lambda_h` for the linear space of directions satisfying (5).
If the gradient is an automorphism, all Nambu directions are transported
constant target derivatives, so `Lambda_h=C^n`.

### In HC4 the only dimensions are two and four

Suppose `h` has four variables and nonautomorphic gradient. The
retained reproduced theorem HCTANG2, together with HC4MYHF0, supplies
the unique constant jointly affine two-plane `W`. In constant affine
source coordinates write

```
h=D(z)+w^T K(z),  J(K) in C*,
eta=grad D+DK^T w,
grad h=(eta,K).                                    (6)
```

The change `(z,w)->(z,eta)` is a polynomial source automorphism;
its inverse uses the polynomial matrix `DK^(-T)`. We use it to
write actual derivations, not to claim Hessian invariance under a
nonlinear source change. In these coordinates the Nambu direction
`c=(a,b)` is

```
D_c=a dot partial_eta + nu_b,
nu_b=(DK^(-1)b) dot partial_z, holding eta fixed.     (7)
```

For `b=0`, this is locally nilpotent. If `b!=0` and it were
locally finite, its restriction to the invariant subring `C[z]`
would make `nu_b` locally finite. The plane result (3) would make
`K` invertible, and (6) would then give the gradient inverse by
recovering `z` from `K` and `w` from `eta`. This contradicts the
assumption. Consequently

```
Lambda_h=W^perp,
dim Lambda_h=2 under nonautomorphy;
dim Lambda_h=4 under automorphy.                    (8)
```

Directions are naturally source covectors, with the same coefficients
as constant translation directions of the gradient target. This
identification is preserved by constant linear source changes. Under
a source translation `x=A u+b0`, the energy changes by
`-b0 dot F`, whose second derivative in any `D_c` is zero.
Thus finite orbit and the quadratic identity in (4) transform as
claimed under the affine coordinates used in (6).

In particular **three linearly independent directions with finite
energy orbit force a polynomial gradient inverse**, in all degrees
and for every actual HC4 graph embedding, with no ambient-observable
or generic-class hypothesis.

There is also a one-direction test. In (6),

```
G=z dot eta-D(z).
```

Every direction in the conditional plane (8) has `D_c^2 G=0`.
Therefore a direction satisfying (5) but `q_h(c)!=0` forces the
inverse. More generally two finite-orbit directions with
`c^T H^(-1)d!=0` suffice, since their sum has nonzero square
when needed. These are sufficient certificates, not proofs that
such extra directions always exist.

## 3. A quadratic cone and mixed trace certificates (HC4MYIT2)

Define the algebraic cone of constant covectors

```
Z_h={c in C^4:q_h(c)=0 as a polynomial in x}.        (9)
```

Although `H^(-1)` may have arbitrarily large source degree, the
conditions are homogeneous quadrics in four direction coordinates:
expand it into its finitely many source monomial coefficient matrices
and require `c^T H_alpha c=0` for every coefficient matrix.

Every direction in (9) has a particularly strong certificate.
Commutativity and (4) imply

```
D_c^2 x_i=delta_i(D_c^2 G)=delta_i(q_h(c))=0.
```

Thus `D_c` is locally nilpotent and its exponential is the explicit
polynomial quasi-translation

```
x -> x+t H^(-1)(x)c,                               (10)
```

with inverse obtained by replacing `t` by `-t`. It translates the
gradient target by `tc`, since `D_c F=c`. This proves
`Z_h subset Lambda_h` without Tang or action classification.
Conversely, (8) and the affine energy in (6) show that under
nonautomorphy every direction of `Lambda_h` lies in (9). Hence

```
Z_h=Lambda_h=W^perp
under hypothetical HC4 nonautomorphy.               (11)
```

### A zero-square plane reconstructs the actual coefficient pair

Conversely, any supplied linear two-plane contained in `Z_h`
recovers an actual cotangent form, even without assuming
nonautomorphy or invoking Tang. Use a constant linear source change
so this covector plane consists of the first two coordinate
covectors. Polarization makes the upper-left two-by-two block of
the actual inverse Hessian zero, so write

```
H^(-1)=[[0,A],[A^T,B]].
```

Its determinant is `det(A)^2`. Because `det H` is a nonzero
constant, `det A` is a nonzero constant too; hence `A^(-1)` is
polynomial. Block inversion gives

```
H=[[-A^(-T) B A^(-1), A^(-T)],
                         [A^(-1),0]].              (12)
```

Thus the actual potential is jointly affine in the last two
source variables:

```
h=D(z)+w_1 K_1(z)+w_2 K_2(z),
D=h(z,0), K_i=partial_(w_i)h(z,0), J(K) in C*.
```

This is a literal reconstruction using constant coordinates and
the given polynomial, not an assertion about the Hessian after
a nonlinear change. Therefore the remaining linear-plane case
has a directly recoverable intrinsic Keller pair; it need not
be approached through larger derivative windows or arbitrary
ambient extensions. Its inverse problem is exactly that pair's
inverse problem. When the entire cone is that plane the plane
is unique, while the recovered pair can still be automorphic.

It follows that **unless the exact quadratic cone (9) is a linear
two-plane, the gradient is a polynomial automorphism**. This includes
the zero cone, a line, a nonlinear cone, or a cone whose linear span
has dimension at least three. For a concrete nonlinearity certificate,
`c,d in Z_h` but `c+d notin Z_h` suffice. The whole space
`Lambda_h` is linear, so this supplies a finite-orbit direction with
nonzero inverse-Hessian square. No component count or degree bound
in the source variables is required.

The earlier HC4MYHF1 pure-trace criterion is compatible with the same
intrinsic covector plane. Under nonautomorphy the nonzero covectors
with one pure whole Hodge trace are exactly `W^perp minus {0}`.
Therefore directions certified by any mixture of

- an exact pure whole Hodge trace as in HC4MYHF1;
- an exact finite recurrence or termination of the actual energy orbit;
- the polynomial quadratic identity `q_h(c)=0`;

all lie in the same two-plane under nonautomorphy. If their combined
linear span has dimension at least three, the gradient is invertible.
For example two independent pure-trace covectors and one transverse
finite-energy direction suffice. These tests can be mixed without
assuming they classify all trace directions for automorphic gradients.

### The linear two-plane is a real unresolved case

For

```
h=s^4+t^4+s u+t v,
```

the gradient has the explicit inverse
`(a,b,c,d)->(c,d,a-4c^3,b-4d^3)`, while

```
q_h(a1,a2,b1,b2)=2(a1 b1+a2 b2)-12s^2 b1^2-12t^2 b2^2.
```

Its exact zero cone is just `b1=b2=0`, a linear two-plane.
Thus (11) is a necessary condition for nonautomorphy, not a converse.
The transverse direction `b1=1` has nonzero square `-12s^2`,
but its energy orbit terminates at order five. The finite-energy
criterion closes it while the zero-cone shape alone does not.

For the quadratic `h=su+tv`, the cone is the nonlinear quadric
`a1 b1+a2 b2=0`. The individually zero directions `a1=1`
and `b1=1` have a sum of nonzero square. This is an exact
quasi-translation and nonlinearity control, not a nonautomorphic sample.

Local finiteness must not be inferred from a polynomial slice alone:
`partial_x+y partial_y` has slice `x` but fixes every derivative
of `y` equal to `y`, so is not locally nilpotent. It is locally
finite; Section 1 excludes this eigenvalue specifically for Keller
Nambu derivations by the field norm, not for arbitrary derivations.
The canonical relative derivation from the IS family is not an
automatic third intrinsic certificate either. Explicitly, for
`M=f-T`, `V=M Y+g` and `J=J(f,g)`, put

```
N=((g_2+Y f_2)/J)partial_X1
                  -((g_1+Y f_1)/J)partial_X2.
```

It fixes `T,Y,Z,V` and has `N(M)=1`. If `K` is invertible,
true Keller coordinates make it `partial_f-Y partial_g`, an LND.
Conversely if it is locally nilpotent, extend it to the geometric
generic domain `S` of IS1. The slice proof above gives
`S=(ker N)[M]`, so `S/(M)` is a domain. IS1 identifies this
quotient as the product of `d` polynomial rings. Hence `d=1`,
and the retained birational Keller theorem gives the inverse.
Local nilpotence here is exactly the original inverse obligation;
the mere relative slice does not provide it.

The [source record](sources/intrinsic-nambu-context.json),
[checker](verify_intrinsic_nambu.py),
[receipt](receipt-intrinsic-nambu.json) and
[validation](validation-intrinsic-nambu.json) distinguish the written
general proof from exact controls. The retained Tang classification
and regular-parameter theorem remain explicit inputs. No new
unproved classification, formal verification, independent replay,
external review or novelty priority is claimed. A remaining linear
two-plane cone, without a further finite-energy or trace certificate,
does not close either plane endpoint or unrestricted HC4/JC2.
