# Reduced pencils preserve rational-parameter monodromy

The previous [all-factor obstruction](FACTOR_RESULTANT_BARRIER.md) left
two unequal factors and a nonnormal integral order. This note takes a
different all-degree route: certify that a rational target parameter
remains relatively algebraically closed after Keller pullback. It closes
scaled families with arbitrary root groups and arbitrary leading
multiplicities, including cases outside the leading triple gcd.

The registered claims are `HC4MYIG0/1/2` in
[MATH_STATUS.json](../../MATH_STATUS.json). Neither plane endpoint nor
unrestricted HC4/JC2 is resolved.

## 1. Pencil tests for a closed rational parameter (HC4MYIG0)

Work over `C`, let `B=C[x_1,...,x_n]`, `L=Frac(B)`, and let
`q=N/M` be nonconstant, with `N,M` coprime and `M!=0`.

**Reduced-pencil test.** If every finite member

```
N-lambda M,   lambda in C,
```

is nonconstant and square-free, then `C(q)` is relatively algebraically
closed in `L`.

**One exceptional polynomial fibre.** If `q=H` is polynomial and there
is `lambda_0` such that every `H-lambda`, `lambda!=lambda_0`, is
square-free, while the gcd of the multiplicities of the irreducible
factors of `H-lambda_0` is one, the same conclusion holds.

**One exceptional rational fibre with independent numerator and denominator.**
The same one-exception test works for `q=N/M` if `N,M` are algebraically
independent: all other finite pencil members are square-free, and the
factor multiplicity gcd of `N-lambda_0 M` is one.

These are sufficient tests, not characterizations of all closed rational
functions. Critical points, singular reduced fibres and reducible fibres
are permitted. The nonconstant condition in the first test is essential.

We explicitly import Petravchuk--Iena, *On closed rational functions in
several variables*, arXiv:math/0701588v1, Lemma 2(1), equivalence (a)--(c):
a nonconstant rational function generates a relatively algebraically
closed subfield if and only if it has no rational composition of degree
greater than one. The proof below establishes this noncomposition. The
[source manifest](sources/reduced-pencil-import.json) distinguishes this
import from the direct pencil argument and the retained Riemann--Hurwitz
input.

### 1.1 At most one parameter value can have no source divisor

Suppose `q=R(h)`, where `h=H_0/H_1` is reduced and nonconstant and
`R:P1->P1` has degree `d>1`. A point `c` of the `h` line corresponds to
a nonzero linear form `ell_c(H_0,H_1)`, well-defined up to a scalar.
If this is nonconstant, its zero set has a prime divisor. At most one
point `c` can give a nonzero constant: two independent linear forms
constant would make both `H_0,H_1` constant.

Write `R=P/Q` using coprime homogeneous binary forms of degree `d`.
The substitutions `P(H_0,H_1)` and `Q(H_0,H_1)` are coprime in `B`:
a common prime would yield a common projective zero of `P,Q`, since
`H_0,H_1` are not both zero at that prime. Thus they are the reduced
numerator and denominator of `q`, up to one common nonzero scalar.

If `c` is ramified with index `e_c>=2` and `R(c)=lambda` is finite,
the polynomial `N-lambda M` contains
`ell_c(H_0,H_1)^e_c`. Distinct point forms remain pairwise coprime after
substitution. Hence in the reduced-pencil test every such `c` must be
the one possible point whose point form is a nonzero constant.

### 1.2 Ramification forces a forbidden constant pencil member

Let `r>=1` be the number of distinct poles of `R`. Ramification over
infinity contributes `d-r`. Riemann--Hurwitz for `P1->P1` says that
total ramification is `2d-2`, so the contribution over finite values is

```
d+r-2.
```

By Section 1.1 all this contribution must occur at at most one point,
where it is at most `d-1`. Thus `r=1`, and there is one totally ramified
point `c_0` over a finite value `lambda_0`. The homogeneous form
`P-lambda_0 Q` is a nonzero scalar times `ell_(c_0)^d`. Its substitution
therefore makes `N-lambda_0 M` a nonzero constant. This contradicts the
nonconstant hypothesis. No rational composition of degree greater than
one exists, proving the first test by the explicit imported equivalence.

The Riemann--Hurwitz input is the retained Stacks `0C1B`; its projective,
finite, separable hypotheses hold for the rational map `R` in
characteristic zero. No affine etale curve map is asserted finite.

The example `q=1/x_1^d`, `d>1`, explains the constant-member boundary:
`1-lambda x_1^d` is square-free for every nonzero `lambda`, but its
`lambda=0` member is constant. This `q` is composite and is not closed.

### 1.3 The polynomial one-exception test

For polynomial `H=R(h)`, every pole point of `R` must have constant
point form after substitution, since `H` has no source prime pole.
There is at most one such point. A projective change moving it to
infinity makes `h` polynomial and `R` polynomial, so we may write
`H=R(h)` with both polynomial and `deg R=d>1`.

A finite critical value of `R` different from `lambda_0` would make
`H-lambda` non-square-free; a nonconstant polynomial `h` has a
nonconstant fibre polynomial at every finite value. Hence every root
of `R'` lies among the roots of `R-lambda_0`. If that polynomial has
distinct roots `c_i` with multiplicities `e_i`, the roots there account
for `sum(e_i-1)=d-r` roots of `R'`, counting multiplicity. Its degree is
`d-1`, so `r=1`. Consequently

```
R(t)-lambda_0=c(t-c_0)^d,
H-lambda_0=c(h-c_0)^d.
```

Every source factor multiplicity is divisible by `d`, contradicting
their gcd one. This proves the second test. Its last counting step is
elementary polynomial differentiation, with no new classification input.

### 1.4 One exceptional rational fibre and algebraic independence

If some point form `ell_c(H_0,H_1)` were a nonzero constant, both
`H_0,H_1` would be affine polynomials in one of them. Thus both
homogeneous substitutions `P(H_0,H_1),Q(H_0,H_1)`, and hence `N,M`,
would be polynomials in a single element, contradicting their algebraic
independence. So no projective parameter point has a constant point form.

In a nontrivial composition `q=R(h)`, all finite ramification of `R`
must therefore lie over `lambda_0`. If `r_0,r_infty` count its distinct
preimages of `lambda_0,infinity`, Riemann--Hurwitz gives
`(d-r_0)+(d-r_infty)=2d-2`. Hence both counts are one. The homogeneous
form `P-lambda_0 Q` is a scalar times one point form to the `d`th power.
Every factor multiplicity of `N-lambda_0 M` is divisible by `d`, again
contradicting gcd one. This proves the third test.

Independence cannot be dropped. Take
`q=(s^3-6s+4)/s^3=4(-1/2+1/s)^3-3(-1/2+1/s)`.
The only nonreduced finite member is at `lambda=-1`, with factorization
`2(s-1)^2(s+2)` and multiplicity gcd one. The member at `lambda=1` is
the reduced linear polynomial `4-6s`, and all other finite members are
reduced. Nevertheless `q` is composite. Its numerator and denominator
are algebraically dependent, so the third test correctly excludes it.

## 2. Keller-stable rational parameters and arbitrary finite covers (HC4MYIG1)

Let `K:A^n->A^n` be any complex Keller map. Suppose a rational target
parameter `q=N/M`, in reduced form, passes one of the tests in Section 1.
Then `q(K)` passes the same test on the source. Indeed the retained
[HC4MYHW0 divisor pullback proof](ABELIAN_CLASSGROUP_OUTPUT_OBSTRUCTION.md)
shows that nonconstant irreducible target polynomials have nonconstant
square-free pullbacks, and distinct irreducibles have coprime pullbacks.
Thus reduced finite pencil members remain reduced, numerator and
denominator remain coprime, and factor multiplicity gcds are unchanged.
Dominance prevents any nonconstant target member from becoming constant.
It also preserves algebraic independence of numerator and denominator.

Consequently

```
k=C(q(K)) is relatively algebraically closed in L=C(x_1,...,x_n).
```

For every finite Galois extension `E/k` in a common algebraic closure,
`E intersect L=k`. The finite Galois extension `EL/L` has restriction
image in `Gal(E/k)` whose fixed field is this intersection. That image
is the entire group. Therefore

```
[EL:L]=[E:k],     Gal(EL/L)=Gal(E/k).          (1)
```

This preserves the full group, not merely root transitivity. Any finite
separable extension is handled by its finite Galois closure. There is no
claim that arbitrary finite covers of the full multi-parameter target
stay disjoint, nor that every fixed parameter fibre is irreducible.

### 2.1 Useful parameters with unrestricted coefficients

In two target coordinates `U,V`, the reduced-pencil test applies to

```
q=(V+B(U))/A(U),   A!=0, A,B in C[U],         (2)
```

since `V+B(U)-lambda A(U)` is a polynomial coordinate and is coprime
to `A(U)`. No condition on the zeros or multiplicities of `A` is needed.
The same conclusion follows directly on a Keller source: each
`g+B(f)-lambda A(f)` is a nowhere-critical Keller component.

The polynomial one-exception test applies to

```
q=Q(U)(V+B(U)),   Q!=0, Q,B in C[U].         (3)
```

For nonzero `lambda` the member `Q(U)(V+B(U))-lambda` is primitive,
linear in `V` and irreducible, hence square-free. At zero, the factor
`V+B(U)` occurs with multiplicity one and is coprime to `Q(U)`, so
the multiplicity gcd is one. Arbitrary repeated zeros of `Q` are allowed.
For example `U^2 V` has a critical line, yet every finite cover of its
one-variable field keeps its group after Keller pullback.

More generally the first test applies to

```
q=(Q(U)V+C(U))/D(U),
Q square-free and nonzero, D nonzero,
gcd(Q,C,D)=1.                                (4)
```

Each finite numerator has square-free content dividing `Q`. After
removing that content the polynomial linear in `V` is primitive and
irreducible, and shares no factor with the content. It is nonconstant.
This includes `UV`, which has a reduced reducible zero fibre and a
critical point. The polynomial and rational tests are not restricted to
parameters with no critical points.

An important extension combines scaling and division:

```
q=Q(U)(V+B(U))/D(U),
Q,D nonzero and coprime, B polynomial.        (4a)
```

For nonzero `lambda` its finite numerator is primitive and linear in `V`,
since `gcd(Q,QB-lambda D)=1`; it is irreducible and reduced. At zero
the factor `V+B(U)` has multiplicity one. If `D` is constant, the
polynomial one-exception test applies. If `D` is nonconstant, the
numerator varies with `V` while the denominator is a nonconstant function
of `U`, so they are algebraically independent, and Section 1.4 applies.
Thus `(V+B(U))/A(U)` is admitted for every nonzero rational `A in C(U)`,
with `B` polynomial, not just for polynomial or reciprocal-polynomial
`A`.

### 2.2 The filtered class-monodromy criterion allows these rational parameters

Keep all the hypotheses of
[HC4MYHX1](REGULAR_PARAMETER_MONODROMY.md#2-nonabelian-monodromy-from-one-regular-target-parameter-hc4myhx1):
the common free nonnegative filtration, generic graded domain with scalar
degree zero, geometric integrality of every affine height-one target
fibre, and normal geometric generic ambient ring with torsion-free class
group. Its one-parameter certificate may now use a rational target `q`
passing a test in Section 1. If a finite Galois cover of `C(q)`
carries the entire finite class action, then every carried closed stable
hypersurface embedding forces the plane Keller pair invertible.

The proof is the retained HX1 argument with (1): the source and target
absolute Galois groups have the same full image on that finite cover.
The class of a source component, fixed by the source group, is therefore
fixed by the whole target group. Its conjugates sum to a principal
divisor, so torsion-freeness makes each class zero. The retained global
filtered component pencil then gives the inverse. The filtration and
all affine height-one fibre hypotheses still hold over the full
polynomial target; poles of `q` do not authorize dropping those fibres.

## 3. Arbitrary scaled one-variable outputs close (HC4MYIG2)

Let `R in C[z]` be any nonconstant polynomial. Let `A,c` be nonzero
polynomials in `C[T]`, and let `beta,B` be arbitrary polynomials. Then

```
P(T,z)=A(T) R(c(T)z+beta(T))+B(T)             (5)
```

stays irreducible in the equation `P(f_1,Z)+f_2` under every
complex Keller base change in any dimension `n>=2`. Its full splitting
group is preserved. All source, parameter, leading multiplicity,
`R` degree and affine-substitution degrees are unrestricted.

Indeed with `q=(f_2+B(f_1))/A(f_1)`, the equation divided by `A(f_1)`
is

```
R(c(f_1)Z+beta(f_1))+q=0.
```

The polynomial `R(Y)+q` is irreducible over `C(q)`: it is primitive
and linear in the transcendental parameter `q` before passing to the
fraction field. Its finite splitting field keeps its full group by
(1). The source-field affine change in `Z` is invertible. On the
target, adjoining `T` and any other target coordinates is purely
transcendental over `C(q)`, so the same is the full target splitting
group as well.

There is a complementary reciprocal-scaling family. If `Q,c` are
nonzero, `beta,B` are polynomial, and

```
P(T,z)=R(c(T)z+beta(T))/Q(T)+B(T) in C[T,z],  (6)
```

then the same conclusions hold. The displayed membership is required:
the division must yield a polynomial. Multiply the source equation by
`Q(f_1)` and use parameter (3),
`q=Q(f_1)(f_2+B(f_1))`. Again the equation is `R(Y)+q=0` after an
invertible affine change, and (1) applies.

For every plane instance of (5) or (6), all closed stable embeddings
carrying a plane Keller pair through `(T,xy-P(T,z))`, `r>=1`, force the
pair to be an automorphism. Full root transitivity invokes the retained
[HC4MYHZ2 filtered output proof](LEADING_FIBRE_UNIT_GATE.md), whose
filtration and every affine height-one geometric fibre conditions apply
to every polynomial `P` of positive `z` degree. No ambient or actual
source Hessian invariance is inferred from the affine change over the
function field.

Two unbounded examples show what this adds:

- `P=T^(3k)(z^3+z)+B(T)`, every `k>=1`, has a repeated leading level,
  permits the intermediate pole order `k`, and is closed by (5).
  The cubic cover `Y^3+Y+q` has nonsquare discriminant `-4-27q^2`
  and full group `S3`.
- `P=T^k z^m+z^(m-1)+B(T)`, every `k>=1,m>=2`, is (6) with
  `R(Y)=Y^m+Y^(m-1)`, `c=T^k`, `Q=T^(k(m-1))`. It is fully
  irreducible after Keller pullback in all degrees. The earlier
  empty-support and all-factor gates alone left possible unequal
  higher factors when `m>=5`.

### 3.1 Rational scalings and rational affine input

In fact (5) and (6) are contained in the following larger family:

```
A,c in C(T)*, beta in C(T), B in C[T],
P(T,z)=A(T)R(c(T)z+beta(T))+B(T) in C[T,z].    (7)
```

The final polynomiality is required, but `A,c,beta` may be rational.
Write `A=D/Q` in lowest terms. The parameter is precisely (4a),
`q=Q(f_1)(f_2+B(f_1))/D(f_1)`. Its field is regular by Section 2.1,
and `c(f_1)Z+beta(f_1)` is an invertible affine change over the source
field. The proof of full irreducibility, full group preservation and
the carried plane inverse is unchanged.

In particular every sparse polynomial

```
P=a(T)z^m+b(T)z^(m-1)+d(T), a!=0, m>=2,     (8)
```

closes, with arbitrary `a,b,d` and no gcd or multiplicity restriction.
If `b!=0`, use `R(Y)=Y^m+Y^(m-1)`, `c=a/b`,
`A=b^m/a^(m-1)`, `beta=0`, `B=d`. If `b=0`, use `R(Y)=Y^m`,
`A=a,c=1`. This upgrades the earlier sparse empty-support rational-root
gate to full irreducibility in every degree and also covers common
repeated leading zeros. It is not an assertion about arbitrary families
whose rational output translation `B` fails the parameter tests.

## 4. Remaining boundary, failed shortcuts and verification

This is a field-regularity argument, not an assertion that all parameter
fibres are irreducible or that all affine etale maps are finite.
`UV` supplies a reduced reducible fibre. `1/x^d` supplies the forbidden
constant pencil member, and `(x/y)^d` supplies a multiple finite member.
Each test requires its stated exception and independence conditions.

The factorization-surface and Darboux-cofactor routes did not supply an
additional properness theorem. In particular a derivation with a
polynomial slice is not assumed locally nilpotent. The broad closure
here instead comes from the reduced pencil and its rational-parameter
monodromy.

General coefficient variation not reducible to (7), class actions not
certified to descend from one of these regular rational parameters,
two unequal factors outside these closures, other joint repeated-leading
cases, uncovered output geometries, singleton-zero and unit-ideal plane
endpoints and unrestricted HC4/JC2 remain open. No necessary pattern is
asserted realizable by a Keller pair.

Run [the exact replay](verify_reduced_pencils.py):

```sh
timeout 60s .venv/bin/python research/jacobian-followup/hc4-triangular-graphs/verify_reduced_pencils.py
```

The [receipt](receipt-reduced-pencils.json) contains finite controls for
the general written arguments, and the
[validation](validation-reduced-pencils.json) separates scoped checks
from unrelated whole-repository missing-artifact failures. The existing
automorphic actual-unit-Hessian controls for arbitrary `P` are preserved.
Formal verification, independent replay, external review and novelty
priority are not claimed.
