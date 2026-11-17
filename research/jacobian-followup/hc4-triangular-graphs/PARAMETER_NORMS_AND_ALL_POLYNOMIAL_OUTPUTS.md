# Parameter norms and every polynomial triangular output

15 November 2026. This note closes the entire plane output class

```
Pi=(T,xy-P(T,z)), P in C[T,z], deg_z(P)>0.     (1)
```

For every complex plane Keller pair `(f,g)`, the polynomial
`P(f,Z)+g` is irreducible over `C(s,t)`. Consequently every carried
closed stable embedding through (1), with `r>=1`, forces the Keller
inverse by the retained global filtered theorem HC4MYHZ2. Every
source, embedding, coefficient, leading-multiplicity and `z` degree is
unrestricted. Full splitting-group preservation is not asserted.

The proof combines four unrestricted statements, recorded as
`HC4MYIK0/1/2/3` in [MATH_STATUS.json](../../MATH_STATUS.json):
normal monogenic covers after algebraic parameter base change;
integral algebraic affine monic models; integrality of centered
monic coefficients from integral critical values; and a plane
branch-inertia block argument. The first three work in all Keller
dimensions. The last, and the unrestricted conclusion above, use
the plane nonproperness theorem and remain plane-only.

This does not resolve either plane endpoint or unrestricted HC4/JC2.
It retires the remaining plane coefficient-family calculations for
(1). Other output geometries and the actual endpoint source orbits
are separate questions.

## 1. Normal monogenic covers over an algebraic parameter curve (HC4MYIK0)

Let `K:A^n->A^n` be complex Keller, `n>=2`, with components
`f_1,...,f_n`. Write

```
k=C(T), A=C[T,V_2,...,V_n], F=Frac(A),
B=C[x_1,...,x_n], L=Frac(B), T -> f_1.
```

Let `k'/k` be finite Galois, let `D=Spec A_D` be the finite normal
curve over `A1_T` in `k'`, and put

```
A'=A_D[V_2,...,V_n], B'=B tensor_C[T] A_D,
F'=k'(V_2,...,V_n), L'=L k'.                  (2)
```

Suppose `R in A'[Y]` is monic, irreducible over `F'`, and
`Spec(A'[Y]/(R))` is normal. Then its source pullback is irreducible
over `L'`. All source, cover, parameter and polynomial degrees are
unrestricted. The parameter curve may be ramified over `A1`, need
not be an affine line and may have nonconstant units.

### 1.1 The actual source ring after base change

The parameter `f_1:A^n->A1` is smooth: its gradient never vanishes,
as follows from the invertible Keller differential. The normal
complex curve `D` is nonsingular. Smoothness is preserved by base
change and composition, so `Spec B'` is smooth over `C`, hence normal.
The standard structural inputs are explicitly recorded in
[the source audit](sources/parameter-norm-import.json).

The module `A_D` is finite torsion-free over the PID `C[T]`, so is
finite free. Thus `B'` is finite free over `B`. The retained
[HC4MYHX0](REGULAR_PARAMETER_MONODROMY.md#1-nowhere-critical-polynomial-parameters-are-regular-hc4myhx0)
makes `L` linearly disjoint from `k'/k`. Therefore `B'` embeds in
the field `L'=L k'`, is a domain and has that fraction field. The
map `Spec B'->Spec A'` is the base change of `K` and is etale.
The normality input Stacks `025P` makes the root cover's pullback
normal as well.

If `R` factored into two proper monic factors over `L'`, normality
of `B'` would put their coefficients in `B'`. Their root components
are disjoint in the normal pullback. Their degrees cannot drop,
so their cross-resultant is a unit

```
r in (B')*.                                  (3)
```

This is precisely the retained HC4MYHY0 component argument, now
over the actual smooth source `Spec B'`.

### 1.2 The norm replaces the false constant-unit assertion

Multiplication by `r` is an invertible endomorphism of the finite
free `B`-module `B'`. Its determinant is therefore a unit of `B`,
hence a nonzero complex constant. On fraction fields this gives

```
Norm_(L'/L)(r)=product_(sigma in Gal(k'/k)) sigma(r)
             =c in C*.                       (4)
```

The determinant construction and its compatibility with field
extension are elementary here; Stacks `0BD2` is retained as a
primary reference for the norm convention. Equation (4) does not
make `r` constant.

To locate every conjugate in a common target field, take the
compositum `M` over `F'` of the splitting fields of all parameter
conjugates `R^sigma`. This is finite Galois over `F` as well as
over `F'`: the finite parameter Galois group permutes these
polynomials. The cross-resultant (3), a product of root differences,
belongs to `M intersect L'`. Every `sigma in Gal(L'/L)` extends
to the algebraic closure fixing `L`; it stabilizes `M` because
`M/F` is Galois. Its `sigma(r)` is a product of differences of
roots of the monic polynomial `R^sigma`.

All those roots are integral over `A'`. In particular every
`sigma(r)` is integral over `A'`. Equation (4) writes

```
1/r=c^(-1) product_(sigma!=1) sigma(r),        (5)
```

so `1/r` is integral too. Hence `r` is a unit in the finite
normalization of `A'` in `M`. Every `Gal(M/F')` conjugate of `r`
has the same property. Notice the role of all parameter-conjugate
monic models; integrality of just one conjugate would not justify (5).

### 1.3 All finite inertia fixes every conjugate root partition

Let `G=Gal(M/F')`, let `H` be the source image of `Gal(ML'/L')`
in `G`, and let `S` be the proper subset of the roots of `R`
selected by the presumed source factor. Then `H` fixes `S`.
Let `G_S` be its set stabilizer and

```
N=intersection_(gamma in G) gamma G_S gamma^(-1). (6)
```

At every height-one prime of the target normalization, all roots
of `R` are integral and all conjugate cross-resultants are units.
An inertia element taking a root in `S` outside `S` would make
one cross-difference have positive valuation, while all others
have nonnegative valuation. The resultant could not be a unit.
Applying this to every conjugate subset shows that all actual
finite divisorial inertia lies in `N`.

The retained Stacks `0BSD` inertia description and `0BMB` purity
therefore make the finite quotient normalization in `M^N` etale
over `D times A^(n-1)`. Restrict this algebraic cover to the
section with all `V_i=0`. Its connected restriction is an algebraic
finite etale cover of `D`; its pullback has the same topological
monodromy because `C^(n-1)` is simply connected. The retained
full faithfulness of Milne Theorem 3.4 makes the topological
isomorphism algebraic. Both covers already exist algebraically;
essential surjectivity is not used. Thus

```
M^N=E(V_2,...,V_n), E/k' finite Galois.        (7)
```

HC4MYHX0 still preserves this parameter cover after passing to
`L'`: take a finite Galois closure of `E/k` containing `k'`, and
use its linear disjointness from `L`. It follows that the source
image `H` surjects onto `G/N`. Hence `HN=G`. But both `H` and
`N` lie in `G_S`, so `G_S=G`. This contradicts transitivity on
the roots of the irreducible `R`. It proves the theorem.

No full splitting-group equality is inferred from this argument.
It proves irreducibility, even after the specified algebraic
parameter base change.

### A control that forbids the earlier shortcut

Take the nowhere-critical polynomial `f=t(2+s^2 t)` and the
smooth parameter-cover source `u^2=f`. Then

```
(1+s^2 t+s u)(1+s^2 t-s u)=1,
(1-s^2 t) f_t/2+s^3 f_s/2=1.                (8)
```

The first factor is a nonconstant unit of norm one. The second
identity certifies smoothness of `f`. No second Keller component
for this example is asserted. The control shows exactly why a
parameter-cover source cannot be replaced by `A^n`, or its units
declared constant; the proof uses their norms instead.

## 2. Integral algebraic affine monic models (HC4MYIK1)

Let `P in C(T)[z]` have positive degree. Suppose that after a
finite algebraic parameter extension there are `c!=0,beta` such that

```
Q(T,Y)=P(T,(Y-beta)/c)
```

is monic and every coefficient is integral over `C[T]`. Then
`P(f_1,Z)+f_2` is irreducible for every complex Keller map in
every `n>=2`. Enlarge the parameter field to a finite Galois
extension and use its curve `D` above. Now

```
R(T,V_2,Y)=Q(T,Y)+V_2
```

is monic and irreducible over `F'`, since it is linear in the
independent `V_2`. Its total space is `D times A^(n-1)`, after
eliminating `V_2`, hence is smooth and normal. Section 1 proves
irreducibility over `L'`, and therefore over `L`. The invertible
field change `Y=c(f_1)Z+beta(f_1)` gives the stated conclusion.
For polynomial `P` in the plane, HC4MYHZ2 then gives the carried
inverse in (1), with its unchanged global filtration hypotheses.

An explicit sufficient test without centering is the following.
Write `P=a z^m+sum_(i=0)^(m-1) b_i z^i`, with coefficients in
`C(T)` and `a!=0`. If

```
b_i(T)^m/a(T)^i belongs to C[T] for every i<m, (9)
```

take `c=a^(1/m)`, `beta=0`. The monic coefficients
`b_i/c^i` satisfy monic equations `X^m-b_i^m/a^i=0`, so are
integral. A rational input translation may be made first when
it improves (9). This is a sufficient certificate, not a claim
that every polynomial already meets it.

### Every binomial is now covered in the plane

For `P=a z^m+b z^q+d`, `a,b` nonzero polynomials and `1<=q<m`,
the no-finite-pole case `a^q divides b^m` satisfies (9), so is
closed in all Keller dimensions by this section. The complementary
case is closed in the plane by HC4MYIJ3. Thus every such binomial
is irreducible after plane Keller pullback, with arbitrary
coefficients, common zeros, multiplicities and exponent gaps.
The combined conclusion is irreducibility; the finite-pole branch's
stronger full-group theorem is retained in its own scope.

## 3. Integral critical values force an integral centered model (HC4MYIK2)

The following valuation lemma is direct and works in every degree.
Let `v` be a nontrivial nonarchimedean valuation with residue
characteristic zero, extended to an algebraic closure. If a monic
centered polynomial

```
Q(Y)=Y^m+sum_(i=0)^(m-2) a_i Y^i, m>=2,     (10)
```

has `v(Q(alpha))>=0` at every critical point `Q'(alpha)=0`, then
every `a_i` has nonnegative valuation.

Suppose otherwise. Put

```
mu=min_(a_i!=0) v(a_i)/(m-i)<0,
choose c with v(c)=mu,
Q_0(Z)=Q(cZ)/c^m.                           (11)
```

All coefficients of `Q_0` are integral, at least one lower
coefficient is a unit, and it is still monic and centered. The
derivative has unit leading coefficient `m`, so all its critical
points are integral. Their critical values have positive valuation,
by the assumed nonnegative valuations for `Q` and `mu<0`.
After a finite extension splitting the derivative, reduction of its
factorization shows that every critical point of the reduction
`q_0` has value zero. No derivative root disappears because the
residue characteristic is zero.

Over an algebraically closed characteristic-zero residue field,
a degree-`m` polynomial all of whose critical values are zero is
a pure power. Indeed, if it has `r` distinct zeros with
multiplicities `e_1,...,e_r`, its derivative has `sum(e_i-1)=m-r`
zeros there, counted with multiplicity. If all `m-1` critical
points are there, `r=1`. Monicity and centeredness make
`q_0=Z^m`. This contradicts the surviving lower unit coefficient.
It proves the lemma. The same argument may be performed after
extending the residue field; it does not require completeness.

Now let `P in C(T)[z]`, degree `m>=2`, and assume every finite
critical value, as an algebraic function over `C(T)`, is integral
over `C[T]`. First translate the input to remove its degree-`m-1`
coefficient; then scale the input by an algebraic `m`th root of
the leading coefficient. This gives a monic centered `Q`, with
unchanged critical values, over a finite parameter extension.
At every place above every finite `T` level its critical values
are integral, so the lemma makes all coefficients integral.
They belong to the finite normal parameter curve ring. Section 2
therefore proves source irreducibility in all dimensions.

Degree one is immediate by an affine input change, and requires
no critical-value condition. No general dynamical good-reduction
theorem is imported; (10)--(11) are the full coefficient argument.

## 4. Every polynomial output in the plane (HC4MYIK3)

Let `k=C(T)`, and initially allow any nonconstant `P in k[z]`.
Let `M/k(V)` be the finite splitting field of `P(z)+V`, with
group `G`, and let `H` be its source group under a plane Keller
pair `(f,g)`. It is enough to prove that `H` acts transitively
on the `m` roots. Let `N_bad` be the normal subgroup generated
by actual divisorial inertia over curves in `A2_(T,V)` with no
nonconstant polynomial `A1` parametrization. The retained plane
[HC4MYII0/3](BRANCH_CURVES_AND_RATIONAL_TRIANGULAR_PARAMETERS.md)
gives

```
N_bad subset H.                              (12)
```

Its orbits on the roots form a `G`-invariant block system. If
there is just one orbit, (12) already proves transitivity. Suppose
there are `r>1` orbits. We show that `H` permutes these blocks
transitively as well.

### 4.1 The block quotient is again a polynomial cover

Fix a root `z_0`. Its root field is `k(z_0)` with
`V=-P(z_0)`, and its point stabilizer is `G_0`. The stabilizer
of the `N_bad` orbit containing `z_0` is `N_bad G_0`. Thus

```
k(P(z_0)) subset E=M^(N_bad G_0) subset k(z_0). (13)
```

The explicitly imported ordinary one-variable Luroth theorem
over the field `k` writes `E=k(S(z_0))` for rational `S`.
Here `S` can be chosen polynomial, with

```
P=R composed with S, R,S in k[z], deg R=r.    (14)
```

For clarity the polynomial refinement needs no Ritt classification.
The rational outer function initially has all its poles lying
under the unique pole `z_0=infinity` of `P`. Surjectivity of the
finite projective map defined by `S` makes this exactly one
geometric pole. Its point is `S(infinity)`, which is `k`-rational.
A fractional linear change of `S` over `k` sends it to infinity.
The adjusted `R` then has only its infinity pole, so is polynomial;
the adjusted `S` has only the pole `z_0=infinity`, so is polynomial
too. This proves (14) from the ordinary Luroth input.

The `r` distinct values `S(z_i)` are exactly the roots of
`R(Y)+V` and correspond to the blocks. Their splitting field
is a Galois subfield of `M`. Since `N_bad` fixes every block,
it acts trivially on this outer splitting field.

### 4.2 Every remaining critical value is integral

Every horizontal finite branch curve of `R(Y)+V` is an actual
branch curve of `M`, and its inertia acts nontrivially on the
blocks. It therefore cannot be one of the nonparametrizable curves:
their entire inertia was included in `N_bad`. Thus every such
outer branch curve admits a nonconstant polynomial `A1` map.

Each finite critical value `R(alpha)` of the polynomial `R`
defines one of these horizontal branch curves through
`V=-R(alpha)`. This is genuine ramification even for a multiple
critical point: locally the first nonzero higher derivative gives
ramification index at least two. At the generic `T` point the
curve cover is constant-extended by its algebraic critical-value
field, so that ramification persists. Poles at particular finite
parameter levels do not alter this horizontal divisorial assertion.

Suppose such a curve has a polynomial parametrization
`T=T(x), V=V(x)`. It is horizontal, so `T(x)` is nonconstant.
The parameter `x` is integral over `C[T(x)]`: divide the equation
`T(X)-T(x)=0` by the constant leading coefficient of `T(X)`.
Therefore `V(x)` is integral over `C[T(x)]`. The dominant map
embeds the curve's field into `C(x)` over `C(T)`, so its value
`V`, and every conjugate, is integral over `C[T]`.

All critical values of `R` are consequently integral over `C[T]`.
Section 3 followed by Section 2 proves `R(f,Y)+g` irreducible.
Hence `H` acts transitively on the blocks. Since `N_bad subset H`
is transitive inside each of its blocks, `H` is transitive on the
original roots. This proves

```
P(f,Z)+g irreducible over C(s,t)
for every nonconstant P in C(T)[Z].          (15)
```

Degree one is included; no decomposition step is needed there.
The argument allows reducible critical loci, repeated critical
points and values, imprimitive generic groups, leading poles and
arbitrary leading multiplicities. No generic-group classification,
simple branching assumption or bounded computation enters it.

### 4.3 The global output consequence and its limits

When `P in C[T,z]`, the full
[HC4MYHZ2 proof](LEADING_FIBRE_UNIT_GATE.md) applies to (15).
Its common global filtered basis, torsion-free geometric root
class lattice, every affine height-one target fibre and denominator
saturation conditions are already proved there for every positive
`z`-degree `P`. A carried closed stable embedding through (1)
therefore forces the plane Keller inverse, for every `r>=1`.

This is the whole coefficient class, not just the binomial,
square-free-leading, bounded-discriminant or fixed-polynomial
subclasses. Their earlier proofs and full-group statements remain
valid in their narrower or higher-dimensional scopes. Equation (15)
itself asserts irreducibility, not preservation of every splitting
group or every target finite cover. It does not make a parameter
base-change source an affine plane, and uses no nonlinear source
Hessian invariance.

## 5. Verification and remaining problem

The [checker](verify_parameter_norms.py) and
[receipt](receipt-parameter-norms.json) give exact controls for these
written proofs, separate from any claim of finite exhaustion. The
[source record](sources/parameter-norm-import.json) identifies
the retained purity, normality, parameter-regularity and comparison
inputs, the ordinary Luroth import, and the structural smoothness
and norm references. The
[validation](validation-parameter-norms.json) records scoped checks
separately from unrelated repository failures.

No unrestricted higher-dimensional version of (15) is claimed:
the branch-curve containment in (12) is plane-specific. Other
ambient output geometries, obtaining compatible ambient observables
for arbitrary graph embeddings, actual unit-Hessian representatives
in the two endpoint source orbits, both plane endpoints and
unrestricted HC4/JC2 remain open. Formal verification, independent
replay, external review and novelty priority are not claimed.
