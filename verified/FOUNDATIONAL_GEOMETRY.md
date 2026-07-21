# Foundational geometry: C01--C04

This note isolates the complete foundational proof chain.  It uses no claim
after C04.  The four inputs are elementary polynomial identities, two explicit
root charts, normalization plus Zariski's Main Theorem, and a standard
algebraic inertia argument for symmetric monodromy.

## Theorem A: the original map and collision

Put `u=1+xy` and

\[
F(x,y,z)=\left(
u^3z+y^2u(4+3xy),
y+3xu^2z+3xy^2(4+3xy),
2x-3x^2y-x^3z
\right).
\]

Then

\[
\det DF=-2                                                   \tag{1}
\]

and

\[
F(0,0,-1/4)=F(1,-3/2,13/2)=F(-1,3/2,13/2)=(-1/4,0,0).       \tag{2}
\]

### Proof

On the dense open `x ne 0`, set

\[
t=y+1/x,\qquad r=2/x,\qquad c=F_3(x,y,z).
\]

Direct substitution gives

\[
a=t^2+{rt\over2}-ct^3,\qquad b=r+4t-3ct^2.                  \tag{3}
\]

The two coordinate Jacobians are

\[
\det{\partial(a,b,c)\over\partial(t,r,c)}={r\over2},
\qquad
\det{\partial(t,r,c)\over\partial(x,y,z)}=-2x.
\]

Since `r=2/x`, their product is `-2`.  Both sides of (1) are
polynomials, so equality on the dense open proves it everywhere.  Equation
(2) is exact substitution.  The three source points are visibly distinct.
Thus `F` is an everywhere-etale noninjective polynomial map in characteristic
zero.  Appending identity coordinates proves C01 in every dimension at least
three.

There are three independent certificates for this finite calculation: two
local implementations and Dean Cureton's separately authored pinned Lean 4
formalization.  The Lean scope and attribution are recorded in
[LEAN_C01.md](LEAN_C01.md).

A complementary positive covariance test is
[C01_INVARIANCE_REGRESSION.md](C01_INVARIANCE_REGRESSION.md).  It applies
nonlinear source and target automorphisms, replaces the primitive element,
changes the projective root coordinate, and stabilizes, then verifies that
the canonical normalization-boundary object is transported rather than
changed.

## Theorem B: the cubic marked-root isomorphism

Intrinsically, let

\[
\pi:\mathbb P^1\times\operatorname{Sym}^2(\mathbb P^1)
\longrightarrow\operatorname{Sym}^3(\mathbb P^1),\qquad
(p,\{q,r\})\longmapsto\{p,q,r\},
\]

let `R` be its ramification divisor, and let `H` be a hyperplane tangent but
not osculating to the small diagonal `{3p}`.  With

\[
X=(\mathbb P^1\times\operatorname{Sym}^2(\mathbb P^1))
\setminus(R\cup\pi^{-1}(H)),\qquad
Y=\operatorname{Sym}^3(\mathbb P^1)\setminus H,
\]

one has `Y ~= A^3`.  After a change of coordinate, take `H` to be the
vanishing of the `U^2V` coefficient and normalize that coefficient to `-2`.
Then `pi^(-1)(Y)` is the incidence below, `R` is its repeated-marked-root
locus, and hence `X=I^simp`.  The two-chart proof will show `X ~= A^3` and
identify `pi|X` with `F`.

For a target `(a,b,c)`, let

\[
Q_{a,b,c}(U,V)=cU^3-2U^2V+bUV^2-2aV^3
\]

and let

\[
I=\{((a,b,c),[U:V]):Q_{a,b,c}(U,V)=0\}
\subset\mathbb A^3\times\mathbb P^1.
\]

No fiber of `I -> A^3` is the whole projective line because the coefficient
of `U^2V` is `-2`.  Hence every fiber is an effective Cartier divisor of
degree three and the projection is finite flat of degree three.  Let
`I^simp` be the locus where the marked root has multiplicity one.

The morphism

\[
\Phi:\mathbb A^3\longrightarrow I^{\rm simp},\qquad
(x,y,z)\longmapsto(F(x,y,z),[1+xy:x])                         \tag{4}
\]

is an isomorphism.  Under it, `F` is the projection which forgets the marked
simple root.

### Proof on the chart `V ne 0`

Set `u=U/V` and

\[
P(T)=cT^3-2T^2+bT-2a.
\]

The chart ring is

\[
k[a,b,c,u,1/P'(u)]/(P(u)).                                    \tag{5}
\]

The inverse morphism is

\[
v={P'(u)\over2},\qquad x=1/v,qquad y=u-v,qquad
z=5v^2-3uv-cv^3.                                               \tag{6}
\]

Substitution proves both compositions.  On the source this is exactly the
open `x ne 0`, since `u=y+1/x` and `P'(u)=2/x`.

### Proof on the chart `U ne 0`

Set `s=V/U` and

\[
R(s)=c-2s+bs^2-2as^3,qquad d=1-bs+3as^2=-R'(s)/2.
\]

The chart ring is

\[
k[a,b,c,s,1/d]/(R(s)).                                        \tag{7}
\]

The inverse is

\[
x={s\over d},\qquad y=b-3as,
\]

\[
z=d\left(a-4b^2+(b^3+22ab)s-(30a^2+8ab^2)s^2
+21a^2bs^3-18a^3s^4\right).                                  \tag{8}
\]

Formula (8) is regular at `s=0`; there it gives
`(x,y,z)=(0,b,a-4b^2)`.  Direct substitution again proves both
compositions.  On the overlap, `u=1/s` and `v=d/s`, so (6) and (8) agree.
The two source opens `x ne 0` and `1+xy ne 0` cover `A^3`, and the two root
charts cover `P^1`; therefore the chartwise inverses glue and prove (4).

This is a scheme-theoretic isomorphism, not merely a bijection on algebraic
points.  In particular, a target fiber of `F` is canonically the set of simple
projective roots of its binary cubic.

Under the standard identification `Sym^3(P^1)=P^3`, the chart where the
`U^2V` coefficient is nonzero is the complement of the tangent
nonosculating hyperplane used in C02; normalizing that coefficient to `-2`
gives the target coordinates `(a,b,c)`.  Marking one root identifies `I` with
the restriction of

\[
\mathbb P^1\times\operatorname{Sym}^2(\mathbb P^1)
\longrightarrow\operatorname{Sym}^3(\mathbb P^1).
\]

The repeated-marked-root divisor is precisely its ramification divisor.
Thus `I^simp` is the abstract open `X` in C02, and Theorem B proves both
$X\cong\mathbb A^3$ and the advertised coordinate realization of the
addition map.

## C03 as a direct corollary

Over `C`, a nonzero binary cubic has respectively three, one, or zero simple
roots according as its multiplicity type is `(1,1,1)`, `(2,1)`, or `(3)`.
Its discriminant is `-4Q`, where

\[
Q=27a^2c^2-18abc+16a+b^3c-b^2,
\]

and the triple-root locus is

\[
\Gamma=V(3bc-4,12a-b^2).
\]

Theorem B therefore gives the `3/1/0` fiber table and
$F(\mathbb C^3)=\mathbb C^3\setminus\Gamma$.  A repeated affine root has
reconstruction
coordinate `x=2/P'(t)` and supplies an escaping branch as `P'(t)->0`.
Conversely, off `V(Q)` every finite root remains simple, while the only root
which can tend to projective infinity is the regular `s`-chart branch (8).
Hence the nonproperness set is exactly `V(Q)`.  This proves C03 without using
the weighted theory.

## Theorem C: the weighted normalized marked-root space

Let `k` have characteristic zero.  Let `H in k[W]` have degree `n>=3` and

\[
H(0)=H(1)=0,\quad H'(0)=0,\quad H'(1)=-c\ne0,
\quad \kappa=H''(1)/c\ne-2.
\]

Choose `b_0 ne 0`, put

\[
a_0=-{1+\kappa\over2+\kappa},\quad
v=xy,\quad S=x^2z,\quad u=1+v,
\]

\[
\gamma=1+a_0v+b_0S,\qquad W=u\gamma.
\]

With `p=H'` and `q=(Wp-H)/c`, define

\[
C=x\gamma,\qquad
B={c+p(W)/\gamma\over x},\qquad
A={u+q(W)/\gamma^2\over x^2}.                    \tag{9}
\]

The endpoint conditions make (9) polynomial and

\[
\det DG_H=b_0c.                                   \tag{10}
\]

For

\[
I_H=V(H(W)-BCW+cAC^2)\subset\mathbb A^3\times\mathbb A^1,
\]

let `tilde(I)_H` be its normalization.  On this normalization define `R_H`
to be the simultaneous regularity locus of

\[
\gamma={BC-H'(W)\over c},\quad
x={C\over\gamma},\quad
y={W-\gamma\over C},
\]

\[
z={\gamma-1-a_0(W/\gamma-1)\over b_0x^2}.         \tag{11}
\]

Then the marked-root morphism lifts to an isomorphism

\[
\mathbb A^3\xrightarrow{\sim}R_H,                 \tag{12}
\]

and `G_H` is the restriction of the finite degree-`n` root projection.  Over
`C ne 0`, `R_H` is exactly the simple-root locus.

### Polynomiality and determinant

The divisibilities `W|p(W)` and `W^2|q(W)` remove the displayed powers of
`gamma`.  At `x=0`, the numerator of `B` vanishes because `p(1)=-c`.  The
numerator of `A` and its first `x`-derivative vanish because

\[
q(1)=-1,\qquad q'(1)=\kappa,\qquad
1+\kappa+a_0(2+\kappa)=0.
\]

Thus the remaining quotients by `x` and `x^2` are polynomial.  For (10), put

\[
s=BC=c\gamma+p(W),\qquad
t=cAC^2=cW\gamma+Wp(W)-H(W).
\]

The determinant chain

\[
x^3\cdot b_0\gamma^2\cdot(-c^2\gamma)
=(-cC^3)\det DG_H
\]

and `C=x gamma` give (10) on `C ne 0`, hence everywhere.

### Normalization argument

The incidence polynomial is primitive and linear in `A`: its two
coefficients `cC^2` and `H(W)-BCW` are coprime.  Thus `I_H` is integral.  It
is finite flat of degree `n` over target space because its defining equation
is monic in `W` after division by the leading coefficient of `H`.

On `C ne 0`, equations (11) invert the marking exactly when
`E_W=H'(W)-BC ne 0`; this gives a dense birational isomorphism with the
simple-root incidence.  Since `A^3` is normal, the marking factors uniquely
through `tilde(I)_H`.  Equation (10) makes the lift quasi-finite.  Zariski's
Main Theorem makes this birational lift an open immersion into the normal
scheme `tilde(I)_H`.

The functions (11) are the inverse on the dense simple-root open.  Wherever
they are regular, the same polynomial identities extend and make them inverse
to the open immersion.  Its image is therefore exactly `R_H`, proving (12).
This argument includes all normalized branches over `C=0` and uses no
classification of later dicritical or exceptional divisors.

## Theorem D: irreducibility and full symmetric monodromy

Let `k` have characteristic zero, let `H in k[W]` have degree `n>=2`, and set

\[
E(W)=H(W)-sW+t.
\]

Then `E` is geometrically irreducible over `k(s,t)`, its repeated-root
discriminant is irreducible with birational normalization

\[
r\longmapsto(H'(r),rH'(r)-H(r)),                  \tag{13}
\]

and both its geometric and arithmetic Galois/monodromy groups are `S_n`.

### Irreducibility and discriminant normalization

In a factorization of `E` in `bar(k)[s,W,t]`, one factor has `t`-degree zero.
Comparison of the coefficient of `t` forces that factor to be a unit.
Gauss's lemma gives geometric irreducibility over `bar(k)(s,t)`.

A repeated root satisfies `s=H'(r)` and
`t=rH'(r)-H(r)`, proving (13).  Its two coordinates have pole orders `n-1`
and `n` at infinity.  The generic degree of (13) divides both orders and is
therefore one.  The discriminant is irreducible and its generic polynomial
has one double root, because `H''(r) ne 0` generically.

### Algebraic monodromy

Let `G` be the geometric Galois group and let `N` be the normal subgroup
generated by codimension-one inertia over the discriminant.  Generic inertia
is a transposition by the ordinary double-root local equation.  The fixed
cover for `N` is unramified in codimension one over `A^2_(s,t)`; purity makes
it finite etale everywhere.  Since

\[
\pi_1^{\rm et}(\mathbb A^2_{\bar k})=1,
\]

the fixed cover is trivial and `N=G`.  Thus `G` is generated by
transpositions.  Geometric irreducibility makes its action transitive; the
graph of its generating transpositions is connected, so those edge
transpositions generate `S_n`.  Hence `G=S_n`.  The arithmetic group contains
`G` and lies in `S_n`, so it too is `S_n`.

The only external inputs in this last paragraph are
[Zariski--Nagata purity](https://stacks.math.columbia.edu/tag/0BMB) and the
triviality of finite etale covers of affine space over an algebraically closed
characteristic-zero field.  For the latter, a finite cover descends to a
finitely generated characteristic-zero subfield, which embeds in `C`; the
[Riemann existence comparison](https://firmaprim.github.io/sga/sga-1/)
(SGA 1, Expose XII, Theorem 5.1 and Corollary 5.2) reduces it to the ordinary
topological fact that complex affine space is contractible.

## Reproduction and proof boundary

Run the complete local foundation suite with

```text
make verify-foundations
```

Add Dean Cureton's pinned Lean C01 build with

```text
make verify-foundations-formal
```

The exact scripts prove all displayed finite identities and both chart
compositions.  The all-`n` quantifiers in Theorems C and D are proved by the
written arguments above; bounded-degree symbolic examples are regression
tests only.  No later claim, and in particular no affine-rigidity
conjecture, is used.
