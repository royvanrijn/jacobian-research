# Danielewski multi-dicritical maps with permutation reductions

The smooth balanced triple-root map extends to a complete family with
arbitrarily many geometric dicritical divisors.  Every map is an open
immersion from affine three-space into a smooth Danielewski threefold, has
constant residue Jacobian, and has an exact boundary-complement
decomposition.

Arithmetic twisting of the missing boundary roots then produces bijective
reductions over finite fields.  In particular,

\[
P(x)=x(x^2+1)
\]

gives two geometric dicritical divisors and bijections on rational points
for every prime `p=3 mod 4`.

The remaining gap is that the smooth target is not affine three-space when
any dicritical divisor is present.

The exact next obstruction is isolated in the
[Poisson-contraction frontier](DANIELEWSKI_POISSON_CONTRACTION_FRONTIER.md):
an \(a\)-preserving affine target would require a constant-bracket pair in
the Danielewski coordinate ring, but the nonzero periods of the residue
form exclude every such pair when \(\deg P>1\).  Any remaining contraction
must mix \(a\) into at least two target coordinates.

## 1. General construction

Let `k` have characteristic zero, and let

\[
P(x)=xQ(x)\in k[x]
\]

be squarefree of degree `n>=1`.  Define

\[
T_P=
\{cw=P(x)\}
\subset\mathbb A^4_{a,c,x,w}.                       \tag{1}
\]

The map

\[
\boxed{
\Phi_P:\mathbb A^3_{a,b,c}\longrightarrow T_P,
\qquad
(a,b,c)\longmapsto
\left(a,c,bc,bQ(bc)\right)
}                                                     \tag{2}
\]

is polynomial because

\[
c\,bQ(bc)=bc\,Q(bc)=P(bc).
\]

## 2. Smoothness and residue Jacobian

For `H=cw-P(x)`, the gradient in `(c,w,x)` is

\[
(w,c,-P'(x)).
\]

If it vanished on `H=0`, then `c=w=0`, `P(x)=P'(x)=0`, contradicting
squarefreeness.  Therefore `T_P` is smooth.

Use

\[
\Omega_P=\frac{da\wedge dc\wedge dx}{c}.            \tag{3}
\]

Since

\[
\det\frac{\partial(a,c,bc)}{\partial(a,b,c)}=-c,
\]

one has

\[
\boxed{
\Phi_P^*\Omega_P=-da\wedge db\wedge dc.
}                                                     \tag{4}
\]

Every map in the family has constant residue Jacobian `-1`.

## 3. Exact image and dicritical divisors

On `c!=0`, the inverse is `b=x/c`.  On `c=0`, equation (1) forces `P(x)=0`.

At the root `x=0`, formula (2) gives

\[
w=bQ(0)=bP'(0),
\]

so this component is also in the image and reconstructs

\[
b=\frac{w}{P'(0)}.
\]

For a nonzero geometric root `alpha` of `P`, the divisor

\[
D_\alpha=\{c=0,\ x=\alpha\}\simeq\mathbb A^2_{a,w}  \tag{5}
\]

has no affine-source preimage: `x=bc` would vanish when `c=0`.  At its
generic point, `c` is a uniformizer and `b=x/c` has valuation `-1`.
Thus `D_alpha` is dicritical.

Consequently

\[
\boxed{
\Phi_P(\mathbb A^3)
=T_P\setminus
\coprod_{\substack{P(\alpha)=0\\\alpha\ne0}}D_\alpha.
}                                                     \tag{6}
\]

The map is an open immersion, and it has exactly `n-1` geometric
dicritical divisors.

## 4. Motivic ledger

Over an algebraic closure, stratify the Danielewski surface
`D_P={cw=P(x)}` by the `n` simple roots of `P`:

\[
[D_P]
=(\mathbb L-n)(\mathbb L-1)+n(2\mathbb L-1)
=\mathbb L^2+(n-1)\mathbb L.                        \tag{7}
\]

Therefore

\[
\boxed{
[T_P]
=\mathbb L^3+(n-1)\mathbb L^2
=[\mathbb A^3]+
\sum_{\alpha\ne0}[D_\alpha].
}                                                     \tag{8}
\]

The motivic excess is exactly one affine plane per dicritical divisor.
In particular,

\[
T_P\simeq\mathbb A^3
\quad\Longrightarrow\quad n=1,
\]

and for `n=1` the map is the identity presentation
`cw=x`, `w=b`.

Within this family, affine target and nonzero dicritical boundary are
therefore mutually exclusive.

## 5. Finite fields

Let `P` have coefficients in `F_q`, remain squarefree, and let

\[
r_q=\#\{\alpha\in\mathbb F_q:P(\alpha)=0\}.
\]

The same fiber stratification gives

\[
\boxed{
\#T_P(\mathbb F_q)
=q^3+(r_q-1)q^2.
}                                                     \tag{9}
\]

The map (2) is injective, and its complement consists of one rational
affine plane for every nonzero rational root of `P`.  Hence

\[
\boxed{
\Phi_P:\mathbb A^3(\mathbb F_q)\longrightarrow
T_P(\mathbb F_q)
\text{ is bijective}
\iff r_q=1.
}                                                     \tag{10}
\]

Only the distinguished root `0` may be rational.

## 6. Two dicriticals and infinitely many bijective primes

Take

\[
P(x)=x(x^2+1).                                      \tag{11}
\]

Over an algebraic closure, the nonzero roots `i,-i` give two distinct
dicritical divisors.  For an odd prime `p`,

\[
r_p=1
\iff -1\text{ is not a square in }\mathbb F_p
\iff p\equiv3\pmod4.                                \tag{12}
\]

Thus (2) becomes

\[
(a,b,c)\longmapsto
\left(a,c,bc,b((bc)^2+1)\right),                    \tag{13}
\]

and induces a bijection

\[
\boxed{
\mathbb A^3(\mathbb F_p)
\overset{\sim}{\longrightarrow}
T_P(\mathbb F_p)
\quad\text{for every }p\equiv3\pmod4.
}                                                     \tag{14}
\]

There are infinitely many such primes.  This realizes the proposed
arithmetic strategy: geometric dicritical components can be retained while
their Galois orbit contributes no rational boundary points over selected
finite fields.

This is a bijection between the affine source and a smooth Cox target, not
a permutation polynomial self-map of affine three-space.

## 7. Reproduction

Run

```bash
.venv/bin/python scripts/verify_danielewski_multi_dicritical_family.py
```

The checker verifies smoothness, residue Jacobian, exact open-immersion
complement, motivic and finite-field ledgers, and the bijections for the
quadratic twist over several primes `p=3 mod 4`.
