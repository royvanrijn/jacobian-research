# Fixed-map local solubility, global points, and height densities

These are written deductions for the single target line of the fixed map
`Phi` in the [fixed-map manuscript](../../unpublished/fixed-map-hasse-failures/main.tex).
They use its exact rational classification and rational-height counting theorem;
the analytic input is not reproved here. They do not classify the remaining
targets of `A^3`, and make no new formal-verification or external-review claim.

For `a in Q`, write

\[
 y_a=(-1,32a/9,(8a+1)/3),\qquad K=\mathbb Q(\zeta_3).
\]

The manuscript identifies the literal fibre with the derivative-localized root
scheme of `f_a(X)=(X^3-a)(X^2+X+1)`. For `a != 0,1`, this polynomial is
separable and the fibre has full rank five. Localization removes multiple roots.

## 1. Exact arithmetic classification on the line

**Theorem 1 (all local and global fibres).** The fibre at `y_a` has a point
over every completion of `Q` precisely when `a != 0` and

\[
 3\mid v_p(a)\quad(p\equiv2\pmod3),\qquad
 a\in\mathbb Q_3^{\times3}.                                      \tag{1}
\]

It has a rational point precisely when `a` is a nonzero rational cube.
Thus its Hasse failures are precisely the locally soluble noncubes.
The abstract fibre algebras are

| Parameter | Literal fibre algebra | Rank | Rational point |
|---|---|---:|---|
| `a=0` | `K` | 2 | no; also no real point |
| `a=1` | `Q` | 1 | yes |
| `a=c^3`, `c in Q^*`, `c != 1` | `Q x K x K` | 5 | yes |
| `a` a noncube | `Q(cuberoot(a)) x K` | 5 | no |

In the last row, (1) decides local solubility. All *globally soluble rank-five*
fibres on this line therefore have the same abstract algebra.

**Proof.** Away from `0,1`, apply the manuscript's rational classification.
For `a=c^3 != 0,1`, factor

\[
 X^3-c^3=(X-c)(X^2+cX+c^2).
\]

The quadratic has roots `c zeta_3, c zeta_3^2`, so its quotient is `K`.
All three factors are coprime, giving the displayed algebra by the Chinese
remainder theorem. For a noncube, the cubic is irreducible. At `a=0`,
derivative localization removes `X^3`; at `a=1`, the factorization is
`(X-1)(X^2+X+1)^2`, and only `X-1` survives. The latter fibre is locally
soluble and the former is not. Every nonzero rational cube satisfies (1).
This also checks the exceptional parameters. \(\square\)

**Theorem 1.1 (isomorphism classes).** For noncubes `a,b in Q^*`,

\[
 \mathbb Q(\sqrt[3]{a})\simeq_{\mathbb Q}\mathbb Q(\sqrt[3]{b})
 \quad\Longleftrightarrow\quad
 b/a\in\mathbb Q^{\times3}\ \hbox{or}\ b/a^2\in\mathbb Q^{\times3}. \tag{1a}
\]

Consequently their rank-five fibre algebras are isomorphic precisely under
this condition. In particular the Hasse-failing fibre algebras are classified
by the nonzero locally admissible classes in `Q^*/Q^{*3}`, modulo inversion.
Here locally admissible means satisfying (1), which is invariant under
multiplication by a rational cube and under inversion.

**Proof.** Put `z=cuberoot(a)` and `M_a=K(z)`, the splitting field over `Q`.
A noncube in `Q` cannot become a cube in the quadratic field `K`: if `q=c^3`
with `q in Q^*` and `c in K^*`, then

\[
 N_{K/\mathbb Q}(c)^3=q^2,\qquad
 q=\left(q/N_{K/\mathbb Q}(c)\right)^3.
\]

Thus `M_a/K` is cyclic of degree three, with an automorphism `sigma(z)=zeta_3 z`.
An isomorphism of the cubic fields extends to their normal closures. Inside
a fixed algebraic closure these normal closures must be equal, since a normal
extension is preserved as a set by every `Q`-automorphism. If `w^3=b` in this
common field, `sigma(w)=zeta_3^j w` with `j=1` or `2`; `j=0` is excluded by
the preceding cube-injectivity argument. The eigenspaces in the basis
`1,z,z^2` are `Kz` and `Kz^2`, hence `w=c z^j` for `c in K^*`. Therefore
`b/a^j=c^3` is rational, and the same norm argument makes it a rational cube.
Conversely `b=q^3a` or `b=q^3a^2` gives the same cubic field, since
`Q(z^2)=Q(z)` using `z=(z^2)^2/a`. Finally the product algebra `L x K` has a
unique field component of dimension three, so its isomorphisms preserve the
cubic field component. This proves the fibre-algebra assertion. \(\square\)

## 2. All targets in projective height

Let `U(B), L(B), G(B), N(B)` count all targets, everywhere locally soluble
fibres, globally soluble fibres, and Hasse failures on this line, respectively,
with projective height at most `B`. Since `a -> y_a` is injective, Theorem 1 gives

\[
 L(B)=N(B)+G(B)\quad\hbox{for every }B.                            \tag{2}
\]

For reduced `a=m/n`, `n>0`, `m != 0`, the manuscript proves

\[
 H(y_{m/n})=
 \frac{\max(9n,32|m|,|24m+3n|)}{\gcd(n,8)\gcd(m,3)}.             \tag{3}
\]

The parameter `a=0` has height three and contributes only `O(1)` to `U(B)`.
Put

\[
 \Omega=\{(x,y):y>0,\ \max(9y,32|x|,|24x+3y|)\le1\}.
\]

Its negative and positive parts have areas `72/20736` and `71/20736`.

**Theorem 2 (all target count).**

\[
 U(B)=C_{\rm all}B^2+O(B\log(2B)),\qquad
 C_{\rm all}=\frac{143}{864\zeta(2)}.                            \tag{4}
\]

**Proof.** Use congruence classes modulo 24, which determine the content in
(3) and primitivity at 2 and 3. Multiplication by an integer coprime to 6
preserves these conditions. The remaining Mobius sieve has Euler factor

\[
 \prod_{p\ne2,3}(1-p^{-2})=\frac{3}{2\zeta(2)}.
\]

At 2 the disjoint cases `n` odd, `v_2(n)=1`, `v_2(n)=2`, and `8|n` have
measures `1/2,1/8,1/16,1/16` after imposing primitive parity. Their content
squares are `1,4,16,64`. At 3, the cases `3` not dividing `m`, and
`3|m, 3` not dividing `n`, have measures `2/3,2/9` and weights `1,9`. Thus

\[
 D_2=\tfrac12+\tfrac48+\tfrac{16}{16}+\tfrac{64}{16}=6,
 \qquad D_3=\tfrac23+9\tfrac29=\tfrac83.
\]

The coefficient is consequently

\[
 \operatorname{area}(\Omega)D_2D_3\frac{3}{2\zeta(2)}
 =\frac{143}{20736}\frac{24}{\zeta(2)}=C_{\rm all}.
\]

For each congruence class and content `d`, the polygon lattice error after
dividing by a sieve integer `k` is `O(B/k+1)`. The nonzero-numerator sieve
truncates at `k=O(B)` because numerator and denominator are `O(B)`.
Summing errors costs `O(B log(2B))`, and completing the main Euler sum costs
`O(B)`. Handle the zero numerator separately. \(\square\)

## 3. Global cubes with an error term

Define

\[
 \Gamma=\{(u,v):v>0,\ \max(9v^3,32|u|^3,|24u^3+3v^3|)\le1\},
 \qquad A=\operatorname{area}(\Gamma).
\]

Its negative part is a rectangle. The positive part is a rectangle and a
strip under a cubic arc. Their areas are

\[
 \begin{aligned}
 A_-&=288^{-1/3},\\
 A_+&=324^{-1/3}
      +\int_{36^{-1/3}}^{32^{-1/3}}
                  \left(\frac{1-24u^3}{3}\right)^{1/3}\,du,\\
 A&=A_-+A_+\simeq0.302598873573506.
 \end{aligned}                                                   \tag{5}
\]

The horizontal boundary switches to the arc at `u^3=1/36`; on this arc `v`
is bounded away from zero. The region has a piecewise smooth boundary,
including the axis segments of its closure.

**Theorem 3 (globally soluble fibres).**

\[
 G(B)=C_{\rm glob}B^{2/3}+O(B^{1/3}\log(2B)),\qquad
 C_{\rm glob}=\frac{(3+3^{2/3})A}{2\zeta(2)}
              \simeq0.467261172801867.                           \tag{6}
\]

The count includes `a=1` and excludes `a=0`. Restricting by sign gives the
same error term and constants obtained by replacing `A` with `A_+` or `A_-`.

**Proof.** Write a nonzero rational cube uniquely as `a=(r/s)^3`, with
`r != 0`, `s>0`, `gcd(r,s)=1`. Its content in (3) is

\[
 d=8^e3^f,\qquad e=\mathbf1_{2\mid s},\quad f=\mathbf1_{3\mid r}.
\]

Put `T=B^(1/3)`. The height condition is `(r,s) in d^(1/3) T Gamma`.
Primitive-at-2,3 residue pairs modulo 6 have counts

\[
 R_{00}=12,\quad R_{01}=4,\quad R_{10}=6,\quad R_{11}=2.
\]

Each allowed set is invariant under simultaneous multiplication by an integer
coprime to 6. For each such `k`, the lattice estimate is therefore

\[
 \#\bigl((d^{1/3}T/k)\Gamma\cap\mathbb Z^2
           \text{ in the allowed classes}\bigr)
 =\frac{R_{ef}}{36}A\frac{d^{2/3}T^2}{k^2}+O(T/k+1).
\]

Retain `r != 0,s>0`. Mobius inversion truncates at `k=O(T)`, with total
error `O(T log(2T))`; completing the main Euler sum costs `O(T)`.
The Euler factor is `3/(2 zeta(2))` and the finite weighted residue sum is

\[
 \sum_{e,f}\frac{R_{ef}}{36}(8^e3^f)^{2/3}
 =1+\frac{3^{2/3}}3.
\]

Their product proves (6). The same proof in each half of the region gives
the signed constants. \(\square\)

## 4. Relative densities and Hasse failures

Write

\[
 \mathcal P=\prod_{p\equiv2\ (3)}\frac{1+p^{-3}}{1-p^{-3}},\qquad
 C_\Phi=\frac{11\mathcal P}{1728\sqrt3\,\zeta(2)}.
\]

The product converges absolutely. The manuscript's rational-height theorem is

\[
 N(B)\sim C_\Phi\frac{B^2}{\log B}.                              \tag{7}
\]

**Corollary 4 (density comparisons).**

\[
 \begin{aligned}
 L(B)&\sim C_\Phi B^2/\log B,\\
 \frac{L(B)}{U(B)}&\sim\frac{\mathcal P}{26\sqrt3\log B},\\
 \frac{G(B)}{L(B)}&\sim
       \frac{C_{\rm glob}}{C_\Phi}\frac{\log B}{B^{4/3}},\\
 \frac{N(B)}{L(B)}&=
       1-\frac{C_{\rm glob}}{C_\Phi}
                     \frac{\log B}{B^{4/3}}(1+o(1)).
 \end{aligned}                                                   \tag{8}
\]

Locally soluble targets have density zero among all rational targets on this
line. Hasse failures have density one among its locally soluble targets.
These statements concern the specified height ordering and line only.

**Proof.** Combine (2), (4), (6), and (7), and use
`C_Phi/C_all = P/(26 sqrt(3))`. \(\square\)

The analytic error in (7) has not been proved to be `o(B^(2/3))`. One therefore
cannot append a justified second term `-C_glob B^(2/3)` to its single leading
term. The valid refinement, using (2), is

\[
 N(B)=L(B)-C_{\rm glob}B^{2/3}+O(B^{1/3}\log(2B)).                \tag{9}
\]

**Corollary 5 (sign distribution).** With subscripts restricting the sign of `a`,

\[
 N_+(B)\sim\frac{71}{143}C_\Phi\frac{B^2}{\log B},\qquad
 N_-(B)\sim\frac{72}{143}C_\Phi\frac{B^2}{\log B}.                \tag{10}
\]

The same proportions hold for `L` and `U`. For `G`, the proportions are instead
`A_+/A` and `A_-/A`, approximately `0.499579379523` and `0.500420620477`.

**Proof.** In the manuscript's polygon-counting proof, retain one sign region
at every valuation stratum. All arithmetic weights, including the modulo-nine
cube condition, are unchanged by sign because `-1` is a cube. Only the area
changes. The areas of `Omega_+,Omega_-` are `71/20736,72/20736`. This proves
(10); the all-target proof has the same property. The global count is negligible
for `L`, and its own signed version follows from Theorem 3. \(\square\)

## 5. Reproducible finite checks and assurance

Run from the repository root:

```sh
python3 research/arithmetic/verify_fixed_map_local_densities.py --output research/arithmetic/fixed_map_local_densities.verification.json
```

The checker uses Python's standard library. It compares projective-coordinate
GCD height against the content formula, direct reduced-parameter cube counts
against independent primitive-cube enumeration, exact local and exceptional
controls, and the finite local factors. Direct enumeration ends at height 300;
cube enumeration also reaches height `10^6`. The compact
[receipt](fixed_map_local_densities.verification.json) records the actual domain,
counts, source hashes, and numerical constants.

These finite checks do not establish the asymptotic theorems, their error terms,
the unbounded local criterion, or the imported Selberg-Delange application.
Those claims depend on the written arguments and the manuscript's identified
inputs. Independent mathematical and external specialist review remain absent.
