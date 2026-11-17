# Poisson-square rigidity for the three-layer band pattern

## Result

Consider

\[
 P=X^3A(t)+X^2B(t),\qquad
 Q=X^2C(t)+XD(t)
\]

with

\[
 \deg A\le3,\quad \deg B\le4,\quad
 \deg C\le2,\quad \deg D\le3
\]

over a characteristic-zero field.  The exact band compiler in
[`cas/poisson_square_rigidity.py`](cas/poisson_square_rigidity.py) expands

\[
 [P,Q]_{X,t}=X^2
\]

into precisely

\[
\begin{aligned}
3AC'-2A'C&=0, &(K4)\\
3AD'-A'D+2BC'-2B'C&=0, &(K3)\\
2BD'-B'D&=1. &(K2)
\end{aligned}
\]

There are 17 scalar coefficient equations before using this triangular
structure.

The entire geometric reduced solution set has three irreducible components.
On the component where `AC` is nonzero, every point is exactly

\[
\boxed{
\begin{aligned}
h&=t-r,\\
A&=ah^3,\\
C&=-\frac32ad^2h^2,\\
D&=d+eh,\\
B&=-\frac{h}{d}-\frac{e}{2d^2}h^2,
\end{aligned}}
\qquad a,d\ne0.                                      \tag{1}
\]

The other two components are `C=0` and `A=0`; explicit charts for them are
given in Section 3.  Thus the weighted tangent model is forced, up to its
natural affine and scaling parameters, at every field-valued point with
`AC != 0`, and every remaining field-valued architecture is classified.

## 1. Highest-layer rigidity

Factor `A` and `C` over an algebraic closure.  At a root where their
multiplicities are `m,n`, equation `(K4)` gives

\[
 3n-2m=0.
\]

Hence `m=3k`, `n=2k`.  The exact degree bounds leave one common root and
`k=1`, proving

\[
 A=ah^3,\qquad C=ch^2,\qquad h=t-r.
\]

This is the same forced cube/square seen in the weighted tangent pair.  It is
now a consequence of the band support and top Poisson layer, not of the
three-dimensional construction.

## 2. The remaining classification

Write

\[
B=b_0+b_1h+\cdots+b_4h^4,\qquad
D=d_0+d_1h+d_2h^2+d_3h^3.
\]

The middle layer gives

\[
b_0=0,\quad
b_1=\frac{3ad_0}{2c},\quad
b_3=\frac{3ad_2}{2c},\quad
b_4=\frac{3ad_3}{2c},
\]

leaving only `b2,d0,d1,d2,d3`.  After substitution, the nonzero coefficients
of `(K2)-1`, from degree six down to degree zero, are

\[
\begin{array}{c|l}
6&3ad_3^2/c\\
5&9ad_2d_3/(2c)\\
4&-(6ad_1d_3-3ad_2^2-8b_2cd_3)/(2c)\\
3&(3ad_0d_3-3ad_1d_2+4b_2cd_2)/(2c)\\
1&d_0(3ad_1-4b_2c)/(2c)\\
0&-3ad_0^2/(2c)-1.
\end{array}
\]

At a field-valued point, the degree-six equation gives `d3=0`, the
degree-four equation then gives `d2=0`, and the constant equation makes
`d0` nonzero.  The last two equations give

\[
 c=-\frac32ad_0^2,\qquad
 b_2=\frac{3ad_1}{4c}.
\]

Renaming `d0=d` and `d1=e` yields (1).

## 3. Complete degenerate classification

The bottom equation alone is already rigid:

\[
 2BD'-B'D=1.                                           \tag{2}
\]

It makes `B` and `D` coprime.  Comparing highest coefficients excludes
`deg(D)=3` and `deg(D)=2`: after making `D` monic and centered, their
coefficient equations successively kill all available coefficients of `B`
and leave the constant contradiction `-1=0`.  Consequently `deg(D)<=1`.
All solutions of (2) are

\[
\begin{array}{c|c}
D=d,\ d\ne0 & B=b-t/d,\\[2mm]
D=d_0+\delta t,\ \delta\ne0
  &B=kD^2+\dfrac1{2\delta}.
\end{array}                                             \tag{3}
\]

If `C=0`, the middle layer is `3AD'-A'D=0`.  It gives

\[
\begin{array}{c|c}
D=d & A=\alpha,\\
\deg D=1 & A=\alpha D^3.
\end{array}                                             \tag{4}
\]

If `A=0`, the middle layer is `2(BC'-B'C)=0`.  Coprimality of
`B,D` and (3) give

\[
 C=\lambda B.                                          \tag{5}
\]

These are four affine charts, with the constant and linear `D` charts joining
in the closure.  They form two irreducible components: one with `C=0`, and
one with `A=0`.

Finally suppose both `A` and `C` are nonzero.  The root-multiplicity argument
for `(K4)` gives

\[
 A=aH^3,\qquad C=cH^2
\]

for a common polynomial `H`.  The degree bounds give `deg(H)<=1`.  If `H`
were constant, `(K3)` would be incompatible with (2); hence `H` is affine
and Section 2 yields the tangent family (1).  This proves that (1), (3)--(5)
cover every field-valued point in the support box.

An independent exact radical computation agrees with the proof.  Singular
finds dimension four and exactly three dimension-four minimal primes: the
tangent-component closure, `C=0`, and `A=0`.  The executable audit is
[`cas/poisson_square_radical.sing`](cas/poisson_square_radical.sing).

## 4. Exact tangent-pencil identification

Put

\[
 W=Xh,\qquad \gamma=dX,\qquad
 H(W)=\frac e2W^2-\frac{ad^2}{2}W^3.
\]

Then the family (1) satisfies

\[
\boxed{
 Q=H'(W)+\gamma,\qquad
 -d^2P=W\bigl(H'(W)+\gamma\bigr)-H(W).}               \tag{6}
\]

Therefore `(Q,-d^2P)` is exactly the tangent map `Phi_H`.  If `e` is
nonzero, affine rescaling normalizes the quadratic and cubic coefficients of
`H` to the foundational form `W^2(1-W)`.  If `e=0`, the result is the
triple-root cubic degeneration; it is still a tangent pencil but has smaller
actual lower-band support.

The foundational point is

\[
 r=1,\qquad a=\frac12,\qquad d=2,\qquad e=-4,
\]

which recovers

\[
A=\tfrac12(t-1)^3,\quad
B=\tfrac12(t-2)(t-1),\quad
C=-3(t-1)^2,\quad
D=-2(2t-3).
\]

## 5. Scheme-theoretic caveat

The classification above completely describes the geometric reduced locus,
but not its full scheme structure.  On the tangent component, the constant
and linear equations make `d0` a unit and impose the displayed relation for
`b2`; the cubic equation then forces `d3=0`.  The degree-four equation
reduces to a nonzero multiple of

\[
 d_2^2=0.
\]

This nilpotent is genuine.  For example, over
`k[epsilon]/(epsilon^2)`, take

\[
 a=1,\quad c=-3/2,\quad d_0=1,\quad
 d_1=b_2=d_3=0,\quad d_2=\varepsilon.
\]

Every residual coefficient vanishes, while `d2` need not.  Higher-band
infinitesimal directions therefore survive even though every field-valued
point is tangent type.

The `C=0` component has a longer generic curvilinear thickening.  In its
affine `D` chart from (3)--(4), write `z=epsilon*c`, keep `B,D` fixed, and set

\[
\begin{aligned}
C_\varepsilon&=zD^2-\frac{2z^2}{3\alpha\delta},\\
A_\varepsilon&=\alpha D^3-
\left(\frac z\delta+\frac{4kz^2}{3\alpha\delta}\right)D,
\qquad \varepsilon^3=0.
\end{aligned}                                          \tag{7}
\]

All three layers still equal `(0,0,1)` modulo `epsilon^3`, but the displayed
residual is nonzero modulo `epsilon^4`.

This length is exact generically.  Put `s=D`,
`B=ks^2+1/(2delta)`, and write
`C=c0+c1*s+c2*s^2`.  The middle layer determines

\[
A=-\frac{c_1}{3\delta}
 +\left(2kc_0-\frac{c_2}{\delta}\right)s
 +2kc_1s^2+\alpha s^3.
\]

On the open set where `alpha` and
`L=-3alpha+4k*c2` are units, the top layer is locally equivalent to

\[
c_1=0,\qquad
\delta Lc_0=2c_2^2,\qquad
c_2^3=0.                                               \tag{8}
\]

Thus its transverse local algebra is exactly `k[c2]/(c2^3)`.

The Jacobian of the original 17 coefficient equations gives the following
exact ranks at rational points in the generic charts:

| Reduced component | Jacobian rank | Scheme tangent dimension | Reduced dimension |
| --- | ---: | ---: | ---: |
| tangent closure | 11 | 5 | 4 |
| `C=0` | 11 | 5 | 4 |
| `A=0` | 12 | 4 | 4 |

Together with the tangent slice `d2^2=0`, this gives generic minimal-component
multiplicities

\[
 2\quad\text{(tangent)},\qquad
 3\quad\text{(`C=0`)},\qquad
 1\quad\text{(`A=0`)}.                                 \tag{9}
\]

The checker verifies (7)--(8), constructs the sample points, and computes the
tangent ranks exactly.

This is not an incidental nuisance: the normalized sixteen-monomial weighted
ansatz has exactly such a dual-number direction, detected by the square
`(B01-3)^2`.  The new result says:

- no unclassified field-valued architecture occurs anywhere in the
  three-layer band box;
- its reduced scheme has exactly three irreducible components;
- nilpotent thickening data can still differ and must be retained when a
  coefficient scheme, rather than its point set, matters.

The exact generic multiplicities are therefore known.  What remains of the
global primary-decomposition gap is the presence or absence of embedded
components and the complete primary structure where the three minimal
components meet.

The reduced incidence is now exact.  Write the three components as `T`
(tangent closure), `C0` (`C=0`), and `A0` (`A=0`).  There are three
irreducible intersection branches, with dense charts

\[
\begin{array}{c|c}
S&
A=C=0,\quad D=d_0+\delta t,\quad
B=kD^2+\frac1{2\delta}\\[2mm]
S_C&
C=0,\quad D=d,\quad B=b-\frac td,\quad A=\alpha\\[2mm]
S_A&
A=0,\quad D=d_0+\delta t,\quad
B=\frac1{2\delta},\quad C=c .
\end{array}
\]

The exact radical identities are

\[
 T\cap C0=S\cup S_C,\qquad
 T\cap A0=S\cup S_A,\qquad
 C0\cap A0=T\cap C0\cap A0=S.                       \tag{10}
\]

All three branches have dimension three.  Exact coefficient-Jacobian ranks
on their dense charts are respectively `8,9,10`, so the scheme tangent
dimensions are `8` on `S`, `7` on `S_C`, and `6` on `S_A`.  The Singular
audit proves the radical equalities and branch counts `(2,2,1)`; the Python
audit supplies the displayed charts and tangent ranks.

Chosen tangent-kernel slices distinguish the two non-core branches even at
the nilpotent level.  After fixing three rational parameters on `S_C`, the
four-dimensional tangent-kernel slice is

\[
 k[u_0,u_1,u_2,u_3]/(u_0,u_1,u_2,u_3)^2,             \tag{11}
\]

of length five, Hilbert vector `(1,4)`, and socle dimension four.  On `S_A`
the corresponding three-dimensional kernel slice is

\[
\frac{k[u_0,u_1,u_2]}{
(u_2^3,u_0^2,u_0u_1-\frac9{14}u_2^2,u_1^2,
 u_0u_2,u_1u_2-3u_2^2)},                             \tag{12}
\]

also of length five, but with Hilbert vector `(1,3,1)` and one-dimensional
socle.  Thus the two boundary collisions have different exact nilpotent
normal signatures.  These are deliberately described as chosen
tangent-kernel slices: they do not identify the full completed local rings.

The tangent-kernel slices do not by themselves determine the completed local
rings.  Section 6 now goes further in a different direction: it gives the
complete associated-prime set and a two-level primary filtration of the
global coefficient algebra.

## 6. Exact primary filtration on the principal chart

The bottom constant equation is stronger than a convenient chart choice:

\[
 2b_0d_1-b_1d_0=1.                                   \tag{13}
\]

It first shows that `D(d0)` and `D(d1)` cover the coefficient scheme.  The
new characteristic-zero saturation audit proves the sharper identity

\[
 I:d_0^\infty=I.                                      \tag{14}
\]

Thus multiplication by `d0` on the full coefficient algebra is injective.
In particular, no associated prime—minimal or embedded—contains `d0`.
Every possible embedded component is therefore detected after localizing at
`d0`; there is no separate primary phenomenon hidden on the divisor
`d0=0`.

The equations carry the exact `G_m` action

\[
 D\mapsto\lambda D,\qquad B\mapsto\lambda^{-1}B,\qquad
 A\mapsto\lambda^{-2}A,\qquad C\mapsto C.
\]

On `D(d0)`, choosing `lambda=d0^{-1}` identifies the primary problem with
the normalized slice `d0=1`, times `G_m`.  The original scheme has dimension
four and standard-basis size 743; after eliminating `d0-1`, the normalized
slice has dimension three and standard-basis size 511 (512 if `d0-1` is
retained as a generator in the sixteen-variable ring).

Write `R` for the normalized fifteen-variable polynomial ring and `I_0` for
the normalized coefficient ideal.  Separator saturations for the tangent,
`C=0`, and `A=0` components have an intersection strictly larger than
`I_0`.  Thus the scheme is not unmixed.  More precisely, the top coefficient
`d3` gives the exact sequence

\[
0\longrightarrow R/(I_0:d_3)
 \mathop{\longrightarrow}^{\cdot d_3} R/I_0
 \longrightarrow R/I_1\longrightarrow0,
\qquad I_1=I_0+(d_3).                                \tag{15}
\]

The colon ideal `I0:d3` has an exact five-component primary decomposition.
Its associated primes are

\[
 S,\quad S_C,\quad S_A,\quad
 K_C=S\cap S_C,\quad K_A=S\cap S_A,                  \tag{16}
\]

with normalized dimensions `2,2,2,1,1`.  Here the intersections in (16)
mean their irreducible intersection branches, not scheme-theoretic ideal
sums.  The two curves are disjoint: `K_C` has `d1=0`, whereas `K_A` has
`2*b0*d1=1`.  The annihilator is nonradical.  If `K` is its reduced support
ideal, then

\[
 K^4d_3\subset I_0,\qquad K^3d_3\not\subset I_0.      \tag{17}
\]

The three surface-primary closures do not yet reconstruct `I0:d3`; their
defect is supported exactly on `K_C union K_A`.  Its annihilator has two
primary components, one on each curve.  This is the complete second gluing
layer, not merely a reduced support calculation.

Continue with `I2=I1+(d2)`.  A second exact sequence uses `d2`.  The colon
`I1:d2` has four primary components: the tangent minimal prime `T` and the
three surface primes `S,S_C,S_A`.  Finally, `I2` already contains `b4` and
`b3` and has exactly three primary components, all of dimension three, with
radicals `T,C0,A0`.

Applying the associated-prime inclusions for (15) and its `d2` analogue
gives the complete answer:

\[
\boxed{
\operatorname{Ass}(R/I_0)=
\{T,C0,A0;\ S,S_C,S_A;\ K_C,K_A\}.}                  \tag{18}
\]

After restoring the `G_m` factor, the full sixteen-variable scheme therefore
has exactly three dimension-four minimal primes, three dimension-three
embedded primes, and two dimension-two embedded primes.  There are no other
associated primes.  A displayed primary decomposition of `I_0` itself is
not canonical at its embedded primes; the colon filtration (15)--(18) is a
canonical and faster certificate of the complete associated-prime set.

### Transverse normal-module fibers

Both multiplication kernels in the filtration are cyclic, so their zeroth
Fitting ideals are exactly `I0:d3` and `I1:d2`.  Exact transverse slices,
computed first at rational points and then replayed over the rational
function fields of the branch parameters, give the following generic local
Hilbert vectors:

| stratum | dim. | `d3` Hilbert | `d3` socle | `d2` Hilbert | `d2` socle |
| --- | ---: | --- | ---: | --- | ---: |
| `T` | 3 | — | — | `(1)` | 1 |
| `S` | 2 | `(1,2)` | 2 | `(1,3,3,1)` | 2 |
| `S_C` | 2 | `(1,2,2)` | 2 | `(1,2,1)` | 1 |
| `S_A` | 2 | `(1,2,3,1)` | 2 | `(1,1)` | 1 |
| `K_C` | 1 | `(1,3,5,7,6,3)` | 4 | — | — |
| `K_A` | 1 | `(1,4,8,10,6,1)` | 6 | — | — |

Every corresponding reduced generic slice has length one, so the totals
`1; 3,5,7,25,30; 8,4,2` measure normal nilpotent thickness rather than
multiple reduced points in the parameter fiber.  In particular, the `d2`
surface fibers have the binomial Hilbert vectors `(1+t)^3`, `(1+t)^2`, and
`(1+t)`, while the deeper `d3` curve fibers are substantially thicker.
The symmetric vector on `S` is not Gorenstein: its socle has dimension two.
Elimination gives an exact presentation by three quadrics and one necessary
cubic.  The `S_C` fiber has an exact four-relation presentation with socle
one, and the `S_A` fiber is exactly the dual-number algebra.  Explicit
presentations are also finite on the `d3` layer:

- `S` is the square-zero algebra in two generators;
- `S_C` is the two-generator length-five algebra
  \[
  K[c_0,d_1]/(4c_0^2+15u c_0d_1+9u^2d_1^2,\ d_1^3,\ c_0d_1^2);
  \]
- `S_A` has two generators and a four-element Gröbner presentation;
- `K_C` has three generators and a nine-element Gröbner presentation;
- `K_A` has four generators and an eighteen-element Gröbner presentation.

Here `K` is the appropriate branch rational-function field and `u` is the
generic `a0` coordinate on `S_C`.  These presentations and their standard
monomial bases determine multiplication tables algorithmically.  A simpler
canonical normal form for the larger three fibers remains open.

The search-facing module
[`cas/poisson_square_filtered_modules.py`](cas/poisson_square_filtered_modules.py)
packages all eight normalized parametrizations and these certified Hilbert
vectors.  Given proposed lower coefficient equations, it substitutes them
on every dense chart and reports:

- `preserved` when they vanish identically;
- `eliminated` when their localized restriction ideal is the unit ideal;
- `cut` when only a proper closed sublocus can remain.

This is the first executable route by which the primary filtration feeds
back into later Newton layers before a global coefficient ideal is built.

## 7. Executable checks and search use

Run:

```bash
.venv/bin/python plane-jc/cas/test_poisson_square_rigidity.py
.venv/bin/python plane-jc/cas/poisson_square_rigidity.py
Singular -q plane-jc/cas/poisson_square_radical.sing
Singular -q plane-jc/cas/poisson_square_primary_charts.sing
Singular -q plane-jc/cas/poisson_square_separator_primary.sing
Singular -q plane-jc/cas/poisson_square_normalized_defect.sing
.venv/bin/python plane-jc/cas/test_poisson_square_filtered_modules.py
```

The Python checks verify the generic support compiler, all three layers,
the tangent parametrization, the four degenerate charts, the exclusion of
quadratic and cubic `D`, the tangent identity (6), the foundational
specialization, the length-two and length-three transverse families, and the
three component tangent ranks.  It also checks dense charts and tangent
dimensions for the three reduced intersection branches and the two
length-five kernel-slice algebras (11)--(12).  The Singular audit
computes the radical and its three minimal primes directly, proves (10), and
checks the intersection branch counts `(2,2,1)`.
The chart audit proves (14), verifies that the two saturated chart closures
reconstruct `I`, and records the exact dimensions and standard-basis sizes.
The separator audit proves that the three generic component closures do not
reconstruct `I` and supplies explicit separator-torsion witnesses.  The
normalized filtration audit proves (16)--(18), including the five- and
four-component colon decompositions, the final three-component quotient, and
the transverse Hilbert vectors in the table.  The Python filter verifies all
eight parametrizations against the original 17 equations and tests exact
localized lower-band survival.

For a proposed plane Newton normal form, the useful prefilter is now:

1. compile its top bands;
2. compare them with the support box above;
3. replace the full reduced 17-equation locus by the tangent component and
   the two explicit degenerate components above;
4. retain the nonreduced structure explicitly if infinitesimal or
   multiplicity data matter;
5. pass only genuinely new lower bands to Gröbner or de Rham elimination.

This does not solve `JC(2)`: the equation has the intentionally vanishing
Jacobian `X^2`, and the theorem covers one band architecture.  It does
eliminate this box as a source of unclassified reduced Poisson-square
models—the weighted tangent pair is forced off the two degenerate
components.
