# Complete scalar degree-nine closure for `HC4`

## Status and scope

This note continues
[`HC4_HIGHER_DEGREE_PENCIL_OBSTRUCTIONS.md`](HC4_HIGHER_DEGREE_PENCIL_OBSTRUCTIONS.md)
and the degree-eight theorem `HC4RSD41`.

> **Theorem HC4RSD43 — complete scalar degree-nine closure.**  In the
> synchronized scalar reverse-Schur packet of `HC4RSD20`, every border
> coefficient whose leading binary form has degree nine has a fixed ruling.
> Hence every such packet reduces to `HC2` or to the exact `JC2` cotangent
> endpoint.

Consequently the scalar reverse-Schur branch is closed through leading degree
nine.  This is not unrestricted `HC4`: degree ten and above, polynomially
moving matrix flags, non-scalar/coisotropic pivots, direct four-variable
constructions, and the `JC2` endpoint remain open.

## 1. Local root values

Let a root of the degree-\(d\) binary top \(f\) have multiplicity \(m\), and
let the first transverse coefficient \(g\) have degree \(e\).  Away from a
resonance, put \(n=\lceil m/2\rceil\).  Comparing the first nonzero local
terms in

\[
 q\det B_f=b_g^{\mathsf T}\operatorname{adj}(B_f)b_g
\]

gives the root value

\[
 \kappa_{d,e}(m)
 =-\frac{C_{d,m,e,n}}{dm(d-m)},                     \tag{1.1}
\]

where \(C_{d,m,e,n}\) is the resonance polynomial from `HC4RSD27`.

Write

\[
 f=A^2B,
\]

where \(B\) is the product of the roots having odd multiplicity.  Then

\[
 g=ABH,
 \qquad q=BC,
\]

and (1.1) becomes the interpolation rule

\[
 C(r)=\kappa_{d,e}(m_r)H(r)^2                     \tag{1.2}
\]

at every distinct root \(r\).

For \(d=9\), the decisive rows are

\[
\begin{array}{c|rrrrrrrr}
 m&1&2&3&4&5&6&7&8\\ \hline
 \kappa_{9,6}(m)&23/72&47/126&4/9&11/45&71/180&-1/18&10/63&-14/9\\
 \kappa_{9,7}(m)&11/18&59/126&11/18&1/9&14/45&-13/18&-85/126&-44/9.
\end{array}                                         \tag{1.3}
\]

The equality

\[
 \kappa_{9,7}(1)=\kappa_{9,7}(3)=11/18             \tag{1.4}
\]

explains every nonresonant exceptional ray that survives the local filter.

## 2. Two roots, three roots, and the pure power

The complete two-root partitions

\[
 8+1,\quad7+2,\quad6+3,\quad5+4
\]

are closed for every admissible \(e=5,6,7\).  The complete three-root
partitions are also closed.  Four isolated \(e=7\) rays have standard basis

\[
 (a^3,r_0,r_1,r_2,r_3,s_0,s_1).
\]

The resonant `(3,3,3)` rows have respectively

\[
 (u^3,uv(u+v),v^3)
\]

at \(e=5\), and

\[
 (a^2r,\ a^3-81r/2,\ r^2)
\]

at \(e=6\).

For the pure ninth power, every coefficient involving a second passive
linear form is removed successively by an explicit square coefficient.  The
remainder depends on only two linear forms and is a fixed cylinder.

## 3. Three odd multiplicities

At \(e=6\), both \(H\) and \(C\) are constant.  Every relevant partition
contains two unequal values in the first row of (1.3), so \(H=0\).

At \(e=7\), interpolation closes three partitions immediately.  Three
one-ray packets reach the complete face:

\[
 (6,1,1,1),\qquad(4,3,1,1),\qquad(3,3,2,1).
\]

Over the corresponding cross-ratio function field, each has exact basis

\[
 (a^2,r_0,r_1,r_2,r_3,s_0,s_1).                    \tag{3.1}
\]

The remaining `(5,2,1,1)` projective Schur system is the unit ideal on both
charts.  Thus every three-odd packet closes.

## 4. Five odd multiplicities

Here \(H\) and \(C\) are constant at \(e=7\).  Equation (1.3) eliminates
every partition except

\[
 (3,3,1,1,1),
\]

where all root values equal \(11/18\).  Normalize five roots to
\(0,\infty,1,\lambda,\mu\), put

\[
 B=xy(x-y)(x-\lambda y)(x-\mu y),
\]

and use

\[
 f=x^2y^2B,
 \qquad g=a xyB,
 \qquad q=\frac{11}{18}a^2B.
\]

The complete weight-nine face again has basis (3.1).

## 5. One odd multiplicity

The remaining nonresonant root partitions are

\[
 (4,2,2,1),\qquad(3,2,2,2),\qquad(2,2,2,2,1).
\]

The `(4,2,2,1)` projective Schur cover is exact.  At \(e=6\), both charts
are unit ideals.  At \(e=7\), the leading projective chart is supported only
at \(\lambda=0,1\), where roots collide, and the other two charts are units.

For the last two partitions, send the odd root to infinity.  The affine
problem becomes

\[
 f=A(t)^2,
 \qquad g=A(t)H(t).
\]

At \(e=6\), interpolation gives \(C=47H^2/126\) and the identity

\[
 21A(H')^2+3H^2A''-14HA'H'=0.                     \tag{5.1}
\]

After an affine normalization, \(H=1\) or \(H=t\).  The first forces
\(A''=0\).  On \(H=t\), the coefficient of \(t^j\) is

\[
 3j^2-17j+21,
\]

which is nonzero for \(0\le j\le4\).

At \(e=7\), root interpolation leaves

\[
 126ALA''-98L(A')^2-63A(H')^2-25H^2A''+70HA'H'=0. \tag{5.2}
\]

For `(3,2,2,2)`, \(A\) is cubic and \(L\) is linear.  For
`(2,2,2,2,1)`, \(A\) is quartic and \(L\) is constant.  Classify the
quadratic \(H\) under affine changes as

\[
 1,\qquad t,\qquad t^2,\qquad t(t-1).
\]

The first, second, and fourth charts have unit coefficient ideals.  The
third leaves only

\[
 A=t^3,\ L=t/7
 \qquad\text{or}\qquad
 A=t^4,\ L=1/7,
\]

which collapses all finite roots and is outside the declared squarefree
packet.  This closes the final nonresonant row.

## 6. The exceptional `(m,e,n)=(3,5,1)` resonance

The only proper-root resonance below \(\lceil m/2\rceil\) in degree nine is

\[
 (m,e,n)=(3,5,1).
\]

The `(4,3,1,1)`, `(3,2,2,1,1)`, and `(3,2,2,2)` Schur systems are unit
ideals.  The projective `(3,3,2,1)` system is supported only on the root
collision \(\lambda=1\).

The sole distinct-root Schur survivor is `(3,3,1,1,1)`, with

\[
 \lambda+\mu+1=0,
 \qquad \mu^2+\mu+1=0.                             \tag{6.1}
\]

No higher same-weight tail is available at \(e=5\).  In the quadratic
residue field, the complete bordered polynomial contains

\[
 [x^{16}y^4z]J(c)=36a^3.                            \tag{6.2}
\]

Thus the equianharmonic survivor is empty as well.

Sections 2--6 exhaust every partition of nine and prove `HC4RSD43`.

## 7. Reproduction

Run the independent exact checkers

```bash
.venv/bin/python scripts/verify_hc4_degree_nine_two_root.py
.venv/bin/python scripts/verify_hc4_degree_nine_three_root_survivors.py
.venv/bin/python scripts/verify_hc4_degree_nine_pure_power.py
.venv/bin/python scripts/verify_hc4_degree_nine_three_odd.py
.venv/bin/python scripts/verify_hc4_degree_nine_five_odd.py
.venv/bin/python scripts/verify_hc4_degree_nine_one_odd.py
.venv/bin/python scripts/verify_hc4_degree_nine_resonance.py
```

They write the corresponding JSON records under
`artifacts/generated-results/`.  Every promoted result is over
characteristic zero; no modular or numerical survivor is used.
