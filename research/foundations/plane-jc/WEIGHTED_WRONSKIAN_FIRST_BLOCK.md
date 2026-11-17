# The weighted-Wronskian structure of the first plane block

The first equation in the audited `(72,108)` plane reduction is

\[
 2AD'-3A'D=t^2,
 \qquad
 A=t+a_2t^2+\cdots+a_7t^7+t^8,
 \qquad
 D=d_2t^2+\cdots+d_{12}t^{12}.                 \tag{1}
\]

The archived calculation solves the eleven coefficients of `D` triangularly
and obtains six equations in `a_2,...,a_7`.  This note gives a conceptual
description of that linear stage.  It also shows that the archived degree-35
eliminant is naturally a degree-5 equation on the residual-scaling quotient.

## 1. Weighted Euler form

Put

\[
 A=tU,\qquad D=t^2V,
\]

where `deg(U)=7`, `U(0)=1`, and `deg(V)<=10`.  Equation (1) becomes

\[
 \boxed{T_U(V):=UV+2tUV'-3tU'V=1.}             \tag{2}
\]

Modulo `U`,

\[
 -3tU'V\equiv1\pmod U.                         \tag{3}
\]

Thus the values of `V` at the roots of `U` are interpolation data, and the
degree-`<7` representative is the modular inverse

\[
 V_0=(-3tU')^{-1}\pmod U.                       \tag{4}
\]

No extra squarefree branch is being discarded.  If `U(beta)=U'(beta)=0`,
then (2) evaluated at `beta` gives `0=1`.  Hence every solution of (2) has
`gcd(U,U')=1`, and (4) exists on the entire solution locus.

Write

\[
 V=V_0+US,\qquad \deg S\le3.
\]

The remaining equation is

\[
 \boxed{
 US+2tUS'-tU'S=
 {1-T_U(V_0)\over U}.
 }                                                   \tag{5}
\]

The correction operator on the left maps `k[t]_{<=3}` to `k[t]_{<=9}`.
The apparent coefficient in degree ten is

\[
 1+2\cdot3-7=0.
\]

It has generic rank four.  Its cokernel therefore has dimension `10-4=6`,
which explains the six residual equations without expanding the eleven
coefficients of `D`.

An exact subresultant implementation of the first block can consequently use:

1. the extended resultant of `U` and `-3tU'` to construct `V_0`;
2. a `10 x 4` coefficient matrix for (5);
3. six left-kernel pairings, or equivalently the rank condition on the
   augmented `10 x 5` matrix.

Clearing the resultant denominator gives polynomial compatibility equations.
Because the solution locus is automatically squarefree, saturation by the
discriminant does not remove a genuine solution.

## 2. Hyperelliptic Hermite reduction

Let

\[
 C_A:\quad y^2=A(t).
\]

For squarefree `A` of degree eight, `C_A` is a genus-three hyperelliptic
curve.  A direct calculation gives

\[
 d\left({D\over y^3}\right)
 = {2AD'-3A'D\over2y^5}\,dt.                    \tag{6}
\]

Thus (1) asks whether the second-kind differential

\[
 \omega_A={t^2\,dt\over2y^5}                    \tag{7}
\]

is exact with a primitive having the prescribed behavior at infinity.
Equation (4) is precisely the first Hermite-reduction step: it fixes the
principal parts at the finite branch points.  Equation (5) removes the next
pole layer.  What remains is the compact-curve de Rham class.

The ordinary residues of (7) at the branch points vanish automatically, so
ordinary residue tests alone do not give the compatibility equations.  The
correct residue formulation is the second-kind/de Rham one.  Its obstruction
space has dimension

\[
 \dim H^1_{\mathrm{dR}}(C_A)=2g=6.              \tag{8}
\]

Equivalently, a standard even-degree hyperelliptic Hermite reduction produces
a remainder `R(t)dt/y` with `deg(R)<=6`.  Its two residues at infinity are
opposite, and the second-kind condition kills their one-dimensional
logarithmic component.  The remaining six coordinates are exactly the six
compatibility conditions.

This six-dimensional space has a concrete basis.  Write

\[
 R=r_0+r_1t+\cdots+r_6t^6.
\]

Expanding at either point over infinity gives the second-kind relation

\[
\begin{aligned}
0={}&r_3-{a_7\over2}r_4
 +\left({3a_7^2\over8}-{a_6\over2}\right)r_5\\
&+\left(-{a_5\over2}+{3a_7a_6\over4}-{5a_7^3\over16}\right)r_6. \tag{9}
\end{aligned}
\]

Consequently one convenient compact-curve de Rham basis is represented by

\[
\begin{gathered}
{dt\over y},\quad {t\,dt\over y},\quad {t^2\,dt\over y},\\
\left(t^4+{a_7\over2}t^3\right){dt\over y},\\
\left(t^5+\left({a_6\over2}-{3a_7^2\over8}\right)t^3\right){dt\over y},\\
\left(t^6+\left({a_5\over2}-{3a_7a_6\over4}+{5a_7^3\over16}\right)t^3\right){dt\over y}.
\end{gathered}                                      \tag{10}
\]

The omitted `t^3dt/y` direction is the logarithmic class measuring the two
opposite residues at infinity.  Formula (10) turns the abstract dimension
count into six canonical obstruction coordinates, up to the displayed choice
of infinity coordinate.

## 3. Binary transvectant

Homogenize `A,D` to binary forms `mathcal A,mathcal D` of degrees eight and
twelve.  With the displayed Jacobian convention,

\[
 \mathcal A_Y\mathcal D_X-\mathcal A_X\mathcal D_Y
 =8AD'-12A'D
 =4(2AD'-3A'D).                                  \tag{11}
\]

Hence the first block is a restricted first-transvectant equation.  The
ratio `8:12=2:3` is exactly what cancels the nominal top coefficient.  This is
the representation-theoretic bridge to the relative `(2,3)` scaling in the
quadratic-cubic factorization space.  It does not identify the two moduli
problems, but it identifies their weight character.

## 4. The degree-35 field is a sevenfold cover of a quintic field

The normalization retains the action

\[
 A(t)\longmapsto\zeta^{-1}A(\zeta t),\qquad
 D(t)\longmapsto\zeta^{-2}D(\zeta t),
 \qquad \zeta^7=1.                                \tag{12}
\]

Thus

\[
 a_k\longmapsto\zeta^{k-1}a_k.
\]

On the chart `a_7 != 0`, invariant coordinates are

\[
 q=a_7^7,\qquad x_k=a_ka_7^{k-1}\quad(2\le k\le6). \tag{13}
\]

The preserved lex basis visibly respects this grading.  Its eliminant has
only the powers `35,28,21,14,7,0`, so

\[
 H(a_7)=h(a_7^7)=h(q),\qquad \deg h=5.            \tag{14}
\]

Explicitly,

\[
\begin{aligned}
h(q)={}&9374377445732q^5+62410476400737833472q^4\\
&+265472843532245531128968765q^3\\
&+591414847960503971284831143987840q^2\\
&+586529490054134032292876680565455306752q\\
&-1888043347611739526396142670327809715470336.
\end{aligned}                                      \tag{15}
\]

After multiplying the other five lex relations by the appropriate power of
`a_7`, they express `x_2,...,x_6` polynomially in `q`.  The quotient
first-block scheme is therefore a graph over the quintic scheme `h(q)=0`.
The archived irreducible degree-35 field has the tower

\[
 \mathbb Q\subset\mathbb Q(q)\subset\mathbb Q(a_7),
 \qquad [\mathbb Q(q):\mathbb Q]=5,
 \qquad [\mathbb Q(a_7):\mathbb Q(q)]=7.          \tag{16}
\]

This is the most economical interpretation of the degree 35: five quotient
solutions, each carrying the residual sevenfold scaling orbit.

The quotient does not conceal a radical parametrization.  Exact computation
shows that `h` is irreducible over `Q` and has Galois group `S_5`.  Thus the
arithmetic separates into a genuinely nonsolvable quintic quotient followed
by the degree-seven Kummer step `a_7^7=q`.

## 5. Reproduction

Run

```bash
python3 plane-jc/cas/weighted_wronskian_first_block.py
```

in an environment providing SymPy.  The checker verifies (2), the generic
rank-four correction operator and its six-dimensional cokernel, (9), the
`Z/7` grading of all six archived lex relations, and the quotient-quintic
form (14), including irreducibility and Galois group `S_5`.  It reads the
preserved exact lex output and does not rerun the large Groebner-basis
calculation.
