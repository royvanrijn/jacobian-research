# Endpoint reduction for the cancellation contact resultant

For positive integers `m,r`, put

\[
 p(v)=v(1-v)^m,
 \quad
 K_{m,r}(w)=w^{-r-1}\int_0^w p(v)^r\,dv,
\]

and

\[
 L_{m,r}(w)=w^{-r-1}\int_0^w\{p(w)-p(v)\}^r\,dv.
\]

The full cancellation prime-intersection diagram requires
`Res_w(K_{m,r},L_{m,r}) != 0`.  This note records a uniform endpoint
reduction and proves the four complete columns `r=1,2,3,4`.  It does **not**
claim the remaining all-parameter theorem for `r>=5`.

## 1. Triangular beta-primitive identity

Define `M_k(w)=K_{m,k}(w)`, including `M_0=1`.  Scaling `v=wt` and expanding
the binomial gives

\[
 \boxed{
 L_{m,r}(w)=\sum_{k=0}^r(-1)^k\binom rk
 (1-w)^{m(r-k)}M_k(w).}                               \tag{1}
\]

The last summand is `(-1)^r K_{m,r}`.  Thus modulo `K_{m,r}`, the contact
coefficient is a triangular combination of lower incomplete-beta
primitives.  Also

\[
 K_{m,r}(w)=\frac1{r+1}
 {}_2F_1(-mr,r+1;r+2;w),                              \tag{2}
\]

so (1) is the useful contiguous structure behind the proposed
hypergeometric search.

## 2. Fixed-`r` endpoint equations

Put

\[
 y=1-w,\qquad z=y^m.
\]

After `u=1-wt`,

\[
 M_k(w)=\frac{1}{(1-y)^{k+1}}
 \left(\beta_k-yz^kT_k(y)\right),                    \tag{3}
\]

where

\[
 \beta_k=\frac{k!}{\prod_{j=1}^{k+1}(mk+j)},\qquad
 T_k(y)=\sum_{j=0}^k\frac{(-1)^j\binom kj}{mk+j+1}y^j.
                                                                    \tag{4}
\]

Equations (1), (3), and (4) turn a degree-`mr` resultant in `w` into two
degree-`r` equations in the auxiliary endpoint variable `z`, together with
the single exponential relation `z=y^m`.  For each fixed `r`, eliminating
`z` has degree independent of the growing degree `mr`.  This is a sharper
symbolic-search formulation than interpolating the raw Sylvester resultant.

The points `y=1` and `y=0` never give a common root: they are respectively
`w=0`, where `K_{m,r}=1/(r+1)`, and `w=1`, where
`K_{m,r}=Beta(r+1,mr+1)`.

There is also a uniform root disk which will be decisive for `r=3`.  Put
`n=mr`.  The complementary incomplete-beta identity gives

\[
 K_{m,r}(1-y)=\beta_r B_{n,r}(y),\qquad
 B_{n,r}(y)=\sum_{j=0}^{n}\binom{j+r}{r}y^j.          \tag{5}
\]

The coefficients are positive, and their consecutive ratios are

\[
 \frac{\binom{j+r-1}{r}}{\binom{j+r}{r}}=\frac{j}{j+r}.
\]

Enestrom--Kakeya therefore puts every zero of `K_{m,r}(1-y)` in

\[
 \boxed{|y|\le \frac{n}{n+r}=\frac{m}{m+1}.}          \tag{6}
\]

## 3. The complete column `r=1`

Equation (1) becomes

\[
 L_{m,1}+K_{m,1}=(1-w)^m.                             \tag{7}
\]

A common root would therefore be `w=1`, but

\[
 K_{m,1}(1)=B(2,m+1)=\frac1{(m+1)(m+2)}\ne0.
\]

More precisely,

\[
 \boxed{
 \operatorname{Res}_w(K_{m,1},L_{m,1})
 =\big((m+1)(m+2)\big)^{-m}.}                         \tag{8}
\]

This also gives a literal Bézout certificate.  Let `U_m(w)` be the truncation
through degree `m-1` in `1-w` of the Taylor series `1/K_{m,1}(w)` at `w=1`.
Then `U_mK_{m,1}-1` is divisible by `(1-w)^m=K_{m,1}+L_{m,1}`.  Dividing it
produces explicit polynomials `A_m,B_m` with

\[
 A_mK_{m,1}+B_mL_{m,1}=1.                             \tag{9}
\]

## 4. The complete column `r=2`

At a common root, `M_2=K_{m,2}=0`.  Equation (1) reads

\[
 L_{m,2}=z^2-2zM_1+M_2.
\]

Neither `y` nor `z=y^m` is zero, so `L_{m,2}=0` forces `2M_1=z`.
Using (3)--(4), this is

\[
 zA_m(y)=2,                                           \tag{10}
\]

with

\[
 A_m(y)=m(m+1)y^2-2m(m+2)y+(m+1)(m+2).
\]

Substitute `z=2/A_m(y)` into `M_2=0` and clear the displayed nonzero
denominators.  Direct simplification gives

\[
 (m+1)^2(y-1)^3\bigl(m^2y-(m+2)^2\bigr)=0.            \tag{11}
\]

The first candidate was excluded above.  The second is

\[
 y=\frac{(m+2)^2}{m^2}>1.
\]

At this value

\[
 A_m(y)=\frac{2(m+2)(5m^2+10m+4)}{m^3}>2,
\]

while `z=y^m>1`; hence `zA_m(y)>2`, contradicting (10).  Therefore

\[
 \boxed{\operatorname{Res}_w(K_{m,2},L_{m,2})\ne0
 \quad\text{for every }m\ge1.}                       \tag{12}
\]

## 5. The complete column `r=3`

At a common root, `M_3=0`, and (1) gives

\[
 L_{m,3}=z^3-3z^2M_1+3zM_2-M_3.
\]

Since `z` is nonzero, `L_{m,3}=0` reduces to

\[
 z^2-3zM_1+3M_2=0.                                  \tag{13}
\]

Put `D=1-y`.  After inserting (3), clearing `D^3` in (13), and clearing
`D^4` in `M_3=0`, the two equations have the sparse form

\[
 az^2+bz+c=0,\qquad dz^3+e=0,                         \tag{14}
\]

where

\[
 \begin{aligned}
 a&=D^3+3yDT_1-3yT_2,& b&=-3\beta_1D,& c&=3\beta_2,\\
 d&=-yT_3,&& e=\beta_3.
 \end{aligned}                                       \tag{15}
\]

The quadratic--binomial-cubic resultant is

\[
 \operatorname{Res}_z(az^2+bz+c,dz^3+e)
 =a^3e^2+3abcde-b^3de+c^3d^2.                        \tag{16}
\]

After substitution from (15), (16) factors as the excluded endpoint
`(y-1)^3` times a degree-six polynomial `H_m(y)` in `Q[m][y]`.  There is no
need to factor this sextic or to solve the additional equation `z=y^m`.
Instead put

\[
 \rho=\frac{m}{m+1},\qquad Q_m(u)=u^6H_m(\rho/u).     \tag{17}
\]

Write
`Q_m(u)=q_0u^6+q_1u^5+...+q_6`.  Let `A_m` and `B_m` be the lower triangular
Toeplitz matrices

\[
 (A_m)_{ij}=q_{i-j},\qquad (B_m)_{ij}=q_{6-i+j}
 \quad(0\le j\le i\le5).                             \tag{18}
\]

The Schur--Cohn criterion says that every zero of `Q_m` lies in the open
unit disk precisely when

\[
 C_m=A_mA_m^{\mathsf T}-B_mB_m^{\mathsf T}
\]

is positive definite.  After clearing a harmless common denominator from
`Q_m`, the six leading principal minors of `C_m` have degrees

\[
 29,\ 56,\ 81,\ 104,\ 125,\ 144                     \tag{19}
\]

in `m`.  Every coefficient of all six minors is strictly positive.  Thus
Sylvester's criterion gives `C_m>0` for every positive `m`, and hence

\[
 H_m(y)=0\quad\Longrightarrow\quad
 \left|\frac{\rho}{y}\right|<1
 \quad\Longrightarrow\quad |y|>\frac{m}{m+1}.        \tag{20}
\]

This is disjoint from the `K_{m,3}` root disk (6).  A common root of
`K_{m,3}` and `L_{m,3}` is therefore impossible, proving

\[
 \boxed{\operatorname{Res}_w(K_{m,3},L_{m,3})\ne0
 \quad\text{for every }m\ge1.}                       \tag{21}
\]

## 6. The complete column `r=4`

At a common root, `M_4=0`, and division of (1) by the nonzero `z` gives

\[
 z^3-4z^2M_1+6zM_2-4M_3=0.                         \tag{22}
\]

Clearing `D^4` in (22) and `D^5` in `M_4=0` produces a cubic `E_m(y,z)`
and a quartic `F_m(y,z)` in `z`.  Their subresultant chain has `z`-degrees

\[
 4,3,2,1,0.
\]

The last member factors as

\[
 \operatorname{Res}_z(E_m,F_m)=(y-1)^5H_m(y),       \tag{23}
\]

up to a nonzero rational function of `m`, where `H_m` has degree eleven.
The penultimate member is linear:

\[
 A_m(y)z+B_m(y).                                    \tag{24}
\]

For `1<=m<=21`, direct substitution `z=y^m` gives the exact polynomial gcd

\[
 \gcd\big(E_m(y,y^m),F_m(y,y^m)\big)=(y-1)^5.       \tag{25}
\]

The only candidate in (25) is the excluded endpoint `y=1`, where
`K_{m,4}(0)=1/5`.

For the uniform tail `m>=22`, put `rho=m/(m+1)` and

\[
 Q_m(u)=u^{11}H_m(\rho/u).
\]

The eleven leading minors of the Schur--Cohn form of `Q_m` have degrees

\[
 65,128,189,248,305,360,413,464,513,560,605.        \tag{26}
\]

After `m=n+22`, every coefficient in these minors has the sign

\[
 +,+,-,+,+,+,+,+,+,+,+.                            \tag{27}
\]

Fraction-free `LDL^T` elimination therefore gives nine positive and two
negative pivots.  By the inertia form of Schur--Cohn, `Q_m` has nine roots
inside the unit disk and two outside.  Equivalently, `H_m` has exactly two
roots in `|y|<rho` and nine outside it.  Thus the pure modulus separation
used for `r=3` fails in exactly one conjugate pair.

That pair can nevertheless be localized uniformly.  Set

\[
 t=1/m,\qquad x=m(y-1),\qquad
 c=-\frac{523}{200}+\frac{41}{8}i,qquad R=\frac12. \tag{28}
\]

Then `t^14 H_{1/t}(1+tx)` is a polynomial of degree eleven in `x` and degree
twenty-four in `t`.  A rational Rouche certificate on the 228 intervals of
width at most `1/5000` covering `0<=t<=1/22` proves that it has exactly one
root in `|x-c|<R`.  The conjugate disk contains the other root.  Each cell
uses a rational degree-eleven comparison polynomial; exact root-distance,
coefficient-error, and parameter-derivative bounds prove the strict Rouche
inequality.  Numerical roots only propose those rational comparison
polynomials and are not used in the certificate.

The upper disk lies strictly inside the transformed `K`-disk.  Indeed, after
squaring the triangle-inequality bound, its numerator is

\[
 44600m^3-487877m^2-1209554m-657077,
\]

whose coefficients become strictly positive after `m=n+22`.  Moreover every
point in the upper disk satisfies

\[
 \pi<m\arg(y)<2\pi.                                 \tag{29}
\]

The lower inequality follows from
`atan(q)>q/(1+q^2)` and `pi<22/7`; the resulting shifted quadratic has
coefficients `415000,19243475,199736262`.  For the upper inequality, compare
the disk with the ray of angle `2pi/m`, use the alternating cubic/quartic
Taylor bounds for sine and cosine and
`333/106<pi<355/113`.  The shifted quartic has coefficients

\[
 361647199066050, 22794972938293518,
 487119829014994663, 3797665752647679636,
 5602580029456243464,
\]

all positive.  Consequently `y^m` lies in the lower half-plane.

It remains to compare the actual endpoint variable `z`.  On a root of
`E_m,F_m`, (24) gives `z=-B_m/A_m` provided `A_m` is nonzero.  Substitute
`y=1+tx`, write this rational function as `N/D`, and put `x=a+ib`.  The
identity

\[
 S(t,a,b)=\operatorname{Im}(N\overline D)=|D|^2\operatorname{Im}(z)
\]

is certified positive on the rectangular box containing the disk:

\[
 0\le t\le1/22,\quad -623/200\le a\le-423/200,\quad
 37/8\le b\le45/8.                                  \tag{30}
\]

After mapping (30) to the unit cube, every tensor-product Bernstein
coefficient of `-S/(4b(3t+2))` is strictly negative.  Hence `S>0`; this also
proves `D!=0`.  Thus the recovered `z` lies in the upper half-plane, whereas
(29) puts `y^m` in the lower half-plane.  The conjugate disk gives the
conjugate contradiction.  Therefore

\[
 \boxed{\operatorname{Res}_w(K_{m,4},L_{m,4})\ne0
 \quad\text{for every }m\ge1.}                      \tag{31}
\]

## 7. Remaining step

The all-parameter problem now begins at `r=5`.  Formula (6) still gives a
fixed comparison disk, but the `r=4` analysis shows that demanding every
endpoint-eliminant root lie outside it is too strong.  The reusable strategy
is instead: compute Schur--Cohn inertia, localize the exceptional branches,
and separate `z` from `y^m` by modulus or argument.  Subresultant recurrences
may keep the fixed-`r` eliminants and their Schur transforms smaller than a
direct expansion.

Raw interpolation of `Res_w(K,L)` still obscures this structure because its
degree grows with `mr`.

The general reduction and `r<=3` certificate are in
[`verify_contact_resultant_endpoint_reduction.py`](../scripts/verify_contact_resultant_endpoint_reduction.py).
The exact inertia, Rouche, angle, and Bernstein certificates for `r=4` are in
[`verify_contact_resultant_r4.py`](../scripts/verify_contact_resultant_r4.py).
