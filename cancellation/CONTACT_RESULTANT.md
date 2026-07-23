# Endpoint reduction for the cancellation contact resultant

This note supplies refined projective-branch intersection geometry.  It is
not a dependency of the degreewise multiplicity theorem, whose nilpotency
index follows unconditionally from the two proofs in
the unconditional thick-contact theorem retained in the
[degreewise audit](../DEGREEWISE_MULTIPLICITY_AUDIT.md).
The concise paper appendix is
the refined boundary calculations retained below.

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
reduction and proves the seven complete columns `r=1,2,3,4,5,6,7`.  It does
**not** claim the remaining all-parameter theorem for `r>=8`.

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

## 7. The complete column `r=5`

At a common root, division of (1) by the nonzero `z` gives

\[
 z^4-5z^3M_1+10z^2M_2-10zM_3+5M_4=0.
\]

After endpoint denominators are cleared, this is a quartic
`E_m(y,z)`, while `M_5=0` is a binomial quintic `F_m(y,z)`.  Their
subresultant chain has `z`-degrees

\[
 5,4,3,2,1,0.
\]

The sparse 26-term quartic--binomial-quintic resultant factors as

\[
 \operatorname{Res}_z(E_m,F_m)=(y-1)^5H_m(y),       \tag{32}
\]

up to a nonzero rational function of `m`; the primitive eliminant `H_m` has
bidegree `(49,20)` in `(m,y)`.  The penultimate subresultant is again linear
in `z`, and its two coefficients have explicit sparse formulas.  Direct exact
substitution `z=y^m` gives

\[
 \gcd(E_m(y,y^m),F_m(y,y^m))=(y-1)^5
 \qquad(1\le m\le19).
\]

This leaves only the excluded endpoint `y=1`, where `K_{m,5}(0)=1/6`.

For `m>=20`, put `t=1/m`, `x=m(y-1)`, and

\[
 P_5(t,x)=t^{29}H_{1/t}(1+tx).
\]

This is a polynomial of degree twenty in `x`, degree forty-nine in `t`, with
840 terms.  The exact Routh table of `P_5(0,-4/5-s)` has two sign changes.
On the line `Re(x)=-4/5`, write

\[
 P_5(t,-4/5+i\omega)=R(t,\omega^2)+i\omega S(t,\omega^2).
\]

The crossing resultant `Res_{omega^2}(R,S)` has degree 931 in `t`.  After
`t=T/(20(1+T))`, all 932 coefficients are strictly negative; the omitted
zero-frequency factor has its own one-sign certificate.  Thus exactly two
roots remain to the left of this line for `0<=t<=1/20`.  The `K`-disk is
contained in that half-plane because

\[
 m-\frac45-\frac{m^2}{m+1}=\frac{m-4}{5(m+1)}>0.
\]

The upper exceptional member remains in the fixed rectangle

\[
 -\frac{13}{4}<\operatorname{Re}x<-\frac{29}{10},
 \qquad
 \frac95<\operatorname{Im}x<\frac{11}{5}.           \tag{33}
\]

At `t=0`, an exact complex root count gives one root in (33).  A root could
leave only through one of its four boundary lines.  Each line-crossing
resultant has degree 1911 in `t`, and its transformed coefficients have one
strict sign.  The whole rectangle is inside the `K`-disk: after `m=n+20`,
the numerator of the squared-radius difference has coefficients
`380,22535,443330,2891275`.

The binomial equation gives the branch-independent identity

\[
 z^5=\frac{\beta_5}{yT_5(y)}.
\]

Writing the right side as `N/D` after the `(t,x)` substitution, every
tensor-product Bernstein coefficient of

\[
 9^{10}|N|^2-|D|^2
\]

is strictly positive on the full rectangle and `0<=t<=1/20`.  Thus every
possible `F_m` branch has `|z|>1/9`, including a hypothetical degeneration
of the linear subresultant.  Meanwhile

\[
 \log|y^m|
 \le -\frac{29}{10}
 +\frac1{40}\left(\frac{169}{16}+\frac{121}{25}\right)
 <-\frac52.
\]

The degree-seven Taylor lower bound already gives `exp(5/2)>12`, so
`|y^m|<1/12`.  Therefore `z` cannot equal `y^m`.
Complex conjugation handles the lower member, while all other eliminant roots
are outside the `K`-disk.  Consequently

\[
 \boxed{\operatorname{Res}_w(K_{m,5},L_{m,5})\ne0
 \quad\text{for every }m\ge1.}                      \tag{34}
\]

## 8. A contiguous recurrence in `r`

The endpoint tails in (4) have the moment representation

\[
 T_r(y)=\int_0^1x^{mr}(1-yx)^r\,dx.                 \tag{35}
\]

There is a first-order, inhomogeneous recurrence in `r`.  For
`0<=j<=m`, set

\[
 a_j=mr+j+1,\qquad b_j=(m+1)r+j+2,
\]

and use the empty-product convention.  Define

\[
 A_{m,r}=\prod_{j=0}^{m-1}a_j,
 \qquad B_{i;m,r}=\prod_{j=0}^{i-1}b_j,             \tag{36}
\]

and

\[
 P_{m,r}(y)=y^mB_{m;m,r}
 -(r+1)\sum_{i=0}^{m-1}y^iB_{i;m,r}
                  \prod_{j=i+1}^{m-1}a_j.           \tag{37}
\]

Then

\[
 \boxed{
 y^mB_{m+1;m,r}T_{r+1}(y)
 =(r+1)A_{m,r}T_r(y)+(1-y)^{r+1}P_{m,r}(y).}         \tag{38}
\]

Indeed, if

\[
 U_j=\int_0^1x^{mr+j}(1-yx)^r\,dx,
\]

integration of the derivative of
`x^{mr+j+1}(1-yx)^{r+1}` gives

\[
 a_jU_j-yb_jU_{j+1}=(1-y)^{r+1}.                   \tag{39}
\]

Iterating (39) through `j=m` and using
`T_{r+1}=U_m-yU_{m+1}` yields (38).  At `y=1`, it specializes to the
first-order hypergeometric recurrence

\[
 B_{m+1;m,r}\beta_{r+1}=(r+1)A_{m,r}\beta_r.       \tag{40}
\]

On the exponential surface `z=y^m`, the cleared moment equation
`F_r=\beta_r-yz^rT_r=(1-y)^{r+1}M_r` therefore obeys

\[
B_{m+1;m,r}F_{r+1}
 =(r+1)A_{m,r}F_r-yz^r(1-y)^{r+1}P_{m,r}(y).        \tag{41}
\]

The contact equation has a companion recurrence.  Normalize

\[
 N_k=M_k/z^k,\qquad
 G_r=L_{m,r}/z^r=\sum_{k=0}^r(-1)^k\binom rkN_k.     \tag{42}
\]

Thus `G` is the alternating binomial transform of `N`.  If `S` is the forward
shift on sequences and

\[
 \Theta=r(1-S^{-1}),
\]

then alternating binomial transform carries multiplication by `k` to
`Theta`, and carries the shift `N_k -> N_{k+1}` to `1-S`.  Put

\[
 C_m(k)=\prod_{j=0}^{m}\bigl((m+1)k+j+2\bigr),
 \qquad
 Q_m(k)=(k+1)\prod_{j=0}^{m-1}(mk+j+1).             \tag{43}
\]

Dividing the moment recurrence corresponding to (41) by `z^k` and taking the
alternating binomial transform gives

\[
 \boxed{
 z(1-y)C_m(\Theta)(1-S)G-Q_m(\Theta)G
 =-y\,\widehat P,}                                  \tag{44}
\]

where

\[
 \widehat P_r=\sum_{k=0}^r(-1)^k\binom rkP_{m,k}(y).
\]

For fixed `m`, (37) is a polynomial in `k` of degree at most `m`.  Hence
`\widehat P_r=0` for every `r>m`, and (44) is then a homogeneous recurrence of
order at most `m+2` for the second endpoint equation.  Equations (41) and
(44) give the requested contiguous system in `r` without constructing a new
Sylvester matrix in each column.

The remaining uniform step is now precise: propagate the polynomial
subresultant chain of the two endpoint equations under (41) and (44), and
prove that its terminal member cannot vanish on `z=y^m`.  The recurrences
alone do not yet establish that nonvanishing.

## 9. Arithmetic all-`r` ranges

There is a second uniform route which does not require an endpoint eliminant.
Both contact polynomials have degree `mr`, and their first coefficients are

\[
 \begin{aligned}
 K_{m,r}(w)&=\frac1{r+1}-\frac{mr}{r+2}w+O(w^2),\\
 L_{m,r}(w)&=\frac1{r+1}
  -\frac{mr(r+3)}{(r+1)(r+2)}w+O(w^2).              \tag{45}
 \end{aligned}
\]

Thus they have the same nonzero constant coefficient but are not associates.
If `K_(m,r)` is irreducible and the resultant vanished, then `K_(m,r)` would
divide `L_(m,r)`.  Equal degrees and constant coefficients would force
`K_(m,r)=L_(m,r)`, contradicting (45).  Hence

\[
 \boxed{K_{m,r}\text{ irreducible over }\mathbb Q
 \quad\Longrightarrow\quad
 \operatorname{Res}_w(K_{m,r},L_{m,r})\ne0.}        \tag{46}
\]

The fractional-linear identity checked in the boundary diagram identifies
`K_(m,r)` with the cancellation parameter polynomial
`P_((m+1)r+1,mr)` up to a nonzero scalar.  The
[two-prime diagonal theorem](DIAGONAL_TWO_PRIME_IRREDUCIBILITY.md) therefore
gives the following contact-resultant theorem immediately:

\[
 \boxed{
 \begin{gathered}
  \operatorname{Res}(K_{m,r},L_{m,r})\ne0
       \quad(1\le m\le1000,\ r\ge1),\\
  \text{and for every }m\ge2,\text{ the same holds whenever }r>X_m/m,
  \\
  \text{and whenever }mr\ge K_0,\quad r\ge5(mr)^{21/40}.
 \end{gathered}}                                    \tag{47}
\]

Here `X_m` is the explicit Dusart threshold (8f), while `K_0` is the
non-numerical Baker--Harman--Pintz threshold, in the linked arithmetic
theorem.

More generally, (46) transfers every proved parameter-irreducibility range,
including the prime and two-prime interval criteria.  This is independent of
the six complete fixed-`r` columns: it supplies 1000 complete all-`r`
columns in the transverse direction and an explicit effective tail in every
fixed-`m` column.

Consequently a still-unresolved pair in `OP-CR` must have `m>=1001`, `r>=8`,
satisfy `r<=X_m/m`, have composite `N=(m+1)r+1`, and
have at most one prime in
`(mr,(m+1)r+1)`, unless another recorded irreducibility criterion applies.
If `mr>=K_0`, it must additionally satisfy `r<5(mr)^(21/40)`.
If that interval contains exactly one prime `N-u`, the only possible factor
degrees of `K` are `u` and `mr-u`; hence any common factor with `L` must have
one of those degrees.  This factor-degree restriction is the useful bridge
between the arithmetic and endpoint approaches.

## 10. The exact first-open-column reduction `r=6`

The bounded-degree endpoint formulation also makes the first column beyond
the five complete fixed-`r` theorems computationally exact.  At a common
root with `r=6`, divide the triangular contact equation by the nonzero
`z=y^m` and clear `(1-y)^6`.  Together with the moment equation this gives

\[
 \begin{aligned}
 E_6(m,y,z)&=\sum_{k=0}^{5}(-1)^k\binom6k
 z^{5-k}\bigl(\beta_k-yz^kT_k(y)\bigr)(1-y)^{5-k}=0,\\
 F_6(m,y,z)&=\beta_6-yz^6T_6(y)=0.                 \tag{48}
 \end{aligned}
\]

Thus `E_6` has degree five and `F_6` degree six in `z`, independently of
`m`.  Clear their parameter denominators over `Q[m,y]` and form the exact
`z`-resultant.  Direct polynomial arithmetic gives

\[
 \boxed{
  \operatorname{Res}_z(E_6,F_6)=(y-1)^7H_6(m,y),
  \qquad \deg_yH_6=29,\quad\deg_mH_6=90,
  \quad H_6(m,1)\ne0.}                             \tag{49}
\]

The factor at `y=1` is the already excluded endpoint.  Hence every genuine
common root in the `r=6` column projects to one of only 29 roots of the
explicit bounded-degree eliminant `H_6`, rather than to a root of a
degree-`6m` Sylvester resultant.  This is a rigorous reduction, not yet the
uniform nonvanishing theorem: one must reconstruct the corresponding `z`
and rule out `z=y^m` on all 29 branches.  The exact construction and degree
certificate are replayed by
[`verify_contact_resultant_r6_reduction.py`](../scripts/verify_contact_resultant_r6_reduction.py).

## 11. The 29 branches at `m=infinity`

The reduction (49) changes the first open column into a fixed-degree
algebraic-function problem.  Put `t=1/m` and use the boundary-layer scale

\[
 y=1+ct.
\]

For `k>=1`, define

\[
 b_k=\frac{k!}{k^{k+1}},\qquad
 U_k(c)=\sum_{j=0}^k\binom kj(-c)^{k-j}\frac{j!}{k^{j+1}}
       =\int_0^\infty e^{-ks}(s-c)^k\,ds.             \tag{50}
\]

These are exactly the limits

\[
 m^{k+1}\beta_k\longrightarrow b_k,\qquad
 m^{k+1}T_k(1+c/m)\longrightarrow U_k(c).
\]

Consequently `m^6 E_6` and `m^7 F_6` tend to the fixed equations

\[
 \begin{aligned}
 e_6(c,z)&=c^6z^5+\sum_{k=1}^5(-1)^k\binom6k
 z^{5-k}\bigl(b_k-z^kU_k(c)\bigr)(-c)^{5-k},\\
 f_6(c,z)&=b_6-z^6U_6(c).                            \tag{51}
 \end{aligned}
\]

Exact elimination gives

\[
 \boxed{\operatorname{Res}_z(e_6,f_6)=\lambda c^7P_{29}(c)},\qquad
 \lambda\in\mathbb Q^\*,                            \tag{52}
\]

where `P_29` is squarefree of degree 29, `P_29(0)!=0`, and
`gcd(P_29,U_6)=1`.  The last equality excludes a limiting solution at
`z=infinity`.  This is not merely a formal limit of the endpoint equations.
If

\[
 H_6(m,1+x)=\sum h_{a,b}x^am^b,
\]

then `b-a<=61` for every nonzero coefficient, and the complete top edge is

\[
 \sum_{b-a=61}h_{a,b}c^a=\lambda'P_{29}(c),
 \qquad \lambda'\in\mathbb Q^\*.                    \tag{53}
\]

Thus all 29 branches, counted with multiplicity, are ordinary analytic
branches at `t=0`:

\[
 \boxed{y_i(m)=1+\frac{c_i}{m}+O(m^{-2}),\qquad
        P_{29}(c_i)=0.}                              \tag{54}
\]

The penultimate subresultant of (51) is `A(c)z+B(c)`, and exact arithmetic
gives `gcd(P_29,A)=1`.  It therefore reconstructs one algebraic value

\[
 d_i=-B(c_i)/A(c_i),\qquad z_i(m)=d_i+O(m^{-1})       \tag{55}
\]

on every branch.  On the other hand,

\[
 y_i(m)^m\longrightarrow e^{c_i}.                   \tag{56}
\]

Every `c_i` is nonzero algebraic, so Lindemann--Weierstrass makes `e^(c_i)`
transcendental, whereas `d_i` is algebraic.  Hence `d_i!=e^(c_i)` on all
29 branches.  Finiteness of the branch set now gives

\[
 \boxed{\operatorname{Res}_w(K_{m,6},L_{m,6})\ne0
        \quad\hbox{for every sufficiently large integer }m.} \tag{57}
\]

This is an eventual theorem, not yet the complete `r=6` column: the proof
does not record an explicit threshold.  An effective root-isolation and
remainder estimate for (54)--(56), followed by finite modular gcds below
that threshold, is the direct remaining task.  The exact Newton edge,
squarefreeness, and linear reconstruction are replayed by
[`verify_contact_resultant_r6_asymptotic.py`](../scripts/verify_contact_resultant_r6_asymptotic.py).

A stronger certificate is available at the limiting edge.  The
29 roots consist of one positive real root and fourteen conjugate pairs.
Their expansions `y=1+x_0/m+x_1/m^2+O(m^-3)` and
`z=z_0+z_1/m+O(m^-2)` are reconstructed explicitly, and rational Rouche
disks plus exact exponential bounds prove `|z_0|!=|exp(x_0)|` on every
branch.  This removes the transcendence input from the limiting separation,
but still does not supply a uniform positive-`t` threshold.  The 15
complex-conjugation classes and replay are in
[`R6_BRANCH_ATLAS.md`](R6_BRANCH_ATLAS.md) and
[`explore_contact_resultant_r6_branch_atlas.py`](../scripts/explore_contact_resultant_r6_branch_atlas.py).

## 12. Effective completion of the `r=6` column

The limiting atlas can be continued effectively all the way to the existing
finite certificates.  In the chart `t=1/m`, `y=1+tx`, cancellation of the
common denominator in the normalized sixth-moment equation gives

\[
 z^6={10\over (1+tx)Q(t,x)},                         \tag{58}
\]

where

\[
\begin{aligned}
Q={}&10t^6x^6+147t^5x^6-10t^5x^5
 +812t^4x^6-137t^4x^5+10t^4x^4\\
&+2205t^3x^6-675t^3x^5+125t^3x^4-10t^3x^3\\
&+3150t^2x^6-1530t^2x^5+525t^2x^4-110t^2x^3+10t^2x^2\\
&+2268tx^6-1620tx^5+900tx^4-360tx^3+90tx^2-10tx\\
&+648x^6-648x^5+540x^4-360x^3+180x^2-60x+10.
\end{aligned}                                       \tag{59}
\]

Thus a common root with `z=y^m` would satisfy the single scalar identity

\[
 10=Q(t,x)(1+tx)^{6/t+1}.                            \tag{60}
\]

Partition `0<=t<=1/41` into 256 rational cells.  On each cell, numerical
continuation proposes centers for the 29 roots of `P_6(t,x)`.  These centers
are not trusted.  With 320-bit Arb balls, Taylor expansion in `x` verifies
the strict Rouche inequality

\[
 |a_1(t)|r>|a_0(t)|+\sum_{j=2}^{29}|a_j(t)|r^j       \tag{61}
\]

uniformly on the cell, after `P_6(t,c+w)=sum a_j(t)w^j`.  The 29 resulting
disks are pairwise disjoint.  Hence every cell has 29 certified root tubes,
one root per tube, and these exhaust the eliminant.  There are
`256*29=7424` such tubes.

For every nonreal tube put

\[
 \Delta(t,x)=\log(Q(t,x)/10)
 +(6+t){\log(1+tx)\over t}.                          \tag{62}
\]

The quotient in (62) is evaluated regularly at `t=0` by its convergent
series

\[
 {\log(1+tx)\over t}
 =x\sum_{k\ge0}{(-tx)^k\over k+1},                  \tag{63}
\]

with an explicit geometric tail.  All nonreal tubes have `|tx|<1`.
Equation (60) would force `Re Delta=0` and `Im Delta` to be an integral
multiple of `2*pi`.  Arb enclosures exclude the first condition on 5920
tubes and exclude the second on the remaining 1504.  The unique positive
real tube is invariant under conjugation and therefore contains a real root;
direct real logarithmic bounds give `Delta>0` there.  This proves

\[
 \boxed{\operatorname{Res}_w(K_{m,6},L_{m,6})\ne0
        \quad\text{for every integer }m\ge41.}       \tag{64}
\]

Finally, degree-preserving reduction modulo `1,000,003` gives monic gcd one
for every `1<=m<=40`.  Combining the finite and tail certificates yields

\[
 \boxed{\operatorname{Res}_w(K_{m,6},L_{m,6})\ne0
        \quad\text{for every integer }m\ge1.}        \tag{65}
\]

The complete replay is
[`verify_contact_resultant_r6_effective.py`](../scripts/verify_contact_resultant_r6_effective.py).
It requires `python-flint`; floating-point root finding supplies candidates
only, while every accepted tube and inequality is certified by Arb.

## 13. Finite single-prime certificates

As finite evidence, direct endpoint gcds modulo `1,000,003` give no common
factor in the following grid:

| `r` | certified `m` |
|---:|---:|
| 5 | `1<=m<=50` |
| 6 | `1<=m<=40` |
| 7 | `1<=m<=30` |
| 8 | `1<=m<=25` |
| 9 | `1<=m<=20` |
| 10 | `1<=m<=16` |
| 11 | `1<=m<=12` |
| 12 | `1<=m<=10` |

All endpoint denominators and both leading coefficients are units modulo the
chosen prime.  Degree is therefore preserved, and a modular monic gcd equal
to one certifies coprimality over `Q`.  These 203 certificates are reproduced
by
[`verify_contact_resultant_modular_grid.py`](../scripts/verify_contact_resultant_modular_grid.py).
For `r=6` and `r=7` these are subsets of the finite parts of the complete
column certificates.  The rows with `r>=8` are evidence only and do not
replace a uniform argument.

## 14. Scope boundary (`OP-CR`)

The all-parameter problem is now confined to the finite staircase in each
column given by `m>=1001`, `8<=r<=X_m/m`, after removing pairs covered by
another parameter-irreducibility theorem.  Above the non-numerical threshold
`K_0`, it is further confined by `r<5(mr)^(21/40)`.  Formula (6) still gives
a fixed comparison disk, but the `r=4` analysis shows that demanding every
endpoint-eliminant root lie outside it is too strong.  For `r=6`, (58)--(65)
complete the whole column without a whole-eliminant Schur--Cohn attack.

The branch mechanism is structural for a fixed `r`, but it is not presently
uniform in `r`.  For arbitrary fixed `r` the boundary-layer limit has the
same form

\[
 \begin{aligned}
 e_r(c,z)&=(-c)^rz^{r-1}+\sum_{k=1}^{r-1}(-1)^k\binom rk
 z^{r-k-1}\bigl(b_k-z^kU_k(c)\bigr)(-c)^{r-k-1},\\
 f_r(c,z)&=b_r-z^rU_r(c).
 \end{aligned}                                      \tag{66}
\]

As a bounded pilot, exact elimination at `r=7` gives a residual endpoint
eliminant `H_7` of bidegree `(42,126)` in `(y,m)`.  Its complete top Newton
edge is

\[
 \operatorname{Res}_z(e_7,f_7)=\lambda c^7P_{42}(c),
 \qquad \lambda\in\mathbb Q^*,                     \tag{67}
\]

where `P_42` is squarefree, `P_42(0)!=0`,
`gcd(P_42,U_7)=1`, and the penultimate subresultant is linear in `z` with
leading coefficient coprime to `P_42`.  Thus the limiting separation and
unique finite `z` reconstruction used for `r=6` survive structurally in the
next column.  The Newton-edge identity proves that all 42 branches of `H_7`
have `y=1+c/m+O(m^-2)` and a unique algebraic limiting `z`.  Since `c` is
nonzero algebraic, Lindemann--Weierstrass separates that value from `exp(c)`.
Consequently the `r=7` contact resultant is nonzero for every sufficiently
large `m`.

The positive-`t` continuation is also effective in this column.  Put
`t=1/m`, `y=1+tx`, and let `P_7(t,x)` be the normalized degree-42 residual
eliminant.  The seventh-moment equation reduces to

\[
 720=Q_7(t,x)(1+tx)^{7/t+1},
\]

where `Q_7` is an explicit polynomial of bidegree `(7,7)`.  On the 512
rational cells partitioning `0<=t<=1/101`, 320-bit Arb arithmetic certifies
42 pairwise-disjoint Rouche tubes per cell.  These 21,504 tubes exhaust the
eliminant.  There are no real branches.  Rigorous logarithmic enclosures
exclude the scalar identity by its real part on 19,456 tubes and by its phase
on the remaining 2,048.  Thus the contact resultant is nonzero for every
`m>=101`.  Degree-preserving modular gcds modulo `1,000,003` certify
`1<=m<=100`, and therefore

\[
 \boxed{\operatorname{Res}_w(K_{m,7},L_{m,7})\ne0
        \quad\text{for every integer }m\ge1.}
\]

What does not yet generalize structurally is a bound uniform in `r`; the
branch count already rises from 29 to 42.  Repeating the atlas column by
column is therefore a valid fixed-`r` computation, not a solution of the
residual wedge.  The full bounded pilot and the smaller cross-column limiting
audit are replayed by
[`verify_contact_resultant_r7_asymptotic.py`](../scripts/verify_contact_resultant_r7_asymptotic.py)
and
[`verify_contact_resultant_fixed_r_branch_schema.py`](../scripts/verify_contact_resultant_fixed_r_branch_schema.py),
respectively.  The effective continuation and finite completion are replayed
by
[`verify_contact_resultant_r7_effective.py`](../scripts/verify_contact_resultant_r7_effective.py).

Two structural routes remain useful beyond that first column.  First,
combine the factor-degree restriction after (47) with endpoint geometry.  A
bound, depending only on `r`, for the degree of a common factor would exclude
the `mr-u` option immediately and leave the small degree `u` option for a
finite recurrence or moment-rank classification.  Valuations at `w=0,1` and
`infinity`, and the linear recurrences forced by a degree-`d` divisor, are
natural inputs to such a bound.

Second, the hypergeometric equation for `K=K_(m,r)` is

\[
 w(1-w)K''+\{r+2-(r+2-mr)w\}K'+mr(r+1)K=0.          \tag{68}
\]

Its roots avoid `0,1` and are simple.  Hence any differential-operator
formula for `L` reduces modulo (68) to `A_(m,r)(w)K'+B_(m,r)(w)K`; at a root
of `K`, coprimality becomes the local nonvanishing of `A_(m,r)`.  The useful
open step is not merely low differential order, which (68) supplies
automatically, but a contiguous-relation formula keeping `A_(m,r)` of
controlled algebraic complexity.

This continuation is tracked only as `OP-CR` in
[`STATUS.md`](../STATUS.md).  Raw interpolation of `Res_w(K,L)` still obscures this structure because its
degree grows with `mr`.

The general reduction and `r<=3` certificate are in
[`verify_contact_resultant_endpoint_reduction.py`](../scripts/verify_contact_resultant_endpoint_reduction.py).
The exact inertia, Rouche, angle, and Bernstein certificates for `r=4` are in
[`verify_contact_resultant_r4.py`](../scripts/verify_contact_resultant_r4.py).
The exact sparse resultant, subresultant, Routh, boundary-resultant, and
Bernstein certificates for `r=5` are in
[`verify_contact_resultant_r5.py`](../scripts/verify_contact_resultant_r5.py).
The exact degree-29 endpoint-eliminant reduction for the former first open column
`r=6` is in
[`verify_contact_resultant_r6_reduction.py`](../scripts/verify_contact_resultant_r6_reduction.py).
Its exact Newton edge, limiting subresultant, and eventual nonvanishing
certificate are in
[`verify_contact_resultant_r6_asymptotic.py`](../scripts/verify_contact_resultant_r6_asymptotic.py).
The effective 29-tube continuation and finite-range completion are in
[`verify_contact_resultant_r6_effective.py`](../scripts/verify_contact_resultant_r6_effective.py).
The bounded `r=7` limiting-edge generalization audit is in
[`verify_contact_resultant_fixed_r_branch_schema.py`](../scripts/verify_contact_resultant_fixed_r_branch_schema.py).
The exact full-eliminant connection and eventual `r=7` theorem are in
[`verify_contact_resultant_r7_asymptotic.py`](../scripts/verify_contact_resultant_r7_asymptotic.py).
The effective 42-tube continuation and finite-range completion are in
[`verify_contact_resultant_r7_effective.py`](../scripts/verify_contact_resultant_r7_effective.py).
The coefficient comparison and fractional-linear transfer behind (46)--(47)
are checked by
[`verify_contact_resultant_irreducible_ranges.py`](../scripts/verify_contact_resultant_irreducible_ranges.py).
