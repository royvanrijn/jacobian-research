# The first `r=6` branch atlas at `m=infinity`

This note resolves the asymptotic branch search left open by the degree-29
endpoint eliminant in [CONTACT_RESULTANT.md](CONTACT_RESULTANT.md).  It is a
tail theorem, not yet a complete proof for every integer `m`: the argument
proves separation for all sufficiently large `m`, but does not extract an
explicit threshold to connect to the existing certificates through `m=40`.

## 1. The unique Newton chart

Write

\[
 \operatorname{Res}_z(E_6,F_6)=(y-1)^7H_6(m,y),
 \qquad \deg_yH_6=29,\quad\deg_mH_6=90.
\]

Exact coefficient extraction gives the staircase

\[
 [m^{90-k}]H_6(m,y)\quad\hbox{is divisible by}\quad
 (y-1)^{29-k}\qquad(0\le k\le29),                 \tag{1}
\]

and the top coefficient is a nonzero scalar times `(y-1)^29`.  Every
coefficient of `y^j`, `0<=j<=29`, has exact `m`-degree 90.  Consequently no
root escapes to `y=infinity`, every root tends to `y=1`, and the only Newton
chart is

\[
 t=\frac1m,\qquad y=1+tx,
 \qquad P_6(t,x)=t^{61}H_6(1/t,1+tx).              \tag{2}
\]

The transform is a polynomial of bidegree `(29,89)` in `(x,t)`.  Its edge
polynomial

\[
 P_0(x)=P_6(0,x)                                  \tag{3}
\]

has degree 29 and is squarefree and irreducible over `Q`.  It has one
positive real root, no negative real root, and fourteen nonreal conjugate
pairs.  Thus all 29 branches are ordinary power-series branches in `t`; no
fractional Puiseux exponents occur.

## 2. Limiting endpoint equations and reconstruction of `z`

After dividing `E_6` and `F_6` by their exact chart valuations `t^6` and
`t^7`, respectively, their limits are, up to nonzero rational scalars,

\[
\begin{aligned}
 \bar E_0={}&A(x)z^5-64800000x^4z^4-40500000x^3z^3\\
 &-16000000x^2z^2-3796875xz-497664,\\
 A(x)={}&10800000x^6-28260000x^5+28260000x^4+12507000x^3\\
 &+7033300x^2+1308555x+497664,\\
 \bar F_0={}&
 (324x^6-324x^5+270x^4-180x^3+90x^2-30x+5)z^6-5.
                                                               \tag{4}
\end{aligned}
\]

Their exact resultant is `c*x^7*P_0(x)` for a nonzero rational `c`.  The
penultimate subresultant is linear in `z`, so on every root of `P_0`

\[
 z_0=-\frac{B_0(x_0)}{A_0(x_0)}.                  \tag{5}
\]

The denominator is coprime to `P_0`, and the `(x,z)` Jacobian of (4) is
nonzero at all 29 reconstructed points.  The implicit-function theorem
therefore gives

\[
\begin{aligned}
 x(t)&=x_0+x_1t+O(t^2),\\
 y(t)&=1+x_0t+x_1t^2+O(t^3),\\
 z(t)&=z_0+z_1t+O(t^2).                            \tag{6}
\end{aligned}
\]

The checker computes `x_1=-P_{6,t}(0,x_0)/P_0'(x_0)` and obtains `z_1`
from the first-order part of (4).

## 3. Atlas

Only one member of each complex-conjugate pair is shown.  The column `gap` records
the rigorously certified sign of `|z_0|-|exp(x_0)|`.  The argument is included
as a numerical diagnostic, not as part of the certificate.

| class | `x_0` | `x_1` | `z_0` | `z_1` | `|z_0/e^x_0|` | arg | gap |
|---:|---:|---:|---:|---:|---:|---:|:---:|
| 1 | 48.368549 | 903.731871 | 0.010351278 | -0.283583124 | 1.020469558e-23 | 0 | `<` |
| 2 | -13.1520778-17.0668996i | 65.9686419+167.062650i | -0.0228261718-0.00319835376i | -0.201899430-0.148990691i | 11872.16212 | 1.49814788 | `>` |
| 3 | -1.62984302-1.44686568i | 0.993274680+4.07522474i | -0.0299871617-0.213853107i | 0.248190880-0.286180489i | 1.101984856 | -0.263245468 | `>` |
| 4 | -0.839494281-0.469646016i | -1.41418811+0.798766253i | 0.335960667+0.267330595i | -0.624898908-0.00326015632i | 0.9940126023 | 1.14177170 | `<` |
| 5 | -0.427076043-0.524942379i | 0.00170935041-0.173361362i | 0.531489632-0.321283549i | -0.265034261-0.166185210i | 0.9519278562 | -0.0187767129 | `<` |
| 6 | -0.233788321-0.133075746i | 0.139220065-0.0622738224i | 0.301271998-0.731718571i | -0.0326668649-0.115961908i | 0.9997271656 | -1.04714145 | `<` |
| 7 | -0.165771077-0.187763349i | -0.0542740006+0.111734340i | 0.279045710-0.800021055i | 0.0496003632+0.0629283794i | 1.000058718 | -1.04742942 | `>` |
| 8 | -0.146535495-0.421325453i | 0.328534624+0.181411200i | -0.709016355-0.515567824i | -0.202887783-0.361164254i | 1.014998859 | -2.09154512 | `>` |
| 9 | -0.0623593508-0.167748393i | 0.0170839783+0.0835507440i | -0.599008812-0.723843256i | 0.0354392813-0.0649717144i | 1.000008742 | -2.09435708 | `>` |
| 10 | -0.0521604812-0.462988750i | -0.214675637+0.269295199i | 0.861653125-0.396259999i | 0.0745621725+0.244875554i | 0.9991849649 | 0.0319462212 | `<` |
| 11 | -0.000608941107-0.414348669i | -0.317260800+0.0623309583i | 0.117055706-0.977703754i | 0.126501758+0.233552220i | 0.9852858744 | -1.03728971 | `<` |
| 12 | 0.0302941249-0.376778506i | -0.285716475-0.0211155670i | -0.799730087-0.632078579i | 0.235090291+0.130412291i | 0.9889408887 | -2.09597595 | `<` |
| 13 | 0.0333891310-0.181709236i | -0.0319927248+0.0830464199i | -1.01694111+0.186929077i | -0.000346161055-0.0935155301i | 1.000024931 | 3.14151611 | `>` |
| 14 | 0.144166585-0.365692831i | -0.0490835111+0.0498968430i | -1.08671292+0.446460456i | -0.0740664278+0.0194736596i | 1.017118541 | 3.11747303 | `>` |
| 15 | 0.167590431-0.192435507i | -0.0604342573+0.0412638078i | -0.385162593+1.11648991i | -0.0624139799-0.0948609852i | 0.9988216053 | 2.09542423 | `<` |

## 4. Exact limiting-modulus certificate

The numerical roots used to seed the calculation are not trusted.  Around
each 30-digit Gaussian-rational center the checker takes the disk of radius
`10^-24`.  If

\[
 P_0(c+w)=\sum_{j=0}^{29}a_jw^j,
\]

it verifies with rational arithmetic that

\[
 |a_1|_{\rm lower}r>
 |a_0|_{\rm upper}+\sum_{j=2}^{29}|a_j|_{\rm upper}r^j.       \tag{7}
\]

Rouche's theorem puts exactly one root in each disk.  The 29 disks are
pairwise disjoint, so they exhaust all roots.  Taylor bounds for the exact
rational function (5) enclose `|z_0|`; a 220-term rational Taylor sum with a
geometric tail encloses `exp(Re(x_0))`.  All 29 interval pairs are disjoint.
The smallest displayed relative gap is about `8.74e-6` (class 9), and it is
still certified with room to spare.

Finally,

\[
 y(t)^{1/t}=(1+tx(t))^{1/t}\longrightarrow e^{x_0}.           \tag{8}
\]

The strict inequalities above and the finiteness of the atlas imply a common
punctured neighborhood of `t=0` on which `|z(t)| != |y(t)^(1/t)|`.  Hence no
`r=6` branch can satisfy `z=y^m` for all sufficiently large positive integer
`m`.

## 5. Reproduction and remaining task

Run

```bash
.venv/bin/python scripts/explore_contact_resultant_r6_branch_atlas.py
```

The remaining useful computation is quantitative continuation in `t`: find
an explicit rational `t_0` for which all 29 modulus separations persist on
`0<t<=t_0`.  If `1/t_0<=41`, the existing exact finite grid closes the whole
`r=6` column.  The closest limiting classes by modulus are 9, 13, 7, and 6;
their large argument mismatches suggest that a mixed modulus/argument proof
may give a much better cutoff than modulus alone.
