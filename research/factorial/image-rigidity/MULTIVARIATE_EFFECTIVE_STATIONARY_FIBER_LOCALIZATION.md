# Effective localization of stationary fibers on interior field rectangles

Delete the equal rows of an ordinary profile and write
\(d_i=b_i-a_i\), \(r_i=d_i/(a_i+b_i)\), \(w_i=2|d_i||r_i|\),
and \(W=\sum_i w_i\). Let \(z=2q-1\). The centered stationary
map in the upper half-plane is
\[
 X(z)=S(z)+\log(1-z)-\log(1+z),\qquad
 S(z)=\sum_i d_i\log(1+r_i z).
\]
Here all logarithms use their principal determinations.

**Theorem.** Fix \(0<W_-\le W_+<\infty\), \(M\ge0\) and
\(0<\epsilon\le\pi/2\). For any ordinary profile with
\(W_-\le W\le W_+\), every preimage \(z\) in the upper half-plane
of the closed field rectangle
\[
              |\Re X|\le M,\qquad -\pi+\epsilon\le\Im X\le-\epsilon
\]
satisfies \(|z|\le R\) and \(\Im z\ge h>0\), for the explicit
constants below. The same bounds hold for locally uniform limits
of the stationary maps. This gives a numerical version of the
[stationary-fiber compactness theorem](MULTIVARIATE_UNIFORM_SADDLE_MAP_COMPACTNESS.md)
away from the two horizontal boundaries of the field strip.

## The radius bound

Set
\[
 \eta=\min(\epsilon/12,1/8),\quad L=\max(4,8/\eta,8/\epsilon),\quad
 C_1=\pi(1+L)/\sin\eta,\quad C_2=C_1/\sin\eta,\quad
 r_*=\frac{W_-\eta}{8\pi}.
\]
A valid radius is
\[
 R=\max\left(
 2,\ 2e^{32(M+2)/W_-},\ 4C_1/W_-,\
 \frac{1+e^{M+2+C_2}}{r_*}
 \right).                                                     \tag{1}
\]
Every summand of \(\Im S\) is nonnegative, since \(d_i\) and \(r_i\)
have the same sign. The imaginary part of
\(\log(1-z)-\log(1+z)\) is in \((-\pi,0)\). Indeed its exponential
is in the lower half-plane and the difference of logarithms agrees
with its lower-half-plane logarithm near \(z=0\) and hence throughout.
Therefore on the preimage one has \(\Im S\le\pi-\epsilon\).
For \(|z|\ge2\), the real part of the logarithm difference has
absolute value at most \(\log3\), so \(|\Re S|\le M+2\).

First take \(0<\arg z\le\eta\) and write \(z=x(1+i\kappa)\).
Then \(0<\kappa\le2\eta\). For a row of negative sign set
\(B=|r_i|x\). If \(B\ge3/2\), then \(\arg(1-B-i\kappa B)\)
is within \(\arctan(3\kappa)\le\epsilon/2\) of \(-\pi\).
The integer gap \(|d_i|\ge1\) would then make its individual
imaginary contribution exceed \(\pi-\epsilon\).
For every negative row we thus have \(B<3/2\), and
\[
    1-|1-B-i\kappa B|^2=B[2-(1+\kappa^2)B]\ge B/8 .
\]
The real contribution of this row is at least \(|d_i|B/16=w_i x/32\),
by \(-\log u\ge1-u\) applied to the squared modulus.
A positive row contributes at least
\(d_i\log(1+r_i x)\ge w_i\log(1+x)/2\) by concavity.
Consequently \(\Re S\ge W_-\log(1+x)/32\), which gives
\(|z|\le2e^{32(M+2)/W_-}\).
For \(\pi-\eta\le\arg z<\pi\), replace \(z\) by \(-\bar z\)
and send every \((d_i,r_i)\) to \((-d_i,-r_i)\).
This sends \(S\) to \(-\bar S\), preserving the bounds used above.

Now let \(\eta\le\arg z\le\pi-\eta\). For a row with
\(|r_i z|\ge L\), the arguments of \(r_i z\) and \(1+r_i z\)
differ by at most \(\arcsin(1/L)\le2/L\). This follows by
drawing the tangent angles at zero to the radius-one disk
around \(r_i z\). A large row of each sign would together
contribute at least \(\pi-4/L>\pi-\epsilon\), since their integer
gaps have magnitude at least one. Thus all large rows have one
sign. Each corresponding factor has absolute argument at least
\(\eta/2\), and the sum of the absolute gaps in this group is
at most \(2\pi/\eta\).

For a small row \(v=r_i z\) with \(|v|<L\), the bound
\[
 |\arg(1+v)|\ge\frac{\sin\eta}{1+L}|v|
\]
follows from \(|\arg(1+v)|\ge|\Im(1+v)|/|1+v|\).
Also \(|\log(1+v)|\le |v|/\sin\eta\). For this, integrate
\(v/(1+t v)\) from zero to one and use \(|1+t v|\ge\sin\eta\),
since the ray of \(v\) makes angle at least \(\eta\) with either
real direction. Thus the small rows contribute at most \(C_2\)
in absolute real part and their \(w_i\)-mass is at most
\(2C_1/|z|\). If \(|z|>4C_1/W_-\), the large rows have
mass at least \(W_-/2\), so at least one satisfies \(|r_i|\ge r_*\).
Their real logarithms are all positive and their signed contributions
have the same sign. Therefore
\[
        |\Re S|\ge\log(r_*|z|-1)-C_2
\]
whenever \(r_*|z|>1\). Together with \(|\Re S|\le M+2\)
this proves the last radius bound in (1).

## A lower bound on the preimage height

Fix \(R\) from (1) and put
\[
 \rho=\frac1{2(R+1)},\qquad J=2+\lceil W_+(R+1)\rceil,\qquad
 b=\frac1{4J^2},\qquad C=-\log(\rho b),
 \qquad d=\min\bigl(b,e^{-M-W_+R-JC}\bigr).
\]
A valid height is
\[
                       h=\min\!\left(\frac d2,\
                              \frac{\epsilon}{W_++2J/d}\right). \tag{2}
\]

Split off the small rows \(|r_i|<\rho\). Their combined logarithm
has real part of modulus at most \(W_+R\) on \(|z|\le R\), using
\(|\log(1+r_i z)|\le2|r_i z|\). Its imaginary part is bounded in
modulus by \(W_+\Im z\), since each summand is real at \(\Re z\)
and its derivative along the intervening segment is bounded by
\(2|d_i r_i|\).

For the remaining rows, include the two base factors
\(\log(1-z)-\log(1+z)\), collect repeated factors exactly, and
remove any factors with zero resulting coefficient. This keeps
integer coefficients \(c_j\) with \(\sum_j|c_j|\le J\).
The poles \(p_j=-1/r_j\) are distinct real rationals of
denominator at most \(J\). For an ordinary row this follows
from \(p_j=-(a_i+b_i)/d_i\) and
\(|d_i|\le W_+/(2\rho)\), and the base poles are \(\pm1\).
Every two distinct poles are separated by at least \(1/J^2\).

If \(|z-p_j|<b\) for one pole, then all other poles are
at distance at least \(b\). Their real logarithms have absolute
value at most \(C\), since
\(\rho b\le |1+r_i z|\le1+R\) for their factors.
Hence the distinguished pole satisfies
\[
   |c_j|\,|\log|1+r_j z||\le M+W_+R+JC.
\]
Its integer coefficient has absolute value at least one.
As \(|r_j|\le1\), this implies \(|z-p_j|\ge d\).
Poles at distance at least \(b\) likewise have distance at least
\(d\), so the bound holds for all poles.

Write \(t=\Im z\). If \(t<d/2\), each real distance
\(|\Re z-p_j|\) is at least \(d/2\).
Thus the argument of \(1+r_j z\) is within \(2t/d\)
of an integer multiple of \(\pi\), by the elementary arctangent
bound. The combined retained argument with coefficients \(c_j\)
is within \(2Jt/d\) of a multiple of \(\pi\).
Including the small rows bounds the distance of \(\Im X\)
from the nearest integer multiple of \(\pi\) by
\((W_++2J/d)t\). The target rectangle requires this distance
to be at least \(\epsilon\), proving (2).

For locally uniform limits, use the retained rows and diffuse term
\(\tau z/2\) from the stationary-fiber compactness theorem.
The diffuse term contributes at least \(\tau\log(1+x)/32\)
to the real lower bound near the positive real axis.
In the central angles, its imaginary contribution is at least
\(\tau|z|\sin\eta/2\) and its modulus is \(\tau|z|/2\).
It therefore satisfies the same small-row bounds there, counting
\(\tau\) as part of the small-row mass. On the finite disk
its bounds are immediate from the linear formula. The large
retained poles and their integer coefficients are unchanged.
This proves (1) and (2) for the limiting maps as well.

## Finite analytic approximation of the profile family

On any disk \(|z|\le R+1\), choose \(0<\rho'\le1/(2(R+1))\).
There are finitely many enumerably bounded row lists with
\(|r_i|\ge\rho'\), because each such row has
\[
  |d_i|\le W_+/(2\rho'),\qquad a_i+b_i\le W_+/(2\rho'^2),
  \qquad w_i\ge2\rho'.
\]
The total number of retained rows is at most \(W_+/(2\rho')\).
Replace all omitted rows by the linear term \(\tau z/2\),
where \(\tau\) is their total \(w_i\)-mass. Then
\[
 \left|\sum_{|r_i|<\rho'}d_i\log(1+r_i z)-\tau z/2\right|
                      \le \frac12 W_+\rho'(R+1)^2 .           \tag{3}
\]
Indeed \(|\log(1+u)-u|\le |u|^2\) for \(|u|\le1/2\),
by integration of \(u/(1+tu)-u\) for \(0\le t\le1\).
Quantizing \(\tau\) within a mesh size \(\delta\) adds at most
\(\delta(R+1)/2\) to (3). Since the omitted-factor error is
analytic across the full disk, Cauchy's estimate bounds its
\(k\)-th derivative on \(|z|\le R\) by \(k!\) times the bound
on the larger disk. This gives finite rational profile nets for
the phase map and its derivatives. Quantitative contour activation
and the estimates near the endpoints of the field strip remain
separate steps in extracting a terminal cutoff.
