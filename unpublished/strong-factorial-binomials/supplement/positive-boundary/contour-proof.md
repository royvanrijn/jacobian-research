# Contour estimates for the positive-exponent boundary

Write $Q_{m,\gamma,d}(x)=\sum_{k=0}^m\binom mk(\gamma)_{dk}x^{m-k}$.
We prove three estimates for the article's large relative gaps.

| Support | Onset | Required conclusion |
| --- | --- | --- |
| $a=1,\ d=6$ | $N\ge10000$ | $Q_{N,N,d}$ is squarefree |
| $a=1,\ d\ge10$ | $N\ge\max(5000,128(d+1))$ | $Q_{N,N,d}$ is squarefree |
| $a\ge2,\ d\ge3a$ | $N\ge\max(4000,128(a+d))$ | The actual two consecutive moments have no common root |

The article gives the derivative identity that converts each
squarefreeness conclusion into the assertion at orders $N-1,N$.
The third row uses two integrals with consecutive exponents.

## Contour and Gaussian estimates

All the integrals below start as an entire polynomial in $s$
times $e^{-Ns}$. On a local horizontal segment through a saddle
$s_*$, write its phase as
$$
 F(s_*+u)=F(s_*)-Au^2/2+E(u),\qquad
 |E(u)|\le B|u|^3,\quad |u|\le r .
$$
Here $A=-F''(s_*)$. If $\Re A/2-Br\ge b>0$, the identity
$e^v-1=v\int_0^1 e^{tv}\,dt$ bounds the difference of the
two exponentials, after removing $e^{mF(s_*)}$, by
$mB|u|^3e^{-mbu^2}$. The Gaussian tail outside $[-r,r]$
is at most $e^{-mbr^2}/(mbr)$, obtained by replacing $1$
by $u/r$ on the positive tail and doubling.

For amplitude $1/s$, put $v=|s_*|$. In units
of $e^{m\Re F(s_*)}/m$, the local error including Gaussian
completion is bounded by
$$
 \frac{B}{(v-r)b^2}
 +\frac{1}{v(v-r)b}
 +\frac{e^{-mbr^2}}{vbr}.
$$
The first term integrates $mB|u|^3e^{-mbu^2}/(v-r)$.
The second integrates the amplitude difference, bounded by
$|u|/(v(v-r))$. With amplitude one the corresponding terms are
$B/b^2$ and $e^{-mbr^2}/(br)$. With normalized amplitude
$e^{-u}$, multiply the phase-error integrand by $e^{|u|}$
and add $|u|e^{|u|}e^{-mbu^2}$ for the amplitude difference.
The Gaussian completion term is unchanged.

The formulas used for the complex phase estimates come from
$(\log(1+t))''=t''/(1+t)-(t')^2/(1+t)^2$ and
$$
 (\log(1+t))'''=\frac{t'''}{1+t}
  -\frac{3t't''}{(1+t)^2}+\frac{2(t')^3}{(1+t)^3}.
$$
For $t$ a scalar times $s^k$ and $|t|\le q_0<1$, these
give bounds by $|s|^{-2}$ and $|s|^{-3}$ times
$$
 |k(k-1)|y+k^2y^2,\qquad
 |k(k-1)(k-2)|y+3k^2|k-1|y^2+2|k|^3y^3,
 \qquad y=\frac{q_0}{1-q_0}.
$$
They apply on the indicated local disks or horizontal strips,
where the phase logarithms exist. Elsewhere we estimate the modulus
of the original entire integrand. Horizontal deformations use
Cauchy's theorem. The vertical closing segment at infinity tends to
zero by exponential decay.

## Gap six

For $d=6$, the Gamma integral and the rescaling of its variable
give
$$
 Q_{N,N,6}(N^6w)=\frac{N^{7N}}{(N-1)!}I_N(w),\qquad
 I_N(w)=\int_0^\infty [s(w+s^6)]^N e^{-Ns}\frac{ds}{s}.
$$
The apparent pole at zero cancels in the original integrand.
Split $I_N=K_N+J_N$ at $s=3$, with $K_N$ on the short
interval.

### The two saddle approximations

For $J_N$ use the phase
$$
 \Phi(s,w)=7\log s+\log(1+w/s^6)-s.
$$
For $w^{-N}K_N(w)$ use
$$
 \Psi(s,w)=\log s+\log(1+s^6/w)-s.
$$
In either case the saddle equation is $w(s-1)=s^6(7-s)$. Write $s_J,s_K$ for the following solutions.

- For $|w|\le3180$, Rouché's theorem on $|s-7|=1/5$ gives exactly one solution there, counted with multiplicity. The inequality is $(1/5)(34/5)^6/(31/5)>3180$.
- For $|w|\ge1000$, write the equation as $s-1=w^{-1}s^6(7-s)$. On $|s-1|=1/100$ the right side has modulus at most $(101/100)^6(601/100)/1000<1/100$. Thus there is exactly one solution in that disk. This solution is analytic in $1/w$, including $1/w=0$.

The logarithms at the saddles are the branches continued from $w=0$ for $\Phi$, and $1/w=0$ for $\Psi$. They can also be taken by the convergent local logarithmic series. At either saddle the negative second derivative is $a(s)=(s^2-7s+7)/s^2$. Denote the real parts of the two saddle heights by $S_J,S_K$. The elementary logarithmic bound $|\log(1+v)-v|\le |v|^2/(2(1-|v|))$ gives
$$
 S_J>7\log7-7-\tfrac1{25}>6.58,\qquad S_K>-1-\tfrac1{500}.
$$
For the first bound, at $s_J=7+v$ the height difference from $7\log7-7$ is $7\log(1+v/7)-\log(1+v/6)-v$. For the second use $|s_K-1|<1/100$ and $|s_K^6/w|\le(101/100)^6/1000$. The curvature satisfies
$$
 \Re a(s_J)>\tfrac3{25},\quad |a(s_J)|<\tfrac15,
 \qquad \Re a(s_K)>\tfrac9{10},\quad |a(s_K)|<\tfrac{11}{10}.
$$
For example $a(s)-1/7=(s-7)(6s-7)/(7s^2)$ for the first disk and $a(s)-1=-7(s-1)/s^2$ for the second.

The saddle estimates on their full stated domains are
$$
\begin{split}
 J_N(w)&=\frac{e^{N\Phi(s_J,w)}}{s_J}\sqrt{\frac{2\pi}{N a(s_J)}}(1+\epsilon_J(w)),\\
 K_N(w)&=w^N\frac{e^{N\Psi(s_K,w)}}{s_K}\sqrt{\frac{2\pi}{N a(s_K)}}(1+\epsilon_K(w)),
 \qquad |\epsilon_J|,|\epsilon_K|<\frac{25}{\sqrt N}.       \tag{2}
\end{split}
$$
The square roots use $\Re a>0$. The proof retains all portions of the deformed contours.

For the local integral at $s_*$, integrate along $s_*+t$ for $|t|\le \rho$. The choices for $J$ are $(\rho,B,c)=(3/4,1/30,7/200)$, and for $K$ they are $(1/5,5/6,1/4)$. Here $B$ bounds the coefficient in the phase Taylor remainder $B|t|^3$, and $c\le\Re a/2-B\rho$. To check $B$, let $r$ be an upper bound on $|w/s^6|$ or $|s^6/w|$ along the local disk and let $y=r/(1-r)$. Differentiating the phases bounds the absolute third derivatives by
$$
 \frac{14+336y+756y^2+432y^3}{|s|^3},
 \qquad
 \frac{2+120y+540y^2+432y^3}{|s|^3},
$$
respectively. At the radii just specified these are less than $1/5$ and $5$. If $s_{\min}$ is $34/5$ or $99/100$, the local Gaussian replacement and its completion to the whole real line cost, after division by $e^{NS_*}$,
$$
 \frac1N\left[
 \frac{B}{(s_{\min}-\rho)c^2}
 +\frac{1}{s_{\min}(s_{\min}-\rho)c}
 +\frac{1}{s_{\min}c\rho}\right].                          \tag{3}
$$
In fact the exponential Taylor error is bounded by $NB|t|^3 e^{-Nc t^2}$, the amplitude difference by $|t|/(s_{\min}(s_{\min}-\rho))$, and the two Gaussian tails by $e^{-Nc\rho^2}/(Nc\rho)$. The displayed constants in (3) are less than $11$ and $43$.

For $J_N$ deform to the vertical segment $[3,3+i\Im s_J]$ and the horizontal ray from its endpoint. The right vertical side at infinity vanishes. The following bounds control the nonlocal part.

- For $3\le x\le 28/5$, the modulus of the exponential base on the deformed path is at most $(1+1/450)^7 f_{3200}(x)$, where $f_R(x)=x(R+x^6)e^{-x}$. The same bound at $x=3$ applies on the vertical segment. This follows from $|x+iy|\le x(1+y^2/(2x^2))$ and $|y|\le1/5$. On $3\le x\le28/5$, the sign of $f_R'(x)$ is the sign of $x^6(7-x)/(x-1)-R$. The ratio is increasing here, since its derivative has the sign of $-x^2+7x-7$ between 1 and 7. Thus the maximum of $f_{3200}$ is at an endpoint. Both endpoint values are below $705$, and $705(1+1/450)^7<e^{6.575}$. The base loses at least $1/200$ in log modulus from the saddle height.
- For $28/5\le x\le9$, the second derivative of $\Re\Phi(x+i\Im s_J,w)$ is bounded above by $-3/(2x^2)\le-1/54$. Indeed if $r=3200/(28/5)^6$, the correction to the log curvature has modulus at most $(42r/(1-r)+36r^2/(1-r)^2)/x^2$; the term from $7\log s$ has negative magnitude at least $7(1-1/28^2)/(1+1/28^2)^2x^{-2}$. Using stationarity, outside the local interval the log loss is at least $\rho^2/108=1/192$.
- For $x\ge9$, the real phase derivative is at most $7/x+6r_x/(x(1-r_x))-1<-1/5$ with $r_x=3200/x^6$. The log loss at nine is at least $(9-36/5)^2/108=3/100$.

Including $1/|s|$, these three parts cost at most $e^{NS_J}[3e^{-N/200}+(5/(9N))e^{-3N/100}]$. The total error including (3) is below $13e^{NS_J}/N$ at $N\ge10000$. Here $3e^{-N/200}<1/N$ follows from the fourth-order exponential lower bound. The Gaussian leading factor has modulus greater than $(3/4)e^{NS_J}/\sqrt N$ by $|s_J|<36/5,\ |a|<1/5,\ \pi>3$.

For $K_N/w^N$, keep $[0,1/2]$ on the real axis and deform the rest to the horizontal segment at height $\Im s_K$, with vertical segments at $1/2$ and 3. For $1/2\le x\le3/2$ the second derivative of $\Re\Psi(x+i\Im s_K,w)$ is at most $-1/(2x^2)\le-2/9$. In fact the $\log s$ term contributes negative magnitude at least $(1-(1/50)^2)/(1+(1/50)^2)^2x^{-2}$, and the correction has modulus at most $(30r/(1-r)+36r^2/(1-r)^2)/x^2$ for $r=(151/100)^6/1000$. Off the local interval this gives log loss at least $\rho^2/9=1/225$.

On the horizontal interval $[3/2,3]$ the modulus of the base is less than $17/50$. Indeed the real envelope is $f_{1000}(x)/1000$, with its maximum at the endpoints by the same monotonicity argument, and the modulus correction is at most $(1+(1/100)^2/(2(3/2)^2))^7$. The same bound holds on the vertical segment at three. At $1/2$ use $|s|\le 1/2+1/100$ to get an upper bound below $17/50$. The real segment $[0,1/2]$ costs at most $(17/50)^N$: the positive envelope $x^{N-1}(1+x^6/1000)^Ne^{-Nx}$ increases there for $N\ge2$, and its endpoint bound times the interval length is the base bound to the $N$-th power. Since $17/50<e^{-1-1/500-1/225}$, all complements including the concavity interval cost at most $5e^{NS_K-N/225}$, with the amplitudes included. For $N\ge10000$ this is less than $e^{NS_K}/N$, again by a fourth-order exponential bound. The total error is below $44e^{NS_K}/N$, and the leading factor has modulus greater than $2e^{NS_K}/\sqrt N$. This proves (2) with room in both constants. The errors in (2) are analytic since the original pieces and the nonzero saddle approximants are analytic on neighborhoods of the stated domains.

### Cancellation on the overlap has a simple zero

The region $|w|\le1100$ has no zero. On the real short interval $x(1100+x^6)e^{-x}\le1100/e+3^7e^{-3}<520$. Splitting at $1/2$ to handle $dx/x$ as above gives $|K_N|\le6(520)^N$. By (2), $|J_N|>e^{NS_J}/(2\sqrt N)>(720)^N/(2\sqrt N)$. The ratio is bounded by $12\sqrt N(3/4)^N<1$.

The region $|w|\ge3080$ also has no zero. On $[3,4]$ bound the real exponential base of $J_N/w^N$ by $3e^{-3}+4^7e^{-4}/3080<8/25$. On $[4,\infty)$ use $4e^{-4}+7^7e^{-7}/3080<8/25$. The envelope $x(1+x^6/3080)e^{-x}$ has logarithmic derivative at most $7/x-1$, hence at most $-3/10$ for $x\ge10$. Integrating the envelope gives $|J_N/w^N|<4(8/25)^N$. Equation (2) gives $|K_N/w^N|>(9/25)^N/\sqrt N$, and $4\sqrt N(8/9)^N<1$. Both ratio inequalities follow already from $e^v\ge v^2/2$.

On the remaining annulus $1100\le|w|\le3080$, a disk of radius 100 around $w$ lies in both saddle domains. Cauchy's estimate applied to (2) yields $|\epsilon_*'/(1+\epsilon_*)|\le1/300$. Also
$$
 |s_*'(w)|\le\frac{|s_*-1|^2}{6|s_*|^7|a(s_*)|},
 \qquad
 \left|\frac{d}{dw}\log\big(s_*^{-1}a(s_*)^{-1/2}\big)\right|<1/50 .
$$
The first bound is implicit differentiation of the saddle equation, using $w+s_*^6=6s_*^6/(s_*-1)$; the second uses $|a_s|\le7/|s|^2+14/|s|^3$ and the disk bounds. Differentiation of the stationary phases now gives
$$
 \left|\frac{J_N'}{J_N}-\frac{N}{w+s_J^6}\right|<\tfrac1{20},
 \qquad
 \left|\frac{K_N'}{K_N}-\frac{N}{w+s_K^6}\right|<\tfrac1{20}.
$$
The two main terms are separated by more than
$$
 N\left(\frac1{3102}-\frac1{(34/5)^6-3100}\right)>\frac N{3300}.
$$
Here $|s_K|^6<2$ and $|w|\le3080<3100$. At a zero of $I_N$ in this annulus, $K_N=-J_N\ne0$, so $I_N'=J_N(J_N'/J_N-K_N'/K_N)\ne0$. This excludes multiple roots uniformly at $N\ge10000$.

## Lower exponent one and gaps at least ten

Set $D=d+1$, $H=D^D e^{-D}$, $W=eH$, and $c=3D/8$.
The same Gamma expansion gives
$$
 Q_{N,N,d}(N^d w)=\frac{N^{DN}}{(N-1)!}(J_N(w)+K_N(w)).
$$
Here $J_N$ and $K_N$ integrate
$[s(w+s^d)]^N e^{-Ns}/s$ on $[c,\infty)$ and $[0,c]$.
On $|w|\le3W$ use $\Phi=D\log s+\log(1+w/s^d)-s$
for $J_N$. On $|w|\ge W/3$ use
$\Psi=\log s+\log(1+s^d/w)-s$ for $K_N/w^N$.

### Stationary points and local errors

Both saddle equations are $w(s-1)=s^d(D-s)$.
Write $q_*=3/8$, $q_1=2/5$. We use
$$
 \frac{163}{60}<e<\frac{87}{32},\quad W>D4^d>400D^2,\quad
 40D^2q_*^d<1,\quad 14D^3q_*^d<2,\quad
 300Dq_*^d<1,\quad 800Dq_1^d<1,\quad 8D^2q_1^d<1 .
$$
The exponential comparisons follow from Taylor bounds.
For the ones with $q_*^d$ and $q_1^d$, check $d=10$
and compare successive terms using
$(12/11)^3q_*<1$ and $(12/11)^2q_1<1$.
The lower bound on $W=D(D/e)^d$ uses $D\ge11$.
The bound $4^d>400D$ holds at $d=10$ and its ratio grows.

On $|s-D|=1/10$, $|w(s-1)/s^d|<4D^2q_*^d<1/10$.
On $|s-1|=1/(100D)$, $|s^d(D-s)/w|<4D/W<1/(100D)$.
The powers of $|s|/D$ and $|s|$ in these two comparisons
are bounded by $e^{1/10}<10/9$ and $e^{1/100}<4/3$,
respectively. Rouché's theorem gives one simple saddle $s_J,s_K$
in each disk, analytic in $w$ and $1/w$, respectively.
The curvature at either saddle is
$A=(s^2-Ds+D)/s^2$. The identities
$$
 A-1/D=(s-D)((D-1)s-D)/(Ds^2),\qquad
 A-1=-D(s-1)/s^2
$$
give $\Re A_J>4/(5D), |A_J|<6/(5D)$ and
$\Re A_K>9/10, |A_K|<11/10$. For the first disk
$D|A-1/D|\le(1/10)(1+1/(10D))/(1-1/(10D))^2<1/5$.
For the second disk $|A-1|\le(1/100)/(1-1/(100D))^2<1/10$.
The real saddle heights $S_J,S_K$ satisfy
$$
 |S_J-\log H|<1/100,\qquad |S_K+1|<1/100 .
$$
In the first comparison, the Taylor error of $D\log s-s$
is less than $3/(500D)$, and the logarithmic correction
has modulus less than $5Dq_*^d$.
In the second these bounds are $1/(10000D^2)$ and $5/W$.
They follow by summing the series for $\log(1+t)$.
The Gaussian prefactors
$\sqrt{2\pi/(NA_*)}/s_*$ have moduli between $2$ and $3$
times $1/\sqrt{ND}$ and $1/\sqrt N$, respectively,
using the saddle-disk and curvature bounds and $3<\pi<22/7$.

Use $(r,B,b)=(1/2,1/D^2,1/(3D))$ and
$(1/10,1,1/3)$ at $J,K$, respectively.
In the upper local disk $|s-D|<3/5$,
$|w/s^d|\le6Dq_*^d$, $y\le7Dq_*^d<1/100$,
and $|s|>47D/50$.
The third-derivative estimate above bounds its main term by
$3/D^2$ and its correction by $14Dq_*^d<2/D^2$.
For the correction, divide by $y$ and bound the remaining factor
by $(12/11+3/100+2/10000)(50/47)^3<2$.
On the lower local disk $|s-1|\le 101/1000$,
$|s^d/w|\le3q_1^d/D$, $y\le4q_1^d/D$.
Use $|s|>89/100$ to bound the base third derivative by $3$
and the correction by $8D^2q_1^d<1$.
The Gaussian decay inequalities are
$2/(5D)-1/(2D^2)>1/(3D)$ and $9/20-1/10>1/3$.
Keeping the completion term with its exponential replaced by one,
the upper local bracket is bounded by
$$
 \frac9{D-3/5}
 +\frac{3D}{(D-1/10)(D-3/5)}
 +\frac{6D}{D-1/10}
 =6+\frac{63}{5D-3}<\frac{29}{4}.
$$
The lower bracket at $v=999/1000$ is less than $46$.
These give relative errors $ (29/8)\sqrt{D/N}$
and $23/\sqrt N$ on the local parts with Gaussian completion.

### Remainder of the upper contour

Deform the upper ray from $c$ to height $\Im s_J$.
On the central interval $[3D/4,5D/4]$, $q_0=3D\,2^{-d}$
bounds $|w/s^d|$. The second-derivative formula above
gives a correction of modulus at most $D^2(y+y^2)/x^2$
for $s=x+i\Im s_J$. The quantity $D(y+y^2)$ decreases
with $d$, since $D^2 2^{-d}$ and $D^3 4^{-d}$ decrease
and the denominator $1-q_0$ increases. At $d=10$ it is
$371712/982081<2/5$. The base curvature has negative real
part of magnitude exceeding $(99/100)D/x^2$ because
$|\Im s_J|/x<1/80$. Thus $\Re\Phi''<-1/(3D)$.
The central Gaussian tails beyond $r=1/2$, including
$1/|s|\le4/(3D)$, cost at most
$(16/N)e^{NS_J-N/(24D)}$.
On the far ray $x\ge5D/4$ the sharper ratio
$q_0=3D(3/10)^d$ gives $y<1/100$.
Hence $\Re\Phi'\le(1+y)D/x-1<-1/10$.
The central bound gives a loss at its right endpoint of at least
$(D/4-1/10)^2/(6D)>D/200$, so the far-ray cost
is at most $10e^{NS_J-ND/200}/(DN)$.

On $[3D/8,3D/4]$ and its left vertical segment, the modulus
of $s(w+s^d)e^{-s}$, divided by $H$, is at most
$$
 e^{1/(25D)}
 \left[(x/D)^De^{D-x}+3ex e^{-x}\right].
$$
The correction follows from $|s|/x\le e^{(\Im s_J)^2/(2x^2)}$.
Split the interval at $D/2$: the bracket is bounded by
$1/8+9/16$ on the first piece and $2/3+3/16$ on the second.
For the first terms use their values at $D/2,3D/4$, respectively,
and $e<87/32$. For the last terms use $x=3D/8,D/2$,
respectively, and $e>163/60$.
These bounds hold at $D=11$; each exponential base
for the first terms is below one and $D e^{-c_0D}$ decreases
for $c_0=3/8,1/2$ on $D\ge11$.
The bracket is less than $7/8$, so
$(7/8)e^{11/100+1/(25D)}<1$ puts it at least $1/10$
below the saddle in log modulus. The integrated amplitude
is less than $\log 2+(1/10)/c<1$.

### Remainder of the lower contour and simplicity

Keep $[0,1/2]$ real and move $[1/2,c]$ to height
$\Im s_K$, retaining the endpoint verticals.
On $[1/2,3/2]$, $q_0=3q_1^d/D$ bounds $|s^d/w|$
by $|s|<1501/1000$ and $W>D4^d$.
The curvature correction is at most $8Dq_1^d/x^2<1/(100x^2)$.
The base negative curvature has magnitude at least $99/(100x^2)$
since $|\Im s_K|/x<1/100$. This gives
$\Re\Psi''<-1/3$. Both central tails cost at most
$(120/N)e^{NS_K-N/600}$, using $r=1/10$ and amplitude
at most $2$.

On $[0,1/2]$ the real integrable envelope
$x^{N-1}(1+3x^d/W)^N e^{-Nx}$ is increasing.
Its endpoint value times $1/2$ is less than $(8/25)^N$.
On the vertical at $1/2$ the exponential base is less than
$33/100$, with integrated amplitude less than
$2|\Im s_K|<1$. Both cost at most $e^{NS_K-N/100}$.
On $[3/2,c]$ and its right vertical segment the real envelope
is at most $xe^{-x}+(3/e)(x/D)^D e^{D-x}$.
The first term is less than $27/80$ at $x=3/2$.
The second is bounded at $x=c$ by
$(9/8)[(3/8)e^{5/8}]^D<9/400$. At $D=11$ this last
comparison follows from
$(9/8)^8(3/8)^{88}(87/32)^{55}<(9/400)^8$,
and the base is below one.
The imaginary correction is at most $e^{1/(45000D)}$.
Thus the envelope is at least $1/100$ below $S_K$ in
log modulus, since $(9/25)e^{1+21/1000}<1$.
All lower pieces outside the central interval together cost at most
$4D e^{NS_K-N/100}$ by integrating $1/x$ on the far
horizontal interval and including its vertical end.

Divide by the Gaussian leading moduli. The total relative errors
$\epsilon_J,\epsilon_K$ are bounded respectively by
$$
 \sqrt{D/N}\left(29/8+8e^{-N/(24D)}\right)
  +\tfrac12\sqrt{ND}\,e^{-N/10}
  +\frac{5}{\sqrt{ND}}e^{-ND/200},
\qquad
 \frac{23+60e^{-N/600}}{\sqrt N}+2D\sqrt N\,e^{-N/100}.
$$
At $N\ge5000,\ N/D\ge128$ the three upper terms are
less than $33/100,1/1000,1/1000$.
For the first use $e^{-N/(24D)}<(3/8)^5$; for the
other two use $e^x\ge x^4/24$ and $e^x\ge x^2/2$.
The lower terms are less than $33/100,1/2000,1/1000$,
using $N/600>8$ and the decrease of
$N^{3/2}e^{-N/100}$ on $N\ge5000$.
Both errors are less than $1/3$.

In particular the magnitude prefactors after including
$1+\epsilon_*$ lie between $1$ and $4$ times
$1/\sqrt{ND}$ and $1/\sqrt N$. The height errors are at
most $1/100$. All coefficients of $K_N$ in $w$, and of
$J_N/w^N$ in $1/w$, are nonnegative.
For $|w|\le W/2$ compare the upper saddle lower bound
with $K_N(W/2)$; for $|w|\ge2W$ compare the lower
saddle lower bound with $J_N(2W)/(2W)^N$.
The smaller-to-larger ratio in the two comparisons is bounded by
$4\sqrt D\,2^{-N}e^{N/50}<1$ and
$4\,2^{-N}e^{N/50}/\sqrt D<1$.
Indeed $4\sqrt D<N$ and $N(2/3)^N<1$ here.
There are no zeros outside $W/2<|w|<2W$.

On that annulus the disk of radius $W/12$ stays within
both saddle domains. The relative errors are analytic because
the leading terms are analytic and nonvanishing.
Cauchy's estimate gives $|\epsilon_*'/(1+\epsilon_*)|<6/W$.
A disk of radius $W/2$ for $s_J-D$ and one of radius
$W/12$ for $s_K-1$ give
$|s_J'|\le1/(5W)$, $|s_K'|\le12/(100DW)$.
From $A_s=D/s^2-2D/s^3$ the logarithmic derivatives with respect
to $s$ of $s^{-1}A^{-1/2}$ have moduli less than
$3/2$ and $3D$, respectively. Their derivatives with respect
to $w$ are less than $1/W$.
At a stationary point, differentiating the full exponential
part, including $w^N$ on the lower packet, gives $N/(w+s_*^d)$.
Thus these are the main terms of $J_N'/J_N,K_N'/K_N$,
with errors less than $7/W$. Since $|s_K|^d<2<W$ and
$|s_J|^d>D^d/2$, their unscaled moduli satisfy
$$
 |(w+s_K^d)^{-1}|>1/(3W),\qquad
 |(w+s_J^d)^{-1}|<3/D^d<1/(100W).
$$
For the second bound use $2W<D^d/6$ and $300Dq_*^d<1$.
Their difference exceeds $1/(4W)$.
At a zero $K_N=-J_N\ne0$ so
$I_N'=J_N(J_N'/J_N-K_N'/K_N)\ne0$, since $N>56$.
This proves squarefreeness at the claimed onset.

## Lower exponent at least two and gap at least three times it

Let $D=a+d$, $\rho=D/a\ge4$, $n=N-1$, $H=D^De^{-D}$,
$W=He^a/a^a$, $c=5D/8$, and $g(s)=s^a(w+s^d)$.
The actual consecutive moments, after division by their nonzero
monomial scaling factors, are
$$
 I=\int_0^\infty g(s)^N e^{-Ns}\,ds,\qquad
 J=\int_0^\infty g(s)^n e^{-Ns}\,ds.
$$
In both, set $w=A/(BN^d)$. Split at $c$, with lower and
upper parts denoted by $L,U$. The shared phase is
$\Phi=\log g-s$, and the amplitude for the $n$-th power
is $e^{-s}$. Upper estimates apply for $|w|\le6W/5$,
lower estimates for $|w|\ge4W/5$.

### Uniformity in the support

For an upper ratio at $x=\beta D$ and a lower ratio
at $x=\beta a$, use $q_U=(6/5)e^{aF_U}$ and
$q_L=(5/4)e^{aF_L}$, where
$$
 F_U=\log\rho-(\rho-1)(1+\log\beta),\qquad
 F_L=(\rho-1)(1+\log\beta)-\rho\log\rho .
$$
In particular $|w/s^d|\le q_U$ for $\Re s=x=\beta D$.
The bound for $|s^d/w|$ is $q_L$ if $|s|\le\beta a$;
the right lower strip below includes the imaginary correction explicitly.
The positive majorants can be treated as functions of real $a,\rho$.

- For $\beta\ge5/8$, $F_U<0$ and
  $F_{U,\rho}=1/\rho-1-\log\beta<-1/4$, using
  $e^{1/2}>8/5$ and $4<(8/3)^3(5/8)^3$.
  On $[5/8,7/8]$, $q_U<1$ by its value at the corner
  $a=2,\rho=4,\beta=5/8$.
  The logarithmic $\rho$-derivative of
  $(1-1/\rho)q_U/(1+q_U)$ is
  $1/[\rho(\rho-1)]+aF_{U,\rho}/(1+q_U)<1/12-1/4$.
- For $\beta\ge7/8$, the products
  $d(d+1)q_U/D$ and $d^2q_U^2/D$ also decrease in $a,\rho$.
  The first polynomial factor has logarithmic derivatives
  at most $1/2,\,2/3$, respectively. The opposing contributions
  are $F_U\le\log 4-18/7<-1/2$ and
  $aF_{U,\rho}\le-17/14$, by $\log(7/8)>-1/7$.
  The second product has no larger polynomial derivatives
  and doubles the opposing contributions.
- For $\beta\le7/4$, $F_L<0$ and
  $F_{L,\rho}=\log(\beta/\rho)\le\log(7/16)<-1/2$.
  The sign of $F_L(4,7/4)$ follows from
  $(7/4)^3(87/32)^3<4^4$.
  When $r_0=(100/99)q_L\le1/4$, the logarithmic
  $\rho$-derivative of $(\rho-1)r_0/(1+r_0)$
  is at most $1/3-4/5<0$.
- At $\beta=9/8$, the polynomial factor $d^3/a$
  has logarithmic derivatives at most $1,1$;
  $F_L<-47/24$ and $aF_{L,\rho}<-29/12$.
  At $\beta=251/200$ the corresponding bounds for
  $d(d-1)/a$ are $3/5,11/15$, and
  $F_L<-941/600,\ aF_{L,\rho}<-647/300$.
  These follow by substituting $\log4>4/3$ and
  $\log\beta<\beta-1$. For $d^2/a$ the polynomial
  derivatives are at most $1/2,2/3$.
  These comparisons prove decrease for $(d^3/a)q_L$,
  $(d(d-1)/a)q_L$, and $(d^2/a)q_L^2$ at their stated
  values of $\beta$.

Replacing powers of a decreasing $q_0<1$ by the same powers
of $q_0/(1-q_0)$ preserves all the product comparisons,
since the extra factor $(1-q_0)^{-1}$ decreases.
The corner in these comparisons is $a=2,\rho=4,D=8,d=6$.

### Stationary points and local bounds

Stationarity is $w(s-a)=s^d(D-s)$.
On $|s-D|=1/2$ and $|s-a|=1/(2\rho)$ the two
Rouché ratios are at most
$$
 R_U=\frac{12}{5}(d+1/2)\rho^a e^{-d}(1-1/(2D))^{-d},
 \qquad
 R_L=\frac52\rho(d+1/(2\rho))\rho^{-\rho a}e^{d+1/2}.
$$
Both decrease with $a,\rho$. For $R_U$, put
$t=-\log(1-1/(2D))$; then $t<1/(2D-1)$ and
$\partial_a(dt)=(\rho-1)[t-1/(2D-1)]<0$.
The remaining logarithmic derivative bounds are
$1/a+\log\rho-(\rho-1)<0$ and
$1/(\rho-1)+a/\rho-a+a/(2D-1)<0$.
For $R_L$ they are
$1/a+\rho-1-\rho\log\rho<0$ and
$1/\rho+1/(\rho-1)-a\log\rho<0$.
The corner rational bounds use $e>163/60$ for $R_U$
and $e<87/32,\ e^{1/2}<5/3$ for $R_L$:
$$
 R_U\le\frac{85761906966528}{93776847890045}<1,\qquad
 R_L\le \frac52\cdot4\cdot\frac{49}8\cdot\frac53
            \frac{(87/32)^6}{4^8}<1 .
$$
Thus the selected saddles $s_U,s_L$ are unique, simple and
analytic in their respective domains.

Their curvature is $A(s)=[s^2-(D+a-1)s+aD]/s^2$.
To find real-part minima, the inversion $v=1/s$ sends a
disk of center $c_0$ and radius $r$ to one of center
$c_0/(c_0^2-r^2)$ and radius $r/(c_0^2-r^2)$.
On the boundary the derivative of
$\Re[1-(D+a-1)v+aDv^2]$ with respect to the boundary
cosine has the sign of
$$
 -(D+a-1)+\frac{2aD(c_0+2r\cos\theta)}{c_0^2-r^2}.
$$
For the upper disk the positive term is at most $16a/7$
since $D\ge8$; for the lower disk it is greater than
$2D-2$. The signs are negative and positive, respectively,
as $D\ge4a$. The minima are consequently at
$s=D-1/2$, $s=a+a/(2D)$. Evaluation gives
$$
 \Re A_U\ge\frac{(D+a)/2-1/4}{(D-1/2)^2}>\frac1{2D},\qquad
 a\,\Re A_L\ge
 \frac{1/2+(a+1)/(2D)+a/(4D^2)}{(1+1/(2D))^2}>\frac12 .
$$
The identities
$A-1/D=(s-D)((D-1)s-aD)/(Ds^2)$ and
$A-1/a=(s-a)((a-1)s-aD)/(as^2)$ give
$|A_U|<3/D$, $|A_L|<2/a$ by the disk radii.

For the heights $S_*=\Re\Phi(s_*)$, radial integration
of the upper stationary value, whose derivative is
$1/(w+s_U^d)$, gives an error at most $q_0/(1-q_0)<7/80$
around $\log H$, using $q_0=32(60/163)^6$.
Indeed at the upper saddle the ratio $|w/s_U^d|$ is less
than $(6/5)\rho^ae^{-d}e^{1/2}$; this is bounded at the corner
since $\log\rho-(\rho-1)<0$ and decreases in $\rho$.
At the lower saddle the Taylor error of $a\log s-s$
is at most $a/[4D(2D-1)]\le1/240$.
Stationarity gives $|s_L^d/w|\le1/(2\rho d-1)\le1/47$.
Thus
$$
 |S_U-\log H|<7/80,\qquad
 |S_L-\log|w|-a\log a+a|<1/32.
$$

Use local triples $(r,B,b)=(1/2,3/D^2,1/(7D))$ and
$(a/16,4/(3a^2),1/(6a))$.
On $|s-D|\le1$, the ratio $|w/s^d|$ is bounded by
$(6e/5)\rho^a e^{-d}$ because
$(1-1/D)^{-d}<e$. This majorant and its product with $D^2$
decrease in $a,\rho$. The logarithmic derivatives of the
product are at most $1+\log 4-3$ and $1/2+1/2-2$.
Take its corner bound $q_0=(96/5)(60/163)^5$.
The third derivative is bounded, after multiplication by
$D^2$, by $2(8/7)^3+64(y+3y^2+2y^3)<18$.
For the upper decay inequality, insert the sharper curvature
minimum above with $a=2$. After multiplication by
$56D^2(D-1/2)^2$, the required difference is
$6(D-8)^3+89(D-8)^2+354(D-8)+187>0$.
The lower disk lies in $|s-a|\le a/8$.
The third derivative times $a^2$ is bounded by
$(8/7)^3[2+(d^3/a)(y+3y^2+2y^3)]<8$
using the lower ratio at $\beta=9/8$ and its corner bound.
The lower decay follows from $1/(4a)-(4/(3a^2))(a/16)=1/(6a)$.

### Complex strip and endpoint bounds

Move the upper half-line to height $\Im s_U$, including
the vertical at $c$. Retain $[0,a/2]$ on the real axis
for the lower contour, and move $[a/2,c]$ to height
$\Im s_L$, including both endpoint verticals.

On the upper central $[7D/8,5D/4]$, the base log curvature
has negative magnitude exceeding $(39/40)D/x^2$.
The perturbation is bounded by
$[d(d+1)y+d^2y^2]/x^2$, so the monotonicity above reduces
$\Re\Phi''<-1/(4D)$ to the corner at each fixed $\beta=x/D$.
On the upper left $[5D/8,7D/8]$, the map
$t\mapsto(1+t)^{-1}$ takes the disk $|t|\le r_0$ to
the disk of center $1/(1-r_0^2)$, radius $r_0/(1-r_0^2)$.
Using it on $t=w/s^d$ gives, with $|\Im s|\le1/2$,
$$
 \Re\Phi'\ge D/x-1-\frac{d r_0}{x(1+r_0)}
                       -\frac{D}{4x^3(1-r_0^2)} .
$$
Here $\Re(1/s)\ge1/x-1/(4x^3)$ and $|1/s|\le1/x$.
Both corrections decrease at fixed $\beta$, so
$\Re\Phi'>1/100$ again reduces to the corner.
On the far ray $x\ge5D/4$ the same ratio majorant is
below $1/50$. Then $\Re\Phi'\le(1+y)D/x-1<-1/10$.

On the lower central $[3a/4,5a/4]$ one has
$|s|\le(251/200)a$ and $|\Im s|/x\le1/12$.
The base curvature has magnitude at least $(39/40)a/x^2$.
Using the two lower curvature correction comparisons at
$\beta=251/200$, the corner inequality needed is
$$
 39/20-30y-36y^2>
 (3/20)(251/100)^2=189003/200000.
$$
It gives $\Re\Phi''<-3/(10a)$.
On the right strip $[5a/4,7a/4]$, put $\beta=x/a$ and
take $r_0=(100/99)q_L$, which is less than $1/4$ at
the corner. This absorbs $(|s|/x)^d\le e^{1/(8D\beta^2)}<100/99$.
Apply the disk bound to $t/(1+t)=1-(1+t)^{-1}$ with
$t=s^d/w$. Using $|1/s|\le1/x$ and
$\Re(1/s)\ge1/x-(\Im s)^2/x^3$ gives
$$
 \Re\Phi'\le a/x-1+\frac{d r_0}{x(1+r_0)}
                 +\frac{d r_0^2}{4\rho^2 x^3(1-r_0^2)}.
$$
The correction terms decrease at fixed $\beta$ by the
third comparison under uniformity, and by decrease of
$(\rho-1)/(4\rho^2a^2\beta^3)$ in the last term.
Thus $\Re\Phi'<-1/20$ reduces to the corner.

Here are the three corner polynomials for the strip comparisons.
Put
$A_0/B_0=(3/10)8^8(60/163)^6$ and
$A_1/B_1=6022586125125/49539595901075456$.
These are positive rational numbers with the displayed numerators
and denominators, or any equivalent common scaling.
Clearing positive denominators in the upper central, upper left
and lower right inequalities gives respectively
$$
\begin{split}
 &-5B_0^2x^{14}+1248B_0^2x^{12}
     +10A_0B_0x^8-9216A_0B_0x^6-5A_0^2x^2+2208A_0^2,\\
 &-101B_0^2x^{15}+800B_0^2x^{14}-200B_0^2x^{12}
     -600A_0B_0x^8+101A_0^2x^3-200A_0^2x^2,\\
 &-152A_1^2x^{15}+1280A_1^2x^{14}-15A_1^2x^{12}
     -960A_1B_1x^8+152B_1^2x^3-320B_1^2x^2 .
\end{split}
$$
They are positive on $[7,10],[5,7],[5/2,7/2]$,
respectively. For a polynomial $P(l+(h-l)t)=\sum c_i t^i$
of degree $m$ the Bernstein coefficients are
$b_j=\sum_{i\le j}c_i\binom ji/\binom mi$.
The identity
$P(l+(h-l)t)=\sum_j b_j\binom mj t^j(1-t)^{m-j}$
proves positivity on the entire interval when all $b_j>0$.
The exact checker independently clears the denominators and
recomputes all $15+16+16$ positive coefficients.

The central curvatures yield phase losses $u^2/(8D)$ and
$3u^2/(20a)$ at horizontal displacement $u$ from the two
saddles. At the upper endpoints $7D/8,5D/4$ these are at
least $D/5120,D/1280$, and at the lower right endpoint
$5a/4$ at least $9a/2560$, by the saddle radii.
On the upper vertical at $c$, the real high contribution
divided by $H$ is at most $[(5/8)e^{3/8}]^D<15/32$.
The real low contribution is at most
$(6/5)[(5\rho/8)e^{1-5\rho/8}]^a<3/8$.
Both decrease from the corner. With imaginary correction at most
$e^{8/(25D)}\le e^{1/25}$, the loss from $S_U$
is at least $1/50$, since $(27/32)e^{59/400}<1$.
The vertical length is at most $1/2$.

The lower left pieces outside $[3a/4,5a/4]$ have total
length less than $3a/2$. For fixed permitted imaginary part
the base $|s|^a e^{-x}$ increases on $[a/2,3a/4]$:
its log derivative is $ax/(x^2+(\Im s)^2)-1>0$.
The largest imaginary correction is obtained with
$|\Im s|=a/(2D)$. Hence $|s|^2\le (145/256)a^2$
at the endpoint used for the envelope. Also $|s^d/w|<1/500$
by the corner bound at $\beta=151/200$.
Dividing the envelope by $|w|a^a e^{-a}$ gives
$[(\sqrt{145}/16)e^{1/4}]^a(501/500)$.
The bracket is less than one. With $a\ge2$, multiplying
its upper bound at $a=2$ by $e^{1/32+1/100}$ still gives
less than one. This includes the real origin segment by
monotonicity on $[0,a/2]$ and gives a loss $1/100$.

Finally write $x=av$ on the lower far part, $7/4\le v\le5\rho/8$.
The real envelope divided by $|w|a^a e^{-a}$ is at most
$$
 [ve^{1-v}]^a+\tfrac54[(v/\rho)^\rho e^{\rho-v}]^a.
$$
Both bases are below one. The high base decreases in $\rho$
for fixed $v$, since its logarithmic derivative is
$\log(v/\rho)$. On $[7/4,5/2]$, take the corner and bound the
low and high terms at the left and right endpoints on each of
$[7/4,2],[2,9/4],[9/4,19/8],[19/8,5/2]$.
For $v\ge5/2$ bound the low term at $5/2$ and the high term at
$5\rho/8$, where it decreases with $\rho$.
The five endpoint envelopes, including the imaginary factor
$e^{1/196}$ and the height adjustment $e^{1/32}$, are each
less than $e^{-1/200}$. Their arithmetic is checked with rational
Taylor bounds. The far part with its right vertical has length $<2D$.

### Errors for the actual two integrals

The Gaussian leading terms for each saddle are
$$
 L_{I,*}=e^{N\Phi(s_*)}\sqrt{2\pi/(NA_*)},\qquad
 L_{J,*}=e^{n\Phi(s_*)-s_*}\sqrt{2\pi/(nA_*)}.
$$
After factoring out the displayed exponentials their magnitudes
exceed $\sqrt{2D/N}$ at $U$ and $\sqrt{3a/N}$ at $L$
for $I$, with $N$ replaced by $n$ for $J$.

To control the amplitude $e^{-u}$ uniformly as the local
interval grows with $a$, use $|u|\le nbu^2/24+6/(nb)$.
In the local estimate the phase and amplitude terms then gain
multipliers $(24/23)^2e^{6/(nb)}$ and $(24/23)e^{6/(nb)}$.
At the lower saddle $n/a>511$, so bound them by
$C_p=(24/23)^2(511/475)$, $C_a=(24/23)(511/475)$.
The resulting local relative error before completion is bounded by
$48C_p/\sqrt{3na}+6C_a\sqrt{a/(3n)}<9/20$.
For $2\le a\le125/16$, use $n\ge3999$ and convexity in
$\sqrt a$, checking the two endpoints. For $a\ge125/16$,
use $n\ge512a-1$; both terms then decrease.
At the upper saddle $n/D>127968/1000$ gives multipliers
$<5/3,8/5$. With the weaker leading bound $\sqrt{D/n}$
the upper weighted local error before completion is less than
$(5/3)(1000/999)(33/40)+(8/5)(5/8)<5/2$.

The following rows give upper bounds for the separate **relative**
errors, in the order local before completion, Gaussian completion,
central tail, left part, vertical or right strip, and far part.
Write $E(t)=e^{-t}$.
$$
\begin{array}{c|l}
 I_U&
 147/\sqrt{ND},\
 14\sqrt{D/N}E(N/(28D)),\
 16\sqrt{D/N}E(N/(32D)),\
 100E(ND/5120)/\sqrt{ND},\
 \tfrac12\sqrt{N/D}E(N/50),\
 10E(ND/1280)/\sqrt{ND}\\
 J_U&
 5/2,\
 14\sqrt{D/n}E(n/(28D)),\
 32\sqrt{D/N}E(N/(64D)),\
 200E(ND/10240)/\sqrt{ND},\
 \tfrac12\sqrt{N/D}E(N/100),\
 10E(ND/2560)/\sqrt{ND}\\
 I_L&
 24/\sqrt N,\
 48E(N/768)/\sqrt N,\
 160E(3N/2560)/(3\sqrt N),\
 (N/16)E(N/100),\
 10E(9N/1280)/\sqrt N,\
 D\sqrt N E(N/200)\\
 J_L&
 9/20,\
 48E(n/768)/\sqrt n,\
 320E(3N/5120)/(3\sqrt N),\
 (N/16)E(N/128),\
 10E(9n/1280)/\sqrt n,\
 D\sqrt N E(N/200).
\end{array}
$$

For the unweighted terms integrate the Gaussian bounds from
the central curvatures and the exponential bounds from the
left and right derivative strips. The two left $L$ terms use
$(3a/2)/\sqrt{3a/N}=(\sqrt3/2)\sqrt{aN}<N/16$ because
$N\ge512a$. In the weighted estimate one may use $\sqrt{3a/n}$
instead. The worst amplitude factor there is $e^{a+1/8}$;
$-n/100+a+1/8\le-N/128$ by $a\le N/512$.
For the weighted central $U$ tails, use
$nu^2/(8D)-|u|\ge Nu^2/(16D)$ for $|u|\ge1/2$
and $N\ge80D+2$. For the weighted central $L$ tails,
$3nu^2/(20a)-|u|\ge3Nu^2/(40a)$ for $|u|\ge a/16$
and $N\ge216$. The upper left weighted derivative is at
least $n/100-1\ge N/200$; at its right endpoint
$x_U-7D/8\le3D/16$ is absorbed by half the phase loss
because $N\ge1922$.
On the upper vertical $x_U-c\le7D/16$ so
$-n/50+x_U-c\le-N/100$ when $N\ge50D$.
On the upper far ray the amplitude only improves the estimate
and $n\ge N/2$. On the lower right strip the amplitude only
improves the estimate because $x>x_L=\Re s_L$.
On the lower far segment $x_L-x+1/200<0$, so its phase
loss remains $N/200$ in the weighted bound.

The rational budgets for the rows sum to less than
$5/6,\,3,\,2/5,\,18/25$, respectively.
For the upper unweighted local term use $ND\ge32000$
with $D\le17$, or $ND\ge72000$ with $D\ge18$.
The two upper Gaussian tails depend only on $N/D$,
bounded below by $4000/17$ and $128$ in these cases.
For the weighted terms use $n>999N/1000$.
All growing square-root factors multiplying decaying
exponentials in the other terms decrease above $N=4000$,
after using $D\le N/128$ for the last lower terms.
Thus the numerical proof is by endpoint bounds on the listed
monotone quantities, with separate exact Bernstein positivity
on the three strips.

### Excluding a common zero

The saddle height bounds differ by at most $19/160$ from
their baseline difference $\log(W/|w|)$.
At $|w|=4W/5$ the upper height exceeds the lower by more
than $1/10$, and at $|w|=6W/5$ the lower exceeds the upper
by more than $1/20$. These follow from
$e^{7/32}<5/4$ and $e^{27/160}<6/5$.
The relative error bounds and Gaussian prefactors give
$|I_L/I_U|\le54 e^{-N/10}$ at the inner circle and
$|I_U/I_L|\le10\sqrt D e^{-N/20}$ at the outer circle.
For the inner disk, $I_L$ has nonnegative coefficients
in $w$ and is bounded by its value at $w=4W/5$;
the lower bound for $I_U$ holds throughout the disk.
For the exterior region use the nonnegative coefficients of
$I_U/w^N$ in $1/w$, and the lower estimate
for $|I_L/w^N|$ there. These are the same bounds,
and both are less than one at the stated onset.

At a zero of $I$ both saddle estimates apply and
$I_U=-I_L\ne0$. The leading ratio of the two integrals at each
saddle is $\sqrt{N/n}/g(s_*)$. The four relative error bounds
give $|J_U/I_U|<24\sqrt{N/n}/|g(s_U)|$ and
$|J_L/I_L|>\sqrt{N/n}/(5|g(s_L)|)$.
But $\log|g(s_U)/g(s_L)|=S_U-S_L+\Re(s_U-s_L)$.
The height bounds and radii give a lower bound
$d-119/160-\log(6/5)$, whose exponential exceeds $150$
already at $d=6$.
Thus $|J_U/J_L|<120/150<1$ and $J\ne0$.
This proves the required assertion at the actual consecutive orders.
