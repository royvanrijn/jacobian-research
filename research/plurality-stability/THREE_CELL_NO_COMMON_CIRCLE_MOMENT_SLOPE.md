# A slope obstruction for three-cell partitions without a common circle

Let \(A\) be an ordinary balanced partition of the Gaussian plane
such that almost every centered circle has positive angular mass
in at most two labels. If \(A\) ties the regular simplex \(T\) at
\(0<r\le73/100\), then
\[
                          S'_r(A)>S'_r(T).
\]
In particular a maximizing partition at such a correlation must
have a common winner circle if simplex-optimal correlations
accumulate at \(r\) from the right. The statement constrains
this class of right ties. It supplies neither global optimality
on a larger correlation interval nor an exclusion of left ties.

## The first Hermite energy

Let \(s\) be radius, \(q=s^2/2\), and let \(d\mu(s)=s e^{-s^2/2}ds\)
be the radial measure. The variable \(q\) is exponential with
mean one. Write \(p_i(s)\) for the angular fraction of label \(i\),
and \(c_{j,i}(s)=\int 1_{A_i}(s,\theta)e^{ij\theta}\,d\theta/(2\pi)\)
for its angular Fourier coefficients. For a subset of angular
fraction \(p\), the magnitude of the first coefficient is at
most \(\sin(\pi p)/\pi\): after projecting along the phase of
that coefficient, the integral of cosine is maximized by the
centered arc of that fraction. For a binary slice, the two
nonzero labels' first coefficients are opposite. A pure
slice contributes zero.

Write \(M=\int s\,d\mu=\sqrt{\pi/2}\), \(A_0=M/\pi\), and
\(\alpha=17/20\). Tag each genuinely binary radius by its
active pair. Orient the three pair differences \(d_e=e_i-e_j\)
cyclically. For an active pair \(e=(i,j)\), let \(v_e\) be the
integral, over its tagged radii, of \(s c_{1,i}\), read as a
vector in the plane. Let \(\omega_e\) be the corresponding
integral of \(s|c_{1,i}|\), so \(|v_e|\le\omega_e\).
Then the first Hermite energy satisfies
\[
  W_1(A)=\left\|\sum_e d_e\otimes v_e\right\|^2
       \le \left(\sum_e\omega_e\right)^2+\sum_e\omega_e^2 .
                                                               \tag{NS1}
\]
Indeed \(|d_e|^2=2\) and \(|d_e\cdot d_{e'}|=1\) for distinct
pairs. Also \(\sum_e\omega_e\le A_0\). Each tagged set has
\(\mu\)-mass at most \(2/3\), since it omits a label of
Gaussian mass \(1/3\). Maximizing the integral of the radius
over such a set uses the outer \(2/3\) of the radial measure.
For \(a=\sqrt{2\log(3/2)}\) this gives
\[
 \omega_e\le \frac1\pi\left(M-\int_0^a s^2e^{-s^2/2}\,ds\right)
               <\alpha A_0.
\]
To check the last constant, \(a>9/10\) follows from the first
two positive terms of
\(\log(3/2)=2\sum_{j\ge0}(1/5)^{2j+1}/(2j+1)\).
The alternating exponential polynomial through degree three
in \(s^2/2\) gives
\[
 \int_0^a s^2 e^{-s^2/2}\,ds
  >\frac{(9/10)^3}{3}-\frac{(9/10)^5}{10}
       +\frac{(9/10)^7}{56}-\frac{(9/10)^9}{432}
  =\frac{21458663451}{112000000000}>\frac{189}{1000}
  >\frac{3M}{20}.
\]
Here \(M<63/50\) follows from \(\pi<22/7\).

Set \(u=\sum_e\omega_e/A_0\) and \(L=1-u\). If \(u\le\alpha\),
(NS1) gives \(W_1/A_0^2\le 2u^2\). If \(u\ge\alpha\), the
largest possible sum of squared weights with this sum and
individual cap \(\alpha A_0\) has one weight \(\alpha A_0\)
and another \((u-\alpha)A_0\). Indeed if two coordinates are
strictly between zero and the cap, transferring weight from
the smaller to the larger increases the sum of squares
until a boundary is reached. As \(u\le1<2\alpha\), the
bound is \(W_1/A_0^2\le u^2+\alpha^2+(u-\alpha)^2\).
The first bound divided by \(u>0\) is at most \(2\alpha\).
The second divided by \(u\) is increasing on \([\alpha,1]\),
since it equals \(2u-2\alpha+2\alpha^2/u\). Its value at
one is \(1+\alpha^2+(1-\alpha)^2>2\alpha\).
The case \(u=0\) has \(W_1=0\). Consequently
\[
       W_1(A)\le \overline W(1-L),\qquad
       \overline W=\frac{349}{400\pi}.
                                                               \tag{NS2}
\]

## The second Hermite energy

At each radius designate one absent label measurably and
let \(E_s\) be the slice of one of the other two labels,
including a pure or empty slice when applicable. Put
\(\ell(s)=1-\pi|c_1(E_s)|\), which gives
\(\int s\ell(s)\,d\mu=M L\). For angular fraction \(p\)
of \(E_s\) the first coefficient bound gives
\[
        |p-1/2|\le \tfrac12\sqrt{\ell(s)}.
\]
In fact \(1-\sin(\pi p)=1-\cos(\pi|p-1/2|)
\ge 4(p-1/2)^2\), by the chord bound for sine on
\([0,\pi/4]\).

Let \(H\) be the hemisphere centered at the phase of
\(c_1(E_s)\), or any hemisphere if the coefficient is zero.
For the direction of its center, the cosine has the correct
sign on both sides of the hemisphere. Thus
\(\ell=\pi\int_{E_s\triangle H}|\cos\theta|\,d\theta/(2\pi)\)
in coordinates centered on that direction. Let \(d_+\) and \(d_-\) be the angular masses of
\(E_s\setminus H\) and \(H\setminus E_s\) respectively.
For given masses the smallest integral of \(|\cos\theta|\)
on each side uses two intervals next to the zeros of
cosine, giving
\[
 \ell\ge (1-\cos(\pi d_+))+(1-\cos(\pi d_-))
       \ge4(d_+^2+d_-^2)
       =2\bigl(d^2+(p-1/2)^2\bigr),
\]
where \(d=d_++d_-\) and the inequality uses the same
sine chord bound since \(d_\pm\le1/2\). The second
coefficient of the hemisphere is zero, so
\(|c_2(E_s)|\le d\). Put
\[
 L_p=\frac2M\int s(p-1/2)^2\,d\mu,\qquad
 L_d=\frac2M\int s d^2\,d\mu.
\]
Then \(L_p+L_d\le L\) and the first coefficient
bound gives separately \(L_p\le L/2\).

An orthonormal basis of the planar degree-two space is
\(q-1,q\cos(2\theta),q\sin(2\theta)\).
For the last two basis functions the vector of slice
coefficients is supported on the two active labels with
opposite values. Its Euclidean norm is at most
\(\sqrt2 qd\). The squared norm of the integrated angular
part is at most
\[
             \left(2\int s d^2\,d\mu\right)
                       \int q^2/s\,d\mu=\frac{3M^2}{4}L_d,
\]
where \(\int q^2/s\,d\mu=3M/4\).

For the radial degree-two basis function, let \(M_i\) be
the set of radii where label \(i\) is designated absent.
Replace both other label fractions by \(1/2\), giving
\(\widetilde p_i=(1-1_{M_i})/2\).
The baseline radial vector has entries
\(-\frac12\int_{M_i}(q-1)\,d\mu\). The centered exponential
variable \(q-1\) has positive-part and negative-part integrals
both equal to \(1/e\). Thus the sum of positive entries
of the baseline vector and the absolute sum of its negative
entries are each at most \(1/(2e)\). Its squared norm is
at most \(1/(2e^2)\). The difference between the actual
and baseline fraction vectors has squared norm
\(2(p-1/2)^2\) at each radius. The identity
\(\int (q-1)^2/s\,d\mu=3M/4\) bounds its contribution
to the integrated radial vector by \(M\sqrt{3L_p/4}\)
in norm. Adding the radial and angular squares and using
\(L_p\le L/2\) in the cross term gives
\[
 W_2(A)\le B(L):=\frac1{2e^2}
        +\frac{\sqrt3\,M}{2e}\sqrt L
        +\frac{3M^2}{4}L .
                                                               \tag{NS3}
\]

## Comparing derivatives at a tie

For a balanced partition,
\(S_r=1/3+\sum_{j\ge1}r^j W_j\) with \(W_j\ge0\). Put
\(F(r)=S_r(T)-1/3\) and define
\[
          w(r)=2F(r)/r-S'_r(T),\qquad
          R(r)=S'_r(T)/r-F(r)/r^2 .
\]
In the simplex series \(w=W_1(T)-\sum_{j\ge3}(j-2)W_j(T)r^{j-1}\)
is nonincreasing, while
\(R=\sum_{j\ge2}(j-1)W_j(T)r^{j-2}\) is nondecreasing.

Suppose \(S_r(A)=S_r(T)\) and \(S'_r(A)\le S'_r(T)\).
The coefficient inequalities \(j\ge2\) for \(j\ge2\),
and \(j\ge3\) for \(j\ge3\), respectively, give
\[
             W_1(A)\ge w(r),\qquad
             2W_1(A)+r W_2(A)\ge 3F(r)/r-S'_r(T).
                                                               \tag{NS4}
\]
The explicit bounds below give \(w(7/10)>\overline W\).
Thus (NS2) and (NS4) already contradict each other
when \(0<r\le7/10\).

For \(7/10\le r\le73/100\), if \(w(r)>\overline W\)
the same immediate contradiction applies. Otherwise
\(L\le L_*:=1-w(r)/\overline W\). By (NS2) and (NS3),
\[
     2W_1+rW_2\le2\overline W(1-L)+r B(L)
                         \le 2w(r)+r B(L_*).
\]
The second inequality holds since \(r\,3M^2/4>2\overline W\)
makes the preceding expression increasing in \(L\).
As \(2w(r)+rR(r)=3F(r)/r-S'_r(T)\), it suffices to show
\(B(L_*)<R(r)\).

On \([7/10,18/25]\) the explicit bounds below give
\(w(r)\ge w(18/25)>(24/25)\overline W\) and
\(R(r)\ge R(7/10)>9/40\). Hence \(L_*<1/25\).
Using \(e>8/3\), \(M<63/50\), \(\sqrt3<7/4\),
and \(\pi<22/7\), we have
\[
       B(L_*)<\frac9{128}+\frac{1323}{16000}+\frac{66}{1400}
                  <\frac9{40}<R(r).
\]
On \([18/25,73/100]\) the bounds are
\(w(r)\ge w(73/100)>(939/1000)\overline W\) and
\(R(r)\ge R(18/25)>243/1000\). Thus \(L_*<61/1000\).
Using \(e>27/10\) and \(\sqrt{61/1000}<247/1000\),
in addition to the same bounds for \(M,\sqrt3,\pi\),
gives
\[
       B(L_*)<\frac{49261973}{204120000}
                       <\frac{243}{1000}<R(r).
\]
This contradicts (NS4) throughout the asserted range.

If simplex-optimal correlations approach the tie from the
right, the fixed-partition deficit \(S_s(T)-S_s(A)\) is
nonnegative along that sequence and zero at \(r\). Its
derivative at \(r\) is then nonnegative. The strict slope
inequality proved above excludes the two-label-per-circle
condition for \(A\). For an actual maximizer, absence of
a common strict-winner circle implies that condition
since the ordinary cells agree with the strict-winner sets
up to null sets.

## Rational bounds at the endpoints

Use the exact simplex value
\[
             F(r)=\frac{3}{4\pi^2}
                \big[\arccos(-r)^2-\arccos(r/2)^2\big].
\]
It is three times the centered sector probability in
Appendix 1 of de Klerk, Pasechnik and Warners, minus \(1/3\).
The form stated there for negative correlations extends to
\(|r|<1\) by real analyticity of both the Hermite series and
the formula.
[Primary formula, Appendix 1, p. 288](https://repository.tilburguniversity.edu/bitstreams/66666d4a-dfaa-4d38-852a-03cbd8bb034a/download)
[Canonical simplex formula and derivative, D0--D2](THREE_CELL_DEFICIT_AND_PERIMETER.md#exact-correlation-derivative-for-the-three-cell-simplex)

Set \(U=\arcsin r\), \(V=\arcsin(r/2)\),
\(C=(1-r^2)^{-1/2}\), \(D=(1-r^2/4)^{-1/2}\), and
\(a=\pi/2+U,\ b=\pi/2-V\), with lower and upper
bounds denoted below by subscripts \(l,u\). The following
interval endpoints are in units of \(10^{-6}\).
\[
\begin{array}{c|rrrr}
 r & U & V & C & D\\ \hline
 7/10 &[775397,775405]&[357571,357572]&[1400278,1400290]&[1067521,1067522]\\
 18/25 &[803802,803824]&[368267,368268]&[1440972,1441002]&[1071866,1071867]\\
 73/100 &[818321,818332]&[373632,373633]&[1463169,1463183]&[1074105,1074106]
\end{array}
\]
These rational bounds follow by summing terms indexed by
\(0\le j\le N\), with \(N=16\) for the first two rows and
\(N=18\) for the last, of
\[
 \arcsin x=\sum_{j\ge0}\frac{\binom{2j}{j}}{4^j(2j+1)}x^{2j+1},
 \qquad
 (1-x^2)^{-1/2}=\sum_{j\ge0}\frac{\binom{2j}{j}}{4^j}x^{2j}.
\]
All coefficients are in \([0,1]\). The respective remainders
are between zero and \(x^{2N+3}/(1-x^2)\),
\(x^{2N+2}/(1-x^2)\). Use \(\pi_l=333/106<\pi<355/113=\pi_u\),
which also follows from Machin's identity and the alternating
arctangent series through degree nine at \(1/5\) and degree
three at \(1/239\).

To detail the last substitutions, write \(a=\pi/2+U\) and
\(b=\pi/2-V\), so \(S'_r(T)=3(2aC+bD)/(4\pi^2)\).
Define \(a_l=\pi_l/2+U_l,\ a_u=\pi_u/2+U_u,\)
\(b_l=\pi_l/2-V_u,\ b_u=\pi_u/2-V_l\).
A lower bound on \((4\pi^2/3)(w(r)-k\overline W)\) is
\[
            \frac2r(a_l^2-b_u^2)-2a_u C_u-b_u D_u
                          -k\frac{349}{300}\pi_u .
\]
For \((r,k)=(7/10,1)\) this rational number is
\(46966261271047523/53802181500000000000>0\).
For \((r,k)=(18/25,24/25)\) it is
\(365208685002674513/80703272250000000000>0\).
For \((r,k)=(73/100,939/1000)\) it is
\(14536279161076607/8182415103125000000>0\).
The table gives the lower bound
\[
 R(r)\ \ge\ \frac3{4\pi_u^2}
      \left[\frac{2a_l C_l+b_l D_l}{r}
                      -\frac{a_u^2-b_l^2}{r^2}\right].
\]
The bracket is positive at both evaluations.
At \(r=7/10\) the bound is
\(158510923963996304769/693848281000000000000>9/40\).
At \(r=18/25\) it is
\(931656359960267687/3823245630000000000>243/1000\).

## A uniform lower bound on simultaneous angular mass

The slope obstruction is stable under a small amount of
third-label angular mass. On each compact interval
\(J\subset(0,73/100]\) there exists \(\epsilon_J>0\) such
that every ordinary balanced \(A\) with
\(S_r(A)=S_r(T)\) and \(S'_r(A)\le S'_r(T)\), \(r\in J\),
satisfies
\[
                      \int\min_i p_i(s)\,d\mu(s)\ge\epsilon_J .
                                                               \tag{NS5}
\]
Consequently the set of radii where all three labels
occupy at least \(\epsilon_J/2\) of the circle has
radial measure at least \(3\epsilon_J/2\), since
\(\min_i p_i\le1/3\). This statement does not require
stationarity or compactness of actual maximizers.

Write \(\eta=\int\min_i p_i\,d\mu\). On each circle assign
the smallest label's angular slice to one of the other
labels, breaking ties by a fixed ordering. The resulting
partition is binary on almost every circle and differs
from \(A\) on a set of Gaussian mass \(\eta\).
Its masses \(\mu_i'\) satisfy
\(|\mu_i'-1/3|\le\eta\). Hereafter assume
\(0<\eta<1/12\), so \(\mu_i'\ge1/4\).
The case \(\eta=0\) is already excluded, and larger
\(\eta\) needs no small-distance argument.

Choose a measurable set of radii \(E\) of radial mass
\(t=6\eta\) on which the integrated mass of each label
is within \(\eta/4\) of \(t\mu_i'\).
For existence approximate the three radial fractions
in \(L^1(\mu)\) by simple functions on a common finite
measurable partition, with each error below \(\eta/8\).
In each atom take a subset of relative measure \(t\)
by nonatomicity. The approximate label masses are
then exactly proportional, and the stated estimates
follow by the two \(L^1\) errors.
Every actual label mass on \(E\) exceeds its possible
surplus \(\mu_i'-1/3\), since
\(t\mu_i'-\eta/4\ge5\eta/4\).
Repaint all radii in \(E\) as pure circles, assigning
their radial masses so each label finishes at \(1/3\).
The required three masses on \(E\) are nonnegative
and sum to \(t\), so nonatomicity gives such an
assignment. The final balanced binary-slice partition
\(A_{\rm bin}\) differs from \(A\) on Gaussian mass at most
\(7\eta\).

For \(j=1,2\), contraction of the orthogonal Hermite
projection and the fixed indicator norms give
\[
                  |W_j(A)-W_j(A_{\rm bin})|
                      \le 2\sqrt{14\eta}.
\]
Indeed the squared \(L^2\) distance between vector
indicators is twice the changed Gaussian mass.
At \(A\) both lower bounds in (NS4) hold by the tie
and slope assumptions. At \(A_{\rm bin}\), (NS2) and (NS3)
bound \(W_1\) by \(X(L)=\overline W(1-L)\) and \(W_2\)
by \(B(L)\) for some \(L\in[0,1]\).
The strict comparisons in the slope proof show that
at each \((r,L)\in J\times[0,1]\) at least one of
\[
  w(r)-X(L),\qquad
  3F(r)/r-S'_r(T)-2X(L)-r B(L)
\]
is positive. Their maximum is continuous and therefore
has a strictly positive minimum on this compact set.
The Hermite errors above are incompatible with this
minimum when \(\eta\) is sufficiently small. This
proves (NS5).
