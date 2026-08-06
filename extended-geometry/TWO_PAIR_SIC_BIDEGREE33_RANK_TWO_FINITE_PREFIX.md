# The isolated rank-two cubic finite-prefix component

## 1. The exact point

Work on the normalized anti-Weyl chart

\[
s_4=-s_2,\quad s_5=s_1,\quad s_6=-s_0,\qquad
t_2=0,\quad t_3=t_1,\quad t_4=-t_0.
\tag{1.1}
\]

All odd moments vanish identically.  On this chart coefficient rank at
most two is cut out by

\[
\begin{aligned}
R={}&-9s_0s_2+9s_1^2-9s_2^2-9s_3^2-8s_3-t_0^2+4t_1^2+1,\\
I={}&-9s_0s_3+s_0+18s_1s_2+9s_2s_3+9s_2-4t_0t_1.
\end{aligned}                                             \tag{1.2}
\]

Let \(c\) be the rational vector represented by the terminating decimals

\[
(-0.8702368803753021,-1.8595344299498098,
-0.3028762025151954,1.5192025917402650,
-2.3287623389711780,-1.5350489449499409)
\tag{1.3}
\]

in the coordinate order \((s_0,s_1,s_2,s_3,t_0,t_1)\), and put

\[
X=\prod_{j=1}^6[c_j-10^{-10},c_j+10^{-10}].              \tag{1.4}
\]

The exact representation of the point is the system

\[
R=I=\mu_2=\mu_4=\mu_6=\mu_8=0                           \tag{1.5}
\]

together with the rational isolating box \(X\).  No floating-point root
acceptance is part of the proof.

## 2. Certificate

The checker constructs the primitive anti-Weyl moments directly from the
contraction formula and performs a Krawczyk test using rational interval
arithmetic.  If \(J\) is the Jacobian of (1.5), it takes the exact rational
matrix \(A=J(c)^{-1}\) and verifies

\[
c-Af(c)+(1-AJ(X))(X-c)\subset\operatorname{int}(X).       \tag{2.1}
\]

Thus (1.5) has a unique real zero \(\xi\) in \(X\).

In the row-major coefficient matrix, all sixteen \(3\)-by-\(3\) minors
reduce to zero modulo \((R,I)\).  The minor on rows and columns
\(\{1,2\}\) is bounded away from zero on \(X\).  Hence

\[
\operatorname{rank}C(\xi)=2.                              \tag{2.2}
\]

On the smooth twelve-dimensional determinantal chart using the same
minor, rational interval evaluation exhibits a nonzero \(8\)-by-\(8\)
minor of the differential of \((\mu_1,\ldots,\mu_8)\).  Its rank is at
least eight.  At a moment-zero point, the three infinitesimal diagonal
\(\mathrm{SL}_2\) directions and radial scaling lie in the kernel.  They
are independent here: the normalized non-null
\(\operatorname{Sym}^2\) component has one-dimensional stabilizer, and
the nonzero-weight sextic coordinates break that stabilizer; radial
scaling cannot be an infinitesimal \(\mathrm{SL}_2\) direction because
the quadratic discriminant is nonzero.  Therefore the rank is exactly
eight and these four directions exhaust the tangent kernel.  The induced
projective GIT point is reduced and isolated.

The same interval evaluation gives, after removal of the positive integer
contents of the moments,

\[
\boxed{\mu_{10}(\xi)>0,\qquad \mu_{12}(\xi)<0,\qquad
\mu_{14}(\xi)>0.}                                        \tag{2.3}
\]

The non-null normalized \(\operatorname{Sym}^2\) discriminant is a
nonzero invariant, so \(\xi\) is semistable.  We have proved:

> **Theorem 2.1.** The rank-at-most-two nine-moment fiber contains an
> exact-rank-two semistable reduced isolated projective point.  Its first
> missing moment is already nonzero.  Consequently this component does
> not lift to either \(Z(\mu_1,\ldots,\mu_8,\mu_{12})\) on the rank-two
> locus or \(Z(\mu_1,\ldots,\mu_{12},\mu_{14})\) in the ambient space.

This realizes the semistable point previously forced only by the
rank-stratified Hilbert-series calculation.  It is a finite-prefix
component, not an all-order moment-zero point.

## 3. Smaller corrected-component charts

For the remaining anti-Weyl classification set

\[
x=s_1,\quad y=s_3,\quad A=s_0^2,\quad B=s_2^2,\quad
C=t_0^2,\quad D=t_1^2,\quad p=s_0s_2.
\]

The equations \(\mu_2=R=0\) give

\[
\begin{aligned}
A={}&-15B+\frac{14}{3}C+\frac{56}{3}D-6x^2-10y^2-\frac{70}{3},\\
p={}&-B-\frac19C+\frac49D+x^2-y^2-\frac89y+\frac19.
\end{aligned}                                             \tag{3.1}
\]

The remaining rank equations are

\[
p^2=AB                                                     \tag{3.2}
\]

and

\[
16CD=(1-9y)^2A+81(2x+y+1)^2B
 +18(1-9y)(2x+y+1)p.                                     \tag{3.3}
\]

After (3.1), the primitive corrected moments have the following sparse
profiles in \((x,y,B,C,D)\):

| moment | terms | degree |
| --- | ---: | ---: |
| \(\mu_4\) | 39 | 4 |
| \(\mu_6\) | 119 | 6 |
| \(\mu_8\) | 294 | 8 |
| \(\mu_{12}\) | 1218 | 12 |

Exact characteristic-zero Gröbner calculation is now immediate on this
quotient:

\[
\bigl(p^2-AB,\ (3.3),\ \mu_4,\mu_6,\mu_8,\mu_{12}\bigr)
=\mathbb Q[x,y,B,C,D].                                   \tag{3.4}
\]

The `msolve` replay reduces 955 pairs and returns the one-element basis
\([1]\).  Because the quotient equations are necessary for every original
anti-Weyl point, possible squaring artifacts do not weaken the exclusion.

> **Proposition 3.1.** The corrected exact-rank-two system
> \(Z(\mu_1,\ldots,\mu_8,\mu_{12})\) has no point on the normalized
> anti-Weyl chart.

Off the exceptional Hurwitz boundary, a rank-two cubic pencil can instead
be normalized to

\[
A_1=X^2(X+T),\qquad A_2=T^2(\lambda X+T),                 \tag{3.5}
\]

with an arbitrary \(2\)-by-\(4\) channel matrix.  Projective normalization
and elimination of \(\mu_1\) leave seven variables.  The coincident-critical-
point and coincident-critical-value divisors must be treated separately;
(3.5) is not a boundary classification.  This generic seven-variable chart
and its exceptional divisors are the remaining rank-two component-
extraction frontier.  Section 4 closes four complete parameter fibres
and one doubly exceptional elimination branch inside this chart.

## 4. Exact exclusions on the generic Hurwitz chart

Dehomogenize (3.5), normalize the first channel coefficient, and write

\[
\begin{aligned}
F(x,y)&=(1+x)B(y)+x^2(\lambda+x)D(y),\\
B(y)&=1+a_1y+a_2y^2+a_3y^3,\\
D(y)&=b_0+b_1y+b_2y^2+d_3y^3.
\end{aligned}                                             \tag{4.1}
\]

The first moment is linear with nonzero constant pivot and gives

\[
d_3=-1-\frac{a_1}{3}-\frac{\lambda b_2}{3}.              \tag{4.2}
\]

The calculation below is localized at the nonzero quadratic
discriminant

\[
\begin{aligned}
\Delta={}&(9+2a_1+\lambda b_2)^2\\
 &+(3+2\lambda b_1+3b_2)
 \bigl((3-\lambda)a_1+2a_2-3\lambda-\lambda^2b_2\bigr)
\end{aligned}                                             \tag{4.3}
\]

and at the channel minor

\[
M_{01}=b_1-a_1b_0.                                       \tag{4.4}
\]

These two localizations are part of the statement: neither their zero
divisors nor the other projective channel charts are silently included.

There are two characteristic-zero unit calculations.

> **Proposition 4.1 (the complete \(\lambda=0\) fibre).**  In
> \(\mathbb Q[a_1,a_2,a_3,b_0,b_1,b_2,\Delta^{-1},M_{01}^{-1}]\), the
> moments \(\mu_2,\ldots,\mu_8\) at \(\lambda=0\) generate the unit
> ideal.

Thus the entire declared \(\lambda=0\) fibre of the generic Hurwitz
chart contains no all-order pure-moment point.

The coefficient of \(a_3\) in the primitive second moment is

\[
P_1=3(3b_0\lambda+3b_0+4b_1).                            \tag{4.5}
\]

On \(P_1=0\), substitute

\[
b_1=-\frac34b_0(\lambda+1).                              \tag{4.6}
\]

The remaining second moment is linear in \(a_2\).  Its primitive pivot
is

\[
P_2=-\frac34
 \bigl(9b_0\lambda^2+2b_0\lambda+9b_0-16b_2-16\bigr).
\tag{4.7}
\]

On the secondary boundary \(P_2=0\), substitute

\[
b_2=\frac{9b_0\lambda^2+2b_0\lambda+9b_0-16}{16}.        \tag{4.8}
\]

Exact rational elimination after (4.6)--(4.8), retaining the equation
\(\mu_2=0\) and the open condition \(\Delta M_{01}\ne0\), gives

\[
\bigl(\mu_2,\mu_3,\ldots,\mu_8\bigr)=(1).                \tag{4.9}
\]

> **Proposition 4.2 (doubly exceptional pivot branch).**  The locus
> \(P_1=P_2=0\) in the declared generic Hurwitz chart has no common zero
> of \(\mu_2,\ldots,\mu_8\), and hence no all-order pure-moment point.

The same exact calculation can be run after specializing the Hurwitz
parameter, without imposing either pivot equation.  At each of

\[
\lambda=-1,\qquad \lambda=1,\qquad \lambda=2,             \tag{4.10}
\]

the localized ideal generated by \(\mu_2,\ldots,\mu_8\) is the unit
ideal over \(\mathbb Q\).

> **Proposition 4.3 (three further complete fibres).**  The complete
> declared fibres \(\lambda=-1,1,2\) of the generic Hurwitz chart contain
> no common zero of \(\mu_2,\ldots,\mu_8\), and hence no all-order
> pure-moment point.

A separate five-parameter slice retains \(\lambda\) and does not require
either localization.  Impose

\[
b_0=0,\qquad b_1=1,                                    \tag{4.11}
\]

so that

\[
 B=1+a_1y+a_2y^2+a_3y^3,\qquad
 D=y+b_2y^2+\left(-1-\frac{a_1}{3}
                       -\frac{\lambda b_2}{3}\right)y^3.  \tag{4.12}
\]

The constant terms of \(B,D\) and the disjoint degree supports of the
two first-factor cubics show that the coefficient matrix has rank exactly
two everywhere on this slice.  Exact elimination over \(\mathbb Q\)
gives a zero-dimensional quotient of length \(687\) for
\((\mu_2,\ldots,\mu_6)\).  Adjoining \(\mu_7\) leaves support equal to
the single rational point

\[
 (\lambda,a_1,a_2,a_3,b_2)=(1,-4,5,-2,-2),             \tag{4.13}
\]

with local length \(26\).  At this point \(d_3=1\), and the affine form is

\[
 F=(1+x)(1-y)^2(1-2y+x^2y).                            \tag{4.14}
\]

Under the contraction-preserving flag change

\[
 W=W',\quad V=V'-W',\quad Z=Z'+Y',\quad Y=Y',
\]

its homogeneous form becomes

\[
 V'Z'^2\bigl(V'^2Y'-2V'W'Y'+W'^2Z'\bigr).             \tag{4.15}
\]

The three nonzero coefficient positions are
\((1,0),(2,1),(3,1)\), all strictly in the chamber \(i>j\).  Therefore
the relative period has strictly negative \(u\)-valuation,

\[
 P(u,t)=u^{-1}(1-t)^3-2u^{-1}t(1-t)^2
        +u^{-2}t(1-t)^2,                               \tag{4.16}
\]

and the exact creative-telescoping recurrence is simply
\(\nu_{m+1}=0\).  The degree-one multiplier with period numerator \(ut\)
has raw mixed value \(-2\) at \(m=1\), so the mixed sequence is nonzero,
but every degree-\(e\) multiplier vanishes for \(m>e\).

> **Proposition 4.4 (simple-root coefficient slice).**  The point (4.13)
> is the only all-order pure-moment point on the exact-rank-two slice
> (4.11)--(4.12), and it is fixed-flag one-sided and SIC-safe.  This is a
> coefficient-slice theorem, not an independent channel-side
> \(\mathrm{SL}_2\) normalization or a diagonal-\(\mathrm{SL}_2\) orbit
> classification.

All calculations in Propositions 4.1--4.4 construct the moments of \(3F\)
over the \(521\)-bit Mersenne prime \(2^{521}-1\).  Every integer coefficient of
the order-\(m\) contraction is bounded in absolute value by
\((3m)!\,52^m\), and the prime is larger than twice this bound through
order eight.  Symmetric representatives therefore recover the exact
integer polynomials, rather than a heuristic modular image.  Scaling
\(F\) by \(3\) does not change their zero sets.  `msolve` over
\(\mathbb Q\) returns the one-element basis \([-1]\) in each exclusion
calculation.

The square system \(\mu_2,\ldots,\mu_8\) in the original seven
parameters has mixed volume \(74\,144\).  This explains the size of the
unreduced elimination, but is only a numerical-algebraic complexity
measurement.  Generic finite-field elimination, native saturation, and
the reduced branches below exceeded the recorded memory limits; those
failures have no mathematical force.

Away from the four closed fibres \(\lambda\in\{-1,0,1,2\}\), the two
generic pivot pieces still remain.  On the tertiary open inside
\(P_1=0,\ P_2\ne0\), eliminating \(a_2\)
and then \(a_3\) reveals exact powers \(P_2^2,P_2^3,P_2^3,P_2^4\) in
the cleared moments \(\mu_4,\ldots,\mu_7\), respectively.  Cancelling
them leaves \(4\,493,11\,332,17\,566,40\,418\) terms over
\(\mathbb F_{1\,073\,741\,827}\).  Even this four-parameter square
system produced no basis under either `msolve`'s exact-sparse or
probabilistic-sparse kernel with a 5 GB address-space limit; the native
saturation kernel also terminated without a basis.  These figures record
the current computational ceiling only.  In particular, they are neither
a finite-field exclusion nor evidence that the branch is nonempty.

The surviving pieces of this one channel chart are now explicit:

1. \(P_1\ne0\), where \(\mu_2\) eliminates \(a_3\);
2. \(P_1=0,\ P_2\ne0\), where \(\mu_2\) eliminates \(a_2\); and
3. the zero divisors of \(\Delta M_{01}\).

On the first piece, eliminate \(a_3\) with \(\mu_2\) and write the
primitive numerator of the reduced third moment as

\[
 Q=Aa_2^2+Ba_2+C,
\]

where

\[
 A=3b_0^2\lambda^2+b_0^2\lambda+6b_0^2
   +b_0b_1\lambda+5b_0b_1-8b_0b_2-8b_0+4b_1^2.       \tag{4.17}
\]

After inverting \(P_1A\), the quotient by \(Q\) is a free rank-two
module with basis \(1,a_2\).  Thus every later reduced moment has a
unique remainder

\[
 \operatorname{rem}_Q(\mu_m)=U_m+V_ma_2.              \tag{4.18}
\]

This gives an exact nonexpanded split.  On \(V_4\ne0\), the common-zero
conditions are

\[
 AU_4^2-BU_4V_4+CV_4^2=0,
 \qquad U_mV_4-V_mU_4=0\quad(m\ge5).                  \tag{4.19}
\]

On \(V_4=0\), one instead imposes \(U_4=0\) and retains the two module
coordinates of each later moment.  The divisor \(A=0\) splits further
into the linear case \(B\ne0\) and the degeneracy \(A=B=0\), with
\(C=0\) required for a root when the quadratic collapses completely.
This module presentation is the next exact component computation; it
does not assert that any of these branches is empty.  Notice that the
slice (4.11) has \(A=4\), which explains why its characteristic-zero
elimination remains tractable.

The linear divisor itself has one more exact triangular reduction.  At
\(b_0=0\), equation (4.17) gives \(A=4b_1^2\), while
\(M_{01}=b_1-a_1b_0=b_1\).  Thus \(A=0\) on the declared minor open
forces \(b_0\ne0\), and \(A=0\) eliminates

\[
 b_2=\frac{3b_0^2\lambda^2+b_0^2\lambda+6b_0^2
       +b_0b_1\lambda+5b_0b_1-8b_0+4b_1^2}{8b_0}.     \tag{4.20}
\]

Put

\[
 q=\frac{b_1}{b_0},\qquad z=b_0a_2.                    \tag{4.21}
\]

Then (4.20) becomes

\[
 b_2=-1+\frac{b_0}{8}
 \bigl(3\lambda^2+\lambda+6+q\lambda+5q+4q^2\bigr).  \tag{4.22}
\]

Retaining \(z\) and the sparse incidence equation
\(\mathcal Bz+\mathcal C=0\), rather than expanding
\(z=-\mathcal C/\mathcal B\), gives an exact \(\mathbb Q\)-system with
six variables and six equations through \(\mu_7\), including the inverse
localizer.  The ratio chart exposes removable powers
\(b_0,b_0^3,b_0^4,b_0^5,b_0^6\) in the cleared moments
\(\mu_3,\ldots,\mu_7\).  After cancelling them, the incidence equation
has 273 terms and the reduced \(\mu_4,\ldots,\mu_7\) have respectively
\(915,2245,4234,8687\) terms.  The complete `msolve` source is 593,624
bytes.

The first component layer can now be certified without computing that
whole standard basis.  On \(\mathcal B\ne0\), eliminate \(z\) between the
incidence and \(\mu_4\).  The primitive exact numerator \(R_4\) has
11,128 terms, total degree \(44\), and degrees

\[
 \deg_{(\lambda,a_1,b_0,q)}R_4=(32,10,12,28).          \tag{4.23}
\]

Its degree-preserving reduction modulo \(29\) has 10,747 terms and is a
nonzero scalar multiple of the polynomial independently constructed from
the modular moments.  Singular factors it as one degree-\(44\) factor.
Therefore \(R_4\) is irreducible over \(\mathbb Q\).  On the declared
\(\mathcal B\)-open, the incidence-plus-\(\mu_4\) quotient is obtained by
adjoining \(z=-\mathcal C/\mathcal B\), so its nonempty localization is
irreducible.

On the complementary degree-drop locus
\(\mathcal B=\mathcal C=0\), the coefficient \(\mathcal B\) is quadratic
in \(a_1\), with a nonzero constant leading coefficient.  Write

\[
 \mathcal B=A_2a_1^2+A_1a_1+A_0,\qquad
 \mathcal C\equiv U+Va_1\pmod{\mathcal B}.             \tag{4.24}
\]

The exact coefficient term counts of \((A_0,A_1,A_2)\) are
\((49,12,1)\), while \(U,V\) have \(160,70\) terms.  The quadratic norm is

\[
 A_2U^2-A_1UV+A_0V^2.                                  \tag{4.25}
\]

Denominator clearing contributes exactly
\((3\lambda+3+4q)^2=(P_1/(3b_0))^2\).  This factor is invertible on the
present chart.  After removing it, the exact norm \(H\) has 714 terms,
total degree \(26\), and degrees

\[
 \deg_{(\lambda,b_0,q)}H=(18,8,18).                    \tag{4.26}
\]

Its degree-preserving reduction modulo \(29\) is irreducible of degree
\(26\), so \(H\) is irreducible over \(\mathbb Q\).  The \(V\ne0\) part
of (4.24) is birational to \(H=0\).  Since multiplication by
\(\mathcal C\) in the free quadratic algebra defined by \(\mathcal B\)
has the nonzero norm (4.25), \(\mathcal B,\mathcal C\) are a regular
sequence.  The complete intersection is unmixed, and the irreducible
norm with its generically linear gcd leaves one irreducible component.

The first later-moment equation on this component is also exact.  If
\(\mu_4\equiv X_4+Y_4a_1\pmod{\mathcal B}\), it is

\[
 N_4=VX_4-UY_4.                                        \tag{4.27}
\]

After removing \(P_1^3\), \(N_4\) has 1,296 exact terms, total degree
\(27\), \(z\)-degree \(3\), and degrees
\((18,3,9,18)\) in \((\lambda,z,b_0,q)\).  Its exact primitive reduction
modulo \(29\) has 1,252 terms and one degree-\(27\) irreducible factor.
Thus \(N_4\) itself is irreducible over \(\mathbb Q\).

The same construction starts a subresultant ladder instead of expanding a
four-variable standard basis.  Let

\[
 L_4=\operatorname{lc}_z(N_4),\qquad
 L_5=\operatorname{lc}_z(N_5),
\]

where the descended \(N_5\) has first had its exact \(P_1^4\) content
removed.  On \(L_4\ne0\), the two cubic equations \(N_4=N_5=0\) are
equivalent to \(N_4=S_5=0\), where

\[
 S_5=L_4N_5-L_5N_4.                                    \tag{4.28}
\]

The \(z^3\) term cancels.  The exact \(S_5\) has 4,332 terms, \(z\)-degree
two, total degree \(40\), and degrees
\((27,2,13,27)\) in \((\lambda,z,b_0,q)\).  Its
degree-preserving reduction modulo \(29\) has 4,195 terms and one
degree-\(40\) irreducible factor.  The exceptional divisor \(L_4=0\) is
itself irreducible: \(L_4\) has 50 exact terms, degree \(9\), and degrees
\((6,3,6)\) in \((\lambda,b_0,q)\).  Finally
\(\operatorname{lc}_z(S_5)\) has 825 exact terms, degree \(28\), degrees
\((19,9,19)\), and irreducible degree-preserving reduction modulo \(29\).
Thus the next generic step is a cubic--quadratic pseudo-remainder, with
the two displayed leading-coefficient divisors retained as separate exact
branches.

> **Proposition 4.5 (first exact module-descent layer).**  On
> \(P_1\ne0,\ A=0,\ b_0\Delta M_{01}\ne0\), the
> \(\mathcal B\ne0\) incidence-plus-\(\mu_4\) projection is irreducible.
> The complementary locus \(\mathcal B=\mathcal C=0\) has one
> irreducible component with projection \(H=0\), and its descended
> \(\mu_4\) equation \(N_4\) is irreducible.  On \(L_4\ne0\), the
> descended \(\mu_5\) condition is the irreducible quadratic-in-\(z\)
> subresultant \(S_5=0\); the leading-coefficient boundaries \(L_4=0\)
> and \(\operatorname{lc}_z(S_5)=0\) are irreducible divisors.  These are
> component and descent statements, not an assertion that the later
> common zero set is empty.

That forced cubic--quadratic step can also be completed without forming
the ambient ideal.  Write

\[
 S_5=S_2z^2+S_1z+S_0,
 \qquad
 S_2^2N_4\equiv R_0+R_1z\pmod {S_5}.                 \tag{4.29}
\]

The exact pseudo-remainder coefficients \(R_0,R_1\) have respectively
16,668 and 13,475 terms, total degrees \(83\) and \(77\), and
multidegrees

\[
 \deg_{(b_0,q,\lambda)}R_0=(27,56,56),\qquad
 \deg_{(b_0,q,\lambda)}R_1=(25,52,52).               \tag{4.30}
\]

Their quadratic norm satisfies

\[
 S_2R_0^2-S_1R_0R_1+S_0R_1^2
     =b_0^2L_4^2S_2^2K.                              \tag{4.31}
\]

The primitive exact residual \(K\) has 48,469 terms, total degree \(118\),
and multidegree \((37,81,81)\) in \((b_0,q,\lambda)\).  Its
degree-preserving reduction modulo \(29\) has 46,755 terms and one
degree-\(118\) factor, so \(K\) is irreducible over \(\mathbb Q\).
On \(b_0L_4S_2\ne0\), equation \(K=0\) is exactly the projection of the
common-root condition \(N_4=S_5=0\).

The same two-dimensional quotient module absorbs the next moments.  For
\(j=6,7\), first remove the exact \(P_1^{j-2}\) content from the descended
moment \(N_j\), and define \(X_j,Y_j\) by the lead-cleared remainder

\[
 S_2^{\deg_zN_j-1}N_j\equiv X_j+Y_jz\pmod {S_5}.
\]

Pairing this remainder with (4.29) gives the base equation

\[
 D_j=R_1X_j-R_0Y_j.                                  \tag{4.32}
\]

Exact division gives, up to a nonzero rational scalar,

\[
 D_6\doteq b_0L_4S_2^2J_6,
 \qquad
 D_7\doteq b_0^2L_4S_2^2J_7.                        \tag{4.33}
\]

The primitive polynomials \(J_6,J_7\in\mathbb Z[b_0,q,\lambda]\) have
profiles

\[
\begin{array}{c|c|c|c}
 &\text{terms}&\text{total degree}&
   \deg_{(b_0,q,\lambda)}\\ \hline
J_6&63{,}791&130&(41,89,89)\\
J_7&121{,}758&162&(51,111,111).
\end{array}                                           \tag{4.34}
\]

Their independently constructed reductions modulo \(29\) have 61,552
and 117,609 terms.  Singular factors each as one factor of the full
degree, proving both \(J_6\) and \(J_7\) irreducible over \(\mathbb Q\).
On \(VR_1\ne0\), \(z=-R_0/R_1\), so after \(N_4=S_5=0\) the equations
\(J_6=0\) and \(J_7=0\) are equivalent to the descended
\(\mu_6,\mu_7\) conditions.  If \(R_1=0\), equation \(K=0\) on the same
leading-coefficient open forces \(R_0=0\); this simultaneous boundary
must be treated separately rather than divided away.

> **Proposition 4.6 (exact module descent through \(\mu_7\)).**  On the
> open of Proposition 4.5 with \(L_4S_2\ne0\), the projection of
> \(N_4=S_5=0\) is the irreducible degree-\(118\) hypersurface \(K=0\).
> On the further open \(VR_1\ne0\), imposing the descended moments
> \(\mu_6,\mu_7\) is equivalent to adjoining the irreducible base
> equations \(J_6=J_7=0\), of degrees \(130\) and \(162\).  Thus the
> generic finite prefix through \(\mu_7\) has been reduced exactly to
> \(H=K=J_6=J_7=0\) in \((b_0,q,\lambda)\).  This proposition neither
> proves that their common zero set is empty nor covers \(V=0\),
> \(R_0=R_1=0\), or a leading-coefficient boundary.

The direct system \(H,N_4,N_5,N_6\) over \(\mathbb F_{29}\) still exceeded
2.7 GB before producing a basis.  That stopped calculation has no
mathematical force.  The certified descent shows that future work should
continue by quotient-module subresultants and exceptional-divisor splits,
not by asking one ambient basis computation to discover the components.

A first three-variable modular basis request for
\(H,K,J_6,J_7\subset\mathbb F_{29}[b_0,q,\lambda]\) was also stopped after
reaching 2.53 GB, and a coefficient-field resultant
\(\operatorname{Res}_\lambda(H,K)\) was stopped after ten minutes.  These
stopped computations likewise have no mathematical force.  Reduction
modulo \(H\) over \(\mathbb F_{29}(b_0,q)\), however, is cheap: \(H\) has
\(\lambda\)-degree \(18\), while \(K\) reduces to a nonzero remainder of
\(\lambda\)-degree \(17\).  The next exact target is therefore a
fraction-free norm or subresultant divisor in \((b_0,q)\), followed by
intersection with \(J_6,J_7\), rather than another simultaneous basis.

A bounded modular fibre scout supplies a more focused routing signal.
It constructs the same quotient-module equations through \(J_9\), then
for every \(q\in\mathbb F_p\) adjoins an inverse to

\[
 b_0P_1L_4S_2VR_1.                                   \tag{4.35}
\]

Over \(\mathbb F_{29}\), 28 of the 29 rational \(q\)-fibres are unit;
the sole survivor is \(q=22\), with the length-two basis

\[
 \lambda+3s-1,\qquad b_0-s-7,\qquad s^2+s-7=(s-14)^2.
\]

Over \(\mathbb F_{31}\), 27 of 31 fibres are unit and four finite
schemes survive, at \(q=1,12,26,29\).  In both primes the survivor bases
are unchanged after adjoining \(J_8,J_9\).  The discovered factor-removal
profiles for \(J_6,J_7,J_8,J_9\) are respectively

\[
 (P_1,b_0,L_4,S_2)=(4,1,1,2),(5,2,1,2),
 (5,2,1,2),(6,2,1,2).                                 \tag{4.36}
\]

These are exact bounded computations over the two finite fields only.
They neither cover extension-field values of \(q\) nor lift a unit ideal
to characteristic zero.  The tempting small lift of the \(p=29\) base
point is not a rational moment-zero point: after reconstructing \(a_3\)
from \(\mu_2\), the exact value of \(\mu_3\) is nonzero.  The scout is
therefore a guide to an evaluation--interpolation eliminant, not an exact
exclusion or counterexample.

The reusable pattern is:

1. quotient by the lowest-degree monic or constant-leading pivot and keep
   its standard basis rather than substituting a rational root;
2. express the next equation by multiplication coordinates, a norm, or a
   first subresultant;
3. split the leading coefficient before the next pseudo-division;
4. remove only factors already inverted on the current chart; and
5. certify an exact factor by comparing its primitive
   characteristic-zero polynomial with a degree-preserving irreducible
   reduction at a good prime.

Each generic step lowers the active polynomial degree by one.  Every
failure of a leading coefficient becomes a named lower-dimensional chart,
so no component is lost to a global saturation.

On the second principal piece, the reduced \(\mu_3\) is linear in
\(a_3\).  These are exact triangular reductions, but their remaining
unspecialized ideals have not been proved unit.  The projective
\(B(0)=0\) chart, the other channel-minor opens, and the exceptional
cubic-pencil divisors also remain separate.  Propositions 4.1--4.4 are
therefore genuine characteristic-zero exclusions, and Propositions
4.5--4.6 are exact component and descent theorems, not a global rank-two
classification.

## 5. Consequence for the Picard--Fuchs work

The fixed integral fiber used in the existing order-eight/order-fourteen
calculation has raw moments

\[
2502,\quad120004752,\quad24426240171840,\quad
13040211749112437760.                                    \tag{5.1}
\]

It is therefore a generic benchmark, not a point of any relevant
finite-prefix zero fiber.  Its exact reconstructed operators, Ore division,
forward-coefficient audit, and fifty-value bridge remain valuable regression
tests.  They do not settle the corrected component locus.

The characteristic-zero gap is specifically the relative boundary identity.
Closed-path Griffiths--Dwork reduction produces Picard--Fuchs operators for
closed periods, while relative-period theory obtains an inhomogeneous
equation from the exact boundary form.  Reduction-based creative telescoping
can intentionally avoid certificates, which is why an exact operator alone
does not identify the interval endpoint trace.  See
[Lairez](https://arxiv.org/abs/1404.5069),
[Bostan--Lairez--Salvy](https://arxiv.org/abs/1301.4313), and
[Li--Lian--Yau](https://arxiv.org/abs/0910.4215).  The recent twisted-period
extension likewise formulates the problem through a Griffiths--Dwork
\(D\)-module
[on relative twisted forms](https://arxiv.org/abs/2604.09129).

For the compact chart, both interval boundaries \(y=0\) and \(y=\infty\)
must be retained.  The next invariant proof object is the
\(14_{\rm interior}+2_{0}+2_{1}\) relative Gauss--Manin connection.  It
should be constructed only after an exact corrected-system component has
survived the missing-moment audit.  At an algebraic survivor the natural
coefficient field is its residue field \(K\); descent to \(\mathbb Q\)
requires a separate argument and restriction of scalars may increase
scalar order.
