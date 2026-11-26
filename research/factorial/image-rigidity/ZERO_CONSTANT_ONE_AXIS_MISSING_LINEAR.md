# Four consecutive factorial moments in one-axis resonance with a missing linear term

Let
\[
 f=x^2/2+\beta xy+\gamma y^2/2+dx+ey,\qquad
 P_i=\Lambda_2(f^i)/i!,
\]
where \(\Lambda_2(x^a y^b)=a!b!\). Work over a valued extension of
\(\mathbb Q_2\), with \(v(2)=1\) and \(v(0)=+\infty\). Suppose
\(\bar\beta=1\) and \(\bar\gamma=\bar d=\bar e=0\).

**Theorem.** If \(d=0\) or \(e=0\), then no four consecutive positive-index
moments \(P_i\) vanish. In particular this gives the support-four window
in the one-axis residue type for either missing linear coefficient.
The proof works in ramified extensions and over arbitrary residue fields
of characteristic two.

## Two exact force equations and the matching congruence

Assume \(P_n=\cdots=P_{n+3}=0\) for \(n\ge1\), and put \(j=n+2\).
Use the axis moments \(A_i,B_i\), axis quotients \(r_i=A_i/A_{i-1}\),
and their differences \(q_i=r_{i+1}-r_i\) from
[unit-pole removal](ZERO_CONSTANT_ONE_AXIS_UNIT_POLE_REMOVAL.md).
In particular \(A_i,r_i\) are units, \(\bar r_i=1\), and
\[
 r_i=(2i-1)+d^2/r_{i-1}\quad(i\ge2),\qquad
 q_1=2-d/(1+d),\qquad
 q_i=2-\frac{d^2 q_{i-1}}{r_{i-1}r_i}\quad(i\ge2).
 \tag{1}
\]
For \(\lambda=v(d)<\infty\), these imply
\[
 \begin{split}
 v(q_i)<1&\ \Longrightarrow\ v(q_i)=(2i-1)\lambda,\\
 (2i-1)\lambda<1&\ \Longrightarrow\ v(q_i)=(2i-1)\lambda,\\
 v(q_i)>1,\quad i\ge2&\ \Longrightarrow\ \lambda=1/(2i-1).
 \end{split}                                                    \tag{2}
\]
For the first line, read (1) backwards to valuation \(v(q_1)<1\).
For the second, apply (1) forward while the displayed valuation stays
below one. For the third, the subtracted term in (1) has valuation one,
so \(v(q_{i-1})=1-2\lambda<1\), and the first line applies.
These arguments include \(q_i=0\) in the third line.

Set \(\rho=(\gamma d-\beta e)/(\gamma-\beta)\) and
\(z_i=(B_i+\rho B_{i-1})/A_{i-1}\). The two consecutive
[interior boundary equations](ZERO_CONSTANT_ONE_AXIS_UNIT_POLE_REMOVAL.md#exact-interior-and-weighted-axis-equations)
give
\[
     (1-\beta)r_i+e-\beta d+(\gamma-\beta)z_i=0
                      \qquad(i=j,j+1).                         \tag{3}
\]
The series \(P(T)=\sum_{i\ge0}P_iT^i\) satisfies
\[
 P(T)\equiv R(T):=
 \frac{1+(\beta+d+e)T+deT^2}
 {(1+T+d^2T^2)(1+\gamma T+e^2T^2)}
                         \pmod{2\mathcal O[[T]]},              \tag{4}
\]
where \(\mathcal O\) is the valuation ring. This is the
[universal matching congruence](ZERO_CONSTANT_BINARY_QUADRATIC_MATCHING_AND_BOUNDARY_LIFTS.md#the-matching-reduction-and-the-unit-pole).
Its error coefficients are \(2\) times integral polynomials in the
natural coefficients. Thus changing each coefficient by valuation at
least \(w\) changes each error coefficient by valuation at least \(1+w\).
Both rational series in such a comparison are expanded at \(T=0\).

Write \(1+T+d^2T^2=(1+aT)(1+uT)\), with \(a\) a unit
and \(v(u)=2v(d)\). The factorization follows by Hensel lifting the
simple residue root, after passage to a completion if needed.
The coefficients of \((1+aT)P(T)\) at indices \(n+1,n+2,n+3\)
vanish. After multiplying (4) by \(1+aT\), consider the remaining pole
parameters, meaning the numbers \(w_\ell\) in factors \(1+w_\ell T\).
If the minimum finite valuation among them is \(s>0\) and the
numerator is \(1+cT\) with \(c\) a unit, the coefficient at \(n+i\)
has valuation at least \((n+i-1)s\). After scaling by a parameter of
value \(s\) to this power, its residue is a unit multiple of the
complete homogeneous coefficient of degree \(n+i-1\) in the scaled
minimum-valuation parameters. This follows by expanding the denominator
factors geometrically. The constant numerator term has higher valuation.
With \(r\) minimum parameters, these homogeneous coefficients cannot
vanish at \(r\) consecutive nonnegative degrees. Their generating
function is the reciprocal of a degree-\(r\) polynomial with unit
leading coefficient, so \(r\) consecutive zeros would propagate back
to its constant coefficient one. This reasoning allows repeated
parameters.

A comparison phase used below is
\[
             f_0=x^2/2+xy+\gamma y^2/2 .
\]
Its series \(R_0\) from (4) equals \(1/(1+\gamma T)\). The
exact interior equation specializes to
\[
    P_i^0=2iP_{i-1}^0+\gamma^i(2i-1)!!,\quad i\ge1,\qquad
    P_0^0=1.                                                   \tag{5}
\]
The formula holds at \(i=1\) by direct evaluation. If
\(0<v(\gamma)<1\), (5) gives \(v(P_i^0)=i\,v(\gamma)\), by
induction and strict dominance of its last term. At \(\gamma=0\)
this comparison phase has \(P_i^0=2^i i!\).

## The missing \(y\)-linear coefficient

Suppose \(e=0\) and first assume \(d\gamma\ne0\).
Write \(\lambda=v(d)\) and \(\sigma=v(\gamma)\). Here
\(B_i=\gamma^i(2i-1)!!\) and \(v(\rho)=\sigma+\lambda\), so
\(v(z_i)=i\sigma\) for every \(i\ge1\). Subtracting the two equations
(3) gives
\[
       v(1-\beta)+v(q_j)=j\sigma .
\]
Since \(v(q_j)>0\), the \(z_j\) term has strictly higher valuation
than \(1-\beta\) in (3). Therefore \(v(1-\beta)=\lambda\) and
\(v(q_j)=j\sigma-\lambda\).

If \(v(q_j)<1\), (2) gives \(\sigma=2\lambda\) and
\((2j-1)\lambda<1\). The two remaining pole parameters \(u,\gamma\)
in (4) both have value \(2\lambda\), and \((n+1)2\lambda<1\).
Hence the two coefficients at \(n+1,n+2\) cannot both vanish after
removing the unit factor as above. If \(v(q_j)>1\), (2) gives
\(\lambda=1/(2j-1)\) and \(\sigma>2\lambda\). The sole
minimum parameter \(u\) gives a nonzero coefficient of valuation
\(2n\lambda<1\) at \(n+1\), again impossible.

Thus \(v(q_j)=1\) and
\[
                     j\sigma=1+\lambda .
\]
By (2), \((2j-1)\lambda\ge1\), so \(\sigma\le2\lambda\).
Equality is excluded by the same two-pole argument since then
\(\lambda=1/(2j-1)\) and \((n+1)2\lambda<1\).
For strict inequality, the force equation at \(j\) implies
\(v((1-\beta)r_j-\beta d)=j\sigma\).
From (1), \(v(r_j-1)\ge\min(1,2\lambda)\). Expanding the preceding
force expression around \(1-\beta-d\) now yields
\[
                      v(\beta+d-1)\ge\min(1+\lambda,2\lambda).
                                                               \tag{6}
\]
The exact rational difference \(R-R_0\) has numerator
\((\beta+d-1)T-d^2T^2\) over
\((1+aT)(1+uT)(1+\gamma T)\). By (4), its error relative to
\(P-P^0\) has coefficient valuations at least \(1+\lambda=j\sigma\).

If \(\lambda\ge1\), (6) shows that \(P-P^0\) itself has coefficient
valuations at least \(j\sigma\). At \(i=n,n+1\), both \(P_i^0\)
would have this valuation bound. The recurrence (5) at \(n+1\) is
then impossible because its forcing has valuation \((n+1)\sigma<j\sigma\).
If \(\lambda<1\), then \(\sigma=(1+\lambda)/j<1\). Multiply the
difference by \(1+aT\). Its remaining denominator parameters have
valuation at least \(\sigma\), and the two numerator valuations are
strictly greater than \(\sigma\), by (6) and \(\sigma<2\lambda\).
Consequently the rational difference at index \(n+1\) has valuation
strictly greater than \(n\sigma\). The transferred matching error does
too. The coefficient \(P_{n+1}^0+aP_n^0\) has valuation exactly
\(n\sigma\) by (5), whereas \(P_{n+1}+aP_n=0\). This is a contradiction.

If \(\gamma=0\) and \(d\ne0\), the two forces at \(e=0\) reduce to
\((1-\beta)r_i-\beta d=0\) for \(i=j,j+1\). They require
\(\beta\ne1\) and \(q_j=0\), so (2) gives \(\lambda=1/(2j-1)\).
The sole small pole \(u\) gives a nonzero coefficient after the unit
factor is removed, of valuation \(2n\lambda<1\) at \(n+1\).
Finally if \(d=e=0\), dividing each force by \(A_i\) gives
\((1-\beta)+(\gamma-\beta)\gamma^i=0\). Two adjacent such equations
are impossible for nonzero \(\gamma\), since \(\gamma-\beta\) and
\(1-\gamma\) are units. When \(\gamma=0\) they require \(\beta=1\),
and (5) gives nonzero moments.

## The missing \(x\)-linear coefficient

Suppose \(d=0\) and \(e\ne0\), and write \(t=v(e)\) and
\(g=v(\gamma)\). Now \(r_i=2i-1,\ q_i=2\), and
\(v(\rho)=t\). If \(g<t\), the axis recurrence gives
\(v(B_i)=ig\) by induction, since the \(e^2 B_{i-2}\) term
has strictly higher valuation than the \(\gamma(2i-1) B_{i-1}\)
term. It follows that \(v(z_i)=ig\). Subtract the two forces (3)
to get \(v(1-\beta)=jg-1\), then use the first force to get
\(v(1-\beta)=t\). This also gives
\[
 \beta+e-1=-2j(1-\beta)-(\gamma-\beta)z_j,\qquad
                         v(\beta+e-1)\ge1+t=jg.                 \tag{7}
\]
The numerator of \(R-R_0\) over the common denominator
\((1+T)(1+\gamma T+e^2T^2)(1+\gamma T)\) is
\[
     (\beta+e-1)(1+\gamma T)T-e^2T^2(1+T).
\]
The error in the comparison with \(P-P^0\) has valuation at least
\(1+t=jg\) by (4). If \(t\ge1\), the displayed numerator also has
this valuation bound coefficientwise. The first two zero moments
then contradict (5), since the forcing at \(n+1\) has valuation
\((n+1)g<jg\). If \(t<1\), multiply the difference by \(1+T\).
The two parameters in \(1+\gamma T+e^2T^2\) have minimum valuation
\(g\): their sum has valuation \(g<t\) and their product has
valuation \(2t\). All remaining parameters therefore have valuation
at least \(g\). The displayed numerator has zero constant term,
valuation at least \(1+t>g\) at degree one, at least
\(\min(1+t+g,2t)>g\) at degree two, and at least \(2t>2g\)
at degree three. Hence the filtered rational difference at \(n+1\)
has valuation strictly greater than \(ng\), as does the filtered
matching error. But \(P_{n+1}^0+P_n^0\) has valuation exactly \(ng\)
by (5), whereas \(P_{n+1}+P_n=0\).

It remains to treat \(g\ge t\), allowing \(\gamma=0\).
The \(B_i\) recurrence gives \(v(B_i)\ge it\) and then
\(v(z_i)\ge it\). Since \(jt>t\), the first force (3) gives
\(v(1-\beta)=t\). The subtraction of the two forces therefore
requires
\[
     v(z_{j+1}-z_j)=1+t,\qquad (j-1)t\le1.
\]
Remove \(1+T\) from the denominator of (4). The two \(y\) parameters
have valuation \(t\), by their sum and product, and the numerator
has the form \(1+cT\) with \(c\) a unit. If \((j-1)t<1\),
the two coefficients at \(n+1,n+2\) contradict (4) as before.
At equality the same argument applies to the coefficient \(n+1\)
whose scale \(nt\) is below one. At \(n+2\) the scale is one.
Compare the error \(P-R\) to the error at
\((\beta,\gamma,d,e)=(1,0,0,0)\), whose rational series in
(4) is \(1\) and whose moment series has \(i\)-th coefficient
\(2^i i!\) by (5). All parameters have changed by valuation at least
\(t>0\), so each error difference has valuation at least \(1+t\).
At indices \(n+2\) and \(n+1\), both at least two, the anchor error
already has valuation strictly above one. Thus the coefficient at
\(n+2\) of \((1+T)(P-R)\) has valuation strictly above one.
The exact zero moments now force zero residues of the complete
homogeneous coefficients at both degrees \(n,n+1\), which is
impossible for the two unit parameters. This completes the proof.
