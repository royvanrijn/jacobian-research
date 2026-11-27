# Binary reduced lower pencils in profiles (6,5,4) and (6,5,5)

**Theorem KellerFiniteBirationalBinaryLeadingExclusion.**
Let a characteristic-zero three-variable Keller map have ordinary
component degrees \(d_1,d_2,d_3\), with \(\min_{i\ne j}(d_i+d_j)>2\).
Suppose its highest forms \(A,B,C\) are binary in the same linear
variables \(x,y\). If \(k[x,y]\) is finite over \(k[A,B,C]\) and
the two rings have the same fraction field, no such map exists.

**Theorem SexticKeller654ReducedLowerPairInversion.**
A characteristic-zero three-variable Keller map of ordinary component
degrees \((6,5,4)\) is invertible if the quintic and quartic highest
forms are coprime and squarefree over the algebraic closure.

**Theorem SexticKeller655ReducedLowerPairExclusion.**
There is no characteristic-zero three-variable Keller map of ordinary
component degrees \((6,5,5)\) whose two quintic highest forms are
coprime over the algebraic closure.

All arguments below are over an algebraically closed field \(k\) of
characteristic zero. Invertibility descends: an inverse over an
extension is unique, and its coefficients can equivalently be found
over the original field by the linear equations for a left inverse of
the same degree.

## Lifting the finite normalization at the first new weight

The argument also works with positive rational weights. Let \(H_i\)
be the leading forms of \(F_i\), of weights \(d_i\), independent of
a source variable \(z\). Increase the weight of \(z\) to the first
value \(w\) at which a monomial containing \(z\) ties with an \(H_i\).
Such a monomial exists, since a Keller map depends on all three
variables. Let \(L_i\) be the new leading forms. Then \(L_i|_{z=0}=H_i\).
If \(N=\sum_i d_i-\sum_j w_j>0\) at the new weights, the highest
determinant vanishes. If the binary image has dimension two, the
image closure of \(L\) has dimension two as well, contains the
binary image, and is irreducible. The image closures therefore agree.
Evaluation identifies \(k[H_1,H_2,H_3]\) with \(k[L_1,L_2,L_3]\).

Every element of the binary image field integral over its image ring
thus lifts to a rational function integral over \(k[x,y,z]\), hence
to a polynomial. Its restriction is the original element, since a
nonzero denominator from the image ring restricts to a nonzero
polynomial. Homogeneous elements lift at their own weights by
equivariance under weighted scaling, after clearing the weight
denominators. In the finite birational ordinary case the two
variables lift at weight one, so the lifted polynomials are \(x,y\)
themselves. This contradicts the defining tie. Indeed \(w\le d_i\)
for an index where the tie occurs, so the bound on the other two
degrees guarantees \(N>0\).

In the applications two coprime binary highest forms already make
\(k[x,y]\) finite over their algebra: their ideal has finite-dimensional
quotient, and homogeneous induction lifts a finite basis to module
generators. The same argument applies with positive weights and can
be reapplied after a triangular change that preserves the current
weights. For that purpose finiteness of the normalization generators
over the image ring is also transported by the isomorphism above.

We use the [leading pencil saturation proof](SEXTIC_KELLER_MIXED_DEGREE_LEADING_REDUCTIONS.md#descending-through-the-lower-terms)
with positive rational weights. The graded basis argument there uses
the regular sequence of the two normalization generators and works
after clearing weight denominators. The determinant step at loss
\(r\) needs \(N>r\). We will only apply it with \(N>\max_i d_i\),
so this holds throughout the descent.

## Binary (6,5,4)

The [nonbinary reduced pair exclusion](SEXTIC_KELLER_654_NONBINARY_EXCLUSION.md)
gives the binary case under the theorem's hypotheses. Write
\(t=y/x\), \(R=B^4/C^5\), \(u=B/C\), and \(V=A/u^6\). Then
\(k(A,B,C)=k(u,R,V)\), and its degree-zero field is \(k(R,V)\),
since \(u\) is an independent variable of degree one over \(k(t)\).
By Lüroth's theorem this degree-zero field is \(k(r)\) for a
rational binary pencil of degree \(m=[k(t):k(r)]\). The outer
degree of \(R\) in \(r\) is \(s=20/m\), by coprimality.
The [fibre valuation table](SEXTIC_KELLER_MIXED_DEGREE_LEADING_REDUCTIONS.md#the-profile-654)
applies verbatim to the rational function \(R\) on the binary
projective line: its member multiplicities use only coprimality
and squarefreeness of \(B,C\). In particular \(m=20,10,5\) are
impossible. For \(m=20\), the field is \(k(B,C)\) and surjectivity
of \((B,C)\) contradicts \(A\ne0\) in degree six. For \(m=10,5\)
the stated pencil contradictions apply. If \(m=1\), adjoining
\(u\) gives the entire binary field and the first theorem applies.

For \(m=2\), after linear changes the table gives
\[
 B=x(K-a_1x^2)(K-a_2x^2),\qquad
 C=(K-b_1x^2)(K-b_2x^2)
\]
with second generator \(K\) quadratic and \(K(0,y)\) a nonzero
multiple of \(y^2\). The degree-zero parameter is \(K/x^2\),
and adjoining \(u\) recovers \(x\), so the image field is
\(k(x,K)\). Its two generators are integral over the image
ring and \((x,K)\) is onto; in particular \(A\in k[x,K]\).
At the first increase of the weight of \(z\), the lift of
\(x\) is \(x\); a nontrivial lift of \(K\) requires \(w=2\)
and is \(K+\gamma z\) with \(\gamma\ne0\).
Use this as a new variable \(Z\) of weight two. Increase the
weight of \(y\), holding \(x,Z\) at weights one and two.
At the first tie the lift of \(Z\) is again nontrivial and
can only be \(Z+\delta y\), with \(y\) of weight two and
\(\delta\ne0\). Use that as the new second variable and
increase the weight of \(y\) again. At a strictly larger
weight the lifted generators of weights one and two cannot
involve \(y\). This contradicts the existence of the next tie.
In all these uses the two fixed source weights sum to at
most three and the weight of the variable being increased
is at most six, so \(N\ge6\) at the tie.

For \(m=4\) the table gives \(B=x(C+\lambda x^4)\) with
\(\lambda\ne0\), after scaling. The degree-zero parameter is
\(C/x^4\), and adjoining \(u\) gives \(k(x,C)\). The pair
\((x,C)\) is onto since \(C(0,y)\) is a nonzero multiple of
\(y^4\), so \(A=a x^6+b x^2C\). The two normalization
generators are integral over the image ring. At the first
increase of the weight of \(z\), only \(C\) can acquire a
nontrivial lift \(H\). The leading tuple is
\[
 (a x^6+b x^2 H,\ xH+\lambda x^5,\ H).
 \tag{1}
\]
The first tie must occur at one of \(4/3,3/2,2,3,4\),
by the exponent of \(z\) in \(H\). If the lift is linear
in \(z\) and its weight is at least two, the actual quartic is
affine in \(z\), since nonlinear terms would otherwise either
tie or exceed weight four. The repository's
[sextic quartic affine-source inversion](SEXTIC_KELLER_SQUAREFREE_CUBIC_QUADRATIC_MATES.md)
then applies in the original linear frame.

For the other ties, the
[weighted rational slice](KELLER_LEADING_PAIR_RATIONAL_SLICE.md)
for the last two entries of (1) forces the relative area
differential of \(H\) on a generic plane \(x=\text{constant}\)
to be exact on the geometric generic level. In fact their
Jacobian derivation is \((H+5\lambda x^4)\det J(x,H,\cdot)\);
the factor is nonzero and constant on the generic fibre.
Over the algebraic closure of \(k(x)\) scale the plane
coordinates using their weights, so the calculation can be
done with \(x=1\) in \(H\). This uses only a nonzero constant
multiple of the plane area.

At weight \(4/3\) the plane polynomial is \(f(Y)+cZ^3\), with
\(\deg f=4\) and \(c\ne0\). The generic level has nonzero
regular differential \(dY/Z^2\) on its smooth complete model.
At the simple finite roots of \(f(Y)-S\) ramification of
index three gives regularity, and at infinity \(Y,Z\) have
pole orders three and four. At weight \(3/2\) the plane
polynomial is \(f(Y)+a_1(Y)Z^2\) with \(\deg f=4\) and
\(0\le\deg a_1\le1\). On the generic level the substitution
\(M=a_1(Y)Z\) gives \(M^2=a_1(Y)(S-f(Y))\) and the
differential is proportional to \(dY/M\). Its defining
polynomial is squarefree of degree four or five, so this
differential is regular and nonzero on the complete model.
A nonzero regular differential on a smooth complete curve
cannot have a rational primitive. A simple pole also
precludes a rational primitive.

It remains to consider weight two with a nonzero square
term in \(z\). Complete the square by a triangular weighted
change to a variable \(Z\) of weight two. The plane polynomial
becomes \(cZ^2+f_*(Y)\). If \(\deg f_*\ge3\), the same
squarefree hyperelliptic calculation gives a regular form;
if \(\deg f_*=2\), it gives simple poles at infinity.
Thus the lifted \(H\) has form
\[
 H=cZ^2+\gamma x^3 y+\delta x^4,\qquad c\ne0 .
\]
When \(\gamma\ne0\) this is a terminal case treated below.
If \(\gamma=0\), repeat the lifting argument with \(x,Z\)
fixed of weights one and two and increase the weight of the
unused variable \(y\) from one to its first tie. The lifted
quartic still restricts to \(cZ^2+\delta x^4\) on \(y=0\).
Its possible tie weights are \(4/3,3/2,2,3,4\), by the
exponent of \(y\) and the remaining weight contributed by
\(x,Z\). At \(4/3,3/2\) the plane polynomial adds a
nonzero cubic or quadratic in \(Y\), respectively, to
\(cZ^2+\delta\); the hyperelliptic check gives the same
contradiction. At weight two any cross term \(ZY\) is
absorbed by completing the square in \(Z\). The result has
degree at most one in \(Y\) by the hyperelliptic check.
If it is constant, the leading tuple is again independent
of \(y\) in the shifted coordinates, and its next tie must
be strictly above two. At weight \(j=3,4\) the \(y\) term
is necessarily \(c_1x^{4-j}y\), since the leftover weight
is too small for a factor \(Z\). All the cases left are
therefore terminal cases with
\[
 H=cZ^2+\gamma x^{4-j}y+\delta x^4,\qquad
 c\gamma\ne0,\qquad
 (w_y,w_Z)\in\{(1,2),(2,2),(3,2),(4,2)\},\quad j=w_y .
 \tag{2}
\]
The determinants at the first ties have positive expected
degree because the fixed weights sum to at most three and
the moving weight is at most six. Each triangular change
used above has constant Jacobian and preserves the current
weights. In the binary phases \(x,H\) remain integral over
the corresponding image by the transport argument, including
after the changes.

### Descent in the terminal quartic cases

For (2), the pair \((x,H)\) has relatively closed field,
since over \(x\ne0\) it reduces to a coordinate pair by
solving for \(y\), leaving \(Z\) free. It is onto since off
\(x=0\) one can solve for \(y\), and on \(x=0\) the term
\(c Z^2\) still ranges freely (if \(j=4\), the \(y\) term
also remains available). Its polynomial kernel is \(k[x,H]\)
by contraction. Also \(x,H\) are coprime, and
\(N=12-j\ge8\). The derivative minors of (1) include
\(H+5\lambda x^4\) and \(2x(bH+3ax^4)\). They are coprime
unless \(b\ne0\) and \(3a=5\lambda b\). The weighted
saturation proof applies when they are coprime.

In the exceptional case scale the sextic so that \(b=1\).
Use reference polynomials \(G_A(x,F_4),G_B(x,F_4)\), beginning
with (1) evaluated at \((x,F_4)\). At first error loss
\(r\), write \(E_A,E_B\) for the homogeneous errors of
weights \(6-r,5-r\). The linear determinant gives
\(E_A-2xE_B=D(x,H)\), since the derivative columns of
(1) are \((H+5\lambda x^4)(2x,1,0)\) and \((x^2,x,1)\).
If \(j\le2\) then \(N\ge10\), and \(2r<N\) for \(r\le4\).
If \(j\ge3\), then \(2r<N\) for \(r\le3\), and for
\(r\ge4\) the error \(E_B\) has weight at most one,
so is already in \(k[x]\). In all cases a constant or
absent \(E_B\) can be absorbed directly together with \(D\).

For the remaining losses use the
[critical-level descent](SEXTIC_KELLER_MIXED_DEGREE_LEADING_REDUCTIONS.md#descent-on-the-nonbinary-discriminant-level).
Here is the calculation with the weights in (2). Rescale by
the source weights, so the determinant is a nonzero constant
times \(\epsilon^N\). Fix the value of \(F_{4,\epsilon}\) to
\(S=-5\lambda\) and use \(x,Z\) on that level, solving for
\(y\) near \(x=1\) by a formal implicit equation with unit
derivative. Write the first two rescaled components there as
\(a(x,\epsilon)+\epsilon^r U\) and \(b(x,\epsilon)+\epsilon^r V\)
in sextic, quintic order, with \(a,b\) the references.
At order zero in \(\epsilon\), \(U=2xV+D(x,S)\).
Now \(b_x(1,0)=0\), \(b_{xx}(1,0)=20\lambda\), and the
formal solution of \(b_x(\xi(\epsilon),\epsilon)=0\)
starts at \(\xi(0)=1\). At this center the two-coordinate
determinant is
\[
 \epsilon^r a_x V_Z +\epsilon^{2r}(U_x V_Z-U_Z V_x)
\]
up to sign and is zero modulo \(\epsilon^N\). Put
\(g(Z)=V(1,Z,0)\). If \(a_x\) at the center has order
\(\ell<r\), order \(r+\ell\) gives \(g'=0\).
Otherwise order \(2r\) gives
\((\alpha+2g+D_x(1,S))g'=0\), with \(\alpha\) the
order-\(r\) coefficient of \(a_x\) at the center. In either
case \(g\) is constant.

The entire plane level \(H(1,y,Z)=S\) is irreducible
because its \(y\)-coefficient is nonzero constant. Hence
\(E_B(1,y,Z)-g\) is divisible by \(H(1,y,Z)-S\).
With \((y,Z)\)-weights \((j,2)\), the divisor has weight
four and the dividend has weight at most \(5-r\le4\).
Thus the quotient is constant or zero. Homogenizing in
\(x\) puts \(E_B\) in \(k[x,H]\) and the relation with
\(D\) puts \(E_A\) there as well. Add these corrections to
the references, using the actual \(F_4\) for \(H\). The
loss increases at each step, until the tuple factors through
\(x,F_4\), contrary to the nonzero determinant.

## Binary (6,5,5)

First choose two coprime squarefree members of the quintic
pencil by an invertible linear change on the lower components.
Indeed, for coprime independent forms \(B,C\) there is a source
partial derivative \(\partial\) for which \(C\partial B-B\partial C\)
is nonzero. Every repeated factor of \(B-tC\) divides that
fixed polynomial, and distinct members have no shared factor
by coprimality. Thus all but finitely many members are squarefree,
and any two distinct such members are coprime.

The [aligned powers exclusion](SEXTIC_KELLER_ALIGNED_LINEAR_POWERS.md)
already makes the leading tuple binary. Write \(R=B/C\),
\(u=A/B\), and \(\eta=B^6/A^5\). The full image field is
\(k(u,R,\eta)\) with degree-zero field \(k(R,\eta)\); here
\(B=u^5\eta\). Its index in \(k(t)\) divides five because
\(R\) has degree five. Index one is the birational exclusion
above. For index five the degree-zero field is \(k(R)\).
The [quintic pencil valuation argument](SEXTIC_KELLER_MIXED_DEGREE_LEADING_REDUCTIONS.md#the-equal-degree-lower-pencils)
then gives a pure fifth-power member: the orders of \(\eta\)
at zero and infinity are \(1,0\) modulo five, respectively,
and another order nonzero modulo five forces its member to
be a pure fifth power. Write the lower pencil as
\(\langle x^5,H\rangle\). The degree-zero field is
\(k(H/x^5)\), and \(A/B\) gives the radial variable \(x\)
over that field. The image field is \(k(x,H)\), both
generators are integral over the image, and \(A=a x^6+b xH\)
by contraction for the onto pair \((x,H)\).

At the first increase of the weight of \(z\), only the
quintic generator can acquire a nontrivial lift \(H'\).
The weights are \((1,1,w)\), \(1<w\le5\). The pair
\((x,H')\) is onto by its restriction to \(z=0\). Its
field is relatively closed: over \(\overline{k(x)}\),
\(H'\) as a plane polynomial has degree five and highest
ordinary plane form a multiple of \(y^5\), because every
\(z\)-monomial has ordinary degree less than five. It depends
on \(z\). If its field were not relatively closed, its
closed generator could be chosen polynomial, and \(H'\)
would be a polynomial of degree five in a linear form.
The leading \(y^5\) would force this form to be independent
of \(z\), a contradiction. For the polynomial generator fact, use
[Arzhantsev and Petravchuk, Theorem 1 and Proposition 1](https://arxiv.org/pdf/math/0608157v2)
and [Petravchuk and Iena, Remark 2](https://arxiv.org/pdf/math/0701588v2).
The relative closedness descends from \(\overline{k(x)}\),
since its rational functions in \(H'\) that are fixed
over \(k(x)\) belong to \(k(x)(H')\).
If \(b\ne0\), the derivative minors include \(5x^4\) and
\(6ax^5+bH'\), which are coprime as \(H'(0,y,0)\) is
a nonzero multiple of \(y^5\). Saturation applies since
\(N=14-w\ge9\). Thus \(b=0\), and in the new weights
the tuple has leading forms \(x^5,x^6,H'\) in that order.

### Unequal transverse weights in the aligned fifth and sixth powers

Denote the components by \(P,Q,T\), in weights \((1,1,w)\).
Use the rescaling in an indeterminate \(\epsilon\) with
rational exponents, or clear all original weight denominators.
Its determinant is \(\kappa\epsilon^N\) with \(N=14-w\ge9\).
Set \(X=P_\epsilon^{1/5}=x+o(1)\), with coefficients in
\(k[x,x^{-1},y,z]\). Then
\[
 Q_\epsilon\equiv \sum_{0\le i<N}\epsilon^i D_i(X,T_\epsilon)
       \pmod{\epsilon^N},\qquad D_0=X^6 .
 \tag{3}
\]
Indeed at each unmatched order \(i<N\), including fractional
orders, the bracket of the residual coefficient with \(x,H'\)
vanishes by the Keller determinant and the previous
subtractions. The relatively closed pair field and polynomial
contraction give that coefficient in \(k[x,x^{-1},H']\).
It is homogeneous of weight \(6-i\). In this ring the two
generators have integral weights one and five, so any nonzero
\(D_i\) has \(i\) integral, and it is a combination of
\(X^{6-i-5j}T_\epsilon^j\) for \(j\ge0\).
This gives (3) by induction. As in the
[aligned powers calculation](SEXTIC_KELLER_ALIGNED_LINEAR_POWERS.md#formal-expansion-at-the-first-weight),
the formal root subtraction resembles
[Guccione, Guccione and Valqui, Proposition 1.7](https://arxiv.org/pdf/1406.0886v3).
Subtract a constant times \(T\) from \(Q\) so that the
coefficient of \(T_\epsilon\) in \(D_1\) in (3) vanishes;
this subtraction has exactly that effect on (3).

Decrease the weight of \(x\) from one to its first tie \(v\)
in either \(P\) or \(Q\). Such a tie occurs with \(0<v<1\)
since \(P,Q\) cannot both be univariate in \(x\). The weight
is at least \(1/6\) since transverse monomials have weight at
least one. Now \(P_*,Q_*\) are monic in \(x\) of degrees
five and six, and \(T_*=h_0=H'(0,y,z)\) has weight five
and a nonzero \(y^5\) term. The determinant has positive
expected degree \(10v+4-w>0\). Thus \(P_*,Q_*\) are
algebraically dependent over \(k(h_0)\).

The [monic parametrization observation](SEXTIC_KELLER_ALIGNED_LINEAR_POWERS.md#the-first-weight-at-which-the-pure-powers-change)
shows that the coefficients of \(P_*,Q_*\) lie in the relative
algebraic closure of \(k(h_0)\) in \(k(y,z)\). It applies
with these weights because the common monic polynomial
parameter has degree dividing \(5,6\) in \(x\), hence
degree one and weight \(v<1\); no transverse monomial has
that weight, so the parameter is \(x\).
For any homogeneous nonzero coefficient \(a\) of weight
\(s>0\) in that closure, Euler's identity and
\(da\wedge dh_0=0\) give \(a^b/h_0^\alpha\in k^\times\)
for \(\alpha/b=s/5\) in lowest terms. The gcd of the
factor multiplicities of \(h_0\) divides five, by restricting
a putative power to \(z=0\). If \(h_0\) is not a fifth
power, that gcd is one. So every positive residual
weight must be a multiple of five. Since \(5v<5\) and
\(6v<6\), the only possible new term is then a nonzero
multiple of \(h_0\) in \(Q_*\), with \(v=5/6\) and
\(P_*=x^5\).

For comparison in (3) use minimum weights \(1-v,1,0,0\)
on \(x,\epsilon,y,z\), respectively. The minimum parts
of \(P_\epsilon,Q_\epsilon,T_\epsilon\) encode \(P_*,Q_*,h_0\)
and have weights \(5(1-v),6(1-v),0\). A monomial
\(\epsilon^i X^{6-i-5j} T_\epsilon^j\) in (3) has
minimum weight \(i+(1-v)(6-i-5j)\), which is at least
\(6(1-v)\) for every nonzero coefficient. For this last
claim induct on the order in \(\epsilon\). The leading
order of each new correction has distinct weights as \(j\)
varies, and the previous corrections have no terms of
smaller weight. At equality the indices satisfy
\(iv=5(1-v)j\). In the nonpower case \(v=5/6\) and the
formal root has pure minimum part \(x\). Comparing the
minimum-weight coefficient at order one in (3) shows
that the nonzero \(h_0\) term would require a nonzero
\(T_\epsilon\) coefficient in \(D_1\). That coefficient
was removed.

If \(h_0\) is a fifth power, its root has weight one,
so it is a multiple of \(y\). Normalize to \(h_0=y^5\).
The same monic parameter argument then gives
\(P_*,Q_*\in k[x,y]\). Write \(v=p/q\) in lowest terms.
A tie in these binary monic polynomials implies \(q\le6\).
Their minimum parts after rescaling have the form
\[
 x^5 f(\tau),\qquad x^6 g(\tau),\qquad
 \tau=\epsilon^{q-p}y^p/x^q,\qquad f(0)=g(0)=1,
 \quad \deg f\le\lfloor5/q\rfloor,\quad
 \deg g\le\lfloor6/q\rfloor .
\]
At least one of \(f,g\) is nonconstant. The equality
\(ip=5(q-p)j\) in (3) gives the binomial identity in the
[ordinary aligned powers proof](SEXTIC_KELLER_ALIGNED_LINEAR_POWERS.md#formal-expansion-at-the-first-weight),
with \(e=5\), exponent \(N=14-w\), and coefficient at
\(\tau^n\) usable for \(n(q-p)<N\). For \(p=5\) one
has \(q=6\), so \(f=1\) and the order-one comparison
gives \(g=1\) because \(D_1\) has no \(T_\epsilon\) term.
Otherwise \(5\) divides \(n\) for any nonzero correction.
For \(q=6\) the coefficient at \(\tau\) suffices, of
\(\epsilon\)-order at most five. For \(q=4,5\) compare
through \(\tau^2\) (order at most eight); for \(q=3\)
compare through \(\tau^3\) (order at most six). For \(q=2\)
compare through \(\tau^6\) (order six), including the
correction at \(\tau^5\). The binomial coefficient
comparisons from that proof give \(f=g=1\) in each case.
All orders used are strictly below \(N\ge9\), giving
the required contradiction.
