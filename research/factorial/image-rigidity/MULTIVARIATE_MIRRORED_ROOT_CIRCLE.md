# Simple unit-circle roots for mirrored binomial moments

For mirrored exponent profiles with equal Gamma shapes within each
paired row, every fixed-order moment polynomial has simple unit-circle
roots. Private exponents are unrestricted. The consecutive pairs
at orders four and five and at orders five and six strictly interlace
for these profiles. A single mirrored unit gap gives strict
interlacing at all orders by a comparison of two Gegenbauer Jacobi
matrices. Applying
Schur–Szegő
composition separately at each degree does not in general preserve
coprimality of two successive polynomials, even when the input pairs
strictly interlace and have positive log-convex binomial-normalized
coefficients.

## Mirrored profiles at every fixed order

**Theorem.** Fix integer profiles \(b_i\ge0\) and \(d_i\ge1\),
\(1\le i\le r\), \(r\ge1\), and positive real Gamma shapes
\(\gamma_i\). Give the two independent variables \(X_i,Y_i\) common
shape \(\gamma_i\) and rate one. For
\[
 U=\prod_{i=1}^r X_i^{b_i+d_i}Y_i^{b_i},\qquad
 V=\prod_{i=1}^r X_i^{b_i}Y_i^{b_i+d_i},
\]
the polynomial \(Q_m(t)=\mathbb E[(U+tV)^m]\) has \(m\) distinct
zeros on \(|t|=1\) for every \(m\ge1\). Any additional independent
Gamma coordinates with equal exponents in both monomials contribute
only a nonzero scalar to \(Q_m\).

The coefficient of \(t^k\), divided by \(\binom mk\), is
\[
 \prod_i (\gamma_i)_{m b_i+d_i k}
          (\gamma_i)_{m b_i+d_i(m-k)}.
\]
Grouping the factors of a rising factorial by their residues modulo
\(d_i\) gives
\[
 (A)_{d_i k}
   =d_i^{d_i k}\prod_{j=0}^{d_i-1}
              \bigl((A+j)/d_i\bigr)_k .
\]
Here \(A=\gamma_i+m b_i>0\). After removing factors independent of
\(k\), the normalized coefficient array of \(Q_m\) is therefore the
entrywise product of the arrays
\((\lambda)_k(\lambda)_{m-k}\), one for each
\(\lambda=(\gamma_i+m b_i+j)/d_i>0\).

For any \(\lambda>0\) use the generating series
\[
 (1-2xv+v^2)^{-\lambda}=\sum_{m\ge0} C_m^\lambda(x)v^m .
\]
The corresponding binomial polynomial
\[
 P_{m,\lambda}(t)
   =\sum_{k=0}^m\binom mk(\lambda)_k(\lambda)_{m-k}t^k
\]
satisfies
\(P_{m,\lambda}(e^{2i\theta})
   =m! e^{im\theta} C_m^\lambda(\cos\theta)\).
This follows by factoring the generating function into the two
factors with arguments \(e^{i\theta}v\) and \(e^{-i\theta}v\).
For completeness, the standard orthogonality argument here works
for all positive \(\lambda\), including \(\lambda<1\). Apply
\((1-x^2)\partial_x^2-(2\lambda+1)x\partial_x+
(v\partial_v)^2+2\lambda v\partial_v\) to the generating function
to get zero. The coefficient \(y=C_m^\lambda\) then satisfies
\[
 \bigl((1-x^2)^{\lambda+1/2} y'\bigr)'
 +m(m+2\lambda)(1-x^2)^{\lambda-1/2}y=0.
\]
The eigenvalues \(m(m+2\lambda)\) are distinct and the boundary
terms vanish at \(\pm1\). Integration by parts gives orthogonality
for the positive integrable weight \((1-x^2)^{\lambda-1/2}\).
Since the polynomials have exact degrees \(m\), orthogonality to
every lower degree forces \(m\) sign changes in \((-1,1)\). Indeed,
multiplying by the linear factors at all sign-changing zeros would
otherwise give a nonzero integrand of one sign with a multiplier
of degree less than \(m\). Hence every \(P_{m,\lambda}\) has simple
unit-circle roots.

Use the degree-\(m\) Schur–Szegő composition
\[
 \Big(\sum_k\binom mk a_k t^k\Big)\star_m
 \Big(\sum_k\binom mk b_k t^k\Big)
 =\sum_k\binom mk a_k b_k t^k .
\]
Here is its disk argument with the hypotheses made explicit.
Suppose the first factor \(F\) has no zero in \(|t|<1\) and the
second is monic with unit-circle zeros \(\zeta_1,\ldots,\zeta_m\).
Then their composition is
\[
 (-1)^m\Big(\prod_{i=1}^m\zeta_i\Big)
 \widetilde F(-t/\zeta_1,\ldots,-t/\zeta_m),\qquad
 \widetilde F(x_1,\ldots,x_m)=\sum_{k=0}^m a_k e_k(x_1,\ldots,x_m),
\]
where \(e_k\) denotes the elementary symmetric polynomial.
This is the symmetric multi-affine polarization of \(F\). For
\(|t|<1\) its displayed arguments all belong to the open unit disk.
The Grace–Walsh–Szegő coincidence theorem for this convex circular
domain makes its value equal to \(F(\eta)\) at some \(|\eta|<1\).
Thus the composition has no interior root. For real palindromic
inputs of exact degree \(m\) the composition remains palindromic
of degree \(m\), so reciprocity excludes exterior roots as well.
The version of the coincidence theorem used here is Theorem 11,
pp. 478–479 in [Borcea and Brändén, Pólya–Schur master theorems for circular domains and their boundaries](https://annals.math.princeton.edu/wp-content/uploads/annals-v170-n1-p14-p.pdf).

Simplicity follows by openness. In the real palindromic space of
exact degree \(m=2l\), write \(F(t)=t^l H(t+t^{-1})\), where \(H\)
is real of degree \(l\). At \(m=2l+1\), first remove \(t+1\) and
use the same expression for its palindromic quotient. The
circle-rooted set corresponds exactly to \(H\) having all roots
in \([-2,2]\), and its relative interior corresponds exactly to
distinct roots in \((-2,2)\). At an endpoint one can move a root
outside the interval, and at a repeated root one can perturb two
equal real factors to a nonreal conjugate pair. In the other
direction simple interior roots persist as real interior roots
under small real coefficient changes. Distinct interior roots of
\(H\) are equivalent to simple unit-circle roots of \(F\), including
the simple structural factor \(t+1\) at odd degree. For any positive
palindromic circle-rooted \(P_{m,\lambda}\), composition by it is
an invertible linear map on the real palindromic space. The disk
argument shows it maps the circle-rooted set into itself. An
invertible linear map is open, so it also maps the interior into
itself. Composing all factors \(P_{m,\lambda}\) now proves the
location and simplicity claimed for \(Q_m\).

For example \(U=X^3Y,\ V=XY^3\) has no private exponents. The
theorem still gives simple circular roots at all orders. The
[private-exponent theorem](SFC_OVERLAPPING_SWAP_PAIRS.md)
proves the stronger statement of strict interlacing for its
specified parameter range.

## Consecutive interlacing for a single mirrored unit gap

**Theorem.** For the profiles in the first theorem with \(r=1\) and
\(d_1=1\), every consecutive pair strictly interlaces on the unit
circle. This holds for every common positive Gamma shape \(\gamma\)
and every integer \(b_1\ge0\). At ordinary shape it gives all-order
Strong Factorial for \(C(Z)(XY)^{b_1}(A X+B Y)\), where \(C\) is a
monomial in independent variables, over every characteristic-zero
field.

More generally, the zeros of \(C_n^\lambda\) and
\(C_{n+1}^{\lambda+h}\) strictly interlace if \(h\ge0\) and
\(\lambda\ge1+h(n-1)\). The case \(h=0\) holds for every \(\lambda>0\).
For the stated moment pair take \(\lambda=\gamma+n b_1\) and \(h=b_1\).
When \(b_1\ge1\) this satisfies the required bound because
\(\lambda-h(n-1)=\gamma+b_1>1\). The map \(x\mapsto
e^{2i\arccos x}\) carries strict alternation on \((-1,1)\) to
strict alternation on the circle. The binomial polynomials differ
from the mapped Gegenbauer polynomials by nonzero factors as in the
fixed-order proof.

Normalize each \(C_j^\lambda\) to be monic. Differentiating its
generating series with respect to the series variable gives
\[
 (j+1)C_{j+1}^\lambda(x)
    =2(j+\lambda)xC_j^\lambda(x)-(j+2\lambda-1)C_{j-1}^\lambda(x).
\]
The leading coefficient is \(2^j(\lambda)_j/j!\). Hence the monic
recurrence is \(p_{j+1}=xp_j-c_j(\lambda)^2p_{j-1}\), with
\[
 c_j(v)^2=\frac{j(j+2v-1)}{4(j+v)(j+v-1)}
          =\frac14\left(1-\frac{v(v-1)}{(j+v)(j+v-1)}\right).
 \tag{G1}
\]
Thus the characteristic polynomial of the size-\(n\) zero-diagonal
Jacobi matrix \(J_n(v)\) with off-diagonal entries \(c_1(v),\ldots,
c_{n-1}(v)\) is \(p_n\) at \(v=\lambda\). Deleting the first or last
row and column gives strict eigenvalue interlacing. Indeed equality
in principal-submatrix interlacing would give a shared eigenvalue,
and the determinant recurrence with positive off-diagonals
propagates its vanishing to the determinant of the empty matrix.
This also proves the case \(h=0\).

We use the following comparison for two zero-diagonal Jacobi matrices
of the same size with positive off-diagonals. If the successive ratios
of the second matrix's off-diagonal entries to the first matrix's
are \(\rho_1\ge\cdots\ge\rho_{n-1}\ge1\), each positive eigenvalue
of the second is at least the corresponding positive eigenvalue of
the first, in decreasing order. To see this, put \(u_n=0\) and
\(u_j=\log\rho_j-u_{j+1}\). Each \(u_j\) is a nonnegative alternating
sum of decreasing nonnegative numbers. The second matrix is the
first multiplied on both sides by \(\operatorname{diag}(e^{u_j})\).
After putting odd vertices first, the positive eigenvalues of a
zero-diagonal path matrix are the positive singular values of its
off-diagonal block. Each diagonal scaling multiplies the appropriate
rows or columns of that block by factors at least one, which can
only increase its singular values. For example, scaling columns
replaces \(B B^\top\) by \(B D^2 B^\top\ge B B^\top\), and scaling
rows uses \(B^\top B\) the same way.

Let \(L=\lambda+h\). The leading \(n\)-by-\(n\) principal block of
\(J_{n+1}(L)\) is \(J_n(L)\). For this block the comparison uses
the ratios \(c_j(\lambda)/c_j(L)\). For \(v\ge1\), the quantity
\[
 -\partial_v\log c_j(v)
  =\frac12\left(\frac1{j+v}+\frac1{j+v-1}\right)
    -\frac1{j+2v-1}
\]
is positive and decreases in \(j\). The same comparison of
denominators proves the derivative sign in \(j\), since the first
two reciprocal squares are at least the third and one is larger.
Integration for \(v\in[\lambda,L]\) proves that the ratios are
at least one and decrease with \(j\). The positive eigenvalues of
the leading block are therefore at most those of \(J_n(\lambda)\).

The trailing \(n\)-by-\(n\) principal block has entries \(c_{j+1}(L)\).
For \(1\le j\le n-1\) we have \(hj\le\lambda-1\), hence
\[
 \frac{L}{j+1}\le\frac{\lambda}{j},\qquad
 \frac{L-1}{j+1}\le\frac{\lambda-1}{j}.
\]
Since \(v(v-1)/[(j+v)(j+v-1)]\) is the product of the two
increasing ratios \(v/(j+v)\) and \((v-1)/(j+v-1)\), (G1) gives
\(c_{j+1}(L)\ge c_j(\lambda)\). The ratios decrease with \(j\).
In fact, for real \(j>0\),
\[
 2\partial_j\log c_j(v)=\frac1j
  \left(1+\frac1{1+x+y}-\frac1{1+x}-\frac1{1+y}\right),\qquad
 x=\frac vj,\quad y=\frac{v-1}{j}.
\]
The expression in parentheses vanishes at \(x=y=0\) and increases
in each nonnegative variable, by differentiating the two terms
involving that variable. Therefore replacing \((j,v)=(j,\lambda)\)
by \((j+1,L)\) decreases the displayed derivative for every real
\(j\) between \(1\) and \(n-1\). This proves the decrease of the
ratios and the opposite eigenvalue bound for the trailing block.

Write \(\mu_i\) for the eigenvalues of \(J_n(\lambda)\) and \(\nu_i\)
for those of \(J_{n+1}(L)\), both in decreasing order. For each
positive \(\mu_i\), the two block comparisons and strict
principal-submatrix interlacing give \(\nu_i>\mu_i>\nu_{i+1}\).
All matrices have symmetric spectra since changing the signs of
the even vertices conjugates a matrix to its negative. In odd
dimension the simple middle eigenvalue is zero. The positive
inequalities, symmetry, and this middle eigenvalue give strict
interlacing for every \(i\), including the case \(n=1\).
The ordinary factorial claim then follows from the rational
Bézout identity of the two polynomials. \(\square\)

## Consecutive interlacing at terminal five

**Theorem.** Every mirrored profile in the first theorem has coprime
order-four and order-five moment polynomials. The two polynomials
have strictly interlacing roots on the unit circle. In particular,
this proves terminal-five separation for ordinary mirrored
supports over every characteristic-zero field.

Normalize the constant coefficients to one. Write the binomial
coefficient arrays at orders four and five as
\[
 (1,u,v,u,1),\qquad (1,w,y,y,w,1).
\]
For each Gamma factor in the grouping above, put
\(T=(\gamma_i+4b_i+j)/d_i\) and \(T'=T+b_i/d_i\). Then
\(T\le T'\le5T/4\), and the four local entries are
\[
 u_i=\frac{T}{T+3},\qquad
 v_i=\frac{T(T+1)}{(T+2)(T+3)},\qquad
 w_i=\frac{T'}{T'+4},\qquad
 y_i=\frac{T'(T'+1)}{(T'+3)(T'+4)}.
\]
The entries of the original arrays are the products of these
local entries over the factors. In particular \(0<u<1\).

For \(0<a<1\) and \(0<t\le1\) define
\[
 F_a(t)=\frac{(1+a)t+1-a}{(1-a)t+1+a}.
\]
We have \(t^a\le F_a(t)\) and \(F_a(t)F_a(s)\le F_a(ts)\).
For the first inequality, the derivative of \(\log F_a(t)-a\log t\)
is nonpositive, since the product of the numerator and denominator
of \(F_a\) is \(4t+(1-a^2)(1-t)^2\). The difference vanishes at
\(t=1\). For the second, put \(\rho(t)=(1-t)/(1+t)\).
Then \(\rho(F_a(t))=a\rho(t)\) and
\(\rho(ts)=(\rho(t)+\rho(s))/(1+\rho(t)\rho(s))\).
Thus \(\rho(F_a(t)F_a(s))\ge\rho(F_a(ts))\), and \(\rho\) is decreasing.

The local identities \(v_i=u_i F_{1/3}(u_i)\) and
\(y_i=w_i F_{1/2}(w_i)\) now give
\[
 u^{4/3}\le v\le u F_{1/3}(u),
 \qquad
 w^{3/2}\le y\le w F_{1/2}(w).
\]
Also, writing \(r_i=3/T\) gives \(w_i=(1+c_i r_i)^{-1}\)
with \(16/15\le c_i\le4/3\) and \(u_i=(1+r_i)^{-1}\).
The inequalities \((1+r_i)^{4/3}\ge1+(4/3)r_i\) and
\[
 \prod_i(1+16r_i/15)\ge
 1+\frac{16}{15}\left(\prod_i(1+r_i)-1\right)
\]
give \(u^{4/3}\le w\le15u/(16-u)\). The product inequality
follows by expanding and comparing all nonconstant coefficients.
Finally \(t^a\ge t/F_{1-a}(t)\) by the first bound on \(F\).
This gives entirely rational envelopes
\[
 \begin{split}
 w_0=\frac{u^2(u+5)}{5u+1}\ \le &\ w\ \le
        w_1=\frac{15u}{16-u},\\
 v_0=\frac{u^2(u+5)}{5u+1}\ \le &\ v\ \le
        v_1=\frac{u(1+2u)}{2+u},\\
 y_0=\frac{w^2(w+3)}{3w+1}\ \le &\ y\ \le
        y_1=\frac{w(1+3w)}{3+w}.
 \end{split}                                                    \tag{1}
\]

Dividing the quartic by \(t^2\) and the quintic by \(t^2(t+1)\)
gives the respective monic quadratics in \(x=t+t^{-1}\)
\[
 A(x)=x^2+4u x+6v-2,\qquad
 B(x)=x^2+(5w-1)x+10y-5w-1 .
\]
Set \(D=1+4u-5w>0\) and \(C=1+10y-5w-6v\), so \(B-A=C-Dx\).
It suffices to show \(A(C/D)<0\), equivalently
\[
 P=C^2+4uCD+(6v-2)D^2<0.                                      \tag{2}
\]
For fixed \(u,w\), this expression is convex in \(v,y\), because the
quadratic part is the square of a linear function and the other
terms are affine. Thus it suffices to check the four corners
\((v_i,y_j)\) of (1). At each corner,
[the exact checker](verify_multivariate_mirrored_terminal_five.py)
clears the two positive envelope denominators by their squared
product. For the resulting polynomial in \(u,w\), it substitutes
\(w=(w_0+z w_1)/(1+z)\) and clears the positive denominator using
the degree in \(w\). It divides out \((1-u)^4\) exactly, then
substitutes \(u=r/(1+r)\) and clears the positive denominator
using the reduced degree in \(u\). At the four corners all
\(161,85,161,80\) coefficients, respectively, are strictly
negative. Thus the resulting polynomial is negative for
\(r>0,z\ge0\), with negative leading coefficient in \(z\) giving
the endpoint \(w=w_1\) as well. All factors dropped for this test
are positive for \(0<u<1\). This proves (2) on the whole domain
without sampling the row parameters.

If \(a_1<a_2\) are the roots of \(A\), (2) places \(C/D\) between
them. Hence \(B(a_1)>0>B(a_2)\), so the roots \(b_1<b_2\) of \(B\)
satisfy \(a_1<b_1<a_2<b_2\). Neither polynomial in \(t\) has zero
constant coefficient, and the additional quintic root \(-1\)
cannot be a root of the quartic, which is a positive even moment
there. The earlier unit-circle theorem places both quadratic
root sets in the interior of \([-2,2]\), so this ordering gives
strict interlacing on the circle. At ordinary shape the moment
polynomials have rational coefficients, and their rational
Bézout identity applies over any characteristic-zero field.

The replay uses SymPy 1.14.0 in the repository virtual environment.
~~~sh
.venv/bin/python research/factorial/image-rigidity/verify_multivariate_mirrored_terminal_five.py
~~~

## Consecutive interlacing at terminal six

**Theorem.** For every mirrored profile in the first theorem, the
order-five and order-six moment polynomials have strictly interlacing
arguments on the unit circle. In particular they are coprime, and
ordinary mirrored supports satisfy terminal-six separation over
every characteristic-zero field.

**Proof.** Normalize the constant and leading coefficients to one.
Write the two binomial-normalized coefficient arrays as
\[
 (1,u,v,v,u,1),\qquad (1,w,y,yh,y,w,1).
\]
For a grouped mirrored factor the shapes at the two orders are
\(T\) and \(T'\), with \(T\le T'\le6T/5\). Its local values satisfy
\[
 u_i=\frac{T}{T+4},\quad \frac{v_i}{u_i}=\frac{T+1}{T+3},
 \qquad
 w_i=\frac{T'}{T'+5},\quad
 \frac{y_i}{w_i}=\frac{T'+1}{T'+4},\quad
 h_i=\frac{T'+2}{T'+3}.
\]
The array is obtained by multiplying the local values. Using the
Möbius function \(F_a\) from the terminal-five proof, the three local
ratios \(v_i/u_i,y_i/w_i,h_i\) are the full factors
\(F_{1/2}(u_i),F_{3/5}(w_i),F_{1/5}(w_i)\), respectively.
Consequently the same product inequalities give
\[
 \begin{gathered}
 u^{3/2}\le v\le uF_{1/2}(u),\qquad
 w^{8/5}\le y\le wF_{3/5}(w),\qquad
 w^{1/5}\le h\le F_{1/5}(w),\\
 u^{5/4}\le w\le\frac{24u}{25-u},\qquad 0<u<1.                 \tag{3}
 \end{gathered}
\]
For the last pair of bounds, write \(r_i=4/T\). Then
\(u_i=(1+r_i)^{-1}\) and \(w_i=(1+c_i r_i)^{-1}\) with
\(25/24\le c_i\le5/4\). The lower bound follows from
\((1+r_i)^{5/4}\ge1+(5/4)r_i\). For the upper bound, expansion
by subsets gives
\[
 \prod_i(1+(25/24)r_i)\ge
 1+\frac{25}{24}\left(\prod_i(1+r_i)-1\right).
\]
The bounds in (3) allow the choices of \(y\) and \(h\) to vary
independently.

With \(x=t+t^{-1}\), the angle polynomials after division by
\(t^2(t+1)\) and \(t^3\), respectively, are
\[
 A(x)=x^2+(5u-1)x+10v-5u-1,\qquad
 B(x)=x^3+6w x^2+(15y-3)x+20yh-12w .
\]
Write \(a_1<a_2\) for the roots of \(A\). The earlier circle theorem
puts them in \((-2,2)\). The exact interval certificate below proves
\[
 B(a_1)>0>B(a_2)                                               \tag{4}
\]
throughout (3). The sextic gives \(B(-2)=-P_6(-1)<0\) because its
value at \(-1\) is a positive even moment after rescaling. It gives
\(B(2)=P_6(1)>0\) by its positive coefficients. Thus the three
roots of \(B\) lie one each in \((-2,a_1),(a_1,a_2),(a_2,2)\).
The structural quintic root \(t=-1\) separates the two conjugate
halves of the circle, giving the claimed interlacing. For ordinary
Gamma shapes use the rational Bézout identity.

### The interval certificate for the signs

For \(0\le u\le13/20\) the verifier uses the displayed formulas for
\(A,B\). At fixed \(u,w,v\), the value \(B(a_i)\) is separately
affine in \(y\) and \(h\), so the four corners of their ranges in
(3) suffice. The quadratic formula encloses both \(a_i\) for the
entire allowed interval of \(v\).

To retain the cancellations near \(u=1\), put \(p=1-u\),
\(q=1-w\), \(\rho=q/p\), and define
\[
 \begin{aligned}
  v&=1-\tfrac32p+Vp^2,\\
  y&=1-\tfrac85q+\tfrac{12}{25}q^2+\gamma q^3,\\
  h&=1-\tfrac15q-\tfrac2{25}q^2+\delta q^3 .
 \end{aligned}
\]
With \(Y=12/25+\gamma q\), substitution and cancellation give
\[
 \begin{aligned}
 \frac{A(pX-2)}{p^2}&=X^2-5X+10V,\\
 \frac{B(pX-2)}{p^3}
   &=X^3-6\rho X^2+15\rho^2 YX+\rho^3 Z,\\
 Z&=-10\gamma+20\left(y\delta-\frac{Y}{5}
                                -\frac{2}{25}(-\tfrac85+Yq)\right).
 \end{aligned}                                                   \tag{5}
\]
The upper endpoints of \(V,\gamma,\delta\) imposed by (3) are
\[
 \frac{3}{2(4-p)},\qquad \frac{12}{25(5-q)},\qquad
 -\frac{4}{25(5-2q)}.
\]
These identities follow by inserting \(u=1-p,w=1-q\) into the
Möbius expressions in (3). The lower endpoints are the respective
Taylor remainders for \((1-p)^{3/2},(1-q)^{8/5},(1-q)^{1/5}\)
after division by \(p^2,q^3,q^3\). The first two normalized
remainders increase with \(p,q\) and the last decreases with \(q\).
This follows from the integral remainder formulas, since the
relevant derivatives at \(1-t\) are
\((3/4)(1-t)^{-1/2}\), \(-(48/125)(1-t)^{-7/5}\), and
\((36/125)(1-t)^{-14/5}\), respectively. For
\(13/20\le u\le9/10\) these endpoints enclose the scaled formulas
in (5). As before, the four corners suffice because the scaled
cubic is separately affine in \(\gamma,\delta\).

For \(0\le p\le13/100\) use the uniform lower bounds
\[
 V\ge\frac38,\qquad \gamma\ge\frac8{125},\qquad
 \delta\ge-\frac{6}{125(1-q)^3}.
\]
The same derivatives and the Taylor remainder give these bounds.
In the last one the exponent \(3\) bounds \(14/5\) on
\(0<1-q\le1\). Equation (3) gives
\(25/(24+p)\le\rho\le5/4\). This part checks (5) directly over
its enclosing intervals in \(p,\rho\), including \(p=0\) in the
scaled expressions. It overlaps the preceding check and completes
the range \(u<1\).

The [verifier](verify_multivariate_mirrored_terminal_six.py) covers
the three compact ranges by bisection of rational rectangles.
On the first two ranges it discards rectangles disjoint from
the cross-order bounds on \(w\). It proves each required strict
sign using outward-rounded Arb intervals at 110-bit precision.
For the cubic evaluation at a root interval \(I\), it uses the
midpoint value plus \(B'(I)(I-\operatorname{mid}(I))\) on the
first two ranges. All bisection endpoints are rational and the
checks cover the entire rectangles. Replay with python-flint 0.9.0:

~~~sh
.venv/bin/python research/factorial/image-rigidity/verify_multivariate_mirrored_terminal_six.py
~~~

The three checks visited \(6971,835,43\) rectangles, respectively,
and certified all required signs. \(\square\)

## Obstruction to varying the composition degree

For \(j=1,2\), form monic pairs
\[
 f_j(t)=(1+t)(t^2+b_jt+1),\qquad
 g_j(t)=(t^2+a_j t+1)(t^2+c_j t+1)
\]
with
\[
 (a_1,b_1,c_1)=(-1/2,11/6,15/8),\qquad
 (a_2,b_2,c_2)=(-1/2,7/4,433/234).
\]
For each \(j\), we have \(-2<a_j<b_j<c_j<2\). The identity
\(t^2+a_jt+1=0\) on \(t=e^{i\theta}\) gives
\(a_j=-2\cos\theta\), which increases strictly on \(0<\theta<\pi\).
Together with the cubic root \(-1\) and conjugate symmetry this
proves strict interlacing in each input pair. All coefficients
are positive. In binomial normalization the quartics have
endpoint and interior arrays
\[
 (1,11/32,17/96,11/32,1),\qquad
 (1,79/234,503/2808,79/234,1).
\]
Each array is log-convex because its middle coefficient \(v\)
satisfies \(u^2\le v\le u\) for the adjacent coefficient \(u\).
The cubics likewise have normalized arrays \((1,u,u,1)\)
with \(0<u<1\).

Exact coefficient multiplication gives
\[
\begin{aligned}
 f_1\star_3 f_2 &= (1+t)(t^2+\tfrac{115}{72}t+1),\\
 g_1\star_4 g_2 &=
     (t^2+\tfrac{115}{72}t+1)(t^2-\tfrac{707}{624}t+1).
\end{aligned}
\]
For the cubic the interior coefficient is
\((17/6)(11/4)/3=187/72\). The quartic's coefficients at
\(t\) and \(t^2\) are \((11/8)(158/117)/4=869/1872\) and
\((17/16)(503/468)/6=8551/44928\), which agree with the
displayed factors. Both outputs remain circle-rooted, but they
share the nonreal pair from the common quadratic. This example
concerns arbitrary input polynomial pairs. The Gamma coefficient
relations across orders would require their own argument.
