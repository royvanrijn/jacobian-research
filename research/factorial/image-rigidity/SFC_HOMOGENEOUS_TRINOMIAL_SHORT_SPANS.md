# The signed remainder for initial homogeneous trinomial windows

**Theorem.** A nonzero homogeneous bivariate trinomial whose exponent
span in one row is at most 80 has a nonzero factorial moment among
orders 1, 2 and 3 over every characteristic-zero field. The exact
coefficient check includes every independent common translation.

Write the offsets as \((0,d),(a,d-a),(d,0)\), with integers \(0<a<d\),
and the common translation as \((L,H)\). At order \(m\) strip the
factor \((mL)!(mH)!\) from the factorial moment. The coefficient of
a count vector \(\alpha\) in the resulting form is
\[
 \binom m\alpha (mL+1)_{U}(mH+1)_{md-U},
 \qquad U=a\alpha_1+d\alpha_2 ,
\]
where \((z)_n\) denotes the rising factorial. Let \(p_0,p_1,p_2\) be
the first-order coefficients. Substituting
\((p_2X,p_2Y,-p_0X-p_1Y)\) spans the first zero plane since \(p_2>0\).
The second and third forms on that plane are
\[
 A X^2+BXY+CY^2,\qquad D X^3+E X^2Y+F XY^2+G Y^3 .
\]
The quadratic is positive definite over the reals because every
nonzero real coefficient direction gives a nonzero polynomial with
strictly positive square integral in the product-exponential measure.
If the forms share a complex projective zero, reality implies that
the cubic also vanishes at the conjugate root of the quadratic. The
quadratic then divides the cubic. Comparing the \(X^3,XY^2,Y^3\)
coefficients gives
\[
 S:=-R=C^2D+ABG-ACF=0.
\]
The exclusion over \(\mathbb C\) transfers to every characteristic-zero
field by the Nullstellensatz on the two rational projective charts.

## Common factorial factors and exact sign check

Put \(b=d-a\). In the coefficient of \(X^{m-k}Y^k\), let \(i\)
count the factors of \(X\) taken from column \(0\), and \(j\) count
the factors of \(Y\) taken from column \(1\), in place of column \(2\).
Put \(v=id+jb\). The unsigned product in this term is
\[
 (mL+1)_{md-v}(mH+1)_v\,
 (L+1)_d^{\,i+j}(L+1)_a^{\,k-j}
 (H+1)_d^{\,m-k-i}(H+1)_b^{\,k-j}.
 \tag{1}
\]
Its weight is
\(\binom{m}{k}\binom{m-k}{i}\binom{k}{j}(-1)^{m-i-j}\).
The greatest common product of linear factors in these summands,
up to a positive scalar, is
\[
 T_{m,k}=(L+1)_a^k(L+1)_d
    \prod_{\substack{1\le r\le ka\\ m\nmid r}}(L+r/m)
    (H+1)_{\max(b,\lfloor(m-k)d/m\rfloor)}.
 \tag{2}
\]

For the \(L\)-factors, first remove \((L+1)_a^k\). The remaining
first-order \(L\)-factors are \((L+1)_d^i(L+a+1)_b^j\).
For an integer shift \(r\le a\), if \(i=0\) then
\(md-v\ge md-kb\ge ma\). Thus the order-\(m\) product supplies
\(L+r\) when the first-order factors lack it. For \(a<r\le d\),
a first-order factor supplies it when \(i+j>0\). When \(i+j=0\),
the order-\(m\) product supplies it. The summand \(i=j=0\)
bounds the multiplicity after removal of \((L+1)_a^k\) by one
and excludes integer shifts above \(d\). The minimum \(md-v=ka\)
over all summands gives precisely the noninteger \(L\)-shifts.
For the \(H\)-factors, the summand \(i=j=0\) has only integer
shifts and the summand \(i=m-k,\ j=k\) has simple shifts.
A factor \(H+r\) with \(r\le b\) comes from an omitted column
if there is one, and otherwise from \(v\ge mb\). If \(r>b\),
the first-order supply requires an omitted \(X\) and \(r\le d\).
With all \(X\)'s kept, the minimum \(v\) is \((m-k)d\).
These conditions give (2).

The [factored checker](verify_sfc_homogeneous_trinomial_factored_sign.py)
uses primitive integer linear factors. It verifies (2) term by term
and also removes the common factors from the three signed products
in \(S\) before expanding. In every one of the \(3160\) oriented pairs
\(2\le d\le80,\ 1\le a<d\), the quotient has nonnegative coefficients
and a positive constant term. Hence \(S>0\) at all nonnegative real
heights in this span range. Checking both orientations matters for
this specific remainder because interchanging the variables changes
the endpoint eliminated by the chart. The
[unfactored checker](verify_sfc_homogeneous_trinomial_short_spans.py)
records the original 36 reflected profiles through span 12. The new
checker compares its fully expanded polynomial with that implementation
for every oriented pair through span 6.

## Highest homogeneous part at every span

For every \(d\ge2\) and \(0<a<d\) with integer offsets, \(S\) has
total degree \(14d-5\) in \(L,H\). Its highest homogeneous part is
\[
 \frac{2^{4d}3^{3d}}{72}\,a^2 d^3(d-a)^4\,
 L^{7d+4a-5}H^{7d-4a-5}(L+H)^4
 \bigl((d-2)L+(d+2)H\bigr).
 \tag{3}
\]
For each fixed support, this proves \(S>0\) once both heights exceed
an effectively computable threshold, with arbitrary height ratio.
Indeed, the summands in (1) bound the bidegree of \(S\) by
\((M_0,M_1)=(7d+4a,7d-4a)\). For \(d>2\), (3) has a positive
coefficient at every exponent pair of total \(M_0+M_1-5\) within
that rectangle. Each lower monomial therefore has both exponents
bounded by those of some top monomial. If \(L,H\ge T\ge1\),
the lower monomial is at most \(1/T\) times its dominating top monomial.
Let \(c_*\) be the smallest positive coefficient in (3), and let
\(B_-\) be the sum of the absolute values of the negative lower
coefficients. The total negative part is then bounded by
\(B_-/(Tc_*)\) times (3). Any \(T>B_-/c_*\) suffices.
The checked case \(d=2\) has \(S>0\) at all nonnegative heights.

To derive (3), put \(L=1/(xt), H=1/(yt)\) and divide the coefficient
of \(X^{m-k}Y^k\) by \(m^{md}L^{md+ka}H^{md-ka}\).
This same divisor applies to each summand of (1). In a single
normalized summand, \((mL+1)_n/(mL)^n\) becomes
\(\prod_{r=1}^n(1+rx t/m)\). Its logarithm has first three
coefficients \(\sigma_1(n)x/m,-\sigma_2(n)x^2/(2m^2),
\sigma_3(n)x^3/(3m^3)\), where
\[
 \sigma_1(n)=n(n+1)/2,\quad
 \sigma_2(n)=n(n+1)(2n+1)/6,\quad
 \sigma_3(n)=\sigma_1(n)^2.
\]
Use \(y\) in the other row and step \(1\) in the first-order factors
of (1). For logarithmic coefficients \(u,v,w\), the normalized
product has coefficients \(1,u,v+u^2/2,w+uv+u^3/6\).
Sum them with the signed weights following (1). Write \(A_r,\ldots,G_r\)
for the resulting normalized coefficients at order \(t^r\).
The common divisor for each product in \(S\) is
\(2^{4d}3^{3d}L^{7d+4a}H^{7d-4a}\).

Set \(s=x+y,\ \delta=x-y\) and \(K=a^2d^3b^4s^4/72\).
The quadratic coefficients start at order one, and the cubic
coefficients at order two. The first terms used in \(S\) are
\[
 (A_1,B_1,C_1)=\tfrac{s}{2}(d^2,2db,b^2),\qquad
 (D_2,F_2,G_2)=\tfrac{4\delta s}{9}(d^3,3db^2,b^3).
\]
The coefficient of \(t^4\) is zero by these formulas. At order five,
the terms from the order-two quadratic coefficients and the
order-three cubic coefficients, respectively, simplify to
\[
 \begin{aligned}
 2C_1C_2D_2+(A_1B_2+A_2B_1)G_2-(A_1C_2+A_2C_1)F_2&=-3K\delta,\\
 C_1^2D_3+A_1B_1G_3-A_1C_1F_3&=K(ds+5\delta).
 \end{aligned}
\]
These two identities follow by the finite sums above, also checked
symbolically in \(a,d,x,y\) by the factored checker. The total
\(K(ds+2\delta)\) gives (3) on restoring the divisor. The checker
also compares (3) with the actual leading homogeneous part at
each enumerated support.

## Coefficient signs in the seven highest degrees

Keep the normalization used for (3), and put \(a_0=a-1\),
\(b_0=d-a-1\). The coefficients of \(t^5,\ldots,t^{11}\) in the
normalized \(S\) have nonnegative coefficients in
\(\mathbb Q[a_0,b_0,x,y]\). Consequently the seven homogeneous parts
of \(S\) in degrees \(14d-5,\ldots,14d-11\) have nonnegative height
coefficients for every admissible support.

The [exact layer checker](verify_sfc_homogeneous_trinomial_top_layers.py)
uses the same signed counts as (1). For a normalized rising factorial,
the log coefficient at order \(r\) is
\((-1)^{r-1}\sigma_r(n)(x/m)^r/r\), where
\(\sigma_r(n)=\sum_{j=1}^n j^r\) is evaluated as a polynomial.
The power-sum formula follows by expanding
\(\sum_{j=1}^n e^{jz}=(e^{nz}-1)/(1-e^{-z})\), using
\(z/(1-e^{-z})\) for the Bernoulli numbers with \(B_1=1/2\).
If \(s_r\) are the resulting logarithmic coefficients of a summand,
the coefficients of its series satisfy
\[
 q_0=1,\qquad q_n=\frac1n\sum_{r=1}^n r s_r q_{n-r}.
\]
Differentiating the exponential of the log series gives this exact
recurrence. The checker expands the signed sums and then the three
products in \(S\) over \(\mathbb Q[a_0,b_0,x,y]\). It compares
the resulting layers with the expanded height polynomials for all
oriented supports through span 6.

At order \(t^{12}\), the coefficient of \(a_0^4 y^{12}\) is
\(-676639794713/5224277606400\). Thus the shifted-gap certificate
fails at the next layer. A fixed-support height coefficient is
obtained after evaluating the whole polynomial in \(a_0,b_0\).
The negative shifted-gap coefficient therefore leaves the sign of
that height coefficient to be decided by other terms.

## Newton expansion in the integer gaps

Write \(A_0=a-1,\ B_0=d-a-1\), and put \(M_0=7d+4a,\ M_1=7d-4a\).
Define the polynomial in reciprocal heights
\[
 Q_{A_0,B_0}(u,v)=
 \frac{u^{M_0}v^{M_1}}{432^d}\,S(1/u,1/v).
\]
For each fixed support this only reverses the height exponents and
divides the coefficients by a positive number. The Newton expansion is
\[
 Q_{A_0,B_0}=\sum_{i=0}^{A_0}\sum_{j=0}^{B_0}
 \binom{A_0}{i}\binom{B_0}{j}
 \bigl(\Delta_A^i\Delta_B^j Q\bigr)_{0,0},
\]
where the differences increase one gap by one. It follows by iterating
the univariate binomial inversion formula. Each difference on the
right is a polynomial obtained from supports of span at most \(i+j+2\).
Equivalently, in the normalized layer calculation, the coefficients
of \(\binom{A_0}{i}\binom{B_0}{j}\) are these same differences taken
one homogeneous part at a time.

The layer checker converts the powers of each gap using
\[
 z\binom zj=j\binom zj+(j+1)\binom z{j+1}.
\]
It verifies nonnegative coefficients in this basis at orders
\(t^5,\ldots,t^{24}\). Thus the twenty highest homogeneous parts of
\(S\) have nonnegative height coefficients for every admissible support.
In particular, the negative shifted-gap monomial at order \(t^{12}\)
is compatible with positivity in the binomial basis.

Independently, the checker forms all 190 complete polynomial
differences \((\Delta_A^i\Delta_B^j Q)_{0,0}\) with \(i+j\le18\).
It reverses the fully expanded remainders and clears a single common
denominator \(432^{20}\) before differencing. Every coefficient in
these complete differences is nonnegative.

A sufficient recurrence target for the uniform height-coefficient
sign is coefficientwise nonnegativity of the first forward differences
\(\Delta_A Q\) and \(\Delta_B Q\) at every gap pair. This weaker
condition than positivity of all Newton coefficients would propagate
the height coefficient signs from \(Q_{0,0}\). In the original heights
the two targets are coefficientwise nonnegativity of
\[
 \begin{aligned}
   S(L,H;a+1,d+1)-432L^{11}H^3 S(L,H;a,d),\\
   S(L,H;a,d+1)-432L^7H^7 S(L,H;a,d).
 \end{aligned}
\]
The normalizing powers change by \((11,3)\) and \((7,7)\),
respectively. It suffices to prove the first inequality on the edge
\(d-a=1\) and the second everywhere. One can first raise \(a\) along
that edge and then raise \(d-a\) with \(a\) fixed.
The origin proposition below supplies strict positivity once
the coefficientwise signs have been proved.

## Common factors in the gap increments

Put \(b=d-a\). The greatest common product of primitive linear
factors across the expanded products of summands in \(S\) is
\[
 K_{a,d}=(L+1)_d^3(L+1)_a^4(H+1)_d(H+1)_b^2
          \prod_{\substack{1\le j\le a\\ j\ {\rm odd}}}(2L+j).
 \tag{5}
\]
Indeed, the \(L\)-integer factors in (2) occur three times with
length \(d\) and four times with length \(a\) in every product.
The remaining second-order \(L\)-factors have half-integer shifts.
Their intersection across the three products comes from \(B\)
in \(ABG\) and is exactly the odd \(j\le a\) in (5).
The third-order noninteger shifts occur in thirds and the term
\(C^2D\) omits them. In \(H\), the lengths in \(C^2D\) are
\(b,b,d\). The other products contain at least the same factors.
This establishes (5) as a termwise common product.

Write \(P_{a,d}=S/K_{a,d}\). For the step increasing \(b\) at fixed \(a\),
the new common product divided by the old is
\[
 (L+d+1)^3(H+d+1)(H+b+1)^2 .
\]
For the step increasing \(a\) at fixed \(b\), the ratio is
\[
 (L+d+1)^3(L+a+1)^4(H+d+1)(2L+a+1)^\epsilon,\qquad
 \epsilon=\mathbf 1_{a\ {\rm even}}.
\]
Retaining their leading monomials gives sufficient
coefficientwise bounds for the two increments,
\[
 \begin{aligned}
 P_{a,d+1}&\succeq 432 L^4H^4 P_{a,d},\\
 P_{a+1,d+1}&\succeq \frac{432}{2^\epsilon}L^{4-\epsilon}H^2 P_{a,d}
                      \qquad (b=1).
 \end{aligned}
 \tag{6}
\]
Here \(U\succeq V\) denotes coefficientwise nonnegativity of \(U-V\).
These inequalities propagate coefficientwise nonnegativity of
the quotients from \(P_{1,2}\), whose coefficients are checked
directly. Multiplication by (5) then gives the original signs.

In the [exact increment check](verify_sfc_homogeneous_trinomial_gap_increments.py)
through source span 34, even the following multipliers in place of
those in (6) work:
\[
 \rho_b=\frac{432(d+1)^3(b+1)^4(d+3)}{d^3b^4(d+2)},\qquad
 \rho_a=\frac{432(a+1)^2(d+1)^3(d+3)}{2^\epsilon a^2d^3(d+2)}.
 \tag{7}
\]
These are the ratios at the least \(L\) exponent of the respective
leading homogeneous parts, after the monomial shifts in (6).
Indeed, the factor \((L+H)^4\) is shared, and
\((d-1)L+(d+3)H\) minus \((d+3)/(d+2)\) times
\((d-2)L+(d+2)H\) is \(4L/(d+2)\).
The least-\(L\) coefficient gives an upper bound on any multiplier
for the full comparison. The checker verifies that (7) attains
this bound on 561 steps in each direction, checking both
increments at every source support in that range.
Uniform comparison in (6) remains open.

## An all-order factor in the formal gap parameters

Each reciprocal-height layer of \(Q_{a-1,b-1}\) is divisible,
as a polynomial in the gap parameters, by \(a^2d^3b^4\),
where \(d=a+b\). To see this, let
\(\Phi_m(z)\) be the formal normalized product for counts whose
total right offset is \(z\). Its factors have lengths \(md-z,z\)
and the power-sum formulas above give polynomial coefficients in
the lengths. Write \(q_s=\Phi_1(s)\), let \(T_s\) denote translation
in \(z\) by \(s\), and put
\[
 {\cal D}_s=q_0 T_s-q_s,\qquad
 \widehat C_{m,k}=\binom mk
   ({\cal D}_d^{\,m-k}{\cal D}_b^k\Phi_m)(0).
\]
The second formula is exactly the normalized signed count sum
in (1). The operators commute. Since
\({\cal D}_s=q_0(T_s-1)+(q_0-q_s)\) is divisible by \(s\),
\(\widehat C_{m,k}\) is divisible by \(d^{m-k}b^k\).
Each of the three products defining \(S\) has total \(X,Y\)
degrees \(3,4\), which gives \(d^3b^4\).

For the remaining factor, put \(u={\cal D}_d\) and
\(w={\cal D}_b-{\cal D}_d\), so the two original operators are
\(u,u+w\) and \(w\) is divisible by \(a=d-b\). Write \(A,U,V\)
for the evaluations of \(u^2,uw,w^2\) on \(\Phi_2\) at zero,
and \(D,P_1,P_2,P_3\) for those of
\(u^3,u^2w,uw^2,w^3\) on \(\Phi_3\).
In these variables the normalized \(S\) is
\[
 \bigl((2U+V)^2-AV\bigr)D
   -6A(U+V)P_1+3A(A-V)P_2+2A(A+U)P_3 .
\]
Each monomial uses at least two factors of \(w\).
The three linear gap polynomials \(a,d,b\) are pairwise
relatively prime, so their powers give the claimed product.
The factor persists after removing (5) from the normalized series,
with the gap parameters considered separately on each parity of
\(a\). In reciprocal heights the normalized product (5) has
constant term one. Its remaining coefficients are polynomial in
the gaps by the power-sum formula, including the odd shifts by
using length \(\lceil a/2\rceil\) on each parity. Its formal
inverse therefore has polynomial coefficients on each parity.

## Highest gap degrees and the finite test for each coefficient

For \(p,q\ge0\) set \(n=p+q\) and
\(F_{p,q}(a,b)=[x^py^q]Q_{a-1,b-1}(x,y)\).
These are polynomials in the two gaps. For \(n<5\) they vanish.
For every \(n\ge5\), \(F_{p,q}\) has gap degree \(2n\), and its
degree-\(2n\) part has a positive coefficient at every monomial
\(a^2b^4 a^i b^{2n-6-i}\) with \(0\le i\le 2n-6\).

To derive this highest part, replace the gaps by \(\lambda a,\lambda b\)
and the reciprocal heights by \(x/\lambda^2,y/\lambda^2\).
In degree \(n\) in the reciprocal heights, the gap degree is at most
\(2n\). Indeed, the \(r\)-th logarithmic coefficient of each product
of normalized rising factorials has gap degree at most \(r+1\),
with equality \(r+1=2r\) only at \(r=1\). Thus only the quadratic
part of the first logarithmic coefficient contributes to the
degree-\(2n\) part. For a row of order \(m\) with total right offset
\(z\), the resulting highest-gap series is
\[
 \Phi_m^{\rm top}(z)=
   \exp\!\left(\frac{(md-z)^2x+z^2y}{2m}\right).
\]
All subsequent identities are in formal power series in \(x,y\).

Use the three right offsets \(z_i\in\{d,b,0\}\). Dividing each
column by its first-order series
\(p_i=\exp(((d-z_i)^2x+z_i^2y)/2)\) gives the kernel
\[
 \frac{\Phi_m^{\rm top}(\sum_{\nu=1}^m z_{i_\nu})}
      {\prod_{\nu=1}^m p_{i_\nu}}
 =\exp\!\left(-\frac{x+y}{2m}
        \sum_{1\le\nu<\mu\le m}(z_{i_\nu}-z_{i_\mu})^2\right).
\]
Write \(s=x+y\). On the first-moment zero plane with coefficients
\((X,Y,-X-Y)\), the two diagonal quadratic coefficients are
\(2(1-e^{-d^2s/4})\) and \(2(1-e^{-b^2s/4})\).
The cubic on that plane is \(-3XY(X+Y)\det M\), where
\(M_{ij}=e^{-s(z_i-z_j)^2/6}\). For example, for coefficients
\(c_0+c_1+c_2=0\), the pure cubes sum to \(3c_0c_1c_2\)
and \(c_i^2c_j+c_ic_j^2=-c_0c_1c_2\); the resulting expression
is \(3c_0c_1c_2(1-\sum_{i<j}M_{ij}^2+2M_{01}M_{02}M_{12})\).
Consequently \(S\) for the normalized columns is
\(12(1-e^{-d^2s/4})(1-e^{-b^2s/4})\det M\).
Restoring the column scales in the eliminated chart multiplies
this by \((p_0p_2)^3(p_1p_2)^4\).

Factor \(e^{-s z_i^2/6}\) from each row and column of \(M\).
The determinant left over reduces at its zero node to
\((e^{sb^2/3}-1)(e^{sd^2/3}-1)-(e^{sbd/3}-1)^2\).
Expanding the two products and pairing distinct indices gives
the complete series of highest gap degrees,
\[
 12e^{\gamma_x x+\gamma_y y}
 (e^{d^2s/4}-1)(e^{b^2s/4}-1)
 \sum_{1\le i<j}
 \frac{(b^id^j-b^jd^i)^2}{i!j!}(s/3)^{i+j},
 \qquad
 \begin{aligned}
 \gamma_x&=(59a^2+70ab+28b^2)/12,\\
 \gamma_y&=(11d^2+17b^2)/12 .
 \end{aligned}
 \tag{8}
\]
In particular this expression uses the height reversal and the
\(432^d\) normalization of \(Q\). All of its coefficients in
\(a,b,x,y\) are nonnegative. After dividing by \(a^2b^4\) they
remain nonnegative because
\(b^id^j-b^jd^i=b^id^i(d^{j-i}-b^{j-i})\) and dividing
\(d^{j-i}-b^{j-i}\) by \(a=d-b\) leaves a positive sum.
The first term of each difference of exponentials and the term
\(i=1,j=2\) yield \(a^2b^4d^4s^5/72\) times
\(e^{\gamma_x x+\gamma_y y}\). Both quadratics \(\gamma_x,\gamma_y\)
have positive coefficients in every monomial of degree two in
the gaps. For every \(p+q=n\ge5\) the displayed contribution
therefore has all gap monomials asserted above after division
by \(a^2b^4\).

This determines a finite exact test for both gap increments of
each fixed reciprocal-height coefficient. Let \(c_{p,q}>0\) be the
smallest coefficient of the degree-\(2n\) part of \(F_{p,q}\)
among its monomials divisible by \(a^2b^4\).
Let \(N_{p,q}\) be the sum of the absolute values of its negative
coefficients of lower gap degree. Every negative monomial
\(a^i b^j\) has \(i\ge2,j\ge4\) by the formal gap factor.
Raise its exponent in the larger gap by \(2n-i-j\) to get a
monomial in the positive top part. Under a forward difference
in either gap, the larger monomial contributes at least
\(\max(a,b)^{2n-i-j}\) times the forward difference of the
smaller monomial at the point \((a,b)\). If the difference acts
on the raised exponent, use
\[
 (z+1)^{k+\ell}-z^{k+\ell}
 \ \ge\ z^\ell\bigl((z+1)^k-z^k\bigr)
 \qquad(z\ge1,\ \ell\ge0).
\]
If the difference acts on the other variable the factor simply
pulls out. Each negative contribution is thus bounded by
its absolute coefficient divided by \(c_{p,q}\max(a,b)\) times
the forward difference of the whole top part. It follows that
\[
 \Delta_a F_{p,q}(a,b)>0,\qquad \Delta_b F_{p,q}(a,b)>0
 \quad\text{if}\quad
 \max(a,b)>N_{p,q}/c_{p,q},\qquad a,b\ge1.
\]
A smaller cutoff follows by summing negative contributions
assigned to the same top monomial. Write \(f_{ij}=[a^ib^j]F_{p,q}\).
For \(a\ge b\), all negative terms with the same \(j\) are
assigned to \(a^{2n-j}b^j\). Both differences are positive if
\[
 \sum_{i:f_{ij}<0}\frac{-f_{ij}}{f_{2n-j,j}}\,
                  a^{-(2n-i-j)}<1
 \quad\text{for every such }j.
\]
For \(b>a\), group by \(i\) and use \(b\) and \(f_{i,2n-i}\)
instead. Each sum decreases to zero, so exact rational bisection
gives a cutoff in each of the two regions. The layer checker can
test both increments on the remaining finite triangles using
`--fixed-deficit`. This gives a terminating exact test of
either sign assertion for any prescribed reciprocal coefficient.
The cutoff depends on the reciprocal degree, leaving the two
uniform inequalities open. The fixed-deficit replay also certifies
both differences in reciprocal degree 25, extending the initial
band of whole-lattice increment signs by one degree. For that degree the maximum matched cutoffs are 109 and
49 in the respective regions. The exact check visits 34,933
coefficient-gap pairs.

## Strict sign at zero heights for every span

**Proposition.** \(S(0,0;a,d)>0\) for every \(0<a<d\) with integer
offsets.

Let \(k=\binom da\) and use the angular functions
\(f_0=(1-t)^d,\ f_1=k t^a(1-t)^b,\ f_2=t^d\) on \([0,1]\).
The contrasts are \(f_0-f_2,f_1-f_2\). At order \(m\) divide
by \(\int_0^1f_2^m=1/(md+1)\). The resulting remainder differs
from the original by a positive factor, since the chart
coefficients are rescaled by \(d!,d!/k\) and the row by
\((md)!\). Reflection gives \(D=0\), so it suffices to prove
\(BG-CF>0\). Put
\[
\begin{gathered}
 t_0=\binom{2d}{d}^{-1},\quad
 p=\frac{k}{\binom{2d}{a}},\quad r=\frac{k}{\binom{2d}{b}},\quad
 s=\frac{k^2}{\binom{2d}{2a}},\\
 x_0=\binom{3d}{d}^{-1},\quad z=\frac{k}{\binom{3d}{b}},\quad
 w_0=\frac{k^2}{\binom{3d}{2a}},\quad w=\frac{k^2}{\binom{3d}{2b}},\quad
 g=\frac{k^3}{\binom{3d}{3a}},\quad h=\frac{k}{\binom{3d}{d+a}} .
\end{gathered}
\]
The beta integral gives
\[
 B=2(1-t_0+p-r),\quad C=1+s-2r,\quad
 F=3(w_0-w-2h+2z+x_0-1),\quad G=g-3w+3z-1 .
\]
In particular \(B>0\) since \(t_0\le 2^{-d}\) and \(r\le 2^{-b}\).
With \(J=1+w-w_0-2z-x_0\), the inequality
\(BG-CF\ge BG+3CJ\) follows from \(h\ge0\) and \(C>0\).

For the large spans use \(d\ge16\) and \(e=1/8192\).
Termwise binomial ratios and Vandermonde's convolution give
\[
 \begin{gathered}
 p\le2^{-a},\quad r\le2^{-b},\quad z\le3^{-b},\quad
 x_0\le3^{-d},\\
 w_0\le s(2/3)^{2a},\quad
 w\le s(2/3)^{2b},\quad s\le1 .
 \end{gathered}
\]
For instance \(\binom{n}{j}/\binom{mn}{j}\le m^{-j}\)
by comparing descending factorial factors.
If \(c=a\) or \(b\), then
\[
 s=\frac{\binom{2c}{c}}{2^c}
        \prod_{i=0}^{c-1}\frac{d-i}{2(d-i)-1}
 \ \ge\ \frac{\binom{2c}{c}}{4^c}.
\]
For a gap at least 13, its mixed second coefficient
and its corresponding \(w_0\) or \(w\) are at most \(e\),
as is \(z\) when \(b\ge13\). Also \(x_0\le e\).
For a gap of one the corresponding \(w_0\) or \(w\) is
\(2d/[3(3d-1)]\) in \([2/9,23/100]\).
When \(b=1\), \(r=1/2,\ z=1/3,\ g\ge2/9\), the last
bound following by comparing each factor of \(\binom{3d}{3}\)
with \(3d\). When \(b=2\), \(r=(d-1)/(2(2d-1))\ge1/5\)
and the product for \(s\) gives
\(w\le(2/3)^4(3/2)(15/29)^2<2/25\).

The following rational bounds use these estimates:
\[
\begin{array}{c|c|c|c|c}
 & C\ {\rm lower} & J\ {\rm lower} & B\ {\rm upper} & -G\ {\rm upper}\\ \hline
 a,b\ge4 & 7/8 &187/200 &17/8 &28/25\\
 a=1 &3/2-2e &77/100-3e &3 &1+3e\\
 a=2,3 &21/16-2e &65/81-3e &5/2 &1+3e\\
 b=1 &1/2 &5/9-2e &1+2e &421/900\\
 b=2 &7/8 &7/9-2e &8/5+2e &31/25\\
 b=3 &17/16 &25/27-2e &2+2e &307/243
\end{array}
\]
In the first row \(w,w_0<1/25\) and \(2z+x_0<1/40\).
The other rows use the distant-gap bounds and the indicated
specific near-gap ratios. In each row, three times the
product of the two lower bounds exceeds the product
of the two upper bounds. The increment checker verifies these
rational margins and the 105 oriented pairs of span \(2\) through
15 by the same beta formulas. This proves the proposition.

## Integral and shift formulas

There is also a pointwise [sign theorem for separated angular exponents](SFC_HOMOGENEOUS_ANGULAR_SEPARATION.md).
In particular, \(R<0\) whenever the smaller gap is at least
\(6\sqrt{L+H+d+2}\).

Let \(N=L+H+d\), \(e_0=0,e_1=a,e_2=d\) and
\[
 \phi_j(t)=t^{L+e_j}(1-t)^{H+d-e_j},\qquad
 Z_j=\int_0^1\phi_j(t)\,dt,\qquad
 f_j=\phi_j/Z_j .
\]
The change of variables \(x=rt,\ y=r(1-t)\) separates the radial
factor \(\Gamma(mN+2)\) from the angular power integral. Liu and
Sun use this angular reduction in
[Theorem 2.6](https://doi.org/10.1017/S0004972719000546),
which assumes vanishing at every positive order. The
[Strong Factorial Conjecture of Edo and van den Essen](https://arxiv.org/abs/1304.3956)
asks about support-sized consecutive windows. The result here
concerns the initial three orders.

Set \(u=f_0-f_2,\ v=f_1-f_2\) and write
\(\alpha=\int u^2,\ \beta=\int uv,\ \gamma=\int v^2\), with integrals
over \([0,1]\). Since \(p_j/p_2=Z_j/Z_2\), the eliminated angular
phase is \(p_2(Z_0Xu+Z_1Yv)\). Thus \(S\) equals the integral
\[
 \int_0^1\bigl(\gamma^2u^3-3\alpha\gamma uv^2+
                    2\alpha\beta v^3\bigr)\,dt
 \tag{4}
\]
times the positive factor
\[
 (p_2 Z_0)^3(p_2 Z_1)^4
 \left(\frac{\Gamma(2N+2)}{\Gamma(2L+1)\Gamma(2H+1)}\right)^2
 \frac{\Gamma(3N+2)}{\Gamma(3L+1)\Gamma(3H+1)}.
\]
Formula (4) applies to real nonnegative heights as well.

Its integrand takes both signs. Write \(z=f_1/f_2\).
Then \(f_0/f_2=c z^{d/b}\) for a constant \(c>0\), and
\(\int f_2 z=1\). Strict Jensen gives
\(c=1/\int f_2 z^{d/b}<1\).
The ratio \(z\) decreases from infinity to zero, so \(v=0\) at a
unique interior point where \(u=f_2(c-1)<0\). The integrand there is
negative. Near \(t=0\), \(u>0\) and \(v/u\) tends to zero, so
the integrand is positive there. For \(L>0\), adjusting the
integrand by subtracting \(c_1u+c_2v\) cannot make it nonnegative.
Dividing the adjustment by \(u\) near zero, where \(u\) tends to
zero, would require \(c_1\le0\). At the interior zero of \(v\) it
would require \(c_1\ge\gamma^2u^2>0\).

The exact kernel \(W_m(U,V)=(mL+1)_U(mH+1)_V\), where
\(U=a\alpha_1+d\alpha_2,\ V=md-U\), has the shift multipliers
\[
 \begin{array}{c|c}
 L\mapsto L+1 &(mL+U+1)_m/(mL+1)_m\\
 H\mapsto H+1 &(mH+V+1)_m/(mH+1)_m\\
 a\mapsto a+1 &(mL+U+1)_{\alpha_1}/(mH+V-\alpha_1+1)_{\alpha_1}\\
 d\mapsto d+1 &(mL+U+1)_{\alpha_2}(mH+V+1)_{m-\alpha_2}
 \end{array}
\]
on each admissible shift. They follow by cancelling the rising
factorials in the shifted count weight. Their dependence on the
count vector gives different multipliers to the signed summands of
\(S\), leaving a positivity-preserving recurrence for \(S\) open.
The unrestricted inequality and the lower-degree coefficient
signs beyond the checked spans also remain open.

## Exact Beta Gram determinants and the fractional-gap obstruction

The Beta kernel does have an exact positive determinant before
taking the highest-gap approximation. For real \(d>0\),
\(L,H\ge0\), and distinct \(e_i\in[0,d]\), put
\[
 K(e,f)=B(2L+e+f+1,\,2H+2d-e-f+1).
\]
With \(\phi_e(t)=t^{L+e}(1-t)^{H+d-e}\), it satisfies
\[
 \det\bigl(K(e_i,e_j)\bigr)_{i,j=0}^2
 =\frac1{3!}\int_{[0,1]^3}
       \det\bigl(\phi_{e_i}(t_j)\bigr)_{i,j=0}^2^{\,2}
            \,dt_0dt_1dt_2>0.                              \tag{9}
\]
Expand both determinants, integrate term by term, and relabel
the three integration variables to obtain the identity.
Strictness follows from linear independence: after division
by the common positive envelope and substitution
\(x=t/(1-t)\), the three functions are distinct powers
\(x^{e_i}\). Their Wronskian is a nonzero Vandermonde
factor times a power of \(x\). All integrals converge.
This is the Gram specialization of Andréief's formula; see
[Forrester, equation (1.7) and Section 2.2](https://arxiv.org/abs/1806.10411),
which states and proves the integral identity. The argument
above supplies its application here. Dividing each function
by its first integral preserves the determinant sign.

However, the positive exact Gram determinant does not itself
give the sign of \(S\). Allow fractional gaps just for this
test, and define \(S_{\rm ang}\) by (4), using
\(f_i=\phi_{e_i}/\int\phi_{e_i}\) and
\((e_0,e_1,e_2)=(0,a,d)\). At
\((L,H,a,d)=(1,0,1/2,1)\), direct Beta integrals give
\[
 S_{\rm ang}=
 \frac{9(143613\pi^4-996096\pi^2-4194304)}{448000\pi^4}<0,
 \qquad
 \det\bigl(\int f_if_j\bigr)
   =\frac{27(256-25\pi^2)}{400\pi^2}>0.                    \tag{10}
\]
The failure persists at integer span two:
\[
 (L,H,a,d)=(19,0,3/2,2),\qquad
                 -4/3<S_{\rm ang}<-13/10.                  \tag{11}
\]
The Gram determinant is positive there as well. Both examples
have a gap less than one, so neither is an admissible ordinary
trinomial support for the theorem at the start of this note.
They also are not assertions of common zeros of the quadratic
and cubic: \(S=0\) is only a necessary common-zero condition.

The [exact Beta checker](verify_sfc_initial_beta_faces.py)
derives these expressions from Gamma values by recurrence and
certifies their signs with rational interval arithmetic.
For (11), half-integer Gamma values put every coefficient in
\(\mathbb Q(\pi)\). Bounds for \(\pi\) follow from
\(\pi=16\arctan(1/5)-4\arctan(1/239)\) and the alternating
arctangent series, using 18 terms and its next-term bound.
The identity follows by the tangent addition formula, with
the angles in the first quadrant. Polynomial interval
evaluation then verifies the stated rational enclosure
without floating point.

Thus an exact positive lift of (8) needs more than positivity
of the Beta Gram kernel on arbitrary ordered real nodes.
The exact mixed cubic moments must also be controlled, and
the admissible gap restrictions must enter that argument.
The coefficientwise gap inequalities (6) and the initial
integer-support sign remain open.

For the alternative translation route, the
[initial fractional-face theorem](SFC_ROW_TRANSLATION_DIVISOR_SIEVE.md#initial-windows-with-arbitrary-deleted-rows)
removes the fractional/Gamma nonvanishing obstacle at deleted
orders \(1,2\), even when the two columns differ in both rows.
Its divisor lists still need to be eliminated. It supplies no
unrestricted three-moment sign theorem.

## Replay

Run from the repository root with Python 3.12.3, python-flint 0.9.0
and SymPy 1.14.0. The internal replay ran spans 2 through 64 and
65 through 80 as two exact batches.

~~~sh
.venv/bin/python research/factorial/image-rigidity/verify_sfc_homogeneous_trinomial_factored_sign.py --max-span 80 --jobs 6
~~~

~~~sh
.venv/bin/python research/factorial/image-rigidity/verify_sfc_homogeneous_trinomial_top_layers.py
~~~

~~~sh
.venv/bin/python research/factorial/image-rigidity/verify_sfc_homogeneous_trinomial_gap_increments.py
~~~

~~~sh
.venv/bin/python factorial/image-rigidity/verify_sfc_homogeneous_trinomial_top_layers.py --fixed-deficit 25
~~~
