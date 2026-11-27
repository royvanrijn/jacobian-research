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
\(t^5,\ldots,t^{16}\). Thus the twelve highest homogeneous parts of
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
the height coefficient signs from \(Q_{0,0}\). Strict positivity of
\(S(0,0)\) at each support would then give strict positivity at every
nonnegative height pair.

## Integral and shift formulas

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
