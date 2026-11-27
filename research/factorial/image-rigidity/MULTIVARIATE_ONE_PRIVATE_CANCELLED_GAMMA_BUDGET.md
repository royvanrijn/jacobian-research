# Exact Gamma cancellation in the half at a private endpoint

Consider an ordinary two-monomial profile with private rows on
exactly one side. Orient the coefficient index \(k\) from the side
having the private variables. Its private rows are then \((e_i,0)\)
with total \(P=\sum_i e_i\). The remaining unequal rows have
\(a_i,b_i>0\), with \(d_i=b_i-a_i\). Define
\[
 W_{\rm eff}=P-1+\sum_{\rm sh}\frac{d_i^2}{a_i},\qquad
 L=\max(1,\max_{\rm sh}|d_i|/a_i),\qquad
 S_{\rm sh}=\sum_{\rm sh}\frac{d_i^2}{\max(a_i,b_i)} .
\]
An empty shared maximum in \(L\) is omitted. These are
orientation-specific quantities at the endpoint with private
exponents.

For real \(m\ge n\ge2\) use the Gamma interpolation
\(c_{m,x}\) of the endpoint-normalized binomial coefficient
array. Put \(B_m=c_{m,1}\) and
\(f_m(x)=\log(c_{m,x}\Gamma(1+x)/B_m^x)\).
For an integer \(2\le K<n\), write \(\eta=n-K-1/2>0\).
For all integral \(r\) with \(2\le r\le K\), the forward
difference \(\nu_r(m)=\Delta_x^r f_m(0)\) obeys
\[
 |\nu_r(m)|\le H_r,\qquad
 |\partial_m\nu_r(m)|\le \frac{r-1}{\eta} H_r,\qquad
 H_r=(r-2)!\left(\frac{L}{\eta}\right)^{r-2}
                          \frac{W_{\rm eff}}{\eta}.
\tag{1}
\]
The same derivative bounds hold for \(\partial_x^r f_m(x)\)
and \(\partial_m\partial_x^r f_m(x)\) on \(0\le x\le K\).
In particular the direct head estimate has the signed form
\[
 -\frac{k(k-1)W_{\rm eff}}{2\eta^2}
       \le \partial_m f_m(k)\le0 \qquad (2\le k\le K).
\tag{2}
\]
If every shared gap is nonpositive, then every \(\nu_r\)
is nonnegative and nonincreasing in \(m\).

## Cancellation and the derivative bounds

Gauss multiplication writes
\[
 \Gamma(1+e_i y)
  =C_i e_i^{e_i y}\prod_{j=1}^{e_i}\Gamma(y+j/e_i)
\]
for a positive \(C_i\) independent of \(y\). Apply this to
\(y=m-x\) for the private rows. One of their shape-one factors
cancels the \(\Gamma(1+m-x)\) in the binomial denominator.
The other binomial denominator \(\Gamma(1+x)\) is already
removed in \(f_m\). Thus
\[
 f_m(x)=\text{an affine function of }x
          +\sum_{\gamma}\log\Gamma(m-x+\gamma)
          +\sum_{\rm sh}\log\Gamma(1+ma_i+d_i x),
\tag{3}
\]
where the first sum has \(P-1\) terms with \(0<\gamma\le1\).
Its affine part comes from normalization and gives
\(f_m(0)=f_m(1)=0\). This identity holds for real degree and real
index on the original coefficient interval.
[Gauss multiplication formula](https://dlmf.nist.gov/5.5.E6)

For \(y>1/2\) and \(r\ge2\), the positive series for the
absolute polygamma value, bounded by the centered integral
sum for a decreasing convex inverse power, gives
\[
 |\psi_{r-1}(y)|
 =(r-1)!\sum_{j\ge0}(y+j)^{-r}
 \le (r-2)!(y-1/2)^{1-r}.
\]
The same bound with \(r\) replaced by \(r+1\) applies to
\(\psi_r\). A shared Gamma argument in (3), minus \(1/2\),
is at least \(a_i\eta\). This uses \(d_i\ge0\) for a
nondecreasing argument and \(|d_i|\le a_i\) when \(d_i<0\).
A private residual argument minus \(1/2\) is at least \(\eta\).
For a shared row the resulting \(r\)-th derivative is bounded
by
\((r-2)!(|d_i|/a_i)^{r-2}(d_i^2/a_i)/\eta^{r-1}\).
Differentiating also in \(m\) multiplies the corresponding
bound by at most \((r-1)/\eta\), since the extra derivative
of the argument is \(a_i\). The private residual argument
has slope \(-1\) in \(x\) and slope \(1\) in \(m\), giving
the same bounds with \(a_i=|d_i|=1\).
Summing and using the iterated fundamental theorem of calculus
for each forward difference proves (1).
[Trigamma series](https://dlmf.nist.gov/5.15.E1)

For \(r=2\) every term of \(\partial_m f_m''\) is
nonpositive by the sign of \(\psi_2\). Thus each second
forward difference of \(\partial_m f_m\) is nonpositive
and has absolute value at most \(W_{\rm eff}/\eta^2\).
The two zero values at \(0,1\) and discrete summation prove
(2). When the shared gaps are all nonpositive, the sign
of each \(\psi_{r-1}\) cancels the sign of the \(r\)-th
power of each slope in (3). The extra degree derivative
has the opposite sign. Integrating gives the asserted
signs of \(\nu_r\) and \(\partial_m\nu_r\).

## Tight endpoint budgets for the interaction terms

Use the oriented endpoint envelope \(Q(w)\) with all shared
rows retained and any nonnegative private chord parameter
\(s_0\). Set \(q=Q(0)\) and \(S'=s_0+S_{\rm sh}\).
For \(h=1,2,3\), put
\[
 M_1=4,\quad M_2=8,\quad M_3=13,\qquad
 C_h=\max(M_h,P-1+S_{\rm sh}).
\]
Then the actual interaction budget satisfies
\[
                       W_{\rm eff}^{\,h}q
                          \le C_h^{\,h} e^{-S'/2}.
\tag{4}
\]
The envelope definition is the same as the one used for
the [private-side entropy theorem](MULTIVARIATE_ONE_PRIVATE_ENTROPY_THIRTEEN.md#matching-the-long-phase-hypotheses).
An actual coefficient bound furnishing the chosen \(s_0\) is
a separate input to the phase argument.

In the [directed endpoint notation](MULTIVARIATE_ABSOLUTE_CURVATURE_13.md#directed-endpoint-estimates),
an increasing shared row with \(\ell=a_i\), \(r=d_i/\ell>0\)
contributes excess \(V_{d,i}-S_i=\ell g(r)\) where
\(g(r)=r^3/(1+r)\), and it pays the excess logarithmic cost
\(\ell E_+(r)\) in \(q\), with
\[
 E_+(r)=(1+r)\log(1+r)-r-\frac{r^2}{2(1+r)} .
\]
A decreasing row has \(V_{d,i}=S_i\) and pays at least its
baseline cost \(S_i/2\). For any \(h\) above,
\[
                E_+'(r)\ge\frac{h\,g'(r)}{M_h+g(r)} .
\tag{5}
\]
For \(h=3\) this is the cited directed cubic bound.
For \(h=1,2\), use \(\log(1+r)\ge2(t+t^3/3)\),
\(t=r/(r+2)\), to obtain
\[
 E_+'(r)\ge
 \frac{r^2(3r+4)}{2(r+2)(r+1)^2}+
 \frac{2r^3}{3(r+2)^3}.
\]
After subtracting the right side of (5) and clearing positive
factors, the respective polynomials are
\[
 \begin{split}
 &48+88r+12r^2-18r^3+38r^4+44r^5+13r^6,\\
 &96+176r+24r^2-84r^3-12r^4+32r^5+13r^6.
 \end{split}
\]
On \(0\le r\le1\), the combined negative terms are absorbed
by the preceding positive linear and quadratic terms. After
\(r=1+s\), the coefficient lists in increasing powers of \(s\)
are \(225,508,821,834,453,122,13\) and
\(245,162,215,448,343,110,13\), respectively. This proves (5).

Since \(C_h\ge M_h\), integrate (5) with \(C_h\) in place
of \(M_h\) to get \(e^{E_+(r)}\ge(1+g(r)/C_h)^h\).
Applying the integer-power Bernoulli inequality on each
increasing row and then multiplying over the rows gives
\[
 q\le e^{-S'/2}
       \left(1+\sum_{\rm incr}\ell_i g(r_i)/C_h\right)^{-h}.
\]
The identity \(W_{\rm eff}=P-1+S_{\rm sh}+\sum_{\rm incr}
\ell_i g(r_i)\) now proves (4).

In the [direct head comparison](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md#direct-estimates-in-one-orientation),
(2) controls the part of \(\partial_m\log c_{m,k}\) beyond
\(k\partial_m\log B_m\). In the
[connected comparison](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md#connected-estimates-in-one-orientation),
(1) supplies the edge and degree-derivative inputs, and
(4) controls the product of their budget with \(D=Nq\).
Both arguments also require an actual coefficient envelope,
a first-factor bound and a final phase budget. For \(P=5\), the linear bound simplifies to
\(W_{\rm eff}q\le4e^{-s_0/2}\), uniformly over all shared rows,
by the decrease of \((4+s)e^{-s/2}\) on \(s\ge0\).
The [private-degree-five theorem](MULTIVARIATE_ONE_PRIVATE_DEGREE_FIVE_ALL_ANGLES.md)
uses (4) at \(s_0=0\) with the exact private first-factor loss,
and derives its coefficient tail from the cancelled interpolation.

