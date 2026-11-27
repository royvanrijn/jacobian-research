# Two consecutive factorial moments near a mixed quadratic anchor

Let \(K\) be a complete nonarchimedean field of characteristic zero
with a real-valued valuation extending \(v_p(p)=1\). Write
\(\Lambda(X^\nu)=\prod_i\nu_i!\) and \(v(0)=+\infty\).
Put \(\lambda=1/(p-1)\). In two variables, suppose
\[
                        F=c(xy+P),\qquad c\ne0,
              \qquad v([x^iy^j]P)\ge 2\lambda
                                  \quad\text{for all }i,j .          \tag{1}
\]
Then at least one of \(\Lambda(F^n)\) and \(\Lambda(F^{n+1})\)
is nonzero for every \(n\ge0\). More precisely,
\[
        \min_{j=0,1}v\left(\frac{\Lambda(F^{n+j})}
                                   {c^{n+j}(n+j)!}\right)
                      \le \lambda n .                               \tag{2}
\]
There is no restriction on the degree or the perturbation support
in (1). The result includes the boundary \(2/(p-1)\). For
coefficients in a noncomplete valued field one may pass to its
completion because the ordinary moments are unchanged under that
embedding.

The proof first gives a bound for a constant plus \(xy\) by a
determinant calculation with the finite truncations of its
contiguous Bessel series. It then shows that adding the same
constant to a phase with suitably normalized analytic moments
preserves the bound. A translation to the small critical point
and the analytic preparation of a mixed product supply those
normalized moments for (1). The same transfer shows that the
cost of a pole in the original monomial expansion can be much
smaller than the valuation of every low jet minor.

Write \(s_p(m)\) for the base-\(p\) digit sum. We use the
factorial valuation \(v(m!)=\lambda(m-s_p(m))\), as in
[Granville, Section 2](https://www.cecm.sfu.ca/organics/papers/granville/paper/binomial/html/node3-an.shtml).
The series and recurrence in Section 1 use the same consecutive
residues as the [constant-plus-product test for towers](SFC_LOGARITHMIC_TOWER_RESIDUES.md#8-a-test-of-the-window-transfer-on-a-product-plus-a-constant).
The determinant comparison also applies at the convergence
boundary for that series. The all-order result for the two-term
phase over unrestricted coefficient fields is already established
in [the constant-plus-product theorem](CONSTANT_PLUS_PRODUCT_OF_TWO_VARIABLES_ALL_ORDERS.md).

## 1. A quantitative reference pair, including the boundary

Fix \(t\ne0\) with \(\alpha=v(t)\ge2\lambda\), and put
\[
    e_n=\frac{t^n}{n!},\qquad
    V_n=\sum_{k=0}^n\frac{t^k}{k!}(n-k)!
       =\frac{\Lambda((t+xy)^n)}{n!}.
\]
For \(n\ge1\), splitting \(n+1-k=(n+1)-k\) gives
\[
                 V_{n+1}=(n+1)V_n-t V_{n-1}+e_{n+1}.                \tag{3}
\]
For \(T\ge0\) and \(s\ge0\) write
\[
 Z_s^{[T]}=\sum_{j=0}^T\frac{(-t)^j}{j!(s+j)!},\qquad
 Y_n^{[T]}=t^{n+1}Z_{n+1}^{[T]} \quad(n\ge0).
\]
Multiplying each summand of \(Z_s^{[T]}\) by \(s+j\) and
splitting the two contributions gives the contiguous recurrence
\[
 Z_{s-1}^{[T]}-s Z_s^{[T]}+t Z_{s+1}^{[T]}
              =\frac{t(-t)^T}{T!(s+1+T)!}\quad (s\ge1).
\]
With \(D_n^{[T]}=V_nY_{n+1}^{[T]}-V_{n+1}Y_n^{[T]}\), (3) gives
\[
 \frac{D_n^{[T]}}{t^{n+1}}
  =\frac{D_{n-1}^{[T]}}{t^n}
      -e_{n+1}Z_{n+1}^{[T]}+
        \frac{V_n t(-t)^T}{T!(n+2+T)!}\quad(n\ge1).                  \tag{4}
\]
The initial value \(D_0^{[T]}/t\) is
\(-Z_0^{[T]}-tZ_1^{[T]}+t(-t)^T/(T!(T+2)!)\).
It is a unit. Indeed \(Z_0^{[T]}\) is congruent to one by
the digit-sum savings for \(j\ge1\). Both other terms have
strictly positive valuation by the same factorial formula.
Every term of \(e_{n+1}Z_{n+1}^{[T]}\) also has positive
valuation, since its total factorial argument is
\(2(n+1+j)\) and \(s_p(n+1)>0\).

The valuation of the last term in (4), if nonzero, is at
least \(v(V_n)-\lambda n+\lambda s_p(T)\). For each fixed
range \(1\le n\le N\), choose \(T\) with digit sum large enough
that all these terms have positive valuation. Such \(T\) exist,
for example \(T=p^a-1\) as \(a\) grows. Thus
\[
                     v(D_n^{[T]})=\alpha(n+1)\quad(0\le n\le N).
\]
The factorial estimates also give, uniformly in \(T\),
\[
             v(Y_n^{[T]})\ge(\alpha-\lambda)(n+1)+\lambda .
\]
Indeed the total factorial argument of each term is
\(n+1+2j\), and \(s_p(n+1+j)\ge1\). Taking the valuation of
the determinant now gives
\[
                         \min(v(V_n),v(V_{n+1}))\le\lambda n .        \tag{5}
\]
The finite truncation is important when \(\alpha=2\lambda\).
If \(\alpha>2\lambda\), the series converge for every fixed
index. Then (4) gives \(v(D_n)=\alpha(n+1)\) for the genuine
series as well, by passing to the limit at each fixed \(n\).

## 2. Adding a constant to analytically normalized moments

Here the number of variables is arbitrary. Suppose a phase \(Q\)
satisfies
\[
              \frac{\Lambda(Q^m)}{(m!)^2}=A(m)+E_m\quad(m\ge0),       \tag{6}
\]
where \(A\) is a restricted power series on the closed unit disk
and the infimum of the coefficient valuations of \(A-1\) is
strictly positive. Suppose
\[
                             v(E_m)>\lambda s_p(m).                 \tag{7}
\]
Then the bound (2) with \(c=1\) holds for \(F=t+Q\) whenever
\(v(t)\ge2\lambda\), allowing \(t=0\). If \(t=0\), (6) and
(7) give \(v(\Lambda(Q^m)/(m!)^2)=0\), so \(v(\Lambda(Q^n)/n!)=
v(n!)\le\lambda n\).

For \(t\ne0\), expand \(A\) in the basis
\(P_h(z)=\prod_{i=1}^h(z+i)\), \(P_0=1\):
\[
                              A(z)=\sum_{h\ge0}d_hP_h(z).
\]
Both changes of basis between the \(P_h\) and monomials have
integer coefficients because the \(P_h\) are monic integral
polynomials of degree \(h\). Applying the change to polynomial
truncations of \(A\) proves convergence in the restricted
power series norm, with \(d_h\to0\), \(d_0\equiv1\) and
\(d_h\equiv0\) for \(h>0\) modulo the maximal ideal.

First use the model \(E_m=0\). The normalized moments of
\(t+Q\) in this model are
\[
   \mathcal W_n=\sum_{k=0}^n\frac{t^k}{k!}(n-k)!A(n-k)
                =\sum_{h\ge0}d_h T_{n,h},\qquad
   T_{n,h}=\sum_{k=0}^n\frac{t^k}{k!}(n+h-k)! .                     \tag{8}
\]
The shift basis has \(T_{n,0}=V_n\),
\(T_{n,1}=V_{n+1}-e_{n+1}\), and for \(h\ge2\),
\[
      T_{n,h}=(n+h)T_{n,h-1}-t T_{n,h-2}+t e_n(h-2)! .
\]
This follows by splitting \(n+h-k\) as in (3), with the
last term accounting for the last removed summand. Thus
\(T_{n,h}\) is an integral polynomial combination of
\(V_n,\ V_{n+1}-e_{n+1},\ t e_n\), with coefficients in
\(\mathbb Z[n,t]\). At \(h=0,1\) the coefficient vectors
are \((1,0,0)\) and \((0,1,0)\). Summing using \(d_h\) gives
restricted series \(u,v_0,h_0\) in the order variable with
\(u\equiv1\), \(v_0,h_0\equiv0\) modulo the maximal ideal
in coefficient norm and
\[
   \mathcal W_n=u(n)V_n+v_0(n)(V_{n+1}-e_{n+1})+t h_0(n)e_n .
\]
Apply (3) to \(\mathcal W_{n+1}\) to obtain the joint identity
\[
 \binom{\mathcal W_n}{\mathcal W_{n+1}}
   =
   \begin{pmatrix}
      u(n)&v_0(n)\\
      -t v_0(n+1)&u(n+1)+(n+2)v_0(n+1)
   \end{pmatrix}\binom{V_n}{V_{n+1}}
   +
   \binom{-v_0(n)e_{n+1}+t h_0(n)e_n}{t h_0(n+1)e_{n+1}} .           \tag{9}
\]
The matrix and its inverse are integral and hence preserve the
minimum of the two valuations. The error has both valuations
strictly greater than \(\lambda n\), since
\(v(e_m)\ge(\alpha-\lambda)m\) and \(\alpha\ge2\lambda\).
Use (5) to deduce that the model pair has the same minimum
valuation as the reference pair.

For the actual \(Q\), the additional terms in component \(n\)
are \(\sum_{k+m=n} t^k m! E_m/k!\), and analogously at \(n+1\).
Every summand has valuation strictly greater than
\(\lambda n\) at order \(n\), because
\(v(m!E_m)>\lambda m\) and \(v(t^k/k!)\ge\lambda k\).
This proves (2) and the exact preservation of the minimum
valuation of the reference pair when \(t\ne0\).

## 3. Analytic factors near separate coordinate anchors

The assumptions (6)–(7) hold in any number of variables for
\[
       Q=(1+P_0)(X_1+P_1)(X_2+P_2),\qquad
       \inf_\nu v([X^\nu]P_i)\ge 2\lambda,                          \tag{10}
\]
with each \(P_i\) a restricted ordinary power series.
For a restricted power series, factorial evaluation is continuous
because all \(\nu!\) are \(p\)-integral.

Expand the three factors separately. Let \(K_i\) be the number
of perturbations from factor \(i\), and \(K=K_0+K_1+K_2\).
Relative to \((m!)^2\), row 1 can lose at most \(K_1\)
consecutive factors from \(m!\), since the only negative
contributions to its exponent come from departures from
the \(X_1\) anchor. These factors all cancel against
\((m)_{\underline{K_1}}\). Row 2 is treated with its own
\(K_2\). Other rows have fixed nonnegative exponents and
contribute integral scalars. After cancellation every count
term is a polynomial of the order variable times the count
coefficient. The polynomial has integral coefficients before
dividing by the count factorials. Thus each nonzero count
term has coefficient valuation at least \(2\lambda K-
\lambda K=\lambda K\).

For fixed \(K\), only finitely many terms can be below any
given precision, because the coefficients of the perturbation
series tend to zero. Therefore the sum of all the count
polynomials, including the all-anchor term one, defines \(A\)
and \(A-1\) has coefficient valuations at least \(\lambda\).
At the ordinary order \(m\), only \(K_i\le m\) for all \(i\)
are used in the actual expansion. All other count terms have
\(K\ge m+1\), giving the stronger error bound
\[
                         v(E_m)\ge\lambda(m+1)>\lambda s_p(m).       \tag{11}
\]
This is the required normalized approximation.
For \(Q=X_1X_2(1+P_0)\), all \(v([X^\nu]P_0)>\lambda\) with
uniform strict margin also suffice, because the falling
factorial is unaltered and excess counts vanish exactly at
each ordinary order.

## 4. Every small bivariate perturbation

Normalize \(c=1\). The map \((x,y)\mapsto
(-\partial_yP(x,y),-\partial_xP(x,y))\) is a contraction
on \(v(x),v(y)\ge2\lambda\). Its images lie in that ball.
The difference estimate follows by subtracting monomials
at two input pairs and using integral exponents, so it
increases the minimum valuation distance by at least
\(2\lambda\). By completeness there is a stationary point
\(h\) in the ball. The value \(t=F(h)\) has valuation at
least \(2\lambda\). In variables translated by \(h\),
\(G(X)=F(X+h)-t\) has zero constant and linear terms, and
\(G-xy\) has coefficient valuations at least \(2\lambda\).

Apply the [prepared factorization of a translated product](SFC_POWERED_RADIAL_FACTOR_CLASSIFICATION.md#7-analytic-factors-of-a-translated-product-of-coordinate-powers)
to \(G\) with \(a_1=a_2=1\) and \(w_{12}=w_{21}=1\).
Its support inequalities both reduce to total degree at
least two. The possible small correction to the coefficient
of \(xy\) satisfies the same inequalities. The preparation
gives \(G=W_0W_1W_2\) with \(W_0-1,W_1-x,W_2-y\) of
coefficient valuation at least \(2\lambda\).
Translation back by \(-h\) preserves those bounds and yields
\(F-t\) in the form (10). Sections 2–3 prove (2).

The argument uses a scalar value at a stationary point in
two variables. With additional independent variables,
Section 3 still applies to phases admitting (10) and a
scalar constant \(t\). Preparing a general perturbation
with additional variables can instead leave a function of
those variables in place of the scalar. The transfer above
does not evaluate that case.

## 5. The moving-center comparison for upward perturbations

For \(F=t+xy(1+P_0)\) in two variables, assume \(P_0\) is a
polynomial with zero constant term and uniform coefficient valuations
strictly greater than \(\lambda\), and
\(\alpha=v(t)>2\lambda\). Write
\(\delta_t=\alpha-2\lambda\). The genuine pole cost (12) in
[the deepest-center criterion](SFC_DEEPEST_RATIONAL_ORDER_WINDOWS.md#7-a-sharper-selection-using-the-support-and-the-slope-costs)
at each nonnegative integer center \(r\) of a window of
width at least two is \(\mu(r)=\delta_t(r+1)\). The only
perturbation that lowers either row is the constant
monomial, so the first genuine pole uses \(r+1\) copies of
that monomial. Centers with a genuine pole are integers
because the two anchor exponents are one. Negative
window-base centers with a genuine pole occur in a finite
range depending on the width.

For \(F_*=t+xy\), let \(A_r^*,B_r^*\) be the leading and
regular Laurent coefficients \(B_0(r),B_1(r)\) of the
pure-normalized moments
\(\Lambda(F_*^n)/(n!)^2\). In other words the two terms in
the Laurent series are \(A_r^*/z+B_r^*\).
Write \(H_\ell=\sum_{i=1}^\ell 1/i\), \(H_0=0\).
The term of count \(k\) is \(t^k/(k! (n)_{\underline{k}})\).
Consequently
\[
 \begin{split}
     r!A_r^*&=Y_r,\\
     r!(B_r^*+H_r A_r^*)&=V_r+\nu_r,\qquad
     \nu_r=\sum_{\ell\ge0}
       \frac{t^{r+1}(-t)^\ell H_\ell}{\ell!(r+1+\ell)!}.             \tag{12}
 \end{split}
\]
Indeed for \(k\le r\) the regular value gives the summands
of \(V_r/r!\). For \(k=r+1+\ell\), stripping the single pole
gives the residue \(t^{r+1}(-t)^\ell/(\ell!(r+1+\ell)!r!)\).
The sum of reciprocals of the nonzero falling-factorial
offsets is \(H_r-H_\ell\), giving the displayed regular
term. The estimates of Section 1 and \(v(H_\ell)\ge
-\log_p(\ell+1)\) give
\[
 v(Y_r)\ge(\alpha-\lambda)(r+1)+\lambda,\qquad
 v(\nu_r)\ge(\alpha-\lambda)(r+1)-C_t .
\]
The two products containing \(\nu_r\) in the determinant of
the columns \(Y\) and \(V+\nu\) at \(r,r+1\) have valuations
at least \(2(\alpha-\lambda)(r+1)-C_t'\).
This exceeds \(\alpha(r+1)\) for large \(r\).
By Section 1 the determinant therefore has valuation
\(\alpha(r+1)\) for large \(r\). Undoing the harmonic shift
in (12) changes the determinant by
\((H_{r+1}-H_r)Y_rY_{r+1}\), whose valuation also exceeds
\(\alpha(r+1)\) for large \(r\).
It follows that the original jet determinant satisfies
\[
 v\!\begin{vmatrix}A_r^*&B_r^*\\ A_{r+1}^*&B_{r+1}^*\end{vmatrix}
   =\alpha(r+1)-v(r!)-v((r+1)!)
   =\delta_t(r+1)+O(1+\log_p(r+2)).                                 \tag{13}
\]

For \(Q=xy(1+P_0)\), Section 3 gives \(E_m=0\). The same
transport below applies when the approximation has error (11).
In that case the error bound in Section 2 improves by a
linear amount in the order: use \(v(E_m)\ge\lambda(m+1)\)
and \(v(t^k/k!)\ge(\alpha-\lambda)k\), then subtract
\(v(n!)\). Divide (9) by \(n!\) and \((n+1)!\) in its two
components to obtain a
matrix \(\mathcal T(n)\) on the pure-normalized reference
pair, where \(\mathcal T(n)=\operatorname{diag}(1,1/(n+1))\)
times the matrix in (9) times \(\operatorname{diag}(1,n+1)\).
The error after this division has component valuations
growing at least linearly in \(n\), by \(\alpha>2\lambda\).
At each fixed integer \(r\ge0\), the shifted tower expansion
therefore transports the leading and regular columns by
\[
                       A=\mathcal T(r)A^*,\qquad
               B=\mathcal T(r)B^*+\mathcal T'(r)A^* .                \tag{14}
\]
Here each column consists of its two components at \(r,r+1\).
This follows by expanding (9) along \(n=r+p^\ell\), where
the error is beyond every fixed jet precision.
The determinant of \(\mathcal T(r)\) is a unit by (9).
Both \(\mathcal T(r)\) and \(\mathcal T'(r)\) have entries
of valuation at least \(-O(1+\log_p(r+2))\), since the
only denominators introduced are powers of \(r+1\) and
the series and their order derivatives have integral
coefficients. The extra determinant in (14) therefore
has valuation at least \(2\delta_t(r+1)-O(1+\log_p(r+2))\).
This is larger than (13) for all sufficiently large \(r\).
Thus (13) also holds for the perturbed phase and supplies
the minor upper bound at large heights in this class.
The all-order bound (2) is already available without
taking any order limit.

## 6. The monomial pole cost can miss every low jet minor

Let \(u,v,t\) be nonzero with \(v(u),v(v)\ge2\lambda\), and put
\[
    F=(x+u)(y+v)+t,\qquad
    \delta_{\rm mon}=v(uv)-2\lambda,\quad
    \delta_t=v(t)-2\lambda .
\]
Take \(\delta_t>4\delta_{\rm mon}\).
The phase has support size four and is absolutely irreducible:
after translation it is \(xy+t\), primitive and linear in \(y\)
over the polynomial ring in \(x\). It has strictly isolated
effective slope at \(xy\). The constant coefficient has
valuation \(v(uv+t)=v(uv)\).

At a nonnegative integer base \(r\) the genuine monomial pole
cost is
\[
                             \mu(r)=\delta_{\rm mon}(r+1).
\]
In fact if the counts of the constant, \(x\), and \(y\) are
\(k_0,k_x,k_y\), then the row deficits are \(k_0+k_y\) and
\(k_0+k_x\), and their linear effective cost is
\[
       (k_0+k_x)\bigl(v(v)-\lambda\bigr)
        +(k_0+k_y)\bigl(v(u)-\lambda\bigr).
\]
Both deficits must exceed \(r+j\) for a genuine pole at
offset \(j\). Taking \(j=0\) and \(k_0=r+1\) realizes the
minimum.

On the other hand \(Q=(x+u)(y+v)\) satisfies (10) and (11),
so the transfer (14) holds with the scalar \(t\).
In particular the adjacent two-column determinant has
valuation \(\delta_t(r+1)+O(1+\log_p(r+2))\), by (13)–(14).
The leading coefficient at every offset has valuation at
least \(\delta_t(r+1)-O(1+\log_p(r+2))\), by the estimate for
\(Y_r\) and (14). Every fixed regularized jet coefficient
has valuation at least \(-O(1+\log_p(r+2))\), by
[the uniform coefficient estimate (8)](SFC_DEEPEST_RATIONAL_ORDER_WINDOWS.md#3-uniform-taylor-expansion-at-the-selected-center)
and the infimum over counts. Consequently for every
\(0\le J<4\), every minor with columns \(0,\ldots,J\) and
\(J+1\) of the four offsets has valuation at least
\[
                       \delta_t(r+1)-O_J(1+\log_p(r+2)).              \tag{15}
\]
In particular an upper bound of the form
\((J+1-\eta)\mu(r)+O_J(1+\log_p(\mu(r)+2))\) with any
\(\eta>0\) fails uniformly for every such choice of \(J\).

This obstruction occurs at actual large-depth selected centers.
For any fixed \(0<\epsilon<1\), set
\[
     L_r=\lfloor(1-\epsilon)\mu(r)\rfloor+1,\qquad N_r=r+p^{L_r}.
\]
The weighted selection in Section 7 of the deepest-center proof
has \(r\) eligible at depth \(L_r\). Every eligible center at
\(N_r\) has height \(O(1+\log_p(N_r+4))=O(1+L_r)\).
Any distinct one of depth at least \(L_r\) would have difference
from \(r\) of valuation at least \(L_r\) and numerator
\(O(1+L_r)\), which is impossible for large \(r\).
Thus the selected center is \(r\). At these starts (15) satisfies
all the necessary determinant inequalities of the deepest-center
proof despite the two-window nonvanishing theorem in Section 4.
For example \(p=2\), \(u=v=4\) and \(t=2^{12}\) give
\(\delta_{\rm mon}=2\) and \(\delta_t=10\).

Here the scalar critical value controls actual jet cancellation.
The cost in (12) of the deepest-center proof sees the cheaper
monomial paths before their cancellation. A quantitative
separation for general isolated-slope phases must account for
this difference. Two options are to use a grouping before
assigning pole costs or to compare higher jet layers after
lower layers are known to be too divisible to contribute
at the selected depth.
