# Two-window nonvanishing with truncated factorial weights

Let \(K\) be a complete nonarchimedean field of characteristic
zero with a real valuation extending \(v_p(p)=1\).
Put \(\lambda=1/(p-1)\), and for an integer \(e\ge0\) write
\[
           L_e(j)=\sum_{s=1}^e\left\lfloor\frac{j}{p^s}\right\rfloor,
           \qquad \lambda_e=\lambda(1-p^{-e}) .
\]
For a two-variable polynomial \(F=c(xy+P)\), \(c\ne0\), suppose
\[
 \begin{gathered}
        v([x^iy^j]P)+L_e(i)+L_e(j)\ \ge D\quad\text{for every }i,j,\\
                        D>2e+\lambda p^{-e},\qquad D\ge e+2\lambda .
                                                                    \tag{1}
 \end{gathered}
\]
Then, for every \(n\ge0\),
\[
         \min_{j=0,1}v\left(\frac{\Lambda(F^{n+j})}
                                    {c^{n+j}(n+j)!}\right)
                            \le \lambda n.                          \tag{2}
\]
Thus every two consecutive orders have a nonzero ordinary
factorial moment, with no degree or support bound. It suffices
that the coefficient field carry the valuation, since completion
preserves every ordinary moment. The case \(e=0\) recovers
[the unweighted coefficient neighborhood of \(xy\)](SFC_MIXED_ANCHOR_CONSTANT_SHIFT.md).
When \(p=2,e=1,D=3\), condition (1) allows a term
\(2^{3-m}x^{2m}\) for every \(m\ge1\), so the valuations of
individual admissible coefficients can be arbitrarily negative.

A sufficient bound using only effective slopes is obtained by
taking \(p^e\) greater than the degree of \(F\).
Then \(L_e(i)+L_e(j)=v(i!j!)\) on its support, and
\[
       v([x^iy^j]P)+L_e(i)+L_e(j)
        =v([x^iy^j]P)+\lambda(i+j-2)
                            +2\lambda-\lambda(s_p(i)+s_p(j)).
\]
Since \(s_p(i)+s_p(j)\le 2e(p-1)\), an effective slope gap
of at least \(4e\) suffices in (1). This bound grows only
logarithmically with the maximum degree for fixed \(p\).

The proof uses the same two-term reference sequence as the
unweighted theorem. Once the normalized moments of an auxiliary
phase are analytic on order disks \(a+p^e\mathbb Z_p\), their
rising factorial expansion can lose \(L_e(h)\) in coefficient
number \(h\). The finite tridiagonal recurrence recovers that
loss, up to a weight on one entry of the adjacent pair.
Truncated factorial weights on the spatial coefficients supply
such order functions by factoring at a stationary point. A
locally constant unit-valued counterexample in Section 5 shows
that a transfer for arbitrary order functions genuinely needs
a relation between the order disks and the scalar.

## 1. Interpolation on order disks

For a function \(A\) on \(\mathbb Z_p\) represented by restricted
series \(A(a+p^e Z)\) for \(0\le a<p^e\), let \(\nu_e(A)\)
be the infimum of the valuations of their coefficients.
The reflected rising factorial expansion is
\[
                  A(z)=\sum_{h\ge0}d_h (z+1)\cdots(z+h),
                  \qquad
 \begin{cases}
                  v(d_h)+L_e(h)\ge\nu_e(A),\\
                  v(d_h)+L_e(h)\longrightarrow+\infty.
 \end{cases}                                                       \tag{3}
\]
The empty product is one. This is the local Newton coefficient
bound of [Amice](https://www.numdam.org/articles/10.24033/bsmf.1606/),
in the reflected variable \(-z-1\). Here is a proof that
also applies over the coefficient fields used below.

Take \(d_h\) to be the divided difference of \(A\) at the
nodes \(-1,\ldots,-h-1\). In its sum over the nodes, group
them into residue classes modulo \(p^e\). For a class \(a\)
of size \(m_a>0\), write \(Q_a(z)\) for the product of
\(z-z_i\) over nodes outside that class. The contribution
of the class is the divided difference of \(A/Q_a\) at
its own nodes. On its order disk each outside factor,
after removing its constant value at \(a\), has form
\(1+bZ\) with \(v(b)\ge1\). Rescaling its inside nodes
to the disk coordinate \(Z\) costs \(e(m_a-1)\) in the
divided difference. The additional loss from \(Q_a\)
is \(\sum_{i\text{ outside}}v(a-z_i)\). Their sum equals
\[
       \sum_{s=1}^e(\#\{i:z_i\equiv a\bmod p^s\}-1)\le L_e(h),
\]
using the consecutive nodes. On the unit disk the divided
difference of a series at integral nodes costs no coefficient
valuation. For the monomial of degree \(b\) the divided
difference at \(m\) nodes is zero for \(b<m-1\) and a
complete homogeneous symmetric polynomial of degree
\(b-m+1\) in the nodes otherwise, as follows by the divided
difference recursion. Thus it uses only series coefficients
from degree \(m-1\) onward.

For convergence in (3), after removing the constant of
\(Q_a\), every inverse factor has coefficient valuation
at least its degree, and so does the entire inverse product.
Approximate \(A\) on each disk by a polynomial of fixed
degree with a tail of arbitrarily high coefficient valuation.
In the product with the inverse, all coefficients of degree
at least \(m_a-1\) then have valuations growing uniformly
as \(m_a\to\infty\). This proves the second assertion in
(3) even though the outside products change with \(h\).
The new series converges uniformly on \(\mathbb Z_p\),
where each product of \(h\) consecutive factors has
valuation at least \(L_e(h)\). It equals \(A\) at the
negative integer nodes by Newton interpolation and thus
on \(\mathbb Z_p\) by continuity.

## 2. Constant shift with only local order analyticity

Suppose a phase \(Q\) in any number of variables has an
approximation
\[
       \frac{\Lambda(Q^m)}{(m!)^2}=A(m)+E_m,\qquad
       \nu_e(A-1)>0,\qquad v(E_m)>\lambda s_p(m),                     \tag{4}
\]
where \(A\) is as in Section 1 and \(s_p(m)\) is the base-\(p\)
digit sum. For \(v(t)\ge e+2\lambda\), the phase \(t+Q\)
satisfies (2) with \(c=1\). The order radius \(p^{-e}\) and
the degree and dimension of \(Q\) are arbitrary.

For \(t=0\), the normalized values in (4) are units so the
bound follows from \(v(m!)\le\lambda m\). Suppose \(t\ne0\),
write \(\alpha=v(t)\), and use
\[
  V_n=\sum_{k=0}^n\frac{t^k}{k!}(n-k)!,\qquad
  b_n=\frac{t^n}{n!},\qquad
  T_{n,h}=\sum_{k=0}^n\frac{t^k}{k!}(n+h-k)! .
\]
As in [the unweighted transfer](SFC_MIXED_ANCHOR_CONSTANT_SHIFT.md#2-adding-a-constant-to-analytically-normalized-moments),
\(T_{n,0}=V_n\), \(T_{n,1}=V_{n+1}-b_{n+1}\), and
\[
          T_{n,h}=(n+h)T_{n,h-1}-tT_{n,h-2}+t b_n(h-2)!
                                               \quad (h\ge2).       \tag{5}
\]
Write \(T_{n,h}=a_h(n)V_n+v_h(n)(V_{n+1}-b_{n+1})+c_h(n)b_n\).
Set \(d(n)=\min(v_p(n+1),e)\). The coefficient bounds are
\[
 \begin{split}
       v(a_h(n))&\ge L_e(h),\\
       v(v_h(n))&\ge L_e(h)-d(n),\\
       v(c_h(n))&\ge L_e(h)-e.                                     \tag{6}
 \end{split}
\]
To prove these, use the polynomial given by the recurrence
\(K_j=(n+j)K_{j-1}-tK_{j-2}\) on a consecutive list of slots.
Its terms correspond to covering each slot by one singleton
with its linear factor or by a pair of adjacent slots with
weight \(-t\). Every term has valuation at least the sum
of the valuations of the full list capped at \(e\).
Indeed two adjacent integers have at most one nonunit, so
their combined capped valuation is at most \(e\le\alpha\).
A list of length \(h\) has capped sum at least \(L_e(h)\).
For \(v_h(n)\) the slots are \(n+2,\ldots,n+h\), so the
missing initial slot can reduce the bound by at most \(d(n)\).
For \(a_h(n)\) with \(h\ge2\) the slots start at \(n+3\)
and carry a prefactor \(-t\). The two missing slots
together cost at most \(e\). The inhomogeneous coefficient
\(c_h(n)\) is a sum over \(2\le i\le h\) of \(t(i-2)!\)
times the recurrence polynomial for the slots
\(n+i+1,\ldots,n+h\). One term therefore has valuation at
least \(\alpha+L_e(i-2)+L_e(h-i)\).
For any nonnegative lengths \(r,s\),
\(L_e(r)+L_e(s)\ge L_e(r+s)-e\). Also
\(L_e(h)-L_e(h-2)\le e\), since the two consecutive
steps in its floor sums concern at most one multiple of
\(p\). These estimates prove the bound for \(c_h\), and
the initial \(h=0,1\) follow from (5).

Use (3) on \(A-1\). Summing (5) with the expansion
coefficients \(d_h\) of \(A\) is legitimate by (3) and (6).
Let \(u(n),v_0(n),c_0(n)\) be the sums of the three
coefficient columns of (5), respectively. The model
moments \(\mathcal W_n\) obtained with \(E_m=0\) satisfy
\[
 \binom{\mathcal W_n}{\mathcal W_{n+1}}
   =
   \begin{pmatrix}
      u(n)&v_0(n)\\
      -t v_0(n+1)&u(n+1)+(n+2)v_0(n+1)
   \end{pmatrix}\binom{V_n}{V_{n+1}}
   +
   \binom{-v_0(n)b_{n+1}+c_0(n)b_n}{c_0(n+1)b_{n+1}}.                \tag{7}
\]
If \(\mu=\nu_e(A-1)>0\), then
\(v(u(n)-1)\ge\mu\), \(v(v_0(n))\ge\mu-d(n)\)
and \(v(c_0(n))\ge\mu-e\).
Scale both pairs by \(\operatorname{diag}(1,p^{-d(n)})\).
The conjugated matrix is the identity modulo the maximal ideal.
For the second diagonal entry use \(v_p(n+2)\ge d(n+1)\).
For the lower off-diagonal entry use
\(\alpha\ge e\ge d(n)+d(n+1)\). Hence the matrix preserves
the minimum valuation of the scaled pair.

The reference pair has the sharper bound in these coordinates
\[
              \min(v(V_n),v(V_{n+1})-d(n))\le \lambda n-d(n).        \tag{8}
\]
Indeed the reference determinant of
[the finite Bessel truncations](SFC_MIXED_ANCHOR_CONSTANT_SHIFT.md#1-a-quantitative-reference-pair-including-the-boundary)
has valuation \(\alpha(n+1)\); the two companion entries
after inverse scaling have valuations at least
\((\alpha-\lambda)(n+1)+\lambda+d(n)\).
Here \(\alpha-\lambda\ge d(n)\), so the increase in the index
of the first companion entry covers the scaling of the second.
This gives (8).

For \(n\ge1\) both entries of the error in (7) after scaling
have valuations strictly greater than \(\lambda n\):
use \(v(b_s)\ge(e+\lambda)s\),
\(v(v_0(n))\ge\mu-d(n)\), and \(v(c_0(n))\ge\mu-e\).
The additional binomial sum from \(E_m\) at order \(s\)
has valuation \(>\lambda s\), since \(v(m!E_m)>\lambda m\)
and \(v(b_k)\ge\lambda k\). After scaling the pair
this last error has valuations \(>\lambda n-d(n)\).
The scaled model and actual pairs therefore have the
same minimum valuation as the scaled reference pair.
Removing the scaling raises their minimum by at most
\(d(n)\), giving (2). At \(n=0\) the first value in (4)
is a unit, which suffices.

There is also a large-order version when \(\alpha\ge e\) and
\(\alpha>2\lambda\). In (8) replace the right side by
\(\lambda n-g(n)\), where \(g(n)=\min(d(n),\alpha-\lambda)\);
the same determinant comparison gives this bound.
The conjugated matrix is still the identity modulo the maximal
ideal, and when \((\alpha-2\lambda)n\ge2e\) the two errors
in (7) have valuations greater than \(\lambda n\) after
scaling, since their coefficient losses are at most \(2e\).
The errors from \(E_m\) are above the reference bound because
\(d(n)-g(n)\le\lambda\) when \(\alpha\ge e\).
Consequently at every such positive start the two unscaled
moments have minimum valuation at most
\(\lambda n+d(n)-g(n)\).

## 3. Count estimates for weighted analytic factors

For a series \(B(X)=\sum_\nu b_\nu X^\nu\), use the spatial
valuation
\[
                 \mu_e(B)=\inf_\nu \bigl(v(b_\nu)+
                                         \sum_i L_e(\nu_i)\bigr).
\]
We use series whose weighted coefficient valuations tend to
infinity. They form a complete algebra in this valuation:
superadditivity of \(L_e\) proves the product bound for
polynomials and the algebra is their completion. Factorial
evaluation extends continuously since
\(v(\prod_i\nu_i!)\ge\sum_iL_e(\nu_i)\).

Suppose
\[
 Q=(1+U)(X_1+B_1)(X_2+B_2),\qquad
          \mu_e(U)>\lambda p^{-e},\qquad
          \mu_e(B_i)\ge2\lambda\quad(i=1,2).                        \tag{9}
\]
Then (4) holds, even when many coefficients of the factors
have negative valuations. Expand the three powers separately
at order \(m\), and let \(K_0,K_1,K_2\) count their
perturbations. Write \(\gamma_i\) for the total exponent
in row \(i\) from all perturbation monomials. At row 1 the
anchor falling factorial exactly combines with the row
factorial ratio to give the polynomial
\[
              (z)_{\underline{K_1}}\,
              \frac{(z-K_1+\gamma_1)!}{z!}
                =\prod_{j=1}^{\gamma_1}(z-K_1+j),                    \tag{10}
\]
where factorial ratios in \(z\) mean finite rational products.
Row 2 has the same formula with \(K_2,\gamma_2\).
The falling factorial of length \(K_0\) remains, and
other rows contribute the constants \(\gamma_i!\).
On an order disk \(z=a+p^e Z\), each product of \(\ell\)
consecutive integer offsets has coefficient valuation
at least \(L_e(\ell)\); the factors count one contribution
for each multiple up to each level \(p^s\), \(1\le s\le e\).

Thus the positive products at rows \(1,2\) and the spectator
rows recover the spatial weights of all the selected
coefficients by superadditivity. For count factorials of
the unit factor use \(\sum_\nu v(k_{0,\nu}!)\le v(K_0!)\).
The remaining falling factorial for that factor recovers
\(L_e(K_0)\), leaving a cost at most
\(v(\lfloor K_0/p^e\rfloor!)\le\lambda p^{-e}K_0\).
For the other two count factorial costs use \(\lambda K_1\)
and \(\lambda K_2\). The resulting coefficient valuation
bound on every order disk is
\[
       (\mu_e(U)-\lambda p^{-e})K_0+
       (\mu_e(B_1)-\lambda)K_1+
       (\mu_e(B_2)-\lambda)K_2 .                                   \tag{11}
\]
A zero perturbation factor can be omitted. All active
coefficients of the \(K_i\) are strictly positive.
For bounded counts only finitely many perturbation choices
can have valuation bounded above at a given precision,
by the weighted coefficient tails. Thus the count sum over
all \(K_i\) defines an order function \(A\) with
\(\nu_e(A-1)>0\). At an ordinary order \(m\), omitted terms
with \(K_0>m\) vanish because their falling factorial from
the unit factor remains. Other omitted terms have either
\(K_1>m\) or \(K_2>m\), giving \(v(E_m)\ge\lambda(m+1)\)
by (9)–(11). This proves (4).

In particular, for \(Q=X_1X_2(1+U)\) one has \(E_m=0\).
If \(U\) is a polynomial and every nonzero coefficient
\(u_\nu\) satisfies \(v(u_\nu)+\sum_i v(\nu_i!)>0\),
then some finite \(e\) gives \(\mu_e(U)>\lambda p^{-e}\).
Thus adding a scalar \(t\) with \(v(t)\ge e+2\lambda\)
gives a two-window theorem in that whole set, in arbitrary
dimension and degree.

## 4. Factoring after the weighted stationary-point shift

Normalize \(c=1\) in (1). A coordinate derivative of any
perturbation monomial evaluated at \(v(x),v(y)\ge\lambda_e\)
has valuation at least \(D-\lambda_e\), since
\(L_e(j)\le\lambda_e j\). The difference of its values at two
such points is at least the minimum valuation of their
coordinate differences plus \(D-2\lambda_e\), by factoring
one coordinate difference. Thus
\((x,y)\mapsto(-\partial_yP(x,y),-\partial_xP(x,y))\)
is a contraction of this closed ball into itself.
Its fixed point \(h\) has coordinate valuations
at least \(D-\lambda_e\). The critical value
\(t=h_xh_y+P(h)\) has valuation at least \(D\).

Translation by any point of coordinate valuations at least
\(\lambda_e\) preserves the spatial valuation bounds.
For a single monomial use
\[
  v_p\binom{a}{b}+L_e(b)+L_e(a-b)\ge L_e(a),
\]
which follows by separating the floor sums above level \(e\)
in the factorial quotient. Multiplication by a translation
coordinate to the power \(a-b\) supplies at least
\(L_e(a-b)\). Apply the bound in each variable and pass
to the completion. Thus \(F(X+h)-t=xy+E\) with
\(\mu_e(E)\ge D\), and \(E\) has no terms of degree below two.

Decompose the space of such perturbations into the mixed
terms \(A\) divisible by \(xy\), the pure \(y\) terms \(B\)
divisible by \(y^2\), and the pure \(x\) terms \(C\) divisible
by \(x^2\). The projections preserve valuation.
For a triple with all valuations at least \(D\), write
\[
     W=1+A/(xy),\qquad S=x+B/y,\qquad T=y+C/x .
\]
Every division by \(x\) or \(y\) loses at most \(e\) because
\(L_e(j)-L_e(j-1)=\min(v_p(j),e)\).
The product \(WST\) has linear part \(xy+A+B+C\) in the
three unknowns. Its remaining terms are
\[
      (B/y)(C/x)+(B+C)A/(xy)+(B/y)(C/x)A/(xy).
\]
Their valuations are at least \(2D-2e,2D-2e,3D-4e\),
respectively. The difference at two input triples increases
the minimum valuation of their differences by at least
\(D-2e>0\). Solve \(WST=xy+E\) by repeatedly projecting
\(E\) minus these remaining terms into its three parts.
This maps the closed triple ball to itself and is a
contraction there. It gives
\(\mu_e(W-1)\ge D-2e\) and
\(\mu_e(S-x),\mu_e(T-y)\ge D-e\).
Translation back by \(-h\) preserves these bounds, including
the terms from the translated \(x,y\) anchors because
\(v(h_i)\ge D-\lambda_e\ge D-e\).
Consequently \(F-t\) satisfies (9).
Finally \(v(t)\ge D\ge e+2\lambda\), so Section 2 proves (2).

The estimates give one way to exploit a positive effective
slope without imposing an unweighted Gauss bound. They still
require the explicit truncated weights in (1) or the local
conditions (4) and (9). A unique least effective slope by
itself does not supply the positive disk norm or the
factorization bounds in this argument.

## 5. The disk size matters for unrestricted order functions

A positive local disk valuation in (4), without a restriction
relating the disk size to \(t\), does not suffice for the scalar
transfer. There is a unit-valued locally constant model with
\(E_m=0\), \(p=2\), \(t=4\) and two consecutive exact zeros.
This example concerns the abstract order functions in Section 2.

Take the reference \(V_n\) in Section 2 at \(t=4\). Define
\[
     x=\frac{2V_{38}-V_{39}}{2\cdot38!},\qquad
     y=\frac{V_{39}-4V_{38}}{8\cdot37!}.
\]
Let \(A\) be \(1+x\) on \(38+64\mathbb Z_2\), \(1+y\) on
\(37+64\mathbb Z_2\), and \(1\) elsewhere. Then
\(v_2(x)=1\), \(v_2(y)=2\), and \(\nu_6(A-1)=1\).
For a short exact verification, put \(M_j=j!V_j\).
The recurrence is
\[
 M_0=1,\quad M_1=5,\qquad
 M_j=j^2 M_{j-1}-4j(j-1)M_{j-2}+4^j\quad(j\ge2).
\]
Modulo \(2^{75}\) this gives
\(M_{38}\equiv30\cdot2^{70}\) and
\(M_{39}\equiv24\cdot2^{70}\).
The two numerators after multiplying by \(39!\)
are \(78M_{38}-M_{39}\) and \(M_{39}-156M_{38}\),
of valuations \(72\) and \(74\), respectively.
The corresponding complete denominators have valuations
\(71\) and \(72\), giving the claimed values.

In the model binomial sums at orders \(38,39\),
the only changed entries are the \(m=38,37\) terms.
At order \(38\) their linear changes are \(38!x\) and
\(4\cdot37!y\); at order \(39\) the respective changes are
four and two times those entries. The definitions of
\(x,y\) make both resulting sums zero.
Here the scalar has the critical valuation \(2/(p-1)\),
and the analytic disks have order depth six. Extra structure
on actual factorial moment sequences would be needed for
a transfer that allows such independent choices of disk values.
