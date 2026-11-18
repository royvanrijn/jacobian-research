# Real sign balance for reciprocal odd Laurent pencils

A real odd Laurent pencil with all its critical points real cannot use
arbitrarily unbalanced reciprocal sign factors. In degree \(4h+2\),
the plus and minus factors differ in their numbers of reciprocal pairs
by at most two. The sign and the possible difference are controlled by
one explicit ratio threshold.

In the normalized degree-26 problem this excludes the entire
\((2,10)\) mixed case from the all-real locus. The \((4,8)\) case can
occur there only as plus degree eight and minus degree four, with the
positive reciprocal square root less than \(11/13\).
These are real restrictions; the global rational systems without the
all-real hypothesis remain open.

## Setting

Use the monic normalization from the
[mixed coupled-ODE note](MIXED_RECIPROCAL_COUPLED_ODE.md):
\(N\) has degree \(2h+1\), \(N(1)=N'(1)=0\), and

\[
 2uN'-N=(u-1)S_N(u),\qquad f(z)=N(z^2)/z,\qquad L=4h+1.
\]

The residual roots are distinct, nonzero, different from 1, and
paired without fixed points by \(u\leftrightarrow\lambda/u\);
\(N(u)^2/u\) is nonzero at them and has exactly one distinct value
per pair. Thus the critical point at \(u=1\) has zero value, and all
remaining critical values are nonzero.

Assume every critical point of \(f\) is real. All the residual roots in
\(u\) are then positive. Thus \(\lambda>0\); choose the positive
\(\rho=\sqrt\lambda\) and set \(r=1/\rho>0\). Write

\[
 F(v)=\rho^{-(2h+1)}N(\rho v),\quad
 \Sigma(v)=\frac{S_N(\rho v)}{L\rho^{2h}}
          =v^hP(w)Q(w),\quad w=v+v^{-1}.
\]

Here \(P,Q\) are real monic polynomials, of degrees \(p,q\), with
\(p+q=h\). Every root of \(P Q\) is real and exceeds 2. The plus sign
means \(F(v)+vF(1/v)=0\); the minus sign means
\(F(v)-vF(1/v)=0\). The fixed-point-free good locus excludes \(w=2\).

With \(C_j,D_j\) defined by the usual recurrence

\[
 C_0=D_0=1,\quad C_1=w-1,\quad D_1=w+1,\quad
 Z_{j+1}=wZ_j-Z_{j-1},
\]

the monic polynomials

\[
 X=\frac{F(v)+vF(1/v)}{v+1},\qquad
 Y=\frac{F(v)-vF(1/v)}{v-1}
\]

have degree \(2h\), satisfy \(P^2\mid X\), \(Q^2\mid Y\), and obey

\[
 \begin{split}
 X+2(w+2)X'&=L\alpha_h(w)P Q,\\
 Y+2(w-2)Y'&=L\beta_h(w)P Q,\\
 \alpha_h&=D_h-rD_{h-1},\qquad
 \beta_h=C_h-rC_{h-1}.                           \tag{1}
 \end{split}
\]

These are the exact identities already proved in the mixed note.

## Lemma 1. Exact outside-root thresholds

On \(w>2\), the polynomial \(\beta_h\) has one root if \(r>1\)
and none if \(r\le1\). The polynomial \(\alpha_h\) has one root
if

\[
                      r>t_h:=\frac{2h+1}{2h-1},
\]

and none if \(r\le t_h\). Each such outside root is simple.
At equality the corresponding root is at 2, outside the open interval
being counted.

Indeed, \(C_j,D_j\) are positive on \(w>2\), with
\(C_j(2)=1\) and \(D_j(2)=2j+1\). Their consecutive ratios obey

\[
 \frac{Z_j}{Z_{j-1}}
   =w-\frac1{Z_{j-1}/Z_{j-2}}\qquad(j\ge2).
\]

Starting from \(C_1/C_0=w-1\), \(D_1/D_0=w+1\), induction shows
both ratios are strictly increasing on \(w>2\) and tend to infinity.
Their initial values are 1 and \(t_h\). This proves the assertions,
including simplicity and the endpoint cases.

## Theorem 2. Uniform necessary real balance

Put

\[
 \epsilon_\beta=\mathbf1_{\{r>1\}},\qquad
 \epsilon_\alpha=\mathbf1_{\{r>t_h\}}.
\]

Then

\[
 \boxed{\quad q\le p+\epsilon_\beta,\qquad
                 p\le q+1+\epsilon_\alpha.\quad}               \tag{2}
\]

For the first inequality, apply Rolle's theorem to
\(Z=\sqrt{w-2}\,Y\) on \(w\ge2\). It has a zero at the boundary
2 and at least double zeros at the \(q\) roots of \(Q\). Between 2
and the first such root, and between each consecutive pair, there
must be a derivative zero. Equation (1) gives

\[
                  Z'=\frac{L\beta_h P Q}{2\sqrt{w-2}}.
\]

None of those open intervals contains a root of \(Q\). There are at
most \(p+\epsilon_\beta\) eligible distinct roots of \(P\beta_h\),
proving the first inequality. Multiplicities or coincidences with the
outside root can only reduce that number.

Similarly, between the \(p\) roots of \(P\) there are at least \(p-1\)
derivative zeros of \(\sqrt{w+2}\,X\). Its derivative is
\(L\alpha_h P Q/(2\sqrt{w+2})\), giving
\(p-1\le q+\epsilon_\alpha\). No lower bound on \(p,q\) was used;
uniform sign cases are covered.

Since \(p-q\) has the parity of \(h\), (2) gives the following complete
list of sign counts that this necessary condition permits:

| Range of \(r\) | \(h\) even | \(h\) odd |
|---|---|---|
| \(0<r\le1\) | \(p=q=h/2\) | \(p=q+1\) |
| \(1<r\le t_h\) | \(p=q=h/2\) | \(p-q\in\{-1,1\}\) |
| \(r>t_h\) | \(p-q\in\{0,2\}\) | \(p-q\in\{-1,1\}\) |

The table is a necessary condition, not an existence assertion.

## Corollary 3. All-real degree twenty-six

For \(h=6\), \(t_h=13/11\). Therefore

\[
 (p,q)=(3,3),\qquad\text{or}\qquad
 (p,q)=(4,2)\ \text{with }r>13/11.               \tag{3}
\]

The raw sign-factor degrees in \(v\) are \(2p,2q\). Thus:

* Both orientations of the \((2,10)\) split are excluded over
  \(\mathbb R\) when all critical points are real.
* The unbalanced allowed split has plus degree eight and minus
  degree four. It requires \(0<\rho<11/13\), equivalently
  \(0<\lambda<121/169\).
* For \(\rho\ge11/13\), only the balanced \((6,6)\) split remains.

For a rational model in the unbalanced case, \(\lambda\) must moreover
be a rational square: quadratic conjugation for a nonsquare parameter
interchanges the sign factors and would require their degrees to agree.
All these inequalities refer to the fixed normalization \(N(1)=N'(1)=0\)
and to the positive square root \(\rho\). Changing its sign interchanges
the labels and the half-line on which the \(w\)-roots lie.

The known
[rational degree-ten model](RATIONAL_DEGREE_TEN_ALTERNATING_PENCIL.md)
shows that the difference two in Theorem 2 can occur: it has \(h=2\),
\(r=15/2>5/3\), and \((p,q)=(2,0)\), with all critical points real.
Consequently the present argument cannot be strengthened to universal
equality \(p=q\).

The companion
[degree-twelve note](RECIPROCAL_DEGREE_TWELVE_OBSTRUCTION.md)
proves the even critical-degree version. In degrees \(4k\), where the
additional factors \(\alpha_h,\beta_h\) are absent, the same boundary
argument forces exact balance for even \(k\), and a difference one for
odd \(k\).

## Replay and limits

~~~text
.venv/bin/python research/arithmetic/verify_reciprocal_real_sign_balance.py --write-receipt research/arithmetic/reciprocal_real_sign_balance_verification.json
~~~

The checker verifies exact outside-root counts, including both threshold
endpoints, in retained finite controls; polynomial positivity of the
ratio derivatives; the permitted integer sign-count table; and the
sharp all-real degree-ten model. The written induction and Rolle
argument provide the uniform proof. There is no claim that the retained
finite controls establish all degrees, or that each permitted count is
realizable.

The general rational degree-26 coupled systems remain unchanged.
Neither rational nonexistence nor a rational construction follows for
models with nonreal critical points. No higher-degree real reciprocal
existence theorem is asserted.
