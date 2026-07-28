# Separable GVC escape obstructions in two variables

## 1. Results and scope

Work over a characteristic-zero field.  Write a nonzero
constant-coefficient operator in two variables as
\[
 \Lambda=\Lambda_r+\Lambda_{r+1}+\cdots+\Lambda_d,
 \tag{1.1}
\]
where every \(\Lambda_j\) is homogeneous of order \(j\) and
\(\Lambda_r\ne0\).  Thus \(r\) is the lowest positive order of
\(\Lambda\), not its usual order \(d\).

> **Theorem 1.1 — lowest-order binary GVC.**  Suppose that \(\Lambda\)
> has no order-zero term and that \(\deg P\leq r\).  If
> \[
>  \Lambda^m(P^m)=0\qquad(m\geq1),
>  \tag{1.2}
> \]
> then, for every fixed polynomial \(Q\),
> \[
>  \Lambda^m(QP^m)=0\qquad(m\gg0).
> \tag{1.3}
> \]
> If \(\Lambda\) has a nonzero order-zero term, (1.2) already forces
> \(P=0\), so the same conclusion is trivial.

This closes a natural nonhomogeneous-operator escape class.  A remaining
nonhomogeneous binary GVC counterexample must satisfy
\[
 \deg P>\min\{j:\Lambda_j\ne0\}.
 \tag{1.4}
\]
The translated-support strengthening of the
[split-symbol theorem](SPLIT_SYMBOL_GVC_THEOREM.md) removes every degree
restriction when \(\Lambda\) is homogeneous.  It also gives:

> **Theorem 1.2 — factor-unit nonhomogeneous class.**  Suppose
> \[
>  \Lambda=\Lambda_0\Gamma,
> \tag{1.5}
> \]
> where \(\Lambda_0\) is homogeneous with split symbol and \(\Gamma\)
> has a nonzero order-zero term.  Then \(\Lambda\) satisfies GVC for
> arbitrary \(P\).

For binary operators, every homogeneous \(\Lambda_0\) splits.  Thus a
remaining counterexample must also avoid a common homogeneous factor of
all operator pieces whose quotient has nonzero constant term.

More strongly, Theorem 3.4 below closes every binary operator with
nonzero linear part, for arbitrary \(P\).  Formal Weierstrass division
straightens its symbol to a separated drift through a locally finite
differential automorphism.  Consequently the surviving nonhomogeneous
frontier has lowest positive order at least two.  Theorem 3.7 further
shows that its first quadratic-leading polynomial degree is at least four.

There is also a conversion obstruction independent of differential
orders.  Let
\[
 F\in\operatorname{Sym}^4 U\otimes\operatorname{Sym}^4(U^*)
 \tag{1.6}
\]
be the bidegree-\((4,4)\) two-pair witness from
[`TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md`](TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md).

> **Theorem 1.3 — separated multiplicative conversion no-go.**  Start
> with one rank-one auxiliary datum
> \[
>  H=A(\zeta,u)P(z,u)
> \tag{1.7}
> \]
> over any commutative auxiliary algebra.  Apply polynomial substitutions,
> dilations, restrictions, or specializations that are algebra
> homomorphisms, preserve the dual/coordinate separation, and preserve the
> two gradings.  Then the bidegree-\((4,4)\) output has coefficient-matrix
> rank at most one.  A sum of \(s\) such channels has rank at most \(s\).
> Consequently no such one-channel conversion produces \(F\), and every
> finite-channel version producing \(F\) has width at least five.

The theorem includes auxiliary-variable evaluation after polarization,
rank-one dilation followed by separated restriction, and nonlinear
substitution on the two sides of the Segre product.  Coefficient
extraction can evade the rank-one conclusion only because it is not
multiplicative; it therefore does not formally transport all powers.
These statements do not classify general SIC rank strata.

## 2. Isolating the lowest operator order

If \(\Lambda\) has a nonzero constant term \(c\), the top homogeneous
part of \(\Lambda^m(P^m)\) is \(c^mP_e^m\), where \(P_e\) is the top
homogeneous part of \(P\).  Hence (1.2) forces \(P=0\).  Assume from now
on that \(r\geq1\).

If \(\deg P<r\), then
\[
 \operatorname{ord}(\Lambda^m)\geq rm>\deg(QP^m)
 \tag{2.1}
\]
for all large \(m\), which proves (1.3).  It remains to take
\(\deg P=r\).

Every term of \(\Lambda^m\) other than \(\Lambda_r^m\) has order greater
than \(rm=\deg(P^m)\).  Therefore the pure premise gives the exact
identity
\[
 0=\Lambda^m(P^m)=\Lambda_r^m(P_r^m),
 \tag{2.2}
\]
where \(P_r\) is the degree-\(r\) part of \(P\).

Over an algebraic closure the binary symbol of \(\Lambda_r\) splits:
\[
 \Lambda_r=D_{v_1}\cdots D_{v_r}.
 \tag{2.3}
\]
Put \(Vt=\sum_i t_iv_i\) and
\[
 H(t)=\frac{P_r(Vt)}{t_1\cdots t_r}.
 \tag{2.4}
\]
Complete polarization turns (2.2) into
\[
 \operatorname{CT}(H^m)=0\qquad(m\geq1).
 \tag{2.5}
\]
If \(P_r(Vt)\ne0\), the Duistermaat--van der Kallen theorem and strict
rational separation give an integral weight \(w\) such that
\[
 w\mathbin{\cdot}\beta\geq1
 \quad\text{for every }\beta\in\operatorname{Supp}(H).
 \tag{2.6}
\]
Equivalently, if \(W_0=\sum_iw_i\), every monomial \(t^\gamma\) in
\(P_r(Vt)\) has weight at least \(W_0+1\).

## 3. Bounded higher-order insertions

Expand the mixed operator:
\[
 \Lambda^m
 =
 \sum_{\substack{n_r+\cdots+n_d=m}}
 \binom{m}{n_r,\ldots,n_d}
 \prod_{j=r}^d\Lambda_j^{n_j}.
 \tag{3.1}
\]
Set
\[
 s=\sum_{j>r}n_j,\qquad
 \delta=\sum_{j>r}(j-r)n_j.
 \tag{3.2}
\]
The corresponding summand has order \(rm+\delta\).  Since
\(\deg(QP^m)\leq rm+\deg Q\), it vanishes unless
\[
 \delta\leq\deg Q.
 \tag{3.3}
\]
In particular \(s\leq\deg Q\).  Thus only finitely many patterns of
higher homogeneous pieces occur, independently of \(m\).  Each surviving
term has the form
\[
 \Lambda_r^{m-s}B(QP^m),
 \tag{3.4}
\]
where \(s\) and the constant-coefficient operator \(B\) belong to a fixed
finite set.

It remains to record the bounded-defect version of the split-symbol
argument.  Apply the derivatives in the fixed operator \(B\) by the
Leibniz rule.  In every resulting term, only a bounded number of the
\(m\) copies of \(P\) are differentiated.  In the translation
\(z\mapsto z+Vt\), a contribution to
\[
 [t_1^{m-s}\cdots t_r^{m-s}]
 \tag{3.5}
\]
can choose a lower-\(t\)-degree term from only a bounded number of the
remaining copies of \(P\): the target total degree differs from \(rm\)
by the fixed amount \(rs\), while \(B\) and \(Q\) have fixed degrees.
Consequently at least \(m-C\) factors contribute monomials of
\(P_r(Vt)\), for one constant \(C\) independent of \(m\).

The separator (2.6) gives a lower bound
\[
 m(W_0+1)-C'
 \tag{3.6}
\]
for the weight of such a contribution.  The target exponent in (3.5),
together with the bounded translated factors arising from \(B\) and
\(Q\), has weight at most
\[
 mW_0+C''.
 \tag{3.7}
\]
Equations (3.6)--(3.7) are incompatible for large \(m\).  Hence every
term (3.4), and therefore (1.3), vanishes eventually.

If \(P_r(Vt)=0\), every monomial of \(P\) has degree at most \(r-1\) in
the span of the directions \(v_i\).  The fixed operator \(B\) does not
increase this degree, while \(\Lambda_r^{m-s}\) differentiates
\(r(m-s)\) times in that span.  The direct degree inequality again kills
(3.4) for large \(m\).  This completes the proof of Theorem 1.1.
Theorem 1.2 follows from the injectivity argument in Section 5 of the
split-symbol theorem: a constant-coefficient operator with nonzero
constant term preserves the top homogeneous part up to a nonzero scalar
and is therefore injective.

The same filtration argument gives a stronger obstruction without the
hypothesis \(\deg P\leq r\).  Let \(e=\deg P\), \(q=\deg Q\), and write
\([R]_j\) for the degree-\(j\) homogeneous part of a polynomial \(R\).

> **Theorem 3.1 — fixed top-layer vanishing.**  Let \(\Lambda\) be an
> arbitrary nonhomogeneous constant-coefficient operator in two variables
> with lowest positive part \(\Lambda_r\), and suppose
> \(\Lambda^m(P^m)=0\) for every \(m\geq1\).  For every fixed \(Q\) and
> every fixed \(K\geq0\),
> \[
>  [\Lambda^m(QP^m)]_{(e-r)m+q-k}=0
>  \qquad(0\leq k\leq K,\ m\gg0).
> \tag{3.8}
> \]
> Homogeneous parts with negative indices are understood to be zero.

Indeed, the top homogeneous part of the pure identity is
\[
 \Lambda_r^m(P_e^m)=0.
 \tag{3.9}
\]
The arbitrary-degree split-symbol theorem gives a Newton separator for
the full generic translated support of \(P_e\) relative to the factors of
\(\Lambda_r\).

Expand a mixed term by homogeneous operator and polynomial pieces.  Its
deficit from the maximum possible output degree
\[
 D_m=(e-r)m+q
 \tag{3.10}
\]
is the sum of three nonnegative integers:

1. the excess of the selected operator orders above \(rm\);
2. the total degree lost by replacing copies of \(P_e\) with lower pieces
   of \(P\); and
3. the degree lost by replacing \(Q_q\) with a lower piece of \(Q\).

At deficit at most \(K\), only a bounded number of higher operator pieces
and lower polynomial pieces can occur.  After the Leibniz rule, every
summand is therefore a bounded defect of the translated
\((\Lambda_r,P_e)\) coefficient.  The same bounded-defect Newton-gap
argument used above kills every such summand for large \(m\), uniformly
over the finitely many defect patterns.  This proves (3.8).

Theorem 1.1 is now also a quick corollary.  If \(e<r\), total order kills
the mixed expression.  If \(e=r\), then
\(\deg\Lambda^m(QP^m)\leq q\); choosing \(K=q\) in (3.8) kills every
possible homogeneous layer.  For \(e>r\), however, the number of output
layers grows linearly with \(m\).  Thus any remaining counterexample must
have a mixed defect whose depth below the leading degree tends to infinity.
No fixed associated-graded face can carry the failure.

The entire binary linear-plus-quadratic island can also be closed.

> **Theorem 3.2 — binary drift--diffusion theorem.**  Let
> \(\Lambda=\Lambda_1+\Lambda_2\) be a binary constant-coefficient
> operator with \(\Lambda_1\ne0\), and let \(P\) be arbitrary.  The first
> two equations
> \[
>  \Lambda(P)=\Lambda^2(P^2)=0
> \tag{3.11}
> \]
> imply the GVC conclusion.

Choose coordinates with \(\Lambda_1=\partial_x\) and write
\[
 \Lambda
 =\partial_x
  +A\partial_x^2+B\partial_x\partial_y+C\partial_y^2.
 \tag{3.12}
\]
If \(C=0\), then
\[
 \Lambda=\partial_x(1+A\partial_x+B\partial_y)
\]
is factor-unit.  Injectivity of the differential unit turns
\(\Lambda(P)=0\) into \(\partial_xP=0\), and
\(\partial_x^m\) kills \(QP^m\) for \(m>\deg_xQ\).

Suppose \(C\ne0\), and write
\[
 P=\sum_{j=0}^N a_j(x)y^j,\qquad a_N\ne0.
 \tag{3.13}
\]
The coefficient of \(y^N\) in \(\Lambda(P)=0\) is
\[
 a_N'+Aa_N''=0.
 \tag{3.14}
\]
A polynomial solution is constant.  The coefficient of \(y^{N-1}\)
then shows that \(a_{N-1}\) is constant, and the next equations imply
\[
 \deg_y(P_{xy})\leq N-3,\qquad
 \deg_y(P_{xx})\leq N-4.
 \tag{3.15}
\]

Let
\[
 \Gamma(F,G)
 =A F_xG_x+\frac B2(F_xG_y+F_yG_x)+C F_yG_y
 \tag{3.16}
\]
be the carré du champ of the quadratic part.  Since
\(\Lambda(P)=0\) and derivatives commute with \(\Lambda\),
\[
 \Lambda(P^2)=2\Gamma(P,P).
\]
A second application of the product rule gives
\[
\begin{aligned}
 \Lambda^2(P^2)=4\bigl(
 &A^2P_{xx}^2+2ABP_{xx}P_{xy}
 +\frac{B^2}{2}P_{xx}P_{yy}
 +\left(\frac{B^2}{2}+2AC\right)P_{xy}^2\\
 &+2BCP_{xy}P_{yy}+C^2P_{yy}^2
 \bigr).
 \tag{3.17}
\end{aligned}
\]
If \(N\geq2\), equations (3.15)--(3.17) show that the coefficient of
\(y^{2N-4}\) is uniquely
\[
 4C^2N^2(N-1)^2a_N^2\ne0,
 \tag{3.18}
\]
contradicting the second pure equation.  Hence \(N\leq1\).  Writing
\(P=a(x)y+b(x)\) in the first equation shows successively that the
polynomial solutions \(a,b\) are constant.  Thus
\[
 P=ay+b.
 \tag{3.19}
\]

Put \(q_x=\deg_xQ\) and \(q_y=\deg_yQ\).  Every factor of \(\Lambda\)
other than \(C\partial_y^2\) contains at least one \(\partial_x\), so at
most \(q_x\) such factors can contribute to
\(\Lambda^m(QP^m)\).  The remaining factors require at least
\(2(m-q_x)\) \(y\)-derivatives, while
\(\deg_y(QP^m)\leq m+q_y\).  Therefore
\[
 \boxed{
 \Lambda^m(QP^m)=0\quad\text{for }m>2q_x+q_y.
 }
 \tag{3.20}
\]

Consequently a remaining binary counterexample whose lowest order is one
must contain an operator piece of order at least three.

The separated higher-order drift class is also safe.

> **Theorem 3.3 — separated drift theorem.**  Let
> \[
>  \Lambda=\partial_x+h(\partial_y),\qquad
>  h(T)=c_rT^r+c_{r+1}T^{r+1}+\cdots\in k[[T]],
>  \quad c_r\ne0,\quad r\ge2.
> \tag{3.21}
> \]
> For arbitrary \(P\), the first two equations
> \(\Lambda(P)=\Lambda^2(P^2)=0\) imply the GVC conclusion.

The formal series causes no convergence issue: on each polynomial only
finitely many derivative terms act.

Put \(p(y)=P(0,y)\).  The first equation has the unique polynomial
solution
\[
 P=e^{-x h(\partial_y)}p(y),
 \tag{3.22}
\]
where the exponential truncates on polynomials.  For a constant-
coefficient one-variable operator \(D^j\), define its product defect
\[
 \delta_j(f,g)
 =D^j(fg)-fD^jg-gD^jf
 =\sum_{i=1}^{j-1}\binom jiD^if\,D^{j-i}g.
 \tag{3.23}
\]
Since every \(y\)-derivative of \(P\) is again annihilated by \(\Lambda\),
evaluating the second moment at \(x=0\) gives a sum of iterated product
defects.  If
\[
 p(y)=a_Ny^N+\text{lower terms},
\]
the highest possible \(y\)-degree is \(2N-2r\).  It comes only from using
the lowest piece \(c_rD^r\) in both defects, and its coefficient is
\[
 c_r^2a_N^2\kappa_{r,N},
 \tag{3.24}
\]
where
\[
 \kappa_{r,N}
 =
 \sum_{i,j=1}^{r-1}
 \binom ri\binom rj
 (N)_{i+j}(N)_{2r-i-j}.
 \tag{3.25}
\]
For \(N\ge r\), the terms are nonnegative integers and the terms with
\(i+j=r\) are strictly positive.  Thus \(\kappa_{r,N}\ne0\) in
characteristic zero.  The second pure equation forces \(N<r\).

Now \(h(\partial_y)p=0\), so (3.22) reduces to \(P=p(y)\).  In a term of
\(\Lambda^m(QP^m)\), at most \(q_x=\deg_xQ\) factors can be
\(\partial_x\).  Each of the remaining \(h(\partial_y)\) factors has
order at least \(r\), whereas
\[
 \deg_y(QP^m)\leq(r-1)m+q_y.
\]
Consequently
\[
 \boxed{
 \Lambda^m(QP^m)=0
 \quad\text{for }m>r q_x+q_y.
 }
 \tag{3.26}
\]

This closes the entire lowest-order-one frontier.  It is useful to state
the formal straightening explicitly.

> **Theorem 3.4 — formal drift straightening.**  Every binary
> constant-coefficient operator whose symbol has nonzero linear part
> satisfies GVC for arbitrary \(P\).  In fact, after the first two pure
> equations the conclusion already follows.

Choose linear coordinates and rescale so that the linear part is \(\xi\).
The formal implicit-function theorem gives a unique
\(q(\eta)\in\eta^2k[[\eta]]\) with
\[
 A(-q(\eta),\eta)=0.
\]
Formal division by the monic series \(\xi+q(\eta)\) then gives
\[
 A(\xi,\eta)=U(\xi,\eta)(\xi+q(\eta)),
 \qquad U(0,0)=1.
\tag{3.27}
\]
Every formal constant-coefficient series acts locally finitely on
\(k[x,y]\).  Moreover, \(U(\partial)\) is an automorphism there: the
formal inverse of a series with constant term one is again locally finite.
Writing \(L=\partial_x+q(\partial_y)\), commutativity gives
\[
 \Lambda^m=U(\partial)^mL^m.
\]
Injectivity of \(U(\partial)^m\) turns the pure equations for \(\Lambda\)
into the corresponding pure equations for \(L\).  The formal-series
version of Theorem 3.3 gives eventual vanishing of \(L^m(QP^m)\) when
\(q\ne0\); when \(q=0\), the same conclusion is the immediate
\(\partial_x\)-degree cutoff.  Hence
\(\Lambda^m(QP^m)=U(\partial)^mL^m(QP^m)=0\) eventually.

The following finite-jet computations are retained as direct checks of
the straightening theorem.  Before straightening, every binary symbol
with nonzero linear part can be written after normalization as
\[
 \sigma_\Lambda(\xi_x,\xi_y)
 =\xi_x\Gamma(\xi_x,\xi_y)+h(\xi_y),
 \qquad\Gamma(0,0)=1.
\tag{3.27a}
\]

> **Proposition 3.5 — quadratic finite-jet check.**  Let
> \(\Lambda\) be any binary constant-coefficient operator with nonzero
> linear part, and let \(\deg P\leq2\).  The first two pure equations
> imply the GVC conclusion.

Use the decomposition (3.27), and let \(s\) be the lowest order of
\(h\), with \(s=\infty\) when \(h=0\).  If \(s\ge3\), then \(h(P)=0\),
so
\[
 0=\Lambda(P)=\Gamma\,\partial_xP.
\]
The differential unit \(\Gamma\) is injective, hence
\(\partial_xP=0\) and \(P=P(y)\).  Every \(h\)-factor has order at least
\(s>\deg P\), while every other factor contains \(\partial_x\); the usual
degree count proves eventual mixed vanishing.

It remains to take \(s=2\), writing the leading transverse term as
\(C\partial_y^2\), \(C\ne0\).  Only the order-at-most-two part of
\(\Lambda\) acts on \(P\), so the first equation gives
\[
 P=c(y^2-2Cx)+ey+f.
 \tag{3.28}
\]
In \(\Lambda^2(P^2)\), every operator pair of total order greater than
four kills \(P^2\).  The only possible new term beyond the
linear-plus-quadratic calculation is
\(2\Lambda_1\Lambda_3(P^2)\).  Writing
\[
 \Lambda_3=\partial_xK_2+c_3\partial_y^3,
\]
one has
\[
 \partial_x^2K_2(P^2)=0,\qquad
 \partial_x\partial_y^3(P^2)=0
\]
on (3.28).  Therefore the universal identity remains
\[
 \Lambda^2(P^2)=16C^2c^2.
 \tag{3.29}
\]
It forces \(c=0\), after which the mixed derivative count applies.

This is the degree-two finite-jet shadow of Theorem 3.4.

> **Proposition 3.6 — cubic finite-jet check.**  Let \(\Lambda\) be
> an arbitrary binary constant-coefficient operator with nonzero linear
> part, and let \(\deg P\le3\).  The first two pure equations imply the
> GVC conclusion.

Normalize the symbol to
\[
\begin{aligned}
 \Lambda={}&\partial_x
 +A\partial_x^2+B\partial_x\partial_y+C\partial_y^2\\
 &+D\partial_x^3+E\partial_x^2\partial_y
   +F\partial_x\partial_y^2+G\partial_y^3.
 \tag{3.30}
\end{aligned}
\]
Terms above order three do not act on \(P\).  Solving the first equation
for a general cubic therefore gives
\[
\begin{aligned}
 P={}&p_0+p_2y+p_5y^2+p_9y^3\\
 &+x(6BCp_9-2Cp_5-6Cp_9y-6Gp_9).
 \tag{3.31}
\end{aligned}
\]
Terms above order five cannot occur in the second moment.  Retaining the
complete order-four and order-five pieces shows that the second moment has
no \(x\)-dependence.  Its leading coefficients remain
\[
\begin{aligned}
 [y^2]\Lambda^2(P^2)&=144C^2p_9^2,\\
 [y]\Lambda^2(P^2)&=
 96Cp_9(-3BCp_9+Cp_5+9Gp_9).
 \tag{3.32}
\end{aligned}
\]
If \(C\ne0\), the first coefficient forces \(p_9=0\), after which the
constant coefficient is \(16C^2p_5^2\).  Hence \(p_5=0\) and
\(P=p_0+p_2y\).

If \(C=0\), the constant coefficient reduces to
\[
 648G^2p_9^2.
 \tag{3.33}
\]
For \(G\ne0\), this gives \(p_9=0\), so \(P\) has \(y\)-degree at most
two, below the transverse order three.  If \(G=0\), the first equation
already leaves \(P\in k[y]\) of degree at most three.  Either every
remaining operator term is divisible by \(\partial_x\), or the least pure
transverse order is at least four.  In both cases the usual eventual
derivative cutoff applies.  More explicitly, the only order-four
coefficient surviving in the unreduced constant term is the coefficient
of \(\partial_y^4\), multiplied by \(Cp_9^2\); every order-five
coefficient drops out.  Neither can alter the branch argument.

The direct quartic calculation continues in exactly the same way: only
the operator \(7\)-jet can enter the second moment, and its successive
branches begin with
\[
 2304C^4p_4^2,\qquad
 15552G^2p_4^2,\qquad
 39168L^2p_4^2,
\tag{3.34}
\]
where \(C,G,L\) are the pure transverse coefficients of orders two,
three, and four after the previous ones vanish.  This is a useful
finite-jet regression, but Theorem 3.4 subsumes it and every higher
polynomial degree.

The first surviving lowest-order-two cell is also safe.

> **Theorem 3.7 — quadratic-leading cubic-polynomial theorem.**  Suppose
> the lowest positive order of a binary constant-coefficient operator
> \(\Lambda\) is two and \(\deg P\le3\).  The first three pure equations
> imply the GVC conclusion.

We use the following elementary weighted cutoff.  Give \(x,y\) positive
weights \(u,v\), and give \(\partial_x,\partial_y\) the same weights.  If
every monomial of \(\Lambda\) has derivative weight strictly greater than
\(\deg_{u,v}P\), then, for every fixed \(Q\), every summand of
\(\Lambda^m(QP^m)\) has derivative weight greater than the weight of its
input once \(m\) is large.  Thus the mixed expression vanishes
eventually.

It is enough to work after scalar extension to an algebraic closure.
There are two nonzero binary quadratic orbits.  First take the double-line
orbit \(\Lambda_2=\partial_x^2\), and write
\[
 \Lambda_3=A\partial_x^3+B\partial_x^2\partial_y
 +C\partial_x\partial_y^2+D\partial_y^3.
\]
Solving \(\Lambda(P)=0\) gives
\[
\begin{aligned}
P={}&p_0+p_1x+p_2y+p_4xy+p_5y^2+p_8xy^2+p_9y^3\\
&-(Cp_8+3Dp_9)x^2.
\end{aligned}
\tag{3.35}
\]
Let \(E\) be the coefficient of \(\partial_y^4\).  The relevant
second-moment coefficients are
\[
 [y]\Lambda^2(P^2)=96Dp_8^2,
\tag{3.36}
\]
and, after \(D=0\),
\[
 \Lambda^2(P^2)=24p_8^2(C^2+4E).
\tag{3.37}
\]
On the only cancellation branch \(p_8\ne0\),
\(D=0\) and \(E=-C^2/4\).  Retaining the complete operator \(5\)-jet,
the third pure moment is exactly
\[
 \Lambda^3(P^3)=-4608C^3p_8^3.
\tag{3.38}
\]
Hence \(C=E=0\).  Choose \(v=1\) and \(2<u<3\).  Then
\(\deg_{u,v}P=u+2\), whereas \(\partial_x^2\), every surviving cubic and
quartic monomial, and every term of order at least five have strictly
larger weight.  The weighted cutoff applies.

If \(p_8=0\) but \(p_9\ne0\), the constant second-moment coefficient is
\(792D^2p_9^2\), so \(D=0\).  Taking \(v=1\) and
\(3/2<u<2\) makes \(\deg_{u,v}P=3\) and again puts every operator
monomial strictly above it.  If \(p_8=p_9=0\), Theorem 1.1 applies.

For the distinct-line orbit normalize
\(\Lambda_2=\partial_x\partial_y\).  Solving the first equation gives
\[
\begin{aligned}
P={}&p_0+p_1x+p_2y+p_3x^2+p_5y^2+p_6x^3+p_9y^3\\
&-6(Ap_6+Dp_9)xy.
\end{aligned}
\tag{3.39}
\]
The second moment begins with
\[
\begin{aligned}
[xy]\Lambda^2(P^2)&=72p_6p_9,\\
[x]\Lambda^2(P^2)&=24p_6(6Cp_9+p_5),\\
[y]\Lambda^2(P^2)&=24p_9(6Bp_6+p_3).
\end{aligned}
\tag{3.40}
\]
Thus the two pure cubic tips cannot coexist.  If \(p_6\ne0\), then
\(p_9=p_5=0\), and the remaining constant coefficient is
\(288A^2p_6^2\), so \(A=0\).  Weights \(u=1\), \(2<v<3\) give a strict
separator.  The branch \(p_9\ne0\) is symmetric, and if both vanish then
\(\deg P\le2\).  This proves the theorem.

## 4. The rank obstruction for natural conversions

Let \(B\) be the auxiliary coefficient algebra and let
\(\phi:B\to k\) be the final specialization.  Evaluation, restriction,
and polynomial substitution are algebra homomorphisms.  If the
substitutions keep dual variables on the dual side and coordinate
variables on the coordinate side, then
\[
 \phi\!\left(A(\zeta,u)P(z,u)\right)
 =
 \phi(A(\zeta,u))\,\phi(P(z,u)).
 \tag{4.1}
\]
Taking the bidegree-\((4,4)\) component gives the outer product
\[
 \phi(A)_4\otimes\phi(P)_4,
 \tag{4.2}
\]
whose coefficient matrix has rank at most one.  For \(s\) channels,
subadditivity of matrix rank gives rank at most \(s\).

The displayed witness has coefficient matrix
\[
 \begin{pmatrix}
 -1&2&0&0&0\\
 -3/2&2&6&0&0\\
 -1/2&3/2&6&6&0\\
 0&1&3/2&2&2\\
 0&0&-1/2&-3/2&-1
 \end{pmatrix},
 \qquad \det=48.
 \tag{4.3}
\]
Its rank is five, proving Theorem 1.3.

For comparison, a coefficient functional such as
\(\ell_j(G)=[u^j]G\) is not an algebra homomorphism:
\[
 \ell_j(u^a u^b)=\mathbf1_{a+b=j},
 \qquad
 \ell_j(u^a)\ell_j(u^b)
 =\mathbf1_{a=j}\mathbf1_{b=j}.
 \tag{4.4}
\]
Thus coefficient extraction may sum several rank-one coefficient
channels, but in general
\[
 \ell_j(H^m)\ne\ell_j(H)^m.
 \tag{4.5}
\]
An auxiliary polarization using such a functional needs a new all-order
identity; specialization alone supplies none.  With one auxiliary
variable,
\[
 [u^j]\bigl(A(u,\zeta)P(u,z)\bigr)
 =\sum_{a=0}^j A_a(\zeta)P_{j-a}(z)
 \tag{4.6}
\]
has rank at most \(j+1\), so reaching (4.3) requires \(j\geq4\).
For a multi-index \(\alpha\), the analogous bound is
\(\prod_i(\alpha_i+1)\geq5\).

## 5. A converse for single-coefficient polarization

There is a second obstruction on the operator side.  Suppose a
constant-coefficient operator \(\Lambda\) admits a linear-translation,
single-coefficient formula
\[
 \Lambda^mR(z)
 =c_m[t_1^{m\alpha_1}\cdots t_s^{m\alpha_s}]
 R\left(z+\sum_{i=1}^s t_iv_i\right)
 \tag{5.1}
\]
for every polynomial \(R\) and every \(m\geq1\), where the multi-index
\(\alpha\) and the directions \(v_i\) are fixed and \(c_m\ne0\).

> **Theorem 5.1 — single-coefficient polarization converse.**  Every
> operator admitting (5.1) has split homogeneous symbol
> \[
>  \sigma_\Lambda(\xi)
>  =c\prod_{i=1}^s\langle v_i,\xi\rangle^{\alpha_i}.
> \tag{5.2}
> \]
> Conversely, every operator with symbol (5.2) admits (5.1), with explicit
> factorial constants.

For necessity, it is enough to take \(m=1\) and apply both sides formally
to \(R(z)=\exp\langle z,\xi\rangle\).  After cancelling the exponential,
\[
 \sigma_\Lambda(\xi)
 =\frac{c_1}{\prod_i\alpha_i!}
  \prod_i\langle v_i,\xi\rangle^{\alpha_i}.
 \tag{5.3}
\]
Conversely, Taylor expansion gives
\[
 [t^{m\alpha}]R(z+Vt)
 =
 \frac{\prod_iD_{v_i}^{m\alpha_i}R(z)}
      {\prod_i(m\alpha_i)!},
 \tag{5.4}
\]
which proves the claim.

Thus the translated Laurent-power proof is not merely one convenient
implementation inside this architecture: it exhausts every fixed linear
translation followed by one diagonal coefficient.  An irreducible
nonhomogeneous symbol such as
\[
 \xi_x+\xi_y^2
 \tag{5.5}
\]
cannot be reached by this construction.  Any continuation for such a
heat-type operator needs multiple coupled coefficient functionals,
nonlinear translation, or a different intertwiner; none inherits a single
Laurent-power Newton separator automatically.

## 6. Consequences for the requested escape routes

1. **Auxiliary polarization followed by evaluation.**  This is
   multiplicative and remains rank one, so it cannot produce \(F\).
   Coefficient extraction must have width at least five and loses formal
   compatibility with powers.  On the operator side, a single diagonal
   coefficient after linear translation is equivalent to symbol splitting
   by Theorem 5.1.
2. **Rank-one dilation.**  Every separated restriction of one rank-one
   dilation remains rank one.  At least five additive channels are needed
   even before imposing moment identities.
3. **Nonlinear Segre projection.**  A nonlinear substitution applied
   separately to the two factors still has the product form (4.1).
   A map that mixes the dual and coordinate algebras is outside this
   theorem and must separately prove that it intertwines contraction.
4. **Nonhomogeneous constant-coefficient operators.**  Theorem 1.1
   excludes every pair with \(\deg P\) at most the lowest positive
   operator order.  Theorem 1.2 excludes arbitrary polynomial degree when
   the operator is a split homogeneous factor times a differential unit.
5. **Degree raising.**  Homogeneous binary operators now satisfy GVC for
   arbitrary \(P\), so \(\deg P>\operatorname{ord}\Lambda\) is closed in
   that class.  The surviving frontier requires a genuinely
   nonhomogeneous operator, \(\deg P>r\), and failure of the factor-unit
   architecture.  Formal drift straightening closes every symbol with
   nonzero linear part in every polynomial degree.  Such a witness still
   cannot be
   obtained from \(F\) by a separated multiplicative conversion, because every
   operator-polynomial symbol remains a single product
   \(A(\zeta)P(z)\).

No two-variable GVC counterexample or improved ordinary-Laplacian endpoint
is claimed.

## 7. Remaining exact continuations

The results leave four sharply separated directions.

1. **Growing-depth defect cone.**  Theorem 3.1 kills every fixed
   filtration depth.  Record operator-order excess and polynomial-degree
   loss as a two-dimensional semigroup and study its normalized limit
   cone.  A separator positive on every nonzero limiting slope would
   promote fixed-layer vanishing to full GVC; a zero-slope face would
   identify the only possible asymptotic counterexample architecture.
2. **Lowest-order-two anchor.**  Theorem 3.4 removes the entire
   lowest-order-one class.  The first unresolved symbols therefore have
   lowest positive order \(r\ge2\), are genuinely nonhomogeneous, have
   \(\deg P>r\), and do not admit the split factor-unit architecture.
   Theorem 3.7 additionally closes \(r=2,\deg P\le3\), so the first
   exact cell is \(r=2,\deg P=4\).  Theorem 3.1 has already removed every
   fixed leading face.
3. **Coupled coefficient polarization.**  A one-coefficient translation
   is exactly the split-symbol class, while coefficient extraction of the
   rank-five witness needs at least five channels.  The next algebraic
   question is to classify finite families of coefficient functionals
   whose convolution under powers intertwines contraction.  Merely
   matching \(F\) at \(m=1\) is insufficient.
4. **Mixed nonlinear intertwiners.**  A substitution mixing dual and
   coordinate variables can evade the separated rank obstruction, but it
   must carry the Weyl contraction map and the coordinate-only multiplier
   class explicitly.  This is the natural nonlinear Segre continuation.

The ordinary-Laplacian quadraticization of the Dvorsky symbol remains a
separate higher-dimensional route.  None of these continuations requires
or supplies a classification of general SIC rank strata.

## Reproduction

Run

```bash
python3 scripts/verify_separable_gvc_escape_obstructions.py
```

The dependency-free checker verifies the determinant in (4.3), the
rank bounds for separated channels, and the failure of multiplicativity
for coefficient extraction.  It also checks the translated polarization
identity with polynomial degree greater than operator order and a
degree-raising example for
\(\partial_x\partial_y(1+\partial_x+\partial_y^2)\).  The all-order proofs
are the Newton-separator and injectivity arguments, not the finite replay.
