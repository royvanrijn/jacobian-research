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
shows that cubic \(P\) is safe, and Theorem 3.9 below pushes the first
quadratic-leading polynomial degree to at least five.  Theorems
3.10--3.12 close all cubic-leading quartic symbols.  Consequently:

> **Corollary 1.4 — binary degree-four GVC.**  Every constant-coefficient
> operator in two variables satisfies GVC for every polynomial \(P\) of
> degree at most four.  In particular, every binary GVC counterexample
> has \(\deg P\ge5\).  If \(\deg P=5\), its lowest positive operator order
> is one of \(2,3,4\).

Theorems 3.14--3.16 below close the \(r=2,\deg P=5\) row as well.
Therefore:

> **Corollary 1.5 — first degree-five reduction.**  If a binary GVC
> counterexample has \(\deg P=5\), then its lowest positive operator
> order is \(3\) or \(4\).

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

The use of the third equation is sharp for this argument.  Explicitly,
\[
 \Lambda=\partial_x^2+\partial_x\partial_y^2
 -\frac14\partial_y^4,\qquad
 P=xy^2-x^2
\tag{3.38a}
\]
satisfies \(\Lambda(P)=\Lambda^2(P^2)=0\), but
\(\Lambda^3(P^3)=-4608\).

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

After three moments, the quartic cell has only one finite residual branch.

> **Proposition 3.8 — quadratic-leading quartic reduction.**  Suppose
> \(r=2\), \(\deg P=4\), and the pure equations hold through order three.
> If \(\Lambda_2\) has two distinct roots, the GVC conclusion follows.
> If \(\Lambda_2\) is a square, every branch except one has the GVC
> conclusion; the remaining branch is cut out by the two equations
> (3.45) below and projects to an explicit finite sextic set.

For the distinct-root orbit put
\(\Lambda_2=\partial_x\partial_y\).  The leading equation and the second
leading pure equation force \(P_4\) to be a single pure tip; take
\(P_4=x^4\).  Let \(A,E\) be the coefficients of
\(\partial_x^3,\partial_x^4\), and let \(c,d\) be the coefficients of
\(y^3,y^2\) in \(P\).  Successive coefficients of the full second moment
are
\[
\begin{aligned}
[x^2y]\Lambda^2(P^2)&=144c,\\
[x^2]\Lambda^2(P^2)&=48(132A^2+d),\\
[y]\Lambda^2(P^2)\big|_{c=0,\ d=-132A^2}&=9216A^3,\\
[1]\Lambda^2(P^2)\big|_{c=d=A=0}&=31104E^2.
\end{aligned}
\tag{3.41}
\]
Thus \(c=d=A=E=0\), and the first equation reduces \(P\) to a polynomial
in \(x\) plus a linear \(y\)-term.  Weights \(u=1\), \(3<v<4\) strictly
separate every operator monomial from \(P\).  The \(P_4=y^4\) branch is
symmetric.

Now put \(\Lambda_2=\partial_x^2\).  Its leading equation gives
\[
 P_4=y^3(ax+by).
\tag{3.42}
\]
The stabilizer of \(\partial_x^2\) reduces this to \(xy^3\) or \(y^4\).
On the \(xy^3\) branch, denote the coefficients of
\[
\partial_x\partial_y^2,\ \partial_y^3,\ \partial_y^4,\
\partial_x\partial_y^3,\ \partial_y^5,\ \partial_y^6
\]
by \(C,D,E,H,J,K\), respectively.  The second moment successively gives
\[
 D=0,\qquad 20E+C^2=0,\qquad EC=0,\qquad
 J=0,\qquad K=-\frac{17}{40}H^2.
\tag{3.43}
\]
Hence \(C=E=0\), and the complete third moment reduces to
\[
 \Lambda^3(P^3)=-3604176H^3.
\tag{3.44}
\]
Thus \(H=K=0\).  Weights \(v=1\), \(3<u<4\) give the strict cutoff.

It remains to take \(P_4=y^4\).  The second moment first forces \(D=0\).
Put \(z=[xy^2]P\).  Its remaining constant coefficient and the third
moment are, up to nonzero scalar factors,
\[
\begin{aligned}
S(C,E,z)={}&1728E^2-48C^2E+112CEz+4Ez^2\\
&-4C^3z+C^2z^2,\\
T(C,E,z)={}&815616E^3-13824C^2E^2+52416CE^2z
 +1152E^2z^2\\
&-1584C^3Ez+936C^2Ez^2+36CEz^3
-36C^4z^2+C^3z^3.
\end{aligned}
\tag{3.45}
\]
If \(Cz=0\), these equations force \(E=0\), and weights
\(v=1\), \(2<u<3\) close the branch.  If \(Cz\ne0\), eliminating \(E\)
gives
\[
\begin{aligned}
0={}&-2583360C^6+1828368C^5z+1514304C^4z^2\\
&+502328C^3z^3+80916C^2z^4+6117Cz^5+92z^6.
\end{aligned}
\tag{3.46}
\]
Thus, modulo the natural scaling, only finitely many algebraic ratios
\(z/C\) remain after the first three moments.

> **Theorem 3.9 — quadratic-leading quartic theorem.**  Suppose the
> lowest positive order of \(\Lambda\) is two and \(\deg P\le4\).  The
> first four pure equations imply the GVC conclusion.

Only the last branch of Proposition 3.8 needs consideration.  Give
\(x,y\) weights \(2,1\).  The first equation and \(D=0\) give
\[
\begin{aligned}
P&=P_{[4]}+P_{<4},&
P_{[4]}&=y^4+zxy^2-(12E+Cz)x^2,\\
\Lambda&=\Lambda_{[4]}+\Lambda_{>4},&
\Lambda_{[4]}&=\partial_x^2+C\partial_x\partial_y^2
 +E\partial_y^4.
\end{aligned}
\tag{3.47}
\]
The subscripts denote these weights.  Every omitted polynomial term has
weight below four, and every omitted operator monomial has weight above
four.  Since \(P^m\) has weight at most \(4m\) while every monomial of
\(\Lambda^m\) has derivative weight at least \(4m\), equality is possible
only by choosing \(P_{[4]}\) and \(\Lambda_{[4]}\) throughout.  Hence,
for every \(m\),
\[
 \Lambda^m(P^m)=\Lambda_{[4]}^m(P_{[4]}^m).
\tag{3.47a}
\]

The second and third right-hand sides are the polynomials \(S,T\) in
(3.45).  The fourth is \(17280M(C,E,z)\).  Eliminating \(E\) from \(S,M\)
and removing a nonzero scalar and the factor \(C^4z^4\) gives
\[
\begin{aligned}
R_8(C,z)={}&-2557094400C^8+639596160C^7z
 +532247424C^6z^2\\
&+236419896C^5z^3+52199245C^4z^4+7972150C^3z^5\\
&+720528C^2z^6+39490Cz^7+550z^8.
\end{aligned}
\tag{3.48}
\]
After setting \(t=z/C\), the degree-six polynomial \(R_6(1,t)\) from
(3.46) and \(R_8(1,t)\) have gcd one.  More explicitly, their resultant is
\[
-22002331580862445954532620608845574194895939575073794373253473026129579212800.
\tag{3.49}
\]
Thus no branch with \(Cz\ne0\) survives all three equations
\(S=T=M=0\); (3.47a) proves that no higher jet or lower coefficient can
alter this obstruction.  The cases \(Cz=0\) were already closed by strict
weights in Proposition 3.8.  Together with Theorem 3.7 and the other
branches of Proposition 3.8, this proves the theorem.

The triple-root part of the next leading order also closes.

> **Theorem 3.10 — triple-root cubic-leading quartic theorem.**  Suppose
> the lowest positive order of \(\Lambda\) is three,
> \(\Lambda_3\) has a triple root, and \(\deg P\le4\).  The first four
> pure equations imply the GVC conclusion.

Normalize \(\Lambda_3=\partial_x^3\).  The leading equation gives
\[
 P_4=y^2(Ax^2+Bxy+Cy^2).
\tag{3.50}
\]
The stabilizer of \(\partial_x^3\) leaves three cases.

If \(A\ne0\), normalize \(A=1,B=0\), leaving
\(P_4=x^2y^2+\rho y^4\).  Denote the coefficients of
\[
\partial_x^2\partial_y^2,\ \partial_x\partial_y^3,\
\partial_y^4,\ \partial_x\partial_y^4,\
\partial_y^5,\ \partial_y^6
\]
by \(U,D,E,V,J,K\), respectively.  Successive coefficients of the second
moment force
\[
 E=0,\qquad D=0,\qquad V=-\frac29U^2.
\tag{3.51}
\]
The third moment then gives
\[
1555200J=0,\qquad
1920(810K-272U^3)=0,
\tag{3.52}
\]
so \(J=0\) and \(K=136U^3/405\).  The fourth moment is
\[
 \Lambda^4(P^4)=3361505280U^4.
\tag{3.53}
\]
Hence \(U=0\).

These identities are full-jet statements.  With weights
\(\operatorname{wt}(x)=2,\operatorname{wt}(y)=1\), the displayed
polynomial and operator terms in (3.51)--(3.52) are precisely their
weight-six faces; every omitted polynomial term has smaller weight and
every omitted operator term has larger weight.  Thus the same face
identity as (3.47a) applies.  Once \(U=0\), the face is the one-sided pair
\((\partial_x^3,x^2y^2)\), whose mixed values have the direct
\(x\)-degree cutoff.

If \(A=0,B\ne0\), a shear gives \(P_4=xy^3\).  Write \(D,E\) for the
coefficients of \(\partial_x\partial_y^3,\partial_y^4\).  The
\(y\)-coefficient of the third moment is \(6531840E^2\), so \(E=0\);
the remaining constant second-moment coefficient is \(1584D^2\).
Weights \((3,2)\) then leave the one-sided face
\((\partial_x^3,xy^3)\).

Finally, for \(P_4=y^4\), the constant second moment is
\(49536E^2\), so \(E=0\).  Weights \((4,3)\) leave
\((\partial_x^3,y^4)\).  Every branch therefore has an eventual mixed
cutoff.

The double-root orbit has only one endpoint left.

> **Theorem 3.11 — double-root cubic-leading quartic theorem.**
> Suppose \(r=3\), \(\deg P=4\), and
> \(\Lambda_3\) has one double and one simple root.  The first three
> pure equations imply the GVC conclusion.

The first two leading-face equations give
\[
 P_4=x^4
 \quad\text{or}\quad
 P_4=y^3(ax+by).
\tag{3.54}
\]
For \(P_4=x^4\), let \(T\) be the coefficient of \(\partial_x^4\).
Solving the first equation and reading the second gives successively
\[
 [y^3]P=0,\qquad [xy^2]P=0,\qquad
 [y^2]P=-132T^2.
\tag{3.55}
\]
With weights \(\operatorname{wt}(x)=1,\operatorname{wt}(y)=2\), the
weight-four face is
\[
\Lambda_{[4]}=\partial_x^2\partial_y+T\partial_x^4,\qquad
P_{[4]}=x^4-12Tx^2y-132T^2y^2.
\tag{3.56}
\]
All omitted operator terms have larger weight and all omitted polynomial
terms have smaller weight.  Its third moment is
\[
 \Lambda_{[4]}^3(P_{[4]}^3)=129392640T^3,
\tag{3.57}
\]
so \(T=0\), leaving the one-sided pair
\((\partial_x^2\partial_y,x^4)\).

Now take \(a\ne0\), normalize \(a=1\), and write
\(P_4=xy^3+by^4\).  Let \(E,H,J\) be the coefficients of
\(\partial_y^4,\partial_x\partial_y^3,\partial_y^5\).
The second equation first gives
\[
 [x^3]P=-10E,\qquad
 20J=-(248E^2b^2+56EHb-20EU+H^2),
\tag{3.58}
\]
where \(U=[\partial_x^2\partial_y^2]\Lambda\).  The \(x\)-coefficient of
the third moment is \(7827840E^2\), hence \(E=0\).  Its remaining
constant coefficient is then a nonzero multiple of \(H^3\), so
\(H=J=0\).  Choosing \(\operatorname{wt}(y)=1\) and
\(2<\operatorname{wt}(x)<3\) strictly separates every operator monomial
from \(P\), proving the mixed cutoff.

It remains to take \(a=0\), so \(P_4=y^4\).  Put
\[
r=[x^3]P,\qquad z=[xy^2]P,
\]
and retain \(E=[\partial_y^4]\Lambda\) and
\(H=[\partial_x\partial_y^3]\Lambda\).  The relevant coefficients are
\[
\Lambda^2(P^2)=96(372E^2+6Hr+rz),
\qquad
[y]\Lambda^3(P^3)=51840r^2.
\tag{3.59}
\]
Thus \(r=0\), followed by \(E=0\); the first equation also gives
\([x^2y]P=-12E=0\).  With \(\operatorname{wt}(y)=1\) and
\(3/2<\operatorname{wt}(x)<2\), \(P\) has weight at most four and every
operator monomial has weight strictly greater than four.  This closes the
last branch and proves the theorem.

The squarefree cubic orbit unexpectedly reduces to the same \(x^4\)
weighted face as (3.56).

> **Theorem 3.12 — squarefree cubic-leading quartic theorem.**
> Suppose \(r=3\), \(\deg P=4\), and \(\Lambda_3\) is squarefree.  The
> first three pure equations imply the GVC conclusion.

After scalar extension and a linear change of variables, normalize
\[
 \Lambda_3=\partial_x^2\partial_y+\partial_x\partial_y^2.
\tag{3.60}
\]
Write
\[
 P_4=A x^4+a_1x^3y+B x^2y^2+a_3xy^3+C y^4.
\]
The first leading moment gives
\[
 a_1=a_3=-\frac23B.
\tag{3.61}
\]
The coefficients of the second leading moment have Gröbner basis
\[
 AB-BC,\qquad B^2-6BC,\qquad 6AC-BC.
\tag{3.62}
\]
Thus \(P_4\) is proportional to exactly one of
\[
 x^4,\qquad y^4,\qquad (x-y)^4.
\tag{3.63}
\]
The stabilizer of the three roots of (3.60) permutes these three
fourth powers transitively, so it suffices to take \(P_4=x^4\).

Let \(T=[\partial_x^4]\Lambda\), and write
\[
 p_{21}=[x^2y]P,\quad p_{12}=[xy^2]P,\quad
 p_{03}=[y^3]P,\quad p_{02}=[y^2]P,\quad p_{30}=[x^3]P.
\]
The first equation is
\[
 24T+2p_{21}+2p_{12}=0.
\tag{3.64}
\]
After substituting \(p_{21}=-12T-p_{12}\), the nonconstant
coefficients of the second moment are
\[
 [x]\Lambda^2(P^2)=96(5p_{12}+6p_{03}),\qquad
 [y]\Lambda^2(P^2)=288p_{03}.
\tag{3.65}
\]
They force \(p_{03}=p_{12}=0\).  The remaining constant coefficient is
\[
 48\bigl(264T^2+2p_{02}\bigr),
\tag{3.66}
\]
so \(p_{02}=-132T^2\) and \(p_{21}=-12T\).

Give \(x,y\) weights \(1,2\).  The weight-four faces are now exactly
\[
\begin{aligned}
 \Lambda_{[4]}&=\partial_x^2\partial_y+T\partial_x^4,\\
 P_{[4]}&=x^4-12Tx^2y-132T^2y^2.
\end{aligned}
\tag{3.67}
\]
The additional squarefree summand
\(\partial_x\partial_y^2\) has weight five; every other omitted operator
term has weight greater than four, while every omitted polynomial term
has weight less than four.  Since \(P^m\) has weight at most \(4m\) and
every monomial of \(\Lambda^m\) has derivative weight at least \(4m\),
only equality can contribute.  Hence, for every \(m\),
\[
 \Lambda^m(P^m)=\Lambda_{[4]}^m(P_{[4]}^m).
\tag{3.67a}
\]
The third moment of this face is
\[
 \Lambda_{[4]}^3(P_{[4]}^3)=129392640T^3.
\tag{3.68}
\]
Thus \(T=0\), and the remaining face
\((\partial_x^2\partial_y,x^4)\) is one-sided.  The strict weighted
degree count gives the eventual mixed cutoff, proving the theorem.

We can now prove Corollary 1.4.  A nonzero order-zero term forces \(P=0\).
Lowest order one is closed in arbitrary degree by Theorem 3.4.  Lowest
order two is closed through degree four by Theorem 3.9.  At lowest order
three, degree below four is covered by Theorem 1.1, while every cubic
leading symbol becomes triple-root, double-root, or squarefree after
scalar extension; Theorems 3.10--3.12 cover these cases in degree four.
Lowest order at least four is again covered by Theorem 1.1.  Vanishing
after scalar extension descends to the ground field.  Finally, in degree
five Theorem 1.1 excludes lowest order at least five and Theorem 3.4
excludes lowest order one, leaving only \(r=2,3,4\).

The first unresolved quadratic-leading cell already has a finite
top-form reduction.

> **Proposition 3.13 — quadratic-leading quintic top forms.**
> Suppose \(r=2\), \(\deg P=5\), and the first two pure equations hold.
> After scalar extension, a linear change of variables, and rescaling,
> the leading pair \((\Lambda_2,P_5)\) is one of
> \[
> (\partial_x\partial_y,x^5),\qquad
> (\partial_x^2,xy^4),\qquad
> (\partial_x^2,y^5).
> \tag{3.69}
> \]

Indeed, if \(\Lambda_2=\partial_x\partial_y\), the first leading equation
gives \(P_5=Ax^5+By^5\).  The second leading moment is
\[
 (\partial_x\partial_y)^2(P_5^2)
 =800ABx^3y^3,
\tag{3.70}
\]
so exactly one of \(A,B\) is nonzero; swapping \(x,y\) gives the first
normal form.  If \(\Lambda_2=\partial_x^2\), the first leading equation
gives \(P_5=y^4(Ax+By)\).  When \(A\ne0\), a shear preserving
\(\partial_x^2\) removes \(B\), giving \(xy^4\); when \(A=0\), one gets
\(y^5\).  This proposition classifies only the leading top forms.  The
lower polynomial terms and higher operator jets are treated in Theorems
3.14--3.16.

The distinct-root tip in (3.69) closes without a third moment.

> **Theorem 3.14 — distinct-root quadratic-leading quintic theorem.**
> Suppose \(r=2\), \(\deg P=5\), and \(\Lambda_2\) has two distinct
> roots.  The first two pure equations imply the GVC conclusion.

Normalize \(\Lambda_2=\partial_x\partial_y\) and \(P_5=x^5\).  Retain
every operator jet through order eight, since higher jets kill \(P^2\).
After solving the first equation, put
\[
\begin{gathered}
 A=[\partial_x^3]\Lambda,\quad
 B=[\partial_x^4]\Lambda,\quad
 C=[\partial_x^5]\Lambda,\\
 u=[y^4]P,\quad v=[y^3]P,\quad w=[y^2]P.
\end{gathered}
\]
Successive coefficients of the full second moment are
\[
\begin{aligned}
 [x^3y^2]\Lambda^2(P^2)&=480u,\\
 [x^4]\Lambda^2(P^2)\big|_{u=0}&=48000A^2,\\
 [x^3y]\Lambda^2(P^2)\big|_{u=A=0}&=240v,\\
 [x^3]\Lambda^2(P^2)\big|_{u=A=v=0}&=80w,\\
 [x^2]\Lambda^2(P^2)\big|_{u=A=v=w=0}&=1296000B^2,\\
 [1]\Lambda^2(P^2)\big|_{u=A=v=w=B=0}&=3340800C^2.
\end{aligned}
\tag{3.71}
\]
Thus all six displayed parameters vanish.  Substitution back into the
first equation leaves
\[
 P=f(x)+ay,\qquad \deg f=5,
\tag{3.72}
\]
and no pure-\(x\) operator term of order below six.  Consequently
\[
 \Lambda=\partial_y\Gamma+H(\partial_x),
\qquad
\operatorname{ord}_{\min}\Gamma\ge1,\quad
\operatorname{ord}_{\min}H\ge6,
\tag{3.73}
\]
where \(H=0\) is allowed.

This normal form gives the all-order cutoff directly.  Fix \(Q\), put
\(q=\deg_yQ\), and expand
\[
 \Lambda^m
 =\sum_{k=0}^m\binom mk
   (\partial_y\Gamma)^kH^{m-k}.
\tag{3.74}
\]
Set \(N=m-k\).  In the Leibniz expansion of
\(\partial_y^k(QP^m)\), at most \(q\) derivatives hit \(Q\), and every
remaining term contains \(P^{N+j}\) with \(0\le j\le q\).  Its
\(x\)-degree is at most
\[
 \deg_xQ+5(N+q).
\tag{3.75}
\]
But every monomial of \(H^N\) has \(x\)-derivative order at least \(6N\).
Hence the term vanishes whenever
\[
 N>\deg_xQ+5q.
\tag{3.76}
\]
Only bounded \(N\) remain.  For those terms, applying
\(\partial_y^kH^N\) leaves a polynomial of degree bounded independently
of \(m\), while \(\Gamma^k\) has differential order at least
\(k=m-N\).  It therefore kills the result for all sufficiently large
\(m\).  This proves the theorem.

The cutoff argument did not use the special value five.

> **Corollary 3.14a — transverse-linear/high-order class.**  Let
> \[
> P=f(x)+ay,\qquad
> \Lambda=\partial_y\Gamma+H(\partial_x),
> \tag{3.76a}
> \]
> where \(\Gamma\) has no order-zero term and every nonzero monomial of
> \(H\) has order strictly greater than \(\deg f\).  Then all pure moments
> vanish and the GVC mixed conclusion holds for every \(Q\).

Indeed, the proof of (3.75)--(3.76), with \(5\) replaced by
\(\deg f\) and \(6\) by \(\operatorname{ord}_{\min}H\), applies verbatim.
For \(Q=1\), every term with \(N>0\) is killed by \(H^N\), while the
\(N=0\) term becomes constant after \(\partial_y^m\) and is killed by
\(\Gamma^m\).

The first double-line tip has two successive weighted faces.

> **Theorem 3.15 — \(xy^4\) quadratic-leading quintic theorem.**
> Suppose \(r=2\), \(\deg P=5\),
> \(\Lambda_2=\partial_x^2\), and \(P_5=xy^4\).  The first three pure
> equations imply the GVC conclusion.

Retain every operator jet through order eight and solve the first
equation.  Write \(a_{ij}=[\partial_x^i\partial_y^j]\Lambda\).  Successive
coefficients of the second moment first give
\[
\begin{aligned}
 [y^5]\Lambda^2(P^2)&=1152a_{03},\\
 [y^4]\Lambda^2(P^2)\big|_{a_{03}=0}
   &=96(68a_{04}+a_{12}^2),\\
 [xy^2]\Lambda^2(P^2)
   \big|_{a_{03}=0,\ a_{04}=-a_{12}^2/68}
   &=-\frac{115200}{17}a_{12}^3.
\end{aligned}
\tag{3.77}
\]
Hence \(a_{03}=a_{12}=a_{04}=0\).  The next three coefficients give
\[
\begin{aligned}
 a_{05}&=0,\\
 a_{06}&=-\frac{23}{70}a_{13}^2,\\
 a_{07}&=\frac{16}{35}a_{13}^2a_{21}
          -\frac{11}{14}a_{13}a_{14}.
\end{aligned}
\tag{3.78}
\]

Give \(x,y\) weights \(3,1\) and put \(H=a_{13}\).  The maximal
polynomial face and minimal operator face are
\[
\begin{aligned}
 P_{[7]}&=xy^4-12Hx^2y,\\
 \Lambda_{[6]}&=\partial_x^2
   +H\partial_x\partial_y^3
   -\frac{23}{70}H^2\partial_y^6.
\end{aligned}
\tag{3.79}
\]
Every omitted polynomial term has weight below seven and every omitted
operator term has weight above six.  Therefore the weight-three component
of the third pure equation is the corresponding component of this face.
Directly,
\[
 [y^3]\Lambda_{[6]}^3(P_{[7]}^3)
 =-553153536H^3.
\tag{3.80}
\]
Thus \(H=0\).  Equations (3.78) give \(a_{06}=a_{07}=0\), and the
remaining constant second-moment coefficient gives
\[
 a_{08}=-\frac{67}{140}a_{14}^2.
\tag{3.81}
\]

Now give \(x,y\) weights \(4,1\) and put \(J=a_{14}\).  The new faces are
\[
\begin{aligned}
 P_{[8]}&=xy^4-12Jx^2,\\
 \Lambda_{[8]}&=\partial_x^2
   +J\partial_x\partial_y^4
   -\frac{67}{140}J^2\partial_y^8.
\end{aligned}
\tag{3.82}
\]
Their third moment is the constant
\[
 \Lambda_{[8]}^3(P_{[8]}^3)=-5430509568J^3.
\tag{3.83}
\]
Hence \(J=0\).  What remains has maximal polynomial face \(xy^4\),
minimal operator face \(\partial_x^2\), and strict inequalities off those
faces.  Since the \(x\)-degree of \((xy^4)^m\) is only \(m<2m\), the
standard weighted one-sided count gives the eventual mixed cutoff.

The last quadratic-leading quintic tip reduces to six weighted ratios,
four of which die at moment three.

> **Theorem 3.16 — \(y^5\) quadratic-leading quintic theorem.**
> Suppose \(r=2\), \(\deg P=5\),
> \(\Lambda_2=\partial_x^2\), and \(P_5=y^5\).  The first three pure
> equations imply the GVC conclusion.

Solve the first equation and retain all jets through order eight.  The
\(y^4\)-coefficient of the second moment is
\[
 122400a_{03}^2,
\tag{3.84}
\]
so \(a_{03}=0\).  Put
\[
 A=a_{04},\qquad B=a_{12},\qquad z=[xy^3]P.
\]
Give \(x,y\) weights \(2,1\).  The maximal polynomial and minimal
operator faces are
\[
\begin{aligned}
 P_{[5]}&=y^5+zxy^3-3(20A+Bz)x^2y,\\
 \Lambda_{[4]}&=\partial_x^2
   +B\partial_x\partial_y^2+A\partial_y^4.
\end{aligned}
\tag{3.85}
\]
The weight-two component of the second moment is
\[
72\left(
16Az(40A+Bz)x+S(A,B,z)y^2
\right),
\tag{3.86}
\]
where
\[
\begin{aligned}
S={}&24000A^2-1200AB^2+880ABz+20Az^2\\
   &-60B^3z+B^2z^2.
\end{aligned}
\tag{3.87}
\]
The equations \(Az(40A+Bz)=S=0\), modulo the scaling
\[
 (A,B,z)\longmapsto(t^2A,tB,tz),
\tag{3.88}
\]
have six nonzero ratios:
\[
\begin{gathered}
(0,0,1),\quad(0,1,0),\quad(0,1,60),\quad
\left(\frac1{20},1,0\right),\\
\left(-\frac{t}{40},1,t\right),
\qquad t^2+12t+60=0.
\end{gathered}
\tag{3.89}
\]
The origin is the uncorrected face.

The weight-three component of the third moment is the third moment of
(3.85).  On the third and fourth ratios in (3.89) it is respectively
\[
\begin{aligned}
-1119744000y(42x+17y^2),\\
-373248y(4x-27y^2),
\end{aligned}
\tag{3.90}
\]
so neither survives.  On the quadratic pair, reduction modulo
\(t^2+12t+60\) gives
\[
124416y\left((51t+1080)x-(259t+2220)y^2\right),
\tag{3.91}
\]
which is nonzero at both roots.  The remaining ratios are one-sided:
\[
\begin{array}{c|c}
(0,0,1)&
(\partial_x^2,\ y^5+xy^3),\\
(0,1,0)&
(\partial_x^2+\partial_x\partial_y^2,\ y^5),
\end{array}
\tag{3.92}
\]
and the origin gives \((\partial_x^2,y^5)\).  Every term off the displayed
faces has strict weight.  The \(x\)-derivative count on each one-sided
face therefore gives the eventual mixed cutoff and proves the theorem.

Together with Proposition 3.13 and Theorems 3.14--3.15, this closes every
\(r=2,\deg P=5\) pair and proves Corollary 1.5.

The following proposition gives the finite leading classification for the
cubic-leading degree-five row.  Its eight nonhomogeneous correction systems
were subsequently closed by the
[binary degree-five frontier theorem](BINARY_DEGREE_FIVE_GVC_FRONTIER.md).

> **Proposition 3.17 — cubic-leading quintic top forms.**  Suppose
> \(r=3\), \(\deg P=5\), and the first two pure equations hold.  After
> scalar extension, linear changes, and rescaling, the leading pair is
> one of the following eight forms:
> \[
> \begin{array}{c|c}
> \Lambda_3&P_5\\ \hline
> \partial_x^3&
> x^2y^3,\ x(x-y)y^3,\ xy^4,\ y^5\\
> \partial_x^2\partial_y&
> x^5,\ xy^4,\ y^5\\
> \partial_x\partial_y(\partial_x+\partial_y)&
> x^5.
> \end{array}
> \tag{3.93}
> \]

For the triple-root orbit, the first equation gives
\[
 P_5=y^3(Ax^2+Bxy+Cy^2).
\tag{3.94}
\]
The affine stabilizer of \(\partial_x^3\) has four orbits on the binary
quadratic: two distinct finite roots, one double finite root, one finite
root together with infinity, or the double root at infinity.  These give
the first row of (3.93).

For the double-root orbit, normalize
\(\Lambda_3=\partial_x^2\partial_y\).  The first equation gives
\[
 P_5=Ax^5+Bxy^4+Cy^5,
\]
and the second leading moment is
\[
 960Ax y^2(9Bx+5Cy).
\tag{3.95}
\]
Thus either \(A\ne0\) and \(B=C=0\), or \(A=0\); the stabilizer reduces
the latter case to \(xy^4\) or \(y^5\).

Finally normalize the squarefree orbit to
\(\partial_x\partial_y(\partial_x+\partial_y)\).  Solving the first
equation writes
\[
P_5=A x^5-Dx^4y+2Dx^3y^2-2Dx^2y^3+Dxy^4+C y^5.
\tag{3.96}
\]
The coefficients of the second leading moment have Gröbner basis
\[
 D(A+C),\qquad C(5A-D),\qquad D(D+5C).
\tag{3.97}
\]
Its projective zero set is \(x^5,y^5,(x-y)^5\).  Root permutation makes
these one orbit, giving the last row of (3.93).  Proposition 3.17 itself
is only a top-form reduction; the later
[degree-five frontier theorem](BINARY_DEGREE_FIVE_GVC_FRONTIER.md) closes
all eight nonhomogeneous correction systems.

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

Formal factorization does not evade this obstruction through an umbral
change of polynomial coordinates.

> **Theorem 5.2 — multiplicative umbral straightening no-go.**  Let
> \(T:k[x_1,\ldots,x_n]\to k[x_1,\ldots,x_n]\) be a \(k\)-algebra
> automorphism.  If
> \[
> T^{-1}D_vT=F(\partial)
> \tag{5.6}
> \]
> is a locally finite formal constant-coefficient operator, then
> \(F\) is a linear form.  In particular no multiplicative conjugacy can
> straighten a genuinely nonlinear delta operator
> \(D_v+q(\partial)\).

Indeed, conjugation by an algebra automorphism preserves the Leibniz rule,
so \(F(\partial)\) is a derivation.  Applying the Leibniz rule formally to
two exponentials gives
\[
 F(\xi+\eta)=F(\xi)+F(\eta).
\tag{5.7}
\]
Over characteristic zero, an additive formal power series is linear.
Thus formal Hensel factorization of a higher-order symbol into smooth
branches is useful at the operator level, but its usual linear umbral
intertwiner cannot preserve products and therefore cannot transport the
powers \(P^m\).  Any successful use of those branches needs a new
multilinear identity rather than a multiplicative coordinate change.

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
   Formal or umbral straightening does not supply a multiplicative escape:
   Theorem 5.2 forces every such conjugated delta operator to remain
   first-order linear.
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
   filtration depth.  The later
   [uniform face-termination theorem](BINARY_GVC_UNIFORM_FACE_TERMINATION.md)
   classifies the leading Hall locus for every \(r<\deg P\) by one root
   multiplicity and proves that every unequal-weight common-threshold
   face is automatically terminal.  Its density refinement also kills
   every path with only \(o(m)\) excursions from one unequal-weight or
   ordinary-homogeneous face.  Thus a remaining path must have depth
   proportional to \(m\), avoid every single common threshold, and use at
   least two distinct faces with positive limiting density in the
   coefficient at \(\lfloor\rho m\rfloor\).  The missing theorem must turn
   that convolution into a one-radial exposed lower face or reduce its
   Hall/jet support.  Generic translation makes the associated balanced
   radial polytope exactly the intersection of the operator Newton polygon
   with the downward derivative polygon of \(P\).  Multiradial prime
   isolation closes an empty intersection or one with a componentwise
   least point, so the unresolved convolution must join incomparable
   Pareto-minimal endpoints.  A further degree row alone does not settle
   this asymptotic alternative.
2. **Lowest-order-two anchor.**  Theorem 3.4 removes the entire
   lowest-order-one class.  The first unresolved symbols therefore have
   lowest positive order \(r\ge2\), are genuinely nonhomogeneous, have
   \(\deg P>r\), and do not admit the split factor-unit architecture.
   Theorems 3.7 and 3.9 close \(r=2,\deg P\le4\), while Theorems
   3.10--3.12 close every \(r=3,\deg P=4\) leading-symbol orbit.
   Corollary 1.4 therefore moves the universal binary frontier to
   \(\deg P\ge5\).  In degree five only \(r=2,3,4\) can occur.
   Proposition 3.13 reduces the first \(r=2\) row to three top forms;
   Theorems 3.14--3.16 close all three forms.  Thus the first
   quadratic-leading target has degree at least six, while a degree-five
   counterexample must have \(r=3\) or \(4\).  Proposition 3.17 reduces
   the \(r=3\) row to eight top forms, and the
   [degree-five frontier theorem](BINARY_DEGREE_FIVE_GVC_FRONTIER.md)
   closes all eight correction systems.  It also closes the \(r=4\)
   squarefree-quartic row uniformly in its cross-ratio.  The
   [quadruple-root quartic theorem](BINARY_QUARTIC_QUADRUPLE_ROOT_GVC.md)
   now closes the partition \((4)\), with arbitrary lower symbol terms
   and arbitrary higher operator jets.  The
   [triple-plus-simple theorem](BINARY_QUARTIC_TRIPLE_SIMPLE_ROOT_GVC.md)
   and [double-root theorem](BINARY_QUARTIC_DOUBLE_ROOT_GVC.md)
   close the other three quartic partitions.  Therefore every binary
   operator satisfies GVC through polynomial degree five, and the next
   nonhomogeneous degree frontier begins at six.  The
   [complete quintic-leading sextic theorem](BINARY_QUINTIC_ALL_ROOT_PARTITIONS_GVC.md)
   closes the entire \((r,\deg P)=(5,6)\) row, including every quintic
   root partition and arbitrary higher jets.  The later quartic- and
   cubic-leading sextic theorems close \(r=4\) and \(r=3\), respectively.
   The
   [complete quadratic-leading sextic theorem](BINARY_QUADRATIC_ALL_ROOT_PARTITIONS_GVC.md)
   closes \(r=2\), including both nested axes of the pure-sixth-power
   endpoint.  Hence every binary operator satisfies GVC through polynomial
   degree six, and the next nonhomogeneous degree frontier begins at seven.
   Theorem 3.1 has already removed every fixed leading face.
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
