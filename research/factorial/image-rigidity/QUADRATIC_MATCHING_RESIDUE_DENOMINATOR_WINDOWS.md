# Residue denominator windows for arbitrary quadratic matching weights

Consider
\[
 f=c+\sum_{i=1}^m(\alpha_i X_i^2/2+d_iX_i)
                   +\sum_{i<j}\beta_{ij} X_iX_j,
 \qquad P_j=\Lambda_m(f^j)/j!
\]
and let \(s\) be the number of nonzero monomials of \(f\).
Assume at least one nonconstant term is present, work at a
place over two with \(v(2)=1\), and scale the entire phase
so that the minimum valuation of the nonconstant natural
coefficients is zero. Assume \(v(c)\ge1\), allowing \(c=0\).
All coefficient residues below are taken in the residue field.

Form the graph on the variables with one edge \(ij\)
when \(\bar\beta_{ij}\ne0\). Define the two polynomials
\[
 D(T)=\prod_{i=1}^m(1+\bar\alpha_iT+\bar d_i^2T^2),\qquad
 N(T)=\sum_{M\ {\rm matching}}T^{|M|}
       \prod_{ij\in M}\bar\beta_{ij}
       \prod_{i\ {\rm unmatched}}(1+\bar d_iT).
 \tag{1}
\]

**Theorem.** If \(N/D\) has a nonconstant reduced denominator
of degree at most \(s\), then \(f\) satisfies its actual-support
Strong Factorial window, in every order. No relations among
the cross coefficients are required.

The criterion also gives a shorter window from the degrees.
Write \(N/D=Q+R/D_{\rm red}\) over the residue field with
\(\deg R<r=\deg D_{\rm red}\), \(r>0\). At every positive start
a window of length \(r+\max(0,\deg Q)\) suffices, with the
second term taken as zero for \(Q=0\).

Two sufficient tests for the required denominator bound
are given below. They allow the cross graph to have vertices
without either a diagonal or a linear unit.

## Proof of the rational-series criterion

The [universal matching identity](QUADRATIC_STRONG_FACTORIAL_WINDOWS_FROM_MATCHINGS.md#2-the-universal-matching-generating-function)
applies to \(f-c\). The formal series of \(P_j\) for \(f\)
is its product with \(\exp(cT)\). For every \(k\ge1\),
\[
                  v(c^k/k!)\ge k-v_2(k!)=s_2(k)>0,
\]
where \(s_2\) is binary digit sum. The product is integral
coefficientwise and has residue exactly \(N/D\).

A nonzero proper series \(R/D_{\rm red}\) with
\(D_{\rm red}(0)=1\) and positive degree \(r\) cannot have
\(r\) consecutive zero coefficients starting at any index
\(i\ge0\). The coefficient recurrence is valid from order
\(r\) on. Its leading denominator coefficient is nonzero,
so successive applications first propagate the block
back to indices \(0,\ldots,r-1\), then forward to all
indices. That would contradict \(R\ne0\).
When \(N/D\) is in reduced form and its denominator is
nonconstant, \(R\ne0\).

A positive-start block contains at most
\(\max(0,\deg Q)\) indices at which the polynomial part
can contribute. The stated window therefore has a
contiguous suffix of \(r\) indices with no polynomial
contribution. It cannot vanish. Finally each monomial
in a summand of \(N\) uses a distinct subset of the
unit cross and linear coefficients and has degree
equal to the number used. Thus \(\deg N\le s\).
After cancellation this also bounds the reduced numerator
degree. Now \(r+\max(0,\deg Q)\) is the maximum of \(r\)
and the reduced numerator degree, so \(r\le s\) gives
the actual-support assertion.

## A graph bound before general polynomial cancellation

Let \(q\) be the number of unit diagonals, \(\ell\) the
number of unit linears, and \(e_0\) the number of edges.
Let \(u\) count vertices with unit linear, nonunit
diagonal and at least one incident edge. At an isolated
vertex with unit linear and nonunit diagonal, (1) cancels
one factor \(1+\bar d_iT\), because its denominator factor
is \((1+\bar d_iT)^2\). Cancel these factors to obtain
\(N^\circ,D^\circ\). The denominator degree is
\[
                 \deg D^\circ=q+\ell+u .
 \tag{2}
\]
Consequently \(q+\ell+u\le s\) suffices for the denominator
bound. In particular this inequality holds if \(u\le e_0\),
since \(q+\ell+e_0\le s\).

The component condition that every nontrivial tree
component has a vertex without a unit linear or with a
unit diagonal implies \(u\le e_0\). A tree on \(v\ge2\)
vertices has \(v-1\) edges and at most \(v-1\) of the
vertices counted by \(u\) under this condition. A connected
component with a cycle has at least \(v\) edges and at
most \(v\) such vertices. Isolated vertices do not
contribute to \(u\). Other components or nonunit terms
can still supply enough support for (2) when this tree
condition fails.

For a simple sufficient condition that \(N/D\) is
nonpolynomial, let \(\nu_0\) be the matching number of
the subgraph induced by vertices without a unit linear,
and let \(\ell_0\) count all vertices with unit linear
and nonunit diagonal, including isolated ones. A matching
contributes numerator degree at most \(\ell+\nu_0\):
relative to the empty matching, an edge meeting no
unit linear can add one, one meeting exactly one adds
zero, and one meeting two subtracts one. Thus
\[
                    \nu_0< q+\ell_0
 \tag{3}
\]
makes \(\deg N<\deg D=q+\ell+\ell_0\).
As \(N(0)=1\), a nonzero proper rational function is
nonpolynomial. Conditions (2) and (3) are sufficient,
with no need for further polynomial cancellation.

For example, take any number of variables \(Y_i\) and
\(f=\alpha X^2/2+dX+\sum_i\beta_i XY_i\), with \(\alpha,d\)
units and the \(\beta_i\) arbitrary integral coefficients.
The residue denominator has degree two and the numerator
degree at most one. This gives a two-window at every
positive start, even though the cross edges may reach
arbitrarily many variables with neither a diagonal nor a
linear term. The same two-window allows \(c\) of valuation
at least one. This example includes independent cross
coefficients outside the fixed principal-minor relations.

## A connected linear tree at the critical denominator degree

Suppose all nonconstant monomials are unit terms and
consist of one linear at each vertex of a tree on \(v\ge2\)
vertices and one cross monomial for each of its edges.
Suppose \(c=0\), so \(s=2v-1\). Here (1) gives
\[
 D(T)=\prod_i(1+\bar d_i^2T^2)
       =D_2(T^2),\qquad \deg D_2=v .
\]
Write \(N(T)=N_{\rm even}(T^2)+T N_{\rm odd}(T^2)\).
Since \(\deg N\le v\), both numerator degrees are
strictly below \(v\). The even numerator has constant
coefficient one. The even series is therefore nonzero
and proper. The odd series is also nonzero and proper
when \(N_{\rm odd}\ne0\). A window of \(2v-1\) indices
contains \(v\) consecutive coefficients of the parity of
its start. It follows that the actual-support window holds
at every even start, and at every start when
\(N_{\rm odd}\ne0\). If a sufficiently divisible
nonzero constant is present, \(s=2v\) and the denominator
criterion alone suffices since the even series is nonzero
and proper.

The results are local sufficient conditions. An algebraic
coefficient profile needs to pass at only one place above
two. In particular, neither the graph degree bound nor a
nonpolynomial residue series is being asserted for every
quadratic profile.
