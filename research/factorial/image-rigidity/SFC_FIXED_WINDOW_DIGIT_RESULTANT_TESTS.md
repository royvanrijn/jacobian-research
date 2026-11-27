# Fixed-window factorial resultant tests on unbounded exponent profiles

Fix the dimension \(d\), support size \(s\), and a set \(I\) of
\(s\) distinct positive orders. The exponent entries of the
\(s\) distinct columns can vary without bound. Fix a prime
\(p\) and a positive precision \(r\).

For each column \(b_j\), absorb its single factorial moment
\(\prod_i b_{ij}!\) into its coefficient. The new order-\(m\)
moment form has integer coefficients
\[
 B_{m,\alpha}=\binom m\alpha
       \prod_i \frac{\bigl(\sum_j \alpha_j b_{ij}\bigr)!}
                     {\prod_j(b_{ij}!)^{\alpha_j}},
 \qquad |\alpha|=m .
 \tag{1}
\]
Each row quotient is a multinomial coefficient, with a repeated
entry for each copy of a column in the product. The rescaling
preserves the common projective zero problem at all orders.

Let \(P\) be a nonzero integer polynomial in these coefficient
slots, written as a sum of distinct monomials, with coefficients
combined first. For a profile, let \(v_{\min}(P)\) be the
least \(p\)-adic valuation of an individual monomial after
substitution of (1), including its integer scalar coefficient.
A profile survives the test if
\(p^{-v_{\min}(P)}P(B)\) is zero modulo \(p^r\).
The homogeneous resultant of the moment forms is one possible
\(P\). For a counterexample its value is zero.

**Theorem.** Given any finite list of coefficient polynomials
\(P\), one can decide whether any exponent profile survives
every indicated test. If there is a surviving profile, one
can produce one. One can also decide whether the surviving
set is finite and, if finite, enumerate it. Fixed linear
constraints on the entries, such as equal total degrees,
or any regular restriction on their simultaneous base-\(p\)
digits can be included.

When the chosen polynomials vanish at every common projective
zero, an empty surviving set proves the factorial window over
every characteristic-zero field. A surviving profile
gives residue information at the specified precision.

## 1. The finite digit graph

By Legendre's formula, after the factor
\(\binom m\alpha\) the valuation of (1) is
\[
 \frac{1}{p-1}\sum_i\left(
       \sum_j\alpha_j s_p(b_{ij})-
                   s_p\!\left(\sum_j\alpha_j b_{ij}\right)\right),
\]
with \(s_p\) the digit sum. Read all the exponent entries
from least to most significant, padding to a common length.
Each output sum uses a bounded carry \(c\) with transition
\(\sum_j\alpha_j d_j+c=d'+pc'\).
The digit sum difference telescopes to \((p-1)\sum c'\),
including the bounded tail that clears the carries.
Thus the valuation of every coefficient in (1) is
its constant multinomial valuation plus a sum of
integer weights on a finite digit graph.

The units modulo \(p^r\) belong to a finite state space.
For an integer \(N\), write \(U_p(N!)=p^{-v_p(N!)}N!\).
The elementary identity
\[
 U_p(N!)=J_p(N)U_p(\lfloor N/p\rfloor!),\qquad
 J_p(N)=\prod_{1\le j\le N,\ p\nmid j}j
\]
gives a digit machine for \(U_p(N!)\) modulo \(p^r\).
Indeed the factor \(J_p(N)\bmod p^r\) is periodic, with
period \(p^r\) times the order of the product of one
block of \(p^r\) units. Every block has the same product.
Read the prefixes of \(N\) from most to least significant
and accumulate their \(J_p\) factors, maintaining the
prefix modulo that period. Reverse this finite graph,
including each specified output state, to read the sum
digits and input digits least significant first alongside
the carry computations. The same construction is used in
the [relative valuation automata](FACTORIAL_RELATIVE_VALUATION_AUTOMATA.md#1-valuation-and-unit-data-on-one-graph).
The units of all coefficients (1) follow by multiplication
and inversion in the unit group. Distinctness of columns
and the profile restrictions can be checked in the product
of the digit graphs. In particular fixed linear equalities
and inequalities use bounded carry computations and
sign tests, and fixed congruences use finite residues.

Use canonical simultaneous input words, with at least one
nonzero entry in the last input digit except for the all-zero
array. Append a fixed number of zero digits sufficient to
clear every sum carry. This convention makes the length
of a canonical word unbounded precisely when the set of
accepted profiles is infinite.

## 2. Finite linear tests on the transition counts

For each terminal state the monomial units of every detector
\(P\) are fixed modulo \(p^r\), including the units of its
integer scalar coefficients. The monomial valuations
are affine integer functions of the transition counts.

For each detector, choose a monomial of minimum valuation.
For every other monomial prescribe its difference from
that minimum as one of \(0,1,\ldots,r-1\), or as at least \(r\).
These are linear constraints on the counts. At the specified
terminal state the normalized residue is now determined:
sum the monomial units times \(p^\delta\) for the specified
differences \(\delta<r\). Retain only the choices giving zero
for every detector. Thus there are finitely many cases of
linear constraints imposed on accepting paths in a finite
graph.

For completeness, prescribe the set of edges with positive
count. Check that this support is connected after forgetting
edge directions and deleting unused vertices, with the
start and terminal vertices included. Impose the usual
degree balances for a directed path from the start to the
terminal vertex. With nonnegative integer counts, these
conditions give an Euler trail in the directed multigraph
with the prescribed multiplicities, hence a path with exactly
those counts. One can obtain the trail by adding an edge from
the terminal to the start and splicing circuits, then deleting
that edge. Allow an empty path separately. This is an
elementary instance of the [Presburger description of regular
path counts](https://doi.org/10.1145/3649855).
Deciding the resulting integer linear conditions proves
emptiness and recovers a witness word when present.
Unbounded path length is also a Presburger question.
If lengths are bounded, a bound can be found by testing
the same conditions with increasing lower bounds on
the path length, and all shorter words can be tried.
This proves the theorem.

For the initial trinomial window there are two small
detectors. The first form is \(Z_0+Z_1+Z_2\).
After setting \(Z_2=-Z_0-Z_1\), write the next two forms
as \(aX^2+bXY+cY^2\) and
\(dX^3+eX^2Y+fXY^2+gY^3\). The real quadratic is positive
definite, since it is the integral of a square of a
nonzero real polynomial. If the forms share a complex
direction then the cubic is divisible by the quadratic,
by conjugation of the two quadratic roots. Equating the
mixed coefficients of this divisibility gives
\[
             ace-bcd-a^2g=0,\qquad acf-c^2d-abg=0. \tag{2}
\]
Conversely, since \(a,c\ne0\), (2) gives the divisibility.
Either expression can be expanded as a polynomial in the
slots (1) for the digit test.

## 3. A limitation of the raw Sylvester expansion

There is arbitrarily deep cancellation in the raw fixed-window
determinant test even on the homogeneous support
\[
                y^D,\quad xy^{D-1},\quad x^2y^{D-2},
                          \qquad D\ge2 .
\]
This persists after eliminating the first equation and removing
the first vanishing power of \(D\) from each restricted form.
It holds at every prime, and simultaneously at any finite
set of primes.

To see the underlying degeneration, consider a same-side
profile with fixed distinct offsets \(0=u_0<u_1<\cdots<u_{s-1}\).
Normalize by (1), then divide order \(m\) by
\((mD)!/(D!)^m\). The resulting coefficients are rational
functions of \(D\), given by
\[
 \binom m\alpha
      \frac{U!}{\prod_i (u_i!)^{\alpha_i}}\,
      \frac{\prod_i (D)_{\underline{u_i}}^{\alpha_i}}
           {(mD)_{\underline U}},\qquad U=\sum_i\alpha_i u_i .
 \tag{3}
\]
They are regular at \(D=0\): the denominator has a simple
zero if \(U>0\), and there is at least one factor \(D\)
in the numerator in that case. At \(D=0\) the form is
\(Z_0^{m-1}\sum_i Z_i\). The order-one equation is exactly
\(\sum_i Z_i=0\) at every \(D\).
Write \(H_u=\sum_{h=1}^u 1/h\).
On the order-one plane the coefficient of \(D\) in order
\(m\ge2\) is
\[
                    (m-1)Z_0^{m-1}\sum_{i\ge1}H_{u_i}Z_i .
 \tag{4}
\]
Indeed, terms with a single copy of a positive offset \(u_i\)
have coefficient
\(1+(m-1)D H_{u_i-1}+O(D^2)\).
The two-copy terms at order \(D\), summed using the two
ordered column indices, are
\[
 -(m-1)D Z_0^{m-2}
        \left(\sum_{i\ge1}Z_i/u_i\right)\left(\sum_{i\ge1}Z_i\right).
\]
This follows from the leading coefficient of
\((D)_{\underline u}=(-1)^{u-1}(u-1)!D+O(D^2)\)
for \(u>0\). The remaining terms have at least three
positive offsets and vanish to second order.
Substitute \(\sum_{i\ge1}Z_i=-Z_0\) to get (4).

For the specific profile \(0,1,2\), set
\(Z_0=X,Z_1=Y,Z_2=-X-Y\), and divide the restricted
second and third forms of (3) by \(D\).
In descending \(X\)-degree their coefficients at \(D=0\)
are \((-3/2,-1/2,0)\) and \((-3,-1,0,0)\).
The constant coefficients as polynomials in \(X\) at \(Y=1\)
start at \(D/6\) and \(-3D^2/20\), respectively.
In the degree-\(2,3\) Sylvester matrix the last column
contains only these two constants. If the term uses
the quadratic constant, the next-to-last column contributes
a further power of \(D\). If it uses the cubic constant,
at least two powers have already occurred. The latter
choice can be completed using only nonzero constant
entries in the other four columns. Thus the minimum
\(D\)-order of the individual determinant monomials is two.
The exact resultant after the two divisions by \(D\) is
\[
 \frac{2D^3(D-1)^2(37D^4-156D^3+197D^2-36D-50)}
 {(2D-3)^3(2D-1)^3(3D-5)^2(3D-4)^2(3D-2)^2(3D-1)^2}.
 \tag{5}
\]
This follows by inserting (3) in the five-row
Sylvester determinant. Its \(D\)-order at zero is three.
For a fixed prime and integer \(D\) with \(v_p(D)\) large,
the difference between the valuation in (5) and the
minimum valuation of the individual determinant terms
is \(v_p(D)+O_p(1)\). Choose \(D\) divisible by high powers
of every prime in a proposed finite test to make the
normalized determinant vanish at all its fixed precisions.

The profile itself has no common zero at any integer
\(D\ge2\): the quartic numerator factor in (5) has
coefficients \(37,140,149,64,10\) in descending powers of
\(D-2\), all positive. Thus this is cancellation in a
particular presentation of the resultant. Algebraic
simplification or a different choice of detector can
remove it. Formula (4) identifies a common leading
factor at the \(p\)-adic boundary for arbitrary same-side
profiles.
