# Real windows and support depths on polynomial repair carriers

Polynomial models for local acceptance can retain the three
radical caps and the correct target repair height while producing
false additive incidences. A short real window imposes another
constraint. The degree of the additive error is bounded below
in terms of the formal radical, with a further contribution from
its contact orders on the original support. Integer line slicing
turns these constraints into uniform aggregate counts for moving
polynomial carriers of bounded degree.

## A degree bound with contact on the original support

Let \(A,B,C\) be pairwise coprime nonconstant polynomials over
\(\mathbb Q\), in any finite number of variables. Write
\[
 E=A+B-C,\qquad D=\max(\deg A,\deg B,\deg C),\qquad
 r=\deg\operatorname{rad}(ABC).
\]
Assume \(r<D\). Then \(E\ne0\) and
\[
                         \deg E\ge D-2r+3.
 \tag{RW1}
\]
If every irreducible factor of each original polynomial
has multiplicity at least two, the stronger bound is
\[
 \deg E\ \ge D-2r+3+\sum_{P\mid ABC}\bigl(v_P(E)-1\bigr)_+\,\deg P .
 \tag{RW2}
\]
In particular these carriers have a factor \(P\mid ABC\)
with \(v_P(E)\le2\).

A linear dependence among \(A,B,C\) would involve all three,
because they are nonconstant and pairwise coprime. In that
case choose a nonzero two-column first-derivative determinant.
It is unchanged up to a nonzero scalar when one uses another
pair of columns. Thus it is divisible at each \(P\mid ABC\)
by \(P^{v_P(ABC)-1}\). Using the two columns other than a
largest-degree one then gives \(D\le r-1\). Hence the three
polynomials are linearly independent. Theorem 4.1
of de Bondt, applied to \(A+B-C-E=0\) with the error term last
and dimension three, already gives (RW1). The following
determinant calculation retains the additional contact term.
[de Bondt, Theorems 1.2 and 4.1](https://arxiv.org/pdf/0707.0434)

There is a nonzero three-column determinant \(W\) of the three
polynomials whose rows consist of the original functions and
partial derivatives of orders \(0,1,1\) or \(0,1,2\). For an
elementary verification put \(u=B/A,v=C/A\). If two entries of
the gradient give a nonzero minor of \((\nabla u,\nabla v)\),
use these two first derivatives. Otherwise \(\nabla v=h\nabla u\)
for a rational function \(h\), and some \(u_i\ne0\).
The function \(h\) is nonconstant, since constant \(h\) makes
\(v-hu\) constant and contradicts independence. Choose \(j\)
with \(h_j\ne0\) and take derivative rows \(\partial_i\)
and \(\partial_j\partial_i\). Dividing out \(A\) from the
three columns after row operations leaves the nonzero
two-by-two minor \(u_i v_{ij}-v_i u_{ij}=u_i^2 h_j\).
The additional terms involving derivatives of \(A\) drop out
against the zeroth row or against the first derivative row,
using proportionality of the gradients.

Write \(s\in\{2,3\}\) for the sum of the derivative orders
and \(b\in\{1,2\}\) for their maximum, so \(s-b=1\).
Replacing the largest-degree column by \(E\) preserves
the determinant up to sign. Its degree is therefore at most
\[
             \deg A+\deg B+\deg C-D+\deg E-s.
\]
For a factor \(P\) of an original column \(T\), put
\(a=v_P(T)\), \(w=v_P(E)\). The original determinant is
divisible by \(P^{(a-b)_+}\). Use instead the equivalent
columns consisting of \(T,E\) and a second original column.
A determinant term differentiates the two factors by
a total of at most \(s\), so it is divisible by
\(P^{(a+w-s)_+}\). If \(a\ge2\), the combined exponent
lower bound is at least \(a-b+(w-1)_+\). Summing these
bounds over distinct \(P\) yields
\[
 \deg E \ge D-br+s+
              \sum_{P\mid ABC}(v_P(E)-1)_+\,\deg P.
\]
The expression \(D-br+s\) is at least \(D-2r+3\).
Dropping the additional term and using \(a-b\) as a lower
bound even when \(a<b\) also proves (RW1) without the
multiplicity condition. If every \(v_P(E)\ge3\), the sum in
(RW2) is at least \(2r\), which contradicts \(\deg E\le D\).
This contact bound differs from the one using
\(\min(v_P(T),v_P(E))\), since it retains larger contacts
of the discrepancy above an original valuation.
[Formal contact bounds (PC1)–(PC2)](ABC_MINIMAL_REPAIR_POLYNOMIAL_THIRD_LAYER.md#a-bound-on-three-simultaneous-contact-orders)

## Aggregate counts with a real window and local depth

Fix the number \(m\ge1\) of variables, a degree cap \(D_0\),
and coefficient height exponent \(B_0\ge0\). Allow the carriers
\(A,B,C\in\mathbb Z[x_1,\ldots,x_m]\) to vary with \(H\), with
degrees at most \(D_0\) and coefficient magnitudes at most
\(H^{B_0}\), as above. Positions are integer points in
\([-H,H]^m\). Assume \(d_E=\deg E>0\), as holds by (RW1) in
the small formal radical regime \(r/D<1/2\).
For every \(0\le\omega<d_E\),
\[
 \#\{\mathbf x: |E(\mathbf x)|\le H^\omega\}
                     \ll_{m,D_0} H^{m-1+\omega/d_E}.
 \tag{RW3}
\]
In particular if \(r/D\le\lambda\) and the allowed error
is \(H^{\theta D}\), the bound has exponent at most
\(m-1+\theta/(1-2\lambda)\) when \(0<\lambda<1/2\)
and \(0\le\theta<1-2\lambda\). This formulation describes
the real window directly in the height scale \(H^D\).
If the coefficient height is relevant to the chosen height
scale, one may instead use \(\omega=\theta(D+B_0)\) whenever
\(\omega<d_E\).

For (RW3), choose a bounded primitive integer direction
where the top homogeneous form of \(E\) is nonzero.
Such a direction exists on a grid of bounded side length
depending only on \(m,D_0\), by the degree bound for zeros
on Cartesian products. Hence \(O_{m,D_0}(H^{m-1})\)
integer lines parallel to that vector meet the box.
On each line the univariate \(E\) has degree \(d_E\) and
nonzero integral leading coefficient. If its value at a
real parameter has absolute value at most \(H^\omega\),
some complex root is within distance \(H^{\omega/d_E}\)
by factoring over the complex numbers. There are
\(O_{D_0}(1+H^{\omega/d_E})\) integer parameters in the
union of these root neighborhoods. This also includes
the exact zeros.

Now fix a prime pool \(S\) with prime product at most
\(H^{C_0}\) for a fixed exponent \(C_0\). Impose the actual
integer tests
\[
 \begin{split}
  &A(\mathbf x)B(\mathbf x)C(\mathbf x)\ne0,\qquad
      0<|E(\mathbf x)|\le H^\omega,\
  &\operatorname{rad}_{\notin S}
       (A(\mathbf x)B(\mathbf x)C(\mathbf x))^k\mid E(\mathbf x),
       \qquad k\ge3 .
 \end{split}
\]
For each primitive irreducible \(P\mid ABC\) with
\(j=v_P(E)<k\), define
\[
 b_P=
 \begin{cases}
  0,&j=0,\[3pt]
  \displaystyle
    \frac{\min\{\deg P+B_0,\omega/j\}}{\lceil k/j\rceil},&j>0 .
 \end{cases}
\]
When at least one such factor exists, the number of accepted
positions is
\[
 \ll_{m,D_0,B_0,C_0,\varepsilon}
      H^{\,m-1+\min\{\omega/d_E,\ \min_{P:\,v_P(E)<k}b_P\}
            +\varepsilon}.
 \tag{RW4}
\]
The factor is guaranteed by (RW2) under its multiplicity
condition, and \(k\) may be chosen at a fixed budget before
the height limit. The estimate also applies to the subset
obeying primitivity, positivity, the three evaluated
radical caps and the least-repair cutoff.

To prove the local contribution, fix one such \(P\) and
write \(E=P^j E_1\) with \(E_1\) integral and coprime to \(P\).
Primitive factor coefficient magnitudes are \(O_{m,D_0}(H^{B_0})\).
Indeed multiplication of polynomials of bounded degrees has
coefficient norm bounded below by a fixed positive constant
times the product of the coefficient norms, by compactness
after normalizing the two norms. The nonzero complementary
integral factors have coefficient norm at least one.

Choose the slicing direction to preserve the top forms of
both \(P\) and \(E\). Complete the primitive direction to
a bounded unimodular coordinate change, possible uniformly
since only finitely many directions are used. In line
coordinates \((\mathbf y,t)\), \(P\) depends on \(t\), is
irreducible over \(\mathbb Q(\mathbf y)\), and its resultant
in \(t\) with \(E_1\) is a nonzero polynomial in \(\mathbf y\).
This follows from coprimality and Gauss's lemma, and also
holds if \(E_1\) is independent of \(t\). Its zero set
contributes \(O_{m,D_0}(H^{m-1})\) positions by counting
anchors and then line parameters. For \(m=1\) the nonzero
resultant is a constant and there are no exceptional
anchors.

At each other anchor include its resultant primes and
\(S\) in a fixed pool of polynomial prime product. For a
prime dividing \(P(\mathbf y,t)\) outside the pool,
\(E_1(\mathbf y,t)\) is a unit and the local test enforces
\(j v_p(P(\mathbf y,t))\ge k\). For \(j=0\) this confines the
whole support of the evaluated \(P\) to the fixed pool.
For \(j>0\) its outside part is \(\lceil k/j\rceil\)-full.
Moreover \(P(\mathbf y,t)^j\) divides the nonzero integer
\(E(\mathbf y,t)\), giving
\(\lvert P(\mathbf y,t)\rvert\ll H^{\min\{\deg P+B_0,\omega/j\}}\). Enumerate its
outside radical at cost \(O(H^{b_P})\) and then its
complete supported values at uniform height-subpower
cost. Each value has at most \(\deg P\) line parameters.
Sum over the lines and take the minimum with (RW3).
[Uniform support enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

## Application and scope

For a fixed polynomial model with formal radical degree
\(o(D)\) along the outer small-budget limit, the short
real test costs nearly one full parameter dimension in
a dense lattice box. The simultaneous-depth estimate can
improve that bound when the degree and height of a deficient
factor are small relative to its forced fullness order.
The primary degree estimate without the contact term is
classical. The counts here use fixed-degree,
coefficient-uniform line slices and actual divisibility of
integer evaluations. They concern false incidences in the
short real and local tests, with any true polynomial zeros
already confined to \(O(H^{m-1})\) positions by (RW3).

This excludes dense polynomial-carrier realizations of the
pointwise obstruction unless the formal radical-degree
hypothesis is lost. A set concentrated on thin parameter
loci, or a carrier whose evaluations have much smaller
radical than its formal radical degree, escapes this count.
The bounded height and the actual radical caps have yet to
force a carrier of the required kind.
[Pointwise models and real-filtered tests](ABC_CAPPED_SIGNED_POPULARITY.md)

For fixed homogeneous binary forms at primitive parameters, a
separate projective multiplicity and rational-spacing argument
can give a stronger real-window count. The bounds here apply
to arbitrary-dimensional boxes, and keep the local support-depth
condition as well as the real one.
[Homogeneous primitive real-window count (PR3)](ABC_POLYNOMIAL_REAL_FILTER_OBSTRUCTIONS.md#primitive-specializations-in-a-near-equality-window)
