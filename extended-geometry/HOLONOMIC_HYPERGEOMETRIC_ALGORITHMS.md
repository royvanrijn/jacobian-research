# Holonomic and hypergeometric algorithms for contraction moments

## 1. Status and scope

This note records two external algorithmic frameworks that can be applied to
the repository's contraction formulas:

1. creative telescoping for finite hypergeometric multiple sums; and
2. holonomic systems attached to Laurent-polynomial support configurations.

The first framework gives a directly applicable route from coefficient
extraction to exact recurrences. The second gives a possible geometric rank
bound only after additional realization, nondegeneracy, and ordinary-point
hypotheses have been checked.

No new all-order contraction identity or finite moment cutoff is proved here.
In particular, holonomicity by itself does not show that a bounded zero prefix
forces an entire moment sequence to vanish.

## 2. Contraction coefficients are finite hypergeometric sums

Write

\[
 f(\zeta,z)=\sum_{j=1}^N c_j\zeta^{a_j}z^{b_j},
 \qquad a_j,b_j\in\mathbb N^r,
\]

and use the contraction

\[
 \mathcal E_r(\zeta^\alpha z^\beta)=
 \begin{cases}
 \displaystyle\frac{\beta!}{(\beta-\alpha)!}z^{\beta-\alpha},
     &\alpha\leq\beta\text{ componentwise},\\[6pt]
 0,&\text{otherwise}.
 \end{cases}
\]

For \(k=(k_1,\ldots,k_N)\in\mathbb N^N\), put

\[
 |k|=\sum_jk_j,\qquad
 A(k)=\sum_jk_ja_j,\qquad
 B(k)=\sum_jk_jb_j,\qquad
 c^k=\prod_jc_j^{k_j}.
\]

Multinomial expansion followed by contraction gives the exact formula

\[
 \boxed{
 [z^\rho]\mathcal E_r(f^m)
 =
 \sum_{\substack{k\in\mathbb N^N\\
                  |k|=m\\
                  B(k)-A(k)=\rho}}
 \frac{m!}{\prod_jk_j!}\,c^k\,
 \frac{(A(k)+\rho)!}{\rho!}.}
 \tag{2.1}
\]

The sum is finite. After solving the linear constraints for a lattice of
free indices, its summand is a product of constants and factorials of affine
linear forms in \(m\) and the free indices. Every admissible unit shift
therefore changes the summand by a rational function. In the standard
terminology, (2.1) is a proper hypergeometric multisum.

If

\[
 Q(\zeta,z)=\sum_{\ell=1}^Lq_\ell\zeta^{u_\ell}z^{v_\ell},
\]

then the analogous mixed coefficient is the finite outer sum over \(\ell\)
of

\[
 \boxed{
 \sum_{\substack{k\in\mathbb N^N\\
                  |k|=m\\
                  B(k)+v_\ell-A(k)-u_\ell=\rho}}
 q_\ell\frac{m!}{\prod_jk_j!}\,c^k\,
 \frac{(B(k)+v_\ell)!}{\rho!}.}
 \tag{2.2}
\]

Thus fixed multipliers preserve the hypergeometric form. In balanced
bihomogeneous problems only \(\rho=0\) occurs, so the contraction is already
a scalar finite multisum. Formula (2.1) is the sparse general form behind,
for example, the one-sum identities (4.4)--(4.5) in
[`TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md`](TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md).

There are two implementation details:

- choose independent lattice coordinates before telescoping, rather than
  retaining Kronecker constraints; and
- encode the finite support with factorial or binomial zero conventions, or
  keep boundary terms explicitly. A telescoping certificate is valid only
  after its boundary contributions have been proved to vanish.

## 3. The Paule--Schneider method

Paule and Schneider consider a term
\(F(m,s_1,\ldots,s_e)\) hypergeometric in every variable and a definite
sum

\[
 S(m)=\sum_{s_1}\cdots\sum_{s_e}F(m,s_1,\ldots,s_e)
 \tag{3.1}
\]

with finite summand support. Their target is a P-finite recurrence

\[
 p_\gamma(m)S(m+\gamma)+\cdots+p_0(m)S(m)=0,
 \qquad p_i\in K[m].
 \tag{3.2}
\]

Finite support guarantees the existence of some homogeneous recurrence by
the Wilf--Zeilberger/Fasenmyer theory. Their faster construction eliminates
the sums recursively:

1. obtain ordinary and hook-type recurrences for an inner sum;
2. use them as rewrite rules in a finite shift basis;
3. reduce the outer creative-telescoping equation to a parameterized linear
   recurrence for a rational certificate; and
4. telescope the outer index.

The main computational bottleneck is the rational solution of the
parameterized recurrence. The paper develops denominator bounds,
Gosper--Petkovšek numerator information, and related preprocessing for this
step.

The contiguous-relation theory explains a broad guaranteed input class.
For \({}_{q+1}F_q\)-type summands, Theorem 1A cited in their Section 3
guarantees telescoping contiguous relations among \(q+1\) distinct
nonnegative integral parameter shifts. It also bounds the degree of the
polynomial certificate. Zeilberger recurrences and the hook recurrences
needed by the recursive multiple-sum method are special cases.

This guarantee should not be overstated. The paper's displayed double-sum
algorithm can still return `Failure` at a prescribed recurrence order or
when its rational-solution ansatz has no solution. The safe conclusions are:

- a finite proper-hypergeometric sum is P-finite;
- the contiguous theory guarantees the inner hook relations for the stated
  \({}_{q+1}F_q\) shift class; and
- the proposed recursive method is efficient on many examples, but is not
  presented as a universal success theorem for every representation and
  every prescribed order.

## 4. What constitutes an all-order certificate

A recurrence found from sample values is evidence, not a certificate. A
creative-telescoping proof should store an identity of the form

\[
 \sum_{i=0}^{\gamma}p_i(m)F(m+i,k)
 =
 \sum_{j=1}^e\Delta_{k_j}G_j(m,k),
 \tag{4.1}
\]

or the nested analogue used by the recursive method. After division by a
common hypergeometric term, (4.1) is usually a rational-function identity
and can be checked independently. One must then verify the finite-support
boundary cancellation.

The recurrence plus initial values proves the claimed sequence identity.
If \(p_\gamma(m)\ne0\) for every integer \(m\geq m_0\), then
\(\gamma\) consecutive values beginning at \(m_0\) determine the tail.
Nonnegative integer zeros of the leading coefficient create singular
steps and require separate bridge values or a different recurrence.
Accordingly, the recurrence order alone is not always the number of
initial values required.

For repository use, a recurrence certificate should record:

- the exact summand and its support convention;
- the recurrence and all telescoping certificates;
- the exceptional integer zeros of the leading coefficient;
- a primitive integral normalization, including any factorial or Borel
  normalization used to pass between the period and the raw moment;
- the numerator/denominator or companion-matrix singularity ledger needed
  for reduction modulo \(p^a\);
- the initial or bridge values used for propagation; and
- an independent exact checker of (4.1), not only a long numerical replay.

### 4.1 Arithmetic postprocessing after telescoping

Once a recurrence is certified, reduction modulo \(p\) is a second
algorithmic stage, not a formal afterthought.  The recurrence must first be
made primitive over an integral coefficient ring.  Scalar multiplication
by factorials can be harmless in characteristic zero while changing every
positive-characteristic singular step.

For an order-one recurrence
\[
 A(m)S(m+1)=B(m)S(m),\qquad
 A,B\in\mathbb Z[m],\quad \gcd(A,B)=1,                  \tag{4.2}
\]
factor \(A\) and \(B\) before reduction.  Then
\[
 v_p(S(m+1))-v_p(S(m))
 =v_p(B(m))-v_p(A(m))                                  \tag{4.3}
\]
is an exact local transition rule whenever \(S(m)\ne0\).  It turns the
recurrence into a base-\(p\) valuation automaton without expanding the
moments themselves.

The recurrence \(p\)-curvature is the \(p\)-step shift multiplier
\[
 \Psi_p(m)=\prod_{i=0}^{p-1}\frac{B(m+i)}{A(m+i)}.       \tag{4.4}
\]
It is useful for checking the reduced difference module and its true
singularity class.  It is not, by itself, a valuation classifier: shift
norms can cancel numerator zeros against denominator poles that occur at
different steps of the same residue orbit.

For an order-\(\gamma\) recurrence, put it in companion form
\[
 Y(m+1)=C(m)Y(m).                                      \tag{4.5}
\]
The corresponding curvature is
\[
 C(m+p-1)\cdots C(m),                                  \tag{4.6}
\]
but prime-power propagation must retain more data: an integral companion
lattice, the Smith valuations of the individual step matrices, and bridge
vectors across singular steps.  Taking the generic characteristic
polynomial of (4.6) first can erase precisely this information.

The SIC2C4 radial family is the calibration.  Its good-prime recurrence
curvature depends only on the degree \(d\), while two degree-eight
recurrences with the same curvature have respectively re-entrant and
monotone valuation phases.  See
[`TWO_PAIR_SIC_FROBENIUS_CURVATURE_BRIDGE.md`](TWO_PAIR_SIC_FROBENIUS_CURVATURE_BRIDGE.md).
Accordingly, the repository pipeline after creative telescoping is:

1. clear parameter and sequence-normalization denominators primitively;
2. factor scalar leading/trailing coefficients and locate singular steps;
3. build the local valuation or Smith ledger over \(\mathbb Z_p\);
4. compute \(p\)-curvature as a module-level consistency check; and
5. certify bridge data rather than canceling singular factors generically.

## 5. The stable GKZ rank theorem

Let \(A\) be an \(n\times N\) integer matrix of rank \(n\), with columns
\(w_1,\ldots,w_N\), and consider the universal Laurent polynomial

\[
 F(t,x)=\sum_{j=1}^Nx_jt^{w_j},
 \qquad
 \Delta_\infty=\operatorname{conv}(0,w_1,\ldots,w_N)
 \subset\mathbb R^n.
 \tag{5.1}
\]

Lei Fu defines the stable GKZ hypergeometric
\(\mathcal D_{\mathbb A^N}\)-module by a direct-image construction from the
exponential connection of \(F\) twisted by a Kummer connection. The main
theorem has two parts:

1. the stable GKZ module is holonomic; and
2. on a Zariski-open set \(U\) parametrizing Laurent polynomials
   nondegenerate with respect to \(\Delta_\infty\), it is an integrable
   connection of rank
   \[
   \boxed{R_A=n!\operatorname{vol}(\Delta_\infty).}
   \tag{5.2}
   \]

Here nondegeneracy means that, for every face \(\Gamma\) of
\(\Delta_\infty\) not containing the origin, the face polynomial
\(F_\Gamma\) has no common critical point in \((\mathbb C^*)^n\):

\[
 \partial_{t_1}F_\Gamma=\cdots=\partial_{t_n}F_\Gamma=0.
 \tag{5.3}
\]

The stable construction has the expected rank for every Kummer parameter.
This should not be conflated with the classical
\(A\)-hypergeometric-module rank statement, where resonance and the lattice
generated by the columns require separate treatment.

## 6. When rank can bound initial data

The Newton-polytope number \(R_A\) is the rank of a connection on the
generic coefficient open set. It is not, without further work, the order of
a recurrence in \(m\).

A valid route from (5.2) to a moment bound would be:

1. package the moments into a generating function, for example
   \[
   \Phi(s)=\sum_{m\geq0}\mu_m\frac{s^m}{m!};
   \]
2. identify \(\Phi\) as a scalar component of the pullback of the stable
   GKZ connection along an explicit coefficient curve \(x=x(s)\);
3. prove that the curve stays in the nondegenerate open set near the
   expansion point and that the point is ordinary for the pulled-back
   connection; and
4. derive a scalar cyclic-vector equation and audit its singularities.

Under those hypotheses, the pullback has rank \(R_A\), and a scalar
component satisfies a local differential equation of order at most \(R_A\).
At an ordinary point, at most \(R_A\) Taylor data determine that scalar
solution. Since \(\Phi^{(m)}(0)=\mu_m\), this would give the proposed bound
on initial moments.

Three obstructions prevent using (5.2) as an automatic cutoff:

- the coefficient-extraction or factorial functional must first be
  realized in the stable GKZ direct-image family;
- the usual scaling curve \(x(s)=s\,c\) passes through the zero polynomial
  at \(s=0\), which is outside the nondegenerate open set; and
- a scalar recurrence obtained at a singular expansion point can require
  more bridge data than the connection rank suggests.

Thus (5.2) is a candidate *generic local* bound. It is not currently a
repository theorem about the first \(R_A\) contraction moments.

## 7. Audit of the two-pair chart

The coefficient-extraction chart (4.3) of the two-pair counterexample uses

\[
 P(x,v)=\frac{1+x}{2x}\left(1-v^2(1+x)^2\right).
 \tag{7.1}
\]

Its exponent support in the Laurent variables \((x,v)\) is

\[
 \{(-1,0),(0,0),(-1,2),(0,2),(1,2),(2,2)\}.
 \tag{7.2}
\]

The Newton polygon at infinity has vertices

\[
 (-1,0),\ (0,0),\ (2,2),\ (-1,2),
\]

Euclidean area \(4\), and normalized volume \(2!\cdot4=8\). Therefore the
generic stable GKZ family with this support has rank \(8\) on its
nondegenerate coefficient open set.

The counterexample itself is not on that open set. Its top face polynomial
is

\[
 P_{\mathrm{top}}(x,v)
 =-\frac12v^2x^{-1}(1+x)^3.
 \tag{7.3}
\]

At \(x=-1\) both partial derivatives vanish in the torus, because the
factor \((1+x)^3\) has a multiple root. Hence the specialized Laurent
polynomial is degenerate with respect to its Newton polygon.

This example shows both the promise and the limitation of the rank formula.
The support gives the small generic number \(8\), but it does not bound the
special moment sequence directly. In fact the mixed moment in the existing
proof satisfies the first-order recurrence

\[
 (2m+1)B_m=2mB_{m-1},
\]

which is much smaller than the generic rank and was obtained by direct
finite telescoping.

<!-- status-consumer: FTI1 836c1443b2e29cd8 -->

## 8. Recommended workflow

For a new parameterized contraction family:

1. derive (2.1) or (2.2) symbolically and reduce its constraint lattice;
2. try direct finite differences, Gosper telescoping, or
   Paule--Schneider-style recursive telescoping;
3. retain the exact certificate and singular-step audit;
4. primitively normalize the recurrence over the relevant integral
   parameter ring and separate period normalization from raw factorial
   normalization;
5. factor any certified first-order hypergeometric components and group their
   shift quotients by the exact divisor-orbit signature of
   [`FACTORIAL_TRACE_INDEPENDENCE.md`](FACTORIAL_TRACE_INDEPENDENCE.md);
   distinct signatures cannot cancel over exponential-rational coefficients,
   while equal signatures must be merged before initial-value tests;
6. compute the local factor/Smith valuation ledger before taking shift
   norms, and use recurrence \(p\)-curvature only as a module-level check;
7. independently compute the Laurent support and normalized Newton volume;
8. test face nondegeneracy at the actual coefficient point, not only
   generically; and
9. use the GKZ rank as an initial-data bound only after constructing the
   generating-function pullback and proving the expansion point ordinary.

The hypergeometric route is the primary proof-producing method for the
present contraction formulas. The GKZ route is a geometric explanation and
potential complexity bound whose application gates remain open.

## References

- P. Paule and C. Schneider,
  [*Creative Telescoping for Hypergeometric Double Sums*](https://arxiv.org/abs/2401.16314),
  arXiv:2401.16314 (2024). See especially Sections 1, 3, 4, and 6.
- L. Fu,
  [*The stable GKZ hypergeometric \(\mathcal D\)-module*](https://arxiv.org/abs/2602.16941v2),
  arXiv:2602.16941v2 (2026), main theorem.
- Y. Zhou and M. van Hoeij,
  [*Desingularization and \(p\)-Curvature of Recurrence Operators*](https://arxiv.org/abs/2202.08931),
  arXiv:2202.08931 (2022), especially Sections 2--3.
