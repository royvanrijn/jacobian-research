# Fixed factorial windows under independent common translations

Fix distinct offset columns \(b_j\in\mathbb Z_{\ge0}^{d}\), with
\(1\le j\le s\), and a set \(I\) of \(s\) distinct positive orders.
Add the same vector \(L=(L_1,\ldots,L_d)\) to every column. The
heights \(L_i\) are independent nonnegative integers. Write
\(d_j=\sum_i b_{ij}\).

**Theorem.** There is a computable \(H\), depending only on the
offsets and \(I\), such that there is no common nonzero coefficient
direction for the orders in \(I\) when every \(L_i>H\), in each of
the following cases.

- All offset total degrees are distinct.
- Exactly one pair of distinct columns has equal total degree,
  with all other degrees distinct, and \(I\) contains an even order.
- There are three bivariate columns of the same degree and \(I\)
  consists of three consecutive orders starting at any even
  order or at an odd start covered by the
  [odd-start Gaussian theorem](SFC_ODD_GAUSSIAN_LOCAL_WINDOW_TESTS.md).
  This includes every odd start through \(12001\) and every
  \(q=(2k+1)p^e\), with \(k\ge0,e\ge1\) and odd prime
  \(p>4k+2\), without a character or resultant restriction.

In two variables, for each fixed offset support and order set in
this theorem, one can decide by finitely many exact calculations
whether any translated support admits a common coefficient zero
over \(\mathbb C\). The test either returns all failing integer
height pairs or detects a fiber on which every height pair fails.
In the latter case a counterexample can be obtained at any
nonnegative integer point on the fiber. A positive decision for
nonvanishing holds over every characteristic-zero field.

This extends the
[common translation analysis](SFC_TERNARY_COMMON_TRANSLATION_CUBIC_OBSTRUCTION.md)
and the
[single doubled degree block proof](SFC_COMMON_TRANSLATION_DOUBLED_DEGREE_BLOCKS.md#a-single-doubled-block)
to independent heights. In particular the ratios of heights may
tend to zero or infinity along an escaping family.

## 1. Exact normalization and the case of distinct degrees

For \(L_i>0\), strip the scalar factorial \(\prod_i(mL_i)!\)
at order \(m\) and set \(C_j=c_j\prod_i L_i^{b_{ij}}\)
as a coefficient change independent of order. Put
\(\tau_i=1/L_i\). The normalized polynomial form is
\[
  F_m(\tau,C)=\sum_{|\alpha|=m}\binom m\alpha C^\alpha
            \prod_i\prod_{\nu=1}^{\sum_j b_{ij}\alpha_j}
                               (m+\nu\tau_i).
 \tag{1}
\]
For real \(L_i>0\), the same formula is the moment
\[
       \mathbb E\Bigl(\sum_j C_j m^{d_j}Y^{b_j}\Bigr)^m,
 \qquad Y_i=\frac{X_i}{m L_i},\quad
        X_i\sim\operatorname{Gamma}(m L_i+1,1)
 \tag{2}
\]
with independent inputs and unit rate. The Gamma moment identity
gives (2) directly from (1). The positive-real extension is used
to certify an effective bound.

At \(\tau=0\), (1) is \((\sum_j C_j m^{d_j})^m\).
If the degrees are distinct, the evaluations at the \(s\) positive
orders are invertible by Descartes' rule of signs. Equivalently,
each maximal square matrix on distinct degrees and equally many
distinct positive orders has nonzero determinant. This applies to
complex coefficients by taking a nonzero real or imaginary part.
The resultant of (1) is therefore nonzero at the origin, so it
stays nonzero when all reciprocals \(\tau_i\) are sufficiently
small.

## 2. A doubled degree at independent scales

Let \(b,b'\) be the unique pair of equal degree \(d_0\) and let
\(v=b'-b\ne0\). Difference estimates must retain the variance
of this particular pair. In fact the rows where \(v_i=0\) can
have much larger fluctuations than any row where the pair differs.

Put \(\epsilon=\max_i L_i^{-1/2}\) and
\(\sigma^2=\sum_i v_i^2/L_i\). Write \(\delta=Y-\mathbf1\).
The Gamma cumulants at fixed order \(m\) give, at each coordinate,
\(\mathbb E\delta_i=1/(mL_i)\) and
\[
 \kappa_k(\delta_i)=(k-1)!\bigl(1+(mL_i)^{-1}\bigr)(mL_i)^{1-k}
                      \qquad(k\ge2).
\]
For every fixed positive integer \(p\), the moment-cumulant
formula at a suitable even order yields
\(\|\delta_i\|_p=O(L_i^{-1/2})\).
All constants here may depend on the fixed support and orders.
Consequently
\[
 \|Y^{b_j}-1\|_m=O(\epsilon),\qquad
 \|Y^{b'}-Y^b-v\cdot\delta\|_m=O(\epsilon\sigma).
 \tag{3}
\]
For the second estimate, factor \(Y^{\min(b,b')}\) out of the
difference. The remaining difference uses only rows with \(v_i\ne0\).
In those rows \(L_i^{-1/2}\le\sigma\), so each term of degree at
least two has norm \(O(\epsilon\sigma)\) by Hölder's inequality.
Multiplying the linear term by the common factor changes it by
\(O(\epsilon\sigma)\) as well. Expanding the fixed monomials gives
these estimates using only finitely many fixed \(L^p\) norms.

Uniformly as all heights tend to infinity,
\((v\cdot\delta)/\sigma\) converges in moments to a centered real
Gaussian of variance \(1/m\). The mean is \(O(\epsilon)\), by
Cauchy–Schwarz, and the variance differs from \(1/m\) by
\(O(\epsilon^2)\). For \(k\ge3\) the normalized cumulants are
\(O(\epsilon^{k-2})\): use
\[
 \sum_i |v_i|^k L_i^{1-k}
       \le\epsilon^{k-2}\Bigl(\sum_i v_i^2/L_i\Bigr)^{k/2}.
\]
Thus no lower bound on the relative heights of the differing
rows is needed.

Group the coefficients into sums \(A_g\) by distinct total degrees
\(d_g\) and a single contrast \(B\) for the pair. Take \(b\)
as the representative of degree \(d_0\), so the pair contributes
\(m^{d_0}(A_0Y^b+B(Y^{b'}-Y^b))\) inside (2).
Let \(p_m=\sum_g A_g m^{d_g}\). The \(L^m\) norm of the difference
between the phase and its value \(p_m\) at \(\mathbf1\) is bounded
by \(C_m(\epsilon\|A\|+\sigma|B|)\).
If the moment vanishes, the binomial expansion and Hölder's
inequality show
\[
 |p_m|\le C'_m(\epsilon\|A\|+\sigma|B|).
\]
Indeed a fluctuation of norm less than
\((2^{1/m}-1)|p_m|\) cannot cancel \(p_m^m\) in the moment.
The evaluation matrix on degree groups has full column rank.
For sufficiently small \(\epsilon\), absorb the \(\epsilon\|A\|\)
term to get \(\|A\|=O(\sigma|B|)\).

Suppose there are common zeros along a sequence of real heights
with \(\min_i L_i\to\infty\). At these zeros \(B\ne0\); normalize
it to one and pass to a limit of \(A/\sigma\). Denote its entries
by \(\widehat A_g\). The divided phases converge in moments, by
(3), to a Gaussian linear phase with mean
\(\widehat p_m=\sum_g\widehat A_g m^{d_g}\) and variance
\(m^{2d_0}/m\). Hence the numbers
\(\widehat p_m/(\mathrm i m^{d_0})\) are roots of
\(\operatorname{He}_m(\sqrt m\,x)\).
They obey a rational linear relation with all coefficients nonzero,
using the one-dimensional left nullspace of
\((m^{d_g-d_0})_{m,g}\) and its nonzero maximal minors.
An even order contributes a nonzero node. The
[scaled Hermite isolation theorem](SFC_SCALED_HERMITE_ROOT_ISOLATION.md)
forbids this relation. This excludes the escaping sequence.

## 3. Three bivariate columns of equal degree

In this case write \(b_j=(\lambda_j,D-\lambda_j)\), with distinct
\(\lambda_j\). Choose a basis \(U_0,U_1,U_2\) of their monomial
span with leading Taylor parts at \(\mathbf1\) of respective
orders zero, one and two:
\[
               1,\qquad \delta_1-\delta_2,\qquad
                         (\delta_1-\delta_2)^2.
\]
The Vandermonde matrix in the \(\lambda_j\) constructs it from the
sum, first and second powers of the \(\lambda_j\). For example,
the quadratic element has coefficient sums against
\(1,\lambda_j,\lambda_j^2\) equal to \(0,0,2\); differentiating
the monomials then gives the stated quadratic part.

Now \(\sigma^2=1/L_1+1/L_2\) is comparable to \(\epsilon^2\).
Write the phase after removing \(m^D\) as \(S U_0+P U_1+Q U_2\).
Normalize \(\max(|S|,\sigma|P|,\sigma^2|Q|)=1\). Along an escaping
sequence, pass to a limit of the three normalized coefficients.
The cumulant estimates of Section 2 and Taylor remainders of the
fixed basis polynomials show convergence in moments to
\(A+B Z_m+C Z_m^2\), with \(Z_m\) centered Gaussian of variance
\(1/m\) and \((A,B,C)\ne0\).
If \(C=0\), a constant phase is excluded directly. A linear phase
would give the same normalized Hermite node at all orders,
which is impossible by isolation and the presence of an even
order. If \(C\ne0\), dividing by \(C\) gives the
[quadratic Gaussian gate at even starts](SFC_COLLINEAR_TERNARY_EVEN_START_AND_DYADIC_REDUCTION.md)
or the
[proved odd-start gates](SFC_ODD_GAUSSIAN_LOCAL_WINDOW_TESTS.md),
including the complete
[dominant initial-prime family](SFC_ODD_GAUSSIAN_LOCAL_WINDOW_TESTS.md#7-large-prime-powers-dividing-the-initial-order).
Those results exclude the common moment zero. This works for
independent heights since the Taylor and cumulant bounds hold
throughout the positive corner, with no fixed direction.

The same implication holds at any other specified start for which
the rational quadratic Gaussian gate has no solution. At a fixed
start the gate can itself be tested over the algebraic numbers.

## 4. Effective search in two independent heights

The proofs give a corner of positive real \(\tau_i\) on which
the rational polynomial forms (1) have no nonzero common zero.
A terminating search finds an integer \(H\) by [real closed
field quantifier elimination](https://doi.org/10.1016/j.apal.2017.03.001),
Theorem 1. For increasing \(H\), test whether
there is a point \(0<\tau_i<1/H\) and complex coefficients with
\(\sum_j|C_j|^2=1\) annihilating all the forms. Use real and
imaginary parts of those forms as polynomial equations.
In the distinct-degree case let \(G(\tau)=\sum_\alpha
g_\alpha\tau^\alpha\) be the resultant of (1). Then
\(g_0\ne0\), and one can simply take an integer
\(H>1+\sum_{\alpha\ne0}|g_\alpha/g_0|\).

In two variables all failures must therefore lie on the finitely
many fibers \(L_1=h\) or \(L_2=h\), \(0\le h\le H\). On such a
fiber strip each nonzero common factorial without using reciprocal
heights. Each coefficient is then a polynomial in the remaining
height \(L\), with integer coefficients. Compute the homogeneous
resultant \(R_h(L)\) of the \(s\) prescribed forms in the \(s\)
coefficient variables. If \(R_h\ne0\), isolate and test its
nonnegative integer zeros. If \(R_h=0\) identically, every ordinary
height on the fiber gives a common zero over
\(\overline{\mathbb Q}\). Evaluate the polynomial system at one
height and solve in the finitely many coefficient charts to
produce an algebraic point. This gives the claimed decision
procedure. For a consecutive window, any nonzero point found
gives a Strong Factorial counterexample even if some coordinates
vanish, because its actual support window is contained in the
prescribed window.

The decision is uniform in the independent heights for the
specified support and order set. Deciding the conjecture still
requires control of the offsets and the unbounded starting
order. The procedure has no evaluated general height cutoff.

## 5. An exact initial window for the translated support \(1,X,Y^2\)

**Theorem.** For every \(u,v\ge0\) integral and every nonzero
coefficient direction, \(X^uY^v(A+B X+C Y^2)\) has a nonzero
factorial moment among its first three orders. This statement
holds over every characteristic-zero field.

The first stripped form is \(A+(u+1)B+(v+1)(v+2)C\).
Use it to eliminate \(A\). For the remaining orders define
\[
 T_m(A,B,C)=\sum_{i+j\le m}
 \frac{m!}{i!\,j!\,(m-i-j)!}
 A^{m-i-j}B^i C^j
            (mu+1)_i(mv+1)_{2j}.
\]
In the chart \(B=1\), substitute \(C=z\) and
\(A=-(u+1)-(v+1)(v+2)z\) into \(T_2,T_3\).
The resulting quadratic and cubic have a resultant of total
degree eighteen in \(u,v\), with 91 strictly positive nonzero
coefficients, including the constant \(382464\). In the
remaining direction \(B=0\), the quadratic evaluates at
\(C=1\) to the leading coefficient of the chart quadratic.
It has positive nonzero coefficients and constant \(20\).
Thus neither chart can contain a common zero at nonnegative
integer \(u,v\). The
[exact checker](verify_sfc_independent_height_initial_example.py)
recomputes these polynomial identities over the integers.
