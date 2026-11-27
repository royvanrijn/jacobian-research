# Strong Factorial for products of perturbed radial powers

Partition the variables into nonempty blocks \(I_g\), and write
\(L_g=\sum_{i\in I_g}X_i\), \(t_g=|I_g|-1\), and
\(d_g(\nu)=\sum_{i\in I_g}\nu_i\). Consider the presentation
\[
 F=c\prod_j\left(\prod_g L_g^{a_{jg}}
                 +\sum_\nu u_{j,\nu}X^\nu\right)^{e_j},
 \qquad
 e_j\ge1,\quad A_g=\sum_j e_j a_{jg},\quad c\ne0 .
 \tag{1}
\]
The \(a_{jg}\) and exponents are nonnegative integers. The
perturbation supports in different factors may overlap and the
terms may have mixed degrees. Work over a characteristic-zero
field carrying a real-valued valuation extending \(v_p(p)=1\)
on the coefficient ratios.

For each group with \(A_g>0\), write \(A_g=p^{b_g}A_g^\flat\)
with \(p\nmid A_g^\flat\). Call a group descending for a factor
if some perturbation column of that factor has block degree
below its anchor degree. The two conditions are:

1. Within each factor no convex combination of its perturbation
   block-degree columns is strictly below its anchor in two groups.
2. For each factor \(j\) and each group \(g\) descending for it,
   \(A_g^\flat\) divides \(e_j\).

These conditions exactly characterize when the independently
expanded terms, after summing powers of the radial sums, can
have a uniform lower valuation bound at the factorial reference
defined below for some choice of nonzero perturbation coefficients.
When both hold, explicit positive valuation bounds give a
nonzero factorial moment in every window of the actual support
size. There is no upper bound on dimension, degree, multiplicity,
support size or starting order. When either condition fails,
every fixed nonzero coefficient choice has unbounded relative
term valuations for this presentation. This is a statement
about individual grouped terms and gives no vanishing of their sum.

With singleton blocks, (1) is an arbitrary product of powers
of monomial-anchored polynomials. The result extends the
[factorwise certificate](SFC_LOGARITHMIC_TOWER_RESIDUES.md#6-factorwise-expansion-can-recover-an-all-order-window)
by allowing arbitrary total anchors whose prime-to-\(p\) parts
are absorbed by the factor multiplicities. Radial blocks also
extend the
[one-block change of reference](SFC_RADIAL_BLOCK_PERTURBATIONS.md#5-one-descending-radial-block-after-changing-the-reference).
The calculation below uses factorial valuations from
[Granville, Section 2](https://www.cecm.sfu.ca/organics/papers/granville/paper/binomial/html/node3-an.shtml).
For comparison, [Adolphson and Sperber, Theorem 7.14](https://arxiv.org/pdf/2001.03296)
restate Landau's integrality criterion using floors of linear forms.
Here the relevant issue is a uniform lower bound at one prime
after removing a chosen factorial normalization. We work directly
with the rational roots in each fixed-count ratio.

## 1. Reference and grouped expansion

If a group descends for at least one factor, let \(s_g\) be the
least nonnegative residue of \(t_g\) modulo \(p^{b_g}\).
For all other groups let \(s_g=t_g\). Use
\[
                        H_n=c^n\prod_g (A_g n+s_g)!/s_g! .
 \tag{2}
\]
For each factor \(j\), let \(k_{j,\nu}\) be the perturbation
counts, with \(K_j=\sum_\nu k_{j,\nu}\le e_j n\).
Put \(\gamma=\sum_{j,\nu}k_{j,\nu}\nu\), let
\(\Gamma_g=d_g(\gamma)\), and set
\[
 D_{j,g}=\sum_\nu k_{j,\nu}(a_{jg}-d_g(\nu)),
 \qquad D_g=\sum_j D_{j,g}.
\]
The factorial identity for a radial block is
\[
 \Lambda_g(X^{\gamma|_{I_g}}L_g^m)
     =\prod_{i\in I_g}\gamma_i!\,\frac{(m+\Gamma_g+t_g)!}{(\Gamma_g+t_g)!}.
\]
It follows by expanding \(L_g^m\) and taking the degree-\(m\)
coefficient in \(\prod_{i\in I_g}(1-z)^{-\gamma_i-1}\).
The normalized contribution of the chosen counts is
\[
 \prod_j\frac{(e_j n)_{\underline{K_j}}}{\prod_\nu k_{j,\nu}!}
       \prod_{j,\nu}u_{j,\nu}^{k_{j,\nu}}
 \prod_g
       \frac{s_g!\,\prod_{i\in I_g}\gamma_i!}{(\Gamma_g+t_g)!}
       \frac{(A_g n+t_g-D_g)!}{(A_g n+s_g)!}.
 \tag{3}
\]
Every ratio in \(n\) is a finite product or reciprocal product
of linear factors.

## 2. Cancel every integral pole of a fixed count term

Suppose both support conditions hold. If
\(U_g=D_g+s_g-t_g>0\), the row in (3) loses the factors
\(A_g n+s_g-\ell\) for \(0\le\ell<U_g\). Some \(D_{j,g}>0\).
Among these contributing factors choose \(j(g)\) with largest
\(K_j/e_j\). A given factor contributes positively to \(D_{j,g}\)
in at most one group for the chosen counts, by condition 1.
Thus different losing groups choose different factors, even
with overlapping possible descents. We have
\[
 U_g\le D_g\le \sum_{j:D_{j,g}>0} a_{jg}K_j
                          \le A_g K_{j(g)}/e_{j(g)} .
 \tag{4}
\]

A lost factor has root \(\rho=(\ell-s_g)/A_g\).
If \(\rho\in\mathbb Z_p\), then \(p^{b_g}\mid\ell-s_g\).
Since \(0\le s_g<p^{b_g}\), this gives \(\rho\ge0\).
Condition 2 implies \(e_{j(g)}\rho\) is an integer. Also
\(\ell-s_g<D_g\), so (4) gives
\(e_{j(g)}\rho<K_{j(g)}\).
The root is therefore present in the selected falling factorial
\((e_{j(g)}n)_{\underline{K_{j(g)}}}\).
Different lost factors in the group choose distinct factors of
this falling factorial. Their linear coefficient ratios have
valuation \(b_g-v_p(e_{j(g)})\le b_g\).
Every other lost factor has constant valuation
\(v_p(\ell-s_g)<b_g\) on the nonnegative integers.
Hence, after accounting for the selected falling factorial, the
valuation cost in the row is at most
\(b_g U_g\le (b_g A_g/e_{j(g)})K_{j(g)}\).

For the converse, perturb only one factor and fix its counts.
If condition 1 fails, choose counts in that factor giving two
positive \(D_{j,g}\), using rational linear feasibility and
clearing denominators. Scale the counts so that both deficits
exceed the corresponding \(t_g\). The two lost rows then
have a pole factor at \(n=0\), whereas the one active falling
factorial has a simple zero there. Other increasing rows have
positive offsets and give no zero there. Taking \(n=p^\ell\)
shows unbounded negative valuation. All counts are valid for
sufficiently large orders.

If condition 2 fails at \(j,g\), set \(\rho=1/A_g^\flat\).
Use enough copies of a perturbation monomial descending there
that \(D_g>t_g+p^{b_g}\). The lost row contains the factor
\(A_g n-p^{b_g}\) with root \(\rho\in\mathbb Z_p\). But
\(e_j\rho\notin\mathbb Z\), so the only active falling factorial
does not vanish there. Other increasing row factors have
negative roots and cannot cancel it. Along sufficiently large
positive integers tending \(p\)-adically to \(\rho\) the
fixed-count valuation tends to \(-\infty\).
Both obstructions apply for every fixed nonzero choice of
coefficient ratios, which contribute only a constant valuation
to the chosen counts.

## 3. Uniform valuation bounds and support-size windows

Define, for each factor \(j\),
\[
 B_j=\max_{g\text{ descending for }j}\frac{b_g A_g}{e_j},
 \qquad V=\sum_g v_p(t_g!/s_g!),
\]
where the maximum over an empty set is zero. For a perturbation
exponent define
\[
 W_\nu=
 \sum_{g:\,|I_g|\ge2}\left(
      \sum_{i\in I_g}s_p(\nu_i)
      + \boldsymbol 1_{d_g(\nu)>0}\max(0,s_p(t_g)-1)\right).
\]
Here \(s_p\) is the base-\(p\) digit sum. A sufficient explicit
bound for all the conclusions on nonvanishing is
\[
                      v(u_{j,\nu})\ \ge\
                   B_j+V+\frac{W_\nu+1}{p-1}.
 \tag{5}
\]

In fact the angular denominators in (3) can be written
\((s_g!/t_g!)C_g(\gamma)^{-1}\), where
\(C_g(\gamma)=\binom{\Gamma_g+t_g}{t_g,(\gamma_i)_{i\in I_g}}\).
For a singleton \(C_g=1\). In any other group the exact
valuation of \(C_g\) is
\[
       \frac{s_p(t_g)+\sum_{i\in I_g}s_p(\gamma_i)
                                   -s_p(\Gamma_g+t_g)}{p-1}.
\]
When \(\Gamma_g=0\) it is zero. Otherwise charge the constant
\(\max(0,s_p(t_g)-1)\) to any term contributing to that group,
and use subadditivity of digit sums and their submultiplicativity
under multiplication. Thus
\(\sum_g v(C_g(\gamma))\le\sum_{j,\nu} W_\nu k_{j,\nu}/(p-1)\).
The unmatched falling factorials and positive row products are
integral, and Section 2 bounds the residual loss in each
assigned factor by \(B_jK_j\). Consequently (5) bounds the
valuation of every nonanchor summand in (3) below by
\[
             V\left(\sum_jK_j-1\right)
                    +\frac{\sum_{j,\nu}s_p(k_{j,\nu})}{p-1}>0 .
\]
Only the all-anchor term remains modulo the maximal ideal.

For each descending group put \(q_g=\lfloor t_g/p^{b_g}\rfloor\).
In that group the all-anchor ratio in (3) is
\[
             \prod_{s_g<l\le t_g}\frac{A_g n+l}{l}
             \equiv \binom{A_g^\flat n+q_g}{q_g}\pmod{\mathfrak m}.
 \tag{6}
\]
Indeed the factors with \(p^{b_g}\nmid l\) reduce to one,
and the other factors make the indicated binomial.
Nondescending groups contribute one.

If every \(q_g=0\), all moments are nonzero at every order.
Otherwise put \(Q=\max_g q_g\), \(k=\lfloor\log_p Q\rfloor\),
and \(T=\sum_g\lfloor q_g/p^k\rfloor\), over descending groups.
At orders \(n=p^k m\), the base-\(p\) binomial rule reduces
each binomial in (6) to
\[
        \binom{A_g^\flat m+\lfloor q_g/p^k\rfloor}
                      {\lfloor q_g/p^k\rfloor}\pmod p .
\]
Their product is a polynomial over \(\mathbb F_p\) of degree
\(T\) and has constant coefficient one. Among any
\(\min(p,T+1)\) consecutive values of \(m\) some residue is
nonzero, by the polynomial root bound if \(T<p\) and the
residue \(m=0\) otherwise. Thus a sufficient window length
for the original orders is
\[
                         p^k\min(p,T+1)\le Q+\sum_g q_g .
 \tag{7}
\]

Here is the comparison with the actual support. In every group
with \(A_g>0\), the pure coordinate powers of \(L_g^{A_g}\)
have unit coefficients. Products of those terms give at least
\(\prod_{g:A_g>0}|I_g|\) unit-coefficient monomials in the
anchor of (1). Every coefficient contribution using a
perturbation has strictly positive valuation by (5), so
these monomials survive in \(F/c\).
If at least two \(q_g>0\), the product of the corresponding
block sizes is at least \(1+\sum_g q_g+Q\) by expanding
\(\prod_{q_g>0}(1+q_g)\) and keeping a product of the maximal
\(q_g\) with another positive \(q_g\). Thus (7) suffices.
If just one \(q_g>0\) and \(b_g>0\), its residue period is at
most \(p^{k+1}\le p q_g\le t_g\). If \(A_g=1\), the
\(t_g\)-th forward difference of \(\binom{n+t_g}{t_g}\)
is one, giving one nonzero residue per \(t_g+1\) orders.

In the remaining case \(p\nmid A_g\) and \(A_g>1\) with
one positive \(q_g=t_g\), (7) is at most \(2t_g\).
The radial power \(L_g^{A_g}\) has at least \(2t_g\) unit
coefficients. For each nonzero base-\(p\) digit \(a\) of
\(A_g\), the Frobenius expansion contributes
\(\binom{a+t_g}{t_g}\) independent choices, and the digit
choices give distinct monomials. If there are two or more
nonzero digits the count is at least \((t_g+1)^2\).
Otherwise the sole digit is \(A_g\ge2\), giving at least
\(\binom{t_g+2}{2}\ge 2t_g\). The corresponding monomials
again cannot be canceled by the perturbations. These cases
prove the support-size bound for every phase under (5).

For a mixed planar example at \(p=2\), take
\[
       F=(x^2y+u y^2)^3(x^6y+v x^7).
\]
The total anchor is \((12,4)\). The factor decreasing in
\(x\) has multiplicity three, divisible by the odd part
of twelve. The factor decreasing in \(y\) has multiplicity
one, and four is a power of two. Formula (5) is satisfied,
for example, if \(v_2(u)\ge10\) and \(v_2(v)\ge12\).
Thus all its moments are nonzero. The expanded support has
a convex combination strictly below the anchor in both
coordinates. It cannot be treated by the uniform direct
pure-term comparison. For a genuine radial example at the
same prime take \(((x+y)^2+u x+v)^3(x+y)^6\). Its block
anchor is twelve, and the only descending factor has
multiplicity three. Sufficiently high coefficient valuations
give all-order nonvanishing here too.

## 4. Analytic comparison after canceling the poles

Let the coefficients lie in a finite extension of \(\mathbb Q_p\).
There is an analytic version under the two support conditions
and the weaker bounds
\[
                         v(u_{j,\nu})>B_j+\frac1{p-1}.
 \tag{8}
\]
Each grouped term (3), viewed as a function of an order variable
\(z\), belongs to the Tate algebra on the closed unit disk.
Use the exact cancellations of Section 2. The unmatched
denominator factors have roots outside the disk, so their
inverses are scalar multiples of convergent geometric series
with integral coefficients. After removing the scalar costs
bounded in Section 2, the remaining numerator factors also
have integral coefficients.

Let \(K=\sum_j K_j\). The same valuation calculation now
bounds the minimum coefficient valuation of (3).
The angular loss is \(O(1+\log(K+1))\) instead of linear,
by bounding \(\gamma_i\) by a fixed multiple of \(K\) in the
digit sums. Thus the coefficient valuation grows at least as
\[
                    \eta K-O(1+\log(K+1)),\qquad \eta>0,
\]
where \(\eta\) is the least margin in (8). Summing (3) over
all count vectors defines a restricted analytic series
\(\Psi(z)\). For ordinary \(n\) the exact normalized moment
uses just \(K_j\le e_j n\). Every extra term of the analytic
series has \(K>n\), giving effective constants \(c_0>0,C_0\)
such that
\[
         v\left(\Psi(n)-\frac{\Lambda(F^n)}{H_n}\right)
                                   \ge c_0 n-C_0.
 \tag{9}
\]
Excess count terms need not be exactly zero at such an \(n\),
because pole cancellation may have removed their falling
factorial zeros.

For any fixed width \(w\), absence of a common zero on
\(\mathbb Z_p\) of \(\Psi(z),\ldots,\Psi(z+w-1)\) is a
finite-precision certificate for all sufficiently late
\(w\)-windows by compactness and (9). If \(\Psi\ne0\),
its zero set in \(\mathbb Z_p\) is finite. At a common zero
\(\rho\), a nonzero component has finite vanishing order
and its valuation is affine in \(v_p(n-\rho)\) sufficiently
near \(\rho\). Hence any sufficiently late exact zero
window starts at \(n\) with \(v_p(n-\rho)\ge c_1 n-C_1\)
for one of finitely many common centers. This is the same
analytic spacing argument as in the
[earlier residue interpolation framework](FACTORIAL_LINEAR_PRECISION_ANALYTIC_INTERPOLATION.md#spacing-of-exact-zero-windows-at-a-fixed-phase),
with the factorial reference (2). Algebraic common centers
would exclude all sufficiently late exact zero windows by
the polynomial height bound at integer \(n\).
The hypotheses on common zeros are independent obligations.

With singleton blocks and a unique least slope, (9) identifies
the full integer-center tower expansions of the
[logarithmic approach](SFC_LOGARITHMIC_TOWER_RESIDUES.md).
At any such center \(r\) the formal expansion is the Taylor
series of \((z-r)^m\Psi(z)\), where \(m\) is that approach's
pole-order bound. Indeed their values differ by (9) on orders
\(r+p^\ell\), faster than every fixed power of \(p^\ell\).
Uniqueness of the asymptotic coefficients proves the identity.
In particular no full germ is flat in this grouped class
when \(\Psi\not\equiv0\). Under (5) with a strict margin in
(8), the singleton case has \(\Psi(n)\equiv1\pmod{\mathfrak m}\)
by (9) and the all-anchor comparison for large \(n\), so this
nonidentity is automatic.

The classification concerns a presented product and fixed
block partition. It leaves general irreducible phases,
unavailable factorizations, and stronger cancellation
between terms outside the pole-canceling class open.

## 5. Local analytic comparison at a strictly isolated slope

The uniform bounds (8) can be weakened if the order disks
are allowed to be smaller. Keep the two support conditions
and assume just that in every factor the radial anchor has
strictly least effective slope, that is,
\[
     \delta=\min_{j,\nu}
        \left(v(u_{j,\nu})+
                    \frac{|\nu|-\sum_g a_{jg}}{p-1}\right)>0 .
 \tag{10}
\]
Then the normalized moments still admit a locally analytic
comparison on all of \(\mathbb Z_p\) with an effective linear
precision bound of the form (9). The coefficient ratios
themselves can have negative valuation. This comparison
includes phases of arbitrarily high degree and multiplicities.

To prove it, fix a disk \(n=r+p^\ell z\) with \(0\le r<p^\ell\).
For \(a>0\), a consecutive interval of \(N\) factors
\(a n+b\) has Gauss valuation on this disk
\[
             \frac{N}{p-1}+O\bigl(Np^{-\ell}+\ell+v_p(a)\bigr).
 \tag{11}
\]
For each factor its valuation in this algebra is
\(\min(v(a r+b),\ell+v_p(a))\). Sum by counting multiples
of \(p^j\) for \(1\le j\le\ell+v_p(a)\). Within the interval
the count is \(N/p^j+O(1)\). This gives (11) with constants
independent of the starting offset of the interval.
Gauss valuation is multiplicative in the Tate algebra,
as seen by normalizing two nonzero series and reducing them
to nonzero polynomials over the residue field.

After the pole cancellation in Section 2, a grouped term is
analytic on the disk. Its Gauss valuation can therefore be
computed by subtracting the valuations of the original
denominator products from those of the original numerator
products, even when an individual original denominator
vanishes at a point of the disk. Indeed the analytic series
times that denominator equals the numerator.
Use (11) in (3). Its main linear contributions from
falling factorials and \(\prod k_{j,\nu}!\) cancel.
The radial products contribute
\((\Gamma_g-\sum_j a_{jg}K_j+t_g-s_g)/(p-1)\).
The same estimate for a row with \(A_g=0\) follows directly
from Legendre's formula for its constant factorial ratio.
The angular loss and factorial digit errors are logarithmic
in \(K+2\). Consequently for some \(C\) independent of \(\ell,r,K\)
and a constant \(C_\ell\) independent of \(r,K\), the minimum
coefficient valuation of the grouped term is at least
\[
                   (\delta-Cp^{-\ell})K-C_\ell\log(K+2).
\]
Choose one \(\ell\) with \(Cp^{-\ell}<\delta\).
On each of its finitely many disks the sum over count vectors
converges in the Tate algebra. Terms outside the exact moment
range again have \(K>n\), giving (9) with one choice of
constants for all disks.

The window-separation test of Section 4 remains valid for
this locally analytic comparison. Its spacing consequence
uses a nonidentity check on each disk for at least one of
the shifted components being compared. In the singleton
case, if the resulting analytic function is nonzero on
the disk of an integer center, its full normalized germ
there is nonzero by the tower comparison. Neither the
pole criterion nor the slope inequality by itself proves
that nonidentity in this weaker valuation region.

## 6. Directional factors extracted from the Laurent logarithm

There is also a bridge from logarithmic cancellation to an
all-order grouping when no polynomial factorization is given.
Use singleton blocks here. Let \(F=cX^a U\), with \(U=1+Z\) a
Laurent polynomial and \(F\) ordinary. Work over a complete
nonarchimedean extension of \(\mathbb Q_p\). Write \(\mu\) for
minimum coefficient valuation, and suppose
\(\tau=\mu(Z)>1/(p-1)\).
The coefficient algebra here is the restricted Laurent algebra:
the coefficient valuations tend to infinity as the exponent vector
leaves finite subsets of \(\mathbb Z^d\).
The exponential and logarithm converge in this algebra on the
stated small disk and satisfy their usual formal identities.
The bounds \(v_p(j!),v_p(j)\le(j-1)/(p-1)\) for \(j\ge1\)
ensure the convergence and preservation of valuation.

Suppose \(L=\log U\) has no Laurent monomial with negative exponents
in two coordinates. Split it as \(L=L_0+\sum_i L_i\), where \(L_0\)
has all exponents nonnegative and \(L_i\) has its only negative
coordinate in position \(i\). Then
\[
                 U=U_0\prod_i U_i,\qquad U_i=\exp L_i.
 \tag{11}
\]
This is the unique factorization by restricted Laurent units close
to one with these support directions: \(U_0\) has nonnegative
exponents, and each nonconstant perturbation term of \(U_i\)
is strictly negative in coordinate \(i\) and nonnegative in the
others. Indeed exponentiation preserves the indicated cones,
and taking logarithms recovers their disjoint coefficient
projections. It also gives \(\mu(U_i-1)\ge\tau\).

In the negative coordinate each \(U_i\) already has a finite
lower exponent bound \(-a_i\). To see it, write
\(U_i=U/(U_0\prod_{j\ne i}U_j)\). Every other factor and its
inverse have nonnegative \(i\)-exponents. The Laurent polynomial
\(U\) has no \(i\)-exponent below \(-a_i\). The product is
computed by convergent coefficient sums and cannot create a
lower \(i\)-exponent. In particular \(a_i=0\) forces \(U_i=1\).
Thus \(U_i\) is a polynomial in \(x_i^{-1}\) of degree at most
\(a_i\), with coefficients restricted analytic in the other
variables. These individual factors can have infinitely many
terms in the other coordinates.

For \(a_i>0\), write \(a_i=d_i p^{b_i}\) with \(p\nmid d_i\),
and consider \(V_i=\exp(L_i/d_i)\). Assume \(V_i\) has no
\(i\)-exponent below \(-p^{b_i}\) whenever \(L_i\ne0\).
This extra condition is automatic if \(a_i\) is a power of
\(p\), by the finite lower bound just proved. In general
\(V_i\) is the unique \(d_i\)-th root of \(U_i\) close to one
in the same directional algebra. The condition asks that
this root be a polynomial in \(x_i^{-1}\) of degree at most
\(p^{b_i}\). Its candidate coefficients up to that degree are
determined recursively by those of \(U_i\), since \(d_i\)
is a unit, and the candidate's \(d_i\)-th power can be
compared with \(U_i\).

Under this condition,
\[
 F=c\,U_0\prod_{a_i>0}(x_i^{p^{b_i}}V_i)^{d_i}.
 \tag{12}
\]
All factors are now restricted ordinary series. If
\(\tau\ge\max_{i:L_i\ne0}p^{b_i}/(p-1)\), with a vacuous
maximum for no descending factors, every normalized moment
\(\Lambda(F^n)/(c^n\prod_i(a_i n)!)\) is a unit congruent to
one. Although (12) uses analytic factors, the argument gives
ordinary polynomial moments of the original \(F\).

For that argument, \(\Lambda\) extends continuously to
restricted ordinary series because the monomial values
\(\alpha!\) have nonnegative valuation. Expand each factor in
(12) at its own order \(d_i n\) and \(U_0\) at order \(n\).
Every nonanchor term of the \(i\)-th factor can descend only
in coordinate \(i\), by at most \(p^{b_i}\).
If it is used \(K_i\) times, row \(i\) has total exponent at
least \(p^{b_i}(d_i n-K_i)\). By prime-power factorial scaling
the factorial loss is at most
\[
     \frac{p^{b_i}-1}{p-1}K_i
            +v_p\big((d_i n)_{\underline{K_i}}\big).
\]
The falling-factorial part is supplied by the same factor.
Nonanchor terms of \(U_0\) do not lower any exponent, even
if its constant coefficient receives a small correction.
The remaining count factorials cost at most the number of
perturbations divided by \(p-1\), with strictly positive
digit-sum savings for every nonzero count. Thus all
perturbation contributions vanish modulo the maximal ideal.
The expansion converges at each fixed \(n\), since there
are only finitely many factors in the powers and their
coefficient tails tend to zero. Dividing by the reference
at that \(n\) is a fixed scalar operation.

With the strict inequality \(\tau>\max_{i:L_i\ne0}
p^{b_i}/(p-1)\), this also proves that the full integer-center
tower germs are nonzero. The anchor is strictly least in
effective slope: in the factor \(x_i^{p^{b_i}}V_i\) the total
degree can drop by at most \(p^{b_i}\) and each perturbation
coefficient has valuation at least \(\tau\). The other factor
cannot lower total degree. The strict margin survives through
the product and through any cancellation among its terms.
Since every normalized moment has valuation zero, the
valuation of the regularized moment on a sequence \(r+p^\ell\)
is exactly \(m\ell\). Formula (2) of the tower expansion
therefore forces \(B_j(r)=0\) for \(j<m\) and \(B_m(r)\) to
be a unit, at every fixed integer center.

This gives a sufficient structural bridge from absence of
mixed negative logarithmic monomials to all-order
nonvanishing. The leading weighted residues in
[the tower formula](SFC_LOGARITHMIC_TOWER_RESIDUES.md#3-the-same-residue-functional-for-adjacent-offsets)
can vanish by cancellations of those monomials after
weighting, and their vanishing alone does not assert the
support condition used in (11). Establishing a separation
result for those weighted cancellations remains open.

For an example where the analytic factor need not be a polynomial
in the other variables, take \(F=x+t x^2+t y\) with \(v(t)>1/(p-1)\).
The logarithm of \(F/x\) has no negative \(y\)-exponents.
The factor \(x U_1\) of (12) is \(x-\alpha(y)\), with
\(\alpha(y)=(\sqrt{1-4t^2y}-1)/(2t)\) the small root of
\(t x^2+x+t y\). This is analytic and nonpolynomial in \(y\)
when \(t\ne0\). The remaining factor is \(1+t x+t\alpha(y)\).
The new argument permits these factors without treating
them as ordinary polynomial factors.

## 7. Analytic factors of a translated product of coordinate powers

Here is a support condition that constructs the directional factors
without a logarithm test or a polynomial factorization. Let
\[
                  G(X)=\prod_i X_i^{a_i}+P(X),
\]
with all \(a_i\) powers of \(p\), including one. For each ordered
pair \(i\ne j\) choose an integer \(w_{ij}\ge0\).
Suppose every perturbation exponent \(\nu\) satisfies
\[
                  w_{ij}(\nu_i-a_i)+\nu_j-a_j\ge0
                  \qquad(i\ne j).
 \tag{13}
\]
Let the perturbation coefficients and the translation coordinates
\(h_i\) all have valuation at least
\(\tau\ge\max_i a_i/(p-1)\). Then for \(F(X)=G(X+h)\) every
ordinary factorial moment is nonzero, and
\(\Lambda(F^n)/\prod_i(a_i n)!\equiv1\pmod{\mathfrak m}\).
This holds for arbitrary dimensions and perturbation degrees
satisfying (13), over a complete valued extension and hence
for phases over its characteristic-zero subfields. Translation
can give \(F\) columns below the anchor in several coordinates.

Let \(M_i=\prod_{j\ne i} X_j^{w_{ij}}\). Substitute \(X_i=M_i T\).
By (13) the resulting polynomial is
\[
 M_i^{a_i}\prod_{j\ne i}X_j^{a_j}\bigl(T^{a_i}+E_i(T,X_{\ne i})\bigr),
 \qquad \mu(E_i)\ge\tau.
\]
For clarity, the required special Weierstrass factorization
and division follow from the same contraction.
For a series in \(T\), split off its powers of degree below
\(a_i\), and divide its remaining part by \(T^{a_i}\).
Solve \(E_i=T^{a_i}U_i+V_i+U_iV_i\) by applying this split to
\(E_i-U_iV_i\), where \(V_i\) has degree below \(a_i\).
Both unknowns stay of valuation at least \(\tau>0\), and the
product difference increases the valuation distance by at
least \(\tau\). Hence
\(T^{a_i}+E_i=(1+U_i)(T^{a_i}+V_i)\).
For division by \(T^{a_i}+V_i\), split the input minus \(V_i\)
times the provisional quotient in the same way. This gives
existence and uniqueness of the quotient and low-degree
remainder for every restricted series.

Undo the substitution in the monic factor by putting
\[
 W_i(X)=X_i^{a_i}+\sum_{k<a_i}
              [T^k]V_i(T,X_{\ne i})\,M_i^{a_i-k}X_i^k .
\]
It is an ordinary restricted series and differs from
\(X_i^{a_i}\) by valuation at least \(\tau\). Divide \(G\)
by \(W_i\) in \(X_i\). After setting \(X_i=M_iT\), both
\(G\) and the quotient term are divisible by \(T^{a_i}+V_i\).
The substituted remainder has degree below \(a_i\), so
uniqueness makes it zero. Since \(M_i\) is nonzero in a
domain, the original remainder was zero. Thus every \(W_i\)
divides \(G\) in the Tate algebra.

Their product divides \(G\) as well. In the quotient by
\(W_j\), use the free low-degree remainder module in \(X_j\)
from the division just proved. Modulo the maximal ideal,
multiplication by \(W_i\) for \(i\ne j\) is multiplication
by \(X_i^{a_i}\) in the polynomial ring with \(X_j^{a_j}=0\).
It is injective there. To lift injectivity to the quotient,
scale any purported nonzero kernel vector so its remainder
has valuation zero and then reduce. This contradicts the
residue-ring statement. Successive division now proves the
product divisibility. Write \(G=W_0\prod_i W_i\). Comparing
with the anchor and using multiplicativity of Gauss
valuation gives \(\mu(W_0-1)\ge\tau\).

After \(X\mapsto X+h\), the same bounds hold, and the
\(i\)-th factor is still \(X_i^{a_i}\) plus small terms that
can descend only in coordinate \(i\). Its leading coefficient
is one. The unit \(W_0(X+h)\) has no descending coordinate
relative to anchor one. Applying the convergent factorwise
expansion of Section 6 proves the asserted congruence.

For example, in two variables with \(a_1=a_2=1\) one can
take \(w_{12}=w_{21}=1\), so (13) allows every perturbation
monomial of total degree at least two. With \(t\ne0\) small,
\(G=xy+t x^3+t y^3\) is absolutely irreducible. If its
projective cubic \(XY Z+tX^3+tY^3\) had a line factor
containing \(Z\), substitution of \(Z\) as a linear form
in \(X,Y\) would leave the two nonzero pure cubics. A line
independent of \(Z\) also fails: away from the axes it leaves
a nonzero term involving \(Z\), and on either axis it leaves
a nonzero pure cubic. This excludes every possible line
factor. The translated polynomial remains irreducible.
For small translations with nonzero resulting constant,
and with the strict slope margin, its support has a
two-coordinate descent from the mixed anchor. It is thus
outside the fixed pure termwise normalization from the
[relative valuation obstruction](FACTORIAL_RELATIVE_VALUATION_AUTOMATA.md#4-an-obstruction-and-an-exact-stable-perturbation-test).
The analytic factors still give all-order nonvanishing.
