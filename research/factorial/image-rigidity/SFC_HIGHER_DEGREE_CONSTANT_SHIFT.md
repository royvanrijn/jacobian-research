# Constant shifts of phases near higher-degree monomials

Let \(K\) be a complete nonarchimedean field of characteristic
zero with a real valuation extending \(v_p(p)=1\). The ordinary
factorial functional is \(\Lambda(X^\nu)=\prod_i\nu_i!\).
Fix a nonconstant monomial \(M=X^a\) of total degree \(D\ge2\),
and write \(\lambda=1/(p-1)\).

For \(v(t)\ge D\lambda\), every window of \(D\) consecutive
orders for \(t+M\) has a nonzero factorial moment. This holds
also at the valuation boundary. There is a quantitative version
at every start, using the adjoint of the order recurrence.
Under strict inequality the estimate extends to phases \(t+Q\)
whose \(M\)-normalized moments are sufficiently close to one
in a specified analytic norm. In particular, one can take
the independently perturbed powers
\[
                Q=(1+P_0)\prod_{a_i>0}(X_i+P_i)^{a_i}              \tag{1}
\]
with arbitrary polynomial perturbation supports of sufficiently
deep coefficient valuation. No dimension, anchor degree,
perturbation degree or starting order is bounded. If \(t+Q\)
has at least \(D\) monomials, the resulting window suffices for
its actual support-size assertion. With two prime-power anchor
exponents, Section 5 reduces the window length to two, regardless
of the total degree, and includes the critical scalar valuation.

The \(D=2\) reference case is the
[adjacent-pair determinant argument](SFC_MIXED_ANCHOR_CONSTANT_SHIFT.md#1-a-quantitative-reference-pair-including-the-boundary).
The higher-order proof uses a single adjoint solution.
The
[fixed-center residue recurrence](SFC_LOGARITHMIC_TOWER_RESIDUES.md#8-a-test-of-the-window-transfer-on-a-product-plus-a-constant)
already gives \(D\)-windows along prime-power approaches for the
squarefree product monomial. The argument below covers arbitrary
ordinary starts and exponent profiles. The
[uniform eventual two-window theorem](MULTIVARIATE_CONSTANT_PLUS_MONOMIAL_TAIL.md)
and the [degree-three two-window theorem](CONSTANT_PLUS_TOTAL_DEGREE_THREE_ALL_ORDERS.md)
have unrestricted scalar coefficients and different scopes.
For higher degrees, the two-window problem for arbitrary
exponent profiles remains open.

## 1. The order recurrence and its adjoint

Use the reference normalized moments
\[
 b_m=\frac{\prod_i(a_i m)!}{m!},\qquad e_m=\frac{t^m}{m!},
 \qquad V_n=\sum_{k=0}^n e_k b_{n-k}.
\]
Set all three sequences to zero at negative indices.
Every \(b_m\) for \(m\ge0\) is a positive integer because at
least one \(a_i\) is positive. For \(m\ge1\), \(b_m=P(m)b_{m-1}\)
with
\[
 P(Z)=\frac1Z\prod_{a_i>0}\prod_{\ell=0}^{a_i-1}(a_i Z-\ell)
           \in\mathbb Z[Z],\qquad
 A=\prod_{a_i>0}a_i^{a_i}.
\]
The polynomial \(P\) has degree \(D-1\) and leading
coefficient \(A\). Expand
\[
           P(Z-k)=\sum_{j=0}^{D-1}q_j(Z)(k)_{\underline j}.
\]
Each \(q_j\in\mathbb Z[Z]\), by the integral change from
monomials in \(k\) to falling factorials. In particular
\(q_{D-1}=(-1)^{D-1}A\). Inserting the recurrence of \(b\)
into the convolution gives
\[
       (\mathcal AV)_n=e_n,\qquad
       (\mathcal AX)_n=X_n-\sum_{j=0}^{D-1}
                          q_j(n)t^j X_{n-1-j}.                     \tag{3}
\]
The missing \(m=0\) recurrence supplies \(e_n\). For
\(t\ne0\), the formal adjoint on the other side of a cut is
\[
 (\mathcal A^\ast Y)_m
       =Y_m-\sum_{j=0}^{D-1}q_j(m+j+1)t^jY_{m+j+1}.
\]
Consider the finite sums
\[
                       \psi_m^{[T]}
                =\sum_{h=0}^T\frac{(-t)^h}{h!\,b_{m+h}}
                                          \quad(m\ge0).              \tag{4}
\]
In \((\mathcal A^\ast\psi^{[T]})_m\), every coefficient of
\(t^h\) for \(0\le h\le T\) cancels. Indeed
\(q_j(Z)=(-1)^j\nabla^jP(Z)/j!\), so
\[
   \sum_{j=0}^{D-1}(-1)^j(h)_{\underline j}q_j(m+j+1)
       =\sum_{j=0}^{D-1}\binom hj\Delta^jP(m+1)=P(m+h+1).
\]
Here \(\nabla,\Delta\) denote the backward and forward
unit differences. Dividing by \(b_{m+h+1}\) gives
\(1/b_{m+h}\), as required in (4).
The remaining summands are signed multiples of
\[
      \frac{q_j(m+j+1)t^h}{(h-j)!\,b_{m+h+1}},
      \qquad T<h\le T+D-1,\qquad h-T\le j\le D-1.                    \tag{5}
\]

## 2. A unit boundary form at every ordinary order

Let \(\alpha=v(t)\ge D\lambda\). When \(t=0\), all normalized
moments \(b_n\) are nonzero, so the window result is immediate.
For \(t\ne0\), write \(s_p(m)\) for the base-\(p\) digit sum.
Legendre's formula and subadditivity of digit sums give
\[
  v\bigl(h! b_{m+h}\bigr)\le
       (D-1)\lambda m+D\lambda h+\lambda s_p(m).
\]
In fact \(v(b_s)=(D-1)\lambda s
  -\lambda\sum_i s_p(a_i s)+\lambda s_p(s)\), and
\(s_p(m+h)\le s_p(m)+s_p(h)\). Hence for all \(T\),
\[
           v(\psi_m^{[T]})\ge -(D-1)\lambda m-\lambda s_p(m).        \tag{6}
\]
The sum \(\sum_{m=0}^N e_m\psi_m^{[T]}\) is a unit for any
nonnegative \(N,T\). Its term \(m=h=0\) is one. Every other
term has positive valuation, since up to its sign it is
\[
     \frac{t^{m+h}}{\prod_i(a_i(m+h))!}\binom{m+h}{m}
\]
and the factorial valuations leave a positive digit sum for
\(m+h>0\) even when \(\alpha=D\lambda\).

For each fixed \(N\), choose \(T=p^L-C\) with fixed integer
\(C>N+D\) and large \(L\). Then every summand of (5) for
\(0\le m\le N\) has valuation tending to infinity as
\(L\to\infty\), even at \(\alpha=D\lambda\). Applying
Legendre to its denominator and using
\(s_p(m+h+1)\le s_p(h-j)+s_p(m+j+1)\) gives a lower bound
\[
 -(D-1)\lambda(m+1)+\lambda j-\lambda s_p(m+j+1)
                          +\lambda\sum_i s_p(a_i(m+h+1)).
\]
The last sum grows without bound:
\(m+h+1=p^L-c\) for an integer \(c\) in a fixed finite
positive range. For each \(a_i>0\),
\(a_i(p^L-c)=(a_i-1)p^L+(p^L-a_i c)\), and the latter
summand contains an increasing string of digits \(p-1\).

Subtracting the two finite pairings across the cut \(N\) in
(3) gives the identity
\[
 \begin{split}
   \sum_{m=0}^N e_m\psi_m^{[T]}
      -\sum_{m=0}^N V_m(\mathcal A^\ast\psi^{[T]})_m
      &=
        \sum_{\substack{N-D+1\le m\le N,\ m\ge0\\
                        1\le d\le D,\ m+d>N}}
           V_m q_{d-1}(m+d)t^{d-1}\psi_{m+d}^{[T]} .                 \tag{7}
 \end{split}
\]
Terms with row and column on the same side of the cut
cancel exactly, and \(V_m=0\) for \(m<0\). At the chosen \(T\)
the left side is a unit, since the adjoint errors times
the finitely many \(V_m\) have positive valuation.
This proves that the window on the right cannot vanish.
For a quantitative form let \(l_p(j)=1+\lfloor\log_p j\rfloor\)
for \(j\ge1\). At a start \(n\ge0\), use \(N=n+D-1\).
By (6) and \(\alpha\ge D\lambda\), every coefficient of
\(V_m\) on the right has valuation at least
\(-(D-1)\lambda(n+D)-l_p(n+2D)\), since \(m+d\le N+d\)
and \(\lambda s_p(m+d)\le l_p(n+2D)\). Thus
\[
       \min_{0\le j<D}v(V_{n+j})\ \le\
       G(n):=(D-1)\lambda(n+D)+l_p(n+2D).                          \tag{8}
\]
No convergence of the infinite adjoint sequence was used at
the boundary. All reference arguments apply to the factorial
moments over a valued characteristic-zero field by completion.

## 3. Stable forward transfer for general anchor exponents

Suppose now \(D\ge2\), \(t\ne0\) and \(\alpha>D\lambda\).
Write \(\delta=\alpha-D\lambda\), \(u=v_p(A)\), and set
\[
      L=(D-2)\big((D-1)\alpha+u\big),\qquad
      H(n)=\max(0,n-D+2).
\]
Let \(\Phi\) be an analytic function given by a restricted
series \(\Phi(p^{-u}Z)\), with the infimum of coefficient
valuations of \(\Phi(p^{-u}Z)-1\) at least \(\mu\).
Assume
\[
             \mu> L+\sup_{n\ge0}\big(G(n)-(\alpha-\lambda)H(n)\big).
                                                                         \tag{9}
\]
The supremum is finite and effectively bounded, since its
linear drift is \(-\delta n\) and the remaining growth is
logarithmic. If
\[
                     \frac{\Lambda(Q^m)}{\prod_i(a_i m)!}
                                      =\Phi(m)                       \tag{10}
\]
at every ordinary \(m\ge0\), then the normalized moments
\(\Lambda((t+Q)^s)/s!\) have the same minimum valuation as
\(V_s\) on every \(D\)-window. A precise error extension of
(10) is given below.

For the proof, use the polynomial basis indexed by
\(h\ge0\), \(0\le\ell<D-1\):
\[
                         Z^\ell\prod_{i=1}^h P(Z+i).
\]
The degrees uniquely run over all nonnegative integers.
After scaling \(Z=p^{-u}W\), normalize each basis term
by its leading coefficient \(A^h p^{-u(h(D-1)+\ell)}\).
It becomes monic integral. Indeed the leading coefficient
of each \(P(p^{-u}W+i)\) has valuation
\(u-u(D-1)\) and each lower coefficient has valuation
at least that much, since all coefficients in \(P(Z+i)\)
are integral. The leading coefficient of the full scaled
basis element has nonpositive valuation. Triangular
changes between the monic integral basis and monomials
have integral coefficients. Applying this to polynomial
truncations shows that the expansion of \(\Phi-1\) in the
original basis has coefficients of valuation at least
\(\mu\), tending to infinity.

For \(s,h\ge0\) put
\[
                         T_{s,h}=\sum_{k=0}^s e_k b_{s+h-k}.
\]
At \(0\le h<D\), it equals \(V_{s+h}-\sum_{i=1}^h e_{s+i}b_{h-i}\).
At \(h\ge D\), the same computation as (3) gives
\[
    T_{s,h}=\sum_{j=0}^{D-1}q_j(s+h)t^j T_{s,h-1-j}
                                                +e_s C_{s,h},       \tag{11}
\]
with \(v(C_{s,h})\ge0\). Indeed for index \(j\) the shifted sum after
applying \((k)_{\underline j}\) ends at \(s-j\); completing it
to \(s\) subtracts terms \(t^j e_{s-b}b_{h-1-j+b}\)
with \(0\le b<j\) and \(b\le s\). Since
\(t^j e_{s-b}=e_s t^{j-b}(s)_{\underline b}\), all the
coefficients in the resulting boundary correction are integral.
If \(j>s\) this subtracts the entire completed sum, with
the same formula. Thus each \(T_{s,h}\) is a linear
combination with integral coefficients of
\(V_s,\ldots,V_{s+D-1}\) and \(e_s,\ldots,e_{s+D-1}\).

After convolution with the basis term of index \(h,\ell\),
the factor \(\prod_{i=1}^h P(m+i)\) replaces \(b_m\) by
\(b_{m+h}\). Expand \((s-k)^\ell\) in falling factorials
in \(k\), with integral coefficients. The contribution is
an integral combination of \(t^i T_{s-i,h}\) for
\(0\le i\le\min(\ell,s)\). Therefore, using the series
coefficients above, the perturbed normalized value at
index \(s\) is \(V_s\) plus coefficients of valuation
at least \(\mu\) times the reference and forcing values
at indices between \(\max(0,s-D+2)\) and \(s+D-1\).
Every sum converges by the bounded coefficients in (11).

For \(n\le s<n+D\), reduce the reference values at indices
outside this window using (3). Forward recurrence is
integral. In a backward step, solve the equation at
index \(m+D\) for \(V_m\), losing at most
\((D-1)\alpha+u\) from the coefficient of \(V_m\).
At most \(D-2\) backward steps are needed, with every
forcing index at least \(H(n)\). This gives the perturbed
vector as a matrix congruent to the identity applied to
\((V_n,\ldots,V_{n+D-1})\), plus an error of valuation
at least \(\mu-L+(\alpha-\lambda)H(n)\) in every entry.
Indeed the matrix correction has coefficient valuations at
least \(\mu-L>0\), and \(v(e_j)\ge(\alpha-\lambda)j\).
The error is above (8) by (9). A matrix congruent to the
identity and its inverse are integral, so the minimum
valuation is preserved. This proves the transfer.

One can include an error \(E_m\) on the right of (10)
with \(v(E_m)\ge \eta m+\xi\), \(\eta>0\).
For example it suffices that
\[
 \xi>\sup_{n\ge0}\left(
        G(n)-(D-1)\lambda n-\min(\eta,\delta)n
                        +\sum_{a_i>0}l_p(a_i(n+D-1)+1)\right).       \tag{12}
\]
This supremum is again finite and effective.
The extra binomial sum has summands \(t^k b_m E_m/k!\)
with \(k+m=s\). Legendre's formula bounds their valuations
below by \((D-1)\lambda s+\min(\eta,\delta)s+\xi-
\sum_{a_i>0} l_p(a_i s+1)\).
On the window this is above (8) by (12).

## 4. High-degree factors and actual support size

For (1), each \(P_i\) can contain arbitrary monomials in
all variables. Use the polynomial exponents as written,
including multiplicity \(a_i\) on the factor near \(X_i\).
For integers \(m\ge0\), the anchor exponent of row \(i\)
before all perturbation contributions is \(a_i m-K_i\),
where \(K_i\) is the count selected from that factor.
Let \(\gamma_i\) be the total row-\(i\) exponent from all
perturbations. Combining the row ratio with the
falling factorial of this factor gives the polynomial
\[
           (a_i z)_{\underline{K_i}}\,
             \frac{(a_i z-K_i+\gamma_i)!}{(a_i z)!}
                    =\prod_{j=1}^{\gamma_i}(a_i z-K_i+j).
\]
The unit factor supplies a separate falling factorial of
length \(K_0\); exponents of rows with \(a_i=0\) contribute
factorials that are integral scalars. On \(z=p^{-u}Z\),
the coefficient loss of each linear factor is at most \(u\).
If \(K=\sum K_i\) and \(B\) is one plus the largest total
degree of a monomial in the perturbations, then the number
of surviving linear factors is at most \(B K\).
The count factorials cost at most \(\lambda K\).

Consequently if all perturbation coefficients have valuation at
least \(uB+\lambda+\tau\) with \(\tau>0\), the Gauss
valuation of a count term on that disk is at least \(\tau K\).
Summing all count polynomials gives a \(\Phi\) with
\(\mu\ge\tau\). At a nonnegative integer \(m\), the omitted
counts have \(K_0>m\) or \(K_i>a_i m\) for some \(i\),
so each such term has valuation at least \(\tau(m+1)\).
This is (10) with an error satisfying (12) for
\(\eta=\xi=\tau\) once \(\tau\) is large enough.
For instance, when \(\tau\ge\delta\), the supremum in
(12) is independent of further increases in \(\tau\).
Thus the size thresholds in (9) and (12) can both be met
effectively in terms of the phase data, the prime and the
valuation of \(t\).

This gives arbitrary-degree ordinary support-size windows
whenever the support size is at least \(D\). In particular
the method includes small high-degree upward perturbations
\(Q=M(1+P_0)\) and small independent polynomial changes of
the coordinate factors. The degree \(D\) in the window
comes from the number of boundary entries in (7).
For sparse phases with fewer than \(D\) monomials, the
same argument gives a longer window. Obtaining a form
involving only support-sized data remains a separate
obligation, already visible for the two-term \(t+M\).

## 5. Two-window reduction for two prime-power anchors

There is an all-order reduction to the actual two-window length
when only two rows are active and their exponents are
\(a=p^b\), \(d=p^c\), with \(b,c\ge0\). More generally let
\[
        \frac{\Lambda(Q^m)}{(a m)!(d m)!}
                              =\Phi_0(m)+E_m ,
\]
where \(\Phi_0\) is restricted analytic on the closed unit
disk, \(\Phi_0-1\) has positive coefficient valuation, and
\(v(E_m)>\lambda s_p(m)\) for every \(m\ge0\).
If \(v(t)\ge(a+d)\lambda\), every two consecutive ordinary
factorial moments of \(t+Q\) contain a nonzero one. It is
enough that the two active rows exist among any number of
variables. An algebraic extension of the valued field in
the proof does not change the ordinary moment zero set.

For the proof write
\[
       U_b(m)=\frac{(p^b m)!}{p^{(p^b-1)\lambda m}m!},\qquad
       B_p(l)=\prod_{\substack{1\le j\le pl\\p\nmid j}}j .
\]
Removing multiples of \(p\) recursively gives
\(U_b(m)=\prod_{\ell=1}^b B_p(p^{\ell-1}m)\).
These are the unsigned factorial unit products underlying
[Morita's \(p\)-adic gamma function, as used by Barsky](https://numdam.org/item/GAU_1977-1978__5__A1_0.pdf).
We need explicit coefficient valuations after removing their
geometric factor and give the corresponding elementary
interpolation here.

If \(p\) is odd, put \(g=(p-1)!\) and
\(J(Z)=\prod_{r=1}^{p-1}(1+pZ/r)\). Then
\[
          B_p(l)=g^l\prod_{j=0}^{l-1}J(j).
\]
The coefficient of \(Z^h\) in \(\log J\), \(h\ge1\),
has valuation at least \(h-v_p(h)\).
The polynomial discrete sum of \(Z^h\) from zero to
\(z-1\) has coefficient valuations at least
\(-\lfloor\log_p(h+1)\rfloor\): express \(Z^h\)
in the integral falling factorial basis and use
\(\sum_{j=0}^{z-1}(j)_{\underline k}=(z)_{\underline{k+1}}/(k+1)\).
It follows that the discrete sum of \(\log J\) is a
restricted series of coefficient valuation at least one.
Indeed \(h-v_p(h)-\lfloor\log_p(h+1)\rfloor\ge1\)
for odd \(p\) and tends to infinity. Exponentiating gives
a restricted series \(A_p\) with \(A_p-1\) of coefficient
valuation at least one and \(\prod_{j=0}^{l-1}J(j)=A_p(l)\).
Thus \(U_b(m)\) is \(g^{(1+p+\cdots+p^{b-1})m}\) times a
restricted order series congruent to one modulo \(p\)
in coefficient norm.

At \(p=2\), pair adjacent odd factors:
\[
          (4Z+1)(4Z+3)=3\,J_2(Z),\qquad
          J_2(Z)=1+\tfrac{16}{3}Z(Z+1).
\]
The logarithmic term of index \(h\ge1\) has degree at most
\(2h\) and coefficient valuations at least \(4h-v_2(h)\).
The same discrete summation loses at most
\(\lfloor\log_2(2h+1)\rfloor\). Hence the exponential of
the discrete sum is a restricted function \(A_2\) with
\(A_2-1\) of coefficient valuation at least three.
For \(l=2q+\epsilon\), \(\epsilon\in\{0,1\}\), this gives
\[
              B_2(l)=3^q A_2(q)(4q+1)^\epsilon .
\]
For the terms \(l=2^{\ell-1}m\) with \(\ell\ge2\), the
factor \(3^q\) is geometric in \(m\) and the remaining
factor is a restricted series close to one.
For \(\ell=1\), adjoin \(\iota\) with \(\iota^2=-1\).
Then \(3^q=\iota^m\iota^{-\epsilon}(-3)^q\).
On each order disk \(m=\epsilon+2Z\), the factor
\((-3)^Z=\exp(Z\log(-3))\) is restricted analytic
close to one, as are \(A_2(Z)\) and \((4Z+1)^\epsilon\).
Since \(v_2(\iota-1)=1/2\), after the single geometric
factor \(\iota^m\) has been removed the remaining order
function has positive disk coefficient valuation relative
to one on both parity disks.

In all cases, there is a scalar \(C\) of valuation
\(\lambda(a+d-2)\) such that
\[
             (a m)!(d m)!=C^m(m!)^2 A_*(m)
\]
where \(A_*\) has positive coefficient valuation relative
to one, on the whole order disk for odd \(p\), and on the
two parity disks for \(p=2\). Products of the corresponding
unit functions retain these estimates. Then \(Q/C\) has
normalized moments \(A_*(m)\Phi_0(m)+A_*(m)E_m\);
the analytic and error assumptions are preserved.
The scalar \(t/C\) has valuation at least \(2\lambda\).

For odd \(p\) the
[unweighted analytic constant shift](SFC_MIXED_ANCHOR_CONSTANT_SHIFT.md#2-adding-a-constant-to-analytically-normalized-moments)
applies. For \(p=2\), the transfer on
[order disks of depth one](SFC_TRUNCATED_FACTORIAL_CONSTANT_SHIFT.md#2-constant-shift-with-only-local-order-analyticity)
also works at scalar valuation \(2\), as can be checked in
its error calculation. In that proof take \(e=1,\lambda=1\)
and \(\alpha\ge2\). The weighted reference bound (8) there
still holds because \(\alpha-\lambda\ge d(n)\). The model
error in (7) there has coefficient losses at most one
in the first entry and two in the second after scaling.
For \(n\ge1\), use \(v(b_n)\ge n+s_2(n)\) and
\(v(b_{n+1})\ge n+1+s_2(n+1)\), and both digit sums are
at least one. This keeps the model errors above \(n\)
for every strictly positive disk margin. The source errors
already satisfy that proof's hypothesis, and order zero
has value one. Thus the critical scalar valuation suffices
on parity disks as well.

The approximation hypothesis for \(Q\) follows, for example,
from the independent factor presentation (1) with two
active rows and small enough coefficients using the count
argument of Section 4 on the closed unit disk \(u=0\).
Another sufficient presentation is
\[
            Q=(1+P_0)(X_1^{p^b}+P_1)(X_2^{p^c}+P_2).
\]
Assume the coefficient valuations in \(P_0\) are at least
\(2\lambda\), in \(P_1\) at least \(b p^b+2\lambda\),
and in \(P_2\) at least \(c p^c+2\lambda\).
If \(K_1\) perturbations are selected from the first
coordinate factor, its row can lose at most \(p^b K_1\)
linear factors from its factorial reference.
Every lost factor with a root in the \(p\)-integral order
disk matches a root of \((z)_{\underline{K_1}}\), because
the anchor exponent \(p^b\) is a prime power.
Matched factors cost \(b\) in coefficient valuation each.
Unmatched factors have valuation less than \(b\) on their
constant term and their inverses have convergent geometric
expansions. This bounds the cost by \(b p^b K_1\), and the
second row is independent. After count factorial costs,
the coefficient valuation of each term is at least
\(\lambda K\), where \(K\) is the total perturbation count.
Summing the series and the omitted terms as in Section 4
proves the approximation, with
\(v(E_m)\ge\lambda(m+1)\).
In particular the two-window conclusion holds for
the pure two-term phase \(t+X_1^{p^b}X_2^{p^c}\) at
every order under the stated scalar condition, even
though the anchor degree can be arbitrarily large.

The same normalization isolates the size of the unresolved scalar
interval. Write \(\omega=\lambda(a+d-2)\) and \(\beta=v(t)-\omega\)
for \(t\ne0\). At an ordinary order \(n\) the moments divided by
\(t^n\) are
\[
             \sum_{m=0}^n (n)_{\underline m} m!
                       (\Phi_0(m)+E_m)A_*(m)(C/t)^m.
\]
Every factor \((\Phi_0(m)+E_m)A_*(m)\) is a unit.
The normalized moment is a unit at each order when \(\beta<0\),
since all nonconstant summands have positive valuation.
At \(p=2\) there is a stronger conclusion when \(\beta<1\):
at each even order \(n\), the falling interval of length \(m\)
has at least \(\lceil m/2\rceil\) even integers, whereas \(m!\)
contains at least \(\lfloor m/2\rfloor\) factors of two.
Thus each summand with \(m>0\) has valuation at least
\((1-\beta)m>0\). This proves two-windows in this range
also. Together with the deep-scalar result, the remaining
range for this method is \(0\le\beta<2\lambda\) for odd \(p\)
and \(1\le\beta<2\) for \(p=2\), independent of \(a,d\).
