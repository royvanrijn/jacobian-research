# An obstruction to a shared scalar saddle correction

Consider a fixed block \(G\) of common total degree \(d\)
and any set of foreign columns. At a saddle near a block
column \(i\), omit monomials with a positive exponent
in an inactive coordinate. This leaves its exact active
saddle value, since those monomials and their active
derivatives vanish on the face. Write its exponential
saddle base as \(U_i(\boldsymbol\epsilon)\) when the
foreign coefficients are \(\epsilon_j\), and write
\(w_i=U_i(0)\). Both bases include the exponential
density at the saddle. The block coefficients can be
chosen independently.

Ask for a common local analytic correction in the
scalar form
\[
       U_i(\boldsymbol\epsilon)=
            w_i\exp G_{\boldsymbol a}
                           (\boldsymbol\epsilon/w_i).
                                                               \tag{1}
\]
The function is allowed to depend arbitrarily on the
block coefficients \(\boldsymbol a\). Equality is
required between the eligible peaks for every choice
of the block coefficients, including coefficient
restrictions to any two labels.

**Proposition.** For sufficiently large dilation,
condition (1) forces each foreign monomial to survive
either at every block peak or at none. If \(C\) is the
union of the supports of the survivors, two block
columns with different restrictions on \(C\) must
have incomparable supports. With these conditions
imposed, (1) for arbitrary block coefficients holds
exactly when the aggregate ratio measures of the
[common aggregate saddle theorem](SFC_MIXED_DEGREE_AGGREGATED_SADDLES.md)
are equal. Thus relaxing the face restrictions while
keeping exactly the scalar phase form (1) and arbitrary
coefficients would require violating a necessary
first-derivative condition.

## The foreign coefficient derivative

It suffices first to restrict the block to two
columns \(b_i,b_k\). We use actual exponents
\(Db\) at an integer dilation \(D\). Normalize the
pure monomials to have values one at their own
peaks and set
\[
 X_h(y)=\frac{y^{Db_h}}{\prod_\nu b_{\nu h}^{D b_{\nu h}}},
        \qquad H_r=X_i+rX_k ,
\]
with the usual convention on zero entries. At the
first peak the normalized block coefficient ratio
is \(r\); at the second it is \(1/r\), after removing
the nonzero scalar \(r\). The saddle at the first
peak, denoted by \(\theta_i(r)\), is holomorphic on
\(|r|<e^{cD}\) for some \(c>0\) independent of \(D\).
Indeed the other monomial quotient is exponentially
small on the complex peak box, including derivatives
on a smaller box, by peak separation. Its product
with \(r\) is still small on this disk with \(c\)
sufficiently small. The active saddle contraction
with the uniform inverse pure Hessian gives a
holomorphic solution throughout. The same reasoning
gives \(\theta_k(s)\) holomorphic on \(|s|<e^{cD}\)
at the second peak. A term absent at the active
face gives a constant saddle. These are the active
constructions in the
[separated-support estimate](SFC_DILATED_SUPPORT_SADDLE_ENVELOPE.md#uniform-local-saddle-and-endpoint-estimate).

Since the retained block on each face is homogeneous,
the active critical equations give
\(\sum_\nu\theta_{i,\nu}=d\). The derivative of a
stationary value with respect to the coefficient
of a foreign monomial \(y^{D b_j}\), multiplied by
the unperturbed base, is therefore
\[
     w_i \partial_{\epsilon_j}\log U_i(0)
                       =e^{-Dd}\theta_i(r)^{D b_j}.
                                                               \tag{2}
\]
Indeed the active critical gradient is zero,
and the partial derivative of the logarithm of
the polynomial is \(y^{D b_j}/H_r(y)\).
Multiplication by \(w_i=e^{-Dd}H_r(\theta_i(r))\)
gives (2). If the foreign monomial is inactive,
the right side is zero. If it is active, it is
nonzero in the peak neighborhood.

By (1) the derivatives (2) must be equal wherever
both peaks are eligible. In particular this applies
to \(1/2<|r|<2\) for sufficiently large \(D\).
Thus activity of a foreign monomial does not
differ between the two peaks. Applying the same
argument to every two-label restriction proves
the common survivor assertion.

## Equality across the two peaks forces constant response

If a foreign monomial survives at both peaks, the
identity required by (2) is
\[
                         \theta_i(r)^{D b_j}
                                =\theta_k(1/r)^{D b_j}
                     \qquad (1/2<|r|<2).
                                                               \tag{3}
\]
The two sides are holomorphic respectively in a
disk about zero and on a disk about infinity.
The identity theorem propagates (3) to the
overlap of those disks, and they glue to a
holomorphic function on the Riemann sphere.
That function is constant.

Suppose \(S_k\subseteq S_i\). The term \(X_k\)
is then active at the \(i\)-th peak. There the
critical equations can be written exactly as
\[
        \theta_i(r)=b_i+\tau(r)(b_k-b_i),\qquad
        \tau(r)=\frac{r X_k(\theta_i(r))}
                       {X_i(\theta_i(r))+rX_k(\theta_i(r))}.
                                                               \tag{4}
\]
This is the coordinatewise equation obtained by
differentiating \(\log H_r(y)-D\sum y_\nu\).
It has \(\tau(0)=0\) and \(\tau'(0)=X_k(b_i)\ne0\).
Consequently the expression in (3) on the first
peak is constant only if the polynomial
\[
       \prod_{\nu\in S_j}
          \left(b_{\nu i}+u(b_{\nu k}-b_{\nu i})\right)^{D b_{\nu j}}
                                                               \tag{5}
\]
is constant in \(u\). Every varying factor in
(5) has a nonzero leading coefficient, and
its exponent is a positive integer. Thus no
varying factor can occur, and \(b_i,b_k\) agree
on \(S_j\). If \(S_i\subseteq S_k\), exchange
the labels and obtain the same conclusion.
Running over every surviving \(j\) gives equality
of the restrictions on \(C\). In particular
comparable same-degree supports cannot evade the
second face condition by letting the shared
function depend on the block coefficients.

## Returning to the aggregate condition

The necessity just proved permits partitioning
the block by its restriction to \(C\), as in the
aggregate theorem. All other classes are inactive
at its peak. Within its own class the block
factors as a fixed monomial on \(C\) times a
homogeneous polynomial in the other coordinates.
The unperturbed saddle therefore restricts to
\(b_{i,C}\) on \(C\). Its physical linear
response to a foreign coefficient in (1) is
\[
                         e^{-Dd}\prod_{\nu\in C}
                                           b_{\nu i}^{D b_{\nu j}} .
\]
Equality of these nonzero constants is forced by
(2). After absorbing them as common scales on the
foreign variables, the correction functions in
(1) are the normalized functions in the aggregate
theorem. The
[critical-correction converse](SFC_MIXED_DEGREE_AGGREGATED_SADDLES.md#what-the-common-critical-correction-determines)
shows that their equality forces the weighted
ratio measures to agree. Conversely equal ratio
measures give the common physical scales and
the shared correction directly by the aggregate
theorem. This completes the characterization.

The argument addresses exact scalar phase sharing
for arbitrary coefficients. Equivariant coefficient
families can instead identify whole local integrals
at sampled orders through coordinate changes, as
in the
[degree-orbit window theorem](SFC_EQUIVARIANT_DEGREE_ORBIT_WINDOWS.md).
