# Mixed-degree factorial windows from common aggregate saddle data

Let \(b_1,\ldots,b_s\) be distinct nonnegative integer
exponent columns, with total degrees \(d_i\) and coordinate
supports \(S_i\). Let \(G\) range over the degree blocks of
size at least two. For each \(G\), suppose there is a
coordinate set \(C\subseteq\bigcap_{i\in G}S_i\) satisfying
the following conditions.

1. Every foreign column is either supported on \(C\),
   or has a positive coordinate outside \(S_i\) for every
   \(i\in G\). Write \(J\) for the foreign columns
   supported on \(C\).
2. If \(i,k\in G\) have different restrictions
   \(b_{i,C}\ne b_{k,C}\), their coordinate supports
   are incomparable.
3. The following finite weighted measure is independent
   of \(i\in G\):
   \[
      \mu_i=\sum_{\nu\in C} b_{\nu i}\,
            \delta_{\,(b_{\nu j}/b_{\nu i})_{j\in J}} .
                                                               \tag{1}
   \]

**Theorem.** For all sufficiently large integers \(D\),
the factorial moments of
\(\sum_i c_iX^{Db_i}\) satisfy the actual support-size
window at every positive starting order, for arbitrary
coefficients over any characteristic-zero field.
The threshold is uniform on compact normalized shape
families satisfying the three conditions with positive
uniform lower bounds on nonzero entries, column distances,
and distinct degree gaps.

For example, take
\[
      X^2Y^2Z^2U,\qquad XYZ^4V,\qquad X^2Y^2Z^4
\]
with independent coefficients and common dilation.
At the repeated degree \(7\) block, the foreign degree
\(8\) monomial survives at both peaks. On \(C=\{X,Y,Z\}\)
the first restrictions are \((2,2,2)\) and \((1,1,4)\),
and both ratio measures are \(4\delta_1+2\delta_2\).
The variable \(U\) or \(V\) supplies the necessary face
separation. In particular the restrictions may have
different coordinate multisets. Equality in (1) also
holds when the restrictions are related by permutations
fixing each surviving foreign exponent vector.

The fixed-coordinate criterion is the special case
where the restrictions on \(C\) coincide throughout
each repeated block. The result below supplies its
common analytic saddle map under the weaker aggregate
condition (1). The analytic filtering argument then
applies to the same family of mode functions.

## A normalized foreign value depends only on the aggregate data

Fix \(G,C,J\), with \(d=d_G\), and put \(V=C^c\).
For each profile \(\mathbf r\) in the common measure
write \(M_{\mathbf r}=\mu_i(\{\mathbf r\})\).
Then \(d_C=\sum_{\mathbf r}M_{\mathbf r}\) is common
and \(d_V=d-d_C>0\). Indeed if \(d_V=0\), every
column in \(G\) has the same support \(C\). By the
second condition all restrictions would be identical,
contrary to distinctness.

For each \(j\in J\) the positive constant
\[
     L_j=\prod_{\nu\in C} b_{\nu i}^{\,b_{\nu j}}
                                                               \tag{2}
\]
is independent of \(i\). To see this, choose for each
nonzero profile \(\mathbf r\) the first index \(h(\mathbf r)\)
at which its entry is positive. For a coordinate with
profile \(\mathbf r\), one can use
\(b_{\nu i}=b_{\nu,h(\mathbf r)}/r_{h(\mathbf r)}\).
The selected index for that coordinate depends only
on which foreign exponents \(b_{\nu j}\) are positive,
so the product of the corresponding numerators in (2)
is fixed. In each profile the denominator is
\[
   r_{h(\mathbf r)}^{\sum_{\nu:\,\mathrm{profile}=\mathbf r} b_{\nu j}}
                      =r_{h(\mathbf r)}^{M_{\mathbf r}r_j},
\]
which is fixed by (1). Zero profiles contribute a
factor of one. This proves (2), also when \(J\) is
empty or contains a constant.

## The reduced critical equations

Work over \(\mathbb C\) and discard zero coefficients.
The hypotheses persist: deleting columns leaves the
face conditions in force, and deleting a foreign
survivor pushes (1) forward under coordinate projection.
Use the normalization of the
[fixed-coordinate theorem](SFC_MIXED_DEGREE_FIXED_COORDINATE_FOREIGN_TERMS.md#the-common-saddle-map-of-a-repeated-block).
Its peak value at start \(m\) from column \(i\) is
\[
 q_i(m)=c_i(mD)^{Dd_i}e^{-Dd_i}
                   \prod_{\nu} b_{\nu i}^{D b_{\nu i}},
      \qquad Q=\max_i |q_i(m)|.
\]
Here \(0^0=1\). The exponential-range argument
reduces to \(m>e^{\eta D}\) for a fixed \(\eta>0\).
The same eligible peaks and truncated boxes as in
the fixed-coordinate theorem apply.

At a peak \(i\in G\), temporarily keep only the
same-degree columns with restriction \(b_{i,C}\),
and the foreign survivors \(J\). Every temporarily
omitted monomial has a positive coordinate outside
\(S_i\), including same-degree columns with a
different restriction. They can be reinserted
with arbitrarily high inverse-power loss relative
to the local saddle contribution by the
[inactive-endpoint estimate](SFC_MIXED_DEGREE_INCOMPARABLE_COORDINATES.md#two-local-estimates).
The estimate uses vanishing in an inactive coordinate,
and also applies to a discarded monomial of the
same degree.

Denote the kept same-degree polynomial at normalized
time \(m\) by \(H_i(y)\). It factors as
\(y_C^{Db_{i,C}}\widetilde H_i(y_V)\), with the second
factor homogeneous of degree \(Dd_V\). Let
\(\theta_i\) be its local saddle at the peak, with
inactive coordinates zero. Then
\(\theta_{i,C}=b_{i,C}\) and
\(\sum_{\nu\in V}\theta_{i,\nu}=d_V\), as follows
from the active critical equations and Euler's
identity. Write \(w_i=e^{-Dd}H_i(\theta_i)\).
Its relative difference from \(q_i(m)/Q\) is
exponentially small by peak separation.

Define the constants and the time-dependent variables
\[
 \kappa_j=e^{-Dd}\frac{c_j(mD)^{Dd_j}}Q L_j^D,\qquad
 B_j=\frac{\kappa_j}{w_i}(t/m)^{D(d_j-d)} .
                                                               \tag{3}
\]
The constants are independent of \(i\) by (2).
For any peak \(i\in G\), their absolute values are
those of \(q_j(m)/Q\) times the normalized foreign
peak profile at \(b_i\), raised to power \(D\).
They are exponentially small by peak separation,
so the \(B_j\) remain exponentially small at eligible
peaks on the usual relative complex time disk.

Use the ansatz \(x_\nu=b_{\nu i}\alpha_{\mathbf r}\)
for \(\nu\in C\) of profile \(\mathbf r\). On \(V\)
multiply the old saddle \(\theta_{i,V}\) by a common
\(\alpha_V\). The exact critical-value change is the
stationary value of the common function
\[
 \begin{split}
 \mathcal L(\alpha;B)
  &= \log\!\left[
       \alpha_V^{Dd_V}
            \prod_{\mathbf r}\alpha_{\mathbf r}^{D M_{\mathbf r}}
       +\sum_{j\in J} B_j
            \prod_{\mathbf r}\alpha_{\mathbf r}^{D M_{\mathbf r}r_j}
               \right] \\
  &\quad-D\left[
       d_V(\alpha_V-1)+
       \sum_{\mathbf r}M_{\mathbf r}(\alpha_{\mathbf r}-1)
            \right].
 \end{split}                                                    \tag{4}
\]
All logarithms here are taken locally near \(B=0,\alpha=1\).
In fact the restricted monomial powers and the linear
exponential density give exactly (4). To verify
stationarity also in the original coordinates, let
\(T\) denote the first monomial inside the logarithm
and \(T_j\) the monomial multiplying \(B_j\).
At a coordinate \(\nu\in C\) of profile \(\mathbf r\),
the unrestricted partial derivative after substitution
is
\[
   \frac{D}{\alpha_{\mathbf r}}
     \frac{T+\sum_j r_j B_j T_j}{T+\sum_j B_jT_j}-D .
\]
It is the same at every coordinate of that profile,
and its product by \(M_{\mathbf r}\) is the reduced
partial derivative of (4). In \(V\) the foreign terms
are independent of the coordinates. Homogeneity and
the old active saddle make the full active derivatives
there equal, and their product by \(d_V\) is the
\(\alpha_V\) derivative. Thus all active critical
equations hold. The inactive coordinates remain at
the endpoints.

At \(B=0,\alpha=1\) the Hessian of \(\mathcal L/D\)
is diagonal with negative entries \(-d_V\) and
\(-M_{\mathbf r}\). This proves local uniqueness of
the critical point. On a polydisc of radius \(a/D\),
fixed derivatives of the displayed powers have their
usual \(O(D^k)\) bounds. The same critical-point
contraction as in the fixed-coordinate proof therefore
gives a common analytic value \(g_G(B)\), for
\(\sum_j|B_j|\le a'/D\), of size \(O(\sum_j|B_j|)\)
with bounded first derivatives. These estimates are
uniform in the allowed shape families. In particular
each nonzero \(M_{\mathbf r}\) is bounded below by one
of the positive exponent entries, and \(d_V\) also
has that lower bound.

## Transfer to the all-order window

At every eligible repeated peak the exact kept saddle
base is
\[
            (t/m)^{Dd} w_i\exp(g_G(B)),
                                                               \tag{5}
\]
with the constants (3) and the function (4) shared
throughout \(G\). This is formula (4) in the
[fixed-coordinate theorem](SFC_MIXED_DEGREE_FIXED_COORDINATE_FOREIGN_TERMS.md#universal-analytic-functions-for-the-filter).
Its subsequent analytic inverse constructions and
annular interpolation use precisely this common
base and these uniform estimates. The local amplitude
is analytic on a relative complex time disk of radius
comparable to \(1/D\), with leading coefficient
exponentially close to its positive pure value.
Indeed the translated active Hessian and the inactive
linear losses are small perturbations of the pure
ones, exactly as in the
[localized fixed-coordinate estimate](SFC_MIXED_DEGREE_FIXED_COORDINATE_FOREIGN_TERMS.md#universal-analytic-functions-for-the-filter).
Reinstating all omitted monomials does not affect
these mode estimates through any fixed inverse-power
order, by the endpoint estimate above. The eligible
singleton modes require no additional hypothesis.

Consequently the annular interpolation, derivative
transfer and final Vandermonde contradiction of the
[fixed-coordinate proof](SFC_MIXED_DEGREE_FIXED_COORDINATE_FOREIGN_TERMS.md#analytic-interpolation-across-degree-annuli)
apply unchanged. They handle arbitrary relative
coefficient magnitudes and possible confluent saddle
frequencies. The orders below the exponential range
are already covered by the
[all-coefficient exponential-range theorem](SFC_MIXED_DEGREE_TOP_BLOCK_TAIL_WINDOWS.md#all-coefficient-vectors-up-to-an-exponential-starting-order).
Embedding the finitely generated coefficient field
into \(\mathbb C\) proves the characteristic-zero
statement.

The window is the actual support-size window of
[Edo and van den Essen, Conjecture 2.8](https://arxiv.org/abs/1304.3956).

## What the common critical correction determines

For this comparison choose \(C\) to be the union of the
supports of the survivors. Under the hypotheses this is
allowed: removing unused coordinates deletes only the
zero-profile mass from (1) and can only merge restriction
classes in the second condition. Now allow (1) to differ
between peaks while keeping the first two conditions.
The derivation of (4) still gives the individual
critical corrections \(g_i(B)\), with the individual
weights in that formula. Add the mass \(d_{V,i}\) at
the zero profile to \(\mu_i\), and denote the extended
measure of total mass \(d\) by \(\widetilde\mu_i\).

At a critical point, put
\[
       \tau_j=\frac{B_j T_j}{T+\sum_h B_h T_h}
                 = B_j\partial_{B_j}g_i(B).
\]
The second equality follows by differentiating the
stationary value. The map \(B\mapsto\tau\) has identity
derivative at zero. The critical equations just proved
give
\[
 \alpha_{\mathbf r}=1+\sum_j\tau_j(r_j-1),\qquad
 \alpha_V=1-\sum_j\tau_j.
\]
Moreover \(\int r_j\,d\widetilde\mu_i=d_j\).
Using these formulas in (4) yields
\[
 \begin{split}
 g_i(B)
   &=D\int\log\bigl(1+\textstyle\sum_j\tau_j(r_j-1)\bigr)
                                  \,d\widetilde\mu_i(\mathbf r)\\
   &\quad-\log\bigl(1-\textstyle\sum_j\tau_j\bigr)
                  -D\sum_j(d_j-d)\tau_j .
 \end{split}                                                    \tag{6}
\]
Thus the normalized analytic germ \(g_i\) determines
the logarithmic integral in (6), after locally
inverting \(B\mapsto\tau\). Its Taylor coefficients,
together with the total mass \(d\), determine every
polynomial moment of the finite measure
\(\widetilde\mu_i\). Such moments determine the measure:
on the union of two finite sets of ratio profiles,
polynomial interpolation separates every individual
point. Therefore the normalized critical correction
is common exactly when these extended measures
coincide. With the present choice of \(C\), every
profile there is nonzero, so this is equivalent to (1).
This characterizes the phase-sharing step under the
first two conditions. Foreign survivors which differ
between peaks lie outside those conditions.
