# Short terms and real height gaps in small-budget abc boxes

A small short term allows one tall letter to vary inside a narrow
valuation-degree interval at a positive saving. A fixed numerical gap
also permits finer local bounds on the degree measures of the two tall
letters. The first argument needs only a one-group congruence, and the
second uses real parametrization to count a ratio involving several
groups and a cofactor of larger height. Grouping on a variable mesh
also converts the height gap and the numerical weight at high degrees
into bounds on the nonatomic masses. [Variable-mesh estimates](ABC_VARIABLE_MESH_DEGREE_GROUPING.md)

## A strict saving when the short term is small

Let \(N_{\lambda,h}(X)\) count primitive positive triples \(a+b=c\le X\)
satisfying
\[
 \operatorname{rad}(abc)<c^\lambda,\qquad
 \min(a,b)\le c^h .
\]
There are absolute constants \(h_*,\sigma,\lambda_0>0\) such that for
\(0<\lambda<\lambda_0\) and every \(\varepsilon>0\),
\[
 N_{\lambda,h_*}(X)\ll_{\lambda,\varepsilon}
      X^{(1/3-\sigma)\lambda+\varepsilon}.
 \tag{SH1}
\]
In particular a limiting profile at \(d=1/3\) has positive short
numerical height. The proof gives existence of \(h_*,\sigma\) via
compactness.

### A movable window on one tall letter

Suppose \(d=1/3\) and some letter \(A\) has weak degree measure
\((1/3)\delta_0\). This includes numerical short height tending to zero.
Use the [scaled-degree measures](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees).
The three letters all have normalized root mass \(1/3\). Since \(A\)
has nonatomic mass zero, the weak degree measures of the two other
letters \(B,C\) are wholly nonatomic, each of mass \(1/3\) on the
positive scaled degrees.

Freeze the factors of \(A\) outside a shrinking zero window, at
vanishing normalized root-mass cost. Choose its upper endpoint tending
slowly to zero so the retained mass tends to \(1/3\). The remaining
factor of \(A\) has numerical height \(X^{o(1)}\), and the frozen
part enters the fixed coefficient of the signed ratio.

Choose \(\alpha,\beta>0\) so that \([\alpha,2\alpha]\) on \(B\) and \([\beta,2\beta]\) on \(C\)
both carry positive limiting mass. On \(B\) choose a positive-mass
interval \(I\subseteq[\alpha,2\alpha]\) of length \(\sigma_I\) and limiting mass
\(t\), with endpoints of measure zero. Such intervals can have
arbitrarily small length and positive mass tending to zero. Freeze all
of \(B\) outside \(I\), at cost \(1/3-t\). Group the roots in \(I\) as
one variable \(z\), with exponent \(q=\lfloor r/L\rfloor\) for the
lower endpoint \(r\) of \(I\). Thus
\[
 H(z)\le X^{(t+o(1))L},\qquad
 B=B_0z^q V_B,\qquad H(V_B)\le X^{\sigma_I t+o(1)},
 \qquad qL\longrightarrow r\ge\alpha.
\]
Here \(B_0\) includes the frozen factors. This decomposition uses
the integer difference between each true exponent and \(q\);
multiplying its scaled absolute value by the corresponding root mass
bounds its numerical contribution.

On \(C\) choose a modulus from a finite union of intervals contained
in \([\beta,2\beta]\) with continuity endpoints and limiting mass
\(\mu\in[t/16,t/8]\). For small \(t\), first choose a subset of the
nonatomic measure in the interior with mass strictly between
\(t/16\) and \(t/8\), then use inner and outer regularity to surround
a compact subset by a finite union of intervals in that range of
masses. The modulus freezes all the coordinate factors in the union.
It has weight \(w\in[\beta t/16,\beta t/4]\) and satisfies the condition of
\(k\)-fullness for \(kL\to\kappa=\beta/2\).

Write \(A=A_0 A_1\) with \(A_0\) frozen. The congruence from the
signed equation counts the pair \((z,V_B/A_1)\), with \(B_0/A_0\)
absorbed in the target. All variables are units by primitivity.
Since the numerical height of \(A_1\) tends to zero, the wider
factor has height
\(X^{E+o(1)}\) for \(E=\sigma_I t\). Use the
[unequal-height modulus estimate (UW1)](ABC_QUANTITATIVE_DEGREE_WINDOWS.md#a-product-bound-for-the-admissible-wider-height)
with \(m=1\), short height constant \(C_1=9>\kappa t/w\), and
amplification order \(\ell=16\). The fixed wider-height condition
there holds when \(\ell E/w<c_1/(\ell C_1)\). It suffices to choose
\(\sigma_I<\beta c_1/73728\), independently of the mass \(t\).
Also choose \(\sigma_I<\alpha/16\). All height conditions are used with
strict margins, allowing the errors from the limiting box passage.

The grouped root \(z\) is an integer. Its support may overlap the
wider factor. The exact-fiber term of (UW1) is allowed in this
configuration and costs at most \(E/r\le t/16\) in the normalized
exponent. The determinant term is at most
\[
 \frac{2w}{16\kappa}\le t/16 .
\]
At each fixed modulus the exact value represented by a pair
determines the complete ratio of the two opposite letters. Their
coprimality determines their signed values, and hence the third
letter. The number of representations of each recovered integer at
fixed truncation is divisor-many. Restoring the frozen choices
therefore gives
\[
 d\le 1/3-t+\mu+\frac{2w}{16\kappa}+\frac Er
     \le 1/3-\frac{3t}{4},
\]
a contradiction. The width can be chosen before the mass \(t\) and
does not depend on the size of that mass. In particular, numerical
height that is carried entirely by escaped first moment cannot
support a measure concentrated at scaled degree zero.

### Positive scales in every endpoint letter

There are absolute \(0<r_0<R_0\) and \(\gamma>0\) such that every
letter of every limiting profile at \(d=1/3\) has weak degree measure
of mass at least \(\gamma\) in \([r_0,R_0]\). Indeed if the masses
away from zero have no uniform positive lower bound on any fixed
positive cutoff, choose a sequence of endpoint profiles with a
letter of mass at most \(1/j\) above scale \(1/j\). After passing to a
subsequence the letter is fixed. Diagonalize the box sequences,
choosing continuity cutoffs at scales tending to zero and taking the
normalized exponents to tend to \(1/3\). The identified letter then
has weak degree measure \((1/3)\delta_0\), which is excluded above.
This gives uniform positive mass above some positive cutoff. The
uniform first-moment bound places at least half that mass below a
sufficiently large absolute \(R_0\).

### Passage to a uniform sector

Suppose (SH1) fails for every choice of positive \(h_*,\sigma,\lambda_0\).
Take \(h_*\to0\), \(\sigma\to0\) and \(\lambda_0\to0\). A diagonal
sequence of counterexamples gives, by the dyadic height shells and
anatomic box passage, numerical short height tending to zero and
normalized exponent \(d\to1/3\). Indeed the lower count exponent is
\((1/3-o(1))\lambda\), the total root mass is at most
\((1+o(1))\lambda\), and the [global limiting bound](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)
bounds \(d\) above by \(1/3+o(1)\). Thus the root mass is also
comparable with \(\lambda\). The anatomic errors can be taken
\(o(\lambda)\), absorbing dyadic choices and representation losses in
the arbitrary \(\varepsilon\). At each fixed budget the shell height
can be taken sufficiently large for its fixed truncation, as in that
limiting passage. The resulting limiting profile is the one just
excluded, proving (SH1).

## A real strip with unequal root heights

Use the [real algebraic parametrization theorem, Theorem 1](https://arxiv.org/pdf/2009.02082v1#page=1) and the narrow-coordinate degree split from the [bounded-root estimate (GW)](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-bounded-number-of-small-roots).

Fix \(m\ge1\), \(s,\tau>0\), and individual positive heights \(C_i\). Put
\[
 D_i=\max(1,C_i/s),\qquad P=\prod_{i=1}^m D_i,\qquad B_m=m^2+m.
\]
There is an absolute \(c_m>0\), depending only on \(m\), with the following property. Let \(L\) be sufficiently small, and at fixed \(L\) let the \(a_i\) be fixed signed integers. Consider positive rational tuples satisfying
\[
 H(x_i)\le X^{C_i L},\qquad H(V)\le X^E,\qquad
 \left|\frac{V\prod_i x_i^{a_i}}{K}-1\right|\ll X^{-\tau}
\]
for any \(K>0\), where \(E<c_m\tau/P\). Suppose at most \(F_0\) tuples in the sample have any fixed exact value of \(V\prod_i x_i^{a_i}\). Then, uniformly in \(K\), the count is
\[
 O_{m,s,\tau,\mathbf C,E,L,\mathbf a,\varepsilon}
 \left(F_0 X^{B_m s L+\varepsilon}\right).
 \tag{RUW1}
\]
The threshold in \(L\) depends on the fixed height and precision parameters and can be chosen independently of the \(a_i\). Errors tending to zero in the assumed height exponents can be absorbed using the strict wider-height margin and choices of slightly larger \(C_i\). The limit is taken with \(X\) sufficiently large at each fixed truncation.

### Normalization and the narrow-coordinate induction

Divide each variable by the lower endpoint of its dyadic box to put the whole sample in \([1,2]^{m+1}\). This raises the reduced rational heights by at most a factor of two. The normalization replaces \(K\) by a number bounded above and below in terms of the fixed exponents on any nonempty box. On the enlarged positive box \([1/2,3]^{m+1}\), the resulting monomial and its derivatives of each fixed order are therefore bounded uniformly in the target, with constants allowed to depend on the exponents. The number of dyadic boxes is absorbed in \(X^\varepsilon\).

Induct on the dimension \(j\) of a narrow algebraic locus containing the short-coordinate sample in the branch. Carry the original tuple restrictions along every branch, preserving \(F_0\). As in (GW), split geometric components. For a component not defined over \(\mathbb Q\), descend on its intersection with a distinct conjugate and the conjugates of that intersection. For geometrically integral rational components of sufficiently high degree, use a bounded-height finite birational projection of the projective closure to a hypersurface of the same degree. The projective height of the short coordinates is \(O(X^{L\sum_i C_i})\). For sufficiently high degree in terms of \(\mathbf C,s,m\), the hypersurface covering with bounded-degree auxiliary forms costs \(O(X^{sL+\varepsilon})\), and pulling back the auxiliaries descends the dimension. The lower-dimensional sets have degrees bounded at fixed parameters with unrestricted coefficient heights.

On a bounded-degree component, replace its equations by the bounded-degree rational interpolation kernel of its sample. This yields equations of height \(X^{O_{\mathbf C,s,m}(L)}\) whose common zero set is inside the old component. The coefficient estimate uses maximal minors of the sample evaluation matrix, since only a bounded set of narrow monomials and short rational coordinates are involved. The containment follows because bounded-degree separating equations for the old component belong to the base change of the rational kernel. If dimension drops, descend. Otherwise let \(\mathbf t\) be a dominant projection consisting of \(j\) of the original normalized short coordinates, chosen by algebraic independence. Obtain bounded-degree equations \(A_i\) for the trailing short coordinates by elimination, with heights \(X^{O_{\mathbf C,s,m}(L)}\). Process the exceptional locus where the eliminated factors or trailing derivatives vanish by dimension descent. On the remaining locus the \(A_i\) have independent gradients. This step is performed anew after descent, so the height exponents of equations obtained by earlier mixed projections do not propagate.

### The first mixed determinant

Choose \(M=\lceil A_m P\rceil\) with a sufficiently large fixed \(A_m\), and set
\[
 R=\left\lceil ME/(sL)+L^{-1/2}\right\rceil,\qquad
 R_i=\lfloor R/D_i\rfloor,\qquad \alpha=sLR .
\]
Thus \(ME\le\alpha=ME+o(1)\) and \(C_i L R_i\le\alpha\). On \((\mathbf t,V)\) use the rectangular support of degree at most \(R_i\) in each projected coordinate and at most \(M\) in \(V\). Its number \(D\) of columns satisfies \(D\ge A_m R^j\), since \(Q_i:=R_i+1\ge R/D_i\) and \(\prod_{i\ \mathrm{projected}}D_i\le P\). At a rational sample point clear the evaluation row by multiplying each projected denominator to power \(R_i\) and the denominator of \(V\) to power \(M\). The integer entries have height at most \(X^{(j+1)\alpha+o(1)}\), where fixed factors from dyadic normalization are absorbed for sufficiently large \(X\).

At an original sample point, correct \(V\) to the exact monomial level, keeping all short coordinates and the equations \(A_i=0\). The correction is \(O(X^{-\tau})\). The rank-full exact reference set in the enlarged bounded box has dimension \(j\), since \(\partial_V(V\prod_i x_i^{a_i})\) is nonzero and the \(A_i\) are short-coordinate equations. Use the real algebraic charts and parameter cells of side \(\delta=X^{-sL}\). The number of cells is \(O(X^{jsL})\) at fixed \(L\). Their implied constant is uniform in the coefficients and target by the cited theorem.

For \(N=\lceil A'_m R^j\rceil\) with \(A'_m\) large in terms of \(m\), the sum of the first \(N\) Taylor orders in \(j\) variables is greater than \(2(j+2)RN\), and their largest order is \(O_m(R)\). Choose \(A_m\) large enough for \(N\le D/2\). At each chart cell, expand the monomial evaluation vectors at the corrected points in Taylor orders up to the last needed order. Each repeated Taylor vector contributes zero, and the original-to-corrected error can be given the last Taylor weight if \(O_m(R)sL<\tau\). The displayed hypothesis on \(E\), with \(c_m\) small in terms of \(m\), ensures this strict inequality for small \(L\). Denominator clearing then gives each \(N\)-row integer determinant absolute value less than one for sufficiently large \(X\), using the Taylor sum and the entry height bound. Hence the rank is below \(N\). Pigeonholing the evaluations on a maximal independent row set gives a nonzero integer polynomial \(G(\mathbf t,V)\) on the same rectangular support with coefficient height \(X^{O_m(\alpha)}\). The bound follows because a coefficient box of radius \(T\) has order \(T^D\) elements, while its evaluations on \(r\le D/2\) independent rows have order at most \((O_D(T)X^{(j+1)\alpha+o(1)})^r\) possibilities.

At least one \(V\)-coefficient of \(G\) is a nonzero polynomial in the dominant coordinates \(\mathbf t\). Its zero locus bounds the simultaneous vanishing locus of the coefficients on the narrow component, so the image has smaller dimension. Descend there. If \(G\) is independent of \(V\) the same holds on its vanishing locus. Otherwise split the remaining pairs by the order of the \(V\)-root of \(G\). On that subbranch write \(G_r=\partial_V^{r-1}G\), so \(G_r=0\) and \(\partial_VG_r\ne0\). The rows of \((A_i,G_r)\) have full rank. The order \(r\le M\) introduces only constants fixed before \(X\) tends to infinity.

### Collision-energy optimization test (open)

[Long, *Conjugate fibres, certificate rigidity, and irrationality exponents*, Lemma
“Collision-energy optimization” and the water-filling stability
remark](https://github.com/long-mathematics/openai-math-extensions/blob/bc45a6906e11db8ca737f803058120b1b88fdf3f/papers/04-conjugate-fibres-rigidity/conjugate-fibres-rigidity.tex)
jointly prices approximation loss and repeated Taylor row types. Its
completion of squares retains a squared-distance remainder from the
optimal occupancy profile. The result concerns its own interpolation
matrix; it does not improve (RUW1) without a new row-type accounting for
the present mixed determinant.

The bounded test is to fix one dimension \(j\) and the rectangular
support \((R_i,M)\) above, keep each row's original-to-corrected error
and each repeated Taylor-type saving separately, and optimize their
combined determinant exponent. Compare the resulting admissible
cofactor-width threshold \(E\) with \(c_m\tau/P\) under exactly the same
denominator-clearing exponent, rank-correction condition, and second
narrow determinant. Success means a proved strict increase in \(E\) for
fixed parameters after all those costs; a gain in only the analytic
minor estimate is insufficient. Long's analytic bound holds for every
maximal minor, but arithmetic use of the full vector of minors is a
separate open step. Its aggregate scalar optimization also has a stated
limitation, so no abc consequence is inferred here. This is a targeted
primary-source reading, not a replay of the source proof or a new RUW1
bound.

### The rank correction and the purely narrow determinant

In the positive torus, close the resulting finite-wide-fiber strata algebraically. The closure of each stratum has dimension at most \(j\). A component with constant exact monomial contributes at most \(F_0\) original tuples. On any other component some maximal minor \(J\) of the rows from \(A_i,G_r\) and the cleared logarithmic differential of \(V\prod x_i^{a_i}\) is generically nonzero. This follows because the vanishing of every such minor would make the logarithmic differential annihilate the tangent space of the component, forcing the monomial to be constant in characteristic zero. Choose a fixed bounded-integer linear combination (with bounds allowed to depend on the algebraic degrees) of the minors nonzero generically on all of these components. The \(J=0\) loci descend after projection to the narrow variables, since each proper exceptional set has dimension at most \(j-1\). Their coefficient heights are immaterial at the next narrow-coordinate reset.

At an original point with \(J\ne0\), the numerator and denominator heights of its rational evaluation are
\[
 X^{O_m(\alpha)+O_{\mathbf C,s,m}(L)}.
\]
In fact each occurrence of \(G_r\) or any of its derivatives retains the individual degree bounds, so its denominator after evaluation costs at most \(\sum R_iC_i L+ME\le(j+1)\alpha\). The cleared logarithmic row has bounded degree and the low-degree \(A_i\) add only \(O_{\mathbf C,s,m}(L)\) in the exponent. Coefficients involving the \(a_i\) are fixed at fixed \(L\) and enter the multiplicative constants. It follows that a selected maximal Jacobian minor of the system \((A_i,G_r,\phi-K)\), where \(\phi\) is the normalized monomial, has absolute value at least \(X^{-O_m(\alpha)-O_{\mathbf C,s,m}(L)}\) on the sample. The derivative row for \(\phi-K\) differs from the logarithmic row by a scalar bounded above and below at such points. Hold the unselected coordinates fixed. The derivative inverse and the second derivatives of this system on the enlarged box are bounded by \(X^{O_m(\alpha)+O_{\mathbf C,s,m}(L)}\). These bounds follow from the lower bound on the minor, the adjugate formula, and the polynomial coefficient heights. The absolute residual at the original point is \(O(X^{-\tau})\).

For clarity, with derivative inverse norm at most \(U\), Hessian norm at most \(H\), and residual norm at most \(\varepsilon_0\), the map \(z\mapsto z-A^{-1}F(z)\) sends the ball of radius \(2U\varepsilon_0\) around the original point to itself and is a contraction when \(2U^2H\varepsilon_0<1/2\). This follows from Taylor's theorem and the derivative bound on the ball. Its fixed point is an exact solution and still has full rank. Here the fixed smallness condition on \(E\) ensures the inequality, since \(\alpha=ME+o(1)\) and every multiplicative constant in the exponent from the mixed polynomial depends only on \(m\). The terms \(O_{\mathbf C,s,m}(L)\) are absorbed by taking \(L\) small. The correction distance is \(X^{-\tau+O_m(\alpha)+O_{\mathbf C,s,m}(L)}\).

The exact reference set now has dimension \(j-1\). Apply real algebraic charts and cells of side \(X^{-sL}\) to it, and form the second determinant using only monomials of bounded degree in the \(j\) narrow coordinates \(\mathbf t\). Choose the degree large enough in terms of \(\mathbf C,s,m\) so the chart-dimension Taylor gain exceeds the denominator-cleared integer determinant height. The needed Taylor order is independent of \(L,R\); hence the correction distance is smaller than the needed Taylor error for small \(L\). A nonzero auxiliary in the dominant coordinates \(\mathbf t\) results and descends the narrow dimension, with \(O(X^{(j-1)sL})\) records. When \(j=1\), the full-rank exact references are isolated and bounded in number by Bézout at fixed exponents. Each correction neighborhood contains at most one original rational tuple: its radius is smaller than each spacing \(X^{-2C_iL}\) and \(X^{-2E}\) by the same strict smallness condition. At \(j=0\), for each of boundedly many short points the wide coordinate lies in an interval of length \(O(X^{-\tau})\), which also has at most one such rational.

Each positive dimension \(j\) has a covering cost at most \(X^{sL+\varepsilon}\), or a pair of chart costs totaling \(X^{(2j-1)sL+\varepsilon}\), before descent. Intermediate exceptional branches carry no larger cost. Counts on constant components carry \(F_0\). Summing over dimensions and absorbing the fixed-truncation \(\varepsilon\)-losses proves (RUW1).

## Consequence for the two tall measures

Consider a limiting profile at \(d=1/3\). Let \(U,V\) be the tall
pair, with a positive gap to the numerical height of the short letter,
and let \(\tau>0\) be strictly below that gap. Put
\(\Delta=c_U+c_V-1/3\). Suppose the two tall measures have finitely many
positive atomic scales, totaling \(p\) across both letters, counting a
shared scale once per letter. For \(p\ge2\), the nonatomic restrictions
on compact positive ranges away from the atomic scales obey
\[
 \nu_{T,\mathrm{na}}(I)\le
 \max\left\{2\Delta,\ C\left(\frac{|I|}{\tau}\right)^{1/(p-1)}\right\}
 \qquad (T=U,V).
 \tag{RUW3}
\]
The constant \(C\) may depend on the profile and the compact ranges.
For \(p=0,1\), the bound is \(2\Delta\) on all intervals of
sufficiently small length in those ranges.

For an interval \(I\) with continuity endpoints, mass \(t>0\) and
length \(\sigma\), free the zero windows and all shrinking atomic
windows on the tall pair, and free \(I\) on its letter. Fix the
remaining coordinates at normalized cost \(1/3+\Delta-t\). Group
each freed positive window around a representative exponent
\(\lfloor r_i/L\rfloor\). The wider numerical height exponent is at
most \(\sigma t+o(1)\), since the atomic and zero windows have
vanishing scaled widths. The group roots have height exponents
\((v_i+o(1))L\), with \(v_i\) the respective masses. The complete
tall ratio, divided by its frozen factor, lies in the strip (RUW1).
An exact ratio recovers the primitive triple, with divisor-many
representations at fixed truncation.

Take \(m=p+1\) and \(s=t/(2B_m)\). The \(p\) atomic group height
constants are bounded above, and the interval group has height constant
at most \(2t\). Thus the product in (RUW1) is \(O(t^{-p})\).
If \(\sigma<c\tau t^{p-1}\) for a sufficiently small constant \(c>0\),
the normalized strip cost is \(t/2\), which gives
\(d\le 1/3+\Delta-t/2\). Every \(t>2\Delta\) must therefore satisfy
\(\sigma\ge c\tau t^{p-1}\). This proves (RUW3) for \(p\ge2\),
and for \(p=0,1\) it excludes \(t>2\Delta\) on sufficiently short
intervals. General endpoints follow by enlargement to continuity
endpoints.

The [variable-mesh bound](ABC_VARIABLE_MESH_DEGREE_GROUPING.md#using-a-numerical-height-gap)
gives a uniform lower bound for \(\Delta\) in terms of the gap,
with arbitrary atomic spectra. Estimate (RUW3) also bounds the
mass of individual intervals when \(\Delta>0\).
