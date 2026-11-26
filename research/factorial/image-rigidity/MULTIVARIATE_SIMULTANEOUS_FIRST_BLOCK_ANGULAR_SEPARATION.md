# All-degree angular separation from simultaneous first-block factors

In an ordinary two-monomial profile, omit equal rows and select a
nonempty group of unequal rows. Let \(C\ge1\) be an integer and
assign each selected row an integer multiplier \(c_i\mid C\).
For that row write \(M_i=\max(a_i,b_i)\), \(d_i=|a_i-b_i|\)
and \(g_i=c_i d_i\). Define
\[
 M=\max_i c_iM_i,\qquad X=\max_i g_i,\qquad e_i=M-c_iM_i .
\]
The maxima here are over the selected group. Fix an integer
\(B\ge0\) with \(0\le e_i\le BX\), and let \(J\ge4\) be at least
the [absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
Put
\[
 F=\max_i(e_i+g_i),\qquad
 G=\max_i(e_i+(1-1/J)g_i),\qquad H=J(B+1).
\]
Choose \(G<T\le F\) with \(T-G\ge\gamma X\) for some \(\gamma>0\).
Let \(R\) be the number of selected rows with \(e_i+g_i\ge T\).
One can always use \(T=F\) when \(F-G\ge\gamma X\).
The selected group may have arbitrary gap distributions, row counts
and orientations.

Among the retained rows define the log height \(\mathcal L\) by
summing \(d_i\log(JM_i)\). The sum omits any rows with \(JM_i\le HX\).
It also omits rows for which a multiplier \(c_i\mid C\) puts both
scaled endpoints in either of the intervals
\([M-G,M]\) and \([M-(B+1)X,M-F]\).

**Theorem.** Suppose
\[
 X\ge\max\{144C^4(B+1),\,2C/\gamma\},\qquad
 \log X\ge144CH/\gamma,\qquad
 \mathcal L\le\frac{R\gamma X\log M}{64CJ^3}.                  \tag{1}
\]
Then the original profile has disjoint consecutive root arguments
in every degree. In particular it satisfies SFC(2) over every
characteristic-zero field. The retained budget grows with the number
\(R\) of ramps covering the window. Under total absolute curvature
\(S<12\), one prime yields binomial reductions simultaneously at all
orders \(2\le j<J\).

For example, \(0\le e_i\le X/(2J)\) automatically gives
\(F-G\ge X/(2J)\), since \(F\ge X\) and
\(G\le X/(2J)+(1-1/J)X\). Every group with a common scaled upper
endpoint has \(F-G=X/J\), even when its gaps and multiplicities
are arbitrarily irregular. The divisor \(C\) can be any fixed positive integer. In either case an arbitrary fixed
retained core eventually contributes nothing as \(X\) increases
with \(B,C,J\) fixed. The numerical value of \(J\) is unevaluated.

The window also exists for arbitrary offsets of the form
\[
                        e_i\le b g_i+\delta X,
                       \qquad b\ge0,\quad 0\le\delta<1 .
\]
Indeed \(g_i\ge(e_i+g_i-\delta X)/(b+1)\), so
\[
 G\le \left(1-\frac1{J(b+1)}\right)F+
                               \frac{\delta X}{J(b+1)}.
\]
Thus \(F-G\ge(1-\delta)X/[J(b+1)]\).
This covers bands whose upper offsets are bounded by a multiple
of each gap, with no approximation of the gap counts.

## A first-block factor at all low orders

If \(S\ge12\), use the
[absolute curvature theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md).
For \(S<12\), a row with scaled gap \(X\) gives
\(M>X^2/(12C)\). Take every integer \(u\in[G,T)\) such that
\(C\mid M-u\). Put \(y=M-u\) and \(w=y/C\).
The integers \(w\) form a step-one progression.

At a fixed \(2\le j<J\), strip each row's scalar factorial
\((j\min(a_i,b_i))!\). In selected row \(i\), the candidate factor
is the integer \(jy/c_i\). Write \(v_i=e_i+g_i-u\).
If \(v_i\le0\), the candidate lies at or below the
stripped ramp's lower endpoint. If \(v_i>0\), then
\(jv_i\le g_i\) because \(u\ge G\) and \(j<J\).
Thus \(jy/c_i\) is in the first block of length \(d_i\) above
the smaller endpoint at every order \(2\le j<J\).
At least \(R\) selected rows are active, since \(u<T\).
If \(r_+,r_-\) count the active increasing and decreasing rows,
the number of marked literals present in the coefficient array is
\[
        r_-,\quad r_-+r_+,\ \ldots,\ r_-+r_+,\quad r_+ .
\]
The counts do not depend on \(j\), apart from the number of interior
entries. At \(0<k<j\), the excess of the mark count over its
endpoint chord is
\[
              (1-k/j)r_+ +(k/j)r_-\ge R/j.                 \tag{2}
\]

For a prime \(p>HX\) dividing \(w\), put \(E=v_p(w)\).
The factor \(jy/c_i\), when present, has valuation \(E\).
Every other factor in the full selected ramp has difference from
it, after multiplication by \(c_i\), equal to
\[
                       j(u-e_i)-c_i h,\qquad 0\le h<jd_i ,
\]
with nonzero magnitude at most \(J(B+1)X=HX\).
Thus all nonmarked factors in selected rows are units.
The retained rows below the prime cutoff have unit factors.
A retained row in the first permitted scaled interval has every
scaled ramp factor strictly greater than \(j(M-G)\ge jy\).
One in the second has every scaled factor at most \(j(M-F)<jy\).
The same offset bound shows that these rows also contribute units.

Let \(D_j\) be the product of all full stripped ramps of the
retained rows counted in \(\mathcal L\). The binomial coefficients
are units, since \(p>j\). All residual ordinates of the stripped
degree-\(j\) moment lie in \([0,v_p(D_j)]\). By (2), each
valuation in the interior is strictly above the endpoint chord
if \(RE>jv_p(D_j)\).

Consequently, if \(RE>\max_{2\le j<J}jv_p(D_j)\), every low-order
moment has, after a scalar and a local change of variable,
reduction \(1+z^j\). To obtain the change, take a \(j\)-th root
of the ratio of the endpoint coefficients and divide by the
constant coefficient. All roots are units in this local
coordinate and have distinct reductions, since \(p>j\).
Distinct-root quotients are unchanged by the rescaling and reduce
to nontrivial \(j\)-th roots of unity.
The [binomial quotient comparison](MULTIVARIATE_UNIVERSAL_OPPOSED_PADIC_ANGLE_MARKING.md#binomial-reductions-separate-root-quotients)
then excludes any common distinct-root quotient between coprime
orders in the indicated range. Two nonreal roots with the same
argument would give the same quotient with their complex
conjugates. The even moment has no real zero, so adjacent
arguments are disjoint. The pair of orders one and two follows
by the same even positivity, and all larger terminals follow
by the absolute cutoff.

## Existence of a prime for all low orders

Write \(s\) for the number of marks \(w\).
Since \(F-G\le X/J\) and \(X\ge2C/\gamma\),
\[
                   \gamma X/(2C)\le s\le X.
\]
The bound on \(X\) gives \(M\ge2(B+1)X\) and \(M>(2C)^6\),
so every mark exceeds \(M/(2C)\) and has logarithm
greater than \((5/6)\log M\).
It also gives \(\log s\le\log X\le(2/3)\log M\).

Suppose every prime \(p>HX\) dividing any mark has
\(Rv_p(w)\le\max_{2\le j<J}jv_p(D_j)\).
Since the marks have span less than \(HX\), each such prime
divides at most one mark. The large-prime part of the product
of the marks, raised to the power \(R\), therefore divides
\[
                        \prod_{2\le j<J} (D_j)_{>HX}^{\,j}.
\]
The logarithm of this product is at most \(J^3\mathcal L\),
since degree \(j\) uses \(jd_i\) factors per retained row
and each factor is at most \(JM_i\).

As in the [progression product comparison](MULTIVARIATE_LARGE_FACTORIAL_MARK_PRIME_BUDGET.md#a-product-bound-at-a-possible-common-root),
remove one term of maximal valuation for each smaller prime
that divides a mark. The remaining small-prime part divides
\((s-1)!\). The elementary prime bound
\(\pi(HX)\le9HX/\log(HX)\) and (1) give
\(\pi(HX)\le\gamma X/(16C)\le s/8\).
Taking logarithms on the remaining product gives the
incompatible bounds
\[
   (7/8)(5/6)\,s\log M
       <\log(\text{remaining product})
       \le \log((s-1)!)+J^3\mathcal L/R
       \le (2/3+1/32)\,s\log M.
\]
This proves existence of the required prime.
The argument also yields an exact finite test by factoring the
marks and comparing \(Rv_p(w)\) with the actual retained valuations.
The factorial moment coefficients are rational, so the resulting
coprimality extends by Bézout to any characteristic-zero field.
