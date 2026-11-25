# Heights of rational transfers and the coarse Fermat curve

Two tests concern different changes of variables in the small-budget
count. A nonconstant rational chart on an affordable fibre retains
a positive numerical secant scale after any fixed support is removed.
For generalized Fermat equations, the normalized weighted projective
curve can have genus zero at arbitrarily large exact orders. Its
rational points can be realized through common multiplication. The
minimal-repair count requires the primitive quotient and its prescribed
height cutoff.

## A uniform small-output bound for rational charts

Use the balanced fixed-budget setting of the projective packing gap.
In particular the primitive rows have bounded numerical heights at
most \(X\) and complete-letter radical caps bounded in base \(n=X^L\).
They have distinct targets and fixed affordable factors \(F,G,D\)
as in (PS1)–(PS3), with at least one of \(\log_X F,\log_X G\) bounded
below by a fixed positive constant. Use their primitive residual
pairs \((u_j,v_j)\) as projective coordinates. Let \(f\) be a
nonconstant rational function on these coordinates, with degree
bounded by a fixed integer and represented by homogeneous integer
forms of polynomial coefficient height. Both may vary with \(X\).
For every fixed \(K,\rho>0\) and all sufficiently small fixed \(L\),
uniformly in the stated data, the number of rows with defined finite
output and rational height \(H(f(u_j/v_j))\le n^K\) is at most
\(n^{\rho+o_X(1)}\). The height threshold and inner loss are at
fixed \(L\). The small-budget threshold depends on \(\rho,K\), the
radical caps and the initial affordability margin. The degree bound
and the polynomial coefficient-height power may be chosen at each
fixed \(L\), before the height limit. They affect the inner threshold
and errors. The numerical gap in (PS3) depends on the fibre exponent
and the affordability margin, and the bounded refinement in (PS2)
works for any fixed polynomial support bound at that budget.
[Uniform projective packing gap](ABC_AFFORDABLE_FIBRE_SECANT_VALUATION_WEIGHTS.md#setup-and-the-projective-packing-gap)

In fact, start with any subfamily of exponent at least \(\rho\) on
which the output is defined and finite. Apply (PS2)–(PS3) to that
subfamily after refining to dyadic cofactor intervals at subpower
cost. It gives a retained subfamily of size \(N_*\), common modulus
\(M\), effective angular scale \(h\ge\nu>0\), and a fixed
\(\delta>0\) in the two exceptional fractions in (PS2). To choose
the parameters, take the projective gap for class exponent \(\rho/4\),
use a thinning loss less than \(\rho/4\), and then take the numerical
tolerance \(\eta<\nu/4\). The resulting class exponent remains greater
than \(\rho/4\), so (PS3) uses the uniform numerical bound \(h\ge\nu\),
independently of \(L\). The polynomial support cap in (PS2) can be
chosen to include the fixed degree and coefficient height bounds.
[Angle and projective peeling (PS2)–(PS3)](ABC_AFFORDABLE_FIBRE_PROJECTIVE_SECANT_PEELING.md#simultaneous-peeling)

Discard the bounded number of roots of a chosen raw denominator.
Fix any retained anchor \(i\) and any additional prime pool of
polynomial product height. Write \(f\) as \(A(u,v)/B(u,v)\) with
homogeneous integer forms of the same degree, and choose a
representation with \(B(u_i,v_i)\ne0\). Use as the forbidden pool
the prescribed primes together with the primes of \(M B(u_i,v_i)\).
The numerical cap on this pool is uniform: \(M\) divides every
nonzero retained secant, and the value of \(B\) has polynomial
height. For all but \(O(N_* n^{-\delta})\) queries \(j\), the
angular part of (PS2) gives
\[
 \frac{|u_i v_j-u_j v_i|}{M}\gg X^{h-\eta}.
\]
The other part of (PS2) bounds the numerical weight of this
quotient on the forbidden pool by \(\eta\). The full secant powers
on primes outside the pool therefore have product
\(Q\gg X^{h-2\eta}\), which is at least \(X^{\nu/3}\) for large
\(X\). For each \(p^e\mid Q\), with \(e\) the full valuation,
the two primitive cofactor pairs are congruent projectively modulo
\(p^e\): a unit coordinate of the anchor gives the scaling, which
is a unit by primitivity of the query. Evaluating the two forms
scales their values by the same unit power. The second form is a
unit at the anchor and hence at the query. The two rational output
values are consequently congruent modulo \(p^e\), also when reduced
to lowest terms. This proves that their output cross determinant
has outside-pool numerical weight at least \(\nu/3\) whenever the
outputs are distinct. Query repetitions of the anchor output have
bounded multiplicity by the nonconstant rational function and the
distinct projective input directions.

If both rational heights are at most \(n^K\), the absolute output
cross determinant is at most \(2n^{2K}\). Take \(L\) sufficiently
small that \(2K L<\nu/3\). The congruence forces equality of the
output values on more than the bounded possible number of queries.
This contradiction proves the claimed small-output count. The proof
also gives the outside-pool numerical lower bound for output cross
determinants without imposing any cap on the output height.

The bound is uniform enough for piecewise rational charts. With at
most \(n^{\theta+o_X(1)}\) possible nonconstant charts of the stated
degree and polynomial coefficient-height bounds, the union of
their small-output row sets has size at most
\(n^{\theta+\rho+o_X(1)}\), regardless of how the chart is selected
at each row. An affordable subfamily of exponent \(r>\theta\) thus
cannot be covered in this way at all sufficiently small fixed \(L\),
by choosing \(0<\rho<r-\theta\).

### Tuple charts and the pure-term pattern

Let the fibre size be \(N=n^{r+o_X(1)}\), with \(r\ge r_0>0\), and
let \(P\) be a rational function of a fixed number \(m\ge2\) of row
directions of bounded degree and polynomial coefficient height.
Suppose it depends on some slot. Fix that slot as the query and
fix the complementary directions. A nonzero coefficient of the
numerator of the derivative in the query variable is a nonzero
polynomial in the complementary directions. Outside its zero set
the specialized rational function is defined generically and
nonconstant. In a Cartesian product of distinct input directions,
the zero set has \(O_P(N^{m-2})\) complementary tuples, by induction
using the univariate root bound. The implied constant here depends
only on the degree and number of slots. Undefined points of the
specialization are excluded from the count. A complement on which
the specialization has common polynomial factors may simply keep
the raw forms and restrict to points where their denominator is
nonzero. If the specialized rational function is instead evaluated
after cancellation at removable points, there are at most a
bounded number of such queries per nonexceptional complement.
The original tuple collection can be arbitrary, since we bound the
whole query slice and sum over complements.

After replacing the query direction by \(Gv/(Fu)\) or its
reciprocal, the specialized map has polynomial coefficient height
uniformly. The one-slot bound gives a total of
\[
                         O\bigl(N^{m-1}n^{\rho+o_X(1)}\bigr)
\]
defined outputs of height at most \(n^K\), counted with input
multiplicity. Take \(0<\rho<r_0/2\). In the pure-term two-mixed-factor
pattern, the high-output bound already applies to eligible tuples
at radical cap \(n^b\), \(b\le3r\), for sufficiently large fixed
\(K\). Thus the eligible tuples have a power-smaller exponent
across all output heights. This includes genuinely mixed opposite
orientations and charts in three or more slots. The conclusion
is conditional on the output radical cap.
[High-output two-mixed-factor count](ABC_PURE_TERM_DENOMINATOR_RECOVERY.md#outputs-above-a-fixed-power-of-the-radical-scale)

For example, with \(x=V_1/U_1\) and \(y=U_2/V_2\), the chart
\(P(x,y)=2+xy\) has mixed terms \(2U_1V_2+V_1U_2\) and
\(-U_1V_2-V_1U_2\), and pure term \(U_1V_2\). This example
has a separate direct height proof: refine the actual tuple graph
by the denominator-recovery procedure, using a fixed numerical
threshold instead of the radical threshold and recording both
free letters in each slot. A bad slot is one whose residual free
letters have numerical weight above the threshold on complementary
letter primes. The bounded potential is now the total fixed
numerical weight. The positive residual-height lemma gives one
slot residual weight bounded away from zero on the retained
dyadic cell. Each of its prime powers outside the recorded pool
and complementary supports survives in the numerator or denominator
of \(xy\), hence \(H(xy)\) is bounded below by a positive power
of \(X\). The same argument works for \(g(\prod_i x_i^{e_i})\)
with all \(e_i\ne0\) and fixed nonconstant rational \(g\), by
the usual fixed rational-map height comparison.
[Tuple refinement](ABC_PURE_TERM_DENOMINATOR_RECOVERY.md#a-bounded-refinement-of-the-eligible-tuples)
· [Positive height after division](ABC_AFFORDABLE_FIBRE_DETERMINANT_COMPLEXITY.md#positive-height-after-division)

## Exact orders and the genus of the coarse Fermat curve


The calculation below restates the effective-order formula (EC5)
and applies it to the least-repair question.
[Effective-order cover](ABC_EXACT_POWER_REPAIR_COEFFICIENTS.md#the-effective-orders-determine-the-weighted-fermat-curve)

Fix nonzero rational coefficients \(A,B,C\), positive integers
\(a,b,c\ge2\), and \(J=\operatorname{lcm}(a,b,c)\). Let \(W\) be
the normalization of \(Ax^a+By^b=Cz^c\) in projective weights
\((J/a,J/b,J/c)\). It maps to the \(t\)-line by \(t=Ax^a/(Cz^c)\).
Its geometric degree and three ramification indices are
\[
 D=\frac{abc}{J},\qquad h_a=\gcd(a,\operatorname{lcm}(b,c)),\quad
 h_b=\gcd(b,\operatorname{lcm}(a,c)),\quad
 h_c=\gcd(c,\operatorname{lcm}(a,b)).
\]
Its genus satisfies
\[
          2g(W)-2=D\left(1-h_a^{-1}-h_b^{-1}-h_c^{-1}\right).
\]

Indeed, over \(\mathbb C\) and away from the three branch values,
gauge \(z=1\). The two independent root equations for \(x,y\) have
geometric monodromy \(\mu_a\times\mu_b\). The first generator is
monodromy around zero and the second around one, since the two
equations have separate simple zeros. The residual gauge group is
\(\mu_{J/c}\), acting freely when both coordinates are nonzero because
the three weights have gcd one. Its quotient is geometrically
connected of degree \(D\). In this quotient the order of the first
monodromy generator is \(a\) divided by the size of the residual
gauge subgroup fixing both \(y\) and \(z\). The size of that subgroup
is \(\gcd(J/b,J/c)=J/\operatorname{lcm}(b,c)\); its action on \(x\)
is faithful. This gives \(h_a\). Gauge the other slots in turn for
the other two indices. The cover is Galois over \(\mathbb C\) because
the original monodromy is abelian, so each displayed index is uniform
above its branch value. Riemann–Hurwitz gives the formula.

If \(a,b,c\) are pairwise coprime the degree is one and \(W\) is the
rational line over \(\mathbb Q\). If only one pair has a nontrivial
gcd \(d\), the degree is \(d\) and the indices are \((d,d,1)\) in
the corresponding orientation, again giving genus zero. Both cases
allow orders of arbitrary size.

Pairwise coprime orders make the common-multiplier issue explicit.
For unit coefficients take any positive primitive \(U+V=T\). For
each prime choose the smallest nonnegative \(j_p\) such that
\[
 j_p\equiv-v_p(U)\pmod a,\qquad j_p\equiv-v_p(V)\pmod b,
 \qquad j_p\equiv-v_p(T)\pmod c.
\]
The Chinese remainder theorem gives \(j_p<J\) and \(j_p=0\) off
the three supports. Then \(q=\prod_p p^{j_p}\) makes \(qU,qV,qT\)
exact powers of the respective orders, and is the least such common
multiplier by divisibility. If \(p\) has positive order in \(U\)
that is not a multiple of \(a\), then \(j_p\) is a nonzero multiple
of \(bc\). The analogous assertions hold on the other two letters.
Thus every rational \(t\) in \((0,1)\) has a positive integral lift
in the unit-coefficient case. Pure-power lifting by a common
multiplier can be much more expensive than filling the missing
valuations on the one prescribed letter, and produces no bound on
the admissible primitive quotients.

Darmon and Granville obtain finiteness for a fixed hyperbolic
signature and fixed nonzero coefficients by using a separate
triangle cover whose three branch indices are the full orders. For
proper solutions, prime valuations outside the coefficient supports
give specializations unramified off a fixed set. Their argument uses
finiteness of the resulting number fields and Faltings' theorem on
the associated high-genus curve. The weighted projective map above
has indices \(h_a,h_b,h_c\). The minimal-repair problem requires a
uniform count as the coefficient supports vary with the height.
[Darmon–Granville, Theorem 2 and §3](https://www.math.mcgill.ca/darmon/pub/Articles/Research/12.Granville/pub12.pdf)
· [Minimal target-repair equivalence](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)
