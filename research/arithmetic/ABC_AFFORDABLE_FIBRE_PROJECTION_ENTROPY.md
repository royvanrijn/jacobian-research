# Projection entropy and target templates in divisor fibres

A fixed pair of full divisors constrains collisions of sums through
the arithmetic lattice of the original letters. This identifies
the entropies of the projected weighted sum distributions up to
the square of the ordinary affine-lock scale. The same lattice can
recover a whole tuple from a normalized template on the target alone,
without a small bound on the relative additive indices.

## Weighted sums and projection entropy

Fix coprime positive integers \(F,G\), signs \(\sigma,\tau\), and a
finite nonempty family \(\mathcal R\) of integer triples
\[
                C=\sigma A+\tau B,\qquad F\mid A,\quad G\mid B .
\]
Let \(W_A,W_B\) be the widths of the two free complete values on the
family, and put \(W=\min(W_A,W_B)\).
Choose \(k\) rows independently from any probability law on \(\mathcal R\).
Denote their integer vector sum by \(S\), and its target projection
by \(S_C\). When \(W=0\), the target determines the vector sum.
When \(W>0\), put
\[
    p_*=\min\left(1,\;2\exp\!\left(-\frac{(FG)^2}{8kW^2}\right)\right).
\]
Then for some \(0\le p\le p_*\), with \(h_2\) the binary Shannon
entropy in natural logarithms, one has
\[
 0\le H(S)-H(S_C)
   \le h_2(p)+p\log\!\big((1+k W_A)(1+k W_B)\big).
 \tag{PE1}
\]
For an independent copy \(S'\), the excess of the scalar collision
probability \(\mathbb P(S_C=S'_C)\) over the vector collision probability
\(\mathbb P(S=S')\) is at most \(2p_*\).

Indeed, two vector sums with the same target sum have differences on
the free letters equal up to sign, and each difference is divisible by
\(FG\). Choose the free role of width \(W\). On the event that its sum
is within \(FG/4\) of its expectation, any two such values differ by
at most \(FG/2\). The target therefore determines the vector sum on
this event. Hoeffding's inequality gives probability at most \(p_*\)
for its complement. Condition on the binary indicator of that complement.
The conditional entropy on the central event is zero, and on the
complement it is at most the logarithm of the total vector support.
The support bound follows by choosing both free integer coordinates
in intervals of widths \(kW_A,kW_B\). This proves (PE1).
The collision assertion follows by excluding both complementary events.
[Hoeffding, Theorem 2](https://www.cs.rpi.edu/academics/courses/spring06/random/hoefding.pdf)

Without weighting, the same congruence makes the target projection
injective on the entire \(k\)-fold vector sumset when \(kW<FG\).
For weighted sums it suffices for an error smaller than any negative
power of \(X\) that \(FG/W\ge X^\alpha\) and
\(k\le X^{2\alpha-\eta}\), with fixed \(0<\eta<2\alpha\)
and all original values at most \(X\).
A fixed divisor on the target permits the same argument on either
of the other roles by exchanging the roles. If the corresponding
lock margins hold for two projections, their entropies and weighted
collision probabilities agree to the same error. This concerns
the distribution of sums, whose full unweighted support can still
have collisions outside the central region.

The power two in this range is sharp as a geometric test. Choose
distinct primes \(F,G\asymp H\), both greater than five, and take
the four positive primitive triples formed by
\(A\in\{F,3F\}\), \(B\in\{G,5G\}\), \(C=A+B\).
Now \(FG/W\asymp H\). Under the uniform law the two free sums encode
independent \(\operatorname{Bin}(k,1/2)\) variables, so
\(H(S)\ge\log k-O(1)\) by the largest binomial atom bound.
The target sum has variance \(O(kH^2)\). Adding an independent
uniform variable on \([0,1)\) preserves its discrete entropy as
differential entropy, and the Gaussian entropy bound gives
\(H(S_C)\le\tfrac12\log k+\log H+O(1)\).
For \(k\ge H^{2+\eta}\) the entropy difference is therefore at
least \(\tfrac12\eta\log H-O(1)\).
The example tests the lattice statement on actual primitive triples
without a small-radical constraint.

## Counting with a normalized target template

Continue with the same signed equation, now on positive primitive
triples with caps \(A\le Y_A,B\le Y_B,C\le Y_C\), and also fix \(D\mid C\).
Write \(S_0=\min(Y_A,Y_B)\) and assume \(FG>S_0\).
In particular each complete target determines its original row.
Fix a normalized target tuple \((c_i^*)_{i=1}^m\) of positive
integers with greatest common divisor one. Put \(M_C=\max_i c_i^*\).
Seek ordered tuples in the fibre with
\[
  C_i=g_C c_i^*,\qquad g_C=\gcd_i C_i,\qquad
  g_C/D\le H_0.
\]
The number of such tuples is
\[
 O\!\left(1+\frac{S_0}{FG}
                \min\!\left\{H_0,\frac{Y_C}{D M_C}\right\}\right).
 \tag{PE2}
\]
If the product inside the second term is strictly less than one,
there is at most one. The assertion uses no restriction on the
additive relative indices of the target tuple.

For the proof take \(Y_A=S_0\), exchanging roles if necessary.
Since \(D\) divides every \(C_i\), write \(g_C=D h\).
The first original row gives a positive coprime pair
\((A_1/F,h)\), because \(g_C\) divides \(C_1\) and
\(\gcd(A_1,C_1)=1\). Its bounds are \(A_1/F\le S_0/F\) and
\(h\le H_1=\min\{H_0,Y_C/(D M_C)\}\).
Reducing that row modulo \(G\) places \((A_1/F)/h\) in a single
unit residue class modulo \(G\). Indeed \(F,h,D,c_1^*\) are
units there by primitivity. Distinct such reduced pairs have
nonzero determinant divisible by \(G\). After scaling the two
coordinates by their upper bounds their angular spacing is at
least \(G/(2(S_0/F)H_1)\), giving (PE2). Alternatively an absolute
determinant bound by \((S_0/F)H_1<G\) proves the uniqueness
assertion. Each pair gives the content \(g_C\), hence all the
complete targets. Uniqueness at each target gives the whole tuple.
This uses projective packing on a single role's content, in
contrast to packing the ratio of two contents after reconstructing
two templates.
[Two-template count (AFP7)](ABC_AFFORDABLE_FIBRE_ADDITIVE_PATTERNS.md#two-patterned-letters-and-a-full-divisor-on-the-third)

### A target-only defect count

Fix a relative additive shape \(\ell_1,\ldots,\ell_m\) as in (AFP2),
and charge only the target defect in (AFP3). The one-template
argument of (AFP4) counts the normalized target templates at cost
\(O(\mathcal B X^\varepsilon)\) when the outside defect radical
has cap \(\mathcal B\), for the even square case or the three-term
prime-order cases in that result. Apply (PE2) to each template.
When \(S_0H_0<FG\), the resulting tuple count is
\(O(\mathcal B X^\varepsilon)\). Neither the summand defect records
nor the hypothesis \(2\max_i|\ell_i|\max_i v_i<F\)
from (AFP1) is needed. The improvement allows a dictionary of
long relative shapes, whose per-shape defect bounds can be summed.
[Target-template reconstruction (AFP4)](ABC_AFFORDABLE_FIBRE_ADDITIVE_PATTERNS.md#a-count-by-the-joint-defect-record)
· [Powerful-progression results of Bajpai, Bennett and Chan](https://arxiv.org/abs/2302.03113)

For example, start with an affordable record of numerical margin
\(FG/S_0\ge X^{\eta_0}\). Apply the original-letter support peeling
of (DC3) with tolerance \(0<\eta<\eta_0\). The procedure can enlarge
the three full factors by fixing actual supported parts, and the
margin for the enlarged \(F,G\) persists. At a given anchor take the
pool of its residual target primes. For all but a power-small
fraction of queries, the gcd of the two residual targets divides
a pool-supported part of weight less than \(\eta\). The same then
holds for the common residual content of any fixed-length tuple
containing the pair. Thus in these tuples \(H_0=X^\eta\) meets the
uniqueness criterion.
On any collection of such tuples covering the targets of interest,
a shape dictionary exponent \(s\) and target defect cap exponent \(f\)
give a count of covered rows at most \(n^{s+f+o_X(1)}\) in the fixed
record. A lock record dictionary of exponent \(t\) can be restored
at total exponent \(t+s+f\) with height-subpower losses.
[Original-letter support peeling (DC3)](ABC_AFFORDABLE_FIBRE_DETERMINANT_COMPLEXITY.md#removing-concentrated-support-from-a-large-fibre)

The entropy transfer supplies a common test for the projections of a
fibre. The target-template bound supplies a small count when the actual
shape and defect records compress. An unrestricted use still requires
an estimate forcing such compression or separately counting the
families with large projection entropy and expensive shape records.
