# Numerical width and valuation quantiles on a complete fibre

Let \(n=X^L\), with \(L\) fixed before taking the height limit. Work
with positive primitive \(a+b=c\), \(X/2<c\le X\), and a fixed bound
\(C\) on each normalized radical cap. Fix one complete letter \(T\).
Let \(U,V\) be the moving letters. They are related by
\(V=\epsilon U+B\), where \(\epsilon\in\{-1,1\}\) and \(B\)
is fixed. The position of \(T\) in the additive equation is arbitrary.

Consider a height-first limiting fibre of count exponent \(d>0\).
Let \(h_U,h_V\) be its numerical heights in base \(X\), and put
\(u=\min(h_U,h_V)\). Use exact valuation layers, with one
factor per prime \(p\) in the numerical weights and normalized
mass \(\log p/\log n\) at scaled order \(L v_p\).
Write \(\nu_W\) for the weak limit of the resulting mass measures
in anatomic refinements, \(l_W=\nu_W([0,\infty))\), and
\[
e_W=h_W-\int r\,d\nu_W(r).
\]
All such refinements fix the radical and height statistics up to
arbitrarily small losses in the inner exponent. They can equivalently
be taken as the truncation and box selection of the three-variable
anatomic reduction, restricted to the selected fibre.
[Anatomic conventions and escaped moments](ABC_RESIDUAL_QUANTILE_SPECTRUM.md#the-residual-quantile)
· [Exact valuations and retained subfamilies](ABC_QUANTITATIVE_DEGREE_WINDOWS.md)

## Full-divisor locks leave a uniform amount of numerical width

For each \(s>0\) and \(C\) there is \(\theta(s,C)>0\) with the
following property. Suppose the fibre is restricted to a single
class of prescribed full divisors \(F_U\mid U,F_V\mid V\) with
a limiting count exponent at least \(s\). If \(w\) is a lower
limit of \(\log_X(F_UF_V)\) in that class, then
\[
                         w+\theta(s,C)\le u.                 \tag{FW1}
\]
The divisors can vary with the fibre and height. The estimate
uses their numerical sizes and holds without a restriction on
their supports.

Indeed, \(F_U,F_V\) are coprime by primitivity. The conditions
\(U=0\bmod F_U\) and \(\epsilon U+B=0\bmod F_V\)
fix a single class for \(U\) modulo \(F_UF_V\). The admissible
\(U\) lie in an interval of length at most
\(X^{\min(h_U,h_V)+o(1)}\), since \(U,V\) have the same
absolute differences. Translate the progression index so its
least occurring index is zero. If \(u-w\) is less than
half the rank-one small-feature threshold for exponent \(s\),
then all indices lie in \([0,X^{\theta_{\rm FT}}]\) on
sufficiently late boxes and heights. The rank-one estimate
gives an exponent strictly below \(s\), a contradiction.
We can take \(\theta(s,C)\) to be half that threshold.
[Rank-one small features on a fixed complete fibre](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#a-small-feature-bound-on-a-single-fibre)

In particular, with \(d\ge d_0>0\),
\[
                      e_U+e_V+\theta(d_0/2,C)\le u.           \tag{FW2}
\]
Fix a large finite continuity cutoff \(R\). The product of
full prime powers at \(Lv_p>R\) on either moving letter has
radical at most \(n^{(1+o(1))/R}\), since the sum of its
numerical weights is at most \(1+o(1)\). Enumerating the two
tails loses at most \(2/R+o_X(1)\) in the normalized count
exponent by radical enumeration. Their combined numerical
weight tends to
\[
 h_U+h_V-\int_{[0,R]} r\,d(\nu_U+\nu_V)(r).
\]
Apply (FW1) with \(s=d_0/2\) and then let \(R\) tend to infinity.
The cutoffs used in each application are fixed before the inner
height limit. The gap is independent of the cutoffs.
[Radical and supported-part enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

## Quantile packing with atomic and nonatomic profiles

Let \(\sigma=\nu_U+\nu_V\), \(M=l_U+l_V\), and let
\(F_\sigma(t)\) be the maximum first moment of mass \(t\)
taken from the highest scaled degrees, allowing fractional
mass at the last degree. Put
\[
 g=\max_{W=U,V}\sup_{x\ge0}\nu_W(\{x\}).
\]
Then for \(0<t<M\), \(t+g+s<d\) and \(s>0\),
\[
            e_U+e_V+F_\sigma(t)+\theta(s,C)\le u.              \tag{FW3}
\]
The overshoot \(g\) may instead be bounded by the largest
individual prime mass at the chosen prelimit cutoff.
For the nonatomic parts \(\sigma_{\rm na}=(\nu_U)_{\rm na}
+(\nu_V)_{\rm na}\), the analogous assertion is
\[
           e_U+e_V+F_{\sigma_{\rm na}}(t)+\theta(s,C)\le u
              \qquad (0<t<\sigma_{\rm na}([0,\infty)),
                         \ t+s<d).                            \tag{FW4}
\]
Either formula extends in \(t\) by continuity at the endpoints
where the stated cost condition permits approach from below.

For (FW3), choose the most expensive degrees in the prelimit
pooled list, with fractional total mass \(t\), and round up to
the last entire prime. By the usual quantile identity the
fractional weight tends to \(F_\sigma(t)+e_U+e_V\).
The maximal single-prime cost has limit superior at most \(g\).
A prime of cost bounded below has bounded scaled degree by
the first-moment bound, so on a subsequence it produces an
atom on its own letter of at least that cost. The two
squarefree supports of the selected divisors have at most
\(n^{t+g+o(1)}X^{o_X(1)}\) possibilities, by enumeration
of their radical followed by supported powers. Choices of a
selection rule for individual rows do not require a separate
count. The radical of the selected product, its two support
allocations and the powers suffice to record the divisors.
For large inner heights one class retains exponent at least
\(s\), so (FW1) applies. This uses no assumption that the
divisors determine the row by themselves.
[Quantile identity with escape and single-coordinate rounding](ABC_RESIDUAL_QUANTILE_SPECTRUM.md#one-coordinate-rounding-and-the-resulting-inequalities)

For (FW4) fix a small positive loss allowance \(\eta\).
Use a finite upper cutoff \(R\) large enough that the
nonatomic moment above \(R\) is below \(\eta\) and the two
full tails cost below \(\eta\). Fix these tails as additional
records. Inside \([0,R]\), take finitely many atomic locations
whose omitted atomic mass is below \(\eta/R\). Remove short
neighborhoods of these locations, with boundaries at
continuity points, losing nonatomic mass and moment less
than \(\eta\). The retained part is a finite union of compact
continuity intervals. If atoms of mass at least \(\eta\)
remain there, remove their neighborhoods as well. No retained
individual prime has cost greater than \(2\eta\) in late
boxes, by the single-prime atom argument. A prelimit
quantile of mass at most \(t\) on this set, rounded up,
has weight at least \(F_{\sigma_{\rm na}}(t)-O(\eta)\)
in the limit. To see this, restrict an optimizing nonatomic
selection to the retained set, losing only the removed
nonatomic moment. The optimal quantile on the retained
full measure is at least as large as this restriction.
Its prelimit weights converge by the quantile identity on
the compact set. If the available mass is less than \(t\),
take the entire set. The cutoff tail weights contribute
at least \(e_U+e_V\) in addition. The radical cost of the
combined records is at most \(t+O(\eta)\), leaving exponent
\(s\) once \(\eta\) is sufficiently small. Apply (FW1),
then let \(\eta\) decrease to zero with \(s\) fixed.

## Transfer from counts with prescribed complete-letter multiplicity

Suppose a selected balanced family has exponent \(D\) and
the number of different complete \(T\)-values in it has
exponent at most \(D-r\), where \(r>0\). The profiles can
first be made common by anatomic box selection and the
complete \(T\) can then be fixed, leaving a fibre of exponent
at least \(r\). Thus (FW2)–(FW4) hold on this family with
\(d\) replaced by \(r\). The same restrictions hold if one
directly prescribes a lower bound \(n^{r-o(1)}\) on complete
\(T\)-multiplicity after valuation box selection.

These statements give uniform exponent losses on sectors
with margins. For example, fix \(r_0>0\) and \(t,s>0\)
with \(t+s<r_0\). A balanced family with total exponent
approaching \(1/3\), complete \(T\)-dictionary exponent at
most \(1/3-r_0\), and valuation data forcing
\(e_U+e_V+F_{\sigma_{\rm na}}(t)>u-\theta(s,C)\)
by a fixed positive margin cannot saturate. A finite-box
version is obtained by prescribing numerical and radical
bounds on the tail and finite packets before enumerating
them. The finite-box margin and rank-one bound give constants
at each fixed small budget, with the permitted
\(X^\varepsilon\) height losses. This criterion does not
need minimum dictionaries or an upper bound on the number of
atoms. It applies to any positive complete-target multiplicity
stratum in the least-repaired primitive quotient problem:
after choosing the original height shell, fixing a target
support costs at most a height-subpower factor for its possible
complete values. The least repair is prescribed by the
quotient, so that no common multiplier is counted.
[Least repair and fixed target-support multiplicity](ABC_FIXED_LETTER_RADICAL_CAPS.md#size-of-a-minimally-repaired-target-radical-fibre)

The test is vacuous when the available complete-target
exponent is zero. Even at a positive exponent it need not
exclude finite dispersed spectra. For instance, the formal
test data \(u=1\), \(e_U=e_V=0\), and each moving measure
uniform of mass \(1/3\) on a short symmetric interval
about order \(3\) have quantile weight near \(3t\).
For \(t\le1/6\) this only uses about half of the available
unit width. This observation concerns the inequality on
formal measures and does not construct saturated arithmetic
families.


## Discrete rounding on the full three-letter family

The same prime-level selection improves the rounding alternative
in the three-letter quantile test. Fix a dyadic shell on each
complete letter, with upper height exponents \(h_A,h_B,h_C\),
and let \(H_*=\min_{i\ne j}(h_i+h_j)\). On each row pool the
individual primes of all three letters. The cost of a prime
\(p\) is \(\log_n p\) and its weight is \(\log_X(p^{v_p})\).
Sort in descending order of \(Lv_p\), breaking ties by a fixed
rule. Let \(Q(t)\) be the fractional weight obtained by paying
total cost \(t\), and let \(g(t)\) be the cost of the
single prime at which the last fractional purchase is made.
If the threshold needs no fractional purchase put \(g(t)=0\).
Restrict to rows with sufficient total cost for \(Q(t)\)
to be defined.

For fixed \(t,\eta,\gamma>0\), the number of such rows with
\(Q(t)>H_*+\gamma\) and \(g(t)\le d_*-t-\eta\) is at most
\[
                         n^{d_*-\eta+o_X(1)} .                 \tag{FW5}
\]
Indeed, round just that prime up and enumerate the three
selected divisors at cost at most \(d_*-\eta\). For two rows
sharing a record, the determinant formed by their coordinates
on any two complete roles is divisible by the product of the
three selected divisors. The divisibility from the third
role follows from the signed additive equations. On the
pair attaining \(H_*\), the determinant is at most
\(2X^{H_*}\) in absolute value, whereas the recorded product
has weight above \(H_*+\gamma\). The determinant vanishes.
The additive equation and primitivity recover the row.
[Three-factor projective quantile bound](ABC_RESIDUAL_QUANTILE_SPECTRUM.md#one-coordinate-rounding-and-the-resulting-inequalities)

At a balanced one-third endpoint this gives a discrete
alternative without any assumption on how many limiting
atoms there are. If a pooled top cost \(t<1/3\) violates its
numerical weight threshold by a fixed margin, any surviving
full-exponent cohort requires a single actual pivot prime
of cost at least \(1/3-t-o(1)\) in each row. A weak limiting
atom split among small primes cannot supply that pivot.
The prime then satisfies \(v_p\le1/(L\log_n p)\), so its
exact valuation has only finitely many possibilities at
each fixed budget. The required prime may vary across
rows. The statement supplies a discrete alternative in
every balanced valuation profile. It leaves profiles
with smaller top-quantile weight and profiles with an
actual macroscopic pivot prime for further estimates.


The pivot also has a lower bound on its numerical weight.
Fix \(0<\eta<d_*-t\). If \(f\) is the fractional
cost spent on the pivot and \(m=\log_n p\) its full cost,
discard the rows for which the round-up cost
\(t+(m-f)\) is at most \(d_*-\eta\). The same lock count
as in (FW5) bounds this larger class. Discard as well
the rows for which rounding down already gives weight
exceeding \(H_*+\gamma/2\). Their lock costs at most \(t\).
Writing \(s=Lv_p\), every surviving row satisfies
\[
       sf\ge\gamma/2,\qquad m-f>d_*-t-\eta,\qquad
       sm\ge s(d_*-t-\eta)+\gamma/2.                            \tag{FW6}
\]
These bounds hold on the same actual prime. In particular
if the pivot is on role \(W\), its full weight is at most
\(h_W\), which also bounds the possible order \(s\).
This couples the discrete mass obstruction to the numerical
deficit of the fractional packing.

