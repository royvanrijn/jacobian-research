# Compressed shifts and sparse source codegrees

A pair of representations of one complete target gives two
primitive difference triples. When restoring the multiplier of
one of them, only primes absent from all three quotient entries
need new support records. The three corresponding deductions
from the quotient's radical budget pay for that record exactly.
This yields a bound for complete source pairs before the other
source or the shared target is counted. It gives a multiplicity
saving when those pairs have sufficiently small actual target
codegrees.

## Counting the one-side records

Let \(n=X^L\), with \(L>0\) fixed before taking the height limit.
Consider any family of positive primitive rows \(A+B=C\) in
\(X/2<C\le X\), with caps
\[
        \operatorname{rad}(A)\le n^{r_A},\qquad
        \operatorname{rad}(B)\le n^{r_B},\qquad
        \operatorname{rad}(C)\le n^{r_C}.
\]
The nonnegative bounds \(r_A,r_B,r_C\) are fixed. A pair of
distinct rows on the same complete target has shift \(h>0\).
On either moving role \(W\in\{A,B\}\), its smaller and larger
values are \(U,U+h\). Put
\[
 g=(U,h),\qquad (u,H,v)=(U/g,h/g,(U+h)/g).
\]
Call this side tall if \(v\ge \sqrt X/2\). With \(g_A,g_B\)
defined on the two sides, \((g_A,g_B)=1\) by original
primitivity, so \(g_Ag_B\mid h\). Both larger original
values are at least \(h\), and at least one exceeds \(X/4\).
Thus the product of the two normalized larger values
is at least \(X/4\), and at least one side is tall.

Fix a prescribed pool \(P\le X^B\), with \(B>0\) fixed.
For an integer \(s>0\) write
\(\operatorname{rad}_P(s)=\prod_{p\mid s,\ p\nmid P}p\).
The compressed shift record on the chosen side has cost
\[
 \log_n\operatorname{rad}_P(H)
   +\log_n\!\prod_{\substack{p\mid g\\p\nmid PuvH}}p .
 \tag{SCM1}
\]
The second product is the support of the multiplier that is
missing from all three primitive entries and from the pool.
The two supports in (SCM1) are disjoint. Fix a raw bound
\(\operatorname{rad}(h)\le n^{D_0}\). It is automatic with a
sufficiently large fixed \(D_0\) if the compressed cost is
bounded and \(\operatorname{rad}(P)\le n^{B_1}\) for fixed
\(B_1\). Indeed a prime of \(h\) outside \(P\) is either
in \(H\), missing from the three primitive entries, or in
one of the primitive endpoints.

**SCM2.** Fix \(B,D_0,\eta>0\) and the three raw caps.
At all sufficiently small fixed \(L>0\), uniformly in
\(0\le d\le D_0\) and in \(P\), the number of distinct
ordered pairs of source values \((U,U+h)\) on a given tall
side \(W\) with compressed cost at most \(d\) is
\[
       O\!\left(n^{(2r_W+d)/3+\eta}X^\varepsilon\right).
 \tag{SCM2}
\]
There is no charge for the dictionary of complete targets
with these source pairs. The alternative exponent
\(r_W/2+d+\eta\) also holds.

To prove it let \(\alpha\log n\) be the logarithm of the
second product in (SCM1). Every prime in that product was
on both original endpoints and on \(h\). Hence each
primitive endpoint has raw radical at most
\(n^{r_W-\alpha}\), and the outside-\(P\) radical of \(H\)
is at most \(n^{d-\alpha}\). The quotient's radical outside
\(P\) is therefore at most \(n^{2r_W+d-3\alpha}\).
Its raw cap is bounded uniformly in base \(n\), since
the original endpoints and \(h\) have bounded raw caps.

Bin \(\alpha\) on a fixed mesh of width \(\delta>0\). In
a bin with lower endpoint \(\alpha_0\), apply (EB4) to
the tall primitive quotient with paid cap
\(n^{2r_W+d-3\alpha_0}\), using ordinary prime prices
off \(P\) and zero price on \(P\). Split into dyadic shells
of height \(Y\ge\sqrt X/2\). One fixed small raw budget
\(C_*L\), with \(C_*\) sufficiently large in terms of the
raw bounds, is valid for all these shells. The price
cap at base \(Y\) is
\((2r_W+d-3\alpha_0)L\log X/\log Y\), giving the
counting exponent \((2r_W+d-3\alpha_0)/3\) in base \(n\)
before the allowances. The pool has bounded polynomial
height in \(Y\), and its prime count is
\(O(\log Y/\log\log Y)\), so the price condition is
uniform. Once the quotient is recorded, choose the
missing multiplier support at cost \(n^{\alpha_0+\delta}\).
All multiplier powers are now on a prescribed support,
namely the quotient support, \(P\), and the newly
recorded support. Fixed-support enumeration restores
the multiplier and hence \(U,U+h,h\) at subpower cost.
The loss \(\alpha_0\) cancels exactly. Take \(\delta\)
and the normalized allowance in (EB4) small in terms
of \(\eta\). This proves the uniform estimate with
the height limit first.
[Weighted limiting one-third bound (EB4)](ABC_EFFECTIVE_RADICAL_BUDGET.md#counting-by-a-sublogarithmic-prime-partition-sum)
· [Uniform fixed-support enumeration (EB3)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

For the alternative, enumerate the joint support in (SCM1)
at cost \(n^d\), and divide it between its two disjoint parts.
The allocation costs only a subpower factor. Then record
the exact value of \(H\) on its specified support and the
pool, again at subpower cost. Apply (FL1) on the tall
primitive difference shell with fixed complete \(H\) and
free radical caps \(r_W\) in base \(n\). Its raw fixed
letter cap is bounded by the raw shift cap. This counts
\(O(n^{r_W/2+\eta'}X^{\varepsilon'})\) normalized pairs.
The missing multiplier support has already been recorded,
so restore its powers at subpower cost as before.
This gives the alternative exponent.

## Actual codegrees and a strict multiplicity tail

Put \(n^v\) as an upper bound for the number of complete
targets which use both members of any source pair on
the indicated side \(W\) in the selected family. This is
an actual codegree in the projection on \(W,C\), and
can be replaced by a bound on the degree of either
member of the pair. A target completes a fixed source
pair in at most one way because it fixes the other two
source values.

**SCM3.** In a cohort of targets for which each has at
least \(n^{2u-\rho}\) pairs of distinct rows with a tall
\(W\)-side and compressed shift cap \(d\), the number of
complete targets is at most
\[
       O\!\left(n^{(2r_W+d)/3+v-2u+\rho+\eta}
                          X^\varepsilon\right).
 \tag{SCM3}
\]
The raw shift cap and pools are those of (SCM2). If the
tall side is not prescribed and each target has that
many pairs satisfying the condition on one side or the
other, subdivide the targets according to which side
supplies at least half its pairs, then apply (SCM3).
One may use different source caps, compressed caps
and codegree caps in the two estimates. One may always
replace \((2r_W+d)/3\) by
\(\min\{(2r_W+d)/3,r_W/2+d\}\) from (SCM2).

In the balanced normalization \(r_W=2c_0+o(1)\),
\(c_0=1/6\), the desired comparison at multiplicity
exponent \(u\) is against \(n^{2c_0-u}\). The strict saving
criterion on this cohort is
\[
                  d+3v+3\rho<2c_0+3u
 \tag{SCM4}
\]
with a fixed margin left for the cap and count tolerances.
This includes the sparse case \(v=o(1)\), \(d\le 2c_0\)
and \(\rho=o(1)\) for every fixed \(u>0\).
In particular, the codegree bound extends the range
of \(d\) beyond \(2u\) if
\(3v+3\rho<2c_0+u\). The comparison with \(2u\)
is the fixed-endpoint bound applied to the compressed
shift. Indeed at a fixed original endpoint one can use
the pool \(P\) together with the primes of that endpoint
on \(W\): all primes of \(g\) then lie in the pool,
so the outside radical of \(h\) is bounded by that of
\(H\). The fixed-endpoint estimate on near-full balanced
input gives the \(2u\) threshold without a codegree
hypothesis. The additional gain here comes from
counting the source pair once and then charging only its
actual target completions.
[Fixed-endpoint shifts, including adaptive pools (FE2)](ABC_FIXED_ENDPOINT_SHIFT_RADICAL_COUNTS.md#a-uniform-bound-inside-any-fixed-complete-fibre)

The cap parameters may vary by target. For instance,
if the chosen \(d(C),u(C),v(C),\rho(C)\) stay in fixed
bounded ranges and (SCM4) has a uniform margin, take
a sufficiently fine fixed mesh of the four exponents.
Use the upper ends for \(d,v,\rho\) and the lower end
for \(u\) in (SCM3). The associated weighted target
count also gives a strict saving when the original
multiplicity on each bin is at most \(n^{u+o(1)}\).
This uses a fixed mesh selected after the margin and
before the height limit.

There is a useful converse on source projections. If
a critical cohort has \(n^{2c_0-u-o(1)}\) targets of
multiplicity \(n^{u+o(1)}\), each carrying
\(n^{2u-o(1)}\) qualifying pairs with a prescribed
tall orientation and \(d\le2c_0+o(1)\), then its
qualifying pair count is \(n^{2c_0+u-o(1)}\).
(SCM2) gives at most \(n^{2c_0+o(1)}\) source pairs.
Those touching a source value of degree below
\(n^{u-\delta}\) can account for only
\(n^{2c_0+u-\delta+o(1)}\) of the pairs, for any
fixed \(\delta>0\). The high-degree source dictionary
has size at most \(n^{2c_0-u+\delta+o(1)}\) when the
original cohort has at most \(n^{2c_0+o(1)}\) rows.
The remaining high-degree rows support the
full-exponent pair cohort. This extracts a smaller
complete dictionary if the qualified multiplicity
occurs without the sparse-codegree saving.

## A uniform fixed-shift bound without dictionary hypotheses

**SCM5.** With the same fixed raw caps and uniform
bounds as above, for any one fixed \(h\) with raw
radical at most \(n^{D_0}\), the number of shared-target
pairs of rows with this shift is
\[
    O\!\left(n^{(r_A+r_B)/2+\eta}X^\varepsilon\right).
 \tag{SCM5}
\]
In addition to its pointwise use, this may replace \(v\)
by \(r_{\bar W}/2\) in the upper bound (SCM3), using
the completion bound in the proof below. Neither
upper bound presupposes a small dictionary of source
or target values.

Choose a tall side, then choose its \(g\mid h\) at
divisor cost. Its normalized pair has a fixed
complete shift \(h/g\). Apply (FL1) on dyadic
shells above \(\sqrt X/2\), with the two free caps
both at most \(r_W\) in base \(n\). This gives at
most \(n^{r_W/2+\eta'}X^{\varepsilon'}\) choices
of the primary source pair. For this application
use the same fixed \(L\) on a shell with upper
endpoint \(Y\), and replace a cap \(b\) by
\(b\log X/\log Y\) in base \(Y^L\). The shell
height lower bound keeps all caps uniform. The
fixed shift also has bounded radical in these units.

It remains to count completions of any fixed primary
record, with \(h\) now known and the original primary
source endpoints fixed. Choose the other multiplier
\(g'\mid h\), again at divisor cost. Fix a small
\(\theta>0\) as in the small-cofactor consequence
of (NB1), with \(r=r_{\bar W}\). If that cap is
zero, use an arbitrarily small positive allowance
instead. If the larger normalized opposite endpoint
is at least \(X^\theta\), apply (FL1) at the fixed
complete shift \(h/g'\) in its normalized difference
triple. The caps convert to base \(Y^L\) with
\(\log X/\log Y\le1/\theta\). The count is
\(O(n^{r_{\bar W}/2+\eta'}X^{\varepsilon'})\)
after summing the dyadic shells. If the normalized
pair is smaller, take either of its original rows
and fix its primary source as the complete letter.
The opposite original source is divisible by
\(Q=g'\le X\) and its cofactor is at most
\(X^\theta\). Now (NB1) counts at most
\(O(n^{r_{\bar W}/2+\eta'})\) of those
original rows. Both original radical caps and the
third cap hold in that fibre. The opposite complete
value together with the primary record and the
sign of the shift recovers the entire pair. This
proves the completion bound and (SCM5).
[Uniform fixed-letter counts (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)
· [Small cofactor on an original fibre (NB1)](ABC_SHARED_BOX_PRIVATE_VALUATION_BANDS.md#a-positive-valuation-band-on-a-fixed-letter-fibre)

For each fixed \(\chi>0\), the short branch alone has
completion count \(O(n^\chi X^\varepsilon)\) when its
normalized opposite endpoint is at most \(X^{\theta_\chi}\)
for a sufficiently small fixed \(\theta_\chi>0\). Apply
(NB1) with \(r=2\chi\) at the fixed original source, and
sum over the divisor-many choices of \(g'\mid h\). One
can use \(v=\chi\) in (SCM3) for pairs restricted to this
short branch, even when the unrestricted source codegree
is large.

For comparison, using a direct upper bound on the raw
shift radical \(n^b\) in the symmetric balanced caps
gives two aggregate exponents from (SCM2) and (SCM5):
\(2c_0+b\) by direct enumeration of \(h\), and
\((7c_0+b)/3\) by the multiplier cancellation.
Their minimum does not improve the aggregate raw-shift
bound \(2c_0+\min\{b/2,c_0\}\) obtained from (FE2)
on near-full inputs. The separate pair and codegree
bounds in (SCM3), and the pointwise conclusion in
(SCM5), do not require that aggregate comparison.
[Earlier rich-shift counts with minimum dictionaries](ABC_POPULAR_SHIFT_RADICAL_AMPLIFICATION.md#a-bound-for-rich-shifts-with-a-radical-cap)

## Thinning and the repaired normalization

On one fixed complete target the retained
representations form a matching on the two original
summand parts. Under independent Bernoulli thinning
of the two parts with probabilities of product \(p\),
its representations therefore survive as independent
Bernoulli variables of parameter \(p\). Each of the
qualifying cherries survives with probability \(p^2\).
Thus (SCM3) bounds the expected total surviving
qualifying cherries after multiplying the corresponding
unthinned cherry bound by \(p^2\).
If there are \(m_C\) original representations and
\(e_C\) qualifying unordered pairs, the total
number of their surviving pairs, of mean
\(\mu=p^2e_C\), has second moment at most
\(\mu^2+\mu+2m_C p\mu\). Indeed disjoint pairs have
independent survival and every pair meets at most
\(2m_C\) other pairs. Hence the probability of at least
one such pair is at least
\(\mu/(1+\mu+2m_Cp)\). At \(p=n^{-u}\), a
constant-density pair relation on a target with
multiplicity of order \(n^u\) survives with probability
bounded below. The arithmetic estimate itself
does not require random thinning.

In the minimally repaired target setting let
\(Z=Y^{1/M}\) and assume
\(\operatorname{rad}(C)\ge C^{a_-\lambda}\) in
the balance sector. The least repair bound implies
\(\operatorname{rad}(C)\le Z\). Thus an original
dyadic shell satisfies \(X^L\ll Z^{1/a_-}\) when
\(L=\lambda\). A prescribed radical has only
height-subpower many complete targets by
fixed-support enumeration, also subpower in \(Z\)
at fixed \(\lambda,a_-\). One may therefore use a
maximizing complete representative for each occurring
radical, or replace a sum over those representatives
by its maximum at this cost. For example, if such
targets each supply at least \(Z^{2u_*-\rho_*}\)
qualifying pairs and the normalized caps \(d,v\)
are uniform in each original shell, (SCM3) gives
their support count at exponent at most
\((2r_W+d+3v)/(3a_-)-2u_*+\rho_*\)
in \(Z\), up to the fixed exponent allowances
and dyadic-shell summation. A strict margin below
\(1-u_*\) gives a strict support tail saving on
this cohort. The bound counts original primitive
quotients and uses no arbitrary repair multiplier.
[Fixed target supports in the least-repair count (FL5, FL8)](ABC_FIXED_LETTER_RADICAL_CAPS.md#size-of-a-minimally-repaired-target-radical-fibre)

The savings in (SCM4) concern targets with many
pairs of low compressed shift cost and an appropriate
actual source-pair codegree. For the unrestricted
tail, one still needs to force such pairs or count
the targets whose shifts have larger compressed
cost.


## Endpoint support and original-target discounts

The multiplier restoration argument can charge the prime support
retained in the primitive difference. The original target's radical
also gives a codegree discount on some source pairs. These
discounts use different records.

For the normalized triple and \(g\) in (SCM1), let \(d_W\)
denote its actual compressed cost, and retain the notation
\(\alpha\) for its second summand. In addition set
\[
\begin{split}
 \beta&=\log_n\!\prod_{\substack{p\mid(g,H)\\p\nmid P}}p,\\
 \gamma&=\log_n\!\prod_{\substack{p\mid(g,uv)\\p\nmid P}}p,\qquad
 \delta=\log_n\!\prod_{\substack{p\mid g\\p\mid P}}p,\\
 D_W&=d_W-2\beta-\gamma-2\delta .
\end{split}                                                     \tag{SCM8}
\]
These split the support of \(g\) into four disjoint parts, since
\(u,H,v\) are pairwise coprime. The signed price \(D_W\) may
be negative. On the support counted by \(\beta\), both
original source values have lost the prime from their
normalized endpoints, while \(H\) retains it. On the
support counted by \(\gamma\), exactly one normalized
endpoint retains the prime. On the support counted by
\(\delta\), both original endpoints contain the prime and
both outside-\(P\) endpoints omit it.

**SCM9.** Under the hypotheses of (SCM2), the number of
distinct ordered source pairs on a given tall side with
\(D_W\le D\) is
\[
       O\bigl(n^{(2r_W+D)/3+\eta}X^\varepsilon\bigr).
                                                               \tag{SCM9}
\]
The bound holds uniformly over \(D\) in a fixed bounded
range. The same result allows any tall threshold \(X^\theta\)
with fixed \(\theta>0\); take the small-budget threshold
also in terms of \(\theta\).

Indeed, writing radical costs in base \(n\), the total paid
cost of \(u,v\) is at most
\(2r_W-2\alpha-2\beta-\gamma-2\delta\). The paid cost of \(H\)
is \(d_W-\alpha\). Hence the cap for applying (EB4) is
\(2r_W+D-3\alpha\). On each height shell and each fixed
\(\alpha\)-bin use the same restoration as in (SCM2).
The cost \(n^{\alpha+o(1)}\), with the prescribed mesh
allowance, cancels the deduction \(3\alpha\). No new
support need be recorded for \(\beta,\gamma,\delta\):
those primes are already in \(uvH\) or \(P\).
An \(\alpha\)-bin with negative upper paid cap is empty.
The change of tall threshold only changes the fixed bound
on \(\log X/\log Y\) in the quotient shells.

An exact cap convenient for comparisons is
\[
 \kappa=\log_n\operatorname{rad}_{PU}(h),\qquad
 b=\log_n\operatorname{rad}_{P}(h),\qquad
 d_W=\kappa+\alpha+\beta,\qquad \gamma=b-d_W .
\]
Indeed the primes of \(h\) in \(U\) are precisely those
of \(g\). The primes of \(h\) outside \(P\) omitted from
\(d_W\) are precisely the ones retained in a normalized
endpoint. Thus \(D_W=\kappa+\alpha-\beta-\gamma-2\delta\).
The discount for \(\gamma\) can improve the pair count
even when the raw shift cap is much larger than its
compressed cost.

**SCM10.** For a source pair \((U_1,U_2)=(U,U+h)\) on \(W\)
and a complete target \(C\) using it, write \(V_i=C-U_i\)
for its opposite source values. Let \(V\) denote this
opposite role. Define the actual support overlaps
\[
 z_C(C,h)=\log_n\operatorname{rad}\gcd(C,h),\qquad
 z_V(C,U_1,U_2)=
        \max_{i=1,2}\log_n\operatorname{rad}\gcd(V_i,hU_{3-i}).
\]
Restrict the qualifying pairs in (SCM3) further by
\(D_W\le D\), \(z_C(C,h)\ge z_C\) and
\(z_V(C,U_1,U_2)\ge z_V\), using fixed
\(0\le z_C\le r_C,\ 0\le z_V\le r_V\).
Let \(n^{v_0}\) bound the actual number of eligible
target completions of any source pair if such a bound
is supplied. Then the target count is at most
\[
 O\left(
 n^{(2r_W+D)/3+v_{\rm eff}-2u+\rho+\eta}X^\varepsilon
 \right),\qquad
 v_{\rm eff}=\min\left\{v_0,\frac{r_C-z_C}{2},
                             \frac{r_V-z_V}{2}\right\}.
                                                               \tag{SCM10}
\]
Omit \(v_0\) when no separate codegree is known. Choosing
smaller allowances for the two input bounds gives the
displayed allowance. A lower bound of \(\xi>0\) on the
desired margin lets the parameters vary in fixed bounded
ranges by the same finite-mesh argument as in (SCM3).
The least-repair support transfer uses
\((2r_W+D)/3+v_{\rm eff}\) as the numerator exponent
on each original shell.

To count actual completions, fix the entire source pair,
so \(h\) and both \(U_i\) are known. Split into the two
choices of \(i\) attaining the overlap requirement on \(V\).
Use (FL1) on the original row with complete fixed letter
\(U_i\). The pool on \(C\) is \(h\), and the pool on \(V_i\)
is \(hU_{3-i}\). Their heights are at most \(X^2\).
The remaining caps are \(r_C-z_C\) and \(r_V-z_V\).
This counts qualifying completions with exponent
\(v_{\rm eff}\) inclusive of the optional independent
codegree bound and the fixed allowances. The complete
target recovers both rows from their ordered source pair.
Sum over the distinct source pairs in (SCM9) and divide
by the qualifying multiplicity \(n^{2u-\rho}\).
There is no enumeration of the radical making up either
overlap. The prime pool is determined by the source pair.
Primes used in the \(h\)-part and in the \(U_{3-i}\)-part
of the \(V_i\) pool are disjoint on \(V_i\): a prime of
both \(h\) and \(U_{3-i}\) divides \(U_i\), contradicting
primitivity. Taking \(z_C=z_V=0\) already replaces the
unconditional completion exponent in (SCM5) by
\(\min(r_C,r_V)/2\) on a prescribed tall side. The short
opposite-cofactor argument of (SCM5) still gives an
arbitrarily small fixed completion exponent on that
more restrictive branch by using (NB1) with small \(r\).

In balanced caps \(r_A,r_B,r_C=2c_0+o(1)\), the strict
criterion is \(D+3v_{\rm eff}+3\rho<2c_0+3u\) with a fixed
margin. To compare with (FE2) on an already fixed
original row, use the pool \(PU_iV_iC\). Its cost on the
shift is exactly \(\kappa-z_{C,P}-z_{V,P}'\), where
\(z_{C,P}\) and \(z_{V,P}'\) count the supports of
\((C,h)\) and \((V_i,h)\) outside \(P\). These are
disjoint by primitivity. The threshold from this
application of (FE2) is that this cost be less than
\(2u\) by a fixed margin on the full pair exponent.
One can also apply (FL1) inside the fixed target after
excluding primes of the original summands, including
the shared \(g\)-support, from the corresponding moving
letters. The following comparisons concern the exponents
of these bounds, with arbitrarily small fixed tolerances.

* With no opposite or target overlap, take
  \(u=c_0/2\), \(\rho=0\), \(\alpha=\beta=\delta=0\),
  \(\gamma=3c_0/4\), and \(\kappa=9c_0/8\).
  The adaptive shift cost is \(9c_0/8>2u\).
  Removing the common source support from the original
  moving cap leaves half-cap \(5c_0/8>u\).
  The general codegree bound \(v_{\rm eff}=c_0\)
  and the signed price \(D=3c_0/8\) give target
  exponent \(35c_0/24<2c_0-u=3c_0/2\).
* With \(P=1\), all four multiplier masses zero and no
  opposite or cross overlap, take \(u=c_0/2\),
  \(\rho=0\), \(\kappa=D=21c_0/8\) and
  \(z_C=z_{C,P}=3c_0/2\). Use a raw shift cap allowing this
  \(\kappa\). The adaptive shift cost is again \(9c_0/8>2u\),
  while both original moving caps still give half-cap
  \(c_0>u\). Here the original target pool alone gives
  \(v_{\rm eff}=c_0/4\) and target exponent \(35c_0/24\).
  More generally this target-only gain gives an open
  range beyond the adaptive shift threshold when
  \(z_C>2(c_0-u)\) and \(\kappa\) lies sufficiently close
  above \(2u+z_C\).

These are strict comparisons of counting inequalities on
the stated parameter sectors. They do not assert that
large original fibres populate those sectors. The residual
obligation is to force abundant pairs meeting the signed
shift and completion criterion, or to count targets whose
pairs avoid it.

