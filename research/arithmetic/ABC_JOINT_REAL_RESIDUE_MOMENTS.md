# Joint moments for real intervals and unit residues

On an exact valuation box the real logarithmic depth and
the congruence depth add in a residue-list moment. The
estimate applies to arbitrary coupled lists of residue
and interval conditions. It gives joint real and
congruence nonconcentration for actual complete values
and for the quotients left after recording full factors.

## An arbitrary list of congruence and logarithmic bins

At fixed \(L>0\) put \(n=X^L\) and consider a root box
\[
        U=c\prod_{i=1}^m x_i^{e_i},\qquad
        1\le x_i\le n^{a_i},\qquad r=\sum_i a_i ,
\]
with fixed nonzero integer orders, fixed finite lists
and positive rational coefficient \(c\) of polynomial
height. The represented values may be integers or
rationals. All coordinates and \(c\) are units at
\(D\le X^B\). We permit \(D=1\), of depth zero. Write
\(\sigma=\log_n D\). Partition the real logarithm axis
into intervals of length \(\asymp n^{-\xi}\), where
\(0\le \xi\le B/L\); arbitrary common shifts of the
equal-width bins are allowed. Test \(U\) against a
prescribed collection of \(K\) pairs consisting of a
unit residue and a log bin. Replacing these bins by a
list of \(K\) arbitrary intervals of length
\(O(n^{-\xi})\) per indicated residue has the same
bound, since each interval meets \(O(1)\) grid bins.

For \(s=\sigma+\xi\) let \(g_{i,s}\) interpolate
\(g_{i,s}(0)=0\) and
\(g_{i,s}(1/h)=\min(a_i,s/h)\) for positive integers \(h\),
as in SMI1. For \(1\le v\le 2\) put
\[
  G_v(s;\mathbf a)
  =\max_{\substack{w_i\in[0,1]\\\sum_i w_i\le v}}
                 \sum_i g_{i,s}(w_i).
\]
In particular \(G_1\) is the gain of SMI1.

**JRM1.** The list count is at most
\[
     O\!\left(K^\vartheta
             n^{\,r-G_{2\vartheta}(\sigma+\xi;\mathbf a)/2}
                         X^\varepsilon\right)
       \qquad(1/2\le\vartheta\le1,\quad\varepsilon>0).
\]
The constant is uniform in \(D,c,\xi\), in the allowed
unit root alphabets, and in the list, with all fixed
budget and box data as above. It bounds the larger
Cartesian root box, so also applies to correlated
subsets and to distinct complete values represented
by that box.
[Interpolated exact-layer moments for residues alone (SMI1)](ABC_SHALLOW_MOMENT_INTERPOLATION.md#a-residue-list-bound-from-interpolated-moments)
· [Real logarithmic moments (RSI1)](ABC_OCCURRENCE_REAL_PAIR_SPACING.md#a-uniform-relative-strip-bound-for-actual-pairs)

### The list majorant

Choose a period \(P=C_0\log X\) containing the full
tested log range in an interval of length less than
\(P/2\). Enlarge \(C_0\) as needed for the moment orders
below, and put \(N=\lceil P n^\xi\rceil\). It suffices
to use bins of width \(P/N\); a fixed-comparable
interval convention changes the number of bins by
only a bounded factor in the list. On the torus of
period \(P\) write
\[
   K_J(t)=\left|\frac1J\sum_{b=0}^{J-1}
                       \exp(2\pi i b t/P)\right|^2 .
\]
The translates \(K_N(t-jP/N)\), \(0\le j<N\), sum to
one and each has mean \(1/N\). On the bin centered at
\(jP/N\), its translate is at least \((2/\pi)^2\) by
the sine formula. For a finite list \(\mathcal L\)
of unit-residue and circular-bin pairs, the sum of
the corresponding translates for each unit residue
is therefore a majorant up to a fixed factor. Denote
this majorant by \(F(a,t)\). It satisfies
\[
    0\le F\ll1,\qquad
    \int F\ll \frac K{\varphi(D)N}
\]
under probability Haar measure on the finite unit
group and the log torus. The bins of the original
real range have distinct circular indices. The
Fourier frequencies of \(F\) have \(|k|<N\).

Pair the Fourier coefficients of \(F\) with the
product of the slot sums
\[
        S_i(\chi,k)
          =\sum_{x}\chi(x^{e_i})
                         \exp(2\pi i k e_i\log x/P).
\]
The known \(c\) contributes a unit phase. On the
Fourier side use Hölder with conjugate exponents
\(1/(1-\vartheta)\) and \(1/\vartheta\).
Hausdorff--Young on the compact product group gives
\(\|\widehat F\|_{1/(1-\vartheta)}
 \ll(K/(\varphi(D)N))^\vartheta\).
For the other factor put uniform probability on
characters and use the probability weights
\((2N)^{-1}(1-|k|/(2N))\) for \(|k|<2N\).
These weights dominate \(1/(4N)\) on \(|k|<N\).
Converting the Fourier sum norm to this probability
norm costs \(O((\varphi(D)N)^\vartheta)\).
Thus the count is bounded by
\(O(K^\vartheta)\) times the \(L^{1/\vartheta}\)
norm of the product of the slot sums. At
\(\vartheta=1/2\) this is just Parseval and
Cauchy--Schwarz on the list.

### Joint moments in a slot

For an integer \(h\ge1\), expand
\(\mathbb E_{\chi,k}|S_i(\chi,k)|^{2h}\) with the
preceding probability weights. The two collapsed
products \(u,v\le n^{h a_i}\) must satisfy
\(u^{e_i}\equiv v^{e_i}\bmod D\), by character
orthogonality. Their real factor is
\(K_{2N}(e_i\log(u/v))\). Choose \(P\) at least four
times the maximum of \(|e_i|h a_i\log n\) for every
moment used. Then no alias with another period is
possible. The sine formula gives
\[
 K_{2N}(e_i\log(u/v))
           \ll (1+n^\xi|u-v|/n^{h a_i})^{-2},
\]
where \(|\log(u/v)|\ge |u-v|/n^{h a_i}\).
For fixed \(u\) the eligible \(v\)'s lie in
height-subpower many residue classes modulo \(D\),
from the \(e_i\)-th roots of unity. On each such
progression the sum of the last upper bound is
\(O(1+n^{h a_i}/(Dn^\xi))\). This follows by selecting
a nearest progression point and bounding the two
tails by the integral of the decreasing square
kernel at spacing \(D\). Ordered decompositions of
a collapsed product cost height-subpower at fixed
\(h\). Thus
\[
    \|S_i\|_{2h}
       \ll n^{a_i-\min(a_i,(\sigma+\xi)/h)/2}X^{o_X(1)}.
\]
As in SMI1, interpolate adjacent reciprocal indices,
then use Hölder to bound the product norm when
\(\sum_i w_i\le 2\vartheta\).
If a reciprocal index is too small, replace it by
zero using the sup norm. Since \(g_{i,s}(w)\le sw\),
a fixed small cutoff at each budget makes the lost
exponent fit inside \(X^\varepsilon\), uniformly for
\(\sigma+\xi\le 2B/L\). Now only bounded integer
moment indices are needed at the fixed budget,
so one \(P=O_{L,\varepsilon,\rm box}(\log X)\) works
throughout. This proves JRM1.

A complementary elementary estimate fixes all slots
but \(i\). One real bin places \(x_i\) in an interval
of length \(O(n^{a_i-\xi})\), by the logarithmic
derivative. The powered residue leaves height-subpower
many torsion classes for that root. Summing the
result gives \(O(K n^{r-\min(a_i,\sigma+\xi)}
X^\varepsilon)\), useful when one slot has mass
at least the entire combined depth.

## Joint collision entropy and conditional full factors

**JRM2.** Let an actual complete-value law in the box
have maximum mass at most \(n^{-\alpha}\). Its mixed
residue and log-bin variable \(\mathcal B(U)\) satisfies
\[
 \Pr(\mathcal B(U)=\mathcal B(U'))
      \ll n^{\,2(r-\alpha)-G_1(\sigma+\xi;\mathbf a)}
                                      X^\varepsilon .
\]
Here \(U'\) is an independent draw. In fact apply
JRM1 with \(\vartheta=1/2\) to the \(j\) largest cells
to bound the \(j\)-th mass by
\(j^{-1/2}n^{r-\alpha-G_1/2}X^\varepsilon\).
The number of bins and residues is polynomial in \(X\),
so the resulting harmonic sum is height-subpower.
Shannon entropy is at least the corresponding order-two
Rényi entropy. For an arbitrary law with a record \(R\)
fixing the coefficient, modulus and residual cap
box, conditional dyadic rank binning instead gives
\[
 H(\mathcal B(U)\mid R)\ \ge\
       2H(U\mid R)-2\mathbb E r_R
        +\mathbb E G_1(\sigma_R+\xi_R;\mathbf a_R)-o(1).
\]
Entropies are in base \(n\). At a fixed budget the bin
index costs \(O(\log\log X)/\log n\). Uniformity of JRM1
permits arbitrary bins for the actual full records
within the prescribed cap and height bounds.
[Conditional residue rank argument (LFR3)](ABC_BALANCED_ORIGINAL_FIBRE_LIFT.md#full-factor-prices-in-the-shared-residue-register)

For example, at balance and zero complete multiplicity,
the combined depth \(\sigma+\xi\le d/2\) has its full
entropy, where \(d=1/3\), even with macroscopic exact
layers. With fine layers this holds up to depth \(d\).
For \(d\) below the combined depth the entropy
reaches \(d\) in the fine-layer case, by the arbitrary
cap estimate in SMI2. The upper bound here uses the
actual row count \(n^{d+o(1)}\) as well as the mixed
alphabet cap \(n^{\sigma+\xi+o(1)}\). Each of the two
marginals is obtained by using just one component,
so the two are independent at the normalized
entropy scale when the combined depth is fully priced.

After recording any full factor \(F\) of captured
mass \(f\) at zero multiplicity, \(H(U\mid F)=d-f+o(1)\),
as in SMI4. If necessary, bin the sizes of the residual
exact layers at subpower cost at each budget. The
coefficient translates the log bins and multiplies
the unit residue by a known unit. JRM2 gives conditional
entropy \(\sigma+\xi+o(1)\) up to depth \((d-f)/2\)
for arbitrary residual layers. It gives
\(\min(d-f,\sigma+\xi)+o(1)\) for fine layers.
Thus a raw-tight factor record can lose its entire
entropy from the complete letter and still reveal no
normalized information about the mixed probe at the
indicated depths. Fixed lists of compatible probes
on the same role may be combined by least common
multiple before adding the real depth.
[Conditional depth of a full-factor record (SMI4)](ABC_SHALLOW_MOMENT_INTERPOLATION.md#conditional-depth-after-recording-a-full-factor)

## A power saving for residue-dependent short windows

**JRM3.** Work in a balanced chart of true primitive
rows, keeping all three original radical caps.
Let \(U\) be one free role with raw cap \(r=d+o(1)\),
\(d=1/3\), and maximum complete degree \(n^{u+o(1)}\).
Let the chart cost be \(n^{b+o(1)}\) and let a
fixed unit modulus of depth \(\sigma\) define the
residue label on \(U\). The following bound holds
for arbitrary finite lists of windows depending on
the residue, of width \(O(n^{-\xi})\) in \(\log U\).
If the total number of allowed residue–window pairs
is at most \(n^\omega\), the count exponent, after
summing charts with that record cost, is at most
\[
            b+r+u-\tfrac12
                 \bigl(G_1(\sigma+\xi;\mathbf a)-\omega\bigr)
\]
up to the chart, box and counting allowances.
This follows by applying the list count on \(U\)
and restoring its complete degree in the original
graph. The interval centres may be arbitrary
functions of the residue and other fixed records.
For example, if \(\omega\le\sigma+\xi-\gamma\) for
fixed \(\gamma>0\), the chart has a power saving
for \(b+u<\gamma/2\) by a fixed margin, when
\(\sigma+\xi\le d/2\) for arbitrary layers, or
\(\sigma+\xi\le d\) for fine layers. This covers
arbitrarily large collections of separate windows,
without assuming an affine description of their
centres. It supplies a uniform requirement on any
near-zero-multiplicity surviving profile. It leaves
the question of obtaining a small list on the
actual three-cap charts.

With three locks truncated to depth \(d/2\), take
a complete role whose valuation layers are fine.
At \(u=0\) its log bin of depth \(d/2\) can replace
either of the two opposite residue registers in
a pair determining that complete role at the
exponent scale. For \(u>0\) retain the
\(n^{2u+o(1)}\) pair-multiplicity allowance. Apply JRM2 on each
of the finitely many pairs and delete heavy
cells before the same pair-injective selection as
in SRL5. The resulting mixed information is
within each complete role. A bound on simultaneous
real and congruence lifts of all three complete
roles remains necessary to lower the unrestricted
coefficient.
[Three actual shallow registers (SRL5)](ABC_SHALLOW_RESIDUE_REGISTER_SATURATION.md#a-three-register-partial-latin-square)

## Optimized list prices after full-factor records

Use the original primitive row law and all three radical caps.
The complete degree on the tested role is at most \(n^{u+o(1)}\).
A joint record \(R\), priced at \(n^{b+o(1)}\) states, fixes
the indicated full factors of the tested value \(U\).
It may also fix factors of the other two letters. Include
the divisor choices and any fixed coefficients in that record.
In each state the remaining exact root box has raw caps
\(\mathbf a_R\) of sum \(c_R\). The complete value is a
known positive coefficient times the residual root product.
Use the logarithmic bins of the **remaining numerical value**,
or translate the bins for \(U\) by the log of the fixed
coefficient. Allow a unit modulus \(D_R\) of polynomial
height chosen by the record. A sign prescribed by an
original congruence only permutes its unit labels.

**JRM4.** In a bin of these records with common upper
residual raw cap \(c\), list size at most \(n^k\), and
combined residue and logarithmic depth \(s\), the row
count after summing the records has exponent at most
\[
 b+u+\min\left\{c,\ \inf_{1/2\le\vartheta\le1}
       \left(c-\tfrac12 G_{2\vartheta}(s;\mathbf a)
                        +\vartheta k\right)\right\}
\]
up to the count and bin allowances. The displayed
profile gain is used with a common certified lower
bound on \(G_{2\vartheta}(s_R;\mathbf a_R)\) in the bin.
There is also a version for genuinely different
records. If \(\mathcal B\) is the mixed label,
\[
 H(U\mid R)\le
 \mathbb E\!\left[c_R-\tfrac12G_{2\vartheta_R}
                         (s_R;\mathbf a_R)\right]
       +\sum_{r}\Pr(R=r)\,\vartheta_r
                      H(\mathcal B\mid R=r)+o(1),
\]
for \(\vartheta_R\in[1/2,1]\) chosen by the record on
fixed finite parameter meshes. The original row entropy
exceeds \(H(U\mid R)\) by at most \(H(R)+u+o(1)\).

The first assertion is JRM1 followed by the original
complete degree bound. For the entropy version, rank the
actual conditional label probabilities in each record
state and bin the ranks dyadically. If the rank is \(j\),
the normalized log of the upper size of its bin is
at most \(\log_n j+\log_n 2\), which is at most the
realized conditional surprise plus \(\log_n 2\).
Apply JRM1 on each rank bin and average. The entropy
of the rank bin itself is height-subpower. The last
assertion follows from the chain rule and from the
original degree bound, which remains valid after
conditioning. The argument uses the joint record
once, even if its prices or its moduli depend on
several original letters.

### Explicit loss at fine residual caps

Suppose \(c_R=c+o(1)\), \(s_R=s+o(1)\), with fixed
positive \(c,s\), and \(\max_i a_{R,i}=o(1)\) in the
height-first convention. Then the counting exponent
in JRM4 is at most
\[
 b+u+c-\ell(c,s)\bigl(\min(c,s)-k\bigr)_++o(1),
 \qquad
 \ell(c,s)=\begin{cases}
 1& s\le c/2,\\
 c/(2s)& c/2\le s\le c,\\
 1/2& s\ge c .
 \end{cases}
\]
In fact, if \(a_*=\max_i a_i\) and \(T=(s-a_*)_+>0\), take
\(h_i=\lceil T/a_i\rceil\) on each positive-cap slot.
At \(w_i=1/h_i\), its gain is \(a_i\ge T w_i\), because
\(h_i a_i\le T+a_i\le s\). If the total weight is at
most \(2\vartheta\), use every slot. Otherwise scale
the weights proportionally to total \(2\vartheta\) and
use concavity. This proves \(G_{2\vartheta}(s;\mathbf a)
\ge\min(c,2\vartheta T)\); the reverse bound with \(s\)
uses the total cap and the slope bound. Minimizing
\(\vartheta k-\min(c,2\vartheta s)/2\) on
\([1/2,1]\), together with the trivial raw price,
gives the formula. The choices of depth, count margin
and permitted maximum slot cap precede the height
limit. When \(s\le c/4\) with arbitrary caps, the same
linear loss \(s-k\) holds for \(k\le s\). If a slot
has cap at least \(s\), apply the elementary one-slot
test following JRM1. Otherwise
\(\sum_i 1/\lceil s/a_i\rceil\ge c/(2s)\ge2\);
the full linear moment gain \(G_2=2s\) is then
available at \(\vartheta=1\).

For example, suppose the tested role has original cap
\(d=1/3\), and the joint record captures full mass
\(f=d-c\) on it. In the fine-layer range \(s\le c/2\)
the rows have a fixed power saving if
\[
                         s-k>b-f+u
\]
by a fixed margin. The charge on the right can be
arbitrarily small at near-zero multiplicity for a
raw-tight record. At the exponent scale the simpler
square-root list test would require
\(s-k>2(b-f+u)\). The improvement applies even
when different record states choose different probes.
Equivalently, on a full-exponent law with
\(\Delta=c-H(U\mid R)\) and common fine \(c,s\),
the conditional mixed entropy is at least
\(\min(c,s)-\Delta/\ell(c,s)-o(1)\).
In particular shallow list entropy loses at most
\(\Delta\) for \(s\le c/2\), compared with the
\(2\Delta\) loss from the square-root estimate.
These tests still require an actual list deficit
to exclude a chart. The strict fixed-shift estimates
supply no such deficit solely from a broad valuation
spectrum or from a common full divisor.
[Optimized fixed-shift and raw-floor tests](ABC_BALANCED_FIXED_SHIFT_REVERSE_TESTS.md#sharp-raw-active-floor-envelope)


### Collision and pointwise versions

In one record state, assume the maximum mass of a
complete value is at most \(n^{-a}\). With
\(\beta=G_{2\vartheta}(s;\mathbf a)/2-(c-a)\),
JRM1 applied to the \(j\) largest label probabilities
gives
\[
      p_j\le n^{-\beta+o(1)}j^{\vartheta-1},
           \qquad p_j\le 1/j .
\]
When \(\beta>0\), sum the squares of the first bound
up to \(j=n^{\beta/\vartheta}\) and use the second
bound thereafter. The two sums are
\(O(n^{-\beta/\vartheta+o(1)})\), with a harmonic
factor at \(\vartheta=1/2\). This proves a collision
entropy lower bound \(\max(0,\beta/\vartheta)-o(1)\),
uniformly with the inner counting allowance.
At \(\vartheta=1\) the same lower bound holds for
min-entropy. Thus the explicit fine-layer bounds
on the conditional entropy also hold for collision
entropy with \(\Delta=c-a\) when the pointwise complete
mass bound is available in that state.

In particular for \(s\le c/2\) in the fine-layer
case, a single pair consisting of the full-factor
record and a mixed label carries at most
\(n^{u+c-s+o(1)}\) original rows. Any pooled list of
such pairs of exponent \(\omega\) then has row
exponent at most \(\omega+u+c-s+o(1)\).
This needs no uniform bound on the number of
allowed labels per record. For \(c=d-f\) it saves
if \(\omega<f+s-u\) by a fixed margin. With arbitrary
layer caps the same assertion holds at \(s\le c/4\).
[Full-exponent three-letter chart and degrees](ABC_BALANCED_DICTIONARY_RESIDUE_CELLS.md#regular-balanced-cohorts-and-absolute-residue-pricing)
· [One-label loss at arbitrary shallow caps (SGM1)](ABC_SECANT_GCD_MEAN_OVERLAPS.md#a-one-residue-count-at-low-combined-depth)


## Fine residual caps from small individual primes

**JRM5.** Fix \(a>0\) and take an exact-valuation box at
fixed \(L\). It suffices for the fine-layer estimates
that each varying individual prime is at most \(n^a\).
An exact valuation layer of large mass may be divided
into subroots of caps at most \(2a+o_X(1)\), with total
cap equal to the residual radical cap up to \(o_X(1)\).
The number of subroots is bounded at each fixed budget
and fixed valuation truncation.

Indeed, within each integer order put the distinct
prime factors in increasing order. Greedily multiply
consecutive primes into a squarefree block until the
normalized size first reaches \(a\), then start a new
block. Every completed block has mass in \([a,2a]\)
and there is at most one smaller block per order.
A total raw cap \(r\) therefore bounds the number of
completed blocks by \(r/a\). Preserve the original
order on every block. Fix their number per order and
dyadically bin their actual sizes. The fixed number
of blocks makes this a height-subpower refinement;
the sum of the new upper caps is at most
\(r+O_{L,a}(1/\log n)\). Padding by unit roots
is harmless. For any fixed full-factor record do
this on its remaining exact factors. The largest
original valuation layer imposes no restriction on
the new moment cap.

The factors at orders above \(K/L\) have total raw mass
at most \(H/K\) on values \(\le X^H\). They may be placed
in the full record before this argument, at that
arbitrarily small raw cost. Both the order truncation
and the required size \(a\) are fixed before the height
limit. Different fixed records and coefficient choices
translate the real bins and the unit residues exactly
as in JRM4.
[Exact-layer moments and conditional records](#optimized-list-prices-after-full-factor-records)

### Retaining the nonatomic mass after large-prime fixing

On a balanced sequence with scaled valuation measure
\(\nu\) of mass \(d=1/3\) on the tested role, record
the full powers on primes larger than \(n^a\), with
\(a>0\) fixed. The recorded mass \(f\) satisfies
\(\limsup f\le d-\nu_{\rm na}([0,\infty))\), before
the arbitrarily small extra valuation-tail charge.
Indeed there are at most \((d+o(1))/a\) such primes.
Each has scaled order \(L v_p\le H/a\) on values up
to \(X^H\). Along a subsequence the sizes and scaled
orders of these finitely many marked primes converge,
after padding with zero masses. Their limiting
submeasure is atomic and is dominated by \(\nu\).
This proves the inequality, including when a marked
order tends to zero. The argument applies uniformly
on anatomic size boxes converging to \(\nu\), by
taking a violating sequence of actual rows if needed.

A mass bin \(f\) costs at most \(n^{f+o(1)}\) full
records by radical enumeration followed by
fixed-support enumeration. At zero complete
multiplicity an exponent-full chart therefore has
a subchart of exponent at least \(d-f-o(1)\) after
selecting one record, with residual cap \(d-f+o(1)\).
An initial cheap lock and the extra high-order tail
record are charged separately. The residual cap is
at least the original nonatomic mass up to the stated
errors. The upper bound \(2a\) on block masses is
available even when the original chart has
macroscopic atomic exact layers. At a balanced
endpoint at least two roles have nonatomic mass at
least \(1/6\). Hence arbitrarily fine block caps
with at least that residual mass can be obtained
on a tested role, while a different role supplies
an economical lock. This is a claim about the
retained full-factor chart; selecting its record
can cost a positive exponent in the original family.
[Pairwise nonatomic constraints](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)
· [Uniform enumeration on fixed supports (EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

## Pair mass in large and small mixed cells

**JRM6.** On a fixed full-factor chart of actual rows,
let one complete role have residual cap \(r>0\), maximum
slot cap \(a_*\), and original complete-degree cap \(n^u\).
Assume its uniform row law has size \(N\ge n^{r-\eta_0}\).
Take a unit modulus and logarithmic bins of combined
depth \(0<s<r\), at fixed distance from \(r\). Let
\(p_j\) be the mixed cell probabilities in decreasing
order and put
\[
 \vartheta=\min(1,r/(2s))>1/2,\qquad
 \mu=u+\eta_0+\vartheta a_* .
\]
All displayed exponents may absorb an arbitrary fixed
positive count and box allowance. By JRM1 and the
fine-slot gain preceding JRM4,
\[
                     p_j\le
 n^{-\vartheta s+\mu+o(1)}\,j^{\vartheta-1}.
\]
Indeed each value has probability at most \(n^u/N\)
and \(G_{2\vartheta}\ge 2\vartheta(s-a_*)_+\), with the
inequality weakened to zero if \(a_*\ge s\).

For fixed \(0<\rho<s\), the cells of probability at least
\(n^{-s+\rho}\) together have squared mass at most
\[
                n^{-s+2\mu-(2\vartheta-1)\rho+o(1)}.
\]
There are at most \(n^{s-\rho}\) such cells. Sum the square
of the rank bound up to this rank, using
\(\sum_{j\le J}j^{2\vartheta-2}\ll(1+\log J)J^{2\vartheta-1}\).
The cells of probability at most \(n^{-s-\rho}\)
have total squared mass at most \(n^{-s-\rho}\)
by the trivial \(L^1\) bound. The total squared mass
is at least \(n^{-s-o(1)}\) by the number of possible
joint labels. Thus if
\[
                  2\mu<(2\vartheta-1)\rho
                  \quad\hbox{and}\quad r-\eta_0>s
\]
by fixed margins, all but a power-small fraction
of same-cell ordered distinct-row pairs fall in
cells with between \(N n^{-s-\rho}\) and
\(N n^{-s+\rho}\) rows. The second margin lets one
subtract the diagonal \(N\) from the pair count.
Repeated complete values on the tested role account
for at most \(N n^u\) ordered pairs, which are also
power-small when \(r-\eta_0-u>s\). The upper bound
on the discarded portion follows from the squared
mass and so controls the pair population, rather
than merely the mass of individual rows.

These statements use a fixed modulus and the actual
mixed cells within each record. For the near-occurrence
application take \(\rho\), the depth margin and the
desired count margin first, then take the budgets
and individual prime caps sufficiently small, fixing
large primes as in JRM5 when needed. Every use of
the moment estimate takes the height limit before
the small-budget limit. The saving deteriorates at
\(s=r\), where the reciprocal exponent is \(1/2\)
and the upper-tail argument has no power margin.
