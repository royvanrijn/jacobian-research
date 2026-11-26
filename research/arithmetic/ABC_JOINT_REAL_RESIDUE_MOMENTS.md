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
