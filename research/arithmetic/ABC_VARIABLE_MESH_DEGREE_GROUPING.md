# Grouping valuation degrees on a variable mesh

At the normalized \(1/3\) endpoint, the nonatomic masses of any
two letters sum to at least \(1/3\). A single grouped root detects
positive excess over this threshold from either a numerical height
gap or high valuation degrees. The mesh can depend on the limiting
measures. Its spacing may be arbitrarily small, while the captured
fraction of nonatomic mass has a lower bound independent of that
spacing.

## Bounds from the pair excess

Use the [limiting profile notation](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees).
Let \(\nu_T\) be the weak scaled-degree measure and \(c_T\) its
nonatomic mass on letter \(T\). Let \(h_T\) be the limit of the first
moments of the prelimit \(T\)-measures, equivalently the numerical
height exponent. Pass to a subsequence for these moment limits.
Define
\[
 \Delta_{UV}=c_U+c_V-\tfrac13,\qquad
 e_T=h_T-\int_0^\infty r\,d\nu_T(r)
\]
for distinct letters \(U,V,T\). At \(d=1/3\) the
[existing endpoint conditions](ABC_WIDE_COFACTOR_CONGRUENCES.md#necessary-conditions-at-normalized-exponent-13)
give total root mass \(1/3\) on each letter and \(\Delta_{UV}\ge0\)
for every pair. Truncating the first moments gives \(e_T\ge0\).
Normalize the maximum numerical height to \(1\).

There are absolute constants \(C,R_*>0\) with the following bounds
for every such profile.

- If \(U,V\) are the two tall letters and the short height is \(h_T\),
  then
  \[
   1-h_T\le C\Delta_{UV}.
   \tag{VM1}
  \]
- For every pair \(U,V\) and its complementary letter \(T\),
  \[
   e_T+\int_{[R_*,\infty)}r\,d\nu_T(r)\le C\Delta_{UV}.
   \tag{VM2}
  \]

In particular, saturation on a pair forces the complementary letter
to have bounded scaled-degree support and convergent first moments.
Saturation on the tall pair forces the short numerical height to be
\(1\).

### One grouped root on the free pair

Fix \(0<\theta\le1/8\), and write \(c=c_U+c_V\). We will freeze
coordinates on \(U,V\) at normalized cost at most \(c-t+t/16\), for
some \(t\ge\theta/60\). The remaining ratio, up to the fixed frozen
factor, will have the form \(V_0 Z^q\), with a spacing \(0<b\le1\)
and \(q=\lfloor b/L\rfloor\) such that
\[
 H(Z)\le X^{(3/b+o(1))L},\qquad
 H(V_0)\le X^{2\theta b t+o(1)} .
 \tag{VM3}
\]
Here \(Z,V_0\) are positive rationals. The limit errors in (VM3) refer
to the sequence of boxes, with \(\theta,t,b\) fixed.

At least one of the two letters has nonatomic mass at least \(1/6\).
On that letter choose a compact interval \(J\subset(0,\infty)\) with
continuity endpoints whose nonatomic mass \(m_J\) is at least \(1/12\).
Set
\[
 M=\lceil 4/\theta\rceil,\qquad t_0=m_J/M .
\]
Then \(t_0\ge\theta/60\). Choose finitely many positive atoms on the pair
whose omitted positive atomic mass is at most \(t_0/16\). Approximate
all selected scales \(r_i\) simultaneously on a mesh of spacing \(a\):
\[
 \|r_i/a\|<\theta t_0/(16M).
\]
Here \(\|x\|\) denotes distance to the nearest integer. Such \(a\) can
be arbitrarily small. For example, take \(a=1/n\) with an arbitrarily
large positive return time of the vector of fractional parts
\((nr_i)_i\) to this neighborhood of zero. The pigeonhole proof of
simultaneous approximation gives these return times by considering
differences of arbitrarily widely spaced indices in a single cell of
a finite partition of the unit cube. Take \(a<1\), below half the lower
endpoint of \(J\), and below half every selected atom scale.

For each \(r\in J\), the \(M+1\) fractional parts of \(i r/a\),
\(0\le i\le M\), give a difference \(1\le j\le M\) for which
\(\|j r/a\|\le1/M\). Thus some \(j\le M\) captures nonatomic mass
at least \(t_0\) in the set \(\|j r/a\|\le1/M\) within \(J\).
Set \(b=a/j\). Enlarge the normalized radius of these windows to a
number between \(\theta/2\) and \(\theta\) whose endpoints within
\(J\) carry no mass for the full measure of the selected letter.
Such a choice avoids at most countably many radii. The resulting
union within \(J\) has finitely many intervals, with nonatomic
mass \(t\ge t_0\). Every retained selected atom on either free
letter satisfies \(\|r_i/b\|<\theta t_0/16\).

Keep the coordinates of the window union on the selected letter,
the shrinking neighborhoods of all the selected atoms on both
letters and both shrinking zero windows. For the atom neighborhoods,
take the weak-convergence diagonal that recovers the respective atom
masses. If an atom neighborhood meets the union on the selected letter,
it eventually lies in the same mesh cell and uses the same mesh
point. Consequently the retained normalized mass is at least the
total atomic mass, less \(t_0/16\), plus \(t\). This gives the
asserted frozen cost.

For each retained positive coordinate of original degree \(e\),
write \(e=iq+(e-iq)\), where \(ib\) is its assigned mesh point.
Put the root to power \(i\), with sign determined by the side of
the ratio, in \(Z\). Put its remaining power in \(V_0\), along with
all the zero-window factors. Since \(b\) and the finite collection
of cells are fixed, the error in using \(q=\lfloor b/L\rfloor\) has
vanishing total numerical height. The sum of the weights \(i\) times
their retained root masses is bounded by \(2/b+2/3\le3/b\),
using \(i\le r/b+1\), total pair first moment at most \(2\), and
total pair root mass \(2/3\). This proves the narrow height bound.

In units of \(b\), the window union costs at most \(\theta\)
per unit of nonatomic mass and per unit of unlisted atomic mass.
The contribution from all selected atoms is at most
\(\theta t_0/16\) times their mass, which is at most \(2/3\).
Together with the vanishing zero-window contribution, the error is
at most
\[
 b\bigl(\theta(t+t_0/16)+(\theta t_0/16)(2/3)\bigr)+o(1)
 \le 2\theta b t+o(1).
\]
This gives the second bound in (VM3).

For each fixing of the frozen factors, a fixed exact value
\(V_0Z^q\) gives the complete ratio of the two free letters.
Their coprimality recovers the two magnitudes, and the signed
equation recovers the third. A given triple has divisor-many
coordinate representations at fixed truncation, so the exact
fibers in the grouped sample cost only \(X^\varepsilon\). This
permits the unamplified forms of both counting estimates with
possibly overlapping prime supports for \(Z,V_0\).

### Using a numerical height gap

Suppose \(0<\tau<1-h_T\). The complete tall ratio is
\(1+O(X^{-\tau})\), by the equation and the two tall heights.
After dividing by the frozen ratio this gives the strip target
for the [unequal-height real bound (RUW1)](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-real-strip-with-unequal-root-heights).
Let \(c_{\mathrm r}\) be its constant in dimension one. Set
\[
 \eta_{\mathrm r}=\min(1/8,c_{\mathrm r}/256),\qquad
 \theta=\eta_{\mathrm r}\tau .
\]
Apply (VM3). Use \(m=1\), \(s=t/8\), and narrow height
constant \(C_1=4/b\) in (RUW1), so its height product is
\(32/(bt)\). The wider exponent \(2\theta bt\) is strictly below
\(c_{\mathrm r}\tau bt/32\), leaving room for the height errors.
The count in the strip costs \(2s=t/4\) in the normalized exponent,
and the already frozen factors cost at most \(c-t+t/16\).
Therefore
\[
 \tfrac13=d\le c-\frac{11t}{16},\qquad
 t\ge \eta_{\mathrm r}\tau/60.
\]
Letting \(\tau\) increase to \(1-h_T\) proves (VM1).

### Counts near a fixed ratio

The same real argument counts a subfamily for which the complete
ratio of **any** two letters, divided by a fixed positive target
\(K\), differs from \(1\) by \(O(X^{-\tau})\), for fixed
\(0<\tau<1\). The bound in (RUW1) is uniform in \(K\), and its
application above only uses this strip condition. Thus in boxes
with limiting profile at \(d=1/3\), the subfamily count is
\[
 X^{(c_U+c_V-k_0\tau+o(1))L}
\]
up to arbitrary fixed-truncation exponent losses, with an absolute
\(k_0>0\). One can take \(k_0=11\eta_{\mathrm r}/960\).

At pair saturation this gives a normalized saving from \(1/3\)
for each fixed target, including a target depending on the height
or box. Taking the union over positive rational targets of
height at most \(X^{\gamma L}\) costs at most \(2\gamma\) in the
normalized exponent. In particular the saving persists when
\(2\gamma<k_0\tau\). For equally tall addends the sector
\(|a-b|\le X^{1-\tau}\) is contained in the strip with target \(1\)
and precision any fixed \(\tau'<\tau\), by their numerical heights.

### Locking a high-degree tail

Let \(R>0\) be a continuity point for \(\nu_T\). For the prelimit
\(T\)-measure freeze all whole coordinates of scaled degree at least
\(R\). This costs
\(\mu=\nu_T([R,\infty))\) in normalized mass and gives a modulus with
limiting numerical weight
\[
 w=h_T-\int_{[0,R)}r\,d\nu_T(r)
   =e_T+\int_{[R,\infty)}r\,d\nu_T(r).
\]
The first equality uses convergence of the full moments and of the
bounded truncated moments. In particular \(R\mu\le w\le1\).
If \(w=0\), the desired estimate holds. For \(w>0\), let
\(c_{\mathrm p}\) be the dimension-one constant of the
[unequal-height modulus bound (UW)](ABC_QUANTITATIVE_DEGREE_WINDOWS.md#a-product-bound-for-the-admissible-wider-height).
Take
\[
 \eta_{\mathrm p}=\min(1/8,c_{\mathrm p}/256),\qquad
 \theta=\eta_{\mathrm p}w .
\]
Apply (VM3) with this \(\theta\). Use \(kL\to\kappa=8w/t\)
for the modulus. Since \(t\ge\eta_{\mathrm p}w/60\), this choice
is possible whenever \(R>480/\eta_{\mathrm p}\). In (UW) take
\(m=1\), narrow height constant \(C_1=32/(bt)\), and wider exponent
just above \(2\theta bt/w\) relative to the modulus. The narrow
bound follows from \((3/b)\kappa/w=24/(bt)\). The wider bound
follows from
\[
 (2\theta bt/w)C_1=64\eta_{\mathrm p}<c_{\mathrm p}.
\]
The signed equation gives a unit-target congruence for \(V_0Z^q\),
as the two free letters are coprime to the whole locked modulus.
Its counting cost is \(2w/\kappa=t/4\).

Take \(R\ge960/\eta_{\mathrm p}\) so that \(\mu\le t/16\).
Adding the frozen costs on all three letters gives
\[
 \tfrac13=d\le c-t+t/16+\mu+t/4\le c-\frac{5t}{8},
 \qquad t\ge\eta_{\mathrm p}w/60.
\]
Thus \(w\le (96/\eta_{\mathrm p})\Delta_{UV}\). For (VM2),
fix \(R_*=1+960/\eta_{\mathrm p}\) and let continuity points
\(R<R_*\) increase to \(R_*\). The estimate includes any atom
at \(R_*\) and any escaped first moment.

## Consequences at the minimum total nonatomic mass

If \(\Gamma=c_A+c_B+c_C=1/2\), the three nonnegative pair
excesses sum to \(2\Gamma-1=0\). Thus \(c_T=1/6\) on every
letter. By (VM1) all numerical heights are \(1\); by (VM2)
all three weak measures have first moment \(1\), with no escaped
moment, and are supported in the same bounded interval
\([0,R_*]\). In particular the \(1/3\) endpoint at this nonatomic
total requires three full-height letters with bounded scaled
degrees, independently of the number or arithmetic relations of
their positive atomic scales.

These constraints still admit formal measures. For example, give
each letter an atom of mass \(1/6\) at \(3\) and distribute its
remaining mass \(1/6\) uniformly on \([2,4]\). Each first moment is
\(1\) and each pair excess is zero. This example concerns the
displayed necessary conditions on measures.

## Atomic and nonatomic weight near the minimum mass

Write \(\nu_{T,\mathrm{at}}\) and \(\nu_{T,\mathrm{na}}\) for the
atomic and nonatomic parts of \(\nu_T\), and define their numerical
weights by
\[
 A_T=\int r\,d\nu_{T,\mathrm{at}}(r),\qquad
 B_T=\int r\,d\nu_{T,\mathrm{na}}(r).
\]
There are absolute \(\eta,\gamma_0>0\) such that every endpoint
profile with \(0\le\Gamma-1/2\le\gamma_0\) satisfies
\[
 \min(A_U,B_U)+\min(A_V,B_V)\ge\eta
 \tag{VM7}
\]
for each pair of letters. In particular at least two letters
must carry a uniformly positive numerical weight of each type.
After decreasing \(\gamma_0\), both types on these letters also
carry uniformly positive root mass in one fixed compact positive
degree range.

First consider any endpoint at which \(c_U+c_V=1/3\) and the
only possible atom on either \(U\) or \(V\) is at zero. The third
letter has nonatomic mass at least \(1/6\). Choose a compact band
\([b,2b]\), \(b>0\), carrying positive nonatomic mass of the third
letter. On one of the free letters take a patch of small positive
width \(\sigma\), away from zero, with nonatomic mass \(t>0\).
Such patches exist with \(\sigma\) arbitrarily small and \(t\)
tending to zero. Free both zero windows and the patch, and freeze
the rest on \(U,V\) for cost \(1/3-t\). Group the patch at a
representative exponent, with root height \(X^{(t+o(1))L}\)
and wider cofactor height \(X^{\sigma t+o(1)}\).

Inside the band on the third letter, choose a finite union of
continuity intervals of mass \(\mu\in[t/64,t/32]\). This follows by
first selecting a subset of the nonatomic part away from all
atoms, then using inner and outer regularity to cover a compact
subset by intervals. Freeze its factors as a full modulus of
numerical weight \(w\in[bt/64,bt/16]\), with \(kL\to b/2\).
In (UW) take dimension one and narrow cap \(C_1=64\).
The wider-height condition holds for small enough \(\sigma\),
independently of \(t\), since \(\sigma t/w\le64\sigma/b\).
The unamplified count costs at most \(t/4\), and its exact ratio
recovers the primitive triple. Including the two freezes gives
\(d\le1/3-t+t/32+t/4<1/3\). Thus such a saturated pair cannot
have both atomic parts supported at zero.

If (VM7) fails uniformly as \(\Gamma-1/2\to0\), choose two fixed
letters and one of the two types on each whose combined numerical
weight tends to zero. Pass to a subsequence fixing these choices.
The three pair inequalities imply \(c_T\to1/6\) on every
letter, since any pair has sum at least \(1/3\) and the total
tends to \(1/2\). The chosen submeasure of each of \(U,V\) has
mass tending to \(1/6\), and its first moment tends to zero.
It therefore converges weakly to \((1/6)\delta_0\).
Take a weak limit of the full letter measures. As in the
[compact-scale argument](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#positive-scales-in-every-endpoint-letter),
diagonalizing the underlying sequences of boxes makes this weak
limit another \(d=1/3\) profile; one takes the boxes far enough
along each sequence to approximate successively many bounded
continuous test integrals and the limiting normalized count.
Each of the two letters has at least \(1/6\) mass at zero in
the new profile, so its nonatomic mass is at most \(1/6\).
Their pair inequality forces equality on both. Hence neither
has any positive atom, contradicting the preceding patch count.
This proves (VM7). The argument tracks the chosen submeasures
and does not assume continuity of the operation of taking the
atomic part.

At most one letter has \(\min(A_T,B_T)<\eta/2\). On the other
two letters the weight contributed below \(r_0\) by either type
is at most \(r_0/3\). The tail in either type above \(R_*\) is
bounded by (VM2), which is \(O(\Gamma-1/2)\) since every pair
excess is at most \(2(\Gamma-1/2)\). Choose \(r_0>0\) and
\(\gamma_0>0\) so that these two errors sum to at most \(\eta/4\).
Each type then has numerical weight at least \(\eta/4\) in
\([r_0,R_*]\), hence root mass at least \(\eta/(4R_*)\)
there.

## A lower bound for the phase discrepancy of the whole pair

Let \(d=1/3\), and suppose a modulus letter \(T\) has numerical
height at least \(19/20\) and every atom of its degree measure
has mass strictly less than \(1/5\). Then its complementary
pair satisfies
\[
 \int\|r/a\|\,d(\nu_U+\nu_V)(r)>\frac1{24000}
 \qquad(0<a\le1).
 \tag{VM8}
\]
In particular every balanced one-third endpoint has a pair for
which this holds at every such spacing. If \(\Gamma-1/2\) is
sufficiently small, it holds on all three pairs.

Fix \(a\), and let \(D\) be the displayed integral. Round every
prelimit scaled degree \(eL\) of the free pair to the nearest
multiple \(ia\), using either nearest point at a tie. In the signed
ratio put the \(i\)-th power of its root in \(Z\), and the remaining
power \(e-iq\) in \(V_0\), with \(q=\lfloor a/L\rfloor\).
The exceptional coefficients are fixed. Since \(i\ge0\),
each of the numerator and denominator of \(Z\) has height
exponent bounded using that letter's first moment. Thus
\[
 H(Z)\le X^{(1/a+1/6+o(1))L},\qquad
 H(V_0)\le X^{aD+o(1)} .
\]
Here \(i\le(eL)/a+1/2\), each prelimit normalized root mass
tends to \(1/3\), and each numerical height is at most \(1\).
Changing \(ia\) to \(iqL\) costs at most
\(L\sum i\beta_e=o(1)\) numerically, with \(\beta_e\) the
normalized root masses. This uses the first-moment bound even
for escaping prelimit weight. The remaining error is \(a\)
times the integral of the bounded continuous function \(\|r/a\|\),
which converges by weak convergence.

Every individual prelimit root on \(T\) has normalized mass less
than \(1/5\) for small \(L\). Indeed, an atom of mass at least
\(1/5\) would result from a subsequence of such coordinates,
whose scaled degrees are bounded by their first moments.
Take the \(T\)-coordinates of scaled degree at least
\(\kappa=6/5\). They have total numerical weight at least
\(19/20-(6/5)(1/3)+o(1)=11/20+o(1)\).
Group them greedily into blocks of normalized mass at most
\(1/5\), starting a new block when the next coordinate would
overflow the current one. Each consecutive pair of blocks has
combined mass greater than \(1/5\), so there are at most three
blocks. Freeze one with limiting weight \(w\ge11/60\) at
normalized cost \(\mu\le1/5\), passing to a subsequence.
Use its whole factors as the modulus \(n\) with \(kL\to6/5\).
The residue modulus in (QD1) has normalized logarithmic mass
at least \(\mu\), since each locked root contributes all its
prime powers to it. It also has mass at least \(w/\kappa-\mu\)
by the rounding loss for each prime of \(n\). Thus
\(b\ge w/(2\kappa)\) in (QD2).

Assume \(D\le1/24000\). Apply the amplified
[quantitative congruence estimate (QD2)](ABC_QUANTITATIVE_DEGREE_WINDOWS.md#a-quantitative-congruence-estimate-for-one-grouped-root)
to \((Z,V_0)\) with \(\ell=\lfloor8w/\kappa\rfloor+1\).
This gives \(w/(\ell\kappa)<1/8\). We have \(\ell\ge2\) and
\(w/\ell\ge11/120\): for \(\ell=2\) use \(w\ge11/60\),
while for \(\ell\ge3\) use \(w\ge(\ell-1)\kappa/8\).
Take \(r=a/100\), \(M=\lfloor220/a\rfloor\), and
\(\alpha=1/2\). Then \((M+1)r\in(11/5,221/100]\).
For the amplified heights \(h=\ell(1/a+1/6)\),
\(E=\ell a/24000\) one has
\[
 \frac B\ell\le\frac7{600}+\frac{220}{24000}
   =\frac{120}{5760}
   <\frac{121}{5760}
   \le \frac{b(M+1)r}{4\ell}.
\]
The precision and terminal conditions follow from
\((M+1)r/2\le221/200<6/5\) and \(4B/\ell\le1/12<11/120\).
The counted ratio recovers the primitive triple. Including the
modulus freeze gives \(d<1/5+1/8<1/3\), a contradiction.

Among the two tall letters at least one has nonatomic mass at
least \(1/6\), by their pair inequality. It has numerical
height one and atomic mass at most \(1/6\), so it supplies a
modulus letter for (VM8). Near \(\Gamma=1/2\), (VM1) and the
pair inequalities give height at least \(19/20\) and total
atomic mass less than \(1/5\) on every letter.

### A bound for any positive nonatomic modulus mass

The shallow residue modulus of (QD2) also gives a bound for
any balanced one-third profile, including low numerical height.
For any letter \(T\) with \(c_T>0\), its complementary pair
satisfies
\[
 (1+a/6)\int\|r/a\|\,d(\nu_U+\nu_V)(r)
     >\frac{h_T c_T^2}{44000}\qquad(a>0).
 \tag{VM9}
\]
This form uses the unsigned distance to the mesh.

Write \(c=c_T\), \(H=h_T>0\), and \(J=1+a/6\).
The preceding whole-letter grouping has narrow height constant
\(J/a\) and wider exponent \(aD\), where \(D\) is the integral
in (VM9), for every fixed positive spacing \(a\).
Take \(\kappa=H\). On \(T\), the coordinates above this scaled
degree have numerical weight at least \(2H/3+o(1)\).
Every coordinate has normalized mass at most \(1/3-c+o(1)\)
by the weak atom argument. Greedily split these high coordinates
into blocks of normalized mass at most \(1/3-c/3\). There are
at most three blocks since two consecutive blocks have combined
mass exceeding that cap and four blocks would exceed the
total \(1/3\). Lock one with limiting weight \(w\ge2H/9\)
and cost \(\mu\le1/3-c/3\), with \(kL\to H\).

The available full residue mass is at least
\(\max(\mu,w/\kappa-\mu)\ge1/9\). Select a shallow divisor
\(Q\mid Q_n\) by adding copies of whole roots as in the proof
of (QD3), until its mass \(b_Q\) reaches \(1/10\).
The overshoot is at most \(\mu\) in the limit,
so \(b_Q\le1/10+1/3=13/30\).
Use amplification order \(\ell=\lceil 2/c\rceil\) in (QD2),
and assume \(D\le H/(8000\ell^2 J)\).
Take
\[
 r=\frac{aH}{100\ell J},\qquad
 M=\left\lceil \frac{160\ell J}{a}\right\rceil,\qquad
 \alpha=1/2.
\]
Here \(r<a=\lim qL\) and \(a/J\le6\).
The amplified height cost satisfies
\[
 B\le \left(\frac1{100}+\frac{166}{8000}\right)H
       =\frac{123}{4000}H.
\]
The determinant threshold is at least \(H/25\) because
\((M+1)r\ge(8/5)H\), while the precision uses
\((M+1)r/2\le(43/50)H<\kappa\).
Also \(4B\le(123/1000)H<w\).
Thus the count after the freeze costs at most \(b_Q/\ell\le13c/60\),
giving \(d\le1/3-c/3+13c/60<1/3\).
Since \(\ell\le(7/3)/c\), the stated weaker threshold in (VM9)
follows.

## A restriction on concentration near a mesh

Suppose the modulus letter \(T\) has positive nonatomic mass.
At pair saturation this follows from \(\Gamma\ge1/2\).
On one free letter choose a compact continuity interval in
\((0,\infty)\). Fix a spacing \(a\in(0,1]\), a phase
\(\varphi\in(0,1]\), and a radius \(0<\theta\le1/8\).
Retain the windows given by \(\|r/a-\varphi\|\le\theta\).
Assume the endpoints carry no mass, and write \(t>0\)
for their nonatomic mass. Suppose a finite list of positive
atoms on \(U,V\) omits at most \(t/16\) atomic mass and satisfies
\[
 \sum_{r_i\ \mathrm{selected}}\nu_{T_i}(\{r_i\})\,\|r_i/a\|
 \le \tfrac12\theta t.
 \tag{VM4}
\]
For all sufficiently small \(\theta\), uniformly in \(a,\varphi\)
and in the compact interval on the free letter, one has
\[
 t\le \max\{(8/5)\Delta_{UV},\, C_T\theta\}.
 \tag{VM5}
\]
The constant and smallness threshold can depend on a compact
positive range carrying nonatomic mass on the modulus letter.

### Retaining shifted windows

Retain the shifted windows on the selected letter and shrinking
neighborhoods of the selected atoms and zero on the pair.
Treat all selected atom neighborhoods as separate from the shifted
windows, removing any overlap from the shifted windows. The lost
nonatomic mass tends to zero. As in (VM3), the freeze costs at most
\(c_U+c_V-t+t/16\). Assign each positive shifted coordinate to its
mesh point \((i+\varphi)a\), and each selected atom to the closest
unshifted point \(ia\). Form a single signed rational root \(Z_1\)
using all the integer weights \(i\). Form \(Z_2\) using just the
shifted coordinates, with weight one and the sign of their letter.
The corresponding powers are
\(q_1=\lfloor a/L\rfloor\) and \(q_2=\lfloor \varphi a/L\rfloor\).
Both the common mesh and the phase are fixed before the small-budget
limit. Their height bounds are
\[
 H(Z_1)\le X^{(3/a+o(1))L},\qquad
 H(Z_2)\le X^{((17/16)t+o(1))L}.
\]
Indeed \(|i|\le r/a+1\) for either kind of mesh point, and the
shifted unlisted atomic mass is at most \(t/16\).
The wider rational factor \(V_0\) contains the zero windows and the
remaining powers after replacing the original exponents by
\(iq_1+q_2\) for the shifted windows and \(iq_1\) for the
selected atom windows. The rounding has vanishing numerical cost
because the lists of cells and selected atoms are fixed.
The wider exponent is at most
\[
 a\bigl(\theta(t+t/16)+\tfrac12\theta t\bigr)+o(1)
 <2a\theta t+o(1).
\]
The first term bounds the shifted window errors, and (VM4) bounds
the selected atom errors. At fixed frozen factors, \(V_0 Z_1^{q_1}
Z_2^{q_2}\) as an exact value determines the primitive triple.
Thus exact fibers cost at most \(X^\varepsilon\) at fixed truncation,
as in (VM3).

Choose a compact positive range \([b,B]\) containing nonatomic
mass on \(T\), with continuity endpoints. Let \(c_2\) be the
constant of (UW) in dimension two and fix
\(K>\max(1,18432/c_2)\). For sufficiently small \(\theta\),
choose a finite union of continuity intervals in this range with
weighted mass \(w\in[K\theta,2K\theta]\) for the measure
\(r\,d\nu_T(r)\). Its nonatomic part allows a Borel subset of
mass strictly between these endpoints containing no atoms. Take
a compact subset with mass above the lower endpoint, surround it
by an open set of mass below the upper endpoint, and take a finite
interval cover. The whole factors in the union give a full modulus
of numerical weight \(w\) and normalized cost at most \(w/b\).

If \(t>48K\theta/b\), set \(\kappa=24w/t<b\). In (UW)
take \(m=2\), \(C_1=96/(at)\), and \(C_2=48\). These caps exceed
the required narrow exponents, and the wider hypothesis holds
since
\[
 (2a\theta t/w) C_1 C_2=9216\theta/w\le9216/K<c_2 .
\]
The chart count costs \(6w/\kappa=t/4\) and the modulus costs
less than \(t/16\). Including the pair freeze gives
\(1/3\le c_U+c_V-5t/8\), proving (VM5).

### Bounds on the degree measures modulo the mesh

Let the weighted atomic discrepancy of the pair be
\[
 D(a)=\sum_{W=U,V}\ \sum_{r>0}\nu_W(\{r\})\,\|r/a\| .
\]
Write \(\mu_W^a\) for the image of the nonatomic part of \(\nu_W\)
under \(r\mapsto r/a\bmod 1\) on the circle. For every \(W\in\{U,V\}\)
and every arc \(I\) of the circle, the preceding estimate implies
\[
 \mu_W^a(I)\ \ll_T\
      \Delta_{UV}+|I|+\sqrt{D(a)}
 \qquad(0<a\le1).
 \tag{VM6}
\]
Here the implied constant is uniform in the mesh spacing and arc.
To prove it, first restrict to a compact positive continuity interval
on the free letter. Put \(\eta=\max\{|I|,\sqrt{2D(a)}\}\).
When \(\eta\) is small, enlarge the arc to a centered phase window
of radius \(\theta\in[\eta,2\eta]\) with continuity endpoints.
If its captured mass \(t\) is at most \(\eta\), the estimate holds.
Otherwise \(D(a)\le\theta t/2\), so any finite selection of atoms
omitting at most \(t/16\) satisfies (VM4). Apply (VM5). For larger
\(\eta\) the total mass bound suffices. A singleton has zero mass
under \(\mu_W^a\), since its preimage in the free letter is countable.
The constant in the compact-interval count is uniform, so
exhaustion gives (VM6).

At saturation \(\Delta_{UV}=0\), every weak limit of
\(\mu_W^{a_j}\) along mesh spacings with \(D(a_j)\to0\)
has a bounded density on the circle. For an open arc use weak
convergence and (VM6) to bound its limiting mass by a constant
times its length, then extend by regularity. Such mesh sequences
always exist: first keep finitely many atoms with small omitted
mass and use simultaneous approximation on their scales.
If all positive atoms on the pair lie in a common positive lattice,
choose its spacing or a positive integer subdivision with \(a\le1\).
Then \(D(a)=0\), so (VM6) bounds the folded measures themselves.
An interval of sufficiently small length on the degree axis maps
to an arc of length \(|I|/a\). Consequently both nonatomic degree
measures have bounded densities, also when the lattice has
infinitely many atomic scales and the numerical heights agree.

### An incommensurate ternary example

Let \(C\) be the standard middle-thirds Cantor random variable on \([0,1]\), whose independent base-three digits are zero and two with equal probability. Choose three positive numbers \(r_W\in(2,4)\) such that \(r_W/r_{W'}\notin\mathbb Q\) for distinct letters and such that along an unbounded sequence of integers \(k\), all \(\|3^k r_W\|\) tend to zero. For example, enumerate the letters by \(i=0,1,2\) and set
\[
 r_i=3+\sum_{n\ge1}3^{-(3n+i)!}.
\]
The nonzero fractional digits of the three scales occur at distinct factorial positions, with arbitrarily long common zero blocks. To verify incommensurability, if an integer linear combination of the scales with at least one nonzero coefficient were zero, multiply by \(3^{m!}\) for large \(m\). The preceding digits contribute an integer. The first following factorial digit with nonzero coefficient dominates the rest, since the factorial gaps tend to infinity. The resulting tail is nonzero and smaller than one in absolute value, a contradiction. This argument also proves the scales are individually irrational. Taking \(k=m!\) shows the required common recurrence.

Define the three formal root measures by putting mass \(1/6\) at \(r_W\) and mass \(1/6\) in the law of \(5-r_W+2C\). Each measure has total mass \(1/3\), first moment one, and nonatomic mass \(1/6\). These formal profiles satisfy the basic nonatomic constraints and the compact-support and height-gap conclusions. Each pair has two positive incommensurate atomic scales. The diffuse interval masses are \(O(|I|^\alpha)\) with \(\alpha=\log2/\log3>1/2\), by comparing \(I\) with the equal-length stage cells of the Cantor construction. Thus the measures also meet the [modular interval bound](ABC_QUANTITATIVE_DEGREE_WINDOWS.md#a-restriction-at-equality-of-the-nonatomic-total) with exponent \(1/2\) for two positive atoms.

At frequency \(3^k\), the fractional part of \(3^k C\) has the same law as \(C\), by dropping its first \(k\) ternary digits. The translated Cantor interval on each letter has offset \(5-r_W\), whose phase at \(3^k\) tends to zero along the common sequence. For every small \(\theta>0\), choose \(k\) in this sequence so large that all atomic phase errors are much smaller than \(\theta^{1+\alpha}\), where \(\alpha=\log 2/\log 3\). The Cantor mass of \(C\le\theta/8\) is at least a constant times \(\theta^\alpha\), by specifying its first \(\lceil\log_3(8/\theta)\rceil\) digits to be zero. Thus the phase mesh with base \(3^{-k}\) retains \(t\gg\theta^\alpha\) of the diffuse mass, with \(t/\theta\to\infty\). A fixed continuity interval \(J\) containing the whole diffuse support is sufficient. The positive atoms meet (VM4) by the choice of \(k\), and the remaining letter has positive nonatomic mass. The phase criterion excludes this family.

## Pair surplus at positive normalized exponent

Fix \(d_0>0\). Consider any limiting anatomic box with \(d\ge d_0\)
and any free pair of letters \(U,V\) with \(c_U+c_V\ge d\). Let
\(T\) be the complementary letter, put \(\Delta=c_U+c_V-d\),
and write \(e_T=h_T-\int r\,d\nu_T\). Then there are fixed
\(C,R_*>0\) depending only on \(d_0\) such that
\[
 e_T+\int_{[R_*,\infty)}r\,d\nu_T(r)\le C\Delta .\tag{VQ1}
\]
If \(U,V\) are both tall, then also
\[
 1-h_T\le C\Delta .
 \tag{VQ2}
\]
These conclusions apply to the pairs supplied by (LQ8) whenever
\(l_C<d\). At quarter equality with \(l_C=0\), if \(c_A=1/4\),
the complementary \(B\) is tall, has a bounded scaled-degree
support, and has no escaped first moment.

Use the construction of (VM3), with the sole initial change
that one of the free letters has nonatomic mass at least
\(d_0/2\) and supplies a compact patch of nonatomic mass at
least \(d_0/4\). For \(0<\theta\le1/8\), the same
simultaneous approximation and choice from
\(\lceil4/\theta\rceil\) meshes gives captured mass
\(t\ge d_0\theta/20\). Retain the two zero windows, finitely
many positive atoms missing at most \(t/16\), and the patch
on the chosen mesh. The frozen cost is at most
\(c_U+c_V-t+t/16\). The numerical height bounds are
unchanged:
\[
 H(Z)\le X^{(3/b+o(1))L},\qquad
 H(V_0)\le X^{2\theta bt+o(1)}.
\]
The total root mass on the free pair is now bounded by one
instead of \(2/3\). This still gives narrow height constant
\(2/b+1\le3/b\), and the selected atoms' numerical error is
at most another \(\theta bt/16\). Thus the same wider bound
holds. Exact values still recover the primitive triple.

For the tail, take a continuity cutoff \(R\) for \(\nu_T\).
The full high-degree factor has numerical weight
\(w=e_T+\int_{[R,\infty)}r\,d\nu_T\) and root cost
\(\mu\le w/R\). If \(w>0\), let
\(\eta_p=\min(1/8,c_p/256)\), with \(c_p\) from one-root (UW),
set \(\theta=\eta_pw\), and choose fullness parameter
\(kL\to\kappa=8w/t\). The choice is possible for
\(R>160/(d_0\eta_p)\). Use narrow constant \(32/(bt)\):
the wider relative cap times this constant is at most
\(64\eta_p\) with a strict margin to the (UW) threshold.
The chart cost is \(2w/\kappa=t/4\). If
\(R\ge320/(d_0\eta_p)\), also \(\mu\le t/16\).
The count gives \(d\le c_U+c_V-5t/8\), so
\(w\le32\Delta/(d_0\eta_p)\).
Let continuity cutoffs increase to
\(R_*=1+320/(d_0\eta_p)\), which includes an atom exactly
at \(R_*\).

For the real gap, choose \(0<\tau<1-h_T\) and the
complete tall-ratio strip of width \(X^{-\tau}\).
Set \(\eta_r=\min(1/8,c_r/256)\), with \(c_r\) from
one-root (RUW1), and use \(\theta=\eta_r\tau\).
Take its chart parameter \(s=t/8\) and narrow constant
\(4/b\). The height product is \(32/(bt)\), and the
wider cap fits the strict real-strip threshold.
The strip cost \(2s=t/4\) gives
\(d\le c_U+c_V-11t/16\). Therefore
\(1-h_T\le320\Delta/(11d_0\eta_r)\).

For the asserted quarter consequence, take free pair
\(A,C\) with \(\Delta=0\). If \(B\) were numerically
short, \(A,C\) would both be tall, contradicting the
zero-surplus real-gap estimate. The tail estimate gives
the compact support and moment assertion.

### A quarter-bound saving from near-maximal atomic mass and a high tail

Fix \(0\le\vartheta<1/3\) and
a positive \(w_0\), put \(d_*=(1+\vartheta)/4\), and take
\(d_0=1/8\) in the preceding estimates. Apply the
small-budget radical cap
\(\operatorname{rad}(C)\le c^{\vartheta\lambda}\).
Fix a positive integer \(K\) and \(0<\epsilon<w_0/C\).
Suppose \(A\) captures normalized root mass at least
\(d_*-\epsilon\) in at most \(K\) intervals of sufficiently small
fixed width in scaled degree. Suppose \(B\) has numerical weight
at least \(w_0\) at actual scaled degrees greater than \(R_*+1\).
For small fixed \(\lambda\) there is a fixed \(\delta>0\) giving a
count exponent \((d_*-\delta)\lambda+\varepsilon\).
A failing sequence with the width, \(\delta\) and \(\lambda\)
tending to zero would force quarter equality with
\(s=\vartheta\). The captured finite intervals force
atomic mass of \(A\) at least \(d_*-\epsilon\), so
\(c_A+c_C-d_*\le\epsilon\). The actual high-degree
moment on \(B\) survives the anatomic compression by
choosing the packed degree \(k\) with \(k\lambda\to\infty\)
and losing only the vanishing numerical remainder.
Thus the escaped-plus-high tail weight is at least \(w_0\),
contrary to the tail estimate.
