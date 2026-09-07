# Ordinary SIC in characteristic \(p\) and the quartic phase diagram

## 1. Integral model and statement

Let
\[
 R=\xi _1z_1+\xi _2z_2,\qquad
 Z=\xi _1z_2,\qquad
 W=2\xi _2z_1,\qquad
 T=\xi _1z_1-\xi _2z_2
 \tag{1.1}
\]
and clear the only denominator in the characteristic-zero witness by
\[
 \boxed{\widetilde F=2F
  =(R+Z)\left(2R^2W-(2R+Z)T^2\right).}
 \tag{1.2}
\]
Thus \(\widetilde F\) is defined over \(\mathbb Z\).  Write
\(\mathcal E_{2,p}\) for ordinary two-pair contraction after reduction to a
field \(k\) of characteristic \(p\).

> **Theorem 1.1 (ordinary-contraction phase diagram).** For every prime
> \(p\) and every \(m\geq1\),
> \[
>  \mathcal E_{2,p}(\widetilde F^m)=0,
>  \qquad
>  \mathcal E_{2,p}(Z\widetilde F^m)
>  =\overline{
>    2^m\frac{(4m+2)!\,m!}{(2m+1)!!}}.
>  \tag{1.3}
> \]
> The bar means reduction of the displayed **integer**, not inversion of
> the odd double factorial in \(k\).  The mixed contraction in (1.3) is
> nonzero exactly when
> \[
>  \boxed{p\ \text{is odd and}\ 4m+2<p.}
>  \tag{1.4}
> \]
> Consequently, for a fixed odd prime \(p\), its nonzero orders are exactly
> \[
>  \boxed{m=1,\ldots,\left\lfloor\frac{p-3}{4}\right\rfloor.}
>  \tag{1.5}
> \]

For odd \(p\), the original \(F=\widetilde F/2\) is also defined and
\[
 \mathcal E_{2,p}(F^m)=0,\qquad
 \mathcal E_{2,p}(ZF^m)
 =\overline{\frac{(4m+2)!\,m!}{(2m+1)!!}}.
 \tag{1.6}
\]
The nonvanishing criterion is again (1.4).  There is therefore no periodic
or infinite mixed-nonvanishing regime in a fixed characteristic: one gets
only the finite initial block (1.5).  In particular the block is empty for
\(p=2,3,5\).

> **Theorem 1.2 (geometric exceptional primes).**
>
> 1. The bilinear coordinate change in (1.1), and hence the displayed
>    rank-one quadric, is nondegenerate exactly when \(p\ne2\).
> 2. The \(5\) by \(5\) coefficient matrix of \(\widetilde F\) has full
>    rank exactly when \(p\ne2,3\).
> 3. Thus both the quadric chart and the polynomial tensor are
>    nondegenerate exactly in characteristics \(p\geq5\).
> 4. The reduction of \(\widetilde F\) is in the diagonal
>    \(\mathrm{SL}_2\)-nullcone at \(p=2\), but is outside the nullcone at
>    \(p=3\) and at every \(p\geq5\).

The last statement at \(p=3\) is an exact Gröbner/Hilbert--Mumford
calculation, not an inference from characteristic-zero semisimplicity.

## 2. The two degeneration loci

The four bilinears satisfy the integral identity
\[
 T^2=R^2-2ZW. \tag{2.1}
\]
In the ordered bilinear coordinates
\[
 (\xi _1z_1,\xi _1z_2,\xi _2z_1,\xi _2z_2),
\]
the linear change to \((R,Z,W,T)\) has determinant \(-4\).  Equivalently,
the quadric
\[
 T^2-R^2+2ZW=0 \tag{2.2}
\]
is smooth for odd \(p\), whereas at \(p=2\) it is the doubled hyperplane
\((T-R)^2=0\), and \(W\) itself has collapsed to zero.

The coefficient matrix of \(F\) has determinant \(48\).  Scaling the
whole bidegree-\((4,4)\) tensor by \(2\) scales its \(5\) by \(5\)
determinant by \(2^5\), so
\[
 \det C_{\widetilde F}=2^5\cdot48=1536=2^9\cdot3. \tag{2.3}
\]
Exact row reduction gives rank four at \(p=2,3\), and rank five at every
other prime.  This separates the quadric exception \(2\) from the
additional tensor-rank exception \(3\).

## 3. Ordinary moments over \(\mathbb Z\)

The characteristic-zero integral identity for \(F\) gives
\[
 \mathcal E_2(F^m)=0,\qquad
 \mathcal E_2(ZF^m)=\frac{(4m+2)!\,m!}{(2m+1)!!}. \tag{3.1}
\]
Multiplying \(F^m\) by \(2^m\) proves the integral identities
\[
\begin{aligned}
 \mathcal E_2(\widetilde F^m)&=0,\\
 \mathcal E_2(Z\widetilde F^m)
 &=2^m\frac{(4m+2)!\,m!}{(2m+1)!!}.
\end{aligned}
\tag{3.2}
\]
Formal differentiation has integral structure constants, so reduction of
(3.2) proves (1.3) in every characteristic.  No Gaussian argument and no
division modulo \(p\) is being used here.

It is useful to rewrite the unscaled mixed integer as
\[
\begin{aligned}
 N_m
 &=\frac{(4m+2)!\,m!}{(2m+1)!!}\\
 &=2^m\frac{(4m+2)!(m!)^2}{(2m+1)!}\\
 &=2^m\binom{4m+2}{2m+1}(2m+1)!(m!)^2.
\end{aligned}
\tag{3.3}
\]
The mixed moment of \(\widetilde F\) is \(2^mN_m\).  This factorization
also makes clear why treating \((2m+1)!!\) as an invertible denominator
would give the wrong modular question.

## 4. Legendre, Lucas, and Kummer criteria

Fix an odd prime \(p\).  Since powers of two are units, (3.3) gives
\[
 v_p(N_m)
 =v_p((4m+2)!)+2v_p(m!)-v_p((2m+1)!). \tag{4.1}
\]
Legendre's formula yields the floor-sum form
\[
 \boxed{
 v_p(N_m)=
 \sum_{j\geq1}\left(
 \left\lfloor\frac{4m+2}{p^j}\right\rfloor
+2\left\lfloor\frac m{p^j}\right\rfloor
-\left\lfloor\frac{2m+1}{p^j}\right\rfloor
 \right).}
 \tag{4.2}
\]
If \(s_p(n)\) is the sum of the base-\(p\) digits of \(n\), the equivalent
digit formula is
\[
 \boxed{
 v_p(N_m)=
 \frac{
 4m+1-s_p(4m+2)-2s_p(m)+s_p(2m+1)
 }{p-1}.}
 \tag{4.3}
\]

There is an especially short Lucas--Kummer proof of the zero criterion.
The last line of (3.3) is nonzero modulo \(p\) only if
\[
 m<p,\qquad 2m+1<p. \tag{4.4}
\]
Under (4.4), Kummer's theorem says that
\(\binom{4m+2}{2m+1}\) is nonzero precisely when adding
\((2m+1)+(2m+1)\) produces no base-\(p\) carry.  Because \(2m+1\) is then
a single base-\(p\) digit, this says
\[
 2(2m+1)<p, \tag{4.5}
\]
which is exactly (1.4).  Conversely, (4.5) makes every factorial and the
binomial coefficient in (3.3) nonzero.  This proves Theorem 1.1.

When \(4m+2<p\), an explicit nonzero residue is
\[
 N_m\equiv
 2^m\binom{4m+2}{2m+1}(2m+1)!(m!)^2\pmod p, \tag{4.6}
\]
and the residue for \(\widetilde F\) has \(2^m\) replaced by \(4^m\).
For \(p=2\), the factor \(2^m\) in the integral \(N_m\), and a fortiori
the extra factor for \(\widetilde F\), makes every positive-order mixed
moment zero.

### 4.1 The all-degree propagated phase diagram

The same criterion holds uniformly for the stronger radial-power family.
Let
\[
 d=4r+k,\qquad r\geq1,\quad k\geq0,\qquad
 \widetilde G_{r,k}=R^k\widetilde F^{\,r}. \tag{4.7}
\]
This is \(2^r\) times the characteristic-zero form \(R^kF^r\).  Its
integral moments are
\[
\begin{aligned}
 \mathcal E_2(\widetilde G_{r,k}^m)&=0,\\
 \mathcal E_2(Z\widetilde G_{r,k}^m)
 &=2^{rm}\frac{(dm+2)!(rm)!}{(2rm+1)!!}\\
 &=4^{rm}\frac{(dm+2)!((rm)!)^2}{(2rm+1)!}.
\end{aligned}
\tag{4.8}
\]

> **Theorem 4.1 (all propagated degrees).** For every odd prime \(p\),
> every \(d=4r+k\), and every \(m\geq1\),
> \[
>  \boxed{
>  \mathcal E_{2,p}(Z\widetilde G_{r,k}^m)\ne0
>  \quad\Longleftrightarrow\quad dm+2<p.}
>  \tag{4.9}
> \]
> Thus every propagated degree \(d\) has exactly the initial nonzero block
> \[
>  m=1,\ldots,\left\lfloor\frac{p-3}{d}\right\rfloor. \tag{4.10}
> \]

Put \(s=rm\).  Since \(dm+2\geq4s+2\), the last expression in (4.8) is
\[
 4^s(s!)^2\prod_{j=2s+2}^{dm+2}j. \tag{4.11}
\]
If \(dm+2<p\), every factor is a unit.  Conversely, if \(s\geq p\), the
factor \(s!\) vanishes.  If \(s<p\) and
\(2s+1<p\leq dm+2\), the product contains \(p\).  Finally, if
\(s<p\leq2s+1\), then
\[
 2s+2\leq2p\leq4s+2\leq dm+2,
\]
so the product contains \(2p\).  This proves (4.9) without any digit
casework.

Multiplication by the invariant \(R^k\) does not change torus weights, and
the least weight of a nonzero power is the corresponding multiple of the
least weight of the base.  Hence the nullcone classification propagates
as well: every \(\widetilde G_{r,k}\) is nullcone at \(p=2\), semistable
at \(p=3\), and semistable and tensor-nondegenerate at every good prime
where the relevant coefficient tensor retains its characteristic-zero
rank.  Only semistability, not a uniform tensor-rank claim in higher
degree, is asserted here.

### 4.2 Quantifier order and \(p\)-adic thickness

The phase diagram has a useful uniform reformulation.

> **Corollary 4.2 (arbitrarily long finite-prefix failures).** Fix a
> propagated degree \(d\geq4\) and a cutoff \(M\geq1\).  Reduction at every
> prime
> \[
>  p>dM+2 \tag{4.12}
> \]
> retains all mixed nonvanishings of orders \(1,\ldots,M\).  Conversely,
> characteristic \(p\) retains at most
> \(\lfloor(p-3)/d\rfloor\) consecutive orders.

Thus the two quantifier orders differ:
\[
\begin{aligned}
 &\text{for every }M\text{, all sufficiently large }p
   \text{ retain the first }M\text{ failures},\\
 &\text{for every fixed }p,\text{ only finitely many failures survive}.
\end{aligned}
\tag{4.13}
\]
This is the precise arithmetic shadow of the characteristic-zero
counterexample.  It supplies unbounded finite-prefix counterexamples over
finite fields, but never one asymptotic counterexample over a fixed finite
field.

The same formula records the thickness before complete reduction.  For an
odd prime \(p\), put \(s=rm\).  The mixed moment in (4.8) has valuation
\[
 \boxed{
 v_p=
 \frac{
 dm+1-s_p(dm+2)-2s_p(s)+s_p(2s+1)
 }{p-1}.}
 \tag{4.14}
\]
It is nonzero modulo \(p^a\) exactly when the integer in (4.14) is less
than \(a\).  For fixed \(p\) and \(a\), the valuation tends to infinity
linearly with \(m\), up to logarithmic digit-sum fluctuations, so every
fixed Artinian reduction again retains only finitely many orders.

Over \(\mathbb Z_p\) or \(\mathbb Q_p\), by contrast, the displayed mixed
moments are nonzero integers for every \(m\).  More generally, if \(k\)
is a perfect field of characteristic \(p\), the canonical map
\(\mathbb Z_p\to W(k)\) sends every nonzero mixed moment to a nonzero
element because \(W(k)\) is \(p\)-torsion-free.  Thus the ordinary
contraction witness survives unchanged over
\(\operatorname {Frac}W(k)\); no Hasse or divided-power modification is
needed on the Witt lift.  Its image in the truncated Witt ring
\(W_a(k)=W(k)/p^a\) survives precisely when (4.14) is less than \(a\).

This is the clean lift-level formulation: it is the passage to a fixed
residue characteristic, not \(p\)-adic size, that destroys the asymptotic
witness.  A statement internal to one residue field necessarily sees
only finite prefixes.

There is a sharp distinction between radial powers and the larger
propagated family at prime-power level.  Write
\[
 A_s=4^s\frac{(4s+2)!(s!)^2}{(2s+1)!}\qquad(s\geq1). \tag{4.14a}
\]
Direct cancellation gives the integral quotient
\[
 \boxed{\frac{A_{s+1}}{A_s}
 =16(4s+3)(4s+5)(s+1)^2.} \tag{4.14b}
\]
Consequently \(v_p(A_s)\) is nondecreasing for every prime, and for odd
\(p\)
\[
 v_p(A_{s+1})-v_p(A_s)
 =v_p(4s+3)+v_p(4s+5)+2v_p(s+1). \tag{4.14c}
\]
When \(k=0\), the moment of \(\widetilde F^r\) at order \(m\) is
\(A_{rm}\).  Its surviving orders modulo \(p^a\) therefore form an
initial interval.

This monotonicity does not extend even to the first non-radial lift.
For \(\widetilde G_{1,1}=R\widetilde F\), so \(d=5\), formula (4.8)
gives
\[
\begin{array}{c|c|c}
 m&v_{11}\!\left(
 \mathcal E_2(Z\widetilde G_{1,1}^m)\right)
 &\mathcal E_2(Z\widetilde G_{1,1}^m)\bmod 11^2\\ \hline
 4&2&0\\
 5&1&22.
\end{array}
\tag{4.14d}
\]
Thus a mixed moment can vanish modulo \(p^a\) and reappear at the next
order.  Formula (4.14) remains a complete decision procedure, and the
survival set is still finite, but outside the radial-power subfamily it
cannot be described by one initial-interval cutoff.  Any Witt-level
classification must retain the digit automaton or an equivalent
state-dependent recurrence.

That recurrence can be made completely local.  If \(M_m\) denotes the
integer in (4.8), then
\[
\begin{aligned}
v_p(M_{m+1})-v_p(M_m)
={}&r\,v_p(4)
 +\sum_{j=3}^{d+2}v_p(dm+j)
 +2\sum_{j=1}^{r}v_p(rm+j)\\
&-\sum_{j=2}^{2r+1}v_p(2rm+j).
\end{aligned}
\tag{4.14e}
\]
For \(d=4r\), this signed expression collapses to the sum of the
nonnegative increments (4.14c) from \(s=rm\) through
\(s=r(m+1)-1\).  For \(d>4r\), the uncancelled final block is the source
of negative jumps.  Formula (4.14e), rather than the factorial-sized
integer \(M_m\), is the natural input for a base-\(p\) automaton.

#### 4.2.1 Recurrence curvature and the integral singularity ledger

The local recurrence in (4.14e) has a useful operator formulation.  Cancel
the consecutive quotient
\[
 \frac{M_{m+1}}{M_m}
 =
 4^r
 \frac{
  \prod_{j=3}^{d+2}(dm+j)
  \left(\prod_{j=1}^{r}(rm+j)\right)^2
 }{
  \prod_{j=2}^{2r+1}(2rm+j)
 }
 \tag{4.14f}
\]
over \(\mathbb Q[m]\), and write the result as \(B_{d,r}(m)/A_{d,r}(m)\)
with coprime primitive \(A_{d,r},B_{d,r}\in\mathbb Z[m]\).  Then
\[
 A_{d,r}(m)M_{m+1}=B_{d,r}(m)M_m                     \tag{4.14g}
\]
is the minimal order-one recurrence, and (4.14e) becomes the exact local
transition
\[
 \boxed{
 v_p(M_{m+1})-v_p(M_m)
 =v_p(B_{d,r}(m))-v_p(A_{d,r}(m)).}                    \tag{4.14h}
\]
Thus the primitive factorization, not the factorial-sized moment, is the
canonical input for the digit automaton.

The recurrence also has a characteristic-\(p\) curvature.  If
\(\tau(m)=m+1\), its \(p\)-step multiplier is
\[
 \Psi_{p,d,r}(m)
 =\prod_{i=0}^{p-1}
 \frac{B_{d,r}(m+i)}{A_{d,r}(m+i)}.                    \tag{4.14i}
\]
For \(p>d+2\), the shift norm of every affine factor \(am+b\) is
\(a(m^p-m)\), so
\[
 \boxed{\Psi_{p,d,r}=d^d(m^p-m)^d.}                    \tag{4.14j}
\]
This exact curvature is a global consistency invariant, not a replacement
for (4.14h): taking the shift norm cancels the separate numerator and
denominator factors that control prime-power re-entry.  For example,
\((d,r)=(8,1)\) and \((8,2)\) have the same curvature
\(8^8(m^p-m)^8\).  The first primitive recurrence has forward factor
\(2m+3\) and admits negative valuation jumps; the second has forward
factor \(1\) and is valuation-monotone.

The normalized angular beta period supplies a parallel warning.  Every
first radial lift \(R^k\widetilde F\) has the same Picard--Fuchs operator
up to a dilation of its generating variable, although its factorial
valuation phase depends on \(d=4+k\).  The exact recurrence derivation,
good-prime curvature proof, differential \(p\)-curvature experiment through
\(101\), and degree-eight control are in
[`TWO_PAIR_SIC_FROBENIUS_CURVATURE_BRIDGE.md`](TWO_PAIR_SIC_FROBENIUS_CURVATURE_BRIDGE.md).

### 4.3 Non-power Hopf profiles and numerator-prime holes

The non-power family \(\Phi_h\) has mixed moment
\[
 M_{h,m}=(4hm+2)!\,C_{h,m},\qquad
 C_{h,m}=\int_0^1
 (1-v^2)^m(1+v^2)^{(h-1)m}\,dv. \tag{4.15}
\]
Put \(s=hm\).  Expanding the integrand shows that the reduced denominator
of \(C_{h,m}\) divides
\[
 \operatorname {lcm}(1,3,5,\ldots,2s+1). \tag{4.16}
\]

> **Theorem 4.3 (universal cutoff for every Hopf profile).** For every
> prime \(p\),
> \[
>  p\leq4hm+2\quad\Longrightarrow\quad
>  M_{h,m}\equiv0\pmod p. \tag{4.17}
> \]
> If \(p>4hm+2\), then
> \[
>  M_{h,m}\not\equiv0\pmod p
>  \quad\Longleftrightarrow\quad
>  C_{h,m}\not\equiv0\pmod p, \tag{4.18}
> \]
> where every denominator in the integral expansion is then a unit.

For an odd \(p\), let \(e\) be maximal with \(p^e\leq2s+1\), taking
\(e=0\) if no such power exists.  Formula (4.16) bounds the denominator
valuation by \(e\).  If \(e=0\), the factorial in (4.15) contains \(p\).
If \(e\geq1\), then \(2p^e\leq4s+2\), and
\[
 v_p((4s+2)!)
 \geq\sum_{j=1}^e
 \left\lfloor\frac{2p^e}{p^j}\right\rfloor
 \geq2e>e. \tag{4.19}
\]
Thus at least one factor of \(p\) remains after every possible denominator
cancellation.  Characteristic two is immediate because (4.16) is odd.
This proves (4.17); (4.18) follows because both the factorial and all
integral denominators are units above the cutoff.

Height two admits a complete answer.  Here
\[
 C_{2,m}=\int_0^1(1-v^4)^m\,dv
 =\frac{4^m m!}{\prod_{j=0}^m(4j+1)}. \tag{4.20}
\]
The last equality follows formally from \(C_{2,0}=1\) and the
integration-by-parts recurrence
\[
 (4m+1)C_{2,m}=4mC_{2,m-1}. \tag{4.21}
\]
All numerator primes in (4.20) are at most \(m\).  Combining this with
Theorem 4.3 gives:

> **Corollary 4.4 (complete \(\Phi_2\) phase diagram).**
> \[
>  \boxed{
>  \mathcal E_{2,p}(Z\Phi_2^m)\ne0
>  \quad\Longleftrightarrow\quad 8m+2<p.}
>  \tag{4.22}
> \]

For general height, (4.18) can have additional holes above the factorial
cutoff.  The first small exact examples found by the checker are
\[
\begin{array}{c|c|c|c}
 h&m&p&\text{reason}\\ \hline
 6&1&47&C_{6,1}=3008/1287,\quad3008=2^6\cdot47,\\
 4&5&89&C_{4,5}=6606540505088/1261867452417,\\
 &&&6606540505088=2^{23}\cdot89\cdot8849.
\end{array}
\tag{4.23}
\]
Both primes are strictly larger than \(4hm+2\).  Hence the clean cutoff
is sufficient for \(h=1,2\) but only necessary for arbitrary \(h\).
Classifying the numerator-prime holes is a genuine terminating
hypergeometric Lucas problem.

## 5. The Image kernel does not fail; differentiation becomes nilpotent

Let
\[
 {\cal M}_{r,p}
 =\sum_{i=1}^r(\partial_{z_i}-\xi_i)
 k[\xi_1,\ldots,\xi_r,z_1,\ldots,z_r]. \tag{5.1}
\]
The Image-kernel identity itself remains valid:
\[
 \boxed{{\cal M}_{r,p}=\ker\mathcal E_{r,p}} \tag{5.2}
\]
over every field, in every characteristic.  Indeed, in one pair,
\[
 \begin{aligned}
 &(\partial_z-\xi)
 \sum_{j=0}^{a-1}\xi^{a-1-j}\partial_z^j(z^b)\\
 &\hspace{35mm}=\partial_z^a(z^b)-\xi^az^b.
 \end{aligned}
 \tag{5.3}
\]
This is a division-free telescoping identity.  Applying it successively
in all variables reduces every monomial to its contraction, so the kernel
is the displayed sum of images.  The reverse containment follows directly
from
\(\mathcal E_{r,p}((\partial_{z_i}-\xi_i)H)=0\).

What changes is
\[
 \partial_{z_i}^{\,p}=0 \tag{5.4}
\]
on the polynomial ring.  Let
\[
 \delta_\xi(f)=
 \min\{|\alpha|:[\xi^\alpha z^\beta]f\ne0\}. \tag{5.5}
\]
If \(\delta_\xi(f)>0\), then every monomial of \(g f^m\) has total
\(\xi\)-degree at least \(m\delta_\xi(f)\).  If
\[
 m\delta_\xi(f)>r(p-1), \tag{5.6}
\]
some \(\xi_i\)-exponent is at least \(p\), so that monomial has zero
ordinary contraction.  Hence:

> **Theorem 5.1 (automatic ordinary SIC away from the dual-degree-zero
> face).** For every fixed \(g\) and every polynomial \(f\) with
> \(\delta_\xi(f)>0\),
> \[
>  \mathcal E_{r,p}(g f^m)=0
>  \quad\text{and}\quad
>  gf^m\in{\cal M}_{r,p}
>  \qquad\left(m>\frac{r(p-1)}{\delta_\xi(f)}\right).
>  \tag{5.7}
> \]
> This conclusion does not require any pure-moment premise.

For \(\widetilde F\), \(r=2\) and \(\delta_\xi(\widetilde F)=4\), so every
multiplier has
eventual Image membership once \(m>(p-1)/2\).  Thus the reductions of
\(\widetilde F\) are not positive-characteristic SIC counterexamples.
They are finite-prefix failures of mixed membership when (1.4) holds.
The ordinary derivative formulation makes balanced positive-degree SIC
automatic rather than producing a modular analogue of the
characteristic-zero phenomenon.

### 5.1 Frobenius closes the nonhomogeneous frontier

The apparent dual-degree-zero frontier collapses at the \(p\)-th pure
moment.  Put
\[
 f_0(z)=f(0,z). \tag{5.8}
\]
The Freshman's dream and (5.4) give the exact identity
\[
 \boxed{\mathcal E_{r,p}(f^p)=f_0(z)^p.} \tag{5.9}
\]
Indeed, write
\[
 f=\sum_{\alpha,\beta}c_{\alpha,\beta}\xi^\alpha z^\beta.
\]
Then
\[
 f^p=\sum_{\alpha,\beta}
 c_{\alpha,\beta}^p\xi^{p\alpha}z^{p\beta}.
\]
If \(\alpha\ne0\), some contraction order \(p\alpha_i\) is at least
\(p\), so that term is killed by (5.4).  The terms with \(\alpha=0\)
sum to \(f_0^p\).

> **Theorem 5.2 (full ordinary SIC in characteristic \(p\)).** Let \(k\)
> be any field of characteristic \(p>0\), and let \(r\geq1\).  If
> \[
>  \mathcal E_{r,p}(f^p)=0, \tag{5.10}
> \]
> then, for every fixed multiplier \(g\),
> \[
>  \boxed{
>  \mathcal E_{r,p}(gf^m)=0,\qquad
>  gf^m\in{\cal M}_{r,p}
>  \quad\text{for every }m\geq p.}
>  \tag{5.11}
> \]
> In particular, the ordinary Special Image Conjecture holds in every
> pair dimension over every positive-characteristic field.  Its full
> all-moment premise can be replaced by the single \(p\)-th moment.

By (5.9), condition (5.10) forces \(f_0=0\), since the polynomial ring is
reduced.  Every monomial of \(f\) then has positive dual degree.  Hence
every monomial of
\[
 f^p=\sum_{\alpha\ne0,\beta}
 c_{\alpha,\beta}^p\xi^{p\alpha}z^{p\beta} \tag{5.12}
\]
has some \(\xi_i\)-exponent at least \(p\).  For \(m\geq p\), factor
\(gf^m=g f^p f^{m-p}\).  Multiplication cannot decrease that exponent, so
ordinary contraction kills every monomial.  The Image-kernel identity
(5.2) proves membership.  If the actual minimum dual degree is large,
the independent pigeonhole cutoff (5.7) can still be earlier than \(p\).

The universal cutoff \(m\geq p\) is sharp, even under the full pure-moment
premise.  Take
\[
 f=\xi _1,\qquad g=z_1^{p-1}. \tag{5.13}
\]
Then \(\mathcal E_{r,p}(f^m)=0\) for every \(m\geq1\), while
\[
 \mathcal E_{r,p}(g f^{p-1})=(p-1)!\ne0. \tag{5.14}
\]
All mixed contractions vanish from order \(p\) onward.

This resolves the ordinary positive-characteristic formulation completely:
the reduced SIC2C4 moments illustrate finite-prefix nonmembership, but no
ordinary polynomial in any pair dimension can give an asymptotic
positive-characteristic SIC counterexample.

## 6. A naive divided-power replacement

The most immediate attempt to remove (5.4) is to replace ordinary
derivatives by Hasse derivatives:
\[
 {\cal H}_2(\xi^\alpha z^\beta)
 =\binom{\beta_1}{\alpha_1}
  \binom{\beta_2}{\alpha_2}z^{\beta-\alpha}. \tag{6.1}
\]
For a balanced polynomial, its scalar Hasse contraction is simply the sum
of its diagonal coefficient-matrix entries.  Equivalently,
\[
 {\cal H}_2(H)
 =\operatorname {CT}_u H(1,u,1,u^{-1}). \tag{6.2}
\]
On this chart,
\[
 R=2,\qquad Z=u^{-1},\qquad W=2u,\qquad T=0,
 \qquad
 \widetilde F=16(1+2u). \tag{6.3}
\]
Constant-term extraction therefore gives the all-order identities
\[
 \boxed{
 {\cal H}_2(\widetilde F^m)=16^m,\qquad
 {\cal H}_2(Z\widetilde F^m)=2m\,16^m.}
 \tag{6.4}
\]

Thus the naive divided-power/Hasse contraction is cleaner algebraically
but does **not** preserve the characteristic-zero cancellation: for every
odd \(p\), the pure moments in (6.4) never vanish.  At \(p=2\) both
displayed sequences vanish only because the integral seed itself has
undergone the exceptional degeneration described below.  Moreover,
\({\cal H}_2\) is not the kernel functional for the ordinary operators
\(\partial_{z_i}-\xi_i\).

A nontrivial positive-characteristic Image conjecture would therefore
have to specify a divided-power algebra, its multiplication, and compatible
Image operators simultaneously.  Merely replacing ordinary derivatives
by Hasse derivatives is not a transfer of SIC(2), and the present seed does
not satisfy its pure premise in odd characteristic.

### 6.1 A complete Hasse-compatible Image system

There is nevertheless a canonical kernel theorem if one permits all
divided differential orders.  For one pair define
\[
\begin{aligned}
 \partial_z^{[a]}(z^b)&=\binom ba z^{b-a},\\
 X_\xi^{[a]}(\xi^c)&=\binom{c+a}{a}\xi^{c+a},\\
 D_a&=\partial_z^{[a]}-X_\xi^{[a]}\qquad(a\geq1).
\end{aligned}
\tag{6.5}
\]
The coordinate and dual operators act on their own variables.  In several
pairs use multi-indices, or equivalently products of the one-pair
operators.

> **Theorem 6.1 (complete Hasse Image-kernel identity).** Over every
> field,
> \[
>  \boxed{
>  \ker{\cal H}_r
>  =\sum_{\alpha\in\mathbb N^r\setminus\{0\}}
>    \operatorname {Im}
>    \left(\partial_z^{[\alpha]}-X_\xi^{[\alpha]}\right).}
>  \tag{6.6}
> \]

The termwise intertwining is the elementary trinomial identity
\[
 \binom b a\binom{b-a}c
 =\binom{c+a}a\binom b{c+a}. \tag{6.7}
\]
It proves
\({\cal H}_r(\partial_z^{[\alpha]}G)
={\cal H}_r(X_\xi^{[\alpha]}G)\).
Conversely, on a monomial with no dual factor in the input,
\[
 D_\alpha(z^\beta)
 =\binom\beta\alpha z^{\beta-\alpha}
  -\xi^\alpha z^\beta
 ={\cal H}_r(\xi^\alpha z^\beta)-\xi^\alpha z^\beta.
 \tag{6.8}
\]
Thus every polynomial differs from its Hasse contraction by an element of
the right-hand side of (6.6), proving the reverse inclusion without
division.

In characteristic \(p\), the apparently enormous operator family has a
clean \(p\)-typical reduction.

> **Corollary 6.2 (\(p\)-power generators suffice).**
> \[
>  \boxed{
>  \ker{\cal H}_r
>  =\sum_{i=1}^r\sum_{j\geq0}
>   \operatorname {Im}
>   \left(\partial_{z_i}^{[p^j]}-X_{\xi_i}^{[p^j]}\right).}
>  \tag{6.9}
> \]

Indeed, if \(\binom{a+b}{a}\) is a unit, then
\[
 \binom{a+b}{a}D_{a+b}
 =\partial^{[a]}D_b+D_aX^{[b]}. \tag{6.10}
\]
Build any integer from its base-\(p\) digits by adding copies of \(p^j\)
without carries.  Lucas's theorem makes every binomial coefficient used
in (6.10) a unit.  Products in different coordinates are separated by the
usual telescoping identity for a difference of two products.

This gives a precise positive-characteristic Image problem on the
ordinary polynomial vector space:
\[
 {\cal M}^{\mathrm H}_{r,p}:=\ker{\cal H}_r. \tag{6.11}
\]
It is not the original SIC problem: it uses countably many \(p\)-typical
operators rather than \(r\) first-order operators.  For the present seed,
(6.4) excludes the pure premise at every odd prime.  At \(p=2\),
\(\widetilde F=ZR^2(R+Z)\) is already supported in strictly one-sided
torus weights.  A fixed multiplier cannot offset the weight gap of its
large powers, so
\[
 {\cal H}_2(g\widetilde F^m)=0\qquad(m\gg0). \tag{6.12}
\]
Thus the exceptional seed does not fail this Hasse-SIC analogue either.

### 6.2 The Hasse Image is not Mathieu

The global Hasse analogue nevertheless fails in the smallest possible
pair dimension.

> **Theorem 6.3 (one-pair Hasse-SIC counterexample).** In every
> characteristic \(p>0\), put
> \[
>  f=\xi z^p,\qquad g=z. \tag{6.13}
> \]
> Then
> \[
>  {\cal H}_1(f^m)=0\qquad(m\geq1), \tag{6.14}
> \]
> but, for every \(e\geq1\) and
> \[
>  m_e=1+p+\cdots+p^{e-1}=\frac{p^e-1}{p-1}, \tag{6.15}
> \]
> one has
> \[
>  {\cal H}_1(gf^{m_e})
>  =z^{(p-1)m_e+1}\ne0. \tag{6.16}
> \]
> Hence \({\cal M}^{\mathrm H}_{1,p}\), and therefore the corresponding
> \(p\)-typical Image in every positive pair dimension, is not a Mathieu
> subspace.

Indeed,
\[
 {\cal H}_1(f^m)
 =\binom{pm}{m}z^{(p-1)m}.
\]
If \(p^t\) is the exact power of \(p\) dividing \(m\), the \(t\)-th
base-\(p\) digit of \(m\) is nonzero while the same digit of \(pm\) is
zero.  Lucas's theorem gives \(\binom{pm}{m}=0\).  On the other hand,
\(m_e\) has \(e\) consecutive base-\(p\) digits equal to one, and
\(pm_e+1\) has \(e+1\) such digits.  Lucas now gives
\[
 \binom{pm_e+1}{m_e}=1,
\]
proving (6.16).  This failure is intrinsic to the Hasse replacement and
is unrelated to the semistable quartic seed.

### 6.3 Why a fully divided-power algebra becomes vacuous

One might instead change the multiplication as well and work in the
divided-power algebra \(\Gamma(V)\), which restores the natural
equivariant pairing.  This removes the mismatch between Hasse contraction
and ordinary symmetric-power multiplication, but introduces a decisive
Frobenius obstruction:
\[
 \boxed{x^p=0\quad\text{for every finite positive-degree }
 x\in\Gamma(V)\ \text{over characteristic }p.} \tag{6.17}
\]
For a divided-power monomial \(v^{[n]}\), the coefficient of
\((v^{[n]})^p\) is
\[
 \frac{(pn)!}{(n!)^p},
\]
whose \(p\)-adic valuation is the positive base-\(p\) digit sum \(s_p(n)\).
The multivariable statement follows from a nonzero coordinate, and the
Freshman's dream handles finite sums.  Consequently every positive-degree
candidate has \(f^m=0\) for \(m\geq p\).

There is therefore a three-way obstruction:

1. ordinary differentiation preserves the original Image operators and
   satisfies the sharp full theorem \(\operatorname{SIC}(r)\) for every
   \(r\);
2. Hasse contraction on the ordinary polynomial algebra stays
   nonnilpotent, its exact Image description requires the infinite
   \(p\)-typical system (6.9), and that Image already fails the Mathieu
   property in one pair by Theorem 6.3;
3. the fully equivariant divided-power algebra restores the natural
   pairing, but Frobenius makes every positive-degree element
   \(p\)-nilpotent and the eventual-power question vacuous.

Any genuinely nontrivial modular replacement must evade all three
mechanisms, for example by using a \(p\)-adic or Witt lift before reduction,
or by asking a uniform-in-\(p\) finite-prefix question rather than an
eventual-power question inside one fixed characteristic.

## 7. Exceptional geometry and enlarged moment-zero loci

At \(p=2\),
\[
 W=0,\qquad T=R,\qquad
 \widetilde F=ZR^2(R+Z). \tag{7.1}
\]
Choose the diagonal one-parameter subgroup orientation for which \(Z\)
has positive weight.  Since \(R\) has weight zero, both terms
\(R^3Z+R^2Z^2\) have positive weight.  Hence the displayed reduction lies
in the nullcone.

At \(p=3\),
\[
 \widetilde F=(R+Z)
 \left(R^3-R^2W-R^2Z+RZW-Z^2W\right). \tag{7.2}
\]
The characteristic-zero separator
\[
 I_2(\widetilde F)
 =\operatorname {tr}((DC_{\widetilde F})^2)=4608
 =2^9\cdot3^2 \tag{7.3}
\]
vanishes modulo \(3\), and the usual
\(\operatorname{End}(\operatorname{Sym}^4)\) self-duality matrix \(D\)
is itself singular there.  Semistability therefore cannot be imported
from the characteristic-zero representation decomposition.

For an exact modular Hilbert--Mumford test, write
\[
 g=\begin{pmatrix}a&b\\c&d\end{pmatrix},\qquad ad-bc=1, \tag{7.4}
\]
apply \(z\mapsto gz\) and \(\xi\mapsto g^{-T}\xi\), and set to zero every
coefficient of \(g\widetilde F\) whose standard torus weight is
nonpositive.  In the basis
\(\xi _1^i\xi _2^{4-i}z_1^jz_2^{4-j}\), that means all positions
\(j\leq i\).  Over \(\mathbb F_3[a,b,c,d]\), these fifteen coefficient
equations together with \(ad-bc-1\) have reduced Gröbner basis
\[
 \{1\}. \tag{7.5}
\]
Thus no conjugate over the algebraic closure is strictly one-sided, and
\(\widetilde F\bmod3\) is outside the nullcone.  For every \(p\geq5\),
(7.3) is nonzero, so it directly separates \(\widetilde F\) from the
nullcone.

Coefficient-tensor rank has a different higher-degree phase diagram.
For the first propagated slice \(R^k\widetilde F\), direct exact
determinants give
\[
\begin{array}{c|c|c|c}
k&d&\det C_{R^k\widetilde F}&
 \text{exceptional }p\text{ (rank)}\\ \hline
0&4&2^9\cdot3&2(4),3(4)\\
1&5&-2^{10}\cdot17&2(4),17(5)\\
2&6&2^{14}\cdot23&2(4),23(6)\\
3&7&-2^{12}\cdot3^2\cdot5\cdot83&2(4),3(7),5(7),83(7)\\
4&8&2^{16}\cdot5^2\cdot13\cdot17&2(8),5(7),13(8),17(8).
\end{array}
\tag{7.6}
\]
There is an exact linear-size determinant algorithm for this entire
\(r=1\) slice.  Index rows and columns from zero, and put
\(C^{(k)}=C_{R^k\widetilde F}\).  Its only nonzero diagonals satisfy
\[
\begin{array}{c|c}
\delta&\displaystyle
 C^{(k)}_{i,i+\delta}=[x^i]q_{\delta,k}(x)\\ \hline
1&4(1+x)^{k+3}\\
0&-2(1+x)^{k+2}(x^2-4x+1)\\
-1&-3x(x-1)^2(1+x)^{k+1}\\
-2&-x^2(x-1)^2(1+x)^k.
\end{array}
\tag{7.7}
\]
Indeed, multiplication by
\[
R^k=\sum_{a=0}^k\binom ka
(\xi _1z_1)^a(\xi _2z_2)^{k-a}
\]
convolves every diagonal of \(C_{\widetilde F}\) by the coefficients of
\((1+x)^k\); factoring the four base diagonal polynomials gives (7.7).

Write \(a_i=C^{(k)}_{i,i}\), \(b_i=C^{(k)}_{i,i+1}\),
\(c_i=C^{(k)}_{i,i-1}\), and \(e_i=C^{(k)}_{i,i-2}\).  If
\(\Delta_n\) is the determinant of the leading \(n\)-by-\(n\) block,
with \(\Delta_0=1\), lower-Hessenberg expansion gives
\[
\boxed{
\Delta_n=a_{n-1}\Delta_{n-1}
-c_{n-1}b_{n-2}\Delta_{n-2}
+e_{n-1}b_{n-3}b_{n-2}\Delta_{n-3},}
\tag{7.8}
\]
where terms with negative indices are omitted.  Therefore
\[
\det C_{R^k\widetilde F}=\Delta_{k+5}. \tag{7.9}
\]
This reduces every exact determinant and modular-rank gate in the first
slice to \(O(k)\) binomial arithmetic.  It does not yet factor
\(\Delta_{k+5}\), whose large primitive prime factors explain why a
simple fixed exceptional set is unavailable.

Characteristic two admits a closed rank formula for every \(k\).  Put
\(n=k+2\), and let \(s_2(q)\) denote binary digit sum.  From (7.1),
\[
R^k\widetilde F=ZR^n(R+Z)
=ZR^{n+1}+Z^2R^n. \tag{7.10}
\]
After deleting the forced zero first row and last column, its coefficient
matrix is lower bidiagonal of size \(n+2\), with diagonal
\(\binom{n+1}{a}\) and subdiagonal \(\binom n a\).  Lucas's theorem then
gives
\[
\boxed{
\operatorname {rank}_{\mathbb F_2}C_{R^k\widetilde F}
=2^{\,1+s_2(\lfloor(k+2)/2\rfloor)}.}
\tag{7.11}
\]
For completeness, if \(n=2q\), the nonzero diagonal entries occur in
pairs indexed by the binary submasks of \(q\), giving rank
\(2^{1+s_2(q)}\).  If \(n=2q+1\), the two subdiagonal entries indexed by
each submask of \(q\) give a permutation minor of that size.  Pascal's
identity shows every remaining nonzero diagonal column is either already
among those columns or duplicates the preceding odd column, proving the
matching upper bound.

Thus there is no degree-independent exceptional set for full coefficient
rank: primes \(17,23,83\) already occur in degrees five through seven.
This rank loss does not imply nullcone collapse.  Since \(R\) is
invariant and has torus weight zero, multiplication by \(R^k\) preserves
the nonzero weight components of every conjugate.  Hence
\(R^k\widetilde F\) remains semistable at every odd prime, including all
the new exceptional rank primes in (7.6); only characteristic two
inherits the displayed seed's nullcone collapse.

Small characteristic does create much larger failures of the
**moment--nullcone equality**, but not SIC failures.  On
\[
 V_4=\operatorname{Sym}^4(U^*)\otimes\operatorname{Sym}^4(U),
\]
ordinary contraction of every \(H^m\) is automatically zero as soon as
\(4m>2(p-1)\).  Therefore:

- at \(p=2\), every point of the \(25\)-dimensional \(V_4\) has all pure
  moments zero;
- at \(p=3\), every moment with \(m\geq2\) vanishes universally, while
  the first moment is the single central-coefficient condition
  \(c_{22}=0\), so the all-moment-zero locus is a \(24\)-dimensional
  hyperplane;
- in both cases the invariant point \(R^4\) is a semistable member of the
  all-moment-zero locus.  At \(p=3\), the displayed
  \(\widetilde F\) is another semistable member; at \(p=2\), the displayed
  seed itself has collapsed into the nullcone.

Together with the \(p\geq5\) separator (7.3), this gives a uniform
conclusion.

> **Theorem 7.1 (modular moment--nullcone failure).** The equality between
> the common zero locus of all ordinary contraction moments on \(V_4\)
> and the diagonal \(\mathrm{SL}_2\)-nullcone fails in every positive
> characteristic.  For \(p\geq3\), \(\widetilde F\) itself is a semistable
> all-moment-zero point.  For \(p=2\), the displayed seed is nullcone but
> the invariant point \(R^4\) is semistable and all of \(V_4\) is
> moment-zero.

These enlarged loci reflect factorial nilpotence in ordinary contraction.
They should not be reported as positive-characteristic Image-Mathieu
counterexamples, because Theorem 5.1 forces eventual mixed Image
membership.

## 8. Reproduction and scope

Run
```bash
.venv/bin/python scripts/verify_two_pair_sic_characteristic_p.py
```
from the repository root.  The checker:

- reconstructs \(\widetilde F\), (2.1), (2.3), and the modular ranks;
- replays the ordinary and Hasse moments through order eight;
- checks the valuation, floor-sum, digit-sum, and carry criteria for all
  primes through \(101\);
- checks every radial-power family \(R^k\widetilde F^r\) through degree
  twenty and order \(2p\) against the cutoff \(dm+2<p\);
- proves the integral consecutive-term quotient (4.14b), checks radial
  valuation monotonicity and recurrence (4.14e), and verifies the
  non-radial re-entry (4.14d);
- verifies the full Frobenius collapse
  \(\mathcal E_{r,p}(f^p)=f(0,z)^p\) monomialwise and the sharp factorial
  cutoff example;
- checks the \(\Phi_2\) closed form through order eighty, the universal
  profile cutoff through height eight, and the two numerator-prime holes
  in (4.23);
- audits the Hasse/divided-raising intertwining and the \(p\)-typical
  no-carry generator reduction;
- verifies the explicit characteristic-two destabilization;
- recomputes the characteristic-three Gröbner basis (7.5);
- computes the exact coefficient determinants (7.6) and modular ranks
  through degree eight, and audits the all-\(k\) diagonal symbols and
  continuant (7.7)--(7.9) through \(k=20\);
- checks the all-\(k\) characteristic-two rank formula (7.11) through
  \(k=128\); and
- writes
  `artifacts/generated-results/two_pair_sic_characteristic_p.json`.

The finite cutoffs are regression checks.  The all-order claims are proved
by the integral moment identity, formulas (3.3)--(4.5), the telescoping
identity (5.3), and the constant-term calculation (6.2)--(6.4).

## 9. Remaining continuations

The phase diagram, Theorem 5.2, and Theorem 6.3 settle both the ordinary
and Hasse Image-Mathieu questions over a fixed residue field.  Retaining
finer arithmetic and geometric structure leaves four distinct problems.

1. **Higher non-power Hopf profiles.** Theorem 4.3 gives the universal
   necessary cutoff and Corollary 4.4 completely classifies \(\Phi_2\).
   Starting at higher profiles, numerator-prime holes such as
   \((h,m,p)=(6,1,47)\) and \((4,5,89)\) occur above the cutoff.  Classify
   these holes using a Lucas theory for the terminating hypergeometric
   integral.
2. **Exact modular moment-zero schemes.** At \(p=2\) the scheme is all of
   \(V_4\), and at \(p=3\) its reduced set is the hyperplane \(c_{22}=0\).
   Scheme structure, quotient geometry, and semistable components for
   \(p\geq5\) remain unclassified.
3. **Prime-power and Witt formulations.** Formula (4.14) gives an exact
   decision procedure modulo \(p^a\).  Formula (4.14b) makes radial-power
   survival an initial interval, but (4.14d) shows that non-radial
   survival can disappear and re-enter.  The ordinary witness already
   survives over \(\operatorname {Frac}W(k)\), so the remaining problem
   is not to modify contraction.  For the radial propagations, (4.14h)
   now supplies the exact factorized digit automaton; the next task is to
   minimize its carry states and extract maximal surviving orders without
   returning to factorial expansion.  For non-power profiles and future
   higher-order telescopers, retain an integral companion lattice and its
   stepwise singular factors before taking \(p\)-curvature.  It remains to
   decide whether a useful ring-theoretic Mathieu property commutes with
   the inverse system \(W_a(k)\).
4. **Higher-degree exceptional tensor rank.** Semistability propagates to
   \(R^k\widetilde F^r\).  Equations (7.7)--(7.9) give an exact
   continuant for every \(k\) when \(r=1\), with the initial
   factorizations in (7.6), while (7.11) completely classifies their
   characteristic-two ranks.  Factor or analyze the odd prime divisors
   of this recurrence, and find the corresponding structured determinant
   for \(r>1\).

These questions concern hypergeometric arithmetic, scheme geometry,
arithmetic lifts, and tensor rank rather than new fixed-field SIC
counterexamples.
