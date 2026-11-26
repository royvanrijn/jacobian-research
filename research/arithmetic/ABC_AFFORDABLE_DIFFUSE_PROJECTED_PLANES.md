# A strict count on affordable projected radical planes

A plane on two projected power roots and the original third radical
gives a saving at the balanced endpoint. The third radical has short
numerical height at the scale of the projected roots, even when its
coefficient in the plane is large. Reduce the plane modulo a packet
on one projected letter. A narrow valuation window on the other
letter then has a short remaining multiplier. This uses the original
radicals to recover primitive rows and allows arbitrary exact
valuation spectra on the free letters.

## The projected plane criterion

For a prime pool \(P\), define
\[
 {\cal R}_{q,P}(U)
  =\prod_{p\notin P}p^{\lfloor v_p(U)/q\rfloor},
 \qquad
 {\cal S}_P(T)=\operatorname{rad}_{\bar P}(T).
\]
Here \(q\) is a positive integer. Work on a shell
\(X/2<c\le X\) of original positive primitive rows \(a+b=c\),
with \(\operatorname{rad}(abc)\le X^\lambda\), and write
\(n=X^\lambda\). A record specifies an ordering \(U,V,T\) of
the original letters, an integer \(q\), a pool \(P\) of
prime product at most \(X^{B_0}\), and a nonzero integral
normal \((A,B,C)\) of coefficient height at most \(X^{B_0}\).
It covers only original rows whose actual projected coordinates
\[
 x={\cal R}_{q,P}(U),\qquad
 y={\cal R}_{q,P}(V),\qquad
 z={\cal S}_P(T)
\]
satisfy \(Ax+By+Cz=0\). All three coordinates are pairwise
coprime. The order and the pool may differ between records.

**AP1.** Fix \(B_0>0\). There exist \(\eta,\delta,\lambda_0>0\)
such that at every fixed \(0<\lambda<\lambda_0\), a family covered
by at most \(n^\eta\) such records with \(q\lambda\le\eta\)
has size
\[
 O_{B_0,\lambda,\varepsilon}
       \bigl(X^{(1/3-\delta)\lambda+\varepsilon}\bigr)
                  \qquad(\varepsilon>0).
\]
The bound is uniform in the record values. In particular, an
endpoint test may allow \(q\) to grow as \(\lambda\to0\)
provided \(q\lambda\to0\). The free letters may have
arbitrary exact valuation layers.

The limiting form is useful for the proof. If (AP1) fails, take
\(\eta,\delta,\lambda_0\to0\) to get near-endpoint families on
budgets tending to zero, with \(q\lambda\to0\) and at most
\(n^{o(1)}\) records. Pigeonholing leaves one record with original
row count \(n^{1/3-o(1)}\). At each fixed budget \(q\) ranges
over finitely many integers, so one may fix it on an inner
height subsequence. The known limiting bound supplies the upper
endpoint exponent. Throughout, record costs and exponents are
taken in the height limit at fixed budget, before passage to the
small-budget subsequence.
[Limiting one-third bound and height convention](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)

Use the anatomic factorization chosen as in the proof of BBLT
Lemma 2.2 with coefficient loss \(o(\lambda)\). Each exceptional
initial coefficient and the one bundled root on a letter have
height \(X^{o(\lambda)}\) in this sequence. Fix these at vanishing
normalized cost. Outside their prime supports, the remaining
roots are squarefree and use the actual exact degree \(e\) at
their primes. Enlarge \(P\) by these fixed supports and by the
nonzero plane coefficients. The enlarged pool still has polynomial
prime product uniformly. For every fixed \(\chi>0\), rows with
at least \(\chi\lambda\log X\) of original logarithmic radical
on this enlarged pool have a strict saving by (EB5). Thus one
may remove them for a cutoff \(\chi\to0\) sufficiently slowly,
retaining the endpoint exponent. Fix all original parts supported
on the enlarged pool at subpower cost. The old projected
coordinates are fixed multiples of the new outside coordinates.
Absorb those multiples into \(A,B,C\) and use the outside
coordinates henceforth. The new coefficients are supported on
the enlarged pool. Boxing the outside and full original root
sizes costs a subpower factor. Their normalized degree measures
agree in the limit, since only \(o(\lambda\log X)\) of original
root logarithm has been removed.
[BBLT, Lemma 2.2 and its exact-slot construction](https://arxiv.org/html/2410.12234v2#S2)
· [Fixed-support enumeration and the prescribed-pool discount](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

The resulting original box satisfies the balanced endpoint
constraints. Each original letter has normalized radical mass
\(1/3\), has positive mass in a compact interval of positive
scaled degrees, and each pair of letters has combined nonatomic
mass at least \(1/3\). Normalizing the measures by \(\lambda\)
in place of the box root mass gives the same limits, since the
box root mass divided by \(\lambda\) tends to one at saturation.
These facts apply to the retained original rows regardless of
their projected equation.
[Balance and pairwise diffusion](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)
· [Positive compact scales in every endpoint letter](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#positive-scales-in-every-endpoint-letter)

If one coefficient is zero, the two coordinates with nonzero
coefficients have a fixed ratio and hence are fixed separately
by pairwise coprimality. If one of them is \(z\), this fixes the
original \(T\)-support outside the pool. Otherwise both projected
power roots are fixed. An original positive compact band on
either corresponding complete letter lies in the projected
support, since \(q\lambda\to0\). In either case add the fixed
support to the pool and apply (EB5) with a fixed positive
discount, contradicting saturation. If only one coefficient is
nonzero, there are no rows. It suffices to take all coefficients
nonzero below.

## A short multiplier from the third radical

Interchange \(x,y\) if necessary so that the original \(V\)
measure has nonatomic mass at least \(1/6\). This follows from
the pairwise bound on the two original letters underlying
\(x,y\), so the argument works for any orientation of \(T\).

Set \(Y=X^{1/q}\), \(L'=q\lambda\). Then \(n=Y^{L'}\).
An exact valuation slot \(e\) on \(U\) or \(V\) contributes
its outside squarefree root raised to \(\lfloor e/q\rfloor\)
in \(x\) or \(y\). The projected degree measured in the new
units differs from the original one by at most \(L'\):
\[
                0\le e\lambda
                   -L'\lfloor e/q\rfloor<L'.
\]
Thus positive compact-band masses and the sizes of the
underlying squarefree roots are preserved under this change.
Everything below uses actual powers in the projected
equation. The original third radical is at most \(n\).

Choose \(0<b<b'\) and a compact subset of \((b,b')\)
disjoint from all atoms of the \(V\)-measure carrying mass
greater than \(1/12\). Such a set exists by inner regularity
after deleting the countable atom set and restricting the
nonatomic measure to compact positive scales. Use shrinking
neighborhoods of this set in the original scaled positions
\(e\lambda\), with continuity boundaries chosen to retain its
mass in the limit. The maximum mass of a single coordinate in
these neighborhoods tends to zero. Otherwise a subsequence
would put an atom on the compact set. In particular, for any
fixed \(0<t\le1/3\), a block of these slots can be chosen with
total mass \(t/8+o(1)\), by stopping when its mass first
reaches \(t/8\). Its scaled positions remain in \((b,b')\).
The block depends on the box, and its actual outside root
values can be fixed at cost \(n^{t/8+o(1)}\), using an
arbitrarily small exponent allowance at each fixed budget.

Take \(C_1=16\) and an integer \(\ell\ge32b'/b\).
On the original \(U\)-letter choose a fixed interval
\(I=[r,r+\sigma]\) inside positive compact scales, with
continuity endpoints and positive limiting mass \(t\).
Require its width to satisfy
\[
         \sigma< \min\{r/16,\ c_1 b/(1024\ell^2)\},
\]
where \(c_1\) is the absolute constant in (UW) for one
narrow factor. Such an interval exists because the original
\(U\) measure has positive mass on a compact positive band.
Choose the width bound first using a fixed positive lower
endpoint of a slightly enlarged band. A finite cover by
sufficiently narrow continuity intervals then gives some
positive mass. This argument also allows the interval mass
to arise from an atom of the \(U\) measure.

Use the preceding \(V\)-block at this \(t\), and denote its
fixed product in the projected \(y\) by \(D\). All its prime
valuations are at least \(k=\lfloor b/(2L')\rfloor\) at
sufficiently late budgets. Its numerical weight \(w=\log_Y D\)
has bounds \(bt/8-o(1)\le w\le b't/8+o(1)\). In particular
\(kL'\to b/2\). Freeze the original \(U\) slots outside \(I\),
including any outside root with degree \(e<q\), at cost
\(n^{1/3-t+o(1)}\). Set
\[
 \begin{split}
   j_0 &=\lfloor r/L'\rfloor,\qquad
   R=\prod_{e\lambda\in I} u_e,\qquad
   Q=\prod_{e\lambda\in I}u_e^{\,\lfloor e/q\rfloor-j_0},\\
   x&=x_0 R^{j_0}Q,\qquad
   R\le n^{t+o(1)},\qquad Q\le Y^{\sigma t+o(1)} .
 \end{split}
\]
Here \(u_e\) is the outside squarefree root for exact degree
\(e\) on \(U\), \(x_0\) is fixed, and all powers in \(Q\)
are nonnegative. Rounding the projected orders accounts for
an \(o(1)\) change in the displayed \(Y\)-exponent.

Modulo the actual divisor \(D\mid y\), the projected plane
gives
\[
                  (Q/z)R^{j_0}\equiv -C/(Ax_0)\pmod D .
\]
All variables and coefficients in the congruence are units
by original coprimality and the support convention. The rational
wide multiplier \(Q/z\) has height at most \(Y^{E+o(1)}\),
with \(E=\sigma t\), since \(z\le Y^{L'}\).
The uniform one-root estimate (UW1) applies with
\(C_1=16\), \(E/w\) as the wide exponent relative to \(D\),
and the chosen \(\ell\). Indeed,
\((kL')t/w\le4+o(1)<C_1\), and
\(\ell E/w\le8\ell\sigma/b+o(1)<c_1/(\ell C_1)\).
Both inequalities have fixed margins. Its exponent to base
\(n\) is at most
\[
    \frac{2w}{\ell k L'}+\frac{E}{j_0 L'}+o(1)
          \ \le\ \frac{t}{16}+\frac{t}{16}+o(1).
\]
The first estimate uses \(\ell\ge32b'/b\), and the second
uses the bound on \(\sigma/r\). This is a direct use of
(UW1), including its allowance for overlap between \(Q\)
and \(R\).
[Product-bound modulus estimate (UW1)](ABC_QUANTITATIVE_DEGREE_WINDOWS.md#a-product-bound-for-the-admissible-wider-height)
· [Underlying movable-window argument](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-movable-window-on-one-tall-letter)

Each pair \((R,Q/z)\) gives \(x/z=x_0 R^{j_0}Q/z\) exactly.
Because \(x,z\) are coprime, the reduced fraction gives
both integers. The projected power root \(x\) retains
every outside prime on \(U\) with original valuation at
least \(q\). Every omitted outside prime of \(U\) was
already fixed with its slot outside \(I\). The recovered
\(z\) gives the entire outside support of \(T\).
Consequently the original \(U,T\) have at most subpower
many choices on these now prescribed supports at height
at most \(X\). The original signed sum then determines
\(V\). The omitted low valuations of \(V\) enter only
through this last exact sum. Restoring the frozen costs
gives the normalized bound
\[
                 1/3-t+t/8+t/16+t/16+o(1)
                      =1/3-3t/4+o(1).
\]
The fixed \(t>0\) depends on the limit profile. It
contradicts every putative near-endpoint sequence of
counterexamples to (AP1), which gives the uniform
existence assertion without a uniform interval width
or window mass in advance.

## Polynomial variants and the remaining extraction

**AP2.** The same criterion with all \(A,B,C\ne0\) holds for
\(Ax+By+CF(z)=0\), if each record also specifies a nonconstant
integer polynomial \(F\) of bounded degree and coefficient
height at most \(n^{B_1}\), for a fixed \(B_1\).
Constants may also depend on the degree and \(B_1\).
When extending the pool, the old \(z\) becomes a fixed
multiple \(z_P z\) of the new one, with \(z_P\le n\).
Replace \(F(Z)\) by \(F(z_P Z)\), which retains the
stated type of height bound with \(B_1\) increased by
at most its degree. There are boundedly many possible
\(z\) for \(F(z)=0\); fixing them and applying (EB5) to
the original \(T\)-support excludes that branch at saturation.
For \(F(z)\ne0\), the modulus congruence uses \(Q/F(z)\) in
place of \(Q/z\). Its height bound persists since
\(|F(z)|\le n^{O_{B_1,\deg F}(1)}=Y^{o(1)}\).
The polynomial value is a unit modulo \(D\), as the
projected equation gives \(CF(z)\equiv -Ax\pmod D\).
Any prime dividing both \(x\) and \(F(z)\) must divide
\(By\), hence \(B\) by pairwise coprimality. Such a prime
has already been removed from \(x\). The exact ratio
therefore gives \(x,F(z)\), and the last value has
boundedly many preimages \(z\). This proves (AP2).
It includes affine planes when the constant divided by
the third coefficient is rational of height \(n^{O(1)}\),
by clearing that denominator in the polynomial and its
scalar coefficient.

For the unrestricted endpoint, it suffices to obtain
these records at vanishing normalized cost on
\(n^{1/3-o(1)}\) actual original rows of each putative
balanced saturation sequence, with \(q\lambda\to0\).
One precise source would be a large simultaneous residue
branch of the projected coordinates. If \(D_i\mid T_i\)
are prescribed pairwise coprime moduli and the complementary
unit ratios of the projected coordinates are fixed modulo
them, the determinant criterion (SR2) forces a plane
as soon as the product of the moduli exceeds six times
the product of the three projected coordinate bounds.
More generally a rank-three branch has at most
\(O(H_1H_2H_3/(D_1D_2D_3))\) distinct projected
vectors. For comparison with a count of original rows,
each fixed projected vector has at most
\(n^{\zeta+o_X(1)}\) preimages when the outside
\(e<q\) radical on the cheaper of the two power-root
letters is capped by \(n^\zeta\). To see this, enumerate
that omitted radical and the two original letters
visible from the selected root and \(z\), as in the
recovery above. This allows an occupancy version of the
rank criterion with this multiplicity included.

In particular a partition into simultaneous residue branches
would suffice if branches of size \(n^{1/3-o(1)}\)
carry \(n^{1/3-o(1)}\) rows and fall on the rank-two side.
The global limiting bound then bounds the number of such
branches by \(n^{o(1)}\). Any plane normal can be chosen
at polynomial height by the cross product of two independent
projected vectors. In rank one, a single projected vector
determines a suitable plane at polynomial height. This
formulation only requires economical planes or sufficiently
large rank-two branches. The actual corrected complementary
power-root residues must still be controlled to obtain them.
[Simultaneous residue lattice and the rank criterion](ABC_SHALLOW_REPAIR_RADICAL_PLANES.md#rank-and-occupancy-of-simultaneous-residue-records)
· [Projected roots and their corrected residues](ABC_DIVISIBLE_POWER_ROOT_PLANES.md#costs-and-arbitrary-high-valuations)

## A single mixed residue suffices

**AP3.** Use the retained boxes after the fixed-pool preparation
and support fixation of (AP1). Use its height-first setup and
limiting profiles, with \(q\lambda\to0\), and assume no plane.

For any chosen \(T\), orient \(U,V\) so that \(V\) has
nonatomic mass at least \(1/6\). Choose the compact scales
\(b,b'\), exponent window \(I=[r,r+\sigma]\) of \(U\)-mass
\(t>0\), and integer \(\ell\) as in the preceding proof.
There is a block on \(V\) that defines, for each row,
a divisor \(D\mid{\cal R}_{q,P}(V)\) of numerical weight
\(bt/8-o(1)\le\log_Y D\le b't/8+o(1)\), with every prime
valuation at least \(k=\lfloor b/(2q\lambda)\rfloor\).
Its dictionary has at most \(n^{t/8+o(1)}\) values.
For each prescribed \(D\) and unit \(K\bmod D\), there are
at most
\[
                          n^{1/3-7t/8+o(1)}X^\varepsilon      \tag{AP3}
\]
rows satisfying \(x/z\equiv K\bmod D\), where
\(x={\cal R}_{q,P}(U)\) and \(z={\cal S}_P(T)\).
The \(o(1)\) covers the outer box allowances and all inner
exponent allowances can be taken arbitrarily small at fixed
budget. The same count works for any prescribed
\(k\)-full divisor \(D\mid V\) in the given numerical
range. In that form no construction of \(D\) from projected
valuation slots is required.

The block construction above only uses the measure on
\(V\), so still applies. Freeze the \(U\) slots outside the
window as before at cost \(n^{1/3-t+o(1)}\).
Then \(x=x_0R^{j_0}Q\), and the residue becomes
\((Q/z)R^{j_0}\equiv K/x_0\bmod D\). The unit,
height and fullness conditions of (UW1) hold with the same
strict margins. At fixed \(D,K\) the count of the resulting
wide and narrow pairs costs at most \(n^{t/8+o(1)}\).
Recover \(x,z\) and their original letters exactly at
subpower cost as before. This proves (AP3).

Consequently any collection of \((D,K)\) records of
exponent less than \(7t/8\) by a fixed margin covers a
power-smaller set of rows. If \(D\) comes from the chosen
block and each \(D\) has its own dictionary of residues
of exponent less than \(3t/4\), with a fixed margin,
the same conclusion holds. Equivalently, a full-exponent
family must have at least \(n^{7t/8-o(1)}\) occupied
\((D,K)\) records, and cannot place a full-exponent
subfamily on dictionaries below those thresholds.

The arithmetic object whose dispersion this tests can
be written using the original correction factors. Let
\(U=x^q E_U\) and \(T=z^q E_T\), with \(E_T\) rational
and both corrections units modulo \(D\). The original
signed sum implies
\[
                    (x/z)^q\equiv \pm E_T/E_U\bmod D .
\]
At fixed \(q\) and polynomial-height \(D\), each unit
power residue has subpower many original unit root residues.
Here \(q\) is bounded at each fixed budget, so this cost
is inner subpower throughout the limiting sequence.
A dictionary of the actual corrected ratio modulo this
one divisor therefore suffices. The needed extraction
concerns one complementary residue, even in the absence
of a simultaneous rank-two branch.
