# Two-target counts from common summand divisors

A pair of rows with fixed relative odds has two actual target
values that are linear forms in the same primitive pair.
A divisor common to the first summands, or common to the
second summands, fixes the quotient of those target values
modulo the divisor. This permits reconstruction on the
roots of the target quotient. In particular parity needs
only two target positions when the summand divisibility
has enough numerical weight. Least repair gives a second
application in which the paid record is a residue of the
repair multipliers and the recovered coordinates are
target radicals.

## Relative odds and full summand divisors

Fix \(B_0>0\), \(X\ge2\), signs
\(\sigma,\tau\in\{-1,1\}\), and two ordered positive primitive
triples
\[
 T_i=\sigma U_i+\tau V_i>0,\qquad
   U_i,V_i,T_i\le X^{B_0},\qquad i=1,2.
\]
Prescribe the relative odds
\[
        r=\frac{U_2V_1}{U_1V_2}\ne1
\]
and positive integers \(F,G\) such that \(F\mid U_1,U_2\)
and \(G\mid V_1,V_2\). If a candidate exists then \(F,G\)
are coprime and both targets are units modulo \(Q=FG\).
The signs and the order of the two rows are fixed.

Put \(x=\gcd(U_1,U_2)\), \(y=\gcd(V_1,V_2)\), and write
\(U_i=d_i x,\ V_i=e_i y\), with \(\gcd(x,y)=1\).
Outside the primes in the reduced numerator and denominator
of \(r\), the valuations of \(U_i/V_i\) agree. In each
primitive row at most one of the two valuations can be
positive, so the difference determines both. Thus the four
cofactors \(d_i,e_i\) are supported on a polynomial-height
integer fixed by \(r\), and fixing them costs
\(X^{o_X(1)}\) uniformly in the odds. Here and below these
subpower counts use arbitrary inner exponent allowances.
With the cofactors fixed, set
\[
            L_i(x,y)=\sigma d_i x+\tau e_i y,\qquad
            \Delta=d_1e_2-d_2e_1\ne0 .
\]
The determinant is nonzero by \(r\ne1\). The quotient
\(T_1/T_2=L_1(x,y)/L_2(x,y)\) determines \(x/y\) and hence
the positive coprime pair \(x,y\). Also \(\gcd(T_1,T_2)\)
divides \(|\Delta|\): it divides both \(\Delta x\) and
\(\Delta y\) by eliminating a term in the two equations.
[Fixed odds and common factors](ABC_MULTI_TARGET_MIXED_KERNEL_COUNTS.md#fixed-cross-ratios-and-the-common-factors)
· [Uniform supported-value enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

Since \(F\mid x\) and \(V_i\) is a unit modulo \(F\),
\(y,e_1,e_2\) are units there. Similarly \(x,d_1,d_2\) are
units modulo \(G\). Therefore, in unit group notation,
\[
       \frac{T_1}{T_2}\equiv \frac{e_1}{e_2}\pmod F,\qquad
       \frac{T_1}{T_2}\equiv \frac{d_1}{d_2}\pmod G.        \tag{TL1}
\]
Combining the two fixes one residue modulo \(Q\). No prime
of \(Q\) needs to be removed when it divides \(\Delta\)
or the numerator or denominator of \(r\).

## The mixed target defect count

For a positive rational \(t\) and integer \(q\ge2\), let
\[
       \operatorname{rad}_q(t)=\prod_{p:\ q\nmid v_p(t)}p .
\]
Prescribe a prime pool \(P\) of polynomial prime product,
allowed to depend on the odds and the original locks.
Write a subscript \(\overline P\) for the product outside
this pool. Impose individual target heights \(T_i\le Y_i\)
and
\(\operatorname{rad}_{q,\overline P}(T_1/T_2)\le B\).

**TL2.** Uniformly in the fixed data the pair count is
\[
 \ll_{\varepsilon,B_0,B_P}
  B X^\varepsilon\left(1+\frac{(Y_1Y_2)^{1/q}}{FG}\right),  \tag{TL2}
\]
where the product of the primes in \(P\) is at most \(X^{B_P}\).
For \(q=2\), the defect of \(T_1/T_2\) is exactly the
squarefree kernel of the product \(T_1T_2\), including
at their common primes. At general \(q\) one may use
\(T_1T_2^{q-1}\). Individual \(q\)-defect caps on the two
targets also suffice by multiplying them.

Indeed, put \(T_1/T_2=J(u/v)^q\) with \(u,v\) positive
coprime integers, by extracting powers separately from the
reduced numerator and denominator. The rational \(J=J_1/J_2\)
is reduced and its prime powers have exponents less than
\(q\). The retained prime support of \(J\) has at most \(B\)
integer radical records. Both \(J_i\le X^{B_0}\), so fixing
\(J\) costs \(O(BX^\varepsilon)\), with an appropriately
smaller inner allowance. For each record
\[
             u\le A_1=(Y_1/J_1)^{1/q},\qquad
             v\le A_2=(Y_2/J_2)^{1/q}.
\]
All the quantities in the rational factorization are units
modulo \(Q\). By (TL1), \((u/v)^q\) is a fixed unit there.
At a prime power \(p^a\) the kernel of the \(q\)-th power
map on units has size at most \(2q\), using the cyclic
unit group for odd \(p\) and its two cyclic factors at
\(2^a\). Thus at every fixed \(q\) there are
\((2q)^{\omega(Q)}=X^{o_X(1)}\) possible residues of \(u/v\).
These include all residue primes dividing \(q\).

In any one residue, distinct reduced pairs have a nonzero
cross determinant divisible by \(Q\). After scaling the
two coordinates by \(A_1,A_2\), their angles in the first
quadrant differ by at least \(Q/(2A_1A_2)\), since their
norms are at most \(\sqrt2\) and the sine of the angle
difference is at most that difference. Thus there are at
most \(1+\pi A_1 A_2/Q\) fractions in that residue. The
fraction and \(J\) recover the complete pair. The estimate
is uniform in \(q\) as stated. Above a large fixed cutoff
depending on \(\varepsilon,B_0\), enumerate the at most
\(X^{2B_0/q}\) pairs \((u,v)\) directly and use the modular
argument only for the finitely many smaller orders.
Allowing the choice of \(q\) to vary on a tuple costs only
a subpower: orders above
\(\lceil B_0\log_2 X\rceil+1\) give the same defect and
both roots equal to one. The choices of cofactors and all
fixed-support allowances are included in (TL2).

If the actual values of \(\log(T_1/T_2)\) lie in one
prescribed interval of length \(\ell\), the numerator
\((Y_1Y_2)^{1/q}\) in the second term of (TL2) can further
be multiplied by \(\min(1,\ell/q)\). For fixed \(J\)
the root quotients have logarithmic range of length
\(\ell/q\). The derivative of the angle of the rescaled
pair with respect to \(\log(u/v)\) has absolute value at
most \(1/2\), by differentiating \(\arctan(c e^{-t})\)
for the constant \(c=A_1/A_2\). The same packing in this
shorter interval gives the improved estimate. The
interval centre may depend on the fixed records.

## Repair multiplier residues and target support recovery

For fixed integer \(M\ge2\) write
\[
  R_i=\operatorname{rad}(T_i),\qquad
  H_i=H_M(T_i)=R_i^M W_i,\qquad Z_i=H_i/T_i .
\]
Here \(W_i\) is the excess valuation part
\(\prod_{p\mid T_i}p^{(v_p(T_i)-M)_+}\).
The \(Z_i\) are the actual least-repair multipliers.
Assume \(R_i\le R_i^*\).

**TL3.** Suppose for the given odds and summand locks
the residue of \(Z_1W_2/(Z_2W_1)\) modulo \(Q\) belongs
to a prescribed dictionary of at most \(D\) classes.
The pair count is
\[
      \ll_{M,\varepsilon,B_0}
          D X^\varepsilon\left(1+\frac{R_1^*R_2^*}{FG}\right).
                                                               \tag{TL3}
\]
If the dictionary instead concerns only \(Z_1/Z_2\bmod Q\)
and \(T_i\le Y_i\), multiply this bound by
\((Y_1Y_2)^{1/(M+1)}\). One can use any smaller bound on
the number of possible classes of \(W_1/W_2\bmod Q\)
in place of that last factor.

For the first assertion all quantities are again units.
Equation (TL1) and the recorded class fix a residue of
\[
       (R_1/R_2)^M=
          (T_1/T_2)\,\frac{Z_1W_2}{Z_2W_1}.
\]
The root-class enumeration and angular spacing used in
(TL2) bound the number of reduced radical ratios by
\(X^{o_X(1)}(1+R_1^*R_2^*/Q)\) per class.
The shared radical of the targets divides
\(\operatorname{rad}|\Delta|\). Fix it from a divisor
list at subpower cost. Then the reduced ratio determines
both \(R_1,R_2\). With the two supports fixed, the original
complete target values have only subpower many
possibilities, and they recover the original pairs by
the two independent linear equations. No summation
over common multiples of primitive quotients occurs.

For the second assertion the tails satisfy
\[
          \operatorname{rad}(W_1W_2)
                    \le (Y_1Y_2)^{1/(M+1)}.
\]
Enumerate this radical and then the two tail values on
its fixed support at subpower cost. They specify the
needed tail ratio class. This proves the second bound.
If \(R_i^*\le X^{C_iL}\) and \(Q\ge X^\omega\) with
fixed \(\omega>0\), the projective factor in (TL3) is
bounded at all sufficiently small fixed \(L\).
If \(M=\lceil K_0/L\rceil\) with \(K_0>0\) fixed
and \(Y_i\le X^{h_i}\), the second version costs at
most \(n^{(h_1+h_2)/K_0}\) in tails, where \(n=X^L\).
The order is fixed inside the height limit. The record
of the corrected ratio has no tail cost.
[Least target repair of a primitive quotient](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

## Graph count and a positive sector

Use a simple outside graph with one primitive witness per
edge and with the complete triple determining its outside
labels. Let \(R,C,E\) be the row, column, and edge counts.
For any ordered pair of distinct rows suppose the
relative-odds record and compatible full summand lock
record belong to a joint dictionary of size at most \(K\).
These dictionaries can vary with the ordered row pair.
Use (TL2) with the same upper bound \(B X^\varepsilon A\)
uniformly on the records, or (TL3) with this bound in
place of its right-hand side. Then
\[
         E\ll C+R\bigl(CKBA\bigr)^{1/2}X^\varepsilon .       \tag{TL4}
\]
For \(E\ge2C\), Cauchy–Schwarz gives at least
\(E^2/(2C)\) ordered pairs of distinct rows sharing a
column. Their number is at most \(R^2 K B X^\varepsilon A\)
by the per-record pair count, since each tuple determines
its column. This proves the claim. In the repair version
absorb the repair residue dictionary and the optional
tail charge into \(B\). In the ordinary defect version
\(A=1+(Y_1Y_2)^{1/q}/Q_{\min}\), using a lower bound
\(Q_{\min}\) for the size of the summand lock product.
The repair version uses \(A=1+R_1^*R_2^*/Q_{\min}\).

Here is one outside-profile consequence in the notation
of (OA1), with normalized count exponent \(d\) and
nonatomic masses \(c_U,c_V,c_T\). Assume \(c_T>0\).
Suppose actual whole prime-power blocks can be locked on
the diffuse pieces of \(U,V\) at normalized costs
\(b_U,b_V\) and combined numerical weight \(w\), with
the tolerances tending to zero in the limiting profile.
For example, such blocks can be chosen in compact
positive valuation bands with vanishing individual
root masses after excluding atomic windows. The cost
and the removed outside mass on each letter are the same
for these locks. Write \(s=b_U+b_V\). Fix the two locks
at cost \(s\), and a positive-weight full block on \(T\)
at arbitrarily small cost. The one-projection form of
(TR2)–(TR3) gives pair thinning and relative-odds
dictionaries with arbitrarily small normalized exponent,
as used in (MTK2) and (PT4). The remaining graph has
edge exponent at least \(d-s\) and part exponents at most
\(c_U-b_U,c_V-b_V\).

If the complete targets have height at most
\(X^{h_T+o(1)}\), and actual \(q\)-defect cap
\(n^{f+o(1)}\) on each target, then for
\(w>2h_T/q\) with fixed numerical margin, (TL4) gives
\[
 d\le s+\max\!\left(c_V-b_V,\
                c_U-b_U+\tfrac12(c_V-b_V)+f\right).         \tag{TL5}
\]
The same bound holds with \(f\) replaced by half the
joint pair-defect cap exponent when such a pair cap is
available. Exchange \(U,V\) to choose the cheaper
orientation. For instance at minimum nonatomic mass
\(c_U=c_V=1/6\), the inequalities \(b_U<1/6\) and
\(b_V+2f<1/6\) suffice when the stated weight margin
is available. With \(q=2\) this uses the target parity
condition only and requires no parity condition on
either summand. For \(q\ge3\) the ordinary numerical
weight threshold is lower. The formulation is valid
away from the minimum nonatomic mass whenever the
right side of (TL5) is strictly below the tested
count exponent. Strict margins give power savings
in the corresponding families at sufficiently small
fixed budgets, by taking all box and inner height
allowances below those margins.

For example, take \(q=3\), \(h_T\le1\), and
\(c_U,c_V\le1/6+1/500\). Suppose a diffuse block on
\(U\) of cost \(3/20\) has weight at least \(7/10\).
Take \(b_V=0\) and \(f=1/15\). Formula (TL5) gives
\[
 d\le \max\left(\frac{19}{60}+\frac1{500},\
                \frac{19}{60}+\frac3{1000}\right)<\frac13.
\]
The three-target estimate (PT5) alone does not give
a saving at this kernel cap for equal \(c_U=c_V=1/6\).
For a parity version of the same example, allow
\(h_T\le 2/3\) and use \(q=2\). The same block
weight clears \(h_T\) and the same count follows.
These are sector conditions on a block from an actual
profile. Neither block is assumed to exist on every
balanced profile.

In the least-repair version take any fixed positive
weight \(w\) from the summand locks. With a corrected
residue dictionary of exponent at most \(\rho\)
uniformly for each odds and lock record, (TL5) holds
with \(f=\rho/2\). With the uncorrected \(Z_1/Z_2\)
dictionary use \(f=\rho/2+h_T/K_0\) at order
\(M=\lceil K_0/L\rceil\). Thus a sufficiently small
*residue* dictionary suffices even when the three
original prime valuation profiles are arbitrary.
Put \(u=c_U-b_U\), \(v=c_V-b_V\). Conversely, if
\(d-s>v\), the maximum size of a corrected residue
dictionary over fixed odds and row-pair records must
have normalized exponent at least
\[
                         2(d-s)-2u-v .
\]
Subtract \(2h_T/K_0\) for the uncorrected multiplier
ratio. At the minimum nonatomic endpoint, with locks
of costs tending to zero, the corrected threshold is
\(1/6\). The argument gives no uniform upper estimate
for that complexity.

[Complete relative odds (TR2–TR3)](ABC_DIFFUSE_FACTOR_GRAPH.md#counting-complete-cross-ratios-by-a-full-modulus)
· [One projection and outside labels (OA1)](ABC_DIFFUSE_FACTOR_GRAPH.md#a-single-anchor-in-one-projection)

The mixed-kernel counts (MTK1) and (PT1) use four targets
for parity or three for higher prime orders, with no
common divisor condition on the summands. Here two
target positions suffice because the summand locks
control their quotient. The argument uses unit power
classes and projective spacing, including the primes
where the relative-odds coefficients have bad reduction.
The least-repair version differs from the Pell and
Thue recovery from exact target differences in (WD3):
it records a modular multiplier ratio and recovers the
original targets from their radical ratio. Both counts
leave a uniform source of economical repair residue
records open.
[Mixed squarefree kernels (MTK1)](ABC_MULTI_TARGET_MIXED_KERNEL_COUNTS.md#the-mixed-kernel-count)
· [Three-target prime-order count (PT1)](ABC_THREE_TARGET_PRIME_ORDER_STARS.md#count-for-three-actual-targets-at-prescribed-odds-quotients)
· [Least-repair target difference recovery (WD3)](ABC_REPAIRED_WEIGHTED_TARGET_DIFFERENCES.md#pair-records-and-their-count)

## A direct weighted residue count for primitive quotients

Fix a positive divisor \(D\) on a designated complete
letter \(T=\sigma U+\tau V\) of an original primitive row.
Fix \(M\ge2\) and least-repair the two free letters,
writing \(R_W=\operatorname{rad}(W)\) and
\(WZ_W=H_M(W)=R_W^M E_W\) for \(W\in\{U,V\}\).
Here \(E_W\) is the excess valuation part of the repair.
Work on a shell of height at
most \(X\), with \(n=X^L\), \(R_U R_V\le n^{C_0}\)
for fixed \(C_0>0\), and \(D\ge X^\omega\) for fixed
\(\omega>0\). Take \(C_0L<\omega\).

**TL6.** On the original rows in this divisor record,
the map to the corrected multiplier residue
\[
        \rho_D=\frac{Z_U E_V}{Z_V E_U}
                     \pmod D
\]
has fibers of size \(O_{M,\omega,C_0,L,\varepsilon}
(X^\varepsilon)\). The constants are uniform in \(D\)
and the other restrictions on the rows. The proof counts
primitive quotients and uses no common-multiplier
parameter.

All free letters and all their repair factors are units
modulo \(D\) by original primitivity. The signed
equation gives \(U/V\equiv-\tau/\sigma\bmod D\). Hence
\[
         (R_U/R_V)^M\equiv(-\tau/\sigma)\rho_D\bmod D.
\]
At fixed \(M\) there are \(X^{o_X(1)}\) root residues
modulo \(D\). Use individual dyadic bins for the two radical sizes,
at subpower cost. Within each occupied bin the product
of the two upper size bounds is at most \(4n^{C_0}\).
The cross determinant is divisible by \(D\) and its
absolute value is at most this product. It therefore
vanishes for sufficiently large heights at fixed \(L\).
A reduced fraction gives the two individual
radicals exactly, because \(U,V\) are coprime. Fixing
the two supports leaves \(X^{o_X(1)}\) choices of the
complete \(U,V\), which determine \(T\). All bounds
keep \(L,M\) fixed in the height limit.
[Uniform supported-value enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

Consequently a divisor dictionary of at most \(n^s\)
records, each of height at least \(X^\omega\), together
with at most \(n^\rho\) attained corrected residues
per record gives a count
\[
                   O_{\varepsilon}(n^{s+\rho}X^\varepsilon)
\]
with constants depending on the fixed parameters of
(TL6). If the recorded image is instead that of
\(Z_U/Z_V\bmod D\), the same bound has an extra factor
\(n^{(h_U+h_V)/K_0}\) when \(U\le X^{h_U}\),
\(V\le X^{h_V}\), and \(M=\lceil K_0/L\rceil\).
Indeed the tail union radical is at most
\((X^{h_U+h_V})^{1/(M+1)}\), and each of its recorded
values allows only subpower many tail values.

For any covering divisor dictionary let \(I_D\) be the
number of corrected residues attained by actual rows
with the given divisor and orientation. Then the
covered count has the same exponent as
\(\sum_D I_D\), including the finite set of possible
orientations. One inequality is (TL6). For the other,
a row occurs in at most
\(\tau(a)+\tau(b)+\tau(c)=X^{o_X(1)}\) oriented divisor
records. Arbitrary separate restrictions on the three
original radicals and on the primitive height can be
kept throughout. When \(\operatorname{rad}(D)\le n^s\)
it costs at most \(n^{s+o_X(1)}\) to enumerate the
whole divisor at height at most \(X\). The images
\(I_D\) count genuine simultaneous occurrences in
those records.

At any balanced endpoint, the pairwise diffusion
constraints give a letter of positive nonatomic mass.
Select a positive scaled-degree band of that measure
and lock a full block of arbitrarily small fixed
normalized cost \(s>0\) in it. The resulting
\(\omega>0\) may depend on the profile and the chosen
block. Along a height-first sequence tending to the
profile, apply (TL6) in the inner height limit at each
sufficiently late fixed budget with \(C_0 L<\omega\).
Equivalently, the affordable
whole-coordinate pair locks supply a single divisor
with positive numerical weight by taking the
heavier member. Thus the weighted criterion applies
to the entire balanced endpoint whenever these
attained images have a uniform power saving. Counting
all residue classes, or just independent dictionaries
for the two free radicals, supplies no such estimate.
For every fixed \(D,M\) the corrected image has the
same exponent as the dictionary of original radical
ratios in that record, by the same residue
reconstruction in both directions. Choosing a repair
order changes the modular representation of that
incidence and leaves a separate arithmetic bound to
be proved.
[Pairwise diffusion constraints](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)
· [Affordable whole-coordinate locks](ABC_WIDE_COFACTOR_CONGRUENCES.md#affordable-whole-coordinate-pair-locks-at-every-balanced-endpoint)
· [True support occurrence](ABC_MINIMAL_REPAIR_SUPPORT_OCCURRENCE.md#exact-reduction-to-occurring-support-unions)
