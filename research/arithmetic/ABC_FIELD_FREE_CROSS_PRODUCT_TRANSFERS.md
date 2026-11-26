# Field-free product and comparison transfers for primitive pairs

The determinant of two products of actual scaled operands
is the product of two reduced shifts. Its primitive
additive triple can recover the original pair by a
split quartic count. The only additional record required
in this recovery is the common parity kernel that
disappears from the auxiliary triple. This record is
cheap on typical pairs of an original near-full family.

## An ordinary primitive triple from the product difference

**FP1.** Take two distinct positive primitive rows
\(a_i+b_i=c_i\le X\). Exclude pairs obtained by swapping
the two summands of the same row. Set
\[
 G=(c_1,c_2),\quad u=c_2/G,\quad v=c_1/G,\quad
 (A,B)=(ua_1,ub_1),\quad (C,D)=(va_2,vb_2).
\]
The common sum is \(W=A+B=C+D\).
Write \(h=A-C,\ j=A-D\). Neither is zero, because
an equality would identify or interchange the reduced
primitive proportions. The identity
\[
                 AB-CD=-h j
\]
gives the auxiliary triple with positive values
\[
                  P=AB/d,\quad Q=CD/d,\quad |P-Q|,
                     \qquad d=(AB,CD).
\]
The three absolute values are put in additive order.
Here \(d=(a_1b_1,a_2b_2)\), since \(u,v\) are coprime
and every prime of \(u\) is prime to \(va_2b_2\),
with the reverse assertion for \(v\). Thus \(d\)
is prime to \(uv\).

Assume an actual common divisor \(M\mid G\) has
\(M\ge X^\nu\) for fixed \(\nu>0\) and
\(a_1\equiv a_2\bmod M\). This is an ordinary absolute
residue cell with a shared target lock. The auxiliary
has fixed positive numerical height unless
\(a_1b_1=a_2b_2\). In fact \(u^2\mid P,\ v^2\mid Q\),
\(d\) is a unit at \(M\), and
\[
   M\mid v^2P-u^2Q,\qquad
              |v^2P-u^2Q|\le P Q .
\]
The divisibility follows from both original products
\(a_i b_i\equiv-a_i^2\bmod M\).
If they are unequal the integer on the left is
nonzero, and \(\max(P,Q)\ge\sqrt{M}\).
The equality case consists of at most
\(N X^\varepsilon\) ordered pairs in any original
\(N\)-row family, because fixing \(a_1b_1\) leaves
only divisor-many possibilities for \(a_2,b_2\).
The normalized-lock alternatives \(M\mid h\) or
\(M\mid j\) instead make \(P-Q\) divisible by \(M\)
when the scaled operands are units there. For a
common target lock, this unit condition follows
from that shift divisibility and primitivity.

At the pair-family scale, FP8 below provides a linear auxiliary
height after deleting \(O(NX^\varepsilon)\) pairs, without
requiring a shared lock.

The raw radical mass of the product auxiliary
can be read exactly, including possible prescribed
zero-price support. Write \(\ell(Z)\) for the
normalized radical of the absolute nonzero integer
\(Z\) outside one fixed pool. Put
\[
 \delta=\ell(d),\qquad
 \delta_0=\log_n\operatorname{rad}_{S^c}(d),
 \quad S=\operatorname{supp}(P Q(P-Q))\cup P_0,
\]
where \(P_0\) is the chosen prime pool and
\(\operatorname{rad}_{S^c}\) counts only primes
outside \(S\). Then
\[
 \ell(P Q(P-Q))=
     \ell(AB)+\ell(CD)+\ell(hj)-2\delta-\delta_0 .
                                                               \tag{FP1}
\]
For each prime of \(d\), all three raw values
\(AB,CD,hj\) contain the prime. Normalization
leaves it on exactly one auxiliary role or
removes it altogether. This gives the two
discounts. All other primes have the same
incidence before and after division.

## Subpower recovery from one auxiliary and a parity overlap

**FP2.** With \(P,Q\) prescribed as actual exact auxiliary
values, write \(D_i=\operatorname{sf}_2(a_i b_i)\), the
squarefree kernel of the original product, and
\(D_0=(D_1,D_2)\). Fix a polynomial-height prime pool
\(P_0\). If
\[
     \operatorname{rad}_{(\operatorname{supp}(P Q(P-Q))
                       \,\cup\,P_0)^c}(D_0)\le n^\zeta ,
\]
the original pair count is
\(O(n^\zeta X^\varepsilon)\) at any fixed budget.
The two complete original targets may both vary.
All original caps can be retained in the recovered pairs.

Indeed \(\operatorname{sf}_2(P Q)=\operatorname{sf}_2(D_1D_2)\)
because the scaling and normalization factors in the
product form a rational square. Hence \(P,Q\) give the
symmetric difference of the two parity supports.
The primes on that difference can be allocated to
the two fields at subpower cost. Choose \(D_0\)
outside the prescribed support at cost \(n^\zeta\),
then choose its remaining prime factors on that
support at height-subpower cost. This fixes \(D_1\).

Write \(r=A/C=s/t\) in lowest terms with \(s,t>0\),
and \(x_1=A/W,\ x_2=C/W\). From the auxiliary ratio,
\[
 x_2=\frac{P-Qr}{P-Qr^2},\qquad
 x_1=\frac{s(Pt-Qs)}{Pt^2-Qs^2},
 \qquad
 x_1(1-x_1)=\frac{Pst(Pt-Qs)(t-s)}
                         {(Pt^2-Qs^2)^2}.
\]
The denominator in the first expression cannot
vanish on a genuine pair: if it vanishes, the
numerator vanishes as well, and the two identities
would force \(r=1\). This is excluded by the
nonzero displacement \(h\). All four linear
factors of the numerator have distinct rational
projective zeros since \(P,Q\) are positive
and \(P\ne Q\). The squareclass of
\(x_1(1-x_1)\) is \(D_1\). Thus the
projective direction \((s:t)\) occurs on the
split quartic double cover of the line
\[
                  D_1 z^2 = P s t(Pt-Qs)(t-s).
\]
The affine equation is understood projectively,
with the square coordinate of weight two.
Move a branch point to infinity to obtain a
cubic with three distinct rational roots,
then scale the leading coefficient. It is an
elliptic equation with rational \(2\)-torsion.
The transformation uses only a fixed number of
rational operations and the input degrees are
fixed. Both the equation and the mapped rational
points therefore have polynomial height in \(X\),
uniformly in \(P,Q,D_1,s,t\).
Dujella's rational-torsion bounded-height theorem,
as used for MTK1, now gives \(X^{o_X(1)}\) directions.
Each exact direction determines \(x_1,x_2\)
and hence both original primitive rows. Summing
over the chosen parity record gives the bound.
[Dujella, Theorem 1.1 and Corollary 1.3](https://arxiv.org/html/2312.03655v1)
· [Four distinct linear forms (MTK1)](ABC_MULTI_TARGET_MIXED_KERNEL_COUNTS.md#the-mixed-kernel-count)

In particular an outside radical cap \(n^q\) on
the actual product auxiliary in FP1, at fixed
numerical depth, gives original pair exponent
\[
                             q/3+\zeta+\eta .
                                                               \tag{FP2}
\]
For the auxiliary count, fix its complete powers
on the prescribed pool at height-subpower cost.
These powers give arbitrary fixed coefficients.
Use exact-order singleton layers and the
two-cheapest-letter bound MCB2 on the remaining
radicals. Their three normalized masses sum to
at most \(q\) because the primitive auxiliary
roles are pairwise coprime. Half of the two
cheapest costs at most \(q/3\). FP1 gives
numerical height at least \(X^{\nu/2}\), above
the MCB2 threshold at small fixed budgets.
This argument only uses an outside radical cap
and permits prescribed coefficients of
polynomial height. A positive additive ordering
of the auxiliary leaves boundedly many possible
roles for the prescribed \(P,Q\).
[Arbitrary-coefficient two-cheapest count (MCB2)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#bounded-group-dictionaries-with-arbitrary-fixed-coefficients)
· [Effective radical budgets (EB4–EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

### Recovery from the product difference alone

**FP5.** Use the notation of FP1, and exclude the primitive row
\((1,1,2)\) as well as equal or swapped rows. Write
\(Z=|P-Q|>0\). Fix a prescribed polynomial-height prime pool
\(\mathcal P\), and for each pair let \(S\) be the squarefree product
of its shared input primes outside \(\mathcal P\), where a prime is
shared if it divides both \(a_1b_1c_1\) and \(a_2b_2c_2\). If
\(S\le n^\beta\), and the radical of \(Z\) outside the input primes
and \(\mathcal P\) is at most \(n^s\), then the pair count is
\[
                       O(n^{\beta+s}X^\varepsilon).
\]
In particular an exact prescribed \(Z\) has multiplicity at most
\(O(n^\beta X^\varepsilon)\). There is no shared full-divisor or
numerical-height requirement in this recovery bound.

Indeed \(h+j=u(a_1-b_1)\) and \(h-j=v(b_2-a_2)\) are both nonzero.
With \(d,G\) from FP1, fix their exact values at subpower cost
given \(S,\mathcal P\). Knowing \(Z\), the identity \(|hj|=d Z\)
enumerates signed \(h,j\) at divisor cost. Now \(u\mid h+j\) and
\(v\mid h-j\), so these scaling factors also have at most
divisor-many choices. Then \(c_1=vG,\ c_2=uG\) and the two
summand differences recover the four summands by integer halves.
All enumerated integers have bounded polynomial height in \(X\).
Every prime private to either original row divides exactly one of
\(AB,CD\) and is prime to \(d\), hence is a unit in \(Z\). For a
prime on an original summand this follows directly from the scaling
definitions and primitivity. A private prime on \(c_1\) divides \(v\)
and is prime to \(ua_1b_1\), with the symmetric argument for \(c_2\).
The old prime support of \(Z\) is therefore contained in \(S\) and
\(\mathcal P\). Enumerate the remaining radical at cost \(n^s\)
and fix its full powers together with those on \(S,\mathcal P\)
at subpower cost. The union of old support and normalization
records is paid once.

Here is the implication for the FP2 budget. Suppose the two rows'
original radicals outside \(\mathcal P\) each have price at least
\(1-\gamma\). If \(S\le n^\beta\), the private original primes
supply at least \(2-2\gamma-2\beta\) to the auxiliary radical.
Indeed each lies on exactly one of \(P,Q\) after normalization.
If the whole primitive auxiliary has outside radical at most
\(n^q\), its fresh difference therefore has price at most
\(q-2+2\gamma+2\beta\), and the selected pair exponent is bounded by
\[
                       q-2+2\gamma+3\beta
\]
up to arbitrarily small allowances. In a saturation-size sample
of \(n^{2/3-\rho}\) such pairs this forces
\(q\ge 8/3-\rho-2\gamma-3\beta\).
For generic pairs near balance, EB5 lets \(\gamma,\beta\) be
arbitrarily small after fixing the prescribed pool and then the
budget, with height tending to infinity at that budget. For
\(\gamma\) use the cap on pool mass and the fact that low total raw
mass itself saves by EB4. For \(\beta\), apply EB5 with one input
as an anchor, then sum over anchors. Thus on the exponent-full
generic sector the FP2 auxiliary needs \(q\ge8/3-o(1)\). This
lower bound rules out applying the savings test to typical
low-overlap pairs. A separate selection with substantial shared
support or a different auxiliary is needed.

The argument gives an explicit inequality for any selected pairs
with the stated raw lower bounds. An arbitrary power-thin selection
need not inherit the generic low-overlap property.
[Effective radical budgets (EB4–EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

### Pricing shared support by its actual pair incidence

**FP6.** Let a family \(\mathcal F\) have at most \(n^d X^\varepsilon\)
rows in \(X/2<c\le X\), with original per-role radical caps
\(n^{1/3+\tau}\). Fix a polynomial-height prime pool \(\mathcal P\),
and restrict to a set of FP5 pairs of exponent \(D\).
The outside raw radical of each original row is at least \(1-\gamma\).
Classify pairs into fixed small exponent bins using the following
actual quantities:

- \(m\) is the shared radical price outside \(\mathcal P\), and
  \(m_{i,J}\) is the portion on role \(J\) in row \(i\).
  Thus \(\sum_J m_{i,J}=m\) for both rows, even for primes that change roles.
  Put \(v_i=m-\min_Jm_{i,J}\).
- \(e\le m\) is the price of shared primes absent from both normalized
  operands \(P,Q\). Equivalently the valuations in the two raw scaled
  products \(AB,CD\) are equal there, including equality at zero.
- \(q\) caps the whole outside radical of \(P Q|P-Q|\), and \(b\) is
  the actual price of the joint squarefree shared-support records \(S\)
  used in FP5.

One may retain exponent \(D\) and take
\(b\le\min(m,2d-D)+o(1)\). Indeed each row admits at most subpower
many shared records, since \(S\mid\operatorname{rad}(a b c)\).
A record occurring on \(I\) distinct rows has at most \(I^2\) pairs.
For any fixed small \(\eta>0\), records with
\(I<n^{D-d-\eta}\) account for at most \(n^{D-\eta}X^{o(1)}\)
pairs, by summing incidences. The remaining record count is at most
\(n^{2d-D+\eta}X^{o(1)}\). If \(D-d-\eta\le0\), the count is at most
the full incidence count and gives the same assertion. The
\(o(1)\) denotes fixed arbitrarily small exponent allowances after
stratification and the height-first limiting use.

On the retained set the following bounds hold up to these allowances:
\[
 \begin{aligned}
 D&\le d+1/3+\tau-\tfrac12\max(v_1,v_2),\\
 D&\le 2/3+2\tau-\tfrac12(v_1+v_2)+b,\\
 q&\ge 2-2\gamma-m-e+D-b .
 \end{aligned}
\]
For the first inequality anchor one row. The relevant pool is
its entire prime support together with \(\mathcal P\). Fix the
supported factors of the other row at subpower cost and apply
MCB2 to its two cheapest outside caps. Their sum is at most
\(2/3+2\tau-v_i\). Choose the better anchor direction.
For the second inequality fix the common \(S\) itself and apply
the same argument to each side separately. The original shell
is the numerical height floor in each application, since the fixed
parts enter as coefficients.

For the third inequality, the raw union of the two outside input
radicals has price at least \(2-2\gamma-m\). Every private prime
remains on one of \(P,Q\), by FP5. A shared prime remains on an
operand exactly when its two scaled valuations differ. Thus the
operand radical alone has price at least \(2-2\gamma-m-e\).
In addition, the fresh portion \(s\) on \(|P-Q|\) has \(D\le b+s\)
by the FP5 proof using the actual dictionary of \(S\). The operand
support and the fresh portion are disjoint. This proves the bound.
No coloring of a power-thin selected pair family is required.
[Coefficient-uniform count (MCB2)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#bounded-group-dictionaries-with-arbitrary-fixed-coefficients)

In particular in saturation notation \(d=1/3-o(1)\),
\(D=2/3-\rho-o(1)\), \(\gamma,\tau=o(1)\), the raw FP2 savings
test requires \(q+3\zeta<3D\) by a fixed margin.
Since \(b\le\rho+o(1)\), the inequalities give the necessary
conditions
\[
   m+e>2/3+\rho+3\zeta-o(1),\qquad
   m>1/3+\rho/2-o(1),\qquad m\le3\rho+o(1).
\]
Therefore the test fails for every such selection with
\(\rho\le2/15\) in the limiting inequality. More generally,
\(m-\min_J m_{i,J}\le2\rho+o(1)\) on each original row.
For \(\rho<2/9\) by a fixed margin, a viable pair sector must
have all three roles of both rows contributing shared mass,
with \(\min_Jm_{i,J}>1/3-3\rho/2-o(1)\).
This necessary regime uses the raw budget of one actual primitive
auxiliary. A different auxiliary or a stronger count exploiting its
separate role prices has its own test.

### Exact parity of the missing shared primes

The normalization forces a further restriction in FP6. Let \(t\)
be the shared outside target-to-target radical price, allowing
different exponents on the two targets. Let \(e_{\rm even}\)
and \(e_{\rm odd}\) be the prices of shared primes on summands
in both rows which vanish from the entire normalized auxiliary
and have equal original valuations of the respective parity.
Any missing shared prime belongs to one of these two sets or
to the target-to-target support. Indeed a prime on one
target and on an opposite summand lies on just one of the two
scaled products, so survives normalization. At a prime on
summands on both sides the scaled product has the original
valuation, and only equal valuations can both disappear.

The odd missing primes are exactly the common parity support
of the two original products outside the auxiliary and the pool.
In the raw parity cost convention one must therefore use
\(\zeta\ge e_{\rm odd}\). The actual original union contributes
at least \(2-2\gamma-m-t-e_{\rm even}-e_{\rm odd}\) to the
auxiliary outside support. The fresh part of \(Z\) contributes
at least \(D-b\), as in FP6. Moreover
\(v_i\ge(m+t)/2\), since the target overlap portion on each
row is at least \(t\), and the smaller of the two operand
portions is no more than their average. Thus \(m+t\le4\rho+o(1)\)
for the near-saturation pair exponent \(D=2/3-\rho-o(1)\).
Combining with \(b\le\rho+o(1)\) gives
\[
 q+3\zeta-3D\ \ge\
       2/3-3\rho-e_{\rm even}+2e_{\rm odd}-o(1).
\]
A raw FP2 saving with \(\rho\le2/9\) therefore needs a positive
mass of exactly matching even valuation on original summands
that disappears completely from the auxiliary. Near the first
possible boundary \(\rho=2/15\), its price must be at least
\(4/15-o(1)\). If the parity choices are independently compressed
by an actual dictionary, keep \(3\zeta-e_{\rm odd}\) instead
of \(2e_{\rm odd}\) in the displayed lower bound.

If putatively viable raw-budget strata have
\(\rho\to2/15\), the inequalities force \(m,e\to2/5\).
Indeed \(m+e\ge2/3+\rho-o(1)\), \(e\le m\), and \(m\le3\rho+o(1)\).
Every original role on each side consequently has shared price
\(2/15+o(1)\), by \(m-\min_Jm_{i,J}\le2\rho+o(1)\).
The equality \(e=m-o(1)\) forces essentially all shared target
mass to be target-to-target with equal valuations. The parity
inequality then forces
\(e_{\rm even}=4/15-o(1)\). Thus essentially every shared
summand prime carries the same even valuation in both rows
and vanishes from the auxiliary. This describes a necessary
limiting valuation profile for the earliest remaining raw
FP2 regime, without claiming such a regime occurs.

### An opposite-shift terminal for compatible pairs

**FP7.** Suppose the selected pairs have consistent prime roles:
no prime dividing one letter of the first row divides a different
role in the second row. The following estimate applies directly
to a family after a global role coloring, although coloring need
not preserve an arbitrary preselected pair sector.

Then \(d=\gcd(a_1,a_2)\gcd(b_1,b_2)\) divides \(h\), so \(j\mid Z\).
For exact \(Z\), choose the nonzero signed divisor \(j\).
Both opposite pairs \((A,D)\), \((C,B)\) are primitive, since
\(u,v\) are coprime and cross-role primes are excluded.
They have the same difference \(j=A-D=C-B\). One of the two
pairs has a term at least \(X/4\) by the original height shell
and \(u,v\ge1\). In that pair count both endpoints by MCB2 for
the prescribed shift \(j\), keeping the two actual outside
radical caps after fixing powers in \(\mathcal P\) at subpower
cost. Each endpoint pair determines divisor-many \(d\),
since \(d\mid AD\) or \(d\mid CB\). Then \(|hj|=dZ\) gives
a bounded choice of \(h\). The two known endpoints and the
two signed differences determine all four scaled endpoints.
The primitive rows follow by division by the respective gcds.

In particular, for equal targets the scaling factors are one,
and every fixed \(Z\) has fibre at most
\(n^{1/6+O(\tau)+\eta}X^\varepsilon\), uniformly without a
record for \(G\), \(d\) or the common parity. If the shared
outside summand mass is \(x=m_A+m_B\), and each original row
has outside summand radical mass at least \(2/3-o(1)\), the two
operands retain outside mass at least \(4/3-2x-o(1)\). The
summand lower bounds follow from the original total lower caps
and target upper caps in FP6. The whole
outside \(Z\) cap is at least \(D-1/6-o(1)\) by the fibre bound.
These two contributions add because \(P,Q,Z\) are pairwise
coprime. On the other hand, the count of second rows per fixed
target and first row is at most
\(n^{1/6-\max(m_A,m_B)/2+o(1)}\) by FL1 with the anchor pool.
Hence \(D\le1/2-\max(m_A,m_B)/2+o(1)\), and
\(x\le4\rho-2/3+o(1)\) when \(D=2/3-\rho-o(1)\).
It follows that
\[
             q-3D\ \ge\ 7/6-6\rho-o(1).
\]
Thus the raw FP2 budget fails for compatible equal-target
selections through \(\rho=7/36\) in the limiting inequality,
even if the missing parity record is free. In particular this
applies at the full equal-target pair count for a regular
balanced dictionary of exponent \(\alpha\le7/36\), by forming
the equal-target pairs after global role coloring.
[Fixed-letter count (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

### Linear auxiliary heights and its three separate radical caps

**FP8.** In the original height shell the FP1 product auxiliary
has a numerical height floor without any shared-lock hypothesis.
For a family of \(N\) original rows, all but \(O_\varepsilon(NX^\varepsilon)\)
ordered pairs satisfy
\[
                           Q>c_1v/8,\qquad P>c_2u/8 .
\]
To prove the first bound fix the first row, choose \(G\mid c_1\),
\(d\mid a_1b_1\) and the factorization \(d=d_a d_b\) according to
the primes on the two summands of the second row. These have
divisor-many possibilities uniformly at polynomial height.
Write \(a_2=d_a z_a,\ b_2=d_b z_b\). Primitivity gives
\(\gcd(z_a,z_b)=1\) and makes \(d_a,d_b,z_a,z_b\) units
modulo \(G\), because \(G\mid c_2\). Now \(Q=v^2 z_a z_b\).
If \(Q\le c_1 v/8=v^2G/8\), dyadically bin \(z_a,z_b\).
In any occupied bin the cross-products of two possible pairs are
at most \(4(G/8)\). Since \(d_a z_a+d_b z_b\equiv0\bmod G\),
their reduced ratios coincide modulo \(G\). Their cross determinant
has magnitude less than \(G\), so it vanishes and the reduced
pairs are equal. This gives at most one second row per bin and
proves the exception bound. Reverse the anchor for the bound on
\(P\). In particular exponent-rich selections with pair exponent
\(D>d\) have auxiliary height at least \(X/16\) after removing a
power-smaller set.

Keep the three separate outside radical masses \(p_1,p_2,p_3\) of
\(P,Q,|P-Q|\) in bins and let \(q\ge p_1+p_2+p_3\) cap the total,
up to the bin tolerance. FP2 followed by MCB2 with the proved
height floor gives
\[
             D\le \zeta+\tfrac12
               \bigl(p_1+p_2+p_3-\max_i p_i\bigr)+o(1).
\]
Thus any viable use of the raw test \(q+3\zeta\le 3D+o(1)\)
requires \(p_i\le D-\zeta+o(1)\) for all three roles. In the
saturation setting of FP6, the old support left on the two
operands is exactly its union size less the erased mass \(e\).
Consequently
\[
               p_1+p_2\ \ge\ 2-m-e-o(1)
                            \ \ge\ 2-6\rho-o(1).
\]
Here \(m\le3\rho+o(1)\) is the anchored FP6 bound and we used
\(e\le m\). Since \(D=2/3-\rho-o(1)\), it follows that
\[
                              \rho\ge1/6+\zeta/2-o(1).
\]
In particular the raw FP2 test fails by a fixed margin throughout
\(\rho<1/6\), and \(\rho=1/6\) is also insufficient for a strict
saving in the limiting inequalities. The direct test with the
separate auxiliary caps can be stronger than the raw test.

For a sequence of putatively viable raw strata with
\(\rho\to1/6\), these bounds give \(\zeta\to0\),
\(m,e\to1/2\), every original overlap role mass tending to
\(1/6\), and all three auxiliary radical masses tending to
\(1/2\). Indeed the operand lower bound reaches one, and the
two-cheapest bound forces each auxiliary mass to reach \(1/2\).
The anchored bound forces each role's overlap to reach \(1/6\).
The identity \(m-e=o(1)\) precludes positive target-to-summand
overlap or unequal valuations among the shared target primes.
For raw missing parity, the exact parity inequality also
requires a mass at least \(1/6-o(1)\) of matching even original
summand primes deleted from the whole auxiliary. This is a
necessary limiting profile at the new raw boundary.

### Fixed auxiliary difference without compatible prime roles

**FP9.** At equal original targets \(u=v=1\), the fixed-\(Z\)
terminal can handle arbitrary prime-role overlap. Let \(g_{ij}\)
be the four cross-gcds of the two original summands, and \(p_{ij}\)
their outside radical masses. These gcds are pairwise coprime
and their product is \(d\). Set \(d_h=g_{AA}g_{BB}\) and
\(d_j=g_{AB}g_{BA}\). Then \(d_h\mid h,\ d_j\mid j\), and
\(Z=|(h/d_h)(j/d_j)|\).

For a prescribed \(Z\), choose the two signed quotient factors
at divisor cost. To use the opposite-role shift, enumerate
\(g_{AB},g_{BA}\) at raw price \(p_{AB}+p_{BA}\). This fixes
\(j\). The normalized pairs \((A,D)/g_{AB}\) and
\((C,B)/g_{BA}\) are primitive, with shifts divisible by
\(g_{BA}\) and \(g_{AB}\), respectively. Write \(M_1,M_2\)
for their largest terms. Then \(M_1\ge g_{BA}, M_2\ge g_{AB}\)
and \(M_1M_2\ge\max(A,B,C,D)\ge X/4\). Apply MCB2 to the tall
one, after fixing its complete factors on the corresponding
recorded gcd support and \(\mathcal P\). This fixing also strips
any residual prime powers left after unequal valuations in the
gcd. The resulting exponent is at most
\(1/6-\min(p_{AB},p_{BA})/2+o(1)\).
The fixed tall endpoint pair allows only divisor-many remaining
\(d_h\), since \(d_h\mid AD\) or \(d_h\mid CB\). Recover \(h\)
and the four endpoints. The same-role arm uses \(g_{AA},g_{BB}\)
and \(h\), with \(d_j\mid AC\) or \(d_j\mid DB\).
Thus the fixed \(Z\) fibre exponent is at most \(1/6+F+o(1)\), where
\[
 F=\min\left\{p_{AA}+p_{BB}-\tfrac12\min(p_{AA},p_{BB}),\
              p_{AB}+p_{BA}-\tfrac12\min(p_{AB},p_{BA})\right\}.
\]
The two arms also apply with moving targets if the count of each
normalized fixed shift is priced by its actual scaled endpoint
caps. The four scaled endpoint gcds equal the original cross
gcds, since \((u,v)=1\), \((u,a_2b_2)=1\) and \((v,a_1b_1)=1\).

In the equal-target case, on an exponent-\(D=2/3-\rho-o(1)\)
stratum with outside summand masses \(1/3+o(1)\), FL1 on anchor
fibres bounds every row and column sum of \((p_{ij})\) by
\(s=2\rho-1/3+o(1)\). Write \(x=\sum p_{ij}\).
The expression \(2x+F\) increases with the entries.
Complete the matrix to all marginals \(s\); its diagonal now
has equal entries \(a,a\) and its off-diagonal entries are
\(s-a,s-a\). Thus \(2x+F\le4s+3s/4\).
The operands contribute at least \(4/3-2x-o(1)\) private mass,
and the outside radical of \(Z\) is at least
\(D-1/6-F-o(1)\) by the fibre bound. By adding these disjoint
masses,
\[
                              q-3D\ge17/12-(15/2)\rho-o(1).
\]
So the raw test fails for arbitrary equal-target selections
through \(\rho=17/90\) in the limiting inequalities. FP7
gives the stronger \(7/36\) threshold on compatible pairs.

### A fixed numerical margin for both operands and the difference

**FP10.** Fix a raw cap \(C_0\) and an allowance \(\eta>0\).
There are \(\theta,L_0>0\) such that for \(0<L<L_0\) fixed,
at most \(O_{\eta,C_0,L,\varepsilon}(n^\eta X^\varepsilon)\)
second rows per anchor have
\[
              Q\le v^2G X^\theta
          \quad\hbox{or}\quad Z=|P-Q|\le X^\theta.
\]
All three original caps are at most \(n^{C_0}\), both original
targets belong to \(X/2<c_i\le X\), and the excluded repeated,
swapped and trivial rows are as in FP5. The bound for \(P\)
follows by reversing the anchors. No compatible-prime condition
is needed. The constants are uniform in the anchors and exact
gcd records.

Fix \(G,d,d_a,d_b\) as in FP8 and put \(a_2=d_a z_a,b_2=d_b z_b\).
If \(Q\le v^2 G X^\theta\), dyadic upper endpoints \(M_a,M_b\)
for \(z_a,z_b\) satisfy \(M_aM_b\le4GX^\theta\).
The lattice \(d_a z_a+d_b z_b\equiv0\bmod G\) has index \(G\).
Scale the coordinates by \(M_a,M_b\), so the data are in the
unit square and the resulting lattice has determinant
\(G/(M_aM_b)\). If the points have linear rank at most one,
coprimality gives at most one pair. Otherwise take a shortest
nonzero lattice vector and lift the shortest positive vector
of the lattice projected onto its orthogonal complement.
Reduce the parallel component of that lift by subtracting a
multiple of the shortest vector. The two basis vectors are
bounded in length because of the two independent witnesses in
the unit square. Cramer's rule now gives integer coordinates
of size \(O(M_aM_b/G)\). On these features all three original
letters are exact linear functions. If \(G\ge X^{1/4}\) there
is a fixed tall divisor on the target. Otherwise
\(z_a,z_b\le4X^{1/4+\theta}\); since the larger original summand
is at least \(X/4\), one of \(d_a,d_b\) is at least \(X^{1/4}\)
for sufficiently small \(\theta\). Apply the
[exact-feature lemma with one fixed tall divisor](ABC_AFFORDABLE_FIBRE_AFFINE_RANK.md#exact-small-features-with-one-fixed-tall-divisor).
The dyadic and divisor choices are height-subpower.

For the short difference fix the four cross-gcds of FP9.
The signed integers \(H=h/d_h\), \(J=j/d_j\) both have
magnitude at most \(Z\), since \(Z=|HJ|\) and neither is zero.
Put \(s_1=a_1-b_1\ne0\). The identities
\[
 u=\frac{d_hH+d_jJ}{s_1},\qquad
 (a_2,b_2)=\frac{1}{v s_1}
       (b_1d_hH+a_1d_jJ,\ a_1d_hH+b_1d_jJ)
\]
give exact linear features for the second row. Also
\(c_2=Gu\le 2GdZ\), so one among the three fixed divisors
\(G,d_a,d_b\) is at least \(X^{1/4}\) for small enough
\(\theta\). The same lemma applies.

There is a stronger relative difference statement on a
fixed-heavy-divisor sector. If at least one of \(G,d_a,d_b\)
is at least \(X^\nu\) on the given record for fixed \(\nu>0\),
the same rate bound holds for
\[
                            Z\le v|s_1| X^\theta .
\]
Indeed the linear map from \(a_2,b_2\) to \(h,j\) has matrix
\[
                \frac1G\begin{pmatrix}-b_1&a_1\\a_1&-b_1\end{pmatrix}.
\]
Apply it to \(a_2=d_a z_a,\ b_2=d_b z_b\) on the index-\(G\)
lattice above, then divide the two coordinates by \(d_h,d_j\).
This maps the whole lattice to integer \(H,J\), by the four
fixed gcd divisibilities and the definition of \(u\).
Since \(d_ad_b=d_hd_j\), its image has index \(v|s_1|\)
in the \(H,J\) lattice. Dyadic upper endpoints for \(|H|,|J|\)
have product at most four times the bound on \(Z\).
The same bounded-basis argument gives two small features.
On rank one there is at most one primitive point per positive
original ray, by invertibility of the linear map. This proves
the relative claim with \(\theta\) chosen in terms of \(\nu,\eta\).

In particular any selected pair exponent exceeding the original
row exponent by a fixed amount permits discarding both
small-operand branches and the absolute small-difference branch.
The surviving pairs have
\[
 \log_X d+\log_X G
       \le \min_i(\log_X a_i+\log_X b_i)-\theta,\qquad Z>X^\theta .
\]
This numerical margin holds independently of the separate
auxiliary radical budgets. The relative statement can be
imposed in addition on sectors with a fixed heavy gcd factor.

### Recovering matched square factors from an actual target modulus

**FP11.** Keep \(H,J,d_h,d_j\) as above. Write
\(d_h=\kappa_h x^2,\ d_j=\kappa_j y^2\) with squarefree
positive kernels and \((x,y)=1\).
Give numerical caps \(d_h\le X^{w_h},\ d_j\le X^{w_j}\).
Suppose there are divisors \(D_u\mid u,D_v\mid v\)
of joint price \(b\) satisfying
\[
                              D_uD_v>2X^{(w_h+w_j)/2}.
\]
Here a price may be the raw radical bound or any uniform
conditional catalogue bound for the actual witnesses.
For fixed \(Z\), assume the squarefree kernels \(\kappa_h,\kappa_j\)
are priced at \(\sigma\) by their new support outside
the primes of \(Z\) and a prescribed pool. Then \(d_h,d_j,h,j,u,v\)
and the four summand cross-gcds have a joint catalogue of price
at most \(b+\sigma\), with height-subpower factors.

Indeed choose signed \(H,J\) as divisors of \(Z\), fix the two
kernels and the two actual locks at their stated costs, and use
\[
 (x/y)^2\equiv-\frac{\kappa_jJ}{\kappa_h H}\pmod {D_u},
 \qquad
 (x/y)^2\equiv+\frac{\kappa_jJ}{\kappa_h H}\pmod {D_v}.
\]
The normalizer is prime to \(u,v\). At a prime dividing \(u\)
the two shifts reduce to \(-v a_2,-v b_2\), both units, with the
corresponding assertion at \(v\). All divisions in the two
congruences are therefore legitimate. Unit square roots cost
\(X^{o_X(1)}\) by prime-power root counting. In one CRT class
the cross determinant of any two candidate fractions is divisible
by \(D_uD_v\) and smaller than it, so the two reduced fractions
agree. This fixes \(d_h,d_j,h,j\), and \(u,v\) are divisors of
the nonzero \(h+j,h-j\), respectively. Split \(d_h,d_j\) into
the four cross-gcds at divisor cost.

To recover the original rows, the only further parameter is \(G\).
Write \(F_a=g_{aa}g_{ab}, F_b=g_{ba}g_{bb}\). The equations
\(2a_1=vG+s_1,\ 2b_1=vG-s_1\) determine \(G\) modulo \(d\),
since \(v\) is a unit there and the two moduli \(F_a,F_b\)
are coprime. If \(G<d\) in the selected bin there is at most
one choice, and the fixed-\(Z\) fibre bound is
\(O(n^{b+\sigma}X^\varepsilon)\).
With no numerical condition on \(G/d\), one may instead add
the half-price of the cheapest available pair of residual
summands. On the first row fix all powers on the primes
of \(F_a,F_b\) and of the prescribed pools at subpower cost.
Then MCB2 applies to the fixed shift \(s_1=a_1-b_1\) with
the two residual summand caps and height from the original shell.
The second row gives an independent version of the same bound.
For instance at the tied raw FP8 profile, these residual caps
are \(1/6+o(1)\) on each summand, and the additional fibre
price is at most \(1/12+o(1)\). No numerical dominance of
\(d\) over \(G\) is needed for that weighted terminal.

If \(G/d\le X^\theta\) and at least one of \(F_a,F_b\)
is at least \(X^\nu\), the possible \(G\)'s instead cost
\(n^\eta X^\varepsilon\) for any fixed \(\eta>0\) at sufficiently
small \(\theta,L\), by the exact-feature lemma. Fix an occupied
value of \(G\) as origin for the affine quotient coordinate.
It is of size \(O(X^\theta)\), and all three original first-row
letters are exact affine functions of it.

At a tied raw FP8 profile, \(\sigma=o(1)\) follows under the
raw missing-parity hypothesis. The odd primes retained on
\(Z\) are in its known support, the hidden odd primes have
raw cost \(\zeta=o(1)\), and unequal shared valuations have
raw cost \(m-e=o(1)\). A parity catalogue conditional only
on \(P,Q\) does not supply this \(Z\)-conditional record.

### The actual shared target record at a tied pair boundary

At an exponent-full tied raw FP8 boundary the conditional
\(G\) record after fixing \(d\) has arbitrarily small
normalized price on all but a power-smaller set of pairs.
Indeed FP6's incidence bound lets us use an actual
\((G,d)\) catalogue of exponent at most \(1/6+o(1)\).
One can see this directly by treating the exact \((G,d)\)
as the shared record: at each incident original row both
entries divide its complete product, so there are only
height-subpower many records incident to that row.
The rich-record argument keeps the full pair exponent
and gives the stated catalogue.

Any fixed exact \(d\) has matched residual raw caps at most
\(1/6+o(1)\) on each original summand after fixing the
complete supported parts on \(d\)'s primes at subpower cost.
This follows from \(m-e=o(1)\), the matched summand masses
\(1/6+o(1)\), and the fact that shared target primes are
distinct from the product gcd on both sides.
By MCB2 there are at most \(n^{1/6+o(1)}\) such first
rows. On the selected near-tied pairs each first row
has degree at most \(n^{1/6+o(1)}\) by FP6 anchor stripping.
For any fixed small \(\eta>0\), the \(d\)'s associated with
more than \(n^\eta\) choices of \(G\) have catalogue exponent
at most \(1/6-\eta+o(1)\). They account for at most
\(n^{1/2-\eta+o(1)}\) pairs. After deletion, fixing \(d\)
leaves a \(G\)-price of at most \(\eta\), independently of
their numerical order.

This permits using the \(b+\sigma\) exponent of the square
ratio record for the fixed-\(Z\) fibre at the tied boundary,
with an arbitrarily small additional allowance. The subsequent
dictionary contradiction still needs actual affordable
target-quotient divisors for the square ratio test.

### A close-target dictionary test and the residue ambiguity

Assume an exponent-full tied FP8 family with \(D=1/2-o(1)\)
and auxiliary radical masses \(1/2+o(1)\) each, for the
ordinary auxiliary radical budget. In particular a free support
pool, if used, must add only vanishing raw mass to the auxiliary.
The original input caps remain in force. The auxiliary height
is at least \(X/16\) outside negligible original pairs and at
most \(X^4\). In dyadic height shells it is therefore an
ordinary balanced saturation family, by FP2 and vanishing
raw parity price. Rebasing the budget at that height preserves
the height-first limiting argument. Its three complete
dictionary exponents are equal by RB3, with conversion back
to the original \(n\) scale. FP11 with no additional \(G\)-price
forces the \(Z\) dictionary exponent to be at least \(1/2-b-o(1)\).
Equivalently both joint product–quotient lists
\((a_1 b_1/d,u)\), \((a_2 b_2/d,v)\) have at least that
exponent, as \(u^2\mid P,\ v^2\mid Q\) give the reverse
recovery at divisor cost.

There is an unconditional dictionary upper bound within this
tied profile. Let \(\alpha_C\) be the original complete target
dictionary exponent on the selected pairs. Fix \(G\) and its
prime pool. Because the matched shared target radical has
price \(1/6-o(1)\) on each original row, MCB2 bounds the
first rows with this \(G\) by \(n^{1/4+o(1)}\).
It uses residual target cost at most \(1/6+o(1)\) and one
summand cap \(1/3+o(1)\), with all powers on \(G\)'s support
fixed at subpower cost and the original numerical height
floor intact. Summed over \(G\), the number of pairs
\((G,c_2)\) is at most \(n^{\alpha_C+o(1)}\):
\(G\) is a divisor of \(c_2\). These two records, the first row,
and \(d\mid a_1b_1\) determine \(u\) and \(P\).
Thus the \(P\) dictionary has exponent at most
\(\alpha_C+1/4+o(1)\). The \(Q\) dictionary has the same bound
by reversal, with the actual target dictionary of the
corresponding side. In particular the tied family is
impossible if the square-recovery price, including any
additional \(G\) price, is below \(1/4-\alpha_C\) by a
fixed margin. This conclusion does not require a separate
cheap private target-quotient catalogue.

A different upper bound on those joint lists comes from close
targets, and is independent of the original target dictionary size. If
\(|c_2/c_1-1|\le X^{-\omega}\) for a fixed \(\omega>0\),
then for a fixed first row and its divisor \(G\), the quotient
\(u\) lies in one relative strip of logarithmic width
\(O(X^{-\omega})\) centered at the known \(v\).
Its actual raw cap is \(1/6+o(1)\), because the matched
shared target primes have been removed, with other overlap
cost \(o(1)\). Apply JRM1 with modulus one, one real bin,
and weights summing to at most one. At depth \(\omega/L\)
larger than twice the total cap, choose \(h_i=\lfloor
(\omega/L)/a_i\rfloor\) on each exact-order slot of cap \(a_i\).
Then the gain \(G_1\) reaches the summed slot cap.
The quotient catalogue per anchor has exponent at most
\(1/12+o(1)\), after the usual arbitrarily small tail allowance.
The first row costs \(1/3+o(1)\) and \(d\) is a divisor of
its product. The \(P\) dictionary thus costs at most \(5/12+o(1)\).
The tied family is impossible on a square-recovery sector
of price \(b<1/12\) by a fixed margin when \(G\) has
subpower additional price in FP11.
[Real bin moments (JRM1)](ABC_JOINT_REAL_RESIDUE_MOMENTS.md#an-arbitrary-list-of-congruence-and-logarithmic-bins)
· [Balanced dictionary equality (RB3)](ABC_BALANCED_DICTIONARY_RESIDUE_CELLS.md)

There are actual primitive pairs showing the last residue
can have unbounded \(G\). Take \(u=2,v=3,d_h=1225,d_j=1,H=J=1\),
and the two rows
\[
 ((3G+613)/2,(3G-613)/2,3G),\qquad (G-204,G+204,2G).
\]
For sufficiently large \(G\) with \(G\equiv29\bmod125\),
\(G\equiv41\bmod343\), and \(G\equiv1\bmod(2\cdot3\cdot17\cdot613)\),
the rows are positive and primitive, the two \(a\)'s have
valuation two at 5 and the two \(b\)'s valuation two at 7.
The exact shifts are \(h=1225,j=1\), so the prescribed
overlaps are exact, and \(Z=1\). Both targets lie in the shell
\(X=3G\) and both target radicals can be bounded along the fixed
progression by taking powers of one admissible \(G_0\).
This example imposes no small radical caps on the summands.

### Two shared role locks cannot be tight at small original dictionaries

At the tied FP8 boundary, suppose the original complete
dictionaries all have exponent at most \(\alpha<1/4\) by a
fixed margin. The two original row projections each have
exponent \(1/3-o(1)\), since on each side the individual anchor
degree is at most \(n^{1/6+o(1)}\) by FP6 and the matched
per-role overlap is \(1/6+o(1)\). Refine the shared support
record by its role assignments and the full powers on each
side at subpower cost per incident row. Each incident row has
only height-subpower many such records, by the divisor bound
and uniform enumeration on a fixed support. The rich-record
argument of FP6 can be applied on both sides: delete records
having fewer than \(n^{1/6-\eta}\) rows on either side, by
summing incidences on the other side. The retained records
are rich on both sides with arbitrarily small \(\eta>0\).

A sector where two common full role locks on either side
have product at least \(X^{1-\theta}\) is power-smaller
than the tied pair population, for sufficiently small fixed
\(\theta>0\). To prove this, fix the choice of the two roles
and the side, and suppose the sector has the full tied
exponent. For every retained exact record the two locks
give at most \(n^\eta X^\varepsilon\) original rows per
fixed third complete letter, by JD5. Thus the record has
at least \(n^{1/6-O(\eta)}\) distinct values of the third
letter \(W\). For each distinct third shared full factor
\(F_W\), choose one rich record containing it. The
\((F_W,W)\) incidences across these representatives count
at most \(n^{\alpha+o(1)}\), since \(F_W\mid W\).
Hence the actual \(F_W\) catalogue costs at most
\(\alpha-1/6+O(\eta)\). But fixing \(F_W\), whose raw price
is \(1/6+o(1)\) in the tied sector, leaves at most
\(n^{1/4+o(1)}\) original rows by MCB2 after stripping
the supported powers on that role. Their projection
could then have exponent at most
\(\alpha+1/12+O(\eta)<1/3\), a contradiction.
All inequalities have fixed margins for
\(\alpha\le1/4-\delta\); the usual box and richness
losses can be chosen smaller than those margins.
[Two deep fixed roles (JD5)](ABC_JOINT_DISPLACEMENT_CATALOGUES.md#tight-full-locks-force-almost-complete-register-information)

### The savings test including the usable pair loss

Suppose a chart of \(N\) original rows supplies at least
\(N^2 n^{-\rho}\) ordered usable pairs, up to fixed small exponent
allowances. The pairs must satisfy FP2's actual auxiliary height,
outside radical \(n^q\), and missing parity cap \(n^\zeta\).
Then the FP2 upper bound and this lower bound give
\[
   N\ll n^{q/6+(\zeta+\rho)/2+\eta}X^\varepsilon .
\]
For a cover by charts of price \(b_0\), the strict small-budget
condition is \(q+3\zeta+3\rho+6b_0<2\) by a fixed margin.
One may also apply FP2 to at least \(n^{2/3-\rho}\) pairs extracted
directly from a saturation-sized family, in which case \(\rho\)
already includes the extraction loss and no separate chart cost
should be counted for an estimate uniform on that selection.

Here \(q\) is the radical of the **primitive** auxiliary after
the prescribed free support and actual cancellation in FP1.
Only unrecorded common parity contributes to \(\zeta\). An almost
injective pair recovery does not imply the auxiliary cap. In
particular a new pair record must either lead to a genuine auxiliary
to which FP2 applies or prove a separate counting bound.

The parity record is cheap on typical pairs.
Suppose all three original cap exponents are at
most \(1/3+\tau\). For any fixed first row, a
second row sharing at least \(\sigma\) of
normalized radical with the first row's
fixed complete product support has count
\(O(n^{1/3+\tau-\sigma/3+\eta}X^\varepsilon)\)
by EB4 and EB5. Thus power-full pairs of a
near-full family can be restricted to
\(\zeta\le\sigma\) provided their pair loss
is below this saving, after the actual
chart and cap tolerances. The omitted common
parity mass can be much smaller than either
individual parity field price. If the pairs
have equal targets and the fixed-target
multiplicity approaches \(1/6\), the same
conclusion follows fibre by fibre: sharing
\(\sigma\) on the two free operands forces
one operand to spend at least \(\sigma/2\)
on the first row's pool, so FL1 bounds such
second rows by \(n^{1/6+\tau/2-\sigma/4+\eta}\).
For a chart of exponent \(d_0\) and pair
exponent \(2d_0-\rho\), the sufficient
deletion margin is
\(\sigma/3>1/3-d_0+\tau+\rho\).
For an equal-target selected exponent
\(d_0+u_0-\rho\), it is
\(\sigma/4>1/6-u_0+\tau/2+\rho\).
The strict margins absorb the height allowances.
[Fixed-letter count with prescribed pools (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

Thus on an exponent-full pair sector in a
balanced chart of exponent \(d_0\) near
\(1/3\), at actual lock price and chosen
pair loss \(\rho\), the cross-product auxiliary
needs outside radical cap at least
\(3(2d_0-\rho-\zeta)-o(1)\) if the inequalities
are to coexist. For full-scale generic pairs
on arbitrarily cheap charts this approaches
the two-row cap \(2\). At equal target and
complete multiplicity \(u_0\), a selected
pair exponent \(d_0+u_0-\rho\) instead
requires \(3(d_0+u_0-\rho-\zeta)\).
These are necessary conditions on the
actual product differences, not assumptions
about the radical of an arbitrary new shift.

## Two other field-free recovery records

**FP3.** If the two original complete targets
are independently prescribed, one cross-product
auxiliary can replace the quartic recovery.
Normalize the pair \((a_1b_2,a_2b_1)\) by its
gcd and add the positive difference as the
third entry. Its exact ratio together with
the prescribed \(c_1,c_2\) recovers both
original rows at divisor cost by TR1. With
moving targets, use both \((a_1b_2,a_2b_1)\)
and \((a_1a_2,b_1b_2)\) as the two
normalized comparisons. Their exact ratios
\(\rho,\sigma\) recover the original odds
from the positive square roots of
\(\rho\sigma\) and \(\sigma/\rho\).
Every such comparison is numerically
tall in a common absolute residue cell:
the cross-products agree modulo \(M\),
and their gcd is a unit there. Their
nonzero reduced difference is therefore
a multiple of \(M\).
[Exact pair recovery from a target ratio (TR1)](ABC_DIFFUSE_FACTOR_GRAPH.md#recovery-from-a-matched-target-ratio)

A different second auxiliary is the primitive
comparison of two scaled endpoints, for
example \(A,C,|h|\) divided by
\(\gcd(A,C)\). Its ratio is the \(r\) used in
FP2, so fixing this comparison and the
product auxiliary recovers the two rows
directly without a parity record. An
opposite-role comparison uses \(A,D\) and
then interchanges \(a_2,b_2\) during the
recovery. An absolute same-residue lock
suffices for numerical depth even when
the scale factors are not units mod \(M\).
For example, for normalized endpoints
\(S=A/g_-,T=C/g_-\) with \(g_-=(A,C)\)
one has \(u\mid S,v\mid T\), and
\(M\mid vS-uT\). If \(a_1\ne a_2\),
this nonzero difference is at most \(ST\)
in absolute value. Use the other
same-role comparison if \(a_1=a_2\)
and \(b_1\ne b_2\).
For opposite roles, use instead the
positive \(vS+uT\) divisible by \(M\)
and at most \(2ST\).
The two ratios here can be recovered
from their exact auxiliary additive
triples up to a bounded choice of order.
The product auxiliary must retain the
height floor of FP1. Under the absolute
lock handle equal original operand products
at the separate divisor bound of FP1.
Apply the same MCB2 count used in FP2 to
both actual auxiliary triples with their
own caps and prescribed pools, giving
one-third of the sum of the two
actual auxiliary radical prices.
This alternative trades the parity
overlap price in FP2 for an additional
ordinary radical price. Both recoveries
allow all original caps and genuine
least-repair eligibility.

## A boundary for three-term identities of the same records

**FP4.** In three independent variables
\(B,h,j\), let the seven linear forms be
\[
 B,\quad A=B+h+j,\quad C=B+j,\quad D=B+h,\quad
 h,\quad j,\quad W=2B+h+j .
\]
Among homogeneous polynomial identities given by
a nonzero linear combination of three monomials in
these records, with no common factor and using
both \(h,j\), the only one is
\(AB-CD=-hj\), up to scaling and reordering.
All monomials have nonnegative integer exponents.

Indeed a factor in two monomials must divide
the third, so the supports are pairwise
disjoint, by irreducibility of the seven
distinct linear forms. If \(h,j\) are in
different terms, reduce modulo \(h\).
The remaining term with \(j\) cannot equal
the one without it: at \(h=0\) none of the
other remaining forms is proportional to \(j\).
Hence both are in one term. Reducing at
\(h=0\) requires equal exponent totals in
the other two terms on each class
\(\{A,C\},\{B,D\},\{W\}\). Reducing at
\(j=0\) gives equal totals on
\(\{A,D\},\{B,C\},\{W\}\). Their exponent
difference is consequently
\((k,k,-k,-k,0)\) in the order \(A,B,C,D,W\).
Using disjoint supports and interchanging
terms if needed, the two monomials are
\((AB)^k\) and \((CD)^k\) for \(k>0\),
with opposite coefficients. Their difference is
\[
         -h j\sum_{r=0}^{k-1}(AB)^r(CD)^{k-1-r}.
\]
The sum is nonzero on both \(h=0\) and \(j=0\),
so the two shift exponents of the third
monomial must both be one. Its only other
available factor is \(W\), to exponent
\(2k-2\). At \(h=0\) the sum is
\(k[B(B+j)]^{k-1}\), which is proportional to
\((2B+j)^{2k-2}\) only for \(k=1\).
This classifies the proposed monomial
transfers. The common-field square identities
use square roots of the cross-products, an
additional algebraic structure outside this
record algebra.
[Common-field auxiliary triples (CP4)](ABC_CROSS_PHASE_AND_AUXILIARY_TRIPLES.md#two-auxiliary-abc-triples-in-the-common-field-equal-target-case)

### A deep target divisor on the reduced determinant ratio

**FP12.** Keep the notation and primitive shell assumptions of FP1
and FP5, including the exclusions there. Write \(d_h=g_{aa}g_{bb}\)
and \(d_j=g_{ab}g_{ba}\) for the diagonal and off-diagonal summand
gcd products, as in FP11. Put \(H=h/d_h,\ J=j/d_j\) and \(Z=|HJ|\).
Take the height limit at each fixed \(L>0\), with \(n=X^L\).

Fix one exact common shared-support record, including the role assignments and all exact original powers on either side that are supported on its primes or on any additionally prescribed polynomial-height pool. The fixed complete numbers include \(G,d_h,d_j\). At fixed support, the refinements into complete powers, choices of roles and exact \(G,d_h,d_j\) all have uniform height-subpower cost. Write \(\mathcal P\) for the additional prescribed pool and \(S\) for the shared input primes outside it. The fixed record sets these supports on both rows. Let \(r_Z\) be a cap on \(\log_n\operatorname{rad}_{\mathcal P^c}(Z)\), and \(\ell_{\mathrm{old}}\) a lower bound on the raw of the \(S\)-supported part of \(Z\). All other primes of \(Z\) outside \(\mathcal P\) are fresh, by FP5.

Assume this selected pair family has, per row pair, a divisor \(M\mid u\) or \(M\mid v\), with \(M\ge X^\nu\) for fixed \(\nu>0\), whose actual conditional dictionary price within the fixed record is at most \(\mu\). In particular this holds when \(\operatorname{rad}(M)\le n^\mu\); enumerate the radical then full powers at height-subpower cost. The lock need not be an exact full part, although full divisors of private positive scaled-degree bands are the principal application. All original letter heights and all supported inputs are bounded by fixed powers of \(X\). Keep the three original radical caps.

For every fixed normalized allowance \(\eta>0\), at all sufficiently
small fixed \(L\), the pair count per record is
\[
 O\!\left(n^{(r_Z-\ell_{\mathrm{old}})/2+\mu+\eta} X^\varepsilon\right).
\]
The constant is uniform for fixed original cap bounds, pool height bounds,
numerical depth, budget and allowances. The same bound with \(\mu=0\)
holds without \(M\) if the original summand imbalances have a fixed
positive relative numerical separation, as described below.

#### Exact ratio recovery and the remaining gcd charge

Put \(g=\gcd(|H|,|J|)\) and \(H'=H/g,J'=J/g\). These are signed coprime nonzero integers. At every prime \(p\mid u\), the scaled differences \(h,j\) are congruent to \(-v a_2,-v b_2\) modulo \(p\), both units. At every prime \(p\mid v\), both are congruent to \(u a_1\) modulo \(p\), also a unit. In particular \(d_h,d_j,H,J,g,H',J'\) are units modulo the chosen \(M\). The identities above give the exact unit congruence
\[
        H'/J'\equiv -d_j/d_h\pmod M\quad(M\mid u),
        \qquad
        H'/J'\equiv d_j/d_h\pmod M\quad(M\mid v).
\]

Let \(\tau\) denote the raw mass outside \(S\cup\mathcal P\) on primes of \(g\) whose exact valuations on \(H\) and \(J\) coincide. These and only these gcd primes disappear from both \(H'\) and \(J'\). Hence
\[
 \log_n\operatorname{rad}_{(S\cup\mathcal P)^c}(H'J')
                 \le r_Z-\ell_{\mathrm{old}}-\tau.
\]
Fix the lock at price \(\mu\). Use RX7 with one residue, the two
positive magnitudes, weights \(1,-1\), and the known pool
\(S\cup\mathcal P\). Record the signs and round the individual
outside caps on a fixed fine mesh. The lower depth bound on \(M\)
and the polynomial height bound on the magnitudes give the
exact ratio count with exponent
\((r_Z-\ell_{\mathrm{old}}-\tau)/2+\eta\).

There is also a real version. If
\(\lvert a_1-b_1\rvert\le X^{-\nu}\lvert a_2-b_2\rvert\), then
\(\lvert h+j\rvert\ll X^{-\nu}\lvert h-j\rvert\), since
\(u/v=c_2/c_1\) is bounded by the original height shell.
The signed ratio \(H'/J'\) lies in a relative logarithmic interval
of width \(O(X^{-\nu})\) around \(-d_j/d_h\). The reverse
imbalance gives the centre \(d_j/d_h\). Use JRM1 with modulus one
and one real interval, after fixing the signs, the powers on the
known pool and the high-order tails at arbitrarily small cost.
For root-slot caps \(a_i>0\) of total \(r\), the moment orders
\(h_i=\lceil r/a_i\rceil\) satisfy
\(\sum_i 1/h_i\le1\) and \(h_i a_i\le 2r\), so the gain is the
full summed cap once \(\nu/L>2r\). The real ratio count has the
same exponent, with no lock price.

For each ratio \(H'/J'\), there are two possible signed pairs \(H',J'\), since the pair is reduced. The nonzero integers
\[
            d_h H'+d_j J' \quad\hbox{and}\quad d_h H'-d_j J'
\]
are divisible by \(u\) and \(v\), respectively: the original signed identities hold after multiplying by \(g\), and \(g\) is coprime to \(uv\). Thus \(u,v\) have divisor-many choices, and both original complete targets are now known. Once \(g\) is fixed, the exact signed summand differences and the rows are recovered.

#### Half-price completion on the fixed original target

Fix the oriented \(H',J'\), a divisor choice of \(u,v\), and the
common record. In particular \(c_1=vG\) is fixed. Put
\[
 k=(d_hH'+d_jJ')/u\in\mathbb Z\setminus\{0\},
 \qquad a_1-b_1=kg .
\]
Only the equal-valuation fresh gcd primes can fall outside
\(S\cup\mathcal P\cup\operatorname{supp}(H'J')\), and their raw
cap is the upper endpoint of the \(\tau\) bin. In particular the
prime of any unequal-valuation gcd is already visible in the
reduced ratio. Fix all exact powers of \(g\) on this known
polynomial-height pool at uniform height-subpower cost.

The remaining \(g\)'s have count \(O(n^{\tau/2+\eta}X^\varepsilon)\)
at sufficiently small fixed budgets. Otherwise, for some fixed
positive allowance, take a violating sequence of fixed-budget
height limits at budgets tending to zero. The first row and \(g\)
determine one another on the fixed record. Apply FT1 to these
actual fixed-\(c_1\) fibres after refining the residual summand
measures outside the known pool. Their positive exponent gives
positive nonatomic mass on each residual summand. A compact
positive band on \(a_1\) supplies a whole divisor \(F\mid a_1\)
of fixed positive numerical depth at arbitrarily small fixed
radical price. Choose the band, price and depth margin from the
limiting measure before the inner counts. Now
\[
                kg=2a_1-c_1\equiv -c_1\pmod F.
\]
Primitivity makes \(c_1\) a unit at \(F\), so \(k,g\) are units
there too. The unknown part of \(g\) is in one prescribed unit
class modulo \(F\). RX1 counts it at exponent \(\tau/2\) up to
arbitrarily small allowances. The purchase of \(F\) and the inner
count errors fit within the original allowance, giving a
contradiction. This use of a varying fixed target and growing
pool support is uniform by the height-first formulation of FT1.
If \(\tau=0\) one can enumerate the already supported powers
directly.

Summing over the finitely meshed \(\tau\) range adds its half-price
completion to the half-price of the reduced ratio. Their \(\tau\)
terms cancel and give the claimed per-record exponent. No
matched-even support or parity record is needed for this count.

#### A direct pair test without a parity price

The ratio bound can also be used without a bound for the whole
product auxiliary. Suppose a chart of row exponent \(d\) supplies
at least \(n^{2d-\rho-o(1)}\) usable pairs, each with fresh
difference radical price \(s=r_Z-\ell_{\mathrm{old}}\) or less
after binning, and with one of the two stated ratio constraints
at price \(\mu\). Price the rich shared records by FP6 at
\(\rho+o(1)\). This part of FP6 only uses the divisor incidence
of a shared record on an original row. On a chart with small raw
overlap \(m\) one can use the smaller price \(\min(m,\rho)\).
The direct counting inequality is
\[
 2d-\rho \le \min(m,\rho)+s/2+\mu+o(1).
\]
In particular a sufficient strict saving test for a chart
catalogue of price \(b_0\) is
\[
                      s+4\rho+2\mu+4b_0<4/3
\]
with fixed positive slack. Neither the total auxiliary radical
nor the parity record occurs in this test. For an exponent-full
selection with \(\rho=o(1)\), negligible shared support and
arbitrarily cheap ratio constraints, the fresh shift alone
must have price at least \(4/3-o(1)\). With original outside
radical masses \(1-o(1)\), the private primes retained on the
auxiliary operands then force its total raw price to be at least
\(10/3-o(1)\). In particular this applies to low-overlap typical
pairs in a refined balanced family whose original complete
target has positive nonatomic mass: a compact positive band
outside the vanishing shared support supplies the cheap quotient
lock. For an arbitrary thin pair selection use its actual
\(\rho,m\) and ratio constraints instead of this typical-pair
specialization.

#### Consequence for the earliest raw FP2 budget

Write \(d_0=1/3+o(1)\) and \(D=2/3-\rho-o(1)\) for the original
row and selected pair exponents. Use the hypotheses and allowances
of FP6 and FP8, including the original three caps. In a raw FP2
viable regime \(q+3\zeta\le3D+o(1)\), FP8 gives the auxiliary
height floor and its separate two-cheapest estimate gives
\(r_Z\le D-\zeta+o(1)\). FP6 prices the rich shared record at
\(\rho+o(1)\). Thus FP12 gives
\[
  D\le \rho+\frac{D-\zeta-\ell_{\mathrm{old}}}{2}+\mu+o(1),
  \qquad
  \rho\ge \frac29-\frac23\mu+
                   \frac{\zeta+\ell_{\mathrm{old}}}{3}-o(1).
\]
Arbitrarily cheap positive-depth quotient divisors, or the
power-scale real imbalance, raise the necessary loss from
\(1/6\) to \(2/9\) in this sector. At the first raw FP8 boundary
the per-record exponent is at most \(1/4+\mu+o(1)\), so any fixed
\(\mu<1/12\) excludes it by a fixed margin. The result is uniform
across the original balanced dictionary range.

#### A discrete restriction on the remaining tied profiles

Fix \(0<\kappa<\mu<1/12\). For either of \(u,v\), take the product
of all its full prime powers at primes \(p\le n^\kappa\).
If this product is at least \(X^\beta\) for fixed \(\beta>0\),
it supplies a divisor of height \(X^\nu\) and radical at most
\(n^\mu\), where \(\nu>0\) depends only on the fixed parameters.
Indeed \(\operatorname{rad}(u),\operatorname{rad}(v)\le n^C\)
for fixed \(C\) by the original target caps. The powers with
\(L v_p(u)<a=\beta/(2C)\) supply numerical weight at most
\(aC\), using logarithms to base \(X\). Among the remaining
powers, if their total radical price is at most \(\mu\), use all
of them. Otherwise take powers until their cumulative radical
price first reaches \((\mu-\kappa)/2\). The overshoot is at most
\(\kappa\); the result has radical price below \(\mu\) and
numerical weight at least \(a(\mu-\kappa)/2\). Use the same
choice for \(v\). Each selected power can also be taken whole
in the corresponding original target: any factors of \(G\) at
those primes are already known on the shared record, and the
congruence uses the actual divisor of \(u\) or \(v\).

Consequently in an exponent-full raw FP8 limiting regime
\(\rho\to1/6\), all but power-smaller sectors have the small-prime
parts of **both** quotients of numerical weight tending to zero,
in the height-first sense. More precisely each fixed positive
\(\beta\) has a fixed count gap on the sector above. At that
boundary the radical of either quotient outside the shared record
and the prescribed pool is at most \(n^{1/6+o(1)}\). For any
fixed \(1/18<\kappa<1/12\), at most two primes \(p>n^\kappa\)
from that outside support can carry its remaining numerical
weight. If the prescribed pool has nonvanishing mass, or shared
primes occur on the quotient at unequal original valuations,
their numerical powers are kept as separate prescribed-support
parts in this statement. Thus the **private** quotient has
vanishing small-prime height and at most four remaining
prime-power factors. This condition is stronger than weak
atomicity of its scaled valuation measure. In particular a compact
positive nonatomic band on either private quotient would have
positive numerical weight on small primes at any fixed
\(\kappa>0\); trim its heavy atom windows as in CF1.

No such divisor or other inexpensive source of the fixed ratio
class has been extracted from arbitrary balanced pair families.
The \(2/9\) bound is a conditional obstruction for this
auxiliary budget, across all original dictionary exponents.

[Signed two-slot residue moment (RX7)](ABC_RESIDUE_EXPANSION_AND_JOINT_RAY_COUNTS.md#signed-residue-words-and-diagonal-energy-on-the-cells)
· [Joint real and residue moments (JRM1)](ABC_JOINT_REAL_RESIDUE_MOMENTS.md#an-arbitrary-list-of-congruence-and-logarithmic-bins)
· [Fixed-target nonatomic mass (FT1)](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#nonatomic-mass-outside-prescribed-prime-pools)
· [Compact-band full locks (CF1)](ABC_COMPLETE_LETTER_SMALL_FEATURES.md#limiting-setup-and-cost-of-fixed-prime-supports)
