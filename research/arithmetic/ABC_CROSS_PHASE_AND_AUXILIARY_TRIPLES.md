# Cross phases and auxiliary primitive triples

Fix n = X^L, constants w,B>0 and all normalized cap parameters. The budget L is fixed sufficiently small before the height limit. Let ai+bi=ci≤X be two positive primitive rows and aibi=Di zi² with Di positive squarefree. Repeated and swapped rows are excluded. For a common odd divisor M|gcd(c1,c2), M≥X^w, a lock a1≡a2 modulo M automatically excludes swaps once M>1. The results permit all six original radical caps independently. A prescribed pool P≤X^B can be enlarged by fixed polynomial-height data.

Write G=gcd(c1,c2), u=c2/G, v=c1/G,
Δ−=u a1−v a2, Δ+=u a1+v a2−uvG.
Both are nonzero: zero would make the two reduced primitive proportions identical or complementary, respectively. Define
U−=a1b2, V−=a2b1, U+=a1a2, V+=b1b2,
ν−=−4(U−+V−+2√(U− · V−)),
ν+= 4(U++V+−2√(U+ · V+)). Equivalently, if αi=ai−bi+2zi√(−Di) with positive imaginary part, and C=c1c2, then
ν±=α1α2+conjugate(α1α2)±2C.
In particular ν+−ν−=4C. Both values have height O(X²), and are nonzero. In the real quadratic field K=Q(√sf(D1D2)) when D1≠D2, the identities are
N(ν−)=16(U−−V−)²=16G²Δ−²,
N(ν+)=16(U+−V+)²=16G²Δ+².
For Di equal the same products hold with the formal opposite square-root signs, which are also rational integers.

## Distinct fields and one shift cap

**CP1.** Prescribe the real quadratic squareclass D△=sf(D1D2)>1, with D△≤n^g. Neither Di separately is prescribed. Prescribe M and a list S of cardinality R for a1a2 modulo M (in particular an absolute equal-residue cell has R=1). Define g−=gcd(a1,a2)gcd(b1,b2)=gcd(U−,V−) and g+=gcd(a1,b2)gcd(b1,a2)=gcd(U+,V+). For one choice σ∈{−,+}, assume
rad_{PM}(G Δσ) · rad_{PM}(gσ) ≤ n^q,
where fixed pool factors, and 2, are excluded from the radical. Then the number of pairs is O(R^{1/2} n^{q/2+η} X^ε) for every fixed η,ε>0 and sufficiently small fixed L. A prescribed polynomial-height support of Δσ can instead be incorporated in P. This gives half the remaining outside radical of G at zero cost for the actual powers of Δσ on that support. Alternatively a uniform bound for both shift signs permits taking the smaller price by partitioning into the two cases.

At each odd p|ci, the quadratic field associated to Di splits, because ai and bi are units and −aibi≡ai² modulo p. The factor αi has one active branch of order 2 ord_p(ci) and one unit branch with residue 4ai modulo p^ord_p(ci). For p^j||M the simultaneous unit branches occur in the biquadratic compositum because the imaginary fields are distinct. Their restriction gives one degree-one branch of the real quadratic field. There νσ reduces to 16a1a2 modulo p^j, for both choices of σ. Fix the choices at height-subpower cost since there are at most 2^ω(M) of them. Their product to the stated exponents gives a modulus of norm M in K.

At any odd p|G the same simultaneous branches give a unit, so at most one degree-one prime over p occurs in (νσ). An exclusive target prime divides neither Δσ nor G. At a p outside G dividing Δσ but not gσ, the two cross values Uσ,Vσ are equal and nonzero modulo p. Thus D△ splits, p is unramified, and the conjugate signs of the square root make one reduction of νσ a nonzero unit. Hence the ideal radical there again costs only p, and at all remaining p it costs at most p². The asserted product therefore bounds the absolute ideal-radical cap outside the prescribed pool. Apply AT1 to this integral element. Its field discriminant is ≤4D△, its conjugates are O(X²), and its absolute residue dictionary on the chosen modulus has cardinality ≤R. Per νσ its trace and norm determine the exact unordered Uσ,Vσ through their sum and squared difference. Factor the two known integers to recover all four positive operands at divisor cost, with both targets then deduced. This recovery does not require either target to be fixed.

## A scalar-free quotient at distinct fields

**CP5.** Keep the distinct-field setup of CP1, but
drop the list of actual absolute residues and allow
different targets. Put
\[
        \gamma=\nu_-/\nu_+,\qquad
        h=\Delta_-,\quad j=\Delta_+ .
\]
Prescribe the real field as in CP1. Write \(H\) for the
normalized sum of the off-pool radical masses of \(h,j\),
outside \(2M\), and \(B_1\) for an upper bound on the
outside mass of primes that divide exactly one of \(h,j\)
and also divide both original operand products.
Then the ordered pair count is
\[
                   O\bigl(n^{(H+B_1)/2+\eta}X^\varepsilon\bigr)
\]
up to fixed fine binning of the caps. The common target
radical and the absolute summand residue carry no
additional price. A fixed pool can contain the repeated
operand primes if this is available as a prescribed
record. The actual catalogue cost of such a pool must
be included when it varies.

On the selected unit branch at each \(p^e\parallel M\),
CP1 gives the same \(16a_1a_2\) for both \(\nu_\pm\).
Thus \(\gamma\equiv1\) in the ordinary ideal modulus
of norm \(M\). Its numerator and denominator ideals
have polynomial norms, and all its embeddings are
bounded above and below by powers of \(X\) since the
two original integral values have nonzero integral
norm and polynomial conjugates.

At odd primes of \(G\), the same degree-one branch
is a unit in both factors, so the \(2v_p(G)\)
contribution cancels at the other branch.
In particular a prime on \(G\) contributes to
the ideal of the quotient only if it divides
one of \(h,j\), and then it costs at most one
prime ideal of norm \(p\). At any other prime
outside \(2M\) and the fixed pool, the norm
formulas show that the quotient can have a
nonzero ideal valuation only on \(hj\). If
\(p\) divides exactly one of the shifts and
the corresponding two cross-products are
units, CP1 gives one degree-one prime of norm
\(p\) in the quotient. If those cross-products
are nonunits, the prime must divide both
original operand products. The total weight
of the prime ideals over an individual rational
prime is always at most the square price \(p^2\).
This also covers primes on both shifts.
The radical ideal formed by the union of
numerator and denominator supports consequently
has norm at most \(n^{H+B_1}\).

The absolute ideal moments of AT1 extend here
to fractional ideals as follows. Fix the part
on the prescribed pool and the full prime powers
on numerator or denominator orders above a
cutoff \(J\) at each fixed budget. The tail
ideal radical has norm at most \(X^{O(1)/(J+1)}\),
so these records have price inside any fixed
allowance. Group the remaining squarefree ideals
into slots of exact signed valuation. A slot of
norm cap \(n^\theta\) has at most \(n^\theta
X^{o_X(1)}\) possibilities, by ideal divisors.
The total of the slot caps is at most \(H+B_1\)
up to the finite mesh. Apply the AT1 Fourier
kernel and its ordinary ray group, which already
uses fractional ideals. In the diagonal moment,
products of a positive slot have integral ideals
of short norm, and those of a negative slot are
their inverses. If two products lie in one
torsion and centered-log class, an exact
generator \(\zeta\) of their quotient reduces
to one modulo the ideal modulus. Multiply
\(\zeta-1\) by the integer norm of its denominator
ideal. This gives an integral element of the
modulus. The AT1 conjugate bounds apply with
numerator and denominator norms both short, by
choosing the same arbitrarily small log box
and then sufficiently small fixed \(L\). Thus
any nonzero such element would have absolute
norm simultaneously at least \(M\) and less
than \(M\). This proves the diagonal. Ordered
factorizations of a fixed slot product have
subpower size. Choose the moment multiplier
before the budget threshold to absorb the fixed
class-number exponent, as in AT1. Unit packing
at the claimed embedding height finishes the
fractional element count.

Finally an exact \(\gamma\) determines two
rational squares by norms. Put
\(z=\gamma/(\gamma-1)\). With \(x_i=a_i/c_i\),
the two norm identities give
\[
    N(z)=(x_1-x_2)^2,\qquad
    N(1-z)=(x_1+x_2-1)^2 .
\]
Indeed \(z=-\nu_-/(4c_1c_2)\) and
\(1-z=\nu_+/(4c_1c_2)\).
A choice of signs for the two roots gives
the two proportions, hence at most four
ordered primitive row pairs. This finishes
the claimed count. For a selected pair sector
of full exponent, the expression
\((H+B_1)/2\) can be compared directly with
its row-pair exponent after adding the
actual price of the single real squareclass.
[Absolute ray moment proof (AT1)](ABC_ABSOLUTE_TARGET_RAY_MOMENTS.md#the-two-diagonal-comparisons)

## Common field and two independent shift caps

**CP2.** Suppose instead D1=D2, with no catalogue of possible common fields. The two ν values are rational integers. At every p^j||M record whether the two unit branches agree or are conjugate in the common field. At least one full part of M has size at least X^{w/2}. Retain the stronger part. If this is the conjugate part, swap a2 and b2 in the second row for the duration of the count. This reverses its branch orientation, interchanges Δ− and Δ+ up to signs, and changes the product residue a1a2 to its negative at that part. The hypotheses and count remain equivalent after this relabeling. These choices cost height-subpower, including if the active alignment signs vary over primes on the original family. Use the retained modulus M0 with aligned branch signs. A list of R possible product residues modulo the original M projects to a list of size at most R, with the preceding sign change as needed.

Let G_opp be the full part of G on primes where the two unit branches are opposite AFTER this choice of orientation. At any odd p|G where they agree, both ν values are units at p. At p|G where they are opposite, the formal conjugate expression has simultaneous unit branches, so its residue modulo p is nonzero. Thus the exact orders of ν− and ν+ are
2v_p(G)+2v_p(Δ−), 2v_p(G)+2v_p(Δ+),
respectively. Primitivity at a common target prime makes the two Δ orders unable both to be positive: a1c2/G and a2c1/G are units if depths are equal, and if depths differ either Δ is automatically a unit, and when depths are equal the difference of Δ+ and Δ− is the unit 2v a2 modulo p. At p∤G, either it is exclusive to one target and both ν are units, or at most one ν is divisible by p since ν+−ν−=4c1c2. Therefore the product of rational primes counted once for every distinct positive order in (ν−,ν+) outside P2M0 is bounded by
rad_{P2M0}(G_opp) · rad_{P2M0}(Δ−Δ+).
If the cap on the right is n^q, the pair count is
O(R^{1/2} n^{q/2+η} X^ε).
In particular independent caps α and β on the two outside shift radicals give q at most mass(G_opp)+α+β. A cap for their joint radical is better. In a prescribed absolute residue cell R=1, and the estimate uses no quadratic-field catalogue.

For completeness, the paired integer count follows by using the product of the two ordinary unit groups modulo M0. Both signed ν residues are 16a1a2, or use positive absolute values and adjust one fixed sign. The three possible slot templates are (p,p), (p,1), (1,p), with common minimum order and the exclusive positive excess. Fix supported prime powers and the full tails where any layer exceeds J, choosing J at fixed L; since the heights are O(X²) the radical of the tail is at most X^{4/(J+1)+o(1)}. Its support and powers can be fixed inside a prescribed n^η allowance using fixed-support enumeration. Round the remaining slot charges to a finite mesh with total H≤q+η/4. For a θ-slot use moment order 2m, m=ceil(T H/θ), for fixed T≥1 chosen before L and with zero slots omitted. Each component of a slot product has height at most n^{mθ}≤n^{(T+1)H}. After equality of scalar e-multiples in the group, the second product lies in height-subpower many exact product residue classes from bounded-order torsion at fixed J. For small L both components are below M0, so these classes determine the exact products. The factorization multiplicity is height-subpower. Fourier inversion and Holder, using L² for the permitted residue pairs and the remaining 2m moments for the slots, gives √R and n^{H/2}. The leftover Holder weight goes to the constant one. The small-budget requirement is independent of the subsequent J because the base slot norm uses only p regardless of its scalar power.

To recover from the actual pair (ν−,ν+) in the common-field case, their difference fixes C=c1c2 and
α1α2=(ν++ν−)/4 ± i√(−ν−ν+)/2.
The product is nonreal. Indeed αi/c_i have positive imaginary part, both of unit modulus, so the only possible real product is −1 in this normalization. This forces complementary reduced proportions and therefore swapped rows. Hence the displayed product determines the common imaginary quadratic field. The principal ideal of α1 divides the prescribed product ideal and has at most divisor-many choices uniformly, since its rational norm is bounded polynomially in X. Each such ideal has at most six generators. Dividing recovers α2. Their positive norm square roots and real parts recover the two original targets and summands. No field catalogue is needed.

If only equality of the two residue coordinates is prescribed rather than a small list of actual product residues, work in the quotient by the diagonal instead. All (p,p) layers then have trivial character, and must be fixed at full radical price. Remaining exclusive templates have exact moment diagonals because a small product compared with the identity in the other coordinate is determined by its ratio modulo M0. Thus this variant has one-half price for exclusive layers and full price for common ones. The loss is an explicit arithmetic limitation of this comparison.

## Compressed projective shifts for equal targets

**CP3.** Now let c1=c2=c, which may vary with the pairs. Put h=Δ−=a1−a2, j=Δ+=a1+a2−c, and write j/h=A/B in lowest terms with B>0 and h=Bk, j=Ak, allowing k signed. If the residues are equal modulo an odd M>1 with M|c, then both h and j are nonzero, j±h≠0, and j is an M-unit. The first two facts exclude repeated and swapped rows, and j±h=0 would make a primitive positive row have target 2. The last follows from j≡2a1 modulo M. Thus all of M divides B, and gcd(c,k)|2 since 2a1=c+(A+B)k and gcd(a1,c)=1. Outside 2 the integer k is coprime to every operand: if a prime divides h,j, then at that prime 2ai≡2bi≡c, so a common operand divisor would contradict primitivity.

The original four operands are (c±(A+B)k)/2 and (c±(A−B)k)/2. These are four pairwise nonproportional linear forms in (c,k), since h,j,j−h,j+h are all nonzero. With a fixed ratio A/B, all coefficient heights are O(X), and their pair determinants have polynomial height. Fix the squarefree kernel D△ of the product of those actual four operands. The degree-four product form has equation F(c,k)=D△z² after absorbing 2^4 into the square. Move one of its rational branch directions to infinity by a rational linear fractional transformation of polynomial height. Dividing the square coordinate by the square of the same linear denominator gives a cubic of polynomial height with three distinct rational roots. Scale the leading coefficient to obtain a monic Weierstrass equation. It has rational 2-torsion, and its actual rational x-coordinates have polynomial height, so Dujella's Theorem 1.1 and Corollary 1.3 give X^ε points uniformly. Every x-coordinate fixes the projective direction (c:k), and gcd(c,k)|2 leaves at most two positive scales (either primitiveness fixes the scale directly). This yields X^ε pairs for the fixed ratio and squareclass, including D△=1. If only the outside squarefree kernel has cap n^f, enumerate it at cost n^f and the choices on any prescribed polynomial-height support at height-subpower cost. Projective heights stay polynomial even when c is variable and no numerical h or j is fixed. [Dujella, Corollary 1.3](https://arxiv.org/html/2312.03655v1).

The four cross-gcds are pairwise coprime and coprime to c. Their same-role product G− divides h, and their opposite-role product G+ divides j. Since odd factors of k occur in no operand, MG−,odd|B and G+,odd|A. Thus any prescribed common cross-divisor records give a finite ratio catalogue of size
O((1+X/(M G−,odd))(1+X/G+,odd)),
or the actual joint remaining radical price of A B outside their fixed support, whichever is smaller. Numerical compression follows if the original cross coincidences exhaust the two shifts after the lock. The projective fibre sums directly over such a compressed catalogue.

A further conditional reconstruction uses actual cross residue data. Put ρ=a1b2/(a2b1), σ=a1a2/(b1b2). Then ρ−1=ch/(a2 b1), σ−1=cj/(b1 b2). Fix the full M-supported parts C0 of c and H of h at height-subpower cost. Modulo M the two reduced increments L−=(ρ−1)/(C0H), L+=(σ−1)/C0 are units in the rational local sense, because j and the complementary operands are units at every prime of M. Exact algebra gives
Hj/h = (L+/L−)(b2/a2) ≡ −L+/L− (mod M).
The reduced numerator and denominator of Hj/h are A and B/H. A catalogue of these actual normalized cross-word residues, at any valid available precision Q on which b2/a2 is also known, determines each ratio up to uniqueness when twice the product of the two actual height caps is less than Q. Indeed the cross determinant of two candidate reduced fractions has that upper bound and is divisible by Q. A saving requires the actual joint catalogue of the leading units and the reconstruction height inequality.

## Use on abc pair extraction

For distinct fields the charge requires only the real squareclass D△, and none of the exclusive target primes. The shift should be genuinely compressed: a raw support cost or the larger charge allowed by an inflated shift bound can absorb the improvement. For equal-target pairs the common target is naturally charged once at half price. Pair entropy extraction must overcome the cost of requiring equal targets. In a block with fixed d>0, n^{2d} rows and at most n^d targets, Cauchy gives n^{3d} ordered equal-target pairs but the diagonal repeated pairs and swaps must be excluded. A real-squareclass catalogue of price f and fixed shift-support catalogue of price t give a contradiction if f+t+q/2<3d by a fixed margin, where q bounds the remaining outside target radical. In the common-field sub-sector there is no field price, but the two original shift caps or their joint radical must be cheap enough and the conditional product residue list must be affordable. A fixed ratio catalogue of price r instead combines with its outside symmetric-kernel price f at total r+f by the projective fibre theorem. Extracting such a price from the three original caps and the occupied residue cells remains open.

## Two auxiliary abc triples in the common-field equal-target case

**CP4.** Suppose c is odd and the two original rows have one common parity field. Use the full disjoint target decomposition c=C_A C_O according to whether the two active branches agree or are opposite. After reducing to one strong part of M and relabeling as needed, M0 divides C_A. For the unchanged equal-residue setting, M0 divides h. In the relabeled opposite-residue setting it instead divides j. Let
U=(√(a1 b2)+√(a2 b1))², V=(√(a1 a2)−√(b1 b2))²,
U'=(√(a1 b2)−√(a2 b1))², V'=(√(a1 a2)+√(b1 b2))².
All are positive rational integers (h and j are nonzero), and U+V=U'+V'=c² exactly. At aligned primes over c, U,V are units and U',V' have valuations 2ord_p(c)+2ord_p(h), 2ord_p(c)+2ord_p(j) respectively. At opposite primes the roles of the two pairs are reversed. Since p|c implies h,j cannot both be divisible by p, their exact gcds are C_O² and C_A². At each prime outside c the products UU'=c²h² and VV'=c²j² show that no new operand prime is introduced. Divide by the gcds to obtain primitive rows
(U/C_O²,V/C_O²,C_A²), (U'/C_A²,V'/C_A²,C_O²).
The first operand of each auxiliary row has radical dividing \(\operatorname{rad}(h)\). The second operand of each has radical dividing \(\operatorname{rad}(j)\). The latter two statements include aligned or opposite primes of c by the exact orders just established. The auxiliary target C_A² is divisible by M0², and one normalized operand of the second auxiliary row is divisible by M0². Thus C_O>M0, and both auxiliary target heights are at least X^{2w0} for w0=w/2, or w0=w if the original lock orients entirely in one direction. More strongly, if h_A,j_A are the full C_A-supported parts of h,j, then h_A²+j_A²≤C_O². The corresponding inequality holds with A and O reversed. If c is even, the same statements apply to all odd factors, and power-of-two factors in gcd and auxiliary target must be tracked separately.

The two ordered auxiliary rows recover C_A,C_O by positive square roots, hence c and both ν values and then the original pair at subpower cost. Given the original c, either single auxiliary row suffices for this recovery, since its target determines which normalization was removed (a finite choice of the auxiliary sign is harmless). In particular full normalized radical caps α,β on h,j and ℓ on c give a field-free exponent [2(α+β)+ℓ]/3+η by applying the existing limiting one-third abc count to both auxiliary rows, with separate fixed fine bins for the two target radicals. If the original full target dictionary has a price d, retain the auxiliary row with the smaller target radical and fix c, obtaining exponent d+[α+β+ℓ/2]/3+η. These two numerical conclusions use FULL radical caps on the auxiliary triples: any prescribed pool primes have to be included in the cap, unlike the arithmetic ray counts. To keep the height-first convention, bin the logarithmic auxiliary target height in a finite mesh over [2w0,2], choose fixed upper budget values for each bin, and apply the height-shell estimates on each dyadic shell there. The mesh correction tends to zero at fixed w0, and the dyadic decomposition costs only a logarithmic factor. The height meshes and normalized allowances are fixed before decreasing the budget.

For a powerful-triple transfer, at every odd prime outside the actual cross-gcds G− and G+, at most one of U,U' or V,V' respectively can be divisible after removing the common target factors, if that shift prime survives in an auxiliary operand. Its valuation there is twice the corresponding shift valuation, and the derived targets also square the corresponding original target parts. This shows exactly which three defects must be controlled before applying a strict high-powerful count. The original three radical caps alone do not give such a shift valuation cutoff.

[Absolute ray moments (AT1)](ABC_ABSOLUTE_TARGET_RAY_MOMENTS.md#algebraic-integers-with-a-ray-residue-dictionary)
· [Four distinct form factors (MTK1)](ABC_MULTI_TARGET_MIXED_KERNEL_COUNTS.md#the-mixed-kernel-count)
· [Limiting small-budget count](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)
· [High-powerful transfer](POWERFUL_TRIPLES_ABC_ENDPOINT_TRANSFER.md#theorem-a-strict-saving-from-endpoint-diffusion)


## The common-field transfer for moving and even targets

**CP6.** The two auxiliary triples of CP4 extend to
arbitrary positive \(c_1,c_2\) in the same parity field.
Use the same common absolute lock
\(M\mid(c_1,c_2)\), \(M\) odd and \(M\ge X^\nu\),
with \(a_1\equiv a_2\bmod M\).
Use the nonzero reduced shifts \(h=\Delta_-\),
\(j=\Delta_+\) as before. Apart from repeated
or swapped rows, the two primitive auxiliary triples
are recoverable from the genuine pair and conversely
recover it at height-subpower cost. Both auxiliary
targets have height at least \(X^{\nu/4}/2\)
after the subpower branch orientation choices.

The four positive integer values \(U,V,U',V'\)
are as in CP4, but now
\[
 U+V=U'+V'=c_1c_2,\qquad
 UU'=G^2h^2,\qquad VV'=G^2j^2.
\]
These identities use the unscaled cross products
and \(G=(c_1,c_2)\).
Partition the odd part of \(G\) into its aligned
part \(G_A\) and its opposite part \(G_O\) according
to the same unit branches as CP2. At an odd prime
of \(G\), one pair consists of units and the other
has orders \(2v_p(G)+2v_p(h)\) and
\(2v_p(G)+2v_p(j)\). The shifts cannot both
vanish modulo \(p\): when the two target orders
differ both are units, and with equal orders
their difference is \(2v a_2\bmod p\), a unit.
At an exclusive target prime the original phase
norms are units. Away from both targets a pair
cannot have a common prime divisor since it sums
to \(c_1c_2\).

Consequently \(d_O=(U,V),d_A=(U',V')\) have odd
parts \(G_O^2,G_A^2\) respectively, and both
divide \(G^2\). To check the latter at two,
if \(uv\) is odd then \(d_O,d_A\mid c_1c_2=uvG^2\)
suffices. If \(uv\) is even, one of \(u,v\)
is even and the other unscaled row at the
corresponding even target has odd operands.
Thus \(h,j\) are odd and the two product
identities give the bound at two.
Divide the two pairs by the respective gcds.
The resulting auxiliary target odd parts are
\[
                    (uv)_{\rm odd}G_A^2,\qquad
                    (uv)_{\rm odd}G_O^2.
\]
In each auxiliary the first operand has
radical dividing \(2\operatorname{rad}(h)\)
and the second dividing \(2\operatorname{rad}(j)\).
The odd gcd of the two targets gives
\((uv)_{\rm odd}\), hence the original product
\(c_1c_2\) and the normalizations are determined
away from two. Enumerate the powers of two
at logarithmic-power cost. Restoring the
two phase values recovers the nonreal product
of the two original integral quadratic
elements as in CP2. The common field and
the two original rows then follow at divisor cost.

It remains to verify the positive height floor
for both auxiliary targets. Select the stronger
alignment part \(M_0\ge M^{1/2}\) and swap the
second row if necessary so that this part is
aligned. The original same-residue relation
becomes \(a_2\equiv\pm a_1\bmod M_0\) with one
common sign. The target \(T_A=c_1c_2/d_O\)
is at least \(M_0^2\). Put \(R_0=G^2/d_A\),
so \(T_O=c_1c_2/d_A=uvR_0\).
In the plus residue orientation,
\(U\equiv-4a_1^2\) and \(h\equiv(u-v)a_1\bmod M_0\).
The identity \(U'/d_A=R_0h^2/U\) gives
\[
                 4U'/d_A+R_0(u-v)^2\equiv0\bmod M_0 .
\]
It is a positive integer at most
\(R_0(u+v)^2\le4T_O^2\), using \(U'/d_A\le T_O\)
and \(u,v\ge1\). In the minus residue orientation,
use \(V'\) and \(j\): \(V\equiv-4a_1^2\)
and \(j\equiv(u-v)a_1\). This gives the
same lower bound \(T_O\ge\sqrt{M_0}/2\).
A complete residue class at a fixed odd deep
common divisor is enough for both bounds.

With full normalized shift radical caps
\(b_h,b_j\), and full radical caps \(r_U,r_E\)
on \(c_1c_2\) and \(uv\) respectively, the limiting
one-third estimate thus gives pair exponent
\[
              [2(b_h+b_j)+r_U+r_E]/3+\eta
\]
under the same finite height meshes as CP4.
Indeed the product of the auxiliary odd
target radicals equals the product of the
odd radicals of \(c_1c_2\) and \(uv\).
No common field catalogue enters the count.
If \(c_1c_2\) is fixed, either one of the
two auxiliary triples suffices for recovery
after fixing the bounded choice of its sign,
so its individual actual radical price can
replace the sum. The target orders in each
auxiliary are those of \(uv\) plus twice the
respective branch part of \(G\), at all odd
primes. The operand-order statement of CP4
also extends away from two and the actual
cross-operand gcds. A high-powerful transfer
requires these actual surviving orders.
