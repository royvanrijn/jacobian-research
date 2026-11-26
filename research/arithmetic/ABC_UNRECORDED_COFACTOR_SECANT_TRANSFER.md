# Two-letter counts with a short unrecorded cofactor

A bounded radical cap on the third letter in the two-letter
count can be replaced by a bounded radical cap on a factor,
a fixed coefficient, and a sufficiently short unrecorded
cofactor. The unrecorded part can have arbitrary prime support.
The resulting count applies to a primitive secant coordinate
even when the physical difference contains a common numerical
scaling factor with uncontrolled radical.

## The two-letter inequality with an unrecorded cofactor

Fix \(B,C,\epsilon_0>0\). There are \(\theta,L_0>0\)
depending only on these constants with the following property.
At fixed \(0<L<L_0\), put \(n=X^L\). Count positive pairwise
coprime complete values \(U,V,S\) in a signed equation
\(\pm U\pm V\pm S=0\), with maximum in \((X/2,X]\).
Assume \(0\le a,b\le C\) and
\[
 \operatorname{rad}(U)\le n^a,\qquad
 \operatorname{rad}(V)\le n^b,\qquad
 S=A Z W,\qquad
 \operatorname{rad}(Z)\le n^B,\quad 1\le W\le X^\theta.
\]
Here \(A\) is a fixed positive integer, and \(Z,W\) are positive
integers whose representations may vary with the solution.
The count is
\[
  O_{B,C,\epsilon_0,L,\varepsilon}
        \bigl(X^{((a+b)/2+\epsilon_0)L+\varepsilon}\bigr).
 \tag{UC1}
\]
It is uniform in \(A,a,b\) and the signs. If \(U\) is prescribed
completely, retain its stated raw cap and replace \((a+b)/2\)
by \(b/2\). The same applies with \(U,V\) interchanged.
Representations of \(S\) need not be counted with multiplicity.

### Height-first reduction and cheap blocks in the recorded factor

Use a failure sequence with the candidate parameters
\(\theta,L_0\) tending to zero. At each failing fixed \(L\)
take the height limit, using inner losses smaller than the
fixed normalized failure margin. Use independent anatomic
factorizations of \(U,V,Z\), fixing their initial coefficients
at cost \(o(L)\) in base \(X\). This is the ordinary reduction
for integers of height at most \(X\) with bounds on their
radicals. Choose one representation of each original row.
All root lists are finite before taking the corresponding
height limit, and the total raw root height is \(O_{B,C}(L)\).
If an entire letter is fixed, choose its factorization once
and charge zero for its root dictionaries. Its raw cap is
still in the hypotheses of the moment estimate.
The fixed factor on \(Z\) can be included in \(A\).
Refine into dyadic boxes for all the recorded roots. The
other initial coefficients have numerical height \(X^{o(1)}\).
[Anatomic decomposition for radical caps](https://arxiv.org/html/2410.12234v2#S2)
· [Height-first conventions](ABC_SECANT_ENDPOINT_TRANSFER.md#the-height-first-setting)

Write \(d\) for the limiting normalized row count exponent
and \(w_U,w_V\) for the limiting normalized singleton
dictionary cost totals. In the case of a fixed \(U\) we
have \(w_U=0\). In all cases the strict failure gives
\[
                     \tau=d-\frac{w_U+w_V}{2}>0.
\]
Any block of roots of \(Z\), with orders at least
\(\gamma/L\) for fixed \(\gamma>0\), whose full powers
have positive limiting numerical weight must cost at least
\(\tau\). Indeed, fixing that block gives a full divisor
of \(S\). The two other complete values are units and
are equal up to sign modulo the divisor. The diagonal
moment count behind (GD1) only requires the raw root and
dictionary caps on \(U,V\). It bounds the conditional count
by exponent \((w_U+w_V)/2\) up to vanishing normalized
loss. This argument charges the selected roots of \(Z\)
in their own decomposition. The unrecorded cofactor has
no role in the congruence.
[Full-modulus moments (GD1)](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#diagonal-moments-for-a-joint-dictionary)

Set aside the bounded number of root positions on \(Z\)
with normalized raw cap at least \(\tau/4\). Partition the
remaining roots into greedy blocks with cost at most
\(\tau/2\); every completed block except the last has
cost greater than \(\tau/4\). Thus the number of blocks is
bounded at the small-budget limit. For each fixed \(\gamma\),
their projections onto orders at least \(\gamma/L\)
have vanishing numerical weight, by the preceding lock
bound. The discarded lower orders have numerical weight
\(O_B(\gamma)\). Let \(\gamma\) decrease to zero. The product
of the powers outside the set-aside positions has numerical
height \(X^{o(1)}\). Send the retained small-budget
subsequence further as needed for the finite cost limits.

Consequently \(S\) is its fixed coefficient times a
bounded number of narrow roots, each of height \(X^{O(L)}\),
times a wider integer of height \(X^{o(1)}\).
The latter includes \(W\), since \(\theta\to0\) in the
failure sequence. The narrow roots in the retained list
whose scaled order tends to zero may also be moved into
the wider factor. At a positive-weight full modulus from
another letter, the resulting complete value or the ratio
of two complete values can be released by (GW) at
vanishing normalized exponent cost, as in the bounded-root
argument for (GD3). A fixed coefficient is a unit at
opposite moduli by primitivity. Wide and narrow factors
on the same letter can have common prime support. The
signed amplification counts exact positive output ratios,
whose numerator and denominator recover the primitive pair.
It colours prime supports by numerator and denominator
before amplification, so the product outputs have
divisor-many decompositions. A fixed multiplicative factor
just changes the residue target, or the centre of a real
strip. For the real strip take a single representing tuple
per exact output.
[Wide cofactor count (GW)](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-bounded-number-of-small-roots)
· [Signed amplification with overlap](ABC_WIDE_COFACTOR_CONGRUENCES.md#signed-amplification-with-overlapping-groups)
· [Real strip count (RUW1)](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-real-strip-with-unequal-root-heights)

### Exchanges of locks on the two fully tracked letters

Relabel \(U,V\) so \(w_V\le w_U\); then \(w_V<d\).
First fix \(V\). Any positive numerical weight high
block on \(U\) has cost at least \(d-w_V\), since its
full modulus releases \(S\) at a fixed residue through
(GW). After fixation \(S,V\) determine the original row.
If \(U\) has positive limiting numerical height, find
a high root of cost \(\beta_U\ge d-w_V-\xi\) for any
fixed sufficiently small \(\xi>0\). To do so, pack all
root positions of cost at most \(d-w_V-\xi\) into blocks
of this cap. There are boundedly many blocks, since any
two adjacent ones have total cost above the cap. For
each fixed high-order cutoff their high parts have
vanishing numerical weight. The entire low-order part
has arbitrarily small numerical weight by shrinking the
cutoff. A remaining positive-weight high position has
the required cost. Choose \(\xi<2d-w_U-w_V\), so the
complement of this position has cost strictly below \(d\).

If both \(U,V\) have positive limiting numerical height,
fix this complement on \(U\). Now any positive-weight
high block on \(V\) must cost at least
\(d-w_U+\beta_U>0\); its full modulus releases \(S\)
and the selected root on \(U\) through their ratio.
The same concentration argument gives a positive-weight
high root of cost \(\beta_V\) on \(V\). Fix the complement
of each selected root in turn, then lock the opposite root.
Releasing the remaining root and \(S\) gives
\[
        d\le \beta_U+w_V-\beta_V,\qquad
        d\le \beta_V+w_U-\beta_U .
\]
If two full-weight locks were chosen for different
purposes they are used in separate counts. Adding the
bounds contradicts \(w_U+w_V<2d\).

If one of \(U,V\) has vanishing numerical height, the
other and \(S\) are tall and their signed equation puts
their positive ratio in a relative strip of width
\(X^{-1+o(1)}\) about one. If the tall tracked letter has
cost below \(d\), fix it and release \(S\) with (RUW1).
Otherwise it is \(U\). Fix its complement of cost
\(w_U-\beta_U<d\), and release its one unfixed root and
\(S\) in the same strip. Exact reduced ratios determine
the two primitive tall values and hence the short letter.
Both branches give a contradiction. The leading numerical
heights in these branches are those of the complete
tracked values \(U,V\); neither contains an untracked
large variable coefficient. This proves (UC1), including
the zero-cost specialization for a prescribed tracked
letter. If one uses exact truncation and fixes a high
valuation tail on a tracked letter, its numerical weight
must remain in the argument. A fixed tail of positive
weight is a cost-zero full lock and already contradicts
the positive lower cost in the corresponding exchange. [Earlier two-letter exchanges (GD3)](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#the-two-letter-effective-cost)

## Primitive secant coordinates with an arbitrary scaling gcd

Use the balanced positive primitive shell of (SC1). For
two distinct rows \(x=(a,b,c)\) and \(y=(a',b',c')\), write
\[
 g_*(x,y)=\gcd(|a-a'|,|b-b'|),\qquad
 u_T(x,y)=\frac{T(y)-T(x)}{g_*(x,y)}
\]
on each role \(T=a,b,c\). A zero coordinate is excluded
from its radical test. For a prescribed integer \(Q\ge1\),
write \(R_Q(t)=\prod_{p\mid t,\ p\nmid Q}p\) for \(t\ne0\).

Fix \(B,K,\zeta>0\). At sufficiently small fixed \(L\),
with sufficiently small balance tolerance \(\kappa\), there
is \(\sigma>0\) with the following property. For each
subfamily \(\mathcal B\) let \(F_{\mathcal B,T}\) be its
largest complete-\(T\) fibre. Let \(\mathcal A\) be any
other subfamily of the shell and allow a pool integer
\(Q_x\le X^K\) to be prescribed for each \(x\in\mathcal A\).
Then
\[
 \begin{split}
 &\#\{(x,y)\in\mathcal A\times\mathcal B:
           x\ne y,\ u_T(x,y)\ne0,\
           R_{Q_x}(u_T(x,y))\le n^B\}\\
 &\quad\ll_{B,K,\zeta,L}
 n^{2/3-\sigma}
 +|\mathcal A| F_{\mathcal B,T}\,
                       n^{(1/3+\kappa)/2+\zeta}.
 \end{split}
 \tag{UC2}
\]
The estimate is uniform in the prescribed pools, subfamilies
and role, with all structural parameters fixed before
the height limit. Only one coordinate of the primitive
direction is tested. The raw radical of the physical
shift and that of the common scaling gcd are unrestricted.

First discard a power-small pair class to ensure that
\[
       T(y)/g_T\ge X^{h_0},\qquad
       g_T=\gcd(T(x),T(y)),
 \tag{UC3}
\]
for an absolute \(h_0>0\). This is the compact positive
degree band and first-endpoint prime-pool discard in the
proof of (SC5): for each fixed \(x\) the common prime
support of \(T(y)\) with \(T(x)\) has small radical except
on a power-smaller family, by (EB5). The full prime powers
of \(T(y)\) outside \(T(x)\) then survive division by
\(g_T\) and have positive numerical weight. The global
one-third bound, with a small fixed allowance, makes the
discarded pair class \(O(n^{2/3-\sigma_0})\).
[Height of the normalized coordinate pair (SC5)](ABC_PRIMITIVE_SECANT_RADICAL_AMPLIFICATION.md#the-radical-of-a-physical-coordinate-difference)

Choose the cofactor tolerance \(\theta\) in (UC1) for
the fixed supplied cap \(B\) and the small counting
allowance needed here, allowing for the conversion
from base \(X\) to any output shell \(Y\ge X^{h_0}\).
In particular one may use \(B/h_0\) and \(C=1/h_0\)
with budget unit \(L\) at height \(Y\). Then choose \(0<\eta<h_0\theta/2\). For this \(\eta\),
apply the secant gcd peeling theorem. There is a
pool integer \(P_x\le X^{K_1(\eta)}\) for each original
\(x\), depending only on \(x\) and the original family,
such that off a power-small pair class
\[
             g_*(x,y)=d\,w,\qquad
             d=\gcd(g_*(x,y),P_x),\qquad w\le X^\eta .
 \tag{UC4}
\]
The peeling theorem uses the polynomial curve bound with
degree depending on \(\eta\). The cofactor tolerance was
chosen before fixing this degree and before shrinking the
input budget threshold. Its threshold does not depend on
\(K_1\); fixed-pool enumeration and all inner height
constants can depend on \(K_1\). Discard this second
exceptional pair class, which is independent of \(Q_x\).
[Polynomial curves and the secant gcd](ABC_BOUNDED_DEGREE_CURVES_AND_SECANT_GCD.md)

For fixed \(x\) choose \(g_T\mid T(x)\) and \(d\mid P_x\).
Their choices have height-subpower cost by the divisor
bound. Put \(h_1=\gcd(d,g_T)\) and \(h_2=g_T/h_1\).
Since \(\gcd(d/h_1,h_2)=1\) and \(g_T\mid g_* |u_T|\),
one has \(h_2\mid w|u_T|\). Put \(h_3=\gcd(h_2,w)\).
The normalized difference triple is positive primitive
after moving the largest complete value to the target
position. Its two tracked values are \(T(x)/g_T\) and
\(T(y)/g_T\), and its third value is
\[
      \frac{|T(y)-T(x)|}{g_T}
      =\frac d{h_1}\,
        \frac{|u_T|}{h_2/h_3}\,
        \frac w{h_3}.
\]
The last factor is at most \(X^\eta\). Split off the entire
\(Q_x\)-supported part of the middle factor. It has
height-subpower many possible values uniformly for fixed \(x\).
Fix it and absorb it into the coefficient \(d/h_1\).
The remaining middle factor has radical at most \(n^B\)
by the hypothesis in (UC2). All these steps count
the exact normalized difference triple without requiring
a bound on the raw radical of \(g_*\) or a prescribed
value of \(w\). [Uniform fixed-support count](ABC_EFFECTIVE_RADICAL_BUDGET.md#counting-by-a-sublogarithmic-prime-partition-sum)

Now the first tracked value \(T(x)/g_T\) is fixed.
Apply the fixed-letter version of (UC1) to the two
tracked caps, both at most \(n^{1/3+\kappa}\).
By (UC3) the derived height is bounded below by a fixed
positive power of \(X\). The raw relative cap on the
middle factor in its output shell \(Y\) is at most
\(B L\log X/\log Y\). One can keep \(L\) fixed as the
budget unit and enlarge \(B,C\) by a fixed factor in
the bound (UC1). Alternatively one can round these relative
caps on a finite mesh at the original fixed budget.
Our choice of \(\eta\) keeps the third factor inside the
short numerical allowance at each \(Y\). The dyadic
shell sum and all the preceding subpower choices are
absorbed in \(\zeta\) by choosing the component allowances
smaller. The output and \(g_T\) recover \(T(y)\), which
has at most \(F_{\mathcal B,T}\) preimages. Sum over \(x\)
and add the discarded pair classes. This proves (UC2).

### Consequences for bounded directions and intermediate multiplicity

A relation of exponent \(2/3-o(1)\) on which a given
primitive secant coordinate has bounded normalized
radical forces the corresponding second-endpoint letter
dictionary to have exponent \(1/6+o(1)\) on a full-exponent
projected family. With a raw cap the test is symmetric
and gives the same conclusion on the first side. The
same one-sided conclusion holds with a pool chosen from
the first endpoint. For proof, bin the original
second-endpoint \(T\)-fibre sizes before retaining the
pair relation, at height-subpower cost. Equation (UC2)
forces the retained fibre exponent to be at least
\(1/6-o(1)\), while (FL1) gives its upper bound.
The ambient input family has exponent at most \(1/3+o(1)\),
so the corresponding values have dictionary exponent at
most \(1/6+o(1)\). The projection of the pair class still
has exponent \(1/3-o(1)\), which forces the lower
dictionary bound. Binning all six original degrees yields
both minimum dictionaries on all three roles when the
whole primitive direction radical has a fixed bound
\(n^B\) on a critical pair relation. Pairs with coinciding
complete values in any role have a power-smaller count
by (FL1). No bound on the common scaling gcd is needed.
[Fixed complete-letter caps (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

In the central least-target-repair normalization
\(Z_{\rm rep}=Y^{1/M}=n^{1/3+o(1)}\), a critical stratum of
target multiplicity at most \(Z_{\rm rep}^{u+o(1)}\) for
fixed \(u<1/2\) has at most that many rows per complete
target as well. For every fixed \(B,K\), the primitive
target-coordinate secant \(u_c\) therefore has outside
radical greater than \(n^B\) on all but a power-small
fraction of pairs. The prime pool can depend on the first
primitive quotient. This also applies successively to
a fixed number of affine secants from a central row,
using each earlier secant support in the next pool.
The discarded geometry in (UC2) is independent of this
choice of pools. Consequently the bounded-direction
alternative on a full-exponent pair relation is restricted
to the minimum complete-dictionary regime. Forcing that
alternative, or otherwise estimating the occurrence of
strata with large primitive secant radicals, remains open.
[Repaired target multiplicity range (FL5–FL8)](ABC_FIXED_LETTER_RADICAL_CAPS.md#size-of-a-minimally-repaired-target-radical-fibre)

## Numerically large scaling gcds without a direction cap

Fix \(\tau>0\). In a height-first endpoint sequence,
suppose an ordered pair relation of exponent \(2/3-o(1)\)
has \(g_*(x,y)\ge X^\tau\) throughout. Some full-exponent
second-endpoint subfamily then has minimum complete
dictionaries on at least two roles, even with no radical
bound on the primitive directions. More precisely,
a role whose escaped first moment in the refined endpoint
profile is strictly less than \(\tau\) also has a minimum
dictionary. In particular, with no escaped first moments
all three dictionaries are minimal.

Apply the peeling theorem with \(0<\eta<\tau/4\), discard
its exceptional pairs, and fix an anchor \(x\) with
\(n^{1/3-o(1)}\) surviving neighbors. Formula (UC4) gives
a divisor \(d\mid P_x\) of height at least \(X^{\tau-\eta}\).
Fix \(d\) and fix the full \(P_x\)-supported part \(B_T\)
of the second endpoint on each of its three roles at
height-subpower cost. All retained values satisfy
\(T(y)\equiv T(x)\bmod d\). The quotient \(T(y)/B_T\)
lies in one unit class modulo
\[
                      m_T=\frac{d}{\gcd(d,B_T)}.
\]
Indeed, divide the linear congruence by \(\gcd(d,B_T)\)
and invert the remaining multiplier. The quotient is
a unit since the supported part removed all primes
in \(P_x\). The three divisors \(\gcd(d,B_T)\) of \(d\)
are pairwise coprime. Therefore at most one exceeds
\(\sqrt d\), and at least two \(m_T\) are at least
\(X^{(\tau-\eta)/2}\). Apply (BC1) to each of these two
unit-residue quotients. Each has radical at most
\(n^{1/3+o(1)}\) and represents its original complete
letter bijectively with \(B_T\) fixed. Its dictionary has
exponent at most \(1/6+o(1)\), and (FL1) gives the lower
bound at full row exponent.
[One-letter residue count (BC1)](ABC_BOUNDED_DEGREE_CURVES_AND_SECANT_GCD.md#one-complete-letter-in-a-unit-residue)

For a further role \(T\), suppose its escaped moment
is \(e_T<\tau\). This moment is the limiting numerical
height on degrees going to infinity in scaled degree.
Take \(\eta>0\) so small that \(e_T<\tau-4\eta\). Choose
a fixed cutoff \(R\) at a continuity point so that the
numerical tail of \(T(y)\) on \(Lv_p(T)>R\) is below
\(\tau-3\eta\) in the relevant refined profiles at late
budgets. Discard second endpoints with a fixed positive
amount \(\xi\) of normalized radical on \(P_x\), by (EB5).
Choose \(\xi R<\eta\). Among the survivors, the
\(P_x\)-supported part on degrees at most \(R\) has
numerical height exponent at most \(R\xi\).
It follows that
\(\log_X B_T<\tau-2\eta\), so \(m_T\ge X^\eta\).
Apply the same unit-residue count. When several role
tail bounds are used, choose the parameters for all
three before taking the budget limit.
[Exact valuation profiles and escaped numerical weight](ABC_SECANT_ENDPOINT_TRANSFER.md#the-height-first-setting)
· [Uniform discount on prescribed pools (EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

There is a power-saving form for a specified role.
Fix \(R,\tau,\eta>0\) with \(3\eta<\tau\), and
require the numerical tail above \(R\) on the second
endpoint in that role to be at most \(\tau-3\eta\)
in base \(X\). If its complete-letter fibre cap is
\(n^{1/6-\beta}\) for fixed \(\beta>0\), a power-small
part of the full pair scale has \(g_*(x,y)\ge X^\tau\).
Use the same first-endpoint conditioning on \(d,B_T\),
and sum the unit-residue dictionary bound times the
fibre cap over the first endpoints. Apply (BC1) with
normalized error smaller than \(\beta/4\), and take the
balance and height allowances smaller than this margin.
The power gaps of the peeling and overlap discards
are kept fixed. This gives a numerically small common
gcd on typical pairs in the intermediate target
multiplicity range whenever the target tail is below
the specified numerical cutoff.

### Affordable recorded divisors and the low-order remainder

The same inequality (UC2) holds with a weaker test on the
primitive coordinate. For some fixed \(\theta_1>0\) depending
on the counting allowances, replace
\(R_{Q_x}(u_T)\le n^B\) by the existence of
\[
   |u_T|=S_0 W_0,\qquad R_{Q_x}(S_0)\le n^B,\qquad
                           W_0\le X^{\theta_1}.
\]
To verify the extension, choose \(\eta,\theta_1>0\) with
\(\eta+\theta_1<h_0\theta/2\) in (UC2). In the factorization
of the normalized coordinate difference, replace \(w\)
by \(wW_0\) and \(|u_T|\) by \(S_0\), and redo the division
by \(h_3\). The only change is the upper exponent of the
unrecorded cofactor. It is still within (UC1). The same
preimage argument applies, so the representations of
\(u_T\) need not be canonical.

In particular, consider a critical intermediate-multiplicity
class as above and fix any \(B>0\). In each tested primitive
target difference put all full prime powers of order at least
\(1/(BL)\), and all primes of the pool, into \(S_0\).
The radical of its outside-pool support is at most \(n^B\),
since the primitive difference is at most \(X\).
Consequently a power-small part of the pair scale suffices
to discard every instance in which the remaining
out-of-pool powers of order less than \(1/(BL)\) have
numerical weight at most \(\theta_1\) in base \(X\).
Thus the obstruction in this multiplicity range persists
even after taking the entire high-order part of a primitive
coordinate and all factors on the chosen pool out of the count.
