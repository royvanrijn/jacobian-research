# Small integer parameters in an affordable abc fibre

An affordable pair of full divisors makes distinct targets identify
distinct triples. The interval count extends to windows whose centres
are linear functions of finitely many small integer parameters, with
arbitrary real coefficients. This applies to all balanced endpoint
profiles.
[Affordable records and interval entropy](ABC_WIDE_COFACTOR_CONGRUENCES.md#numerical-interval-entropy)

## Sparsity for finitely many small integer parameters

Write \(n=X^L\) and take the height limit before the small-budget
limit. Fix \(\eta,\epsilon,C>0\) and an integer \(R\ge1\). Consider
a coprime full-divisor record \(F\mid U,\ G\mid V\), with combined
numerical weight exceeding the shorter complete height by at least
\(\eta\). Permute the free letters so that \(V=Gv\) is the shorter
one, write \(U=Fu\), and use the fixed signed relation for \(T\).

There is a fixed \(\delta>0\) with the following property. A
subfamily of this record has at most \(n^{\epsilon+o(1)}\)
distinct targets if it admits integer vectors \(\boldsymbol z\)
with \(R\) coordinates, including the constant coordinate \(1\),
\(|z_i|\le n^C\), and a common real linear function \(\ell\) with
\[
                   |T-\ell(\boldsymbol z)|\le X^\delta .
 \tag{AF1}
\]
The vectors need not distinguish targets. In particular, the
centres can be values of a polynomial of fixed degree in boundedly
many integer parameters of height \(n^{O(1)}\), by using their
monomials as coordinates.

### The exact dependence of the cofactors

Choose independent \(z_i\)-columns of the evaluation matrix
spanning all the feature columns. Include the column \(T\) if it
increases the rank. For any square minor formed by these basis
columns and then \(v\), multiplication of its last column by \(G\)
makes its determinant zero modulo \(F\). Indeed, \(Gv\) is a
signed copy of \(T\) modulo \(F\), and \(T\) lies in the span of
the other columns whether or not it was included. As \(G\) is a
unit modulo \(F\), the original integer determinant is divisible
by \(F\).

If \(T\) was included, subtract \(\ell(\boldsymbol z)\) from it
using the chosen feature basis. Its entries are now bounded by
\(X^\delta\). Writing \(w_G\) for the weight of \(G\), the minor
has absolute value at most
\((R+2)!\,n^{CR}X^{\delta+h_V-w_G+o_X(1)}\).
For \(\delta<\eta/2\) and small enough \(L\) this is less than
\(F\). Thus \(v\) lies exactly in the span of the feature columns
and \(T\), as does \(u\) by the signed relation. All three complete
letters are rational linear functions on the rows
\((\boldsymbol z,T)\).

On any subfamily, again choose the independent \(z_i\)-columns
first, then include \(T\) only if necessary. Any full-rank minor
of this basis is an integer of absolute value at most
\[
 D_X=(R+1)!\,X^\delta\prod_{i=1}^R \max(1,\sup|z_i|),
 \tag{AF2}
\]
by the same subtraction if \(T\) is present.

### Full factors force the rank to drop

For uniformity, suppose the bound fails at every positive width
exponent for one fixed allowance \(\epsilon\). Choose a diagonal
sequence of actual endpoint fibres with \(\delta_j\to0\) and
normalized counts above \(\epsilon\). Take the inner height
limits before passing to the outer sequence. In particular the
numerical exponent of \(D_X\) tends to zero.

Select spanning witness rows of \((\boldsymbol z,T)\) from the
current subfamily. For every moving complete letter there is a
nonzero difference of its values on two witnesses. Otherwise the
letter, as a linear function on the rows, would be constant on
the whole sample, since the span includes the constant coordinate.
For every moving letter strip and fix the parts of its original
roots supported on primes of the chosen witness difference.
Only polynomial-height integers supply these prime supports, so
the [fixed-support count](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound)
makes the cost subpower at fixed truncation. A constant letter
can be fixed outright. Rebox the residual roots and pass to weak
limits of their scaled-degree measures. Witness rows need not
survive the stripping since they will be used only for their
values in the span before stripping.

If a moving letter has positive residual nonatomic mass, choose
a whole block in a bounded positive degree band. The block can
have arbitrarily small fixed normalized cost and fixed positive
numerical weight. Use continuity bands and omit shrinking atom
windows to make individual retained root caps tend to zero.
Freeze the block at its root cost. Its full divisor is coprime
to the chosen witness difference. If the resulting subfamily
still has the old rank, use a nonzero basis minor \(\Delta\)
from (AF2). Integer Cramer's rule expresses \(\Delta\) times
each witness row as an integer combination of the locked rows.
By linearity, the lock divides \(\Delta\) times each witness
letter value and hence their difference times \(\Delta\).
Coprimality forces the lock to divide \(\Delta\), impossible
by its positive weight. The rank drops. If stripping itself
lowered the rank, that drop requires no block.

If every moving letter has purely atomic residual spectrum,
keep shrinking windows around finitely many positive atoms and
at zero. Fix the omitted tails at arbitrarily small root cost,
including any escaped numerical moment. Group the positive
windows at their representative valuation, and put exponent
corrections and the zero window into a rational factor of
height \(X^{o(1)}\). The numerical weight of \(F\) is at least
\(\eta-o(1)\) because \(w_G\le h_V+o(1)\). Discarding its
original positions of sufficiently small scaled degree leaves
a full modulus of positive weight and fullness of order \(1/L\).
The other two letters are coprime to this modulus by primitivity.
Their signed ratio is a fixed unit modulo it. The amplified
[bounded-root count](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-bounded-number-of-small-roots)
releases their grouped variables at arbitrarily small normalized
cost, with fixed coefficient factors absorbed into the twist.
Every exact ratio determines the whole primitive triple. Constant
letters simply have no variables to release.

There are at most \(R\) drops before rank one, when the constant
coordinate makes every complete letter constant. Use a cost less
than \(\epsilon/(4(R+1))\) at each drop and in the atomic terminal.
All reboxing and support losses can be taken arbitrarily small
at fixed truncation in the inner limit. This contradicts a
subfamily above exponent \(\epsilon\) on the diagonal sequence
and proves the stated uniform positive width.



## Rank and root-coordinate size

Suppose two original letters have wholly nonatomic limiting
spectra. Allow \(R=R(L)\) to grow in the outer limit, keeping
it fixed at each inner height limit. Put
\(V_0=\sum_i\log_n\max(1,\sup |z_i|)\).
For a putative fibre exponent \(r>0\), choose a compact positive
degree band \([s,S]\) on both diffuse letters, leaving outside
mass less than \(r/16\) on either. Let \(\omega_L\) be the largest
of their original root-coordinate caps in that band. A sufficient
condition to exclude an exponent \(r\) is
\[
                 R(\delta+LV_0+\omega_L)\longrightarrow0 .
 \tag{AF3}
\]
The argument also works if the left side is bounded by a
sufficiently small positive constant depending on \(r,s\).

Indeed allocate \(r/(4(R+1))\) root cost to each drop, starting
with a family of exponent at least \(r\). While the remaining
count is at least \(n^{r/2}\), the residual root caps of the two
wholly diffuse letters total at least \(r/2-o_X(1)\), since their
complete values determine the triple. This is true after
fixing all current witness support parts, at subpower cost.
Use zero cap for a constant letter. Thus a moving letter
has residual compact-band cap at least \(r/8\). Choose a block
of cap between \(r/(8(R+1))\) and that quantity plus \(\omega_L\).
Its full numerical weight is at least \(sr/(8(R+1))\) up to inner
tolerances, above the minor weight in (AF2) by (AF3). The
overshoot fits the allocation by (AF3). The rank drops while
total costs stay below \(r/2\), so the process cannot reach rank
one. On a generalized progression of parameter volume \(n^v\)
we have \(V_0\le v\). In particular for \(\delta=0\),
\(v=O(1)\), and \(\omega_L=O(L)\), ranks up to a sufficiently
small multiple of \(1/L\) admit a normalized saving.

## Lifting short rational features from two congruences

Work in the iterated balanced endpoint convention, \(n=X^L\), with signed complete letters \(T=U+V\), fixed coprime full divisors \(F\mid U,\ G\mid V\), and numerical height caps \(h_U,h_V\). A negative sign can be absorbed into the signed letter or cofactor. Write \(U=Fu,V=Gv\). Assume the numerical weight of \(FG\) exceeds \(\max(h_U,h_V)\) by a fixed margin.

Suppose each retained edge is assigned two integer vectors \(z_F,z_G\) and positive integer denominators \(d_F,d_G\) with a fixed number of coordinates and entries of magnitude \(\le n^C\) for fixed \(C\). Suppose
\[
 d_FT\equiv a_F\cdot z_F\pmod F,\qquad
 d_GT\equiv a_G\cdot z_G\pmod G
\]
for common integer coefficient vectors \(a_F,a_G\) of arbitrary size. Unit denominator hypotheses are unnecessary in this cleared-congruence formulation. If the relations are instead written as modular fractions, take \(d_F,d_G\) units in their own moduli.

Then the targets in such a subfamily have subpower count in the balanced iterated limit. They need no common Archimedean approximation and the residue coefficients need no height bound.

Indeed, choose independent \(z_F\) columns spanning the integer feature matrix over \(\mathbb Q\). For any augmented full minor using the column \(d_Fv\), multiplication by \(G\) and use of the congruence replaces the latter modulo \(F\) by \(a_F\cdot z_F\). That implies divisibility by \(F\), since the augmented determinant with the full linear combination is identically zero and \(\gcd(F,G)=1\). Its archimedean exponent is at most \(h_V-\log_XG+O(L)\), smaller than \(\log_XF\). Thus \(d_Fv\) is exactly a rational linear function of \(z_F\) on the subfamily. The other side gives \(d_Gu\) exactly as a rational linear function of \(z_G\). There is no reliance on clearing denominators of dependencies modulo \(F\) or \(G\); determinant multilinearity with the original integer \(a\)'s is enough. The argument permits rank zero on either side.

Include \(d=d_Fd_G\) and all coordinates of \(d_Gz_F,d_Fz_G\) in the new integer vector \(z\), so \(|z_i|\le n^{2C}\). Then each of \(dU,dV,dT\) is a common rational linear function of \(z\), and \(d\ne0\). Every bounded-rank minor of this new feature matrix has numerical exponent \(O(L)=o(1)\). The rational form of the rank argument gives subpower occupancy. For clarity, on a current subfamily select spanning witness rows. A moving letter \(W\) differs on two witnesses \(i,j\), since equality on a spanning set would imply the identity of its numerator with a common constant times the denominator coordinate \(d\) on the whole space. Strip and fix the full parts of \(W\) on the support of the nonzero integer \(D_W=d_id_j(W_i-W_j)\), at subpower cost by the fixed-support count. For any arbitrarily cheap full lock \(P\mid W\) on retained rows with positive residual diffuse weight, the lock is coprime to \(D_W\). If rank did not drop, take a nonzero fixed-coordinate minor \(\Delta\) of the retained full rank. By Cramer's rule, \(P\) divides each \(\Delta d_iW_i\), hence also \(\Delta D_W\). This contradicts the subpower numerical bound on \(\Delta\) because \(P\) has fixed positive numerical weight. When no moving letter has any positive residual diffuse mass, the final atomic argument from that proof applies unchanged: a positive numerical full modulus among the original fixed divisors exists by the fixed margin, and the grouped congruence releases the opposite two now purely atomic residual letters, after freezing arbitrarily small remainder and escape root masses. Fixed letters are constants throughout. Allocate arbitrarily small positive cost to each of at most the feature rank number of rank drops. The same diagonal argument as in the exact feature proof supplies uniformity over all coefficients and bounded feature prescriptions.

The one-sided statement needs only that \(FG\) exceed the height of the corresponding free letter by a fixed margin. For example, from the \(F\) congruence one obtains a rational-feature description of \(V\) alone, using \(d_F\) and \(z_F\). If \(V\) is wholly nonatomic, its distinct values have subpower size in this record. Apply only the one-letter rank descent; if the residual root mass vanishes, enumerate the remaining roots and stop. A bounded target dictionary is not required for that conclusion, and no factor constraints on \(U,T\) enter its proof.

For a two-sided proof in the balanced case there is also a short CRT route to the exact rational-feature description: multiply each congruence by its complementary denominator, combine them modulo \(FG\) with a single CRT coefficient for each feature, and choose the coefficients at most \(FG\) in magnitude. Since \(FG\) also exceeds the target's full height by a fixed margin at balance, the added integer quotient by \(FG\) is \(n^{O(1)}\), an extra feature. The determinant route above explains the one-sided assertion and remains valid without the CRT route's auxiliary height comparison on \(T\).

The input is an economical residue feature model, not a deduction from the existence of the records. It can accommodate different modular feature dictionaries at the two moduli.

## Consequences for the interval count

Suppose each affordable pair record of normalized cost
\(t\le 1/3-\eta_0\) can be covered by at most \(n^\rho\)
systems (AF1), with uniform bounds on \(R,C\) and
\(\rho<\eta_0\) by a fixed margin. With
\(\epsilon<\eta_0-\rho\), summing over records gives exponent
at most \(1/3-\eta_0+\rho+\epsilon+o(1)\) on this sector.
The missing step is to obtain sufficiently economical covers
or otherwise bound the separated targets of arbitrary records.

For any \(r>0\) and \(K_0<\infty\), a set of at least \(n^r\)
distinct exact targets on one affordable record cannot have
doubling bounded by \(K_0\) along the iterated endpoint limit.
Freiman's theorem would contain it in a proper generalized
progression of bounded rank and volume \(O_{K_0}(|S|)\).
Its coordinates give (AF1) with zero error. This extends the
previous wholly diffuse conclusion to all balanced profiles.
For the width version, use a subset of the occupied indices
\(\lfloor T/B\rfloor\), where \(B=\lfloor X^\delta\rfloor\).
A progression for the indices approximates the actual targets
with error at most \(B\). The positive width exponent is chosen
after fixing \(r,K_0\). The Balog–Szemerédi–Gowers theorem and
Plünnecke give the corresponding exclusion for energy above any
fixed positive fraction of cubic cardinality, choosing the width
after that fraction. These arguments give no power rate.
[Green–Ruzsa, Theorem 1.1](https://arxiv.org/abs/math/0505198)
· [Reiher–Schoen, Theorem 1.1](https://arxiv.org/abs/2308.10245)
· [Petridis, Theorem 1.1](https://arxiv.org/abs/1101.3507)

Croot–Mao–Yip obtain progression covers from structured
reductions of a whole set at many test primes. A fixed
polynomial-height integer has \(O(\log X/\log y)\) prime
divisors in \([y,2y]\), whereas their Theorem 1.1 needs
a multiple of \(y^{1-2\varepsilon}\log X/\log y\).
In particular at \(y=n^\beta\), \(\beta>0\) fixed, the primes
of one affordable record cannot alone supply these tests.
The high powers on individual targets give divisibilities of
different points, leaving the common structured reductions
needed for that inverse theorem unproved.
[Croot–Mao–Yip, Theorems 1.1 and 1.10–1.12](https://arxiv.org/html/2510.08862v2#S1.SS2)



## High energy despite polylogarithmic intersections

Integer sets can have polylogarithmic intersections with all
the systems (AF1) and still have additive energy with the full
cubic exponent. Fix \(r>0,R,C\), and take \(L>0\) sufficiently
small relative to them. Set
\[
 K=\lceil\log X\rceil,\qquad
 m=\left\lfloor\frac{r\log n}{\log K}\right\rfloor,\qquad
 Y=\lfloor X^{1/4}\rfloor,\qquad
 \mathcal D=\{0,\ldots,K-1\}^m .
\]
Choose \(b_i\) independently and uniformly from \(1,\ldots,Y\).
For any nonzero integral vector \(\boldsymbol d\), the
probability of \(\boldsymbol b\cdot\boldsymbol d=0\) is at
most \(1/Y\) by conditioning on all coordinates but one.
Exclude collisions and nonzero equal-sum differences of
four digit vectors. Also exclude all nonzero combinations
of at most \(R+1\) digit vectors with integer coefficients
bounded by \((R+2)!n^{CR}\). There are at most
\(O_R(n^{2r}+n^{4r}+n^{(R+1)(r+CR)})\) vectors to exclude.
The failure probability tends to zero at fixed \(L\) if
\(L\max(4r,(R+1)(r+CR))<1/4\).
Fix a successful choice and put
\[
 S=\{\lfloor X/2\rfloor+Y\,\boldsymbol b\cdot\boldsymbol a:
                                                \boldsymbol a\in\mathcal D\}.
\]
All values are \(\asymp X\), and their mutual spacing is at
least \(Y\). The map preserves and reflects equal-sum
relations, so
\[
 |S|=K^m=n^{r+o_X(1)},\qquad
 E_+(S)=\bigl((2K^3+K)/3\bigr)^m=n^{3r-o_X(1)}.
\]
The factor \((2/3+O(K^{-2}))^m\) is subpower in the inner
limit and tends to zero.

For a subfamily of \(S\) in (AF1) with \(\delta<1/8\), take
\(k\) independent feature rows. Each other row has an
integer dependence on them with coefficients bounded by
\((R+2)!n^{CR}\), by Cramer's rule. Since the coefficients
sum to zero by the constant coordinate, their combination
of targets is a multiple of \(Y\). The fit (AF1) bounds
its magnitude by \(O_R(n^{CR}X^\delta)<Y\), so it is zero.
The choice of weights implies the same relation in the
digit vectors. These thus lie in the affine span of the
\(k\) witness digits and are determined by at most \(k-1\)
pivot coordinates. The intersection has at most \(K^{R-1}\)
points. The construction uses exact integer addition and
full numerical height. It imposes no target radical or
valuation condition.


## Amplifying full target locks at a fixed shift

Let \(n=X^L\) and take the inner height limit before the small-budget limit. Fix positive integers \(V,F\) of polynomial height in each box, with \(F\mid U\) and \(T=\epsilon V+\omega U\) for fixed signs. The edges are primitive and the complete \(U,T\) are positive. The factor \(F\) is a full record from the \(U\) factorization, so the remaining root box for \(U/F\) has normalized root cap \(a\) and numerical height cap \(H\), in bases \(n,X\) respectively. Assume that this residual factorization has a compact positive degree block allowing a selection of whole positions of normalized cap \(p+o(1)\) and numerical height at most \(S_U p+o(1)\). Such selections exist for every sufficiently small fixed \(p>0\) when a positive residual nonatomic mass is supported in a compact positive band with upper endpoint \(S_U\). With atomic spectra, the proof still applies to any chosen all-row common block of the stated cap and height.

Suppose \(T\) has a band of scaled degrees in \([s,S]\) permitting the selection of whole positions of total root cap \(\beta+o(1)\) and numerical weight at least \(s\beta-o(1)\). For example \(\beta\) can be chosen below any positive nonatomic band mass by omitting shrinking windows of atoms, then selecting whole positions with no overshoot in the limit. Keep also an optional disjoint high tail of root cap \(\zeta+o(1)\) on one row, of weight at least \(e-o(1)\). This can carry any escaped numerical moment \(e\) using arbitrary small fixed \(\zeta>0\) by passing to a sufficiently far continuity cutoff and using actual rooted lower exponents. Take \(e=\zeta=0\) when no tail is used.

For a putative normalized fibre count exponent \(r=a-\delta\), fix an integer \(m\ge2\). The following three strict inequalities are incompatible with that count:
\[
        \delta<p,\qquad
        m\beta+\zeta<a-p-m\delta,\qquad
        e+s m\beta-\binom m2 S_U p>H .                              \tag{FS1}
\]

To prove this write \(U=Fxy\), with \(x\) the selected root block of cap \(p\). Make the simple bipartite graph of the \(x,y\) values. Each solution gives a distinct edge after selecting one factorization per solution, since \(xy\) and the fixed data determine the complete letters. The two caps for graph vertices are \(n^{p+o(1)},n^{a-p+o(1)}\). Choose \(0<\eta<p-\delta\) and discard columns of degree below \(n^{p-\delta-\eta}\). This loses at most \(n^{a-\delta-\eta+o(1)}\) edges. The surviving columns have degrees tending to infinity, so distinct ordered \(m\)-stars count \((1-o(1))\sum_y d_y^m\), where \(d_y\) is the degree at \(y\). By the power mean that sum is at least \(n^{mp+a-p-m\delta+o(1)}\). Thus some distinct fixed \(x_1,\ldots,x_m\) have at least \(n^{a-p-m\delta+o(1)}\) common \(y\) values.

The target in row \(i\) is \(T_i(y)=\epsilon V+\omega Fx_i y\). Record its chosen band roots on each row, and the tail on row one if used. Let \(M_i\mid T_i(y)\) be the product of the complete powers recorded on row \(i\); the band and tail supports are disjoint on the same row. Across rows the locks need not be pairwise coprime. Rather, for \(i\ne j\) any prime of \(\gcd(M_i,M_j)\) is prime to \(Fy\) by primitivity and \(Fy\mid U_i\). From
\[
  T_i-T_j=\omega Fy(x_i-x_j)
\]
it follows with full multiplicities that \(\gcd(M_i,M_j)\mid x_i-x_j\). Thus its numerical weight is at most \(S_U p+o(1)\). For every prime, the sum of its exponents in all \(M_i\) minus their maximum is at most the sum of the pairwise minima. Consequently the least common multiple \(M\) of the recorded \(M_i\) has numerical weight at least
\[
       e+s m\beta-\binom m2 S_U p-o(1).
\]
For each row, \(Fx_i\) is a unit modulo its recorded modulus, again by primitivity. The row congruences specify at most one residue class for \(y\) modulo \(M\): any two solutions agree modulo every \(M_i\). Since \(y\le U/F\le X^{H+o(1)}\), the third inequality gives one \(y\) per record. There are at most \(n^{m\beta+\zeta+o(1)}\) possible records, contradicting the second inequality. All bounds are uniform in \(F,V\) and the fixed row choices. They use the numerical height of \(x_i-x_j\), so large overlaps in the radical are permissible. Root and tail bounds here come from one common dyadic or refined logarithmic root bin for each position; if prelimit indices are unbounded, finer intervals whose ratio is \(1+o(1)\) still leave only \(X^{o_X(1)}\) bins over active roots and give numerical error \(o(1)\) in base \(X\).

In particular if the fixed-shift fibre saturates its residual root bound (\(\delta=0\)) and the residual \(U\) letter has positive compact diffuse mass, every positive nonatomic band of \(T\) with lower endpoint \(s\) must satisfy \(H\ge e+s a\). Otherwise choose \(q<a\) with \(e+s q>H\), choose fixed \(m\ge2\) with \(\beta=q/m\) below the band mass, then choose \(\zeta,p\) small to satisfy (FS1). The margins persist for sufficiently small positive \(\delta\). The result can also be conditioned on an \(n^{o(1)}\) dictionary for \(V\) in a fixed pair record. In that application if the original endpoint count matches the unsaved fixed-pair exponent and the \(G\) root cost is \(g\), the residual \(U\) deficit is \(\delta=g+o(1)\), so check (FS1) quantitatively.

For an exact comparison with just locking a single \(T\) value's factors, set \(a=1/12,\ H=17/60,\ S_U=7/2\), for instance from a residual \(U\) law uniform on \([33/10,7/2]\) with mass \(a\). Let the \(T\) law have mass \(1/300\) uniformly on \([4,41/10]\) and mass \(33/100\) uniformly on the interval of radius \(1/5\) centred at \(1973/660\). Its total mass is \(1/3\) and moment is one. The low interval has upper endpoint below \(16/5\). At any affordable single \(T\) lock of root mass at most \(a\) its numerical weight is no greater than
\[
 \frac1{300}\frac{41}{10}+\left(a-\frac1{300}\right)\frac{16}{5}
 =\frac{809}{3000}<\frac{17}{60}.
\]
Here the row method works with \(m=24,\ s=4,\ \beta=1/320,\ p=1/120000,\ e=\zeta=0\), and \(\delta=1/12000000\). The total recorded band cap is \(3/40\) and the overlap penalty is \(966/120000\). All three inequalities (FS1) are strict. This formal profile demonstrates scope beyond single-target CRT with the given fixed shift and residual box. It does not assert the existence of overflowing fibres. The displayed theorem does not require assumptions about target root counts versus the residual \(U\) cap.

### A uniform gap for comparable residual root caps

Fix \(0<a_0\le a_1<\infty\) and \(H_0,w_0,B_0>0\). Consider fixed-shift primitive pairs \(U,T\) as above, now with full records \(F\mid U\) and \(M\mid T\) both fixed (possibly unit). Assume the residual root boxes of \(U/F\) and \(T/M\) have normalized root caps \(a_U,a_T\in[a_0,a_1]\), numerical heights at most \(H_0\), with prelimit root moments giving the residual numerical heights up to vanishing precision errors. Assume one of \(F,M,V\) has numerical weight at least \(w_0\) and radical at most \(n^{B_0+o(1)}\). There are \(\sigma,\gamma,L_0>0\), depending only on the indicated bounds, such that, for \(L<L_0\) and \(|a_U-a_T|<\gamma\), the normalized fixed-record count is at most
\[
           n^{\min(a_U,a_T)-\sigma+o_X(1)}.
\]
Here the constants and \(o_X(1)\) use the same fixed-budget and fixed-factorization conventions as the surrounding root-box arguments. All numerical data of the shift and records can vary with \(X\) within the given polynomial height assumptions. The count fixes the numerical shift and the product values of the complete recorded factors, not only their supports.

Otherwise choose a diagonal countersequence with \(L\to0\), root caps tending to the same \(a\in[a_0,a_1]\), and normalized count tending to \(a\). Pass to weak limits of the two residual degree measures \(\nu_U,\nu_T\) of mass \(a\), and limiting residual heights \(H_U,H_T\), with escaped moments \(e_U,e_T\). Tightness follows from the bounded numerical first moments. In the fixed integer supplying the hypothesis, retain the
prime powers \(p^k\) with \(kL\ge\kappa=w_0/(4B_0)\).
The discarded weight is at most \(w_0/3\) with sufficiently
small tolerances by the radical bound. The remaining divisor
has positive fixed numerical weight and fullness of order
\(1/L\). If it comes from a free-letter record, exchange the
letters so that this record is \(F\).

When the modulus comes from \(F\mid U\), primitivity gives
a unit congruence for \(T/M\) under it. Indeed the fixed shift and \(M\) are units modulo \(F\). Freeze the residual \(T\) roots outside shrinking windows around finitely many positive atoms and the zero atom, with upper-tail root mass as small as desired. As in the affine rank proof, group the positive windows, with their narrow products as the small variables and their degree discrepancies, and the zero window, as a single \(X^{o(1)}\) correction. Under the given positive full modulus, amplified (GW) releases the exact residual value at arbitrarily small normalized counting cost. This uses positive full weight and \(k_0L\) bounded below, uniform in the target unit and fixed coefficients. If no positive windows remain, use direct rational reconstruction for the small correction. Taking any finite atom list first, then increasing that list, thus gives \(a\le \nu_T^{\rm na}(\mathbb R_{\ge0})\). So \(\nu_T\) is wholly nonatomic.

This gives a compact positive band of \(T/M\) of positive diffuse root mass. It supplies a further complete lock of arbitrarily small positive root cost, and fixed positive numerical weight with normalized fullness bounded below. Freeze it, impose the unit congruence on \(U/F\), and run the same argument. Taking the new lock cost arbitrarily small shows that \(\nu_U\) is wholly nonatomic. All locked factors of the opposite letter are units by primitivity, regardless of whether the original \(M\) also had positive numerical weight.

When the original positive modulus comes from the fixed
shift \(V\), there are two cases. If either residual measure
has positive diffuse mass, its cheap full band supplies a
modulus on that letter. The preceding argument forces the
opposite residual measure to be wholly diffuse. Apply it
back to obtain two wholly diffuse measures. If both are
purely atomic, group the residual windows on both sides.
Their exact ratio is a unit with fixed residue modulo the
positive full divisor of \(V\), including the fixed factors
\(F,M\) in the twist. The amplified count releases their
ratio and hence both complete free values at arbitrarily
small normalized cost, contradicting saturation.

Apply (FS1) with \(U/F\) as residual varied side and the factors of \(T/M\) as the target bands and optional tail, ignoring the already fixed \(M\) in the lcm count. Then apply it in the opposite orientation. If \(s_T\) is the essential supremum of the support of \(\nu_T\), the saturation consequence of (FS1) gives \(H_U\ge e_T+a s_T\), by taking positive compact nonatomic bands approaching \(s_T\). In particular \(s_T<\infty\), and likewise \(s_U<\infty\). Nonatomicity and positive mass give
\[
        H_U\ge e_T+a s_T>e_T+\int r\,d\nu_T(r)=H_T .
\]
Interchanging \(U,T\) gives \(H_T>H_U\), a contradiction. This proves the uniform strict gap by the diagonal argument. The modulus hypothesis matters for this derivation of wholly nonatomic residual profiles. Merely fixing a large integer with no positive weighted power divisor is insufficient for the amplified step.

Several of the same reconstructions work with a varying shift
dictionary when multiplicative rows have enough common columns.
The interval version uses the dictionary's numerical width.
[Common rows with varying shifts](ABC_AGGREGATE_SHIFT_DICTIONARY_STARS.md)

## A modular progression cover implies affine sparsity

Fix \(r\ge1\), \(0<\theta<1\) and
\(0<\varepsilon_0<(1-\theta)/(r+1)\). Put
\(\beta=\theta-(1-\theta)/(r+1)+\varepsilon_0\). In the
setting of (AF1), let \(S\subseteq[1,X]\) be the distinct
targets of the subfamily. At fixed \(L\), take
\[
 (8^{r+3}\log X)^{2/(r\varepsilon_0)}<y\le n^C
\]
with \(C\) fixed and \(y\le X\). Suppose that for each prime
in \([y,2y]\), reduction of \(S\) is contained in a
generalized arithmetic progression of rank at most \(r\)
and cardinality at most \(p^\theta\). In the iterated
small-budget limit the count is at most
\[
                        y^{\max(\beta,0)}+n^{o(1)}.          \tag{AF6}
\]
Here the \(o(1)\) bound expresses an arbitrarily small
fixed normalized allowance for sufficiently small \(L\).

If the threshold \(|S|\ge y^\beta\) is met, Croot–Mao–Yip,
Theorem 1.12, covers \(S\) by a bounded number of translates
of a rank-at-most-\(r\) generalized arithmetic progression
in \(\mathbb Z\) of cardinality \(O_r(y)\). Discard zero
generators. In any generalized arithmetic progression in
\(\mathbb Z\), the number of allowed coefficients of one
nonzero generator is at most the cardinality of the
progression, since all values along that generator alone
are distinct. Thus the integer coefficients are bounded by
\(n^{C+o(1)}\) with no hypothesis on the integer steps or
on uniqueness of representations. Use these coefficients
and the constant feature in (AF1) with zero error on each
translate. Below the threshold the first term of (AF6)
already bounds the size.
[Croot–Mao–Yip, Theorem 1.12](https://arxiv.org/html/2510.08862v2#S1.SS2)

A radical cap for one complete letter cannot furnish this
modular hypothesis by itself. At fixed \(k\ge2\), the set
of \(k\)-th powers in \([X/2,X]\) has radicals at most
\(X^{1/k}\). Its roots occupy a consecutive interval of
length \(\gg_k X^{1/k}\). Repeated residues of these roots
then give the entire set of \(k\)-th powers modulo every
prime \(p\) below this length, with cardinality
\(1+(p-1)/\gcd(k,p-1)\gg_k p\). This example concerns one
letter. The missing input for using (AF6) on actual fibres
is a progression restriction valid at many independent
primes from the three-letter relation.
