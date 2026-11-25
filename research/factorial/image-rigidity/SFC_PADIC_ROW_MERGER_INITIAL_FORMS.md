# A discrete row merger that preserves factorial moment initial forms

Fix a prime \(p\), a positive integer order bound \(J\), and any
finite support of distinct ordinary exponent columns. A pair of
exponent rows can be replaced by one ordinary row, using a power
of \(p\) to separate their digits, so that all moment polynomials
of orders at most \(J\) retain the same weighted initial forms
after one diagonal rescaling of the coefficient variables. This
applies to arbitrary support size.

Let the exponent entries in the two rows be \(A_i,B_i\),
\(1\le i\le m\), and replace them by \(p^e A_i+B_i\), where
\[
                    p^e>J\max_i B_i .
\]
One may choose \(e\ge1\). This map is injective on the distinct
pairs \((A_i,B_i)\) because \(0\le B_i<p^e\). Every support column
therefore remains distinct once the other rows are included.
Write \(M_j(c)=\mathcal L((\sum_i c_i x^{a_i})^j)\) for the
original moment and \(M'_j\) for the merged one. Put
\[
 D_e=1+p+\cdots+p^{e-1},\qquad
                   \lambda_i=(-p)^{D_e A_i}.
\]

**Theorem.** For \(1\le j\le J\), each coefficient of
\(M'_j(\lambda_1^{-1}c_1,\ldots,\lambda_m^{-1}c_m)\) is the
corresponding coefficient of \(M_j\) times a rational \(p\)-adic
unit congruent to one modulo \(p\). All orders use the same
diagonal rescaling. The same conclusion with modulus \(p^r\)
holds if each \(A_i\) is divisible by \(p^{r-1}\),
\(p^{e-r+1}>J\max_i B_i\), and the diagonal is chosen as below.

## The factorial comparison

For a coefficient indexed by \(\nu=(\nu_1,\ldots,\nu_m)\)
with \(\sum_i\nu_i=j\le J\), put \(x=\sum_i\nu_iA_i\) and
\(y=\sum_i\nu_iB_i\). Only the two row factors change. Their
factorial quotient is
\[
 \frac{(p^e x+y)!}{x!y!}
 =\frac{(p^e x)!}{x!}\prod_{h=1}^y(1+p^e x/h).
\]
Since \(0\le y<p^e\), each factor in the product is a \(p\)-adic
unit congruent to one modulo \(p\). Under the sharper restriction
on \(e\) in the theorem it is congruent to one modulo \(p^r\):
each index \(h\) in the product satisfies \(v_p(h)\le e-r\).

Separate the factors divisible by \(p\) in \((p^s x)!\) successively,
for \(s=e,\ldots,1\). This gives the exact identity
\[
 \frac{(p^e x)!}{p^{D_e x}x!}
   =\prod_{s=1}^e\prod_{\substack{1\le h\le p^s x\\p\nmid h}} h .
\]
Each complete block of length \(p\) contributes \(-1\) modulo
\(p\) by Wilson's theorem. There are \(D_e x\) such blocks.
Thus the original factorial quotient is \((-p)^{D_e x}\) times
a \(p\)-adic unit congruent to one. The scalar is exactly
\(\prod_i\lambda_i^{\nu_i}\), proving the claim.

For modulus \(p^r\), let \(\epsilon\) be an integer representative
of the product of all units modulo \(p^r\). Each \(p^s x\) is
divisible by \(p^r\), so the last display is congruent to
\(\epsilon^{D_e x/p^{r-1}}\) by partition into complete blocks.
Use
\(\lambda_i=p^{D_e A_i}\epsilon^{D_e A_i/p^{r-1}}\).
If \(y=0\) the shorter factorial contributes the empty product,
so no additional restriction on \(e\) is needed for that term.
Taking the second row identically zero specializes the comparison
to rowwise inflation by \(p^e\).

## Preservation of root-exclusion certificates

Choose a nonempty set of coefficient indices on which the
coefficient vector is allowed to be nonzero. Set all other
coordinates to zero. On that coordinate torus choose any rational
weight vector \(w\). At each order \(j\), the weighted initial form of
\(M_j=\sum_\nu a_\nu c^\nu\) consists of the terms whose values
\(v_p(a_\nu)+w\cdot\nu\) are smallest, each with coefficient
\(\overline{p^{-v_p(a_\nu)}a_\nu}\) in \(\mathbb F_p\).
The factorial comparison preserves both the values and these units after
rescaling. It thus preserves all these forms simultaneously,
including on each coordinate torus and for each selected order.
Before rescaling, the valuations correspond to translating
\(w\) by the same vector \((D_e A_i)_i\) at every order, and
the unit parts differ by the same diagonal signs.

For a finite set of moment orders, a certificate may exclude
common zeros of the initial forms over
\(\overline{\mathbb F}_p\) for every coordinate torus and every
weight. This implies absence of common projective zeros in
characteristic zero. In fact, if such a zero exists over a
characteristic-zero field, the rational polynomial system also
has an algebraic zero over \(\overline{\mathbb Q}\), by the
Nullstellensatz on a nonempty projective chart. Embed its
coordinate field into \(\overline{\mathbb Q}_p\). Use the
valuations of the nonzero coordinates as \(w\). After adjoining
a suitable root of \(p\), divide each coordinate by the
corresponding power of that root and reduce the resulting
units. In a vanishing moment sum, divide by the smallest
term valuation before reducing. Only the minimizing terms
survive and give a common residue zero of the initial forms.
For fixed support and orders, the list of simultaneous initial
systems is finite, since only finitely many subsets can minimize.

The row merger preserves every certificate of this kind in both
directions. Iterating the merger packs any ordinary finite
support into a single exponent row for the whole finite order
window. In the binary case the theorem preserves the Newton
polygons and all residual polynomials, up to a common shift of
slopes and a unit change of the coefficient ratio. In arbitrary
support it preserves the corresponding systems on all coordinate
tori. This argument transfers the exclusion certificate because
every characteristic-zero zero would give a residue zero. A
residue zero by itself gives no conclusion about lifting to an
exact common zero.
