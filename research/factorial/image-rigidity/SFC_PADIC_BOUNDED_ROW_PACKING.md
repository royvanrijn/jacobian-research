# Bounded-dimensional factorial packing at any fixed prime precision

Fix a prime \(p\), a precision \(\ell\ge1\), a support size \(s\ge2\),
and a finite moment range \(1\le j\le J\). An ordinary profile
is a finite list of nonnegative integer exponent rows with
\(s\) columns. The columns are distinct. Write
\[
 \mathcal M_j(c_1,\ldots,c_s)
  =\sum_{|\nu|=j}\binom{j}{\nu}
        \prod_i\left(\sum_{a=1}^s\nu_a a_{i,a}\right)!\,c^\nu .
\]
The theorem concerns all coefficients and uses the same
positive rational diagonal change of coefficient variables
at every order.

**Theorem.** There is an ordinary profile with at most
\[
         1+(2p^{\ell-1})^s\bigl(\varphi(p^\ell)-1\bigr)        \tag{1}
\]
rows whose moment coefficients, after such a diagonal
change, have exactly the original valuations and exactly
the original unit parts modulo \(p^\ell\) for \(1\le j\le J\).
The row bound is independent of \(J\) and of the original
dimension and exponents. The construction preserves
distinctness of the support columns. At \(\ell=1\), one
row suffices. In the binary case the result preserves
all exact common-argument truncated Newton face comparisons
at the specified precision.

This is a finite-prime reduction. The exponents of the
packed row and the scaling can grow with \(J\) and with
the original profile. The conclusion concerns valuations
and units at the chosen prime.

## The finite defect of inflating a row

Let \(P=p^\ell\) and put
\[
 F_p(t)=\prod_{\substack{1\le u\le t\\p\nmid u}}u,\qquad
 \varepsilon=F_p(P)\bmod P,\qquad D_e=1+p+\cdots+p^{e-1}.
\]
The product of all elements of \((\mathbb Z/P)^\times\)
has square one, by inversion, so \(\varepsilon^2=1\).
The [factorial inflation identity](SFC_PADIC_ROW_MERGER_INITIAL_FORMS.md#the-factorial-comparison)
obtained by successively removing the \(p\)-divisible
factors gives
\[
   \frac{(p^e x)!}{p^{D_e x}x!}
                     =\prod_{h=1}^e F_p(p^h x) .
\]
For \(e\ge\ell\), the terms with \(h\ge\ell\) break into
complete unit blocks of length \(P\). Hence the unit on the
left is congruent to
\[
         \varepsilon^{D_{e-\ell+1}x}\,B_\ell(x)\pmod P,\qquad
                B_\ell(x)=\prod_{h=1}^{\ell-1}F_p(p^h x)
                                      \pmod P.               \tag{2}
\]
The defect \(B_\ell(x)\) is periodic in \(x\) with period
\(2p^{\ell-1}\). In fact increasing \(x\) by that amount
adds \(2p^{h-1}\) complete unit blocks in the \(h\)-th
factor of the defect, each contributing \(\varepsilon\).
All these block counts are even. The empty defect for
\(\ell=1\) is one. This uses the elementary prime-power
block product also underlying
[Granville's factorial unit formula](https://www.sas.rochester.edu/mth/sites/doug-ravenel/otherpapers/granville.pdf).

Group the rows by their \(s\) residues modulo
\(2p^{\ell-1}\). For a fixed multi-index \(\nu\), rows
in the same class give the same defect for their
factorial arguments \(x_i=\sum_a\nu_a a_{i,a}\).
In each class retain only the remainder of its
multiplicity modulo \(\varphi(P)\). All rows removed
from their classes have product of defects equal to
one at every \(\nu\), by Lagrange's theorem in
\((\mathbb Z/P)^\times\). The same set of rows can
therefore be packed for every moment order.
At \(\ell=1\) remove all rows.

This neutralization is needed at higher precision.
For example, modulo \(4\) the single argument defect
\(B_2(x)=F_2(2x)\) is one at \(x=1\) and three at \(x=2\).
It cannot be absorbed into one unit to the power \(x\)
at these two values. This occurs for the row \((0,1)\)
when the factorial argument ranges from zero to two.

## Packing the neutral rows

If the neutral collection is nonempty, replace its
rows by one row with entries \(\sum_i p^{e_i}a_{i,a}\)
in each column, taking \(e_i\ge\ell\) successively
large. When inserting row \(i\), ensure that the
largest possible preceding lower argument \(y\),
over \(|\nu|\le J\), is smaller than \(p^{e_i-\ell+1}\).
Then
\[
       \frac{(p^{e_i}x_i+y)!}{(p^{e_i}x_i)!\,y!}
              =\prod_{u=1}^y(1+p^{e_i}x_i/u)
                    \equiv1\pmod P
\]
as a \(p\)-adic unit. The condition follows from
\(v_p(u)\le e_i-\ell\) for all \(u\le y\).
The integers \(e_i\) also ensure disjoint digit
blocks in each support column. If two columns
differed at a neutral row, they differ in the
packed row. Columns differing in a retained row
still differ there.

In the quotient of coefficients after packing,
(2) contributes just the common positive diagonal
scaling associated to the powers \(D_{e_i}x_i\) of
\(p\) and the powers \(D_{e_i-\ell+1}x_i\) of
\(\varepsilon\). More explicitly, take a positive
integer representative \(\varepsilon_0\) of the unit
\(\varepsilon\), and let
\[
 \lambda_a=p^{\sum_i D_{e_i}a_{i,a}}\,
        \varepsilon_0^{\sum_iD_{e_i-\ell+1}a_{i,a}}
\]
with sums over packed rows. The coefficient quotient
is \(\prod_a\lambda_a^{\nu_a}\) times a \(p\)-adic
unit congruent to one. Thus substituting
\(c_a\mapsto c_a/\lambda_a\) into the new moments
proves the assertion. At most
\((2p^{\ell-1})^s(\varphi(P)-1)\) rows were retained,
giving (1). If the neutral collection is empty,
no packed row is needed.

In the binary case dividing out the constant term
at each order reduces the diagonal change to the
common rational argument scaling
\(z\mapsto(\lambda_1/\lambda_2)z\).
After any common Newton substitution
\(z=\pi^a X,\ \pi^b=p\), and each parent's scalar
normalization, all coefficients above or on its
Newton face therefore agree modulo \(\pi^{b\ell}\).
Every comparison using that finite precision, in
particular the
[truncated face width test](MULTIVARIATE_PHYSICAL_TRUNCATED_NEWTON_FACES.md),
is unchanged.

## Matching the data inside arbitrarily weak binary supports

One can also make the packed binary profile fully shared
with arbitrarily small curvature. Add the same large
integer \(L\) to each entry in each resulting row. If
\(x,y\) are any two original factorial arguments in
that row at one order \(j\le J\), the ratio of the
new and old factorial quotients, for \(x\ge y\), is
\[
           \frac{(jL+x)!/(jL+y)!}{x!/y!}
                      =\prod_{h=y+1}^x(1+jL/h).
\]
For \(x<y\) use the reciprocal. All such factors
are units congruent to one modulo \(P\) once \(L\)
is divisible by a sufficiently large power of \(p\).
Thus the translated profile still reproduces the
original data after the same diagonal rescaling and
separate constant normalizations at each order. Its
gaps stay fixed, while their midpoint curvature
tends to zero. In particular, by the
[uniform weak-curvature theorem](MULTIVARIATE_UNIFORM_WEAK_CURVATURE_ANGLE_SEPARATION.md),
one can realize every such finite-prime record by an
all-degree regular binary profile with the row bound (1).
With no compression this translation argument applies
at any prescribed finite collection of primes by
taking \(L\) divisible by all the required powers.

## Finite recognition of the resulting screens

For binary profiles with \(p,\ell,J\) fixed, the set of
attainable valuations and units after removal of the
common affine \(p\)-slope is effectively definable
in Presburger arithmetic. This also holds for survival
of the common-degree screens from the truncated faces.
The ordinary exponent gaps and the number of rows
need not be bounded.

Extend the [first-layer carry graph](MULTIVARIATE_PADIC_DIGIT_AUTOMATON.md)
as follows. For \(x\ge0\) write \(U_p(x!)=p^{-v_p(x!)}x!\).
Iterating removal of the factors divisible by \(p\)
and then taking complete blocks modulo \(P\) gives
\[
 U_p(x!)\equiv
    \varepsilon^{A_\ell(x)}
       \prod_{t\ge0} F_p\!\left(\lfloor x/p^t\rfloor\bmod P\right)
           \pmod P,\qquad
 A_\ell(x)=\sum_{t\ge\ell}\lfloor x/p^t\rfloor .                \tag{3}
\]
Each argument of \(F_p\) in (3) is a window of
\(\ell\) consecutive digits, starting at position \(t\).
If \(x=\sum_i x_i p^i\), then
\(A_\ell(x)=\sum_{i\ge\ell}x_iD_{i-\ell+1}\).
Only the parity of this exponent is needed.
The weights \(D_m\bmod 2\) have the finite update
\(D_{m+1}=pD_m+1\bmod2\), preceded by the finitely
many positions \(i<\ell\) that have weight zero.

Read the digits of the two entries of a row from
least to most significant, using the finite carry
state for each output \((j-k)a+kb\).
Keep the last \(\ell-1\) digits of each output so that
every completed window in (3) contributes one finite
unit factor. The state also tracks the parity weight
at the current digit and the first \(\ell-1\) positions
where no full window has been read. After the output
terminates, flush the carries and pad by zeros to
complete all final windows. Reset the state before
the next row. Each transition records the same
integer digit-sum differences as in the first-layer
graph, the modulo-\(p^\ell\) unit contributions,
and the signed degree difference modulo \(p-1\).
Net units of the constant-normalized coefficients
are computed by multiplying the corresponding ratios
in (3) by the fixed unit part of \(\binom jk\).
The substitution \(z\mapsto p^{-(\Delta-\delta)/(p-1)}z\)
gives the centered digit-sum valuation formula of
the first-layer graph. Here \(\delta\) is the chosen
representative of \(\Delta\bmod(p-1)\).
This positive rational rescale leaves the
stripped units unchanged.

The weighted walks in this finite graph have a
Presburger description by the balance and connected
edge-support conditions of the first-layer proof.
The integer weights contribute linear sums, and
the finite unit group contributions contribute finite
congruences after choosing a decomposition of the
group. Ordinary distinct columns can be imposed by
requiring a digit mismatch at least once.
For a possible shared exponent valuation \(a/b\),
the edge endpoints and the reduction of its
denominator require only finitely many choices of
indices and congruences of valuations. The scaled
coefficient excesses
\(b\,v_k+ak-\min_h(b\,v_h+ah)\)
are affine functions once the minimizing index is
fixed. At precision \(\pi^{b\ell}\) each excess is
either at least \(b\ell\) or belongs to a finite list
of exact values. The unit information in (3) then
determines the truncated polynomials in the finite
ring \(\mathcal O_{\mathbb Q_p(\pi)}/\pi^{b\ell}\).
The resulting Sylvester width tests can be tabulated
on those finite choices. Linear comparisons and
congruences therefore express whether any profile
survives the chosen finite-precision screen.
The normal form and decidability for these linear
integer conditions can be taken from
[Ginsburg and Spanier](https://msp.org/pjm/1966/16-2/pjm-v16-n2-p09-p.pdf).

This reduction applies one prime at a time. Choosing
a single row packing or a single polynomial
certificate at all primes is a separate problem,
as is an exact common-root exclusion from the
finite-precision surviving data.

