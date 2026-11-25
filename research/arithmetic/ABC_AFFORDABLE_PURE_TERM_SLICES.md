# Marginal costs and polynomial expansion in affordable fibres

A direct transfer with only two mixed factors can have one term
that is a product of original free letters. The two possible
free-letter orientations cannot both occur in a cheap transfer
of this form. In the remaining orientation the pure term is
the denominator in a polynomial chart. Fixing its original
letter values brings the fixed-complete-letter count into play
and permits incidence bounds on a Cartesian grid of actual
fibre slices.

[The two-mixed-factor reduction](ABC_CROSS_RATIO_FRESH_PRIME_OBSTRUCTION.md#two-mixed-factors-and-a-pure-free-letter-term)

## Dyadic regularization of the two free-letter marginals

Write \(n=X^L\), take the height limit at fixed \(L>0\),
and let an affordable record contain \(N=n^{r+o_X(1)}\)
actual rows for some fixed \(r>0\). The signed relation
\(T=\sigma U+\tau V\) makes the graph of \((U,V)\) simple.
Regard the other parameters of the record as fixed. A
rational transfer with \(\ell\ge2\) active slots has
exactly two mixed irreducible factors in different terms
and a pure term proportional to
\(\prod_i W_i^{e_i}\), where \(W_i\) is the chosen original
free letter \(U\) or \(V\) in slot \(i\), and \(e_i\ge1\).

Sort the edges by the dyadic ranges of their original
endpoint degrees on both sides. Some class retains
\(E\ge N/O((1+\log N)^2)\) edges and has endpoint degree
caps \(2D_U,2D_V\), with all its endpoint original degrees
at least \(D_U,D_V\) on the respective sides. Consequently
the retained numbers of distinct endpoints are at most
\(N/D_U,N/D_V\). Simplicity gives
\[
 E\le N^2/(D_UD_V).
\]
Define \(a_W=\log_n(E/(2D_W))\) for \(W=U,V\). A negative
value can be replaced by zero with an inner \(o_X(1)\)
change, since the largest original degree is at most \(N\).
The retained exponents satisfy
\[
                 a_U+a_V\ge r-o_X(1).                 \tag{PT1}
\]
Each \(a_W\) is the retained dictionary exponent up to
\(o_X(1)\), by the two bounds \(E/(2D_W)\) and \(N/D_W\)
on its number of distinct values.

For any fixed polynomial-height auxiliary prime pool and
for nonexceptional anchors, a unary factor \(W\) in one
slot supplies normalized radical mass at least
\(a_W-\epsilon\) exclusive of the other output factors
outside that pool, on all but a power-smaller fraction
of tuples. Indeed fix the other rows, put their unary
values and the relevant resultants with the mixed
factors in the pool, and use the uniform fixed-support
enumeration. A bound \(n^{a_W-\epsilon}\) on the outside
radical allows at most \(n^{a_W-\epsilon+o_X(1)}\)
values. Each occurs on at most \(2D_W\) retained rows.
If the exponent cap is negative the lower-bound assertion
is automatic. The exceptional anchor fraction is
\(O(1/E)\), as in the resultant proof. Retain one prime
contribution per unary factor even if the factor occurs
to a higher power. On taking the finite union of all bad
events, the unary contributions are exclusive of one
another and of the two mixed-factor contributions.
Each survives the common gcd in the primitive output.
[Resultant pools and the mixed-factor bounds](ABC_CROSS_RATIO_FRESH_PRIME_OBSTRUCTION.md#fresh-radical-costs-of-extra-polynomial-factors)

If both free-letter orientations occur, one slot of
each type supplies combined exponent \(r-o_X(1)-2\epsilon\)
by (PT1). The two mixed factors supply
\(2(r-\epsilon)\). Thus the primitive output has
normalized total radical at least \(3r-o_X(1)-4\epsilon\)
on typical tuples of this logarithmically thinned fibre.

## A fixed pure term on uniform marginal slices

Suppose all slots choose the same free letter. Call it \(U\). Thin using only the dyadic
class of the original \(U\)-degrees, retaining \(E\ge
N/O(1+\log N)\) edges. Here the entire fibre of every
retained \(U\)-value is retained and has between \(D\) and
\(2D\) rows. Put \(a=\log_n(E/D)\) and \(q=\log_n D\), so
\(a+q=r+o_X(1)\). The same unary count and the two mixed
factors yield
\[
  \log_n\operatorname{rad}(\text{primitive output})
          \ge 2r+\ell a-o_X(1)-(\ell+2)\epsilon          \tag{PT2}
\]
on typical tuples.

Permute the three output terms by an anharmonic
transformation so that the pure term is the denominator.
The transfer in the coordinate \(x_i=V_i/U_i\) is then a
fixed rational-coefficient polynomial \(P(x_1,\ldots,x_\ell)\).
Indeed each output polynomial is multihomogeneous with
the same multidegree as the pure term, and division by
that product of powers of \(U_i\) yields a polynomial in
the \(x_i\).

Fix an absolute \(b>0\), bounded as \(L\) decreases, and a
positive numerical height lower bound \(h\) for the primitive
output. Call a tuple eligible if its output has all three
terms nonzero, largest absolute value at least \(X^h\), and
total radical at most \(n^b\). Its largest term is at most
\(X^{O(1)}\) because the transfer is fixed. In any collection
obtained by fixing the \(U\)-values for the moving query
slots and all rows in the other slots, the raw pure output
term \(H\) is fixed. After primitive division its absolute
value is a divisor of \(|H|\), so it has \(X^{o_X(1)}\)
possibilities. After reordering absolute values each output
is an actual positive primitive \(abc\) triple and its signed
projective value \(P\) has bounded multiplicity for each
such reordered triple.

Let \(c_*=1/3-\sigma_*\) be a fixed coefficient strictly
less than \(1/3\) from the uniform fixed-complete-letter
estimate (FT1). It follows that the number of eligible
distinct \(P\)-values in the specified slice is at most
\[
                           n^{c_* b+o_X(1)}.             \tag{PT3}
\]
The error means an arbitrarily small inner exponent
allowance. For example, bin outputs by their
\(\log_X\)-heights into intervals \([t,t+\zeta]\) of fixed
small width \(\zeta\) with \(t\ge h-O(\zeta)>0\). Use (FT1)
in that bin with budget just above \(Lb/t\) and height
cap \(X^{t+\zeta}\), then make \(\zeta\) and the exponent
allowance small. The output budgets fall within (FT1)
for sufficiently small fixed \(L\) in terms of \(b,h\).
The bounds are uniform in the fixed nonzero integer
\(H\) and its divisors.
[Fixed-complete-letter saving (FT1)](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-strict-one-third-saving-with-any-one-complete-letter-fixed)

For one moving slot the specialized polynomial is
nonconstant outside a power-smaller fraction of anchor
choices. This follows by choosing a nonzero polynomial
coefficient of its derivative and applying the Cartesian
polynomial zero bound to the distinct anchor directions.
On fixing the corresponding \(U\)-value, the number of
eligible queries is at most \(O_P(1)\) times (PT3), as a
nonconstant polynomial has bounded fibres. Averaging over
fixed rows and the uniform marginal fibres therefore gives
\[
                              c_* b\ge q-o(1)             \tag{PT4}
\]
whenever the eligible fraction has no power saving. In
this and the bounds below the \(o(1)\) represents
arbitrarily small exponent allowances. Precisely, a
fixed strict violation of such an inequality gives a
power-smaller tuple fraction by the same argument with a
small fixed allowance. Together (PT2) and (PT4) exclude a
budget \(b\le 3r\), with a uniform positive margin
proportional to \(r\), if \(\ell>1/(3\sigma_*)\).

## Bivariate and trivariate slices

Fix two or three query slots, then fix the other rows and
the \(U\)-value in each query slot. The query variables
range over real sets of sizes between \(D\) and \(2D\);
their element denominators divide the corresponding
fixed \(U\)-value. A lower bound \(n^{-\eta}\) on the
eligible tuple fraction, for any small fixed \(\eta>0\),
would leave a slice with that density up to a constant
by averaging. If a type of specialized polynomial is
specified, the same assertion applies to the fraction of
tuples that are both eligible and of the specified type.
A distinct-variable restriction in the original tuple
merely restricts the eligible part of each grid. In
particular the grid incidence bounds below do not require
it to be the whole grid.

For a bivariate specialization \(G(x,y)\), use the
Raz–Sharir–Solymosi bound
\[
 \#\{(x,y)\in A_1\times A_2 : G(x,y)\in C\}
       \ll (|A_1|^{2/3}|A_2|^{2/3}+|A_1|+|A_2|)\,|C|^{1/2}
\]
when \(G\) is neither \(g(\phi(x)+\psi(y))\)
nor \(g(\phi(x)\psi(y))\) for real univariate polynomials.
The degrees are fixed and the estimate is uniform in the
coefficients. Take \(C\) to consist of the eligible
values, apply (PT3), and compare with \(n^{-\eta}D^2\).
The necessary inequality as the allowance decreases is
\[
                             c_* b\ge \tfrac43 q-o(1).    \tag{PT5}
\]
A certificate that all but a power-smaller fraction of
anchor specializations have the required nonspeciality
is that \(\partial_x\partial_y\log(G_x/G_y)\) is nonzero
over the field of the anchored variables. Here the
logarithmic derivatives are rational identities. Both
special forms have \(G_x/G_y\) a quotient of separate
one-variable rational functions and hence give zero.
Exclude the zeros of one nonzero coefficient after
clearing denominators.
[Raz, Sharir and Solymosi, Theorem 2](https://arxiv.org/pdf/1401.7419#page=4)

For a trivariate specialization \(G(x,y,z)\) depending
on all three variables, consider the polynomial
\(t-G(x,y,z)\). It is irreducible with all four partial
derivatives nonzero as polynomials. The nonspecial
case of the four-dimensional incidence theorem gives
\[
 \#\{(x,y,z)\in A_1\times A_2\times A_3 : G(x,y,z)\in C\}
                       \ll D^2|C|^{2/3}+D^2+D|C|.
\]
By the proof of the expansion corollary in the same
paper, its exceptional case for a real polynomial graph
implies one of the two compositions by a sum or product
of three univariate polynomials. Comparing with
\(n^{-\eta}D^3\) and using (PT3) gives the necessary
inequality
\[
                             c_* b\ge \tfrac32 q-o(1).    \tag{PT6}
\]
When \(q\) tends to zero this inequality is automatic
up to the allowance. A nonzero derivative
\(\partial_z(G_x/G_y)\) over the field of the anchored
variables is a sufficient generic certificate, provided
all three partial derivatives of \(G\) are nonzero.
Both exceptional compositions give zero. For exactly
three active slots there are no anchored variables and
one can test the nonspeciality directly.
[Raz, Sharir and de Zeeuw, Theorem 1.1 and proof of Corollary 1.2](https://arxiv.org/pdf/1607.03600#page=2)

The multiplicative composition case also gives strong
expansion on these rational slices. It suffices to fix
two nonconstant variables in the composition and all
other rows generically. For a bivariate specialization
\(G=g(\phi(x)\psi(y))\) with nonconstant \(\phi,\psi,g\),
normalize \(\phi,\psi\) to be monic by rescaling \(g\).
If \(e=\deg g\), the leading coefficient of \(G\) in
\(y\), divided by its own leading coefficient in \(x\),
is \(\phi(x)^e\). The analogous formula holds for
\(\psi\). Since the coefficients of \(G\) are rational
and of polynomial height, comparison of coefficients
starting at the leading term shows that \(\phi,\psi\)
have rational coefficients of polynomial height. The
recursive divisions by \(e\) and the number of
operations are bounded in the fixed degree. Each
polynomial evaluated on its query set has a common
denominator of polynomial height, because the
corresponding \(U\)-value is fixed. After clearing
these denominators, a nonzero product has at most
\(X^{o_X(1)}\) factorizations into the two resulting
integer values. The univariate polynomials have
bounded fibres and their zero values account for only
\(O_P(D)\) query pairs. Thus a set of
\(n^{-\eta}D^2\) eligible query pairs gives at least
\(n^{-\eta}D^2/X^{o_X(1)}\) distinct \(G\)-values. Equation (PT3) now
gives
\[
                              c_* b\ge 2q-o(1).          \tag{PT7}
\]

Combining (PT2) with (PT5), (PT6), or (PT7) yields a
strict budget floor above \(3r\) with uniform positive
room in the following respective ranges:
\(\ell\ge4\) and a nonspecial two-slot slice,
\(\ell\ge3\) and a nonspecial three-slot slice, or
\(\ell\ge2\) and a multiplicatively separated two-slot
slice. For instance a putative \(b\le3r\) forces
\(a\le r/\ell\) by (PT2), and the right sides of
(PT5)–(PT7) are then at least \(r\), while
\(c_*b\le(1-3\sigma_*)r\). Sufficiently small fixed
height and density allowances preserve the strict
margin. These conclusions require the specified height
lower bound and radical cap for the transferred outputs.

## Effective radical costs inside the same slices

Use the same \(U\)-degree slices and transfer polynomial.
On a typical input tuple, the unary factors of the pure
term contribute at least \(\ell a-o(1)\) normalized
radical mass to the primitive output, on disjoint
supports. This is the unary part of (PT2).
Discarding its power-small exceptional fraction
does not affect an eligible fraction with no power saving.

Within a slice the raw pure term \(H\) is fixed.
Price its primes at zero and the others at their full
logarithmic cost. The partition sum is sublogarithmic
since \(H\) has polynomial height. The primitive pure
term is supported on these primes, and \(\ell a-o(1)\)
of the output radical budget is carried on them.
Binning the output heights as in (PT3), the
[effective radical bound](ABC_EFFECTIVE_RADICAL_BUDGET.md#counting-by-a-sublogarithmic-prime-partition-sum)
gives at most
\[
                  n^{(b-\ell a)/3+o(1)}                       \tag{PT8}
\]
eligible polynomial values in the slice. The
[zero-cost letter gap](ABC_EFFECTIVE_RADICAL_BUDGET.md#a-strict-effective-saving-when-one-letter-is-free-to-record)
gives a uniform positive subtraction from this
exponent if \(b-\ell a\) has a fixed positive lower
bound. Its budget comparison applies because the
given upper and lower numerical output heights
convert \(Lb\) to a sufficiently small fixed radical
budget comparable to \(L\).

Write \(\gamma=1,4/3,3/2,2\) for the lower exponents
on the output dictionary, in units of \(q=r-a\),
from one-slot bounded fibres, nonspecial bivariate
slices, nonspecial trivariate slices and multiplicative
bivariate slices respectively. The same arguments
that gave (PT4)–(PT7) now imply
\[
                    b\ge \ell a+3\gamma(r-a)-o(1).            \tag{PT9}
\]
At a fixed positive lower bound on \(r\), any candidate
eligible fraction without power saving must have
\(b-\ell a\ge2r-o(1)\) by (PT2), so the strict gap
in (PT8) supplies a positive margin in (PT9).

A nonspecial bivariate slice therefore excludes \(b\le3r\)
for every \(\ell\ge2\). Indeed (PT2) forces
\(a\le r/\ell+o(1)\) at such a budget, and hence
\(\ell a+4(r-a)\ge3r-o(1)\).
The one-slot argument alone excludes \(b\le3r\) for
every \(\ell\ge3\), since
\(\ell a+3(r-a)\ge3r\). The multiplicatively separated
case still excludes \(b\le3r\) for \(\ell=2\).
With the stated numerical output height bound, the
remaining case for this direct argument is a two-slot
additive composition with both pure factors on the
same original free letter.



## Mixed orientations and the two-slot additive case

The dyadic regularization in (PT1) also gives typical slice sizes
\(n^{r-a_U+o(1)}\) and \(n^{r-a_V+o(1)}\) for a row drawn
uniformly from the retained edges. For precision, at a fixed
\(\epsilon>0\) the number of rows at values of \(W\) whose
retained degree is at most \(n^{r-a_W-\epsilon}\) is at most
\[
 \frac{N}{D_W}n^{r-a_W-\epsilon}=E n^{-\epsilon+o_X(1)}.
\]
Use this to restrict each chosen pure-letter slice without
discarding an eligible fraction that lacks a power saving.
The upper slice bounds are the original dyadic degree caps.

Let \(s=\sum_i a_{W_i}\), allowing either orientation, and
suppose an eligible fraction lacks a power saving at \(b\le 3r\).
The mixed-factor bound and the unary argument give
\(b\ge 2r+s-o(1)\). On fixing the pure letters in all moving
slots the raw pure term is constant, so (PT8) holds with \(s\)
in place of \(\ell a\). At fixed \(r\) bounded below, the
zero-cost letter gap gives an upper exponent \(2r/3-\delta\)
for the number of eligible outputs in that slice, with a
fixed positive \(\delta\) after the allowances are taken small.
In fact \(s\ge a_U+a_V\ge r-o(1)\), and
\(b-s\ge2r-o(1)\) as required for that gap.

If both orientations occur and there are at least three active
slots, one of them repeats. Its marginal exponent tends to zero
because \(s\le r+o(1)\). Varying a slot on that side with its
pure letter fixed yields a query set of exponent \(r-o(1)\).
A generic specialization of \(P\) is nonconstant in that slot,
so the bounded polynomial fibre contradicts the output bound.
The Cartesian polynomial zero estimate used for (PT4) applies
to the full retained set of anchor directions before any
pure-letter slice is chosen.

For exactly two mixed slots write the two typical query
exponents as \(q_1=r-a_U\) and \(q_2=r-a_V\). Their sum is
\(r+o(1)\). The one-slot bounded fibres show
\(\max(q_1,q_2)\le2r/3-\delta+o(1)\) at an eligible density
without a power saving. Choose a pair of fixed-pure-letter
slices with that density by averaging, restricting to the
typical slice sizes above. For a nonspecial bivariate map,
the unequal-size Raz–Sharir–Solymosi bound gives incidence
exponent at most
\[
 \max\{2(q_1+q_2)/3,q_1,q_2\}+\frac12(2r/3-\delta)
             \le r-\delta/2+o(1).
\]
This is below the number of eligible pairs. A multiplicatively
separated bivariate map has the full product dictionary exponent
\(q_1+q_2-o(1)\) by the same rational normalization and divisor
argument as before, allowing unequal sizes.

For an additive composition in the two mixed slots, use the
coordinates
\[
                x=V_1/U_1,\qquad y=U_2/V_2 .
\]
Write \(P(x,y)=g(\phi(x)+\psi(y))\) with three nonconstant
real univariate polynomials. This representation may be taken
over \(\mathbb Q\). Translate and scale the inside arguments
to make \(\phi(0)=\psi(0)=0\) and \(\phi\) monic. If \(d=\deg\phi\)
and \(j=\deg g\), the coefficient of \(x^{dj}\) in \(P(x,0)\)
is the leading coefficient \(g_j\) of \(g\). The next \(d-1\)
descending coefficients of \(P(x,0)\) recover the nonconstant
lower coefficients of \(\phi\) recursively: the coefficient
of \(x^{dj-k}\) for \(1\le k<d\) depends linearly on the
coefficient of \(x^{d-k}\) in \(\phi\) with nonzero coefficient
\(j g_j\), with all other contributions already known.
Then division in powers of the known monic \(\phi\) recovers
\(g\) from \(P(x,0)\). The coefficient of \(x^{d(j-1)}\) in
\(P(x,y)-P(x,0)\) is \(j g_j\psi(y)\), which recovers \(\psi\).
For \(d=1\) the intermediate coefficient step is empty.

Fix \(U_1=Fu_0=:A\) and \(V_2=Gv_0=:B\). Suppose \(V\) is the
shorter free letter; otherwise exchange the letters and the
slots. In the first slice put \(V_1=Gv\), and in the second put
\(U_2=Fu\). The affordable margin gives \(v\le o(F)\) uniformly.
For a fixed output of \(P\) there are at most \(j\) possible
values of \(\phi(Gv/A)+\psi(Fu/B)\). For each such value that
occurs, compare to one witness \((v',u')\). Take one fixed
positive integer \(J\) clearing the denominators of both
univariate polynomials and set
\[
                   Q(Z)=JA^d\phi(GZ/A)\in\mathbb Z[Z].
\]
For \(p\mid F\), the factors \(G\) and \(B\) are units at \(p\):
the latter divides the actual second row's letter opposite
\(F\). Likewise \(v,v'\) are units since \(Gv,Gv'\) are
opposite \(F\) in their own rows. The difference of the \(J\psi\)
values is \(p\)-integral, and equality of the sums implies
\[
                 Q(v)-Q(v')\equiv0\pmod{p^{d v_p(F)}}.
\]
Write \(f=v_p(F)\) and \(C_p=v_p(dJ\phi_d)\), with \(\phi_d=1\)
in the monic normalization above. For \(f>2C_p\), every unit
argument has \(Q'\)-valuation exactly \(C_p\): the leading
derivative term is \(dJ\phi_dG^dZ^{d-1}\), while all lower
terms have a factor \(A\), of valuation at least \(f\).
The strong Hensel estimate places a solution of the congruence
in a residue class of an exact root of \(Q(Z)-Q(v')\) modulo
\(p^{df-C_p}\). Thus there are at most \(d\) allowed classes
modulo \(p^{(f-2C_p)_+}\); the statement is trivial at the
discarded powers \(f\le2C_p\). For all primes outside one fixed
finite set \(C_p=0\). Combining the congruences leaves a modulus
at least \(F/\prod_{C_p>0}p^{2C_p}\) with at most
\(d^{\omega(F)}=X^{o_X(1)}\) residue classes. Its polynomial
margin over the range of \(v\) gives at most one \(v\) per
class. Once \(v\) is chosen the fixed value of \(\psi(Fu/B)\)
allows boundedly many \(u\). The possible output values
therefore each have only \(X^{o_X(1)}\) pairs in that slice.

This gives the full product exponent \(q_1+q_2-o(1)\) also
for additively separated two-slot maps with opposite pure
letters. With the stated numerical output height guard, the
only surviving pattern in the direct \(b\le3r\) argument is
the same-orientation two-slot additive composition.

## Denominator recovery on the eligible graph

Support refinement on the eligible graph itself
recovers the original free-letter tags from the reduced
output denominator. This supplies a bound for every
pure-term polynomial chart at outputs exceeding a fixed
power of \(n\). At smaller heights a local root bound
also excludes the two-slot same-orientation case and
the opposite-orientation additive case. These results
hold on arbitrary eligible subgraphs of sufficiently
large exponent, rather than requiring a product of
full slices.
[Denominator recovery and local bounds for small outputs](ABC_PURE_TERM_DENOMINATOR_RECOVERY.md)

## An additive shape without expansion

The affine polynomial \(P(x,y)=1+x+y\) already has
the two-mixed-factor pattern with pure term
\(U_1U_2\). Its other two homogenized terms are
\(U_1U_2+V_1U_2+U_1V_2\) and
\(-V_1U_2-U_1V_2\). Both are irreducible bilinear forms,
since the two coefficient matrices have nonzero
determinants.

An arithmetic example shows why expansion in all
fixed-letter slices cannot be inferred from the
affordable divisor congruences and primitivity. For
an integer \(m\), set \(X=6^m\), \(U=F=5^m\),
\(G=2^m\), and choose
\[
          V=Gv,\qquad T=U+V,\qquad
          v=5j+1,\quad 0\le j<\big\lfloor(3/2)^m/5\big\rfloor .
\]
These are distinct primitive rows with \(V<T\le X\).
The modulus \(FG=10^m\) exceeds the shorter free height
by a fixed numerical margin. The directions \(V/U\)
form an arithmetic progression, and \(P(x,y)\) on the
full square of this progression has \(2N-1\) values
when it has size \(N\). The original targets and
transferred outputs have unrestricted radicals in this
example. The additional arithmetic restriction needed
for the additive compositions is precisely the
simultaneous control of those actual radicals.
