# High free degrees for the sextic tangent-conic scalar

**Theorem SexticTangentConicHighMateRestriction.** Suppose a
characteristic-zero three-variable Keller map of ordinary degree at
most six has a quartic component \(R=g(x)+x^2(xz+y^2)\).
Subtract \(g(0)\) and write \(r=g(x)+x^3v\) after the localized
substitution \(z=v-y^2/x\). After a constant linear combination
of the other two components, their ordered free \(y\)-degrees
and leading coefficients, up to nonzero scalar multiples, must
be one of
\[
\begin{array}{c|c}
(\deg_y P,\deg_y Q)&([y^{\deg_y P}]P,[y^{\deg_y Q}]Q)\\ \hline
(9,6)&(x^{-3},x^{-2})\\
(12,4)&(x^{-6},x^{-2})\\
(12,6)&(x^{-6},x^{-3}).
\end{array}
\]

Work over \(\mathbb C\) by coefficient-field embedding. Write
\(g(0)=0\); submersion makes \(g'(0)\ne0\).
The substitution has determinant one on \(x\ne0\), and the
coefficients belong to
\(\mathbb C[r,x,x^{-1}]=\mathbb C[x,x^{-1},v]\).
The derivation in the bracket is \(D=-x^3\partial_x|_r\), so
\([P,Q]=D(P)Q_y-P_yD(Q)=\kappa\ne0\). The actual values and first
\(y\)-derivatives at \(y=0\) are regular at \(x=0\) in the variables
\((x,v)\). Also \(P(0,y,z),Q(0,y,z)\) in the original affine
coordinates are a plane automorphism, by the Keller equation
and [the retained planar bound](LOW_DEGREE_COORDINATE_KELLER_MAPS_IN_THREE_VARIABLES.md).
The [low-mate theorem](SEXTIC_KELLER_TANGENT_CONIC_QUADRATIC_MATE.md)
already gives \(m=\deg_y P,n=\deg_y Q\ge3\).

## Leading support and the degree pairs

Use the group bounds in the
[localized support calculation](SEXTIC_KELLER_TANGENT_CONIC_QUADRATIC_MATE.md#support-and-the-affine-cases).
The free degrees are at most twelve. In a leading coefficient
\(p\) of free degree \(m\), the \(x\)-exponents at \(v^j\) lie in
\[
              j-\lfloor m/2\rfloor\ \le i\le 6-m-j,
       \qquad j\le\left\lfloor(6-\lceil m/2\rceil)/2\right\rfloor .
 \tag{1}
\]
For any coefficient of \(y^s v^j\), including below the leading
free degree, the same grouping gives the upper bound
\(i\le6-s-j\). At fixed \(r\), (1) gives \(r\)-degree at most
the displayed bound on \(j\), highest \(x\)-exponent at most
\(6-m\) since \(\deg g\le4\), and lowest exponent at least
\(-\lfloor m/2\rfloor-2\lfloor(6-\lceil m/2\rceil)/2\rfloor\).

Write \(k=\mathbb C(r)\). The leading bracket identity is
\(p^n/q^m\in k\), where \(q\) is the other leading coefficient.
If \(m=n\), then \(p/q\in k\). A nonconstant ratio forces
one factor of the form \(r-c\) in one coefficient. For
\(m\ge5\) the \(r\)-degree is at most one, so constant target
combinations would give leading coefficients \(\ell(x),r\ell(x)\).
The first requires lowest exponent at least \(-\lfloor m/2\rfloor\)
and the second highest exponent at most \(2-m\), which are
incompatible. For \(m=3,4\), the \(r\)-degree is at most two.
If the common factor has \(r\)-degree one, constant target
combinations give leading coefficients \(L,rL\), with
\(L=A(x)r+B(x)\). The \(v^2\)-coefficient of \(rL\) in (1)
forces \(A(x)\) to be a scalar times \(x^{-m-2}\). Now the
\(v\)-coefficient \(x^3 A(x)\) of \(L\) violates the lower bound.
Thus the common factor \(\ell\) is univariate. After constant
combinations the maximum invariant degree cannot be two,
since its upper bound on \(\ell\) would be less than the lower
bound from the other coefficient. The leads are therefore
\(\ell,r\ell\), with \(\ell\) proportional to \(x^{2-m}\).
Subtract \(r\) times the first mate. If the second is independent
of \(y\), its original polynomial restricts on \(x=0\) to a
polynomial in \(y^2\), impossible for a plane coordinate.
Indeed its Laurent \(v^j\)-coefficient must be divisible by
\(x^j\) by evaluating the original polynomial at \(z=0\).
If the new free degree is \(n'>0\), the leading ratio gives
\(m\mid (m-2)n'\). For \(m=3\) this is impossible, and for
\(m=4\) it gives \(n'=2\) with odd \(x\)-valuation of the
quadratic coefficient, contrary to the
[quadratic-mate square relation](KELLER_REPEATED_PARALLEL_SCALAR_QUADRATIC_MATE.md).
Thus equal degrees have constant leading ratio and constant
target cancellation leaves unequal degrees within the original
bound.

Now \(m>n\ge3\). Let \(d=\gcd(m,n)\).
Factorization in \(\mathbb C[r,x,x^{-1}]\) and the ratio
give
\[
                 p=e(r)\ell(x)^{m/d},\qquad
                 q=f(r)\ell(x)^{n/d}.
\tag{2}
\]
Indeed a common irreducible factor involving both \(r,x\) would
have exponent at least \(m/d\) in \(p\). For \(m/d\ge3\) this
exceeds the maximal \(r\)-degree two. For \(m/d=2\) one has
\(m\ge6\) and maximal \(r\)-degree one. Thus only invariant
factors and univariate factors occur. The \(x\)-orders of
\(p,q\) have ratio \(m/n\) so their Laurent powers of \(x\)
also take the form (2), with \(e,f\in\mathbb C[r]\).

If \(\ell\) is not a Laurent monomial, there are only the
possibilities
\[
      (m,n)=(6,3)\ \text{or}\ (6,4),\qquad
     \ell\ \text{proportional to }1-a/x,\quad a\ne0 ,
     \qquad e,f\ \text{constant}.
\tag{3}
\]
To see this, use the
[formal residue identity](https://arxiv.org/pdf/2610.36495v1)
of Wang and Yan, Proposition 3.17. In the tangent-conic
coordinates it requires a primitive of \(dx/(x^3\ell^{1/d})\)
over the field adjoining the indicated root and an invariant
root of \(e\). Multiplying that primitive by \(\ell^{1/d}\)
and taking normalized trace gives a rational solution of
\[
                       U'-(\ell'/d\ell)U=x^{-3}.
\tag{4}
\]
A pole of \(U\) at nonzero \(x\) is impossible because the
leading derivative coefficient is a nonzero multiple of
\(-j-b/d\), with pole order \(j>0\) and root order \(b\ge0\).
The lower bound of (1) at fixed \(r\) gives
\({\rm ord}_0(\ell)/d>-2\), so (4) forces \(U\) to have
exactly a double pole at zero.
For \(m\ge5\) the ratio of highest \(x\)-exponent of \(p\)
to \(m\) is strictly between \(-1\) and \(1\).
At infinity (4) therefore permits a highest exponent
above \(-2\) in \(U\) only if that exponent and the ratio
are both zero. Thus either \(U=A/x^2\), forcing a
monomial \(\ell\), or the highest exponent of \(\ell\)
is zero and \(m\le6\). Now the upper bounds in (1) applied
to the highest powers of \(v\) in (2) force \(e,f\) constant.
The lower bounds at \(v^0\) then exclude a nonmonomial
for \(m=5\) and leave exactly (3) for \(m=6\).
For \(m=4,n=3\), the extremal exponents of \(p\) are multiples
of four and lie in \([-6,2]\). Nonmonomial support would
require the pair \((-4,0)\). Highest \(v\)-degree then makes
\(e\) constant, but the \(v^0\) bound forbids the \(-4\) exponent.

For monomial \(\ell\), write \(p=x^a e(r),q=x^b f(r)\).
The \(v^0\) and highest \(v\)-coefficient bounds in (1) give
\[
 \begin{array}{ll}
 a+\operatorname{ord}_0e\ge-\lfloor m/2\rfloor,&
 a+4\deg e\le6-m,\\
 b+\operatorname{ord}_0f\ge-\lfloor n/2\rfloor,&
 b+4\deg f\le6-n,\qquad na=mb .
 \end{array}
\tag{5}
\]
Here \(g'(0)\ne0\) makes the order at zero of \(e(g(x))\)
equal \(\operatorname{ord}_0e\).
The finite comparisons in the
[exact checker](verify_sextic_tangent_conic_high_mates.py)
give the following choices. An entry \((a,b)\) has
constant invariant factors unless otherwise indicated:
\[
\begin{array}{c|l}
(m,n)&(a,b)\\ \hline
(4,3),(5,3),(5,4),(6,5)&(0,0)\\
(6,3)&(0,0),\ (-4,-2)\ (e,f\text{ linear}),\
                       (-2,-1)\ (f\text{ constant or linear})\\
(6,4)&(0,0),\ (-3,-2)\ (f\text{ constant or linear})\\
(8,4)&(-4,-2)\ (f\text{ constant or linear}),\ (-2,-1)\\
(8,6)&(-4,-3)\\
(9,3)&(-3,-1)\ (f\text{ constant or linear})\\
(9,6)&(-3,-2)\\
(10,4)&(-5,-2)\ (f\text{ constant or linear})\\
(10,5),(10,6),(10,8)&(-4,-2),\ (-5,-3),\ (-5,-4)\ \text{respectively}\\
(12,4)&(-6,-2)\ (f\text{ constant or linear})\\
(12,6),(12,8),(12,10)&(-6,-3),\ (-6,-4),\ (-6,-5)\ \text{respectively}.
\end{array}
\]
The table records the degrees of invariant factors, with the
additional vanishing conditions given by (5).

## Eliminations at low orders

First eliminate the entries with \(\gcd(n,b)=1\) for a monomial
leading coefficient \(q=x^bf(r)\). In any subsequent mate the
leading bracket requires its free degree to be divisible by
\(n\), by the \(x\)-valuation at fixed \(r\).
For such a degree, subtract the appropriate power of \(Q\)
with its invariant leading multiplier over \(k\).
Repeated subtraction gives a mate independent of \(y\),
whose bracket with \(Q\) is zero or has degree \(n-1>0\).
This excludes the relevant entries in \((6,3),(8,4),(9,3),(10,5)\).
It also excludes the first nonmonomial case \((6,3)\) in
(3), by using the simple root of \(q\) at \(x=a\).

For coprime \((m,n)\) with constant leads, scale both mates
to monic and center the smaller by translating \(y\).
After dividing the bracket by \(x\), the derivation is
\(-x^2\partial_x|_r\) and the right side is \(\kappa/x\).
The [coprime pole lemma](SEXTIC_KELLER_TRIPLE_LINE_COPRIME_POLE_LEMMA.md)
forces \(Q\) in the centered variable \(Z\) to have
coefficients polynomial in \(r\), except its nonzero scalar
term \(\eta/x\). Indeed the invariant Laurent coefficients
lie in \(\mathbb C[r]\), and the lemma gives
\(\eta^2=-n\kappa/(n+1)\). At the original \(y=0\),
if \(Z\) is regular at \(x=0\) this simple pole persists;
if it has a pole, \(Z^n\) dominates. Each case contradicts
regularity of the actual value \(Q|_{y=0}\).

For \((6,3)\) with constant leads, subtract a constant
multiple of \(Q^2\) to lower the larger positive degree
to at most five. This is a polynomial operation in the
original frame. If the remaining degree is five or four,
its leading coefficient \(A(r)\) is invariant by (1)
and lies in \(\mathbb C[r]\). Center \(Q\), divide the
bracket by \(x\), and divide the larger mate by \(A(r)\).
The [quintic-cubic calculation](SEXTIC_KELLER_TRIPLE_LINE_CUBIC_MATE.md#two-centered-cubic-calculations)
rules out degree five over \(k[x,x^{-1}]\) with right
side \(\kappa/(A x)\). For degree four the same calculation
forces the centered cubic of (4) in that proof.
Its coefficient \(\eta(r)\) is polynomial in \(r\) because
the actual coefficients are Laurent over \(\mathbb C[r]\),
and \(A\eta^2\) is a nonzero scalar. Hence \(\eta\)
is scalar and again yields the preceding contradiction at
the original \(y=0\). A remaining degree three cancels
by an invariant polynomial multiple of \(Q\), since the
multiplier is an invariant in \(\mathbb C[r,x,x^{-1}]\).

If the remainder is quadratic, call it \(Q_2\) and the
constant-leading cubic \(P_3\). The leading coefficient
of \(Q_2\) is \(q_2=e(r)\in\mathbb C[r]\) by the leading
ratio. Use the
[quadratic-mate integration](KELLER_REPEATED_PARALLEL_SCALAR_QUADRATIC_MATE.md#parity-after-centering)
over \(k[x,x^{-1}]\) with \(w=1\):
\(Z=y+B/(2e)\), \(Q_2=eZ^2+K\), and
\(P_3=\Psi(Q_2,r)+\int_0^Z F(Q_2-et^2)\,dt\).
Since \(D=-x^3\partial_x\), one has \(J(K)=\lambda/x^2+U(r)\)
with \(J'=F\). The odd part has degree three, hence \(F\)
degree one and Laurent factorization gives
\[
           K=\tau+u/x,\quad F(T)=2c(T-\tau),\quad
                        cu^2=\lambda\ne0 .
\]
The leading coefficient of the cubic is a nonzero constant,
and equals \(4ce/3\). Thus \(e\) is proportional to \(u^2\),
with \(u\in\mathbb C[r]\).
If \(e(0)\ne0\), \(\tau,u\) are regular at \(r=0\) by the
coefficient localization, and \(u\) is a unit there.
Then \(K\) has a pole at \(x=0\), inconsistent with
\(K=C-B^2/(4e)\) because \(C,B\) are regular there as the
actual value and first derivative at \(y=0\).
Therefore \(u\) has order \(b_0\ge1\) at \(r=0\).
On \(x\ne0,\ r=0\), the cubic retains its nonzero leading
coefficient. Hence \(B\) does not vanish identically there:
otherwise \(Q_2\) would be independent of \(y\) on this
component, contradicting the nonzero bracket.
It follows that \(K\), and thus \(\tau\), have pole
order \(2b_0\) there in the parameter \(r\).

For symbolic \(Q_2\), put \(Z_*^2=(Q_2-\tau)/e\).
The odd integral at \(Z_*\) is \(O_*=(4/3)ce Z_*^3\)
of pole order \(6b_0\), whose Laurent coefficients in \(r\)
are polynomials in \(Q_2\). For the matching branch
of \(Z\) the change is
\[
 O-O_*=-\frac{c(u/x)^2}{2e Z_*}
              \big(1+O((u/x)/(Q_2-\tau))\big).
 \tag{6}
\]
This results from integrating \(-F(K)/(2eZ)\) with respect
to \(K\) at fixed \(Q_2\). The ratio has positive order both
on \(x\ne0,\ r=0\) and at \(x=0\) after setting \(y=0\),
since \(u/x\) is regular at both. Also \(eZ_*\) is a unit
at both. Thus the change is regular on \(x\ne0,\ r=0\).
There \(Q_2\) restricts nonconstantly because \(B\ne0\), so
regularity forces all negative Laurent coefficients of
\(\Psi+O_*\) to vanish identically in symbolic \(Q_2\).
At \(x=0,y=0\), the opposite sign of \(Z_*\) leaves a pole
of order \(6b_0>2\); the matching sign leaves the nonzero
pole of order two from (6). Both contradict polynomiality.

If the remainder is affine instead, write it as \(Q_1=e(r)y+B\),
with invariant leading coefficient by the leading ratio.
The bracket on \(x\ne0,\ r=0\) forces \(e(0)\ne0\); otherwise
\(Q_1\) would lose all \(y\)-dependence there.
In \(Q_1\) coordinates the cubic is
\(\Psi(Q_1,r)+\kappa'/(2x^2e(r))\), with \(\kappa'\ne0\).
The first term has no poles at \(r=0\) by regularity on the
nonvertical component where \(Q_1\) is nonconstant.
The second term then contradicts regularity at the original
\(y=0,x=0\). An independent mate is impossible by the
even-coordinate argument in the equal-degree case.
This excludes every cubic-mate entry.

## Quartic mates and the plane at \(x=0\)

In (6,4) with constant leads, divide \(D\) by \(x\)
after centering and rescaling the two mates. This gives
the square case of the
[sextic-quartic centered equations](SEXTIC_KELLER_TRIPLE_LINE_SEXTIC_QUARTIC_MATES.md#the-square-case)
with \(w=1\), namely the monic Laurent bracket with
right side a nonzero scalar times \(1/x\).
That square-case contradiction uses only the bracket
and applies with this right side.

In the nonmonomial (6,4) case (3), use
\(w^2=1-a/x\) to center and rescale \(Q\) and rescale
the leading coefficient of \(P\). The centered
sextic-quartic bracket has right side \(\kappa/w\)
and derivation \(D=-x^3\partial_x\).
The [nonsquare centered equations](SEXTIC_KELLER_TRIPLE_LINE_SEXTIC_QUARTIC_MATES.md#the-nonsquare-case)
through the identity \(D(d^3/4)=\kappa/w\) use only
the involution \(w\mapsto-w\) and the constant field
of the extension, so apply here.
But \(dx/(x^3 w)=2(1-w^2)\,dw/a^2\).
Its primitive cannot be a cube of a rational function:
the polynomial \(w-w^3/3\) plus any additive constant
would require a linear cube with no quadratic
term and a nonzero linear term.

For the monomial (6,4) entry with \(q=H(r)x^{-2}\)
for a nonconstant linear \(H\),
the same centered argument uses \(w=x^{-1/2}\),
dividing the original quartic mate by \(H\).
It yields \(D(d^3/4)=\kappa/(Hw)\).
Since \(D(w)=1/(2w^3)\), one obtains \(d=\delta w\)
and \(\delta^3\) a nonzero constant times \(1/H\).
This is impossible at the simple root of \(H\).

For use with the other entries, restrict the original
plane automorphism at \(x=0\) to a generic line \(y=y_0\).
It is a polynomial embedding of the line, so one
of its two positive coordinate degrees in \(z\) must
divide the other.
[Abhyankar–Moh degree divisibility as proved by García Barroso–Płoski, Theorem 6.9](sources/garcia-barroso-ploski-1208.0913v1.pdf).
If a mate of free degree \(m\) has leading coefficient
proportional to \(x^{-\lfloor m/2\rfloor}\), the original
\(z\)-degree on that line is exactly \(\lfloor m/2\rfloor\).
Indeed the support grouping bounds the degree from
above. In the group \(N=-\lfloor m/2\rfloor,d=m\)
there is only one possible original monomial contributing
to the given Laurent coefficient at \(v^0\), namely
\(y^{m-2\lfloor m/2\rfloor}z^{\lfloor m/2\rfloor}\)
with no \(x\)-factor. Its coefficient is nonzero.
Thus the degree divisibility excludes the scalar-\(f\)
entries of (6,4), (10,4), and the entries (8,6),
(10,6), (10,8), (12,8), (12,10).

For any other \(n=4\) entry with \(q=(r-c)x^{-2}\)
up to a scalar, restrict to \(r=c,x\ne0\), where the constant bracket
still holds over \(\mathbb C(x)\).
The larger mate retains leading \(x^{-m/2}\), but the
quartic mate drops to a degree \(1\le n_0<4\).
It cannot become independent of \(y\), since the
other mate retains degree \(m>1\) there.
The leading ratio over \(\mathbb C(x)\) forces
\(x\)-valuation \(-n_0/2\) for the new leading
coefficient. Thus \(n_0=2\) with valuation \(-1\),
contradicting the quadratic-mate square relation.
This eliminates the nonconstant-\(f\) entries of
(8,4), (10,4), (12,4).

The entry \( (10,5) \) has valuation \(-2\) in its
smaller coefficient, coprime to five, so repeated
invariant cancellation eliminates it as above.
The same cancellation excludes (8,4) with \(b=-1\).
It remains to eliminate (8,4) with leads \(x^{-4},x^{-2}\).
Rescale these coefficients to equality of leading powers
and subtract \(Q^2\) from \(P\). The new degree is even
and at most six by its leading ratio. Degree four
cancels by an invariant polynomial multiple of \(Q\),
and degree two is incompatible with the odd
leading valuation \(-1\). Degree zero is impossible
by the bracket. In degree six the new coefficient is
\(F(r)x^{-3}\), where \(F\in\mathbb C[r]\) and
\(\deg F\le2\): at the highest \(v\)-degree \(j\), its
\(x\)-exponent is \(3j-3\), whereas the upper
bound for \([y^6](P-Q^2)\) is \(6-j\).

Use the [centered sextic-quartic equations](SEXTIC_KELLER_TRIPLE_LINE_SEXTIC_QUARTIC_MATES.md#the-nonsquare-case)
with \(w=x^{-1/2}\), after dividing the sextic mate
by \(F\). The Galois involution again gives
\(D(d^3/4)=\kappa/(F w)\). Hence \(d=\delta w\) and
\(\delta^3\) is a nonzero constant times \(1/F\).
The degree bound on \(F\) forces it constant, so
\(\delta\) is a nonzero scalar. With the centered
quadratic expression denoted by \(S\), the reconstructed
quartic and sextic shapes are
\[
 Q=S^2+dZ+\text{invariant},\qquad
 P'=S^3+\tfrac32 dZS+\tfrac38 d^2+
                        \text{invariant multiple of }Q+\text{invariant},
\]
where \(S=x^{-1}Y^2+\alpha(r)x\), \(dZ=\delta Y/x\),
\(Y=y+B_3/(4q)\), and the invariants are in \(\mathbb C[r]\)
by Laurent polynomiality.
At \(y=0\), if \(Y\) has a pole or is a unit at \(x=0\),
\(S^2\) gives the unique highest pole of \(Q\).
So \(Y\) vanishes there, and now \(3d^2/8=3\delta^2/(8x)\)
gives an uncancellable pole of \(P'\).
This contradicts the regularity of the actual
values at \(y=0\) and finishes the claimed restriction.

## Exclusion of the tangent-conic quartic form

**Corollary SexticTangentConicQuarticScalarExclusion.**
A characteristic-zero polynomial Keller map of ordinary degree
at most six cannot have a quartic component
\(R=g(x)+x^2(xz+y^2)\).

The leading restriction leaves only the localized free-degree
pairs \((9,6),(12,4),(12,6)\). Their separate exclusions are
proved in the [(9,6)](SEXTIC_KELLER_TANGENT_CONIC_NINE_SIX.md),
[(12,4)](SEXTIC_KELLER_TANGENT_CONIC_TWELVE_FOUR.md), and
[(12,6)](SEXTIC_KELLER_TANGENT_CONIC_TWELVE_SIX.md) calculations.
