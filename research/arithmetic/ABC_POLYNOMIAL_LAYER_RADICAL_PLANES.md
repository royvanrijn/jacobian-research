# A radical plane with several valuation layers on one carrier

A single parameter can vary several exact valuation layers
on one original letter. When its radical, a second radical,
and the third radical lie on a prescribed plane, one
quotient is enough to recover the row at height-subpower
cost. The quotient uses the squarefree part of a
polynomial in the parameter. The proof leaves the third
letter's valuation profile unrestricted.

## The quotient and its completion bound

Fix positive \(B_0\) and a bound for the degrees and
number of the carrier polynomials. Let
\(P_i(t)\in\mathbb Z[t]\) be nonconstant polynomials of
coefficient height at most \(X^{B_0}\), and write
\[
             P(t)=\prod_i P_i(t),\qquad m=\deg P\ge1 .
\]
Consider actual rows in which the first residual radical
is \(x=P(t)\), the second residual radical is \(y\), and
the third is \(z\), with \(t\) integral. Residual means
outside a prescribed prime pool of polynomial prime
product. In particular \(x,y,z\) are positive, pairwise
coprime and squarefree. The signed complete first and
second letters are
\[
       A(t)=a\prod_i P_i(t)^{p_i}=P(t)^q M(t),
          \qquad B y^q,\qquad
       M(t)=a\prod_i P_i(t)^{p_i-q},
\]
where \(a,B\) are fixed nonzero integers of magnitude
at most \(X^{B_0}\), and \(p_i\ge q\ge2\). The individual \(P_i(t)\) on an actual row are positive
exact valuation-layer roots. Assume that at least one
\(p_i>q\), so \(r=\deg M\ge1\). All polynomial factors
of \(M\) are factors of \(P\). The remaining signed
complete value is \(W=A(t)+B y^q\ne0\), with
\(\operatorname{rad}_{\overline{\mathcal P}}(W)=z\).
Assume \(|t|,|A(t)|,|B y^q|,|W|,y,z\le X^{B_0}\).

Prescribe nonzero rational \(\mu,\nu\), of polynomial
height, such that the actual radical plane gives
\[
                              y=\mu P(t)+\nu z .
\]
Put \(L(t)=M(t)+B\mu^q\) and clear denominators.
Take \(L_0\in\mathbb Z[t]\) to be a primitive product
of one copy of each irreducible factor of \(L\) not
dividing \(P\). Either sign can be prescribed
consistently. Write \(r_0=\deg L_0\le r\), including
the possibility \(L_0=1\). Enlarge the fixed pool by
the primes in the denominators and the coefficient
content used to make \(L\) primitive. The primes
of that content lie on the leading coefficient of
the cleared polynomial, hence on \(a\), the leading
coefficients of the \(P_i\), and the denominator of
\(\mu\). Thus the enlarged prime product has polynomial
height even before bounding the orders. Primitive
integer factors with their scalar content give exact
factorization identities in the localized coefficient
ring. If any enlarged-pool primes occur in \(z\),
record their squarefree product in \(z\) at subpower
cost and absorb it as a fixed multiplier in \(\nu\). From now on \(z\) is the third radical off the
enlarged pool. The coefficients of \(L\) and \(L_0\)
do not change in this operation.

Since \(y\equiv\mu P(t)\) at each prime of \(z\) and
\(P(t)\) is a unit there, the original equation makes
\(L(t)\) vanish at every such prime. A factor of \(L\)
that divides \(P\) cannot account for a prime of
\(z\). The removal of repetitions is permitted because
\(z\) is squarefree. Thus there is an exact integer
quotient
\[
                            k=L_0(t)/z .                    \tag{LC1}
\]
The factorization records and the stripping step are
uniform at polynomial height, by the fixed-support
bound and the fixed-degree coefficient bound for
integral factors.
[Uniform supported-value enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)
· [Unequal pure valuation orders](ABC_UNEQUAL_POWER_RADICAL_PLANES.md#the-exact-quotient-from-the-third-support)

**LC2.** For every prescribed \(k\), the number of
actual rows is \(O_\varepsilon(X^\varepsilon)\) in
either of the cases
\[
                        q\ge3,\qquad
                 q=2\ \text{and}\ r_0\le m .
\]
More generally the exact degree test for \(k\ne0\)
suffices whenever
\[
                  \max(qm+r,qe)>m+e+r_0-1,                 \tag{LC2}
\]
with \(e=\deg(\mu P+\nu L_0/k)\). The implicit constants
are uniform at fixed degree, order and polynomial height
bounds. No power condition is imposed on the third
letter. The coefficients of the prescribed carrier
give its valuation layers without fixing their actual
values.

For \(k\ne0\), \(z=L_0(t)/k\) and the second radical
is the single polynomial
\(u(t)=\mu P(t)+\nu L_0(t)/k\).
The polynomials \(A\) and \(u\) are coprime, since
all roots of \(A\) are roots of \(P\) and \(L_0\) is
coprime to \(P\). Thus the two summands in
\[
                            F=A+B u^q
\]
are coprime. Mason–Stothers gives at least
\(\max(qm+r,qe)-m-e+1\) distinct roots of \(F\),
counting finite complex roots. Under (LC2) at least
one of these lies off \(L_0=0\). If \(r_0\le m\),
then \(e\le m\) and \(qm+r\) makes (LC2) strict for
every \(q\ge2\). If \(r_0>m\), then \(e=r_0\) and
\(qe\) makes it strict for \(q\ge3\). This argument
also covers cancellation of the leading terms of
the two summands.
[Browning–Verzobio, Lemma 2.7](https://arxiv.org/html/2608.24512v1#S2.SS1)

Clear the rational denominators in \(F\) and remove
all primitive irreducible factors shared with \(L_0\).
The remaining factor product \(G\) is nonconstant,
coprime to \(L_0\), and has polynomial-height
coefficients. The quotient \(k\) itself has polynomial
height on any nonempty record at fixed degrees.
At an actual \(t\), a prime dividing \(G(t)\) outside
the fixed pool and the cleared denominator must
divide \(W\), hence \(z\), hence \(L_0(t)\). It
therefore divides the nonzero resultant of \(G,L_0\).
Add these resultant primes and the cleared denominator
to the fixed pool. The nonzero value \(G(t)\) is
supported on it. There are only subpower many such
values of polynomial height, each with at most
\(\deg G\) preimages. The parameter and \(k\) then
recover the original row. For \(L_0=1\) the same
argument applies with no common factors and no
resultant exclusions.

If \(k=0\), the parameter is a root of \(L_0\) and
has only boundedly many choices. After fixing it,
\(x=P(t)>0\) and \(A(t)\) are fixed, while
\(y=\mu x+\nu z\). At \(z=0\) the polynomial
\(A(t)+B(\mu x+\nu z)^q\) vanishes with nonzero
derivative. The vanishing follows from \(L_0(t)=0\);
the derivative is \(qB\nu(\mu x)^{q-1}\ne0\).
The quotient by \(z\) is nonconstant with nonzero
constant term. As above, its nonzero value can have
primes outside the given pool only on \(z\), so
such primes divide its constant term after clearing
denominators. The fixed-support count recovers \(z\)
at subpower cost.

The per-quotient statement can also be used with
varying orders. Orders of layer values equal to one
do not affect the complete letters. The other orders
have \(O(\log X)\) possible values apiece. This includes
\(q\) when \(y>1\) or some \(P_i(t)>1\), by the bound
\(q\le p_i\). When all those values are one the complete
letters are already fixed. If \(q\) is large,
\(P(t)\le X^{B_0/q}\) and \(y\le X^{B_0/q}\) count
the two radical parameters at subpower cost. If a
\(p_i\) is large, \(P_i(t)\le X^{B_0/p_i}\) counts \(t\)
at subpower cost, including when that layer value is
one. After fixing \(t\) and the actual \(A(t)\), use
\(F(z)=A(t)+B(\mu P(t)+\nu z)^q\) directly, at a
bounded \(q\). If \(F(0)\ne0\), all primes of \(F(z)\)
outside the fixed pool divide \(z\) and hence the
constant \(F(0)\) after clearing denominators.
If \(F(0)=0\), that root is simple and the fixed
constant-term argument applies to \(F(z)/z\).
Thus this branch needs no reconstruction from \(L_0\).
Only finitely many smaller orders need the factor
and resultant estimates. The degree and height of
\(L_0\) remain explicit in the quotient price below.

### Large excess in order two

When \(q=2\), \(r\ge3m\), and \(L\) is coprime to \(P\),
one may also use the quotient by \(z\) of the primitive
form of the full polynomial \(L\). The degree in its
height price is \(r\), which agrees with \(r_0\) if
\(L\) is squarefree. Denote the new fixed quotient by
\(k\). When it is nonzero write
\(u=\mu P+\theta_0 L\) with the corresponding fixed
nonzero rational \(\theta_0\), accounting for the
primitive scalar in \(L\). Then
\[
  F=L G,\qquad
  G=P^2+2B\mu\theta_0 P+B\theta_0^2 L
    =(P+B\mu\theta_0)^2+B\theta_0^2 M.
\]
The two summands in the last expression are coprime,
since \(P+B\mu\theta_0\) is coprime to \(P\) and all
the factors of \(M\) lie on \(P\). As \(r>2m\), the
degree of \(G\) is \(r\). Mason–Stothers gives at
least \(r-2m+1>m\) distinct roots of \(G\).
At any common root of \(G,L\), the coprimality of
\(P,L\) forces \(P+2B\mu\theta_0=0\). That
polynomial has at most \(m\) distinct roots. Thus a
factor of \(G\) remains off \(L\). Apply the same
resultant and support enumeration. The zero quotient
has the same proof as before.

## Quotient price and a sector with several layers

For \(n=X^L\) at fixed budget, suppose \(L_0\) has
coefficient height at most \(n^\theta\), \(|t|\le n^h\)
and \(z\ge n^{r_z}\) after the fixed-pool stripping.
Then (LC1) gives the normalized quotient price
\[
                               (\theta+r_0h-r_z)_+ .        \tag{LC3}
\]
Here the fixed-degree allowance absorbs the number of
coefficients. Alternatively, if the outside radical of
each nonzero \(k\) is at most \(n^b\) off an additional
prescribed polynomial-height pool, its price is at
most \(b\). Enumerate the outside support and then
its supported numerical values at subpower cost.
These estimates allow the zero quotient.

Thus the original row count is at most
\(X^\varepsilon n^{s+(\theta+r_0h-r_z)_+}\)
for a catalogue of normalized size \(s\) supplying
the carrier, coefficient and plane records, in the
cases of (LC2), or with the full-polynomial price
in the large-excess order-two case. It is enough that this total is
strictly below the tested normalized count exponent.
The catalogue pays for the carrier once and does not
pay separately for the exact layer roots. If
\(r_z\ge r_0h\), the only height price after the
catalogue is the primitive coefficient height
\(\theta\). For example, with \(q=2\) and valuation
orders \(2,3,4\) on three linear carrier forms, one
has \(m=r=3\) and \(r_0\le3\). The bound applies
even at order two and with arbitrary third-letter
valuations. If the three linear layer forms have
sizes \(n^{h+o(1)}\), this condition can use a
third radical of size at least \(n^{3h+o(1)}\).
The corresponding general consequence of a single
grouped root dictionary only prices the three
coordinates at \(h\) instead of \(3h\). The present
bound can absorb larger catalogue prices when
\(L_0\) is cheap.

The restriction to genuine excess \(r>0\) matters in
order two. For instance \(P(t)=t\), \(x=t\),
\(y=g-t\), and \(z=2t-g>0\) give
\(x^2-y^2=gz\) on the plane \(y=x-z\), with \(g\)
fixed in the prime pool. On squarefree coprime
eligible parameters the third value has residual
radical \(z\), and the only nonconstant factor of
the third value is \(z\). One instance is \(g=7\),
\(t=5\), giving the primitive row \(4+21=25\). This is the squarefree residual
branch of equal quadratic orders. For \(r>0\) and
\(q=2,\ r_0>m\), the statement uses the explicit
degree test instead of assuming that this branch
is absent. It makes no extraction claim for an
arbitrary diffuse family.
[Equal quadratic orders on radical planes](ABC_SHALLOW_REPAIR_RADICAL_PLANES.md#consequence-for-exceptional-quadratic-records)
· [Grouped root dictionary estimate (GD3)](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#group-dictionary-count)
