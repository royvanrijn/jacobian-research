# Short algebraic features in the aggregate packet inverse

A collection of diffuse target packets can recover a short feature
of an original summand through a bounded-degree algebraic
correspondence. The feature need not parametrize a rational curve.
The other coordinates of a curve realizing the correspondence can
have larger heights.

## A finite algebraic correspondence

Use the primitive additive graph and the definitions of \(c_0,n,L\)
and the packet parameters \(\beta,\gamma,w,S\) in (AP1)–(AP3).
In particular a packet can be the product of full powers of many
primes whose individual normalized logarithmic sizes tend to zero.
Fix \(k\ge1\), \(u,\sigma,B\ge0\), a degree bound \(D\ge1\), and a
price bound \(\bar s\). A catalogue \(\mathcal F_X\) of at most
\(n^\sigma\) polynomials \(F(Y,T)\in\mathbb Z[Y,T]\) is allowed to
depend on \(X\). Each polynomial is irreducible over \(\mathbb Q\),
has positive degree in both variables, total degree at most \(D\),
and coefficients of absolute value at most \(X^B\).

Suppose each used \(A\) admits \(F\in\mathcal F_X\) and a rational
feature \(t\) with
\[
                    F(A,t)=0,\qquad H(t)\le X^u,
 \qquad H(a/b)=\max(|a|,b)\quad (\gcd(a,b)=1,\ b>0).
\]
If
\[
                       kw>2u,\qquad \sigma+k\bar s<c_0,
                                                               \tag{AL1}
\]
no bin of price \(s\le\bar s\) can meet the near-full edge
threshold in (AP2), with all its tolerances sufficiently small
in terms of the fixed margins. For integer features, \(kw>u\)
suffices. With these assumptions on every near-full bin, the
whole marked sector has \(O_L(n^{2c_0-\delta})\) edges for some
fixed \(\delta>0\), for every sufficiently small fixed \(L>0\).
One can use the default \(\bar s=\beta+\gamma+\eta\) from (AP2)
if (AL1) has a fixed strict margin at \(\beta+\gamma\).
The argument applies to either summand and to signed permutations
of the original three roles.
[Aggregate packet bins and common lists (AP2–AP3)](ABC_AGGREGATE_DEEP_PACKET_INVERSE.md#almost-disjoint-packets-at-a-common-summand)

Choose one witness \(F,t\) on each used \(A\). For fixed \(F\),
its resultant \(\operatorname{Res}_T(F,F_T)\) is a nonzero
polynomial in \(Y\), since \(F\) is irreducible and has positive
\(T\)-degree in characteristic zero. The Sylvester matrix has
bounded size and degree, so there are \(O_D(1)\) possible
exceptional \(A\)'s at which this resultant vanishes. The chosen
witnesses with exceptional \(A\) contribute at most
\(O_D(n^{\sigma+c_0+\kappa})\) edges by the opposite dictionary
bound. Remove those edges. Since \(\sigma<c_0\), this removal
costs only a fixed power-smaller term at the near-full scale.
Perform this deletion in the original bin. For sufficiently small
\(\eta\), a bin originally meeting (AP2) still has at least
\(n^{2c_0-\tau-\eta/50}\) edges. The upper bounds on its packet
dictionary, target multiplicities and residue lists are preserved.
The proof of (AP3) applies with this edge bound, since after dropping
vertices of degree less than \(n^{c_0-\eta/4}\) it still gives at
least \(n^{c_0-\eta/10}\) different \(A\)-values.

For each remaining \(A\), let \(Q_A\) be the product of the
denominator of \(t\) and the absolute numerator of the nonzero
reduced rational \(F_T(A,t)\). This has bounded polynomial height,
uniformly in all witnesses. In fact, writing \(t=a/b\), one
can clear \(F_T(A,t)\) by \(b^D\), and the resulting integer is
\(O_D(X^{B+D(1+u)})\). The bound on \(b\) then gives
\(Q_A\le X^{B_*}\) for a fixed \(B_*\).

Choose \(\xi>0\) such that \(k(w-S\xi)>2u\), or \(>u\) in
the integer case. By (AP3), one list \(M_1,\ldots,M_k\) is
admissible for at least \(n^{c_0-4k\eta}\) distinct \(A\)'s.
Their possible vectors of \(A\)-residues modulo the list
have cardinality at most \(n^{k(s+\kappa+\eta)}\) by (AP2).
For each \(A\), remove from every \(M_i\) all powers of
primes dividing \(Q_A\), leaving \(M_i'\). The vector of
stripped moduli has \(X^{o_X(1)}\) possibilities, by the
divisor bound applied to the fixed original list. For any
such vector used here its least common multiple \(M'\) satisfies
\(\log_X M'\ge k(w-S\xi)\) by (AP3).

Fix \(F\), a vector of stripped moduli and a vector of
\(A\)-residues modulo the original list. For
\(p^e\parallel M'\), choose one member of the stripped list
attaining depth \(e\). The prescribed \(A\)-residue \(r\)
determines the polynomial \(F(r,T)\) modulo \(p^e\).
The feature is \(p\)-integral and its reduction is a simple root
of this polynomial modulo \(p\), because of the support stripped
from \(Q_A\). There are at most \(D\) such roots. Each lifts
uniquely to depth \(e\): at every successive depth the new
base-\(p\) digit solves one linear equation with unit coefficient
\(F_T(r,t)\bmod p\). Thus the possible \(t\)-classes modulo
\(M'\) number at most \(D^{\omega(M')}=X^{o_X(1)}\).
Two rational features of height at most \(X^u\) agreeing
modulo \(M'\) have a cross-difference of absolute value at
most \(2X^{2u}\) divisible by \(M'\), so are equal. For integral
features use their difference and the bound \(2X^u\).

For any fixed rational \(t\), \(F(Y,t)\) is a nonzero polynomial.
Otherwise \(T-t\) divides \(F(Y,T)\), contrary to irreducibility
and positive \(Y\)-degree. Thus each recovered feature yields
at most \(D\) possibilities for \(A\). Summing over the catalogue,
residue lists and stripping choices gives at most
\[
                   n^{\sigma+k(s+\kappa+\eta)+o_X(1)}
\]
distinct \(A\)'s in the common list. The strict price margin
contradicts the (AP3) lower bound, with \(\eta\) chosen small.
All residue and branching estimates are uniform in the
catalogue as it moves with \(X\). The discarded branching
edges and the choice of \(\eta\) can be absorbed by choosing
\(\delta\) smaller than the resulting fixed margins before
taking \(X\) large at each fixed \(L\).

## Scope of the transfer

For example, the short feature may be the \(t\)-coordinate
of \(y^2=P(t)\), where \(P\) is a specified nonconstant squarefree
polynomial of bounded degree. Suppose the original summand is
\(A=\rho\,t^r y^{2j+1}\), for fixed \(r,j\ge0\).
The relation
\[
                  A^2=\rho^2 t^{2r}P(t)^{2j+1}
\]
gives an irreducible bivariate polynomial after clearing
denominators and taking primitive part. Indeed its right side
is not a square in \(\mathbb Q(t)\).
This applies to a Pell norm model or to a curve of positive
genus, with no separate bound on \(H(y)\). The polynomial
and coefficient records still have to fit the size budget
in (AL1). A bounded-degree polynomial image or power class
is the case in which the correspondence is linear in \(Y\).

The input here is the existence of an affordable algebraic
relation to a short feature on an actual original letter.
A collection of deep moduli with few occupied classes does not
by itself construct such a relation. One-letter code dictionaries
can retain diffuse compact valuations, the radical cap and
avoidance of prescribed pools without an empirical packet
discount. To use the transfer in the unrestricted problem
one must force a short feature or decrease the effective
price using the two other original letters.
[Target-only packet price obstruction (AP5)](ABC_AGGREGATE_DEEP_PACKET_INVERSE.md#a-one-letter-obstruction-to-compulsory-packet-compression)
