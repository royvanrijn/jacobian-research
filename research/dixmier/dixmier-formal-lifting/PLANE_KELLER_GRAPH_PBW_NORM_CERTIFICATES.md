# Computing the low filtration of a plane Keller graph

Let \(k\) be a perfect field of characteristic \(p>0\) and \(f=(f_1,f_2)\)
a Jacobian-one plane map with \(\deg f<p\). Put
\(F_i=f_i^{[p]}\), the coefficientwise Frobenius of \(f_i\), \(d_i=\deg F_i\),
and \(Z=k[X,Y]\). The actual splitting graph module \(M\) in
[Kanel-Belov–Fedorov v2](sources/kanel-belov-fedorov-2610.33919v2.pdf),
Propositions 3.2 and 5.2, has division-ring idealizer
\[
 B=Z\langle U,V\rangle/
 ([V,U]-1,\ U^p-F_1(X,Y),\ V^p-F_2(X,Y))
   \simeq\operatorname{End}_{A_1(k)}(M)^{\mathrm{op}} .
\]
Its \(Z\)-basis is \(U^iV^j\), \(0\le i,j<p\).
Write \(F_tB\) for degree at most \(t\) in the intrinsic
reduced-norm filtration in \(X,Y\), and let \(T=\operatorname{gr}B\)
in \(S=k[\xi,\eta]\). The embedding uses the Bernstein
filtration of the source Weyl division ring, with \(X=x^p\), \(Y=y^p\).
This is the filtration in the
[two-degree bound](PLANE_KELLER_GRAPH_BW_FILTRATION.md).

## A uniform central-coefficient bound

**Theorem PlaneGraphCentralPBWCoefficientBound.** For
\(b=\sum_{0\le i,j<p}a_{ij}U^iV^j\in F_t B\), every nonzero
coefficient satisfies
\[
 p\deg a_{ij}\ \le\
 t+(p-1-i)d_1+(p-1-j)d_2-2(p-1).                         \tag{1}
\]
For fixed \(t\) and bounded \(d_1,d_2\), the central coefficient
degrees have a common bound independent of \(p\).

First compute reduced trace over any splitting field of the Weyl
algebra with \(p\)-center fixed. For central \(p\)-roots \(\alpha,\beta\)
of the two powers, represent the generators by
\(\alpha+s\) and \(\beta+\partial_s\) on the truncated polynomial module over the splitting field.
The diagonal of \(s^a\partial_s^b\) is zero unless \(a=b\).
For \(a=b\) its trace equals
\(\sum_{m=0}^{p-1}m(m-1)\cdots(m-a+1)=a!\binom p{a+1}\).
Thus the trace of the PBW basis is zero except at \((p-1,p-1)\),
where it is \(-1\). The trace formula works over any central
base change. By cyclicity, the trace pairing of \(U^{i'}V^{j'}\)
with \(V^{p-1-j}U^{p-1-i}\) is
\(\operatorname{Trd}(U^{p-1-i+i'}V^{p-1-j+j'})\), which is
\(-1\) exactly at \((i',j')=(i,j)\). Indeed reducing either exponent
modulo \(p\) can reach \(p-1\) only if it was already \(p-1\).
Consequently
\[
 a_{ij}=-\operatorname{Trd}\!\left(bV^{p-1-j}U^{p-1-i}\right).  \tag{2}
\]

In the source Weyl algebra \(A_1(k)\), the same trace selects the
coefficient of \(x^{p-1}y^{p-1}\) over \(k[x^p,y^p]\). For any
nonzero trace this gives
\[
 p\deg_{X,Y}\operatorname{Trd}(q)+2(p-1)
          \ \le\ \deg_{\mathrm{Bern}}q.                    \tag{3}
\]
Every source division-ring element can be written with a central
polynomial denominator, since the Weyl algebra is Azumaya and a
domain. Its fraction-center algebra is already the division ring.
Clearing the central denominator in (3) decreases both sides by the
same degree upon division. Reduced trace agrees under the
specified-center identification with the idealizer. Also
\(\deg_{\mathrm{Bern}}U=d_1\), since
\(\operatorname{Nrd}(U)^p=\operatorname{Nrd}(U^p)=F_1^p\);
the corresponding degree of \(V\) is \(d_2\). Applying (3) to
(2) proves (1).

## Exact filtration from an explicit matrix

**Theorem PlaneGraphPBWNormElimination.** For every supplied \(k,f,p\)
as above and \(t\ge0\), the filtered increments through \(t\) can be
computed with exact arithmetic and inverse Frobenius in \(k\) by the following finite
procedure. It uses only the prescribed graph presentation of \(B\)
and determinants of size \(p\).

Introduce two variables \(x_0,y_0\) and send \(X\mapsto x_0^p\),
\(Y\mapsto y_0^p\). On the free rank-\(p\) module with basis
\(1,s,\ldots,s^{p-1}\), use the matrices
\[
 U=f_1(x_0,y_0)I+s,\qquad V=f_2(x_0,y_0)I+\partial_s
 \quad\text{on } k[x_0,y_0][s]/(s^p).                      \tag{4}
\]
These satisfy the three graph relations. They split the pulled-back
algebra. In fact \(s\partial_s\) is diagonal with the distinct
entries \(0,\ldots,p-1\). Its polynomial projectors, together with
\(s\) and \(\partial_s\), generate the matrix units. Thus the
determinant of the matrix of \(b\) in (4) is
\(\operatorname{Nrd}(b)(x_0^p,y_0^p)\).

For \(b\ne0\) in the idealizer, its determinant has total degree
\(p\deg_B b\). Its homogeneous top form is \(\sigma(b)^p\), with
\(\sigma(b)\in S\) evaluated at \((x_0,y_0)\). To verify the second
claim, first take a source Weyl element \(a\). The matrix of its
left multiplication in the center basis \(x^iy^j\) has Bernstein
weights \(i+j\) on the basis. Its highest possible weighted matrix
is multiplication by \(\sigma(a)\) over the graded center. After
localizing that center the determinant is \(\sigma(a)^{p^2}\), the
norm in the purely inseparable extension
\(k(\xi,\eta)/k(\xi^p,\eta^p)\) of degree \(p^2\). This determinant
is nonzero. Compare it with the top term of the identity \(\det L_a=\operatorname{Nrd}(a)^p\) in the original
positive-characteristic center. Comparing unique \(p\)th roots
gives the claim. Division by a central denominator extends it
to the idealizer. The determinant's top-form Frobenius root
therefore gives both degree and the actual principal symbol
using the specified source center coordinates.

Form the finite \(k\)-span \(C_t\) of \(X^rY^sU^iV^j\) with
\(0\le i,j<p\), \(r,s\ge0\), and \(p(r+s)\) at most the right
side of (1). There are no such monomials when that right side is
negative. Equation (1) shows \(F_tB\subset C_t\). Each element of
\(C_t\) has filtered degree at most
\(t+(p-1)(d_1+d_2-2)\) by the degrees of \(X,Y,U,V\).

Use any homogeneous monomial order within each symbol degree.
For each PBW monomial of \(C_t\), use (4) to compute its symbol,
and reduce against a stored element whenever their leading symbol
monomials agree. Match coefficients before subtraction. Each reduction
strictly decreases the ordered pair of degree and leading symbol
monomial. A nonzero unreduced row is stored with its pivot normalized.
The resulting rows span \(C_t\) and have distinct pivots. In a linear
combination the highest pivot cannot cancel, because all other rows
have strictly smaller pivots. Hence the rows of degrees at most \(j\)
are a basis of \(F_jB\) for every \(0\le j\le t\). In particular, one
run at \(t=n+1\) tests both degree increments in the two-degree bound.
A replay certificate records the input graph, ordered PBW monomials,
exact row operations and the determinant top forms of the output
rows. The graph relations and all filtration comparisons occur in the
given finite algebra. The [exact replay](verify_plane_graph_pbw_norm_certificates.py)
uses the graph relations and filtered elimination over \(\mathbb F_2\),
\(\mathbb F_3\), and \(\mathbb F_5\), including a composed shear of degree
four with nontrivial pivot cancellations.

## A single odd degree in sufficiently large characteristic

**Corollary PlaneGraphOddDegreeFullness.** If \(m\ge 3\) is odd, \(p>m\)
and \(T_m=S_m\), then \(T_{2m-2}=S_{2m-2}\) and \(T_t=S_t\) for
\(t\ge (m-1)(2m-3)\). In particular a common such \(m\) for a
collection gives a common bound on BW, and the certificate above
need only reach degree \(m\).

The Bernstein principal-symbol bracket of source fractions is the
usual Poisson bracket, up to the sign determined by \([y,x]=1\).
The commutator drops degree by at least two on Weyl polynomials.
The Leibniz rule and \([q^{-1},r]=-q^{-1}[q,r]q^{-1}\) extend this
to the division ring, so \(T\subset S\) is a graded Poisson
subalgebra. The brackets of \(S_m\) span \(S_{2m-2}\). In fact
the monomials with \(\xi\)-exponents \(i,j\) in degree \(m\) give,
up to sign, \(m(i-j)\xi^{i+j-1}\eta^{2m-i-j-1}\). Every target
exponent occurs with \(0\le i,j\le m\), \(i\ne j\). Its coefficient
is nonzero since \(p>m\). The positive degrees \(m\) and \(2m-2\)
are coprime. For coprime positive integers \(a,b\), every integer
at least \((a-1)(b-1)\) is \(ua+vb\) with \(u,v\ge0\), by taking
\(v\in\{0,\ldots,a-1\}\) in the prescribed residue modulo \(a\).
Products span the resulting homogeneous pieces. Formula
\(\lambda=2\nu\) then bounds BW by half of
\(\sum_{i=1}^{(m-1)(2m-3)-1}(i+1)\).

## Uniformity in each fixed characteristic

**Theorem PlaneGraphFixedCharacteristicUniformity.** Fix \(d\ge1\)
and a prime \(p>d\). Among all Jacobian-one plane maps of degree at
most \(d\) over perfect fields of characteristic \(p\), the BW
invariants of the actual graph modules are uniformly bounded by a
number depending only on \(p,d\).

The graph coefficients range over an affine finite type scheme
\(P_{p,d}/\mathbb F_p\) given by the coefficients of the polynomial
equation \(\det Df=1\). For any fixed degree \(t\), the locus
\(T_t=S_t\) is constructible on this scheme. Use the finite list
of PBW monomials for \(C_t\) with both \(d_1,d_2\) replaced by \(d\).
Take \(t+1\) universal linear combinations of these monomials.
In their determinants (4), require all terms of degree greater
than \(pt\) to vanish. The remaining coefficient matrix, taking
coefficients of \(x_0^{pi}y_0^{p(t-i)}\) for \(0\le i\le t\)
in the degree-\(pt\) parts of these determinants, must have rank
\(t+1\). These are polynomial equations and a nonvanishing minor
on a finite type scheme over \(P_{p,d}\). By
[Chevalley's theorem](https://stacks.math.columbia.edu/tag/00FE)
the image is constructible. The determinant symbol formula proves
that the image consists exactly of the points with a full
degree-\(t\) piece after passage to an algebraic closure of the
residue field. The degrees in the filtered elimination and their
pivot counts persist under extension of perfect fields, since the
same elimination and pivot polynomials work after extension.

For each point of \(P_{p,d}\), including each generic point,
pass to an algebraic closure of its residue field and apply
Kanel-Belov–Fedorov Proposition 5.2. Finite colength gives fullness
in two consecutive positive degrees. Thus the intersections of
the constructible full loci in degrees \(n,n+1\), \(n\ge2\), cover
\(P_{p,d}\). A constructible cover of a Noetherian scheme has a
finite subcover. Indeed one can choose a set covering the generic
point of each irreducible component, retain the dense relatively
open parts within the chosen sets, and apply Noetherian induction
to the remaining proper closed subset. The two-degree bound for
the finitely many indices in this subcover yields the assertion.

**Corollary PlaneGraphLargeCharacteristicSingleDegree.** The uniform
numerical bound of Kanel-Belov–Fedorov v2 Conjecture 1.9 is equivalent
to the following statement. For every \(d\ge1\), there is an odd
\(m\ge3\) such that \(T_m=S_m\) for all degree-at-most-\(d\) Keller
graphs over finite fields of characteristic \(p>\max(d,m)\).

If the uniform BW bound holds, graded fullness above twice the LM
colength gives such an \(m\). Conversely, odd-degree fullness
gives a common BW bound at all primes \(p>\max(d,m)\). The finitely
many intervening primes \(d<p\le m\) have a common bound by the
fixed-characteristic theorem.

The finite certificate works directly with the given graph algebra.
The PBW index set still grows with \(p\). The common odd degree
across all sufficiently large primes remains open, as does
boundedness for one faithful arithmetic family per complex map in
the conditional transfer.
