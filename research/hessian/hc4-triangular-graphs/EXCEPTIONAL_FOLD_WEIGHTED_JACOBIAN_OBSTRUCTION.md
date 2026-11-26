# A weighted Jacobian obstruction for exceptional polar folds

**Theorem (Longitudinal degree of an exceptional polar fold).**
Let a complex polynomial plane Keller map \((P,Q)\) have the
[exceptional polar fold](EXCEPTIONAL_POLAR_FOLD_WHOLE_CHARTS_AND_BIRATIONAL_IMAGES.md)
in its full finite normalization. In the original normalized polynomial
source coordinates \(a,b\), the fold chart is \(a=b^{-2}+bt\), and
\[
\max(\deg_a P,\deg_a Q)\ge18.
\]
There is no restriction on their \(b\)-degrees. If the maximum is
eighteen, an affine target change gives restrictions to the entire
fold line of degrees two and three, with lowest blocks proportional
respectively to the square and cube of \(b^{-3}(ab^2)^5(ab^2-1)\).
This is an admission obstruction in the actual source frame; the
unrestricted exceptional fold remains open.

## 1. A rational equation for the lowest block

Put \(u=ab^2=1+b^3 t\). Polynomial admission in both planes gives
the
[weighted intersection criterion](EXCEPTIONAL_FOLD_POLYNOMIAL_ADMISSION_AND_MARKED_TRACE.md):
in each block \(b^i p_i(u)\), the polynomial \(p_i\) is divisible by
\(u^{\max(0,\lceil-i/2\rceil)}(u-1)^{\max(0,\lceil-i/3\rceil)}\).
In particular both \(0,1\) are roots in every negative-weight block.

Let one actual target coordinate \(P\) have lowest block
\(b^{-m}p(u)\) with \(m>0\). There is a polynomial \(f(u)\) such that
\[
 m f' - \frac{p'}{p}f=J,\qquad
 J=\operatorname{Jac}_{a,b}(P,Q)\in\mathbb C^*.
\tag{1}
\]
Moreover \(f\) has a simple zero at every root of \(p\), and if
\(\deg f\ge2\), then
\[
 \deg p=m\deg f .
\tag{2}
\]
A root of \(p\) cannot have multiplicity exactly \(m\).

For context, the positive-weight Hamiltonian in
[Guccione–Guccione–Valqui, *On the shape of possible counterexamples to the Jacobian Conjecture*](https://arxiv.org/abs/1401.1784v3),
Theorem 2.6 and Proposition 2.11, gives a closely related initial-form
identity. In our coordinates \((a,b)\) its direction is \((2,-1)\),
of positive total weight one, and its Hamiltonian of weight one has
the form \(b^{-1}f(u)\) up to sign and constant normalization.
Here we derive the equation locally in the formal Laurent field,
then combine it with simultaneous polynomial admission and the
actual fold-line map.

The bracket in \(\mathbb C(u)((b))\) is the one induced from the source
Jacobian, namely \(b^2\operatorname{Jac}_{u,b}\). For individual blocks
it satisfies
\[
 \{b^i p,b^j q\}=b^{i+j+1}(j p'q-i p q').
\tag{3}
\]
Write \(P=b^{-m}p\,C\) with \(C\in1+b\mathbb C(u)[[b]]\).
Suppose the lowest remaining block of \(Q\) is \(b^j q\) with \(j<m-1\).
Its bracket with \(b^{-m}p\) is the only contribution at weight
\(j-m+1<0\), so that coefficient is zero. This gives
\(m q'/q=-j p'/p\), and \(q^m p^j\in\mathbb C^*\).
The formal Laurent element
\[
 T=b^j q\,C^{-j/m}
\]
has the same lowest block and satisfies \(T^m=(q^m p^j)P^{-j}\).
It therefore commutes with \(P\). Subtract it from \(Q\) in the formal
Laurent field. This operation raises the lowest weight by at least
one while retaining the exact Jacobian. It need not be a polynomial
target change; it is used solely to compute the first possible
nonzero bracket coefficient. After finitely many such subtractions,
the lowest weight must be exactly \(m-1\), since a larger weight
would give zero constant term in the bracket. With \(r(u)\) its
rational coefficient, (3) gives
\(m p r'+(m-1)p'r=J\), which is (1) for \(f=p r\).

At a finite root \(\alpha\) of \(p\) of multiplicity \(h>0\), if
\(f\) has order \(s\) there, the potential leading term on the left
of (1) has order \(s-1\) and coefficient \(ms-h\) times the local
leading coefficient of \(f\). For \(s\le0\) this coefficient is
nonzero, giving a pole. For \(s\ge2\) the left side vanishes there.
Thus \(s=1\), and the nonzero constant term requires \(h\ne m\).
At any other finite point a pole in \(f\) would give an uncancelled
pole of higher order from \(m f'\). Hence \(f\) is polynomial.
A zero of \(f\) away from the roots of \(p\) must also be simple
by (1). At infinity the potential leading term in (1) has degree
\(\deg f-1\) and coefficient \(m\deg f-\deg p\) times the leading
coefficient of \(f\). For \(\deg f\ge2\) this gives (2).

In this fold chart \(p\) has at least the two distinct roots \(0,1\),
so (2) yields \(\deg p\ge 2m\) for each coordinate with a negative
lowest weight. This is the obstruction that uses the Laurent
Jacobian equation in addition to blockwise admission.

## 2. Applying the equation on the whole fold chart

The fold line is the whole \(b=0\) line of the chart, and the
target restriction generates its function field \(\mathbb C(t)\)
by the whole-chart theorem. Its image is an irreducible component
of the nonproper value set. A restriction of degree one in an affine
target combination would make that image a polynomial graph,
contrary to the no-affine-line component theorem in
[Chau's remark, printed page 43](../../jacobian/affine-complements/sources/global-boundary-topology/chau-2003-two-remarks.pdf).
Neither coordinate restriction can be constant, since birationality
would make the other linear. If both degrees are at most two,
cancelling common quadratic terms gives either a linear restriction
or a field generated by one quadratic. This shows that at least one
coordinate has restriction degree \(j\ge3\).

Every \(t^j\) term on the chart line comes from weight \(-3j\),
because \(u-1=b^3t\). Thus the coordinate just chosen has
\(m\ge3j\ge9\), and by (2) its lowest coefficient has degree
at least \(2m\ge18\). The \(a\)-degree of a polynomial is the
maximum of the \(u\)-degrees of its \(b\)-blocks, giving the bound.

Suppose the maximum \(a\)-degree is eighteen. Equations (1)–(2)
now bound the magnitude of every negative lowest weight by nine,
also after an affine target change. Accordingly all restrictions
have degree at most three. After cancelling a common cubic term
where necessary, take degrees two and three for \(P,Q\).
For \(Q\), the lowest weight equals \(-9\), with order exactly
three at \(u=1\). Its coefficient has degree eighteen and exactly
two distinct roots, by (2). Those roots are \(0,1\), so the
coefficient is proportional to \(u^{15}(u-1)^3\).

Let the lowest block of \(P\) have weight \(-m\) with \(6\le m\le9\).
The bracket of the two lowest blocks has strictly negative weight
and must vanish. Their two coefficient polynomials therefore have
multiplicities in ratio \(m:9\) at every root. At \(u=1\) this
forces an order \(m/3\) for the coefficient of \(P\).
If \(m=9\), that block would contribute a cubic term on the fold
line. Thus \(m=6\), and the coefficient is proportional to
\(u^{10}(u-1)^2\), as asserted.

The proof restricts the lowest blocks in the original frame without
bounding the other weights from above. It supplies no elimination
of all remaining blocks, and the full-normalization
[marked double pole](EXCEPTIONAL_FOLD_POLYNOMIAL_ADMISSION_AND_MARKED_TRACE.md)
still measures the one-velocity obstruction on any actual survivor.
