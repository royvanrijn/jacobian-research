# Counting with two short root progressions and a target divisor

On a generalized Fermat carrier, two original residue lists can
shorten the two root variables simultaneously. A third packet on the
actual target then enters the affine determinant estimate before the
curve count. The variable used in the curve projection is a linear
form in the two short progression indices. The coefficient-uniform
geometry of the original Fermat equation therefore gives a bound in
their shorter widths.

## The simultaneous arithmetic count

Fix \(p,q,r\ge2\), \(B>0\), and \(0<\omega_0<1\).
Use \(X\ge2\). Let \(\alpha,\beta,\gamma\) be nonzero integers of
absolute value at most \(X^B\), and consider positive integral
solutions of
\[
                    \alpha x^p+\beta y^q=\gamma z^r,
             \qquad \gcd(x,y,z)=1.
\]
All three roots are at most \(X^B\). Prescribe intervals \(I_x,I_y\)
for \(x,y\) respectively, one residue class modulo \(Q_x\) in the
first interval and one modulo \(Q_y\) in the second, and \(2\le Z\le
X^B\) with \(z\le Z\). The two moduli are positive integers at most
\(X^B\). Let \(T_0\) be the least integer exceeding \(\max(p,q)/\omega_0\).
Put
\[
 U=T_0+3+\frac{|I_x|}{Q_x},\quad
 V=T_0+3+\frac{|I_y|}{Q_y},\quad H=\max(U,V).
\]
One may intersect the intervals with \([1,X^B]\) before computing
these widths. Impose a further fixed divisor \(M\mid \gamma z^r\),
\(1\le M\le X^B\). Let \(M_0\) be its full part outside the primes
dividing \(\alpha\beta\), the primes dividing \(\gcd(Q_x,Q_y)\),
and the primes at most \(T_0\). Define
\[
 \begin{split}
 S&=\max(p\log U,q\log V),\\
 \mathcal K&=
     \max\left(1,\ \exp\left[
         \sqrt{\frac{\log U\log V}{r}}
                 \left(1-\frac{(1-\omega_0)\log M_0}{2S}\right)
                           \right]\right),\\
 \nu&=\max\left(\frac{2}{\sqrt{\max(p,q)}},\frac1r\right).
 \end{split}
\]
For every \(\varepsilon>0\), the number of solutions is
\[
               O_{p,q,r,B,\omega_0,\varepsilon}
                         \left(X^\varepsilon
                                (\mathcal K^2+\mathcal K H^\nu)\right).
                                                               \tag{DF1}
\]
The estimate is uniform in all coefficients,
root intervals, residue classes and moduli within the indicated
bounds. The widths include the effect of both root residue classes.

### Keeping a good part of the third modulus

The result is vacuous if a prescribed root class misses its interval.
Otherwise choose initial representatives \(x_0,y_0\) from the respective
classes in their intervals. Write
\(x=x_0+Q_x u,\ y=y_0+Q_y v\). Translate each representative by an
additional multiple \(jQ_x\) or \(jQ_y\) with \(0\le j<T_0\), taking
\(T_0\) to be an integer. Use the translated representatives again as
\(x_0,y_0\). The indices still satisfy \(|u|\le U, |v|\le V\).
At a prime dividing \(M_0\), at least one of \(Q_x,Q_y\) is a unit.
As the translation on that side varies, with the other one fixed,
\(\alpha x_0^p+\beta y_0^q\) is a nonzero polynomial of degree at most
\(\max(p,q)\) modulo the prime, by its unit leading coefficient.
The \(T_0\) possible shifts are distinct modulo this prime.
Average over both translations, weighting each prime by its full
\(v_\ell(M_0)\log\ell\). Some choice leaves a full divisor
\(M_*\mid M_0\) with
\[
 \log M_*\ge (1-\omega_0)\log M_0,\qquad
                    \gcd(M_*,\alpha x_0^p+\beta y_0^q)=1.
\]
This step applies to composite divisors and does not require any prime
of positive normalized logarithmic size. If \(M_0=1\), use \(M_*=1\).

### The affine surface and the short curve projection

Apply Browning and Verzobio's affine-surface theorem with side
condition to the surface
\[
  f(u,v,z)=\alpha(x_0+Q_xu)^p+\beta(y_0+Q_yv)^q-\gamma z^r=0
\]
and to \(g(u,v)=\alpha(x_0+Q_xu)^p+\beta(y_0+Q_yv)^q\) modulo \(M_*\).
The surface is absolutely irreducible by the irreducibility of
the original Fermat polynomial and the invertible affine change.
The common divisor condition of that theorem holds by the chosen
translations. Take a lexicographic order with \(z\) first, so that
the leading monomial is \(z^r\) and its box weight is \(Z^r\).
The value \(S\) above is the maximum log box weight on the side
polynomial, since both \(u^p\) and \(v^q\) occur in it.
The determinant parameter of the theorem is then at most
\(\mathcal K\). It provides \(O(X^\varepsilon\mathcal K)\)
auxiliary polynomials, coprime to \(f\), of degrees \(O(\log X)\),
apart from \(O(X^\varepsilon\mathcal K^2)\) exceptional points.
The exponents \(\varepsilon\) can be decreased at this step
to leave room for subsequent height-subpower factors.
[Browning and Verzobio, Counting integer points on affine surfaces with a side condition, Theorem 1.1 and equation (5)](https://arxiv.org/html/2408.11453v2#S1)

For each auxiliary polynomial, take its resultant with \(f\) in
\(z\), or the polynomial itself if independent of \(z\).
The result is nonzero by coprimality and absolute irreducibility.
The sum of the degrees of its irreducible plane factors is
\(O_{p,q,r,\varepsilon}(\log X)\). Count points over each such
factor \(c(u,v)\), of degree \(d\). At every fixed \((u,v)\)
there is at most one positive \(z\). A factor irreducible only over
\(\mathbb Q\) has \(O(d^2)\) rational points, by intersecting
distinct conjugate components. For an absolutely irreducible factor
with \(d\ge\lceil\sqrt{\max(p,q)}/2\rceil\), count its integral
\((u,v)\) by the degree-uniform affine curve theorem. This gives
\(O(d^2 H^{2/\sqrt{\max(p,q)}}(\log H)^{O(1)})\).
[Binyamini, Cluckers and Novikov, Theorem 2](https://arxiv.org/html/2302.04209v2#S5)

For the remaining factors, use their original-coordinate curve
\(c((x-x_0)/Q_x,(y-y_0)/Q_y)=0\). It has the same degree and
absolute irreducibility. If vertical or horizontal, the fixed
original root is positive. The equation on the other root and
\(z\) is an absolutely irreducible shifted power curve because
the binomial with its nonzero constant term has simple roots.
The lopsided integral curve bound uses the pure \(z^r\) term
and gives \(O(X^\varepsilon H^{1/r})\).
For comparable zero-and-pole supports of \(x,y\) on the
normalized plane curve, the divisor argument in Browning and
Verzobio's Lemma 2.10 gives \(O(X^\varepsilon)\) primitive
points in the original root box. The relevant fixed-degree plane
polynomial has coefficients of polynomial height unless there
are already \(O(1)\) integral points, by their Lemma 2.4.
These bounds thus apply without assuming primitive progression
indices.
[Heath-Brown, lopsided curve estimate from Theorem 15](https://ora.ox.ac.uk/objects/uuid%3Ac5ab1264-e082-4664-8c45-ae6d7b6f81a1/files/m73013a42d9d1586d603c5a67baef1e0d)
· [Browning and Verzobio, Sums of three powerful numbers, Lemmas 2.4 and 2.9–2.10](https://arxiv.org/html/2608.24512v1#S2.SS2)

For every other remaining factor, \(4d^2<\max(p,q)\).
Browning and Verzobio's function-field irreducibility lemma
then gives degree \(r\) for the extension adjoining \(z\) to the
function field of the absolutely irreducible plane curve.
This application is to the original \(x,y\) functions in the
function field. The affine change to \(u,v\) leaves that field
and the extension unchanged. Now take the generic linear
projection in the **short** coordinates \(u,v\).
For a suitable integer \(j=O_{p,q,r,d}(1)\), put \(t=u+jv\).
The leading homogeneous part of \(c\) ensures
\([\,\overline{\mathbb Q}(u,v):\overline{\mathbb Q}(t)\,]=d\),
where these are function fields of the curve. The primitive
element argument over \(\overline{\mathbb Q}(z)\) ensures
\(\overline{\mathbb Q}(t,z)=\overline{\mathbb Q}(u,v,z)\).
Both properties exclude only boundedly many \(j\), as in
Browning and Verzobio's Lemma 2.17. Thus \(z\) has degree \(dr\)
over \(\overline{\mathbb Q}(t)\).

The plane relation on \(t,z\) has degree bounded in terms of
\(p,q,r,d\) by elimination in \(v\). If the relation is not
defined over \(\mathbb Q\), it has only boundedly many rational
points by Bézout. Otherwise its coefficient height is
polynomial unless the integral point count is already bounded,
by Lemma 2.4 again. Use \(|t|\ll H\), \(z\le Z\) and the
nonzero term of \(z\)-degree \(dr\) in the lopsided curve
bound. This counts \(O(X^\varepsilon H^{1/(dr)})\) pairs
\(t,z\). Each pair comes from at most \(d\) points \(u,v,z\),
because the chosen projection has constant nonzero leading
coefficient in the degree-\(d\) equation for \(v\).
[Browning and Verzobio, Lemmas 2.11 and 2.17](https://arxiv.org/html/2608.24512v1#S2.SS3)

Summing over all resultant factors and auxiliary polynomials
absorbs their degree and number losses in \(X^\varepsilon\).
The exceptional point count and the per-curve bounds give
(DF1). In the generic projection, performing the affine change
on the surface alone would only shorten the determinant parameter.
Using the short index coordinates for the curve projection
shortens the per-curve bound as well.

## Strict primitive-quotient sectors

Suppose an original abc family in \(X/2<c\le X\) admits
\(a=\alpha x^p,\ b=\beta y^q,\ c=\gamma z^r\) with the input
data of (DF1) from at most \(n^f\) records, \(n=X^L\).
A record includes the coefficient tuple, the root residue
classes together with their two moduli, the target divisor, and
the boxes. If it instead gives the complete-letter residues
\(a\bmod Q_x,b\bmod Q_y\), and these letters are units there,
the root residue classes cost only \(X^{o_X(1)}\) per record.
Indeed each kernel for the \(p\)-th or \(q\)-th power map on
the unit group has size bounded by the corresponding power
raised to \(2\omega(Q_x)\) or \(2\omega(Q_y)\).
Primitivity of the complete row implies primitivity of the roots.
All counts are for such original rows.

As an example, fix \(p=q=16,\ r=2\) and suppose both index
widths are \(U,V\le n^{1/10+o_X(1)}\), \(M_0\ge X^\mu\)
for a fixed \(\mu>0\), and \(f\le1/4\).
Take all sufficiently small fixed \(L>0\).
Here \(S\le(16/10+o_X(1))L\log X\), so the retained
third divisor forces \(\mathcal K=1\). As \(\nu=1/2\),
the combined coefficient is at most \(1/4+1/20=3/10\).
Counting both indices individually costs as much as
\(1/4+1/5\). Even the same reduced-width determinant bound
without the side condition does not give \(1/3\) at the top
of these widths, since its exceptional-point exponent would
include \(2/(10\sqrt2)\) in addition to the record cost.
In this example the two index widths can contain power-many
integers. The count uses the target divisor in addition to the
two shorter progressions.

For a minimally repaired primitive quotient application,
assume also
\(\operatorname{rad}(c)\ge c^{(1/3-\xi)L}\) and that the
least target repair at \(M=\lceil K/L\rceil\) has height at
most \(Y\), for fixed \(K>0\). The same consequence holds if
the least common repair at equal order \(M\) has height at
most \(Y\). In either case \(c^{M(1/3-\xi)L}\le Y\),
since every prime of the original target occurs to exponent
at least \(M\) in the repaired target. If the record
hypothesis and the \(3/10\) bound hold on the dyadic original
shells, summing gives
\[
               O\left(Y^{\, (3/10)/(M(1/3-\xi))+\varepsilon}\right).
\]
For fixed \(0\le\xi<1/30\) this is a power saving under
\(Y^{1/M}\) for primitive quotients. No count is made over
additional common multipliers.
[Least common repair and the primitive quotient exponent](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#the-exponent-of-the-minimally-repaired-primitive-count)

The proof consumes both short root progressions, an affordable
power-carrier record, and a target divisor with positive numerical
weight after removing the support common to both root moduli.
It does not produce those records from arbitrary diffuse degree
profiles.
