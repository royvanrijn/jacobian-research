# Exact parameter exclusions in the even stripped quintic space

Over a characteristic-zero field, let
\[
 f=P(y,t)+xt^2(q+B t^2),\qquad
 P=\sum_{j+r\leq5,\ r\ {\rm even}}p_{jr}y^jt^r,\qquad
 (A,C,D)=(p_{14},p_{32},p_{50}).
\]
Moments are \(M_m=\Lambda(f^m)\) for the
[normalized sphere functional](GLOBAL_SPHERE_PURITY_RECOGNITION.md).
Purity is the vanishing of all positive moments. Eventual purity is
their vanishing for every sufficiently large order.

The even \(D=0,\ BC\ne0\) chart has empty pure-parameter sets for
normalized \(A\) and \(T=A^2/(BC)\), both with monic squarefree
polynomial \(1\). The full \(BD\ne0\) chart also has empty parameter
image, hence unit reduced elimination ideal in the normalized
parameters \(A,C\). In fact each phase in either chart has infinitely
many nonzero moments. This proof uses the leading binary form.
An independent finite computation excludes all nonzero \(P\) in
two specified coefficient slices, with prefixes of orders \(13\)
and \(14\).

## Leading binary form and exact parameter sets

The angular substitution \(x\mapsto\rho x,\ y\mapsto\rho^{-1}y\)
and scalar multiplication by \((\rho B)^{-1}\) preserve moment
vanishing. Choose \(\rho^4=C/B\) when \(D=0,\ BC\ne0\).
This gives \(B=C=1\) and the normalized \(a=A/(B\rho^2)\), so
\(a^2=T\). If \(BD\ne0\), choose \(\rho^6=D/B\) to obtain
\(B=D=1\) and parameters \(a=A/(B\rho^2),\ b=C/(B\rho^4)\).
These are reversible operations over an algebraic extension.

The leading binary form in the original coefficients, obtained by
\(x=u^2,\ y=-v^2/2,\ t=uv\) on the cone, is
\[
 v^4\bigl(Bu^6-\tfrac A2u^4v^2-\tfrac C8u^2v^4-\tfrac D{32}v^6\bigr).
\]
Since \(B\ne0\), the root \(v=0\) has multiplicity exactly four.
When \(D\ne0\), all other roots correspond to nonzero values of
\(u/v\), and the multiplicity of each is at most three because
the sextic is cubic in \((u/v)^2\). If \(D=0,\ C\ne0\), the root
\(u=0\) has multiplicity two and all other roots have multiplicity
at most two. Thus every geometric multiplicity is strictly less
than five in either chart, including repeated-root and
vanishing-parameter strata. The
[stable leading-form theorem](STABLE_LEADING_SPHERE_PHASE_PERIOD_OBSTRUCTION.md)
excludes eventual purity. The geometric pure-parameter sets in the
[effective elimination theorem](EFFECTIVE_PRIME_FACE_PURE_PARAMETER_ELIMINATION.md)
are therefore empty on exactly the indicated families.
The same empty-set conclusion holds at each fixed eventual threshold.

For an even stripped phase of actual degree five with both angular
signs, the [Newton-face reduction](STRIPPED_QUINTIC_NEWTON_FACE_REDUCTION.md)
forces \(B\ne0\) and at least two of \(A,C,D\) nonzero. If \(D=0\)
this forces \(C\ne0\); if \(D\ne0\) the second chart applies.
Consequently the entire mixed-winding even stripped quintic space
is excluded. Pure one-sided phases with zero height core are safe by the
angular cutoff. The individual even proof concerns the stripped
support.

## An independent arithmetic check on repeated-root faces

Here is a second all-order exclusion of \(D=0,\ BC\ne0\) and of the
cubic discriminant-zero chart. The balanced Prym argument in the
\(D=0\) case originated in the pre-existing local working note on
that case. The rank-two ordinary-reduction theorem supplies the
fourfold input for both applications.

Work in the normalized coordinates and put
\(R(z)=1-az/2-bz^2/8-dz^3/32\), where \((b,d)=(1,0)\)
or \(d=1\) according to the chart. Take a number-field specialization
of the finitely generated coefficient algebra preserving the
required nonzero face coefficients and the exact branch-point
multiplicities. A single specialization preserves every tail identity,
as described in [the Newton-face proof](STRIPPED_QUINTIC_NEWTON_FACE_REDUCTION.md#8-coefficient-fields-known-cancellation-and-scope).
Include \(\zeta_5\) in the number field. The top Laurent face,
for weight \(h+i\) on \(x^i t^h\), is \(xt^4 R(t^2/x^2)\).
At a good absolute degree-one prime \(p=10k+1\), multiplying \(M_{2k}\)
by \(p\) leaves exactly \([z^k]R(z)^{2k}\) on reduction because the
top angular-constant height is \(10k=p-1\). All lower even
heights give prime-to-\(p\) denominators.

Put \(\widetilde R(z)=z^3R(1/z)\). On the normalized curve
\(w^5=\widetilde R(u^2)\), the character-one summand is
\(H^1(\mathcal O(-2))\) since all raw homogeneous branch exponents
are less than five. Its class \(w/u\) is anti-invariant for
\(u\mapsto-u\), hence belongs to the Prym of the quotient
\(w^5=\widetilde R(z)\). Frobenius there has coefficient
\[
 [u^{p-1}]\widetilde R(u^2)^{(p-1)/5}
   =[z^{5k}]\widetilde R(z)^{2k}=[z^k]R(z)^{2k}.
\]
Thus ordinary primes of the Prym give nonzero moments. The branch
completion has degree ten on the curve and five on its quotient.
For raw branch exponents \(e\) totaling \(5n\), the
character-\(i\) dimension is
\(in-1-\sum_e\lfloor ie/5\rfloor\), by the normalized cyclic-cover
eigensheaf computation in [the leading binary-form proof](SPHERICAL_LEADING_BINARY_FORM_MOMENT_OBSTRUCTIONS.md).

If \(d=0\) and \(a^2+2\ne0\), these branch exponents are
\((2,4,1,1,1,1)\) and \((1,1,1,2)\), yielding dimensions
\((1,2,2,3)\) and \((0,1,1,2)\). When \(d=1\) and the cubic has
one double root and one simple root, the lists are \((2,2,1,1,4)\)
and \((2,1,2)\), yielding \((1,2,1,2)\) and \((0,1,0,1)\).
Both Pryms therefore have \(\mathbf Q(\zeta_5)\)-signature
\((1,1,1,1)\), with Rosati conjugation from the polarized deck
action. The [rank-two ordinary-reduction proposition](RANK_TWO_CM_ORDINARITY_AND_QUINTIC_FACE_EXCLUSION.md)
applies and supplies density-one ordinary primes after finite
extension. When \(d=0\) and \(a^2+2=0\) the coefficient is
\(\binom{4k}{k}(-a/4)^k\), a unit at all such primes outside a
finite set. A triple root in the cubic case similarly gives
\(\binom{6k}{k}(-\gamma)^k\) for \(R(z)=(1-\gamma z)^3\) with
\(\gamma\ne0\). The cubic discriminant boundary is
\(4a^2b^2+8b^3-16a^3-36ab-27=0\), which equals
\(1024\operatorname{disc}(\widetilde R)\). The stable leading-form proof
above also handles every distinct-root cubic, without an
ordinary-prime assertion about its Prym sixfold.

## Rational moment equations and the triangular system

Normalize \(B=1\) by angular scaling, so \(q\) is the original
ratio \(q/B\). With \(M_1=0\) the general centered phase in this
chart is
\[
 \begin{split}
 f={}&xt^2(q+t^2)+D y^5+a_4y^4+y^3(a_3+bt^2)+y^2(a_2+b_2t^2)\\
     &+y(a_1+b_1t^2+at^4)+b_0(t^2-1/3)+c_0(t^4-1/5).
 \end{split}
\]
Set \(D=1\) or \(D=0,\ b=1\) for the two parameter charts.
The [exact rational generator](compute_even_stripped_quintic_moments.py)
enumerates all balanced angular tuples of monomial multiplicities.
For such a tuple let \(N\) be the number of positive winding
factors and \(H\) the total half \(t\)-degree. Its trace is
\[
 \Lambda((xy)^N t^{2H})=
 \frac{2^N N!(2H)!(N+H)!}{H!(2N+2H+1)!}.
\]
The exact moment is the sum of these traces multiplied by
the corresponding multinomial factors and coefficient monomials.
For instance in either parameter chart,
\[
 M_2=\frac{2a_1 q}{15}+\frac{2a_1}{35}
 +\frac{2b_1q}{35}+\frac{2b_1}{63}
 +\frac{4b_0^2}{45}+\frac{16b_0c_0}{105}+\frac{16c_0^2}{225}
 +\frac{2aq}{63}+\frac{2a}{99}.
\]

The first moment containing a coefficient of \(y^jt^{2r}\)
is \(M_{j+1}\), where it occurs linearly with coefficient
\[
 L_{j,r}=(j+1)\sum_{k=0}^j\binom jk q^{j-k}
                         \Lambda((xy)^j t^{2(j+r+k)}).
\]
This follows by using exactly \(j\) positive factors to balance that
negative factor. For \(j=1,2,3\), \(L_{j,0}\) and \(L_{j,1}\)
are coprime over \(\mathbf Q[q]\) by the checked Euclidean algorithm.
Choose \(U_jL_{j,0}+V_jL_{j,1}=1\) and substitute
\[
   (p_{j0},p_{j2})
    =(U_j C_j+L_{j,1}D_j,\; V_j C_j-L_{j,0}D_j).
\]
This change has determinant \(-1\). Successive substitutions
from \(M_2,M_3,M_4\) then eliminate \(C_1,C_2,C_3\) over
\(\mathbf Q[q]\), including on the vanishing loci of any single
\(L_{j,r}\). The \(j=3\) equation also recovers the parameter \(b\).
If \(b=1\) is prescribed, retain that equation on the transformed
coordinates.

The \(a_4\) coefficient in \(M_5\) is \((8/21879)\ell(q)\), with
\[
 \ell(q)=q^4+\frac{36q^3}{19}+\frac{198q^2}{133}
                   +\frac{1716q}{3059}+\frac{1287}{15295}.
\]
Away from its zeros, \(M_5\) eliminates \(a_4\) rationally and
leaves the explicit \(M_m\), \(m\geq6\), in \(D_1,D_2,D_3,b_0,c_0,q,a\)
on the \(B=D=1\) parameter chart. At a zero of \(\ell(q)\) one
must retain \(a_4\) and impose \(M_5\) after the first three
substitutions. The next section excludes the entire latter chart
using a finite prefix, and also treats the slice \(q=0\).

## Characteristic-zero finite certificates for two slices

Keep \(D\) variable and \(B=1\). The transformed positive
coefficients remain fixed under simultaneous angular and overall
multiplication \(x\mapsto s^{-1}x,\ y\mapsto sy,\ f\mapsto sf\).
This gives weight \(j+1\) to each coefficient of \(y^jt^{2r}\)
and weight zero to \(q\). Consequently \(M_m\) is homogeneous
of weight \(m\) in the centered \(P\) coefficients. This weighting
is essential to the characteristic-zero implication.

The [projective slices checker](verify_even_stripped_quintic_projective_charts.py)
checks these two exact special fibers with msolve's exact sparse
linear algebra.

| Ground condition | Residue field and value | Moment prefix | Eliminated coefficients over the residue field |
|---|---|---|---|
| \(q=0\) | \(\mathbf F_{65521}\), \(q=0\) | \(M_1,\ldots,M_{13}\) | \(a_1,a_2,a_3,a_4\) |
| \(\ell(q)=0\) | \(\mathbf F_{101}\), \(q=69\) | \(M_1,\ldots,M_{14}\) | \(a_1,a_2,a_3\) |

For each row the affine \(D=1\) chart has unit Gröbner ideal.
For \(D=0\) the exact leading ideals contain the following pure
powers of all the respective remaining variables.
\[
 \begin{array}{c|l}
 q=0 & b_2^8,\ b_1^{10},\ b_0^6,\ c_0^9,\ a^{10},\ b^6\\
 \ell(q)=0 & a_4^5,\ b_2^8,\ b_1^8,\ b_0^5,\ c_0^8,\ a^9,\ b^6
 \end{array}
\]
The leading pure powers make the quotient finite-dimensional,
and weighted homogeneity implies the only geometric affine
zero is the origin. The eliminated coefficients also vanish at
that point by their positively weighted equations. Any
geometric affine zero with \(D\ne0\) could be scaled to \(D=1\)
by the weight action after taking a sixth root. These outputs
therefore exclude all nonzero centered \(P\) in each special
fiber.

For the second row, \(15295\ell(q)\) is primitive and irreducible
modulo \(31\). The checker verifies degree four and coprimality
with \(q^{31^2}-q\). This rules out factors of degree at most
two, hence any proper factorization. It also verifies the
simple zero \(q=69\bmod 101\). Thus \(\ell(q)=0\) defines a
single characteristic-zero number-field parameter, and its
number field has a place realizing this second row.

To see the exact characteristic-zero implication, suppose a
phase over a characteristic-zero field kills the stated
prefix in one of these slices and has nonzero \(P\).
Specialize its finite-type coefficient algebra to a number
field, inverting one nonzero coefficient of \(P\). For the
quartic chart, the specialization contains the original
quartic parameter field. Extend the chosen residue-place
valuation to the specialized coefficient field. Write \(w_i>0\)
for the weights of the nonzero coefficients \(c_i\) of \(P\).
After a finite ramified extension choose \(s\) with
\(v(s)=-\min_i v(c_i)/w_i\). The weighted action gives
integral \(s^{w_i}c_i\), at least one a unit. The moments
still vanish, by exact homogeneity. Their coefficients before
reduction are integral at the chosen place. Indeed the
Laurent coefficient expansion divides only by powers of
two for the chart substitution, by \(3,5\) for the centered
core, and by integers at most \(5m+1\) for a moment of order
\(m\). The latter bound follows from \(i+h\leq 5m\) and
angular extraction \(i=0\). Both residue characteristics
exceed \(5m+1\) throughout the corresponding prefix.
Thus reduction produces a nonzero geometric \(P\) in the
computed special fiber, a contradiction. This is a
weighted-projective specialization proof over characteristic
zero. It does not infer a characteristic-zero identity from
a unit ideal on an affine modular chart alone.

### Replay

From the repository root, with python-flint 0.9.0 and msolve
0.6.5, run

~~~sh
python research/foundations/extended-geometry/compute_even_stripped_quintic_moments.py .scratch/even-bd --max-order 15
python research/foundations/extended-geometry/compute_even_stripped_quintic_moments.py .scratch/even-d0 --max-order 13 --d-zero
python research/foundations/extended-geometry/verify_even_stripped_quintic_projective_charts.py
~~~

The generator writes rational moments \(M_m\) and primitive
integer forms \(P_m\). It checks the first four against literal
Laurent expansion and verifies the exact pivot computations.
The projective checker rebuilds the two exact rational inputs,
checks the first four moments again by literal Laurent
expansion and verifies the weighted systems, coefficient
reductions and Gröbner leading powers. These finite
computations use their indicated prefixes. The general
parameter exclusions are given by the leading-form proof.
