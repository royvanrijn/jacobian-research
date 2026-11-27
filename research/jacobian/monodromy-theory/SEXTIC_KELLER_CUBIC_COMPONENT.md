# Sextic Keller maps with a cubic component

**Theorem SexticKellerCubicTargetComponent.** A three-variable
polynomial Keller map of ordinary total degree at most six over
a characteristic-zero field is invertible if a nonconstant linear
combination of its components has degree at most three.

Over \(\mathbb C\), take this combination as \(R\). The
[bounded completion theorem](LOW_DEGREE_COORDINATE_KELLER_MAPS_IN_THREE_VARIABLES.md)
handles degrees at most two and any component with a coordinate
completion of inverse degree at most five, since \(6\cdot5\le100\).
For a cubic \(R\), the
[fibre-only reduction to the exceptional scalar](QUINTIC_KELLER_CUBIC_COMPONENT.md#the-cubic-reduction)
applies unchanged. Specifically, the use of the other two components
in the proof of Wang and Yan's Proposition 3.16 is to make the
residue volume exact on smooth fibres, in Lemmas 3.13 and 3.14.
The remaining eliminations use the coefficients of \(R\)
and its submersion condition. Their alternatives are a completion
of inverse degree at most five and the scalar \(x+x^2v\).
The proof below excludes the exceptional scalar with mates of
degree at most six. The theorem descends to characteristic zero
by embedding the finitely generated coefficient field in
\(\mathbb C\) and using uniqueness of the formal inverse.
[Wang and Yan, Proposition 3.16 and Corollary 3.19](https://arxiv.org/pdf/2610.36495v1).

## The exceptional scalar and three residue restrictions

**Proposition SexticExceptionalCubicScalarExclusion.** There are no
polynomials \(P,Q\in\mathbb C[x,v,y]\) of ordinary degree at most
six with \(\det J(P,Q,x+x^2v)\) a nonzero constant.

Put \(r=x+x^2v\), \(B=\mathbb C[x,v]\), and
\(D=-x^2\partial_x+(1+2xv)\partial_v=-x^2\partial_x|_r\).
The bracket is
\[
               [P,Q]=D(P)Q_y-P_yD(Q)=\kappa\ne0
 \tag{1}
\]
up to the overall sign from the order of source coordinates.
The rational constants of \(D\) are \(\mathbb C(r)\), and
\(B[x^{-1}]=\mathbb C[r,x,x^{-1}]\).
There is no polynomial \(D\)-primitive of a nonzero constant,
as such a primitive would be \(c/x+H(r)\) with \(c\ne0\) and
\(H\in\mathbb C[r]\).

Use the three degree-independent restrictions in the
[exceptional cubic argument](QUINTIC_KELLER_CUBIC_COMPONENT.md#the-exceptional-scalar-through-degree-five).
Both \(y\)-degrees are positive, neither leading \(y\)-coefficient
is a nonzero constant, and two affine mates are impossible.
For further use, the chain-rule identity behind the formal residue
argument is
\[
      D(QP_y S)-\partial_y(QD(P)S)=-\kappa S,\qquad
                              S=P^{-1/m}.
 \tag{2}
\]
Here \(m=\deg_yP>0\), and \(S\) is the binomial series in \(y^{-1}\)
with leading term \(p^{-1/m}y^{-1}\) for \(p=[y^m]P\).
Work over the coefficient ring localized at \(p\) and with a
chosen \(m\)-th root of \(p\) adjoined. Taking the
\(y^{-1}\)-coefficient in (2) shows that \(p^{-1/m}\), times a
nonzero scalar, has a \(D\)-primitive there.
[Wang and Yan, Proposition 3.17](https://arxiv.org/pdf/2610.36495v1).

First suppose \(p=L(r)\) is invariant. The primitive given by
(2), multiplied by \(p^{1/m}\), lies in \(B[L(r)^{-1}]\), since
the binomial coefficients of \(p^{1/m}S\) lie in this ring.
It has \(D\)-derivative a nonzero constant and therefore has the
form \(c/x+H(r)\). If \(L(0)\ne0\), this expression is regular
along each component \(x=0\), \(1+xv=0\) of \(r=0\).
The latter component forces \(H(r)\) regular at \(r=0\),
and the former then cannot cancel \(c/x\). Thus every such
invariant leading coefficient must vanish at \(r=0\).

Next let \(\ell=\alpha x+\beta v+\gamma\) be affine linear
and nonconstant. Whenever (2) requires a primitive of
\(1/\sqrt{\ell}\), it follows that \(\beta=0\).
Indeed, over \(\overline{\mathbb C(r)}((x))\), the substitution
\(v=(r-x)/x^2\) gives \(\ell=\beta r/x^2+O(x^{-1})\)
if \(\beta\ne0\). The square root can be taken in this
Laurent field, and \(dx/(x^2\sqrt\ell)\) has a nonzero simple
pole coefficient. It cannot be a Laurent derivative.
If a primitive of \(1/\ell\) is required in
\(\mathbb C(r,x)\), again \(\beta=0\). When \(\beta\ne0\),
the differential is
\[
                  \frac{dx}{\alpha x^3+\gamma x^2-\beta x+\beta r}.
\]
Over the algebraic closure of \(\mathbb C(r)\) the denominator
has simple roots. Its derivative in \(x\) is nonzero, and any
common root with that derivative would make \(r\) algebraic
over \(\mathbb C\). A simple pole cannot be the derivative of
a rational function.

Finally, for any nonzero \(a\in B\) with no full divisor \(r-c\),
\[
                   B[a^{-1}]\cap\mathbb C(r)=\mathbb C[r].
 \tag{3}
\]
A pole at \(c\) of a reduced fraction in \(r\) would force
\(r-c\) to divide a power of \(a\), by the substituted Bézout
identity for the numerator and denominator.
Since \(r-c\) is squarefree by submersion, this forces
\(r-c\mid a\). In particular (3) holds for any univariate
nonzero \(a(x)\).

## Equal free-variable degrees

Let \(m=\deg_yP\ge n=\deg_yQ\), and let \(p,q\) be their
leading coefficients in \(y\). The degree bounds give
\(\deg p\le6-m,\ \deg q\le6-n,\ m\le5\), where \(m=6\)
would have a constant leading coefficient. The highest
bracket coefficient gives
\[
                             p^n/q^m\in\mathbb C(r).
 \tag{4}
\]
Because each full fibre factor \(r-c\) is squarefree of
degree three, a pole or zero at \(c\) in the reduced fraction
in (4) forces \(r-c\) to divide the denominator or numerator
coefficient, respectively. This also applies to \(p/q\)
when \(m=n\) by relative algebraic closedness of \(\mathbb C(r)\)
in \(\mathbb C(r,x)\).

When \(m=n\ge4\) the degree bounds disallow full fibre factors,
so \(p/q\) is constant. For \(m=n=3\), a coefficient with a
full fibre factor is a scalar multiple of \(r-c\), and a
nonconstant ratio forces both coefficients to be linear
polynomials in \(r\). A constant target combination then has
a nonzero constant leading coefficient, impossible.
For \(m=n=2\), each coefficient can have at most one full
fibre factor. Removing those factors in a nonconstant
ratio leaves proportional coefficients of degree at most
one. Constant target combinations therefore give leading
coefficients \(\ell,r\ell\) with \(\ell\) affine linear.
If \(\ell\) is constant it is forbidden. Otherwise (2) for
the mate with coefficient \(\ell\) requires a primitive of
\(1/\sqrt\ell\), so \(\ell\) is univariate by the restriction
above. Subtract \(r\) times the first mate from the second,
preserving (1). Its new \(y\)-degree must be one, since
zero \(y\)-degree is impossible. If \(q_1\) is its leading
coefficient, (4) now requires \(\ell/q_1^2\in\mathbb C(r)\).
The valuation in \(\mathbb C(r,x)\) at the root of the linear
polynomial \(\ell(x)\) gives an odd valuation on the left
and zero on the right. Thus a nonconstant ratio is impossible
also when \(m=n=2\).
The case \(m=n=1\) is the excluded affine pair. When the ratio
is constant, a constant target shear reduces one \(y\)-degree
and preserves the degree-six bound.

## An affine mate

The [affine-mate theorem for repeated-parallel scalars](KELLER_REPEATED_PARALLEL_SCALAR_AFFINE_MATE.md)
excludes \(n=1\) for \(r=x+x^2v\) without a degree bound.
It uses the two components of \(r=0\).

## Higher unequal degrees

Suppose \(m>n\ge2\). If neither leading coefficient
contains a full fibre factor, (4) is constant. Unique
factorization gives \(p=\lambda\ell^{m/g}\) and
\(q=\mu\ell^{n/g}\) for \(g=\gcd(m,n)\), nonzero scalars
\(\lambda,\mu\), and a nonconstant polynomial \(\ell\).
The inequalities
\[
               (m/g)\deg\ell\le6-m,\qquad
               (n/g)\deg\ell\le6-n
\]
leave only \((m,n)=(3,2),(4,2)\), with \(\ell\) affine
linear. For \((3,2)\), (2) requires a primitive of
\(1/\ell\), and for \((4,2)\) it requires one of
\(1/\sqrt\ell\). Thus \(\ell\) is univariate in \(x\).
In both cases \(q=q(x)\) has degree at most two and
\(Q=qy^2+By+C\).

These univariate cases are covered by the
[centered quadratic pole argument](QUINTIC_KELLER_BINARY_TRIPLE_LINE_SCALAR_EXCLUSION.md#the-cases-with-smaller-degree-two).
Its exact premises are (1), \(\deg_y P=3\) or \(4\),
\(q(x)\) nonconstant of degree at most two, and \(p\)
proportional to \(q^2\) in degree four. In that argument
the odd part of \(P\) after completing the square in \(Q\)
is polynomial over \(B\), because the commuting even part
has coefficients in \(\mathbb C[r]\) by (3).
The constant bracket forces the stationary term of \(Q\)
to have a pole at \(x=0\), and the odd part cannot be
polynomial for either possible order of \(q(x)\) at zero.

If \(p\) has a full fibre factor then \(m=3,\ n=2\) and
\(p\) is a scalar multiple of \(r-c\). Equation (4)
makes \(q\) invariant, with degree at most four, so it is
linear in \(r\). The residue restriction for invariant
leading coefficients makes \(p,q\) both nonzero scalar
multiples of \(r\). Complete the square in \(Q\) over
\(B[r^{-1}]=\mathbb C[r,r^{-1},x,x^{-1}]\), writing
\(Q=qY^2+K\). If \(pY^3+fY\) is the odd part of \(P\),
then comparison at \(Y^0,Y^2\) yields
\[
            -f DK=\kappa,\qquad 2qDf-3p DK=0.
\]
The even part contributes only odd bracket powers.
Consequently \(D(f^2)=-3(p/q)\kappa\) is a nonzero constant.
But \(f\) is a unit of the Laurent polynomial ring by the
first equality, hence \(f=\lambda r^i x^j\), with \(i,j\)
integers. Its square has \(D\)-derivative either zero or
of \(x\)-exponent \(2j+1\), a contradiction.

The remaining possibility is a full fibre factor in \(q\)
but none in \(p\). Write \(q=(r-c) w\) where
\(\deg w\le3-n\). After removing this factor the ratio
\(p^n/w^m\) is constant. The same degree inequalities with
\(6-n\) replaced by \(3-n\) leave only \((m,n)=(4,2)\),
\(p=\lambda\ell^2,\ w=\mu\ell\) with \(\ell\) linear.
The square-root restriction makes \(\ell\) univariate.
Reduce (1) on \(r=c\) with \(x\ne0\), whose coordinate
ring is \(\mathbb C[x,x^{-1}]\) for every \(c\).
The first mate still has degree four in \(y\), since
\(\ell(x)\ne0\) as a Laurent polynomial. If the second
mate lost its entire \(y\)-dependence, the reduced bracket
would force the first to be affine, impossible.
Thus the second mate is affine there with nonzero
coefficient \(b(x)\). The highest bracket equation gives
\(\ell^2/b(x)^4\in\mathbb C^*\), since the induced derivation
is \(-x^2\partial_x\). A nonconstant linear polynomial
\(\ell(x)\) is not a constant times a Laurent square,
by its odd highest exponent. This is the final contradiction.

The [exact checker](verify_sextic_cubic_component.py)
replays the two finite degree comparisons and the
centered bracket equations for the invariant case.
