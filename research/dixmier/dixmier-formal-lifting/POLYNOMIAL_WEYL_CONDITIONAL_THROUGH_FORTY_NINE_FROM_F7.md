# Weyl generation through Bernstein degree and constant-leading order forty-nine

**Theorem PolynomialWeylConditionalThroughFortyNineFromF7.**
A characteristic-zero polynomial Weyl pair with nonzero scalar
commutator generates whenever one member has Bernstein degree at most
\(49\), or has nonzero constant differential leading coefficient and
order at most \(49\).

## Two inputs from polynomial degree pairs

The
[direct degree-forty-two exclusion](POLYNOMIAL_WEYL_FORTY_TWO_DIRECT_F7_CERTIFICATE.md)
derives the exact \(F_7\) certificate system from the Weyl
residual and excludes a proper pair with a member of Bernstein
degree \(42\). The induction below uses only this degree exclusion,
so it also follows from the
[alternative via the planar \(F_7\) theorem](POLYNOMIAL_WEYL_FORTY_TWO_SPLIT_PRIME_F7_TRANSFER.md).
We also use the following consequence of the reduction argument.

**Proposition WeylBernsteinDegreePairTransfersToPlanarKeller.**
If the Bernstein degrees of a characteristic-zero Weyl pair with
nonzero scalar commutator are \(i,j>1\), and neither divides the
other, then there is a plane Keller counterexample over
\(\mathbb C\) of degrees \(i,j\).
Indeed, scale the commutator to one and specialize to algebraic
coefficients using a finite system for the commutator and the two
exact degrees. To encode exact degrees, retain one chosen nonzero
top coefficient for each member. Over each sufficiently large
good prime the center map of the resulting Weyl endomorphism
has nonzero constant Jacobian by
[Lauritzen–Thomsen, Theorem 1.7 and Propositions 1.11–1.12](https://arxiv.org/pdf/1508.02499v1).
In its central variables the degrees are still \(i,j\), since the
Bernstein leading forms of the \(p\)-th powers are the \(p\)-th
powers of the two leading forms. The finite equations for a
commutative Keller pair with these exact degrees then have a
complex solution by the same denominator-clearing argument as in
the split-prime transfer. Such a solution cannot be an
automorphism, by degree divisibility for plane polynomial
automorphisms. This last fact follows from the
Jung–van der Kulk normal form, and is the degree check used
for [Guccione–Guccione–Valqui, Remark 4.2](https://arxiv.org/pdf/1401.1784v3).
In particular \((48,64)\) cannot be a Weyl degree pair:
the corresponding planar pair would contradict
[Moh’s published degree-100 theorem](../../BREAKTHROUGH_TARGET_ASSESSMENT.md#four-variable-hessian-conjecture-and-planar-degree-bounds).

## Bernstein induction

Start with the
[bound through forty-one](POLYNOMIAL_WEYL_BERNSTEIN_DEGREE_THROUGH_FORTY_ONE.md)
and the exclusion of degree \(42\) above.
For each odd degree \(43,45,47,49\) after the preceding even
degree is excluded, the standard Bernstein mate descent gives
a smaller nondividing degree pair with smaller degree at most
that odd degree. No smaller degree at most the preceding even
degree is possible by induction. The odd degree itself has no
proper divisor at least \(16\), contradicting the general
degree-gcd lower bound of Guccione–Guccione–Valqui,
Corollary 7.4, in their
[Weyl shape theorem](sources/guccione-guccione-valqui-1111.6100v3.pdf).
This is the same descent used in the degree-forty-one bound.

For even steps \(44,46,48\), use the quadratic bracket-mate
orientation, leading-degree saturation, and infinity shear in
the degree-forty-one proof. Apply the induction to Bernstein
degree and constant-leading order together, as described below.
In particular the preceding Bernstein bound implies that a member
of a proper pair with nonconstant leading differential coefficient
has order plus leading-coefficient degree at least the current
even degree, by the infinity shear and dividing-mate argument.
After fixed-member differential mate descent, the degrees of the
resulting shaped pair are \(a(g+d),b(g+d)\), with
\[
            a(g+d)=n,\qquad g+d\ge16,\qquad a,b\ge2,
                       \qquad \gcd(a,b)=1 .
\]
Here \(n\) is the current even degree and \((g,d)\) is the
normalized outer corner. The positive branch of the proof
uses the axis with smaller top exponent for differentiation
and the actual infinity shear, so the maximal \(X\)-face
after Fourier orientation is its single corner. The
possibility of a residual dividing smaller mate contradicts
the fixed-member leading-degree bound, as in the earlier
degree-forty-one and forty-two arguments.

For \(n=44,46\) one has \(a=2\).
For \(n=48\) either \(a=2\), or \(a=3\) and \(g+d=16\).
In the latter case the exact shaped gcd-sixteen quotient
restriction from the
[Bernstein degree thirty-three proof](POLYNOMIAL_WEYL_BERNSTEIN_DEGREE_THROUGH_THIRTY_THREE.md)
gives \(b=4\) if \(a=3\). Thus the actual Bernstein degrees
are \(48,64\), excluded by the commutative reduction above,
without its previous field square-class test.

It remains to use \(a=2\) and \(g+d=22,23,24\).
The first-crossing exclusions for normalized sums through
twenty from the degree-forty-one proof apply to the ratio
\(2:b\), irrespective of the outer sum.
At sum twenty-one the only crossing row is
\[
                 C_0=(6,15),\qquad(\rho,\sigma)=(3,-1),
\]
by the
[degree-forty-two crossing calculation](POLYNOMIAL_WEYL_BERNSTEIN_FORTY_TWO_REAL_BOUND.md).
It cannot lie inside any of these larger rectangles.
Indeed write \(s=g+d\in\{22,23,24\}\). Since \(C_0\) is in
the rectangle one has \(g\ge6\), and since the outer
corner is in the actual support,
\(3g-d=4g-s\le 3\), the weight of \(C_0\).
This forces \(g=6\). But \(C_0\) then has maximal
\(X\)-exponent and \(Y\)-exponent \(15<d\),
contradicting the pure maximal \(X\)-face.
Thus no extension of the \(F_7\) arithmetic condition to
an embedded crossing is needed.

Sum twenty-two is excluded by the
[sum-twenty-two successor analysis](POLYNOMIAL_WEYL_BERNSTEIN_FORTY_TWO_REAL_BOUND.md#the-nonsquare-field-bound-through-forty-seven),
and sum twenty-three is prime and fails the first-corner
gcd condition of Guccione–Guccione–Valqui,
Proposition 7.3(3). For sum twenty-four use the
[sum-twenty-four successor analysis](POLYNOMIAL_WEYL_TWO_NONSQUARES_THROUGH_FORTY_NINE.md#the-two-new-normalized-corners-of-sum-twenty-four).
The two successor exclusions at sums twenty-two and
twenty-four impose no field square-class condition:
their tests are the integer crossing and support
lattices, the degree and multiplicity restrictions from
the homogeneous bracket-mate equation, and the rational
residue exclusions at the finitely retained corners.
The nonsquare assumptions in those two sources are used
elsewhere in their inductive steps, and those uses are
replaced here by exclusion of degree forty-two, the
pure-maximal-face argument and the planar exclusion of
\((48,64)\), respectively.
This closes the even Bernstein steps.

## Constant-leading order

After each Bernstein step \(n,n+1\), use the
constant-leading weighted descent in the
[proof through forty-one](POLYNOMIAL_WEYL_BERNSTEIN_DEGREE_THROUGH_FORTY_ONE.md#closing-the-degree-induction)
and the continuation through the next orders in the two
sources above. For orders \(42\) through \(49\), the
smaller-order case after a binomial coordinate change
or a residual dividing mate has order at most \(24\)
and is already excluded. In the intermediate change
the same small-order bound applies if the new leading
coefficient is constant. Otherwise the resulting order
plus leading-coefficient degree is at most the starting
order, contrary to the updated lower bound from the
corresponding Bernstein step by infinity shear and
dividing-mate descent. The weight-one branch has
Bernstein degree equal to the starting order after
depression. This proves the two bounds jointly at each
stage, through \(49\).

The arguments can be made over \(\mathbb C\), as in the
stated \(F_7\) transfer. For any characteristic-zero
coefficient field, the finitely many coefficients of a
given pair lie in a finitely generated field that
embeds into \(\mathbb C\). Generation descends along
these field extensions by faithful flatness on vector
spaces, so the same conditional cutoff holds over the
original field.
