# Weyl generation through Bernstein degree and constant-leading order fifty-nine

**Theorem PolynomialWeylBernsteinThroughFiftyNine.**
A characteristic-zero polynomial Weyl pair with nonzero scalar
commutator generates if a member has Bernstein degree at most \(59\),
or has nonzero constant differential leading coefficient and order
at most \(59\). Every member of a proper pair has order plus
leading-coefficient degree at least \(60\).

**Theorem PolynomialWeylBernsteinFiftySixPlanarSystems.**
If a proper characteristic-zero polynomial Weyl pair has a member
of Bernstein degree \(56\), an equivalent shaped pair has Bernstein
degrees \(56,140\) or \(56,196\). The two possibilities give the explicit
characteristic-zero planar Poisson systems (8) and (10) below.

## First crossings and scalar endpoints

The [bound through fifty-five](POLYNOMIAL_WEYL_BERNSTEIN_THROUGH_FIFTY_FIVE.md)
and the normalization used there give a subrectangular pair with
normalized corner \((g,d)\), nondividing coprime ratio \(2:b\)
with \(b\ge3\) odd, and \(g+d=28\). The maximal-\(X\)
face of either member is its single corner.

The local sum-twenty-five calculation in the
[bound through fifty-one](POLYNOMIAL_WEYL_CONDITIONAL_THROUGH_FIFTY_ONE_FROM_F7.md#degree-fifty-and-the-sum-twenty-five-crossing)
allows only ratio \(2:3\). By the degree-pair transfer and Moh's
bound cited below, this excludes the embedded sum-twenty-five
crossing at actual degrees \((56,84)\) and \((58,87)\) as well.
All lower normalized first-crossing sums are thereby excluded by
the earlier ratio \(2:b\) arguments except the potentially embedded sum
twenty-one. For that row \((u,v)=(6,15)\) has weight \((3,-1)\).
The bounds \(g\ge6\), \(4g-28\le3\) give \(g\le7\), and \(g=6\)
contradicts the pure maximal-\(X\) face. If \(g=7,d=21\), any
intervening normalized integral vertex would have first coordinate
strictly between six and seven. Here both coordinates decrease
strictly as the weight \((1,-s)\) increases in slope from zero up to
the first crossing. Thus the preceding edge has weight \((6,-1)\)
and common form \(X^6Y^{15}L(XY^6)\), with \(L\) linear and \(L(0)\ne0\).
The bracket companion of
[Guccione–Guccione–Valqui, Theorem 4.1(4)](sources/guccione-guccione-valqui-1111.6100v3.pdf)
is \(XYG(XY^6)\). Its equation forces \(21\deg G=14\)
at infinity, although \(G\) vanishes at the root of \(L\).
This excludes the embedded row.

By Proposition 7.3 of the same source, a sum-twenty-eight
first crossing has one of the following rows. Write the common
leading form as \(R=X^hY^jL(z)\) at \(z=X^{-\sigma}Y^\rho\),
with \(L(0)\ne0\), where \(C_1=(h,j)\). Put \(A=\rho h+\sigma j\),
and let \(\mu\) be the largest root multiplicity of
\(E^jL(E^\rho)\). The last column gives the associated scalar
successors before imposing the plane degree bound.
\[
\begin{array}{c|c|c|r|r|r|c|c}
 C_0&C_1&(\rho,\sigma)&A&\deg L&\deg G&\mu&
        (\mu,b,s_{\rm scalar})\\ \hline
 (7,21)&(3,1)&(5,-1)&14&4&1&4&\varnothing\\
 (7,21)&(2,1)&(4,-1)&7&5&2&3,4,5&(4,7,19/28),(5,3,7/12)\\
 (7,21)&(1,0)&(7,-2)&7&3&2&2,3&(2,3,16/21),(3,5,3/5)\\
 (8,20)&(3,0)&(4,-1)&12&5&1&5&\varnothing\\
 (10,18)&(1,0)&(2,-1)&2&9&4&3,\ldots,9&\varnothing .
\end{array}                                                     \tag{1}
\]
The bracket companion \(XYG(z)\) satisfies
\[
 AzG'+(h-j-(\rho+\sigma)zL'/L)G=c\ne0 .
\]
Every nonzero root of \(L\) is a simple zero of \(G\).
At multiplicity \(m\) one must have \(A\ne(\rho+\sigma)m\).
Comparison at infinity gives \(A\deg G=(\rho+\sigma)\deg L-(h-j)\).
In the row \(\mu=4,b=7\), let the quadruple and simple
roots be \(\alpha,\gamma\). Substituting \(G\) proportional
to \((z-\alpha)(z-\gamma)\) gives \(4\gamma=5\alpha\).

Use
[Guccione–Guccione–Valqui II, Proposition 4.3](sources/guccione-guccione-valqui-1205.6827v1.pdf)
and the lattice recurrences in the
[bound through fifty-five](POLYNOMIAL_WEYL_BERNSTEIN_THROUGH_FIFTY_FIVE.md#successor-tests).
At the initial crossing the nonmonomial edge forms have zero
leading bracket, positive weights with nondividing ratio \(2:b\),
and highest-\(Y\) endpoints of negative diagonal difference.
The Weyl commutator is a nonzero scalar. After shifting at a maximal
root, the normalized corner is \(C=((A-\sigma\mu)/\rho,\mu)\).
The actual support denominator is \(l=\rho\) and
\(\beta=l(\mu-C_x)=(\rho+\sigma)\mu-A\).
The initial slope is \(-\sigma/\rho\).
For scalar successors with \(N=2+b\), the reduced denominator
\(D\) of \(ls=l-\beta N/(N\mu-1)\) is coprime to \(\mu\).
Nonmonomiality of the smaller face requires \(D\le2\mu\),
so \(N\le2\beta\). A nonzero scalar bracket requires
\(D\mid2\) or \(D\mid b\), since otherwise both face polynomials
have zero constant term. The table imposes the strict slope bounds.

For a commuting multiple-root successor let \(M,D,j_0,s\)
be as in the same lattice recurrence, with
\(f(\eta)=\eta^{j_0}J(\eta^D)\) of degree \(\mu\).
A partition \((m_i)\) of \((\mu-j_0)/D\) yields
\(D\) nonzero roots of multiplicity \(m_i\) for each part.
The equation \(W'-(M/\mu)(f'/f)W=c\ne0\) requires at most \(M\)
distinct roots and excludes multiplicity \(\mu/M\).
If all root multiplicities times \(M/\mu\) are integers,
then \(f^{M/\mu}\) is polynomial up to scale and
\(W/f^{M/\mu}\) would be a rational primitive of
\(c/f^{M/\mu}\). This is impossible with multiple distinct roots
by the pole-degree argument in the
[bound through fifty-one](POLYNOMIAL_WEYL_CONDITIONAL_THROUGH_FIFTY_ONE_FROM_F7.md).
These tests leave only
\[
\begin{array}{c|c|c|c}
 (l,\mu,\beta)&(M,D,s)&\text{root multiplicities}&
             (l',\mu',\beta')\\ \hline
 (5,4,2)&(2,1,4/5)&3,1&(5,3,1)\\
 (2,6,4)&(9,4,5/8)&2,1,1,1,1&(2,2,1).
\end{array}                                                      \tag{2}
\]
In the first row a shift at the triple root keeps denominator five.
In the second row the double root is zero, so the identity shift
keeps denominator two. Neither output corner has a further
candidate. Proposition 4.3 applies at either input of (2) because
the commuting successors retain positive weight, the same
nondividing ratio, nonmonomiality, and negative highest-\(Y\)
diagonal differences. A commuting single-root shift requires
integral \(ls\), keeps \(C,l\) fixed and strictly increases the
previous slope. This only removes candidates from the tests above,
and can occur finitely many times.

Ratio \(2:3\) would give a planar Keller pair of degrees \(56,84\)
by the
[nondividing Bernstein transfer](POLYNOMIAL_WEYL_CONDITIONAL_THROUGH_FORTY_NINE_FROM_F7.md#two-inputs-from-polynomial-degree-pairs).
The [planar degree-one-hundred theorem](https://doi.org/10.1515/crll.1983.340.140)
excludes it. The remaining rows are
\[
\begin{array}{c|c|c|c|c}
 b&(\rho,\sigma)&R&\text{right slope}&
      \text{optional additional shifts of }Y\\ \hline
 5&(7,-2)&X(X^2Y^7-1)^3&3/5&
                    \lambda_3X^{-3/7},\ \lambda_4X^{-4/7}\\
 7&(4,-1)&X^2Y(XY^4-1)^4(XY^4-5/4)&19/28&
                    \lambda_2X^{-1/2}.
\end{array}                                                      \tag{3}
\]
Scalings make the selected root one in both common forms.
The first right shift is a nonzero multiple of \(X^{\sigma/\rho}\).
The optional single-root shifts have slopes strictly between the
first-crossing and scalar slopes.

## A cusp calculation for the left side

Write \(r=\rho+\sigma\), \(T=X^{1/\rho}\),
\(E=T^{-\sigma}Y\) at the first right weight. Then \(R=T^7f(E)\)
and \(\deg f=21\). In both rows \(f\) has several distinct roots,
and none has multiplicity \(7/r\).
In the Weyl algebra the leading-weight bracket is the Poisson
bracket with \(\{T,E\}=\rho^{-1}T^{1-r}\), up to sign.
In the Poisson algebra use the same weight on commuting variables.

In either algebra put \(K=7(b-2)+r\), equal to \(26\) or \(38\).
Start from the combination of \(Q^2,P^b\) canceling their top
terms, with \(Q^2\) monic in the combination.
Its bracket with \(P\) is a nonzero constant times \(Q\).
For every remaining largest weight \(j>K\), its top coefficient
\(A(E)\) commutes with the leading form of \(P\) and obeys
\(A'=(j/7)(f'/f)A\). Because the root multiplicities of \(f\)
have gcd prime to seven, \(j=7i\) is divisible by seven.
Every such \(i\) belongs to the semigroup generated by \(2,b\)
since \(b-1\le i<2b\). Subtract a multiple of a corresponding
product of \(P,Q\) at each step. Any such term involving \(Q\)
is linear in it and can be written as a symmetrized product
in the Weyl algebra. Its bracket with \(P\) has smaller weight
than that of \(Q\). Completing the square gives
\[
   S=(Q^\ast)^2-F(P),\qquad Q^\ast=Q-D(P)/2,\qquad
   \deg D\le(b-1)/2 .
\]
The bracket of \(P,S\) is a nonzero constant times \(Q^\ast\),
whose first-right-weight leading form is that of \(Q\).
Hence \(S\) has weight \(K\) exactly. If \(A(E)\) is its top
coefficient at that weight, its equation is
\[
 A'-(K/7)(f'/f)A=c f^{b-2},\qquad c\ne0 .
\]
At a root of \(f\) of multiplicity \(m\), the coefficient of the
lowest power shows that \(A\) vanishes to order \((b-2)m+1\),
since \(Km/7\) is nonintegral for the multiplicities in (3).
Thus \(A=f^{b-2}B\), with \(B\) polynomial vanishing at every
root. It satisfies \(B'-(r/7)(f'/f)B=c\). Comparison at
infinity gives \(\deg B=3r\), since its degree exceeds one.
Consequently the highest-\(Y\) endpoint of \(S\) at the
first right weight is \((K,3K)\). The corresponding
endpoints of \(P,Q^\ast\) are \((14,42),(7b,21b)\).

Start the weight sweep immediately to the left of this weight,
toward \((0,1)\) and then through \((-\theta,1)\).
Initially these three monomials are the unique maximizing points.
The weights of \(P,Q^\ast,S\) retain the ratio \(14:7b:K\)
where the leading brackets of \(P,Q^\ast\) and \(P,S\)
must vanish. To see persistence, write a weight \((u,v)\).
The support corner of \(P\) gives \(w(P)\ge14(u+3v)>0\).
The strict inequality \(r(u+3v)>u+v\) holds after the initial
right weight and up to \(\theta<3\). In the first quadrant
\(r>1\), and after \((0,1)\) one has \(r(3-\theta)>1-\theta\).
The ratio \(7b/14\) then forces the leading \(P,Q^\ast\)
bracket to vanish, and the margin
\[
   (K/14-(7b/14-1))w(P)=(r/14)w(P)>u+v
\]
forces the leading \(P,S\) bracket to vanish as well.
At generic weights vanishing brackets make the maximizing
exponent vectors proportional. The ratios are continuous
through each edge and locally constant where the strict
inequalities hold, so they cannot leave the asserted ratio.
This applies in the Weyl algebra up to \((-1,1)\), with
positive contraction loss throughout the open interval.
In the Poisson plane the argument works up to \(\theta<3\).
Also \(Q,Q^\ast\) have equal leading forms on this sweep by
their positive weights and \(\deg D<b/2\).

At an edge in this interval, weighted Euler differentiation
and unique factorization make the three leading forms
proportional to powers of one polynomial (Laurent polynomial
after shearing). The powers are \(14,7b,K\) because these
three integers are coprime. A first left edge at
\((-\theta,1)\) would therefore have a common form
joining \((1,3)\) to \((0,3-\theta)\).
There is no such edge for \(0<\theta<1\).
By continuity the Weyl Euler maxima are \(28,14b\).

Now transfer to the Poisson plane by the
[split-prime argument for degree forty-two](POLYNOMIAL_WEYL_FORTY_TWO_SPLIT_PRIME_F7_TRANSFER.md).
The supports, incoming forms (3), Weyl Euler maxima and actual
right caps are finite coefficient equations with nonzero
conditions. The same construction over a number field, with
the two lowest Euler coefficient polynomials split modulo
the chosen rational primes, makes the central Euler forms
pure corners \(x^{14}y^{42},x^{7b}y^{21b}\).
The actual right caps also transfer. In characteristic \(p\)
prime to \(\rho\), their shift \(Y\mapsto Y+g(X^{1/\rho})\)
induces \(Y^p\mapsto Y^p+g(X^{1/\rho})^p\). The
\((p-1)\)-st \(X\)-derivative of each power in \(g\) vanishes
for sufficiently large \(p\), as in the degree-forty-two proof.
The scalar right filtrations have positive contraction loss,
so taking \(p\)-th powers preserves the caps with weights
divided by \(p\) on the center. The corresponding finite
coefficient equations for the center pair lift to
characteristic zero by the Nullstellensatz and clearing
denominators. This uses Lauritzen–Thomsen,
Propositions 1.11–1.12 and Theorem 1.7(iii)(a),
in the [primary source](https://arxiv.org/pdf/1508.02499v1).

Repeat the cusp calculation with the transferred Poisson pair.
Its pure Euler corner eliminates a first edge at \(\theta=1\).
The only possible first left edge for \(1<\theta<3\) is
\(\theta=2\), with common factor \(xy^3-\gamma y\).
The Laurent shear \(x\mapsto x+\gamma y^{-2}\) takes this
edge to the corner. It preserves the corner terms of \(P,Q,S\)
since it decreases the \(x\)-exponent and their maximal
\(x\)-exponents, obtained at \((1,0)\), are \(14,7b,K\).
The sweep after the shear starts strictly before \(\theta=2\),
where the shear lowers weights. It continues across \(\theta=2\)
with the same ratios. A new first edge for \(2<\theta<3\) is
impossible by integrality. If there is no \(\theta=2\) edge
use \(\gamma=0\).
Letting \(\theta\) tend to three gives \(j\le3i\) on the
supports of the sheared \(P,Q\). In the original variables
they therefore lie in the polynomial subring generated by
\[
                       v=y^{-1},\qquad z=xy^3-\gamma y .
                                                               \tag{4}
\]
The Jacobian in \((z,v)\) is a nonzero scalar times \(v\).

## The two explicit systems

For \(b=5\), put \(T=x^{1/7}\). The cap in the coordinates
\(y=\lambda T^{-2}+\lambda_3T^{-3}+\lambda_4T^{-4}+\eta\),
\(\lambda^7=1\), uses weights \(W(T)=5,W(\eta)=-21\)
and levels \(4,10\). Put
\[
  s=z^2-v-\kappa_1z-\kappa_2,\quad
  \kappa_1=7\lambda^2\lambda_3,\quad
  \kappa_2=7\lambda^4(\lambda\lambda_4-\lambda_3^2).
                                                               \tag{5}
\]
The maximal-weight terms of \(z,s\) are
\(\lambda^3T,7\lambda^5 T^4\eta\), respectively.
Indeed \(z^2-v=x^2y^6-y^{-1}-2\gamma xy^4+\gamma^2y^2\).
Writing \(y=\lambda T^{-2}(1+\varepsilon)\), the first
two terms give \(\lambda^6T^2((1+\varepsilon)^6-
(1+\varepsilon)^{-1})\), whose bracketed factor has
coefficients \(7,14\) in degrees one and two. Subtraction
of \(\kappa_1 z+\kappa_2\) cancels the weights five and zero.
The \(\gamma\)-terms have weights below \(-1\).
Distinct leading monomials for \(z^i s^j\) give the caps
\[
    i+2j\le14,35,\qquad 5i-j\le4,10                                  \tag{6}
\]
for \(P,Q\), respectively. The first caps follow from
the incoming weight \((7,-2)\), where \(z,v\) have
weights \(1,2\) and their prescribed forms are proportional
to \(z^2(z^2-v)^6,z^5(z^2-v)^{15}\).
In \((z,s)\) those become proportional to \(z^2s^6,z^5s^{15}\).
The coefficient \(p_{1,1}\) is nonzero, since otherwise the
bracket of the right leading parts cannot have a \(z^2\) term.

Translate \(z\) to remove the linear term on the right, then
use polynomial scalings and \((z,s)\mapsto(\mu z,\mu^2s)\)
as in [Niedbala Giraudin, Lemma 6.1](https://hexagonmath.org/2611.00717v1).
Subtract constants from \(P,Q\), and multiples of \(P^2,P\)
from \(Q\), to obtain
\[
 p_{1,1}=p_{2,6}=q_{5,15}=1,\quad
 p_{0,0}=q_{0,0}=q_{4,12}=q_{2,6}=0,\quad
 (q_{2,0},q_{3,5},q_{4,10})=(5/16,15/8,5/2).
                                                               \tag{7}
\]
The subtractions preserve the full prescribed highest
\((1,2)\)-homogeneous forms. The normalized equation is
\[
                   \{P,Q\}_{z,s}=-\tfrac58(z^2-\kappa-s).
                                                               \tag{8}
\]

For \(b=7\), put \(T=x^{1/4}\) and
\(y=\lambda T^{-1}+\lambda_2T^{-2}+\eta\), with \(\lambda^4=1\).
The right weights are \(W(T)=7,W(\eta)=-19\), and the
levels for \(P,Q\) are \(2,7\).
Set \(s=z-v-4\lambda^2\lambda_2\). The leading terms of \(v,s\)
are \(\lambda^{-1}T,4\lambda^2T^2\eta\).
Indeed \(z-v=v(xy^4-1)-\gamma y\) and the constant term
is \(4\lambda^2\lambda_2\). The quadratic terms in
\(\lambda_2T^{-1}+\eta T\) and the \(\gamma\)-term have
smaller weight. Distinct leading monomials give the caps
\[
          i+j\le14,49,\qquad 7j-5i\le2,7                         \tag{9}
\]
on monomials \(s^i v^j\) of \(P,Q\), respectively.
The first caps follow from the incoming weight \((4,-1)\),
where \(z,v\) have weights \(1,1\). The prescribed total-degree
tops are proportional to
\(s^8(s+v)^4(s-v/4)^2\) and \(s^{28}(s+v)^{14}(s-v/4)^7\).
Normalize the scales so the coefficients at \((8,6)\),
\((28,21)\) are one, giving factors \(16,-16384\).

With \(w=s^7v^5\), the right leading forms are
\(sv(A+w)\) and \(vG(w)\), with \(G\) monic of degree four.
Their bracket is
\[
        v\bigl((A+8w)G-2w(A+w)G'\bigr).
\]
The nonzero constant factor on the right forces \(A\ne0\).
A common scaling of \(s,v\) and rescaling \(P,Q\)
to retain their total-degree tops sets \(A=1\).
Then the remaining coefficients are determined by
\[
        G=w^4+\tfrac72w^3+\tfrac{35}8w^2+\tfrac{35}{16}w
              +\tfrac{35}{128}.
\]
Subtract constants from both polynomials and powers
\(P^3,P^2,P\) from \(Q\) to set \(q_{3,3},q_{2,2},q_{1,1}\)
to zero in that order. The term \((sv(1+w))^k\) has
coefficient one at \(s^k v^k\), and these subtractions are
strictly below both degree caps of \(Q\).
The normalized equation is
\[
                         \{P,Q\}_{s,v}=\tfrac{35}{128}v .
                                                               \tag{10}
\]

## Exact exclusion of the planar systems

For (8) assign layer \(k\) to monomials of weight \(4-k\)
in \(P\) and \(10-k\) in \(Q\). For (10) use \(2-k\) and
\(7-k\), respectively. In either case the bracket at layer \(k\)
uses a fixed rational linear map on the unknown coefficients
in that layer, plus products of coefficients already determined
at earlier layers. This follows from
\[
 \{x^i y^j,x^{i'}y^{j'}\}
       =(ij'-ji')x^{i+i'-1}y^{j+j'-1}.
\]
Specified terms from the full incoming leading forms appear
at the appropriate layers as fixed values.

The [exact replay](verify_weyl_bernstein_fifty_six.py) assembles
those maps, eliminates all pivot coefficients, and retains the
necessary equations from rows with no pivot. For (8) the layer-eleven
equation is a nonzero rational constant after substitution of
the earlier layers. This proves exclusion directly, using the
stronger incoming form conditions of the transferred pair.

For (10) the layers through 56 leave five free coefficients.
Name them \(A,B,C,D,E\) in the order of the \(Q\)-coefficients
with exponent pairs \((6,5),(5,4),(4,3),(7,5),(3,2)\).
There is one compatibility polynomial at each of layers
\(24,29,34,36,39,41,44,46,48,49,51,53,54,56\).
Give the five variables weights \(2,4,6,7,8\), respectively.
Every monomial of the equation at layer \(k\) has
weight at most \(k\) and congruent to \(k\) modulo twelve.
The replay constructs the 1341 products of these equations
with monomials of weight at most \(72-k\) and congruent to
\(72-k\) modulo twelve. It selects a square system of rank
1258 using modular linear algebra, solves the selected
system over the rationals, then expands the resulting
identity \(1=\sum_k h_k f_k\) over the rationals. The identity
is verified in the full polynomial ring with every coefficient
included. The modular computation is used to select pivots.
This excludes (10) and proves the Bernstein bound through \(56\).

## Degrees fifty-seven through fifty-nine

For degree \(57\) the same shaped reduction gives
\(57=a(g+d)\), \(g+d\ge16\), \(a\ge2\). Thus \(a=3,g+d=19\).
The only possible normalized crossing sums are sixteen and
eighteen. Smaller sums fail the earlier numerical first-crossing
tests, and seventeen and nineteen are prime and fail the
gcd condition. For \(C_0=(u,v)\) of weight \((\rho,\sigma)\)
the outer corner satisfies
\[
 (\rho+|\sigma|)(g-u)\le|\sigma|(g+d-u-v),\qquad g\ge u .          \tag{11}
\]
At sum sixteen the two rows have \(C_0=(4,12),(6,10)\)
with weights \((4,-1),(2,-1)\), respectively.
The first gives \(g=u\). The second gives \(g=u\) or \(g=u+1\).
At sum eighteen the unique row has \(C_0=(6,12)\), weight
\((3,-1)\), giving \(g=u\). In every case \(g=u\) contradicts
the pure maximal-\(X\) face because \(d>v\).
The alternative \(g=u+1\) at sum sixteen is equality in (11)
and the outer corner has larger \(Y\)-exponent. This contradicts
the defining highest-\(Y\) endpoint of the crossing.
Hence degree \(57\) is excluded.

For degree \(58\) the shaped reduction gives \(a=2,g+d=29\).
The sum-twenty-nine table is empty by the gcd test. At embedded
sums twenty-seven and twenty-eight, each first row above
or in the bound through fifty-five has
\(|\sigma|(29-u-v)<\rho+|\sigma|\). By (11) it has \(g=u\),
contradicting the pure maximal-\(X\) face.
For the embedded sum twenty-one, \(g\ge6\) and \(4g-29\le3\).
The case \(g=6\) contradicts the pure maximal-\(X\) face.
The case \(g=8,d=21\) is equality in the weight bound
and contradicts the highest-\(Y\) endpoint \((6,15)\).
For \(g=7,d=22\) the edge preceding \((6,15)\) has weight
\((7,-1)\), by the integral normalized vertex argument
for degree \(56\). The common form is \(X^6Y^{15}L(XY^7)\)
with \(L\) linear and \(L(0)\ne0\).
The bracket companion now requires \(27\deg G=15\),
again impossible. The other lower crossings are excluded
as at degree \(56\).

At degree \(59\) Bernstein mate descent would give a
nondividing pair of smaller degree \(59\), since all degrees
through \(58\) are excluded. This has degree-gcd a proper divisor
of the prime \(59\), contrary to the gcd lower bound sixteen.

The constant-leading weighted descent in the
[bound through fifty-three](POLYNOMIAL_WEYL_BERNSTEIN_THROUGH_FIFTY_THREE.md)
also gives the assertion for starting orders \(56,57,58,59\).
Use induction at each new degree and order together. The smaller
order in the dividing and binomial-change branches is at most
\(\lfloor59/2\rfloor=29\), already excluded. The intermediate
change either gives smaller constant-leading order or a nonconstant
leading coefficient with order plus leading-degree at most the
starting order. The latter contradicts the newly available Bernstein
bound by the infinity shear and dividing-mate induction. The
weight-one branch gives a member of Bernstein degree equal to the
starting order. This proves both assertions and the individual
leading-degree bound. As in the earlier proofs, use the embedding
of a finitely generated characteristic-zero coefficient field
into \(\mathbb C\) and descent of generation by faithful flatness.
