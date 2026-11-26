# Skew-adjoint prime orders and polynomial orders four and six

Let \(K\) be a characteristic-zero one-variable algebraic function field with nonzero derivation
\(\delta\), constant field \(C\), and Ore ring \(K[D;\delta]\). The adjoint
fixes \(K\), sends \(D\) to \(-D\), and reverses products. A scalar mate has
commutator \(\kappa\in C^\times\).

## Prime order

**Theorem SkewAdjointPrimeOrderObstruction.** A skew-adjoint element of odd
prime order in \(K[D;\delta]\) has no scalar mate, including a mate of
order one.

Let \(L^*=-L\), with prime order \(p\), and suppose \([M,L]=\kappa\).
[PrimeOrderMateNormalForm](ALL_COPRIME_ALGEBRAIC_STRING_OBSTRUCTION.md#5-constant-field-descent-and-the-prime-order-conclusion)
gives \(M=A(L)+bD+d\) over the original field, with \(A\in C[T]\) and
\(b,d\in K\). At each descent step of order \(jp\), the leading
commutator coefficient gives \(\delta(\operatorname{lc}(M)/\operatorname{lc}(L)^j)=0\).
Thus subtraction uses an original-field constant. The coprime order
exclusion applies to a residual order at least two and coprime to \(p\).
If the residual is of order zero, the coefficient of \(D^{p-1}\) in its
commutator forces its derivative to vanish. The commutator then vanishes.

Taking the adjoint gives \([M^*,L]=\kappa\). The average has the form
\[
 \frac{M+M^*}{2}=\frac{A(L)+A(-L)}2+d-\frac{\delta(b)}2.
\]
The polynomial commutes with \(L\). A scalar \(h\in K\) cannot give the
required commutator since the \(D^{p-1}\) coefficient in \([h,L]\) is
\(-p\operatorname{lc}(L)\delta(h)\). Its vanishing forces \([h,L]=0\).

**Corollary SelfAdjointSixSkewMateFloor.** If \(L^*=L\) has order six and
\([M,L]\in C^\times\), the skew part \((M-M^*)/2\) is a scalar mate of
order one or odd composite order at least nine.

Indeed a nonzero skew element has odd order. Order three is excluded
by the theorem and order five or seven by the coprime order exclusion.

**Corollary SelfAdjointFourSixExclusion.** A scalar string of orders four
and six cannot have either element self-adjoint.

For a self-adjoint element of order four use
[PowerTwoNormalizedMateDescent](SELFADJOINT_POWER_OF_TWO_STRING_CLASSIFICATION.md).
Its coprime reduction and centralizer argument force mate order one or
a multiple of four. For a self-adjoint element \(L\) of order six, the
skew part of a proposed order-four mate has order at most three, hence
order one. Over a finite coefficient-field extension choose
\(\rho^6=\operatorname{lc}(L)\) with derivation extended algebraically.
Use the skew first-order generator \(E=\rho D+\delta(\rho)/2\). In the
new coordinate \(L\) is monic. The order-six coefficient of the
commutator with a first-order mate forces its leading coefficient in
\(E\) to be constant. A skew first-order mate is therefore a nonzero
constant multiple of \(E\). The scalar equation gives \(L=f(E^2)+w\)
for a constant-coefficient polynomial \(f\) and a coefficient \(w\)
satisfying \([E,w]\in C_{\rm ext}^\times\), where \(C_{\rm ext}\) denotes
the actual extended constant field. The centralizer calculation (A5)–(A6)
of the same note applies to this displayed form in any even order.
Fractional-power descent places the centralizer in \(C_{\rm ext}[w][E]\)
and the Weyl coordinate shift \(w\mapsto w+f(E^2)\) gives centralizer
\(C_{\rm ext}[L]\). The original mate minus the skew mate would be a
commuting element of order four, which is impossible.

**Corollary SelfAdjointSixRegularLongMateFloor.** In the exact-first-image
regular scalar Long chart with classical orders two and three, a
self-adjoint generic order-six string member requires a generic mate
of order at least nine.

If the skew part had order one, the centralizer calculation just given
would place the original mate in the form \(A(L)+M_1\).
Here \(M_1=(M-M^*)/2\) is in the original field, so every positive
coefficient of \(A\) descends by the successive ratios of leading
coefficients, with its constant term absorbed into the remainder.
The valuation argument for
[PrimeOrderSpecializationObstruction](ALL_COPRIME_ALGEBRAIC_STRING_OBSTRUCTION.md#6-a-specialization-obstruction-without-a-unit-leading-coefficient)
now applies to this displayed normal form. The argument from the
minimum coefficient valuation uses only the normal form and a constant
uniformizer. It would force the classical order of \(L\) to divide the
classical order of \(M\), incompatible with the unordered pair
\(\{2,3\}\). The skew part therefore has order at least nine, and taking
the adjoint does not increase the order of \(M\).

The original-field first-order descent uses the constant lower-flow
step in the coprime proof. Its analytic input is Beals–Sattinger,
[arXiv solv-int/9801010](https://arxiv.org/abs/solv-int/9801010), §9,
Theorem 9.1, Lemma 9.3 and Theorem 9.4, with the lower-flow remark on
p. 43 of the [retained source](sources/beals-sattinger-solv-int-9801010v1.pdf).
The adjoint arguments introduce no further analytic input.

## Monic polynomial orders four and six

Let \(k\) be a characteristic-zero field and \(D x-xD=1\). A polynomial
operator of order four with constant nonzero leading coefficient may
be divided by this coefficient and depressed by the polynomial
automorphism \(D\mapsto D+f(x)\), \(x\mapsto x\). Its canonical form is
\[
 P=(D^2+u/2)^2+vD+v'/2+z,\qquad u,v,z\in k[x].
 \tag{1}
\]
The term \(vD+v'/2\) is its skew part. Polynomial admission of (1)
requires an actual polynomial pair and a constant leading coefficient
before depression.

**Theorem PolynomialFourSixConstantSkewExclusion.** In (1), if \(v\) is
constant, then \(P\) has no order-six polynomial scalar mate.

**Proposition PolynomialFourSixCoefficientSystem.** For any \(u,v,z\)
in (1) write \(W=P^{1/4}=D+\sum_{j\ge1}w_jD^{-j}\) in the formal Ore
ring and \(B_i=(W^i)_+\) for its nonnegative-order part. After multiplication
by a nonzero constant and addition of a constant-coefficient multiple
of \(P\) and of \(1\), every order-six polynomial scalar mate has the
form
\[
 Q=B_6+aB_5+bB_3+cB_2+eD,\qquad a,b,c,e\in k.                 \tag{2}
\]
The exact scalar-mate equations are \(r_2=r_1=0\), \(r_0\in k^\times\),
where for \(P=\sum_i p_i D^i,\ Q=\sum_j q_j D^j\),
\[
 r_t=\sum_{i+j-h=t}
   \left(\binom ih p_i q_j^{(h)}-\binom jh q_j p_i^{(h)}\right).
 \tag{3}
\]
All \(r_t\) of order greater than two vanish identically.

To see exhaustion, \([P,(W^i)_+]=-[P,(W^i)_-]\) has order at most two.
A residual of positive order \(j\) whose commutator with \(P\) has order
at most two has constant leading coefficient since the coefficient
in order \(j+3>2\) is four times its derivative. Subtract that multiple
of \(B_j\) and continue, including the order-zero residual whose
derivative is detected in order three. This proves (2). To evaluate
(3), \(4w_j\) occurs in order \(3-j\) of \(W^4\) after the earlier
terms are known. The first five \(w_j\) determine all \(B_i\) in (2).

For the exclusion, \(v=0\) is the self-adjoint case. For \(v\in k^\times\),
the coefficient \(r_2\) equals \(H'\) for
\[
\begin{aligned}
H={}& (u''/2-3z)v-2cv-eu\\
 &+a(-5u^3/32-5uu''/16-5(u')^2/32-5uz/4-5v^2/8-u^{(4)}/16-5z''/4)\\
 &+b(-3u^2/8-u''/4-3z).
\end{aligned}                                                    \tag{4}
\]
Also
\[
\begin{aligned}
r_1-r_2'={}&u z'''+3u'z''/2+u''z'/2+z^{(5)}/2-3zz'
 -(5av/4+2c)z'\\
 &+(5av/16)(uu'+u''')+3v(v+b)u'/4.
\end{aligned}                                                    \tag{5}
\]
These are Ore product identities from (3). The
[exact coefficient checker](verify_order_46_constant_skew.py)
replays them.

If \(u\) is constant, a nonconstant \(z\) gives a unique highest term
from \(zz'\) in (5). Both constant would give a zero commutator.
Put \(n=\deg u\ge1\). If \(a\ne0\), (4) forces \(\deg z=2n\) by the
unique highest terms from \(u^3\) and \(uz\) outside that degree.
Then \(zz'\) gives the unique highest term in (5), a contradiction.
Thus \(a=0\). If \(b=-v\), the highest term in (4) is from \(3vu^2/8\).
If \(b\ne0,-v\), (4) again forces \(\deg z=2n\), contradicted by
(5). Consequently \(b=0\). Equation (4) now gives
\(z=u''/6-eu/(3v)+z_0\). If \(e\ne0\), the unique highest term in
(5) is again from \(zz'\), so \(e=0\).

With \(z=u''/6+z_0\), the remaining identities are
\[
\begin{aligned}
r_1&=u u^{(5)}/6+u'u^{(4)}/4+u^{(7)}/12
 -(c/3+z_0/2)u'''+3v^2u'/4,\\
r_0-r_1'/2&=v u^{(5)}/12+\frac v2(c+3z_0/2)u'.
\end{aligned}                                                    \tag{6}
\]
For \(n>4\) the highest coefficient in \(r_1\) is the square of the
leading coefficient of \(u\) times
\(n(n-1)(n-2)(n-3)(5n-8)/12\). For \(1\le n<4\) the highest term comes
from \(3v^2u'/4\). For \(n=4\), equation \(r_1=0\) first gives
\(u^{(4)}=-3v^2\) and then \(c=-3z_0/2\). Thus \(r_0=0\) in the
only remaining case, proving the result.

## Remaining admissions

A monic admitted polynomial pair of orders four and six must have
nonconstant \(v\) after polynomial depression. Equations (1)–(3)
give the residual within that class. The
[centered rectangle](CENTERED_46_RECTANGLE_NEWTON_FACE_GATE.md) permits
nonconstant order-four leading coefficients. Passing to a coefficient
cover to make the leading coefficient monic can leave the polynomial
ring. Its descent is part of the
[polynomial Ore admission problem](POLYNOMIAL_ORE_INTERTWINER_AND_SUBPRINCIPAL_DESCENT.md).
The calculation does not put a formal position correction in \(A_1\)
or establish a finite matrix realization of a prescribed center map.
