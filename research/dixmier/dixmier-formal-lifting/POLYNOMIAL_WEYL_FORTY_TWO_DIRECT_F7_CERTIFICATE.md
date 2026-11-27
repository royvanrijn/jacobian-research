# Direct reduction of the degree-forty-two Weyl residual to the F7 certificate

**Theorem PolynomialWeylFortyTwoDirectF7Certificate.**
The normalized degree-forty-two polynomial Weyl residual is excluded
by the exact \(F_7\) unit-ideal certificate in Proposition 5.2 of
[Niedbala Giraudin, Hexagon 2611.00896v1](https://hexagonmath.org/2611.00896v1).
The proof below derives its input system directly from the Weyl
residual. In particular it does not use the planar exclusion's
left-propagation input from Guccione–Guccione–Valqui, Corollary 7.4,
its planar successor chain, or a planar census.

## The transferred data

Start from the regular centered Weyl residual used in the
[split-prime transfer](POLYNOMIAL_WEYL_FORTY_TWO_SPLIT_PRIME_F7_TRANSFER.md).
Normalize the fourth-multiplicity root to one by scaling the
coordinates. In the notation \(x,y\) for the transferred planar
coordinates, put \(z=xy^3\) and \(h(z)=(z-1)^4(z-c)\), where
\(c\ne0,1\) and \(2c^2-5c+5=0\).
The transfer gives a characteristic-zero polynomial Jacobian pair
\(P,Q\) with nonzero \(J=\{Q,P\}\), support rectangles
with corners \((12,30),(42,105)\), and initial forms
\[
 \operatorname{in}_{3,-1}P=\alpha_0\bigl(xh(xy^3)\bigr)^2,
 \qquad
 \operatorname{in}_{3,-1}Q=\beta_0\bigl(xh(xy^3)\bigr)^7 .
                                                               \tag{1}
\]
Their maximal-coordinate faces and their \((-1,1)\) faces
are single corners. The root-ratio equation is a finite
coefficient condition and can be kept in that transfer.

We can impose in addition the right cap in the shifted coordinates
\[
       x=T^3,\qquad y=\eta+\omega T^{-1},\qquad \omega^3=1,
       \qquad W(T^i\eta^j)=7i-12j:
       \quad W(P)\le2,\quad W(Q)\le7.                         \tag{2}
\]
Indeed this is the
[actual scalar successor](POLYNOMIAL_WEYL_BERNSTEIN_FORTY_TWO_REAL_BOUND.md)
of the Weyl residual. The finite-pole centering translation
preserves this cap: for \(X_{\rm old}=X+x_0\), one has
\(t_{\rm old}=(t^3+x_0)^{1/3}=t+O(t^{-2})\) and the
shifted \(Y\)-coordinate changes by
\(\omega(t^{-1}-t_{\rm old}^{-1})=O(t^{-4})\).
Thus every correction to \(t\) or to the shifted coordinate
has strictly lower weight. The ramified Ore relation also
has positive contraction loss, so their finite re-expansions
retain the cap. A reciprocal linear scaling in \(X,Y\)
normalizing the nonzero root only rescales the shifted
variables by nonzero constants.
In the ramified Weyl ring with
\(t=X^{1/3}\), the shift in \(Y\) is by \(\omega t^{-1}\).
For the shift on the center in prime characteristic \(p>3\),
the formula
\((Y+g)^p=Y^p+g^p+\partial_X^{p-1}g\) follows from
[Jacobson's formula as in Lauritzen–Thomsen, §2.2](https://arxiv.org/pdf/1508.02499v1):
after the first commutation with a multiplication operator,
a further such commutation vanishes, so the only mixed
term is the \((p-1)\)-st derivative of \(g\).
For \(g=\omega X^{-1/3}\) this derivative vanishes, since
\(-1/3\) is not \(-1\) in the prime field and the derivative
coefficient contains the factors at all integers from
zero through \(p-2\).
Consequently the central shift is
\(Y^p\mapsto Y^p+\omega^p(t^p)^{-1}\).
The weight in (2) has positive contraction loss nine in
the ramified Weyl ring. The cap therefore passes to
the \(p\)-th powers with weights divided by \(p\) on
the central variables. The algebraic and split-prime
arguments in the transfer can include \(\omega\) and
the finite coefficient equations (2), both in the Weyl
specialization and in the return to characteristic zero.
Thus we use (2) directly, with no planar successor theorem.

## A polynomial bracket companion for the left side

Use \(t=x^{1/3}\), \(E=ty\) in commutative symbols.
The \(t\)-degree is \(v_{3,-1}\) on the original polynomial
ring. Starting from \(Q^2-(\beta_0^2/\alpha_0^7)P^7\),
perform the same ordered-degree cancellation as in the
[normalized cusp theorem](POLYNOMIAL_WEYL_FORTY_TWO_NORMALIZED_CUSP_RELATION.md),
now with the Poisson bracket.
For clarity, if the current largest \(t\)-degree is \(K>17\),
the homogeneous top coefficient \(A(E)\) commutes with
\(t^6\alpha_0 h(E^3)^2\), since the bracket with \(P\)
still has \(t\)-degree \(21\).
The equation is \(A'=(K/3)(H'/H)A\), where \(H=h(E^3)\),
and a simple root of \(H\) forces \(K=3j\).
Since \(6\le j\le13\), the power \(H^j\) is a product of
powers of \(H^2,H^7\), and subtraction of the corresponding
polynomial in \(P,Q\) cancels it.
Only terms linear in \(Q\) occur among the subtractions
that involve \(Q\), so completing the square gives
\[
 S=(Q^\ast)^2-F(P),\qquad
 Q^\ast=Q-\tfrac12 D(P),\qquad
 \deg D\le3,\quad \deg F=7.
\]
It satisfies
\[
             \{P,S\}=-2J Q^\ast,\qquad
             \{Q^\ast,S\}=-J F'(P).
\]

The \(t\)-degree of \(S\) is exactly \(17\).
The coefficient at that degree is \(H^5 B(E)\), with
\(B'-(2/3)(H'/H)B\) a nonzero scalar.
At a root of \(H\) of multiplicity \(1\) or \(4\) its
\(t^{17}\)-coefficient has vanishing order exactly
\(5\mu+1\), by the same leading-order equation as
in the normalized cusp theorem, and the comparison
at infinity gives \(\deg B=10\).
Thus \(x^{34}y^{85}\) occurs in \(S\).
At its lowest Euler grade \(a\le-51\), its component
commutes with the lowest component of \(P\),
because the bracket with \(P\) has lowest grade
\(-63\), while \(a-18<-63\).
The latter initial form is the pure monomial
\(\alpha_0 x^{12}y^{30}\) by the split-prime data.
The bracket degree equation gives that the
largest \(y\)-exponent \(e\) of the component of
\(S\) there satisfies \(e=-5a/3\).
Indeed a monomial \(x^{a+e}y^e\) commutes with
\(x^{12}y^{30}\) exactly when
\(30(a+e)=12e\).
The \(t\)-cap \(17\) gives \(-a/3\le17\).
Thus \(a=-51\), and the inequalities \(i-j\ge-51\),
\(3i-j\le17\) show that \(S\) is subrectangular
with corner \((34,85)\).

For a left weight \((-\theta,1)\), \(\theta\ge1\),
apply the same observations also after any shear
\(x\mapsto x+\zeta y^{-2}\).
The corner monomials of \(P,Q^\ast,S\) persist,
since the shear decreases \(x\)-exponents and
none has terms beyond its maximal \(x\)-degree.
As long as \(\theta<5/2\), all three weights are
therefore positive.
Their weight ratios are \(6:21:17\) throughout
this interval. Indeed the leading bracket of
\(P,Q^\ast\) vanishes by its scalar right side
and \(1-\theta\le0\), so their ratio is constant
by continuity and proportional monomial exponents
on the open intervals between edge weights.
Its initial value is \(7/2\).
Write \(v\) for the weight of \(P\).
Whenever the weight of \(S\) exceeds \(5v/2\),
the equation \(\{P,S\}=-2JQ^\ast\) forces its
leading bracket to vanish too. The corresponding
ratio is locally constant in the same way,
and starts at \(17/6>5/2\). It therefore never
leaves that region.
Now the positive weighted leading forms commute,
so by weighted Euler differentiation and unique
factorization they are powers of a common form.
Because the integer ratio \(6:21:17\) is primitive,
the three forms are proportional to its sixth,
twenty-first and seventeenth powers.
The initial form of the original \(Q\) equals
that of \(Q^\ast\) on these weights because
\(\deg D\le3<7/2\).
The same reasoning continues through any further
interval over which the three weights stay positive
after a Laurent shear.

## The necessary left subring

For the original unshifted pair, the first left edge
with \(\theta<5/2\), if present, starts at
\((12,30)\). By the sixth-power conclusion, the
common factor starts at \((2,5)\), and any other
point is \((1,5-\theta)\) or \((0,5-2\theta)\).
Since the first left edge has \(\theta>1\) by the
pure Euler face, it must have slope parameter
\(3/2\) or \(2\).

The case \(\theta=3/2\) would have a common factor
proportional to \(R=x^2y^5+\gamma y^2\) with
\(\gamma\ne0\). We recall the short filtered
Poisson argument in
[Guccione–Guccione–Valqui, Lemma 2.2](https://arxiv.org/pdf/1401.1784v3).
Among homogeneous polynomials of any fixed weight,
those commuting with a homogeneous form of
nonzero weight span at most one dimension.
This follows by applying the weighted Euler
derivation to the proportional gradients.
The span of \(P^iQ^j\) for \(0\le i,j\le n\)
has dimension \((n+1)^2\) and only linearly many
possible leading weights. Thus some element has
a nonzero leading bracket with the top of \(P\).
The operator \(\{\cdot,P\}\) is locally nilpotent
on \(\mathbb C[P,Q]\), so its last nonzero
leading-bracket step gives a homogeneous polynomial
\(A\) with \(\{A,R^6\}\ne0\) and
\(\{\{A,R^6\},R^6\}=0\).
Since \(x^2y^3+\gamma\) is irreducible and occurs
to multiplicity one in \(R\), weighted Euler
differentiation and unique factorization make
\(\{A,R^6\}\) a scalar multiple of an integral
power of \(R\). Thus after scaling, the
polynomial \(\{A,R\}\) equals \(R^i\) for some
integer \(i\ge0\).

Rescale \(\gamma\) to one. Put \(z_0=x^{-2}y^{-3}\)
and \(r=1+z_0\). The homogeneous equation at
weight \((-3,2)\) has
\(A=x^{2i-1}y^{5i-4}F(z_0)\), with \(F\) Laurent,
and gives
\[
      3F r+(4i-5)z_0F r'-4z_0F' r=r^i .
\]
The lowest exponent comparison forces \(F\)
to be polynomial, since a negative exponent
\(k\) gives coefficient \(3-4k\ne0\).
For \(d=\deg F\), the coefficient at \(d+1\)
has factor \(4i-2-4d\ne0\), forcing \(d\le i-1\).
Successively divide the equation by \(r^j\)
and evaluate at \(z_0=-1\) for \(0\le j<i\).
If \(F=r^j F_j\), the surviving coefficient
is a nonzero multiple \(4i-5-4j\) of \(F_j(-1)\),
so \(F_j(-1)=0\). Hence \(r^i\mid F\), a
contradiction. This excludes \(\theta=3/2\).

If \(\theta=2\) is the first edge, its common
factor is proportional to \(y\,g(xy^2)\)
with \(g(T)=(T-\alpha)(T-\beta)\).
If no edge occurs before \(5/2\), set
\(\alpha=\beta=0\); the same expression gives
the monomial at the corner for \(\theta=2\).
Put \(v=y^{-1}\), \(u_0=xy^2-v\) and
\[
                Y_1=u_0-\alpha,\qquad
                Y_2=g(u_0)/v.
\]
The needed assertion is
\[
                    P,Q\in \mathbb C[v,Y_1,Y_2].           \tag{3}
\]

To prove it, let \(\zeta\) be a root of \(g\) of
multiplicity \(e_\zeta\), and write
\(\widetilde x=x-\zeta y^{-2}\). In this coordinate
the pair is obtained by a Laurent shear.
If \(e_\zeta=2\), the common form at \(\theta=2\)
becomes proportional to \(\widetilde x^2 y^5\).
Any first edge before \(5/2\) would require
\(2<\theta<5/2\) with \(\theta\) or \(2\theta\)
integral, which is impossible.
The pair has nonpositive weights at \((-5,2)\)
by continuity. If \(e_\zeta=1\), the translated
common factor at \(\theta=2\) is
\(y(\widetilde x y^2)(\widetilde x y^2+\zeta-\zeta')\),
where \(\zeta'\) is the other root.
Its smaller \(x\)-endpoint is \((1,3)\).
The three actual supports from the respective
powers \(6,21,17\) then stay positive until
\(\theta=3\), so the left sixth-power argument
continues on that interval. A first edge
\(2<\theta<3\) from the normalized endpoint
\((1,3)\) would require an integral \(\theta\),
so no such edge exists. The pair therefore
has nonpositive weights at \((-3,1)\) in this
translated coordinate.

Finally write any polynomial in \(x,y\) as
\(\sum_d v^d\phi_d(u_0)\) with \(d\) integral.
It belongs to the ring in (3) exactly when
\(g(u_0)^a\mid \phi_{-a}(u_0)\) for \(a>0\),
since the ring is generated by \(v,u_0,g(u_0)/v\).
At a root \(\zeta\) this says that all terms
have nonnegative valuation for the monomial
valuation with \(v\)-value \(e_\zeta\) and
\((u_0-\zeta)\)-value one.
For \(e_\zeta=2\) this is the restriction of
the monomial valuation
\((\widetilde x,y)\mapsto(5,-2)\);
indeed \(u_0-\zeta=\widetilde x y^2-v\)
has leading term \(\widetilde x y^2\) there.
For \(e_\zeta=1\) it is the restriction of
\((\widetilde x,y)\mapsto(3,-1)\).
The initial forms of \(v,u_0-\zeta\) are
algebraically independent in both cases,
so the restrictions are the asserted monomial
valuations. The two nonpositivity bounds above
are exactly the needed nonnegative valuations.
This proves (3).

## Exact alignment with the F7 certificate

Put \(\delta=\alpha-\beta\) and
\(s=xy^3-1-\alpha y\). Then
\[
 x=(s+1)v^3+\alpha v^2,\qquad y=v^{-1},\qquad
 Y_1=sv,\qquad Y_2=s(sv+\delta),\qquad
 Y_1^2+\delta Y_1=vY_2.
\]
A basis of the ring in (3) consists of
\(v^iY_2^jY_1^\epsilon\), with \(i,j\ge0\) and
\(\epsilon=0,1\). Its key is
\((b,a)=(2j+\epsilon,i+j+\epsilon)\), so
\(0\le b\le2a\). The basis element at that
key is
\[
 e_{b,a}=\sum_{r=0}^{\lfloor b/2\rfloor}
    \binom{\lfloor b/2\rfloor}{r}\delta^r s^{b-r}v^{a-r}.
                                                               \tag{4}
\]
Spanning follows from the quadratic relation in \(Y_1\);
independence follows from the distinct leading monomials
\(s^b v^a\) in (4).

At weight \((3,-1)\) the \(t\)-degree of a basis element
is exactly \(a\), with initial form \((xy^3-1)^b v^a\).
The symbols at each \(a\) are linearly independent
as functions of \(xy^3\). Thus (1) bounds the support
by \(a\le6\) for \(P\), \(a\le21\) for \(Q\), and
fixes the coefficient polynomials at those two highest
values of \(a\).

For the right weight (2), substitute
\(y=\omega T^{-1}(1+\omega^{-1}T\eta)\).
The leading terms of \(s,v\) are
\(3\omega^{-1}T\eta,\ \omega^{-1}T\), of weights
\(-5,7\). The subtraction of \(\alpha y\) in
\(s\) has smaller weight.
All terms with \(r>0\) in (4) have smaller weight by
\(2r\). The leading monomials are distinct for all
keys. Thus if the coefficient of \(e_{b,a}\) is
nonzero, \(7a-5b\le2\) for \(P\), and \(7a-5b\le7\)
for \(Q\).

In these coordinates the bracket
\(\{P,Q\}_{s,v}\) equals \(-J' v\) for \(J'\ne0\)
by \(\partial(x,y)/\partial(s,v)=-v\).
The possible top-weight keys of \(P\) are \((1,1),(8,6)\);
those of \(Q\) are \((7j,1+5j)\), \(0\le j\le4\).
Set \(w=s^7v^5\) and \(H_1=(1-c)^{-1}\), so
\(2H_1^2+H_1+2=0\).
Scaling the two polynomials makes their top
\(a\)-coefficient polynomials
\[
     s^8(1+s)^2(1+H_1s)^2,\qquad
     s^{28}(1+s)^7(1+H_1s)^7 .
\]
The respective weight-top parts are then
\(sv(A+w)\) and \(vG(w)\), with \(G\) monic of
degree four. The weight-seven bracket equation is
\[
       (A+8w)G-2w(A+w)G'=-J'.
\]
Its constant coefficient shows \(A\ne0\);
a rescaling of \(v,\delta\) and the two polynomial
scales makes \(A=1\) and retains the fixed \(a\)-tops.
The recurrence with \(G\) monic now gives coefficients
\((35/128,35/16,35/8,7/2,1)\) and
\(J'=-35/128\).
Subtract constants from the two polynomials.
Then subtract multiples of \(P^3,P^2,P\) from \(Q\)
to kill its keys \((3,3),(2,2),(1,1)\),
in that order. The top term \((sv(1+w))^k\)
has coefficient one at \(s^kv^k\), and every
\(P^k\) for \(k\le3\) has \(a\le18\), weight at
most six, so these operations preserve all
other top normalizations.

These are exactly the \(F_7\) inputs Psupp,
Qsupp, Pfix, Qfix and J in the machine function
of the pinned ancillary verify_F7F8.py.
Formula (4) gives its expansion function tail.
The bracket of two ordinary monomials is
\[
 \{s^b v^a,s^{b'}v^{a'}\}
       =(ba'-ab')s^{b+b'-1}v^{a+a'-1}.
\]
At each next smaller weight the new unknowns
therefore occur in the bracket with the known top
parts only. All their expansion tails have smaller
weight. The pinned run function performs Gaussian
elimination on precisely these linear columns over
\(\mathbb Q[H_1]/(H_1^2+H_1/2+1)\), and its remaining
rows give necessary compatibility polynomials in the
free coefficients and \(\delta\).
For the initial fixed coefficients, none of the
lower tails from a nonzero fixed key lands on
another nonzero fixed key in \(F_7\), so their
initial expansion in run agrees with (4).
Indeed the \(a\)-top tails lower \(a\), and to land
on a lower key of top weight would require a
source \(b<8\) or \(b<28\) respectively. Between
top-weight keys a tail is impossible because it
strictly lowers weight.

Proposition 5.2 and its ancillary exhibit and
expand the identity \(1=\sum_i h_i f_i\) for eight
of these necessary \(F_7\) conditions, over the
same quadratic number field. The
[targeted F7 checker](verify_direct_f7_certificate.py)
independently assembles the bracket coefficients by expanding
and differentiating the two polynomials, verifies equality of
all seventeen compatibility generators to those of the pinned
script, and checks the exact unit-ideal expansion.
Hence the system derived here has no solution, independently
of the other planar shape reductions in the preprint.

For any characteristic-zero coefficient field, the data descend to a
finitely generated subfield that embeds into \(\mathbb C\). Generation
of the Weyl algebra after this extension descends by faithful flatness
on the vector-space quotient by the generated subalgebra. Thus the
exclusion holds over the original coefficient field.
