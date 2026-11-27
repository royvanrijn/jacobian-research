# Actual coefficient constraints at the degree-forty-two roots

**Theorem PolynomialWeylFortyTwoSecondaryRootConstraints.** Assume the
hypothetical proper degree-\(42\) Weyl pair from the
[degree-forty-two residual](POLYNOMIAL_WEYL_BERNSTEIN_FORTY_TWO_REAL_BOUND.md)
exists. Use its Fourier coordinates \([Y,X]=1\), normalized by
\([Q,P]=1\). In the fractional Weyl algebra write
\[
  t=X^{1/3},\qquad E=tY,\quad [E,t]=1/(3t),\qquad
  P=\sum_{k\le6}t^k P_k(E),\quad Q=\sum_{k\le21}t^k Q_k(E).
\]
The sums are finite and the coefficient functions are polynomials.
Set
\[
 H(E)=(E^3-u)^4(E^3-v),\qquad P_6=\alpha H^2,\quad
 Q_{21}=\beta H^7,\qquad \alpha\beta\ne0 .
\]
Here \(u,v\) are the two roots with distinct multiplicities from
the degree-forty-two residual. Then
\[
        Q_{20}=\frac{7\beta}{2\alpha}H^5P_5,\qquad
        P_5=E(E^3-u)^7(E^3-v)(A_0+A_1E^3).
                                                            \tag{1}
\]
Also \(P_4=E^2(E^3-u)^6 B(E^3)\), with \(\deg B\le3\), and
\(B(v)=(A_0+A_1v)^2/(4\alpha)\).

At each simple root \(\lambda^3=v\) of \(H\), put
\(A=\alpha H'(\lambda)^2\) and
\(\gamma=-P_5'(\lambda)/(2A)\).
The exact shear \(E\mapsto E+\gamma/t\) gives an operator
with coefficients \(\widehat P_k(E)\) in the same normal form.
One has \(\widehat P_3(\lambda)=0\), and
\[
  \left(-2+\frac{3\widehat P_4'(\lambda)}{A}\right)^2
             -\frac{36\widehat P_2(\lambda)}{A}
                         \in\{1,9,25,49\}.                \tag{2}
\]
Equation (2) comes from the full commuting product at a
filtration with zero contraction loss.

At a fourth-multiplicity root \(\lambda^3=u\) used for the
scalar successor, let \(a_0,a_1,c_4\) be the coefficients of
\(t(E-\lambda)\), \(t^6(E-\lambda)^8\), and \(t^{21}(E-\lambda)^{28}\)
in the appropriate scalar faces before writing them in \(E\)
translated by \(\lambda\). Then they are nonzero and satisfy
\[
                            a_0^5=-\frac{384}{35}\frac{a_1^4}{c_4}.
                                                            \tag{3}
\]
These conditions constrain actual Laurent coefficients at
both root orbits. They leave the global Weyl pair open.

## The normal coefficients and the fourth-multiplicity orbit

The subrectangular pair in the degree-forty-two residual
has \(Y\)-degrees at most \(30,105\), and the first
edge is \((3,-1)\). Replacing \(Y\) by \(t^{-1}E\)
gives finite Laurent normal forms. In fact
\([E,t^k]=(k/3)t^{k-2}\) for any integer \(k\), so a
contraction lowers \(t\)-degree by two. This keeps
\(\deg_E P_k\le30,\ \deg_E Q_k\le105\). It also
transports the first edge to \(t\)-degrees \(6,21\),
with leading coefficients \(P_6,Q_{21}\) as displayed.
The cubic deck transformation \(t\mapsto\zeta t\),
\(E\mapsto\zeta E\) fixes the original \(X,Y\).
Consequently every power \(E^j\) appearing in \(P_k\)
or \(Q_k\) satisfies \(j+k\equiv0\pmod3\).

The scalar successor of slope \(4/7\) in the \(X,Y\)
weights is obtained by the shift \(E\mapsto E+\lambda\)
at a fourth-multiplicity root of \(H\).
Equivalently it uses weights \(7,-5\) on \(t,E-\lambda\).
The loss of a contraction is \(9\), since
\([E-\lambda,t]=1/(3t)\). The positive weights at this
successor are \(2,7\). Thus all coefficients of the
original \(P,Q\), expanded about \(E=\lambda\), satisfy
\[
 \operatorname{ord}_{E=\lambda}P_k\ge
            \left\lceil\frac{7k-2}{5}\right\rceil,\qquad
 \operatorname{ord}_{E=\lambda}Q_k\ge
            \left\lceil\frac{7k-7}{5}\right\rceil
\]
when the lower bound is positive. This applies at all
three roots of \(E^3-u\) by deck covariance.
In particular \(P_5\) is divisible by \((E^3-u)^7\),
and \(P_4\) by \((E^3-u)^6\).

For the first equation in (1), compare the coefficient
of \(t^{24}\) in the commutator. Only the pairs of grades
\((6,20)\) and \((5,21)\) contribute, because the top
grades sum to \(27\) and contractions lose even amounts.
The first-contraction expression for two terms
\(t^i A_1(E),t^j B_1(E)\) is
\(t^{i+j-2}(j A_1' B_1-i A_1B_1')/3\).
Hence, with \(V=Q_{20}\), the equation is
\[
 \alpha(40 HH' V-6H^2V')
         +\beta(21P_5'H^7-35P_5 H^6H')=0 .
\]
The stated \(V\) is a particular solution.
The difference satisfies \(3V'=20(H'/H)V\).
A nonzero polynomial difference would have order
\(20/3\) at a simple root of \(H\), which is impossible.

## The three simple roots and their polynomial jets

Fix a simple root \(\lambda^3=v\) of \(H\), put
\(w=E-\lambda\), and give \(t,w\) weights \(1,-s\),
respectively. The loss of one contraction is \(2-s\).
For \(0<s<2\) a potential highest bracket of positive
weight must vanish, since the actual commutator is
scalar. Write the top symbols as \(t^A f(\eta)\),
\(t^B g(\eta)\) with \(\eta=w t^s\). Vanishing of the
bracket gives \(B f'g-A f g'=0\).
Thus if \(f,g\) are nonconstant their degrees are in
the ratio \(A:B\).

First take \(s=1/2\). The leading weight of \(P\) is \(5\)
and its symbol is
\[
                t^5\big(\alpha H'(\lambda)^2\eta^2+P_5(\lambda)\big).
\]
The term from \(Q_{21}\) gives weight \(35/2\) in \(Q\),
with a nonzero \(\eta^7\) coefficient.
Any larger weight \(B>35/2\) would make the degree in
\(\eta\) at most six, because \(B+s\deg g\le21\).
But the vanishing Poisson bracket with the quadratic
symbol requires \(\deg g=2B/5>7\). A constant \(g\)
also cannot commute with that quadratic at \(B>0\),
by the same bracket equation. Hence the actual
weight is \(35/2\). The top bracket, of weight
\(5+35/2-(2-1/2)>0\), is zero and gives \(f^7\)
proportional to \(g^2\). Therefore \(f\) is a square
up to scale, forcing \(P_5(\lambda)=0\).

With \(s=1\), the same argument gives weights \(4,14\)
and makes the quadratic
\[
                A\eta^2+P_5'(\lambda)\eta+P_4(\lambda)
\]
a square. It follows that
\(P_5'(\lambda)^2=4A P_4(\lambda)\).
Now deck covariance and \(\deg_E P_5,\deg_E P_4\le30\)
give the forms in the statement: \(P_5\) also has
the simple factors \(E^3-v\) just obtained,
and its exponents are \(1\) modulo three.
The exponents of \(P_4\) are \(2\) modulo three.
Substituting \(\lambda^3=v\) in the discriminant
equation gives \(B(v)=(A_0+A_1v)^2/(4\alpha)\).

Apply the exact shear \(E\mapsto E+\gamma/t\). It
preserves the scalar commutator and the highest
\(t\)-degree and coefficient of both operators.
It makes \(\widehat P_5\) vanish at least twice at
\(\lambda\), and \(\widehat P_4\) vanish there.
Indeed its changes at degrees five and four are
the ordinary first and second Taylor terms from
\(P_6\) and \(P_5\). The first noncommutative correction
uses \([E,\gamma/t]=-\gamma/(3t^3)\), which lowers
degree by three. At \(s=3/2\) the same weight comparison
now gives a quadratic with no linear term:
\[
                    t^3(A\eta^2+\widehat P_3(\lambda)).
\]
As before \(Q\) has weight exactly \(21/2\), since
a larger weight contradicts the \(t\)-degree cap.
The quadratic must be a square. Thus
\(\widehat P_3(\lambda)=0\), which gives the explicit
jet relation
\[
 P_3(\lambda)+\gamma P_4'(\lambda)
 +\frac{\gamma^2}{2}P_5''(\lambda)
 +\frac{\gamma^3}{6}P_6'''(\lambda)
 -\frac{\gamma}{6}P_6''(\lambda)=0 .
\]
For the last term, in a monomial \(E^j\) its
single-shear-factor contributions have
\(\binom j2[E,\gamma/t] E^{j-2}\) at that degree.

## The Euler condition at a simple root

In the sheared pair, for \(3/2<s<2\) the highest
term of \(P\) has weight \(6-2s\) with the sole
top monomial from \(t^6P_6\). The vanishing orders
just proved make every lower \(t\)-degree term
strictly lower in that interval. The same
Poisson degree argument gives weight \(21-7s\)
for \(Q\). Continuity at \(s=2\) gives weights
\(2,7\), respectively. Now the relation
\([w,t]=1/(3t)\) is homogeneous of weight \(-1\);
the ring is graded at this weight and the
highest components must commute exactly, since
the actual commutator has weight zero.

Put \(C=3t^2w\), so \([C,t]=t\). The two components
have the forms \(t^2 f(C),t^7 g(C)\) with
\(\deg f=2,\deg g=7\). Their degrees follow from
the terms in \(P_6,Q_{21}\) at the simple root
and the bounds on \(t\)-degree. In particular,
using \(w^2=t^{-4}(C-2)C/9\), one has
\[
       f(C)=\frac A9 C(C-2)
                   +\frac{\widehat P_4'(\lambda)}3 C
                   +\widehat P_2(\lambda).
\]
Commutation gives \(f(C+7)g(C)=f(C)g(C+2)\).
In an orbit of roots under translation by integers,
write \(F(z),G(z)\) for the Laurent polynomials
whose coefficients count their multiplicities for
\(f,g\), respectively. The relation on roots is
\[
                  (1+z+\cdots+z^6)F(z)
                               =z^5(1+z)G(z).
\]
Thus \(1+z\) divides \(F\) on every orbit, so the
two roots of \(f\) lie on one orbit and are distinct
with odd gap \(m\). After a monomial shift
\(F=1+z^m\). If \(m\ge9\), the coefficient of \(z^7\)
in \((1+z+\cdots+z^6)(1-z+\cdots+z^{m-1})\) is
\(-1\), contradicting the nonnegative multiplicities
for \(G\). Therefore \(m\in\{1,3,5,7\}\).
The discriminant of \(f\) divided by the square
of its leading coefficient equals \(m^2\). This
gives (2).

## The amplitude at the scalar successor

Use the fourth-multiplicity root of the scalar
successor and write \(w=E-\lambda\).
The lattice slopes force the two leading symbols
to have the form
\[
    t w(a_0+a_1 z),\qquad
    t\sum_{j=0}^4 c_j z^j,\qquad z=t^5w^7 .
\]
Indeed the monomials in \(P\)'s face satisfy
\(7k-5j=2\) and \(k\le6\), and those in \(Q\)'s
face satisfy \(7k-5j=7\) and \(k\le21\).
A nonzero scalar bracket makes
\(a_0 c_0\ne0\); the highest \(t\)-degree terms
give \(a_1c_4\ne0\). Set \(U=tw\), so
\(\{U,t\}=1/3\) and \(z=U^7/t^2\).
With \(A_1(z)=a_0+a_1z\), \(B_1(z)=\sum c_jz^j\),
the scalar condition \([Q,P]=1\) becomes
\[
            (A_1+7zA_1')B_1-2zA_1B_1'=-3 .
\]
It gives \(a_0c_0=-3\) and
\[
        (1-2j)a_0 c_j+(10-2j)a_1c_{j-1}=0
                         \quad (1\le j\le4).
\]
In particular \(c_4=(128/35)(a_1/a_0)^4 c_0\),
which is (3). The resulting face coefficients
are consistent with a scalar Poisson bracket.
The relations in the preceding sections arise
from matching this face with other layers and
root orbits of the same finite operators.

## Source scope

The fractional coordinate algebra and leading-bracket
convention are those of
[Guccione–Guccione–Valqui II, Definition 2.9 and Proposition 4.3](https://arxiv.org/pdf/1205.6827v1).
The leading face and the first scalar successor are
the hypotheses supplied by the degree-forty-two residual.
The simple-root jet and Euler calculations use the
finite normal forms of that same pair.
