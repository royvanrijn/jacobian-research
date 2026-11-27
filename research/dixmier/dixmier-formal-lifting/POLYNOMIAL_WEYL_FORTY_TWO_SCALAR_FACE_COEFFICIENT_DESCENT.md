# Lower scalar weights of the degree-forty-two Weyl residual

**Theorem PolynomialWeylFortyTwoScalarFaceCoefficientDescent.** Use the
hypothetical actual pair and its Laurent coefficients from
[the degree-forty-two root constraints](POLYNOMIAL_WEYL_FORTY_TWO_SECONDARY_ROOT_CONSTRAINTS.md).
Write \(D_u(E)=E^3-u\) for the factor associated to the
multiplicity-four roots. Then
\[
 \begin{aligned}
 P_3(E)&=D_u(E)^5 C_3(E^3),& \deg C_3&\le5,\\
 P_2(E)&=E D_u(E)^4 C_2(E^3),& \deg C_2&\le5,\\
 P_1(E)&=E^2 D_u(E) C_1(E^3),& \deg C_1&\le8,\\
 P_0(E)&=C_0(E^3),& \deg C_0&\le10.
 \end{aligned}                                                   \tag{1}
\]
The latter two polynomials satisfy
\[
 C_1(u)+uC_1'(u)=0,\qquad C_0'(u)=0.                          \tag{2}
\]
Given the higher \(P\) coefficients, the simple-root condition
restricts \(C_2(v)\) to the four values in (9) below. The result
leaves the degree-forty-two residual open.

## Scalar weights and support

Fix a root \(\lambda\) with \(\lambda^3=u\). In the fractional
Weyl algebra put
\[
 w=E-\lambda,\qquad U=t w,\qquad [U,t]=1/3,\qquad z=U^7/t^2.
\]
The variable \(z\) in the formulas below is in the commutative symbol
ring \(\mathbb C[t,t^{-1},U]\). The positive weights of \(t,U\) are
\(7,2\) and each contraction loses weight \(9\). The first
contraction between two normal-form components with symbols \(f,g\)
is \( (f_Ug_t-f_tg_U)/3\). Products can be ordered with \(t\)
on the left, including negative powers of \(t\), by
\([U,t^i]=(i/3)t^{i-1}\).

The highest symbols are
\[
 p_0=U A(z),\qquad q_0=t B(z),\qquad
 A=a_0+a_1z,\quad \deg B=4,\quad a_0a_1 B(0)\ne 0 .
\]
Their weights are \(2,7\). Their scalar amplitude equation is
\((A+7zA')B-2zAB'=-3\), from \([P,Q]=-1\).
Let \(p_d,q_d\) denote the symbols of the normal-order components
at respective weights \(2-d,7-d\).
For \(j\ge0\), replacing \(t^k w^j\) with \(t^{k-j}U^j\)
has this same symbol at weight \(7k-5j\).
Every correction has \(t\)-degree index \(k-2i\) and \(w\)-index
\(j-i\) when rewritten in \(t,w\) order, so has weight
\(7k-5j-9i\) with \(i\ge1\).
Indeed rewriting a product with \(w=t^{-1}U\) preserves the
leading exponents of \(t,U\), and each contraction lowers
both exponents by one. Consequently for \(0\le d\le8\), the
contribution of \(t^k w^j\) to \(p_d,q_d\) uses only this first
symbol, since the top weights are \(2,7\).

Solving \(7k-5j=2-d\) with \(k\le6\), \(j\ge0\) and the
analogous \(7k-5j=7-d\) with \(k\le21\) gives the following
forms for the components used in the proof. The polynomial letters
in this table denote separate polynomials at their indicated weights.
The constants \(\epsilon_1,\epsilon_3\) are the Taylor coefficients
of \(P_3,P_2\) at order \(4,3\) at \(E=\lambda\).
\[
\begin{array}{c|c|c}
 d & p_d & q_d\\ \hline
 1 & \epsilon_1 U^4/t & U^3 G_1(z)\\
 2 & T(z) & U^6/t\,G_2(z)\\
 3 & \epsilon_3 U^3/t & U^2 G_3(z)\\
 5 & U^2/t\,F_5(z) & U G_5(z)\\
 7 & U/t\,F_7(z) & G_7(z)
\end{array}                                                     \tag{3}
\]
Here \(T,F_5,F_7\) have degree at most one and every \(G_d\) is
a polynomial. The constant terms \(F_5(0),F_7(0)\) are the Taylor
coefficients of \(P_1,P_0\) at orders \(2,1\), respectively.
For all weights \(-d\) used below, the exact commutator equation
gives
\[
 \sum_{i+j=d}\big((p_i)_U(q_j)_t-(p_i)_t(q_j)_U\big)=0.       \tag{4}
\]
Indeed the top weights sum to \(9\), and a second contraction
has weight at most \(-9\).

## The coefficients at weights minus one, three, five and seven

At \(d=1\), divide (4) by \(U^3/t\) and use (3) to obtain
\[
 z(-2A G_1'+6A'G_1)+\epsilon_1(4B-z B')=0.
\]
Evaluation at \(z=0\) gives \(\epsilon_1=0\). The remainder is
\((G_1/A^3)'=0\), so \(q_1=c p_0^3\) for some scalar \(c\).
Replace \(Q\) by \(Q-cP^3\). This kills \(q_1\), preserves
\([P,Q]\) and \(q_0\), and preserves the \(t\)-degree cap since
\(P^3\) has \(t\)-degree at most \(18\).
The bound at the multiplicity-four root already gives
\(\operatorname{ord}_\lambda P_3\ge4\). Hence
\(\operatorname{ord}_\lambda P_3\ge5\).

At \(d=3\), the two middle terms in (4) vanish by the first step.
After dividing by \(U^2/t\), the remaining equation is
\[
 z(-2A G_3'+4A'G_3)+\epsilon_3(3B+z B')=0.
\]
Again \(\epsilon_3=0\), improving
\(\operatorname{ord}_\lambda P_2\ge3\) to
\(\operatorname{ord}_\lambda P_2\ge4\).
Now \(q_3=c' p_0^2\), since \((G_3/A^2)'=0\).
Replacing \(Q\) by \(Q-c'P^2\) kills \(q_3\) while preserving
\(q_0,q_1\), the commutator and the \(t\)-degree cap.

At \(d=5\), (4) consists of the two terms with \(p_0\) or \(q_0\).
After dividing by \(U/t\), it is
\[
 2z(A'G_5-A G_5')+(2B+3z B')F_5+7z B F_5'=0.              \tag{5}
\]
It forces \(F_5(0)=0\), which is \(P_1''(\lambda)=0\).

At \(d=7\), (4) has four potentially nonzero terms.
Divide it by \(1/t\). The term from \(p_0,q_7\) is
\(-2z A G_7'\), and the one from \(p_2,q_5\) is \(2z T'G_5\).
The term from \(p_5,q_2\) before division is
\(U^7/t^3\) times a polynomial in \(z\), so it also has a
factor \(z\) after division. Finally the term from \(p_7,q_0\)
becomes \((B+5z B')F_7+7z B F_7'\).
Its constant term forces \(F_7(0)=0\), or \(P_0'(\lambda)=0\).

The argument applies at each of the three fourth-multiplicity roots.
The deck congruence on \(t^kE^j\) is \(k+j\equiv0\bmod3\).
It gives the forms in (1) from the vanishing orders and
\(\deg_E P_k\le30\). In particular the original scalar weight
gives the factor \(D_u\) in \(P_1\). At \(\lambda^3=u\), the
coefficient of \((E-\lambda)^2\) in
\(E^2D_u C_1(E^3)\) equals
\(9u(C_1(u)+uC_1'(u))\).
Also \(P_0'(\lambda)=3\lambda^2 C_0'(u)\).
The root \(u\) is nonzero, so these are exactly (2).
In particular the same normalized \(Q\) works at all three
fourth-multiplicity roots. The deck transformation sends the
bracket computation at one root to that at another, and fixes
both the pair and the scalars in \(q_1=c p_0^3\) and \(q_3=c' p_0^2\).

## Joint condition at the simple roots

Fix a simple root \(\delta^3=v\), set \(A_\delta=\alpha H'(\delta)^2\),
and use \(\gamma=-P_5'(\delta)/(2A_\delta)\) as in
[the simple-root shear](POLYNOMIAL_WEYL_FORTY_TWO_SECONDARY_ROOT_CONSTRAINTS.md#the-three-simple-roots-and-their-polynomial-jets).
All derivatives in the next display are evaluated at \(\delta\).
Put
\[
 \begin{aligned}
 J_4 &=
     P_4'+\gamma P_5''+\frac{\gamma^2}{2}P_6''',\\
 J_2 &=
     \gamma P_3'+\frac{\gamma^2}{2}P_4''
       +\frac{\gamma^3}{6}P_5'''+\frac{\gamma^4}{24}P_6''''
       -\frac{\gamma}{6}P_5''-\frac{\gamma^2}{6}P_6''' .
 \end{aligned}                                                   \tag{6}
\]
For the shear \(E\mapsto E+h\) with \(h=\gamma/t\),
ordinary Taylor contributions lower \(t\)-degree by one per
factor of \(h\).
The correction with one commutation from an \(E^m\) term is
\[
 \binom m2 [E,h](E+h)^{m-2}
 \quad\text{at first contraction order},\qquad
 [E,h]=-\gamma/(3t^3).
\]
To count it, choose the two positions of an \(E\) followed by
an \(h\) in the expansion, then choose the other \(m-2\) factors
of \(E+h\). The term with one \(h\) altogether lowers the
\(t\)-degree by three and contributes the second derivative
with coefficient \(-\gamma/6\). The term with two \(h\)'s
altogether lowers it by four and contributes the third
derivative with coefficient \(-\gamma^2/6\). Terms with two
contractions lower the degree by at least five.
This proves the exact evaluations
\[
  \widehat P_4'(\delta)=J_4,\qquad
  \widehat P_2(\delta)=P_2(\delta)+J_2.                     \tag{7}
\]
By (1),
\[
 P_2(\delta)=\delta(v-u)^4 C_2(v).
                                                               \tag{8}
\]
Now the per-orbit Euler condition from the previous result gives
\[
 C_2(v)\in
 \left\{
 \frac{\tfrac{A_\delta}{36}
       \bigl((-2+3J_4/A_\delta)^2-m^2\bigr)-J_2}
      {\delta(v-u)^4}
       \;:\;m\in\{1,3,5,7\}
 \right\}.                                                     \tag{9}
\]
Each choice in (9) is the same at the three simple roots.
Indeed under \(\delta\mapsto\zeta\delta\) for a cube root of unity
\(\zeta\), the deck congruence multiplies \(A_\delta,J_4,J_2\)
by \(\zeta\) and multiplies \(\gamma\) by \(\zeta^{-1}\).
The right side of (9) is then unchanged.
