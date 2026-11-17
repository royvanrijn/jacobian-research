# Twisted adjoints exclude the mixed quartic square face

**QAS1.** Let F be a characteristic-zero differential field, with
nonzero derivation D and constant field C. In the finite scalar Ore
algebra F[A;D], put
\[
 q=aA^2+bA+c,\qquad a\ne0,\qquad
 P=gq^2+d,\qquad g\in C^\times,\quad b,c,d\in F.             \tag{1}
\]
If a is not a nonzero constant times a square in F, then P has no
finite differential operator mate [P,T]=kappa with kappa in C^times.
For F=k(V) the conclusion also excludes every finite algebraic
coefficient-cover mate, by normalized coefficient trace.

This is a finite-operator theorem, not a statement about arbitrary
elements of the quotient division ring. The proof uses both the
nonnegative order of a finite differential operator and the leading
coefficient equation above order zero.

The mixed critical face
\[
 [A,V]=\eta/V,\qquad q=A+\beta V^{-3}A^2,\qquad
 P=gq^2+\mu V,\qquad g\beta\mu\eta\ne0                      \tag{2}
\]
is therefore excluded in every finite partner order, despite its
perfect-square top symbol. The direct globally admitted mixed-square
family is excluded as well. A nonsymmetric first-order perturbation
escapes the adjoint argument, but a separate universal Laurent-residue
certificate below excludes it in every finite order as well.

## A rational anti-involution making every quadratic q selfadjoint

Put
\[
                         t=(D(a)-b)/a.
\]
Define f^dagger=f for f in F and
\[
                         A^\dagger=-A+t.                    \tag{3}
\]
Extend this map anti-multiplicatively. It preserves the Ore relation:
the adjoint of Af-fA equals f(-A+t)-(-A+t)f=D(f). Applying it
twice fixes A, so it is an anti-involution of the finite algebra.

The coefficient of A in q^dagger is 2D(a)-2ta-b=b. The constant
correction to c is
\[
 D^2(a)-2tD(a)+(t^2-D(t))a-D(b)+tb=0,                     \tag{4}
\]
after using b=D(a)-ta. Thus q^dagger=q, and consequently
P^dagger=P for every position potential d in (1).

If [P,T]=kappa, applying dagger gives [P,T^dagger]=-kappa. Hence
\[
                         T_-=(T-T^\dagger)/2
\]
is another finite mate and T_-^dagger=-T_-. It is nonzero. Its
leading coefficient acquires the sign (-1)^M under dagger, so its
differential order M must be odd.

On the other hand P has leading coefficient p_4=g*a^2. The highest
commutator coefficient, at A^(M+3), is
\[
                  4p_4D(t_M)-M t_M D(p_4)=0.
\]
Its order is positive for every finite M>=0; it cannot come from
the scalar right side. Thus
\[
                         D(t_M^2/a^M)=0.                    \tag{5}
\]
If M is odd, (5) writes a as a constant times the square of
t_M/a^((M-1)/2), contrary to the hypothesis. Every rational mate
must have even order, while the nonzero skew-adjoint mate has odd
order. This proves QAS1.

More generally, any coefficient-fixed rational anti-involution making
an even-order operator selfadjoint reduces a finite scalar mate to
odd order. A leading-coefficient divisor forcing all rational mate
orders to be even therefore gives the same obstruction. The quadratic
square in (1) supplies both conditions explicitly.

For (2), a=beta*V^-3 is nonsquare up to constants because its
V=0 valuation is odd. The twisting coefficient is exactly
\[
                      t=-3\eta V^{-2}-\beta^{-1}V^3.        \tag{6}
\]
No algebraic square-root position change or infinite gauge is used.
Finite algebraic covers cannot help: trace of a hypothetical cover
mate is a rational mate with the same nonzero scalar commutator.

The leading-square hypothesis is material. For a=1 and q=A^2,
the operator A^4+V^2/2 has the finite mate -A, with commutator eta.
Here a is a square, the mate is odd, and (5) gives no contradiction.

## Exhaustive first orders, and a genuinely nonsymmetric perturbation

Now let
\[
 P_\sigma=gq^2+\mu V+\sigma V^{-2}A,
              \quad q=A+\beta V^{-3}A^2,
              \quad \sigma\in k,\quad g\beta\mu\eta\ne0.    \tag{7}
\]
The new term breaks the adjoint symmetry when sigma is nonzero:
the quartic and cubic coefficients already fix the twist (6),
while the added term changes sign in its first-order part. Its
exclusion below therefore uses a separate argument.

Nevertheless every rational mate through order nine is excluded.
The leading coefficient g*beta^2*V^-6 makes every rational mate
order even. After a constant square-root extension choose a_0
with a_0^2=g, and let R=P_sigma^(1/2) have leading term a_0*q.
Both B_2=R_+ and B_6=(P_sigma*R)_+ have commutators of order at
most two. Finite descending subtraction therefore gives the
exhaustive normal form at order at most eight:
\[
                  T=c_8P_\sigma^2+c_6B_6+c_4P_\sigma
                                      +c_2B_2+c_0,
                       \qquad c_j\in k(a_0).               \tag{8}
\]
This covers arbitrary rational position coefficients, not a bounded
coefficient ansatz. Odd leading orders cannot appear at any step
because their leading coefficient would be proportional to
V^(-3M/2). Order zero is constant.

The exact positive parts have the short form
\[
 \begin{split}
 B_2={}&a_0q,\\
 B_6={}&a_0^3q^3+a_0\mu T_\mu+a_0\sigma T_\sigma
                                  +\frac{3\sigma^2}{8a_0\beta}V^{-1},\\
 T_\mu={}&\tfrac32Vq+\tfrac32\beta\eta V^{-4}A
                       +\tfrac34\eta V^{-1}-2\beta\eta^2V^{-6},\\
 T_\sigma={}&\tfrac32\beta V^{-5}A^3
       +(\tfrac32V^{-2}-\tfrac{21}4\beta\eta V^{-7})A^2\\
       &+(-\tfrac32\eta V^{-4}+7\beta\eta^2V^{-9})A
                   +\tfrac98\eta^2V^{-6}-\tfrac{21}2\beta\eta^3V^{-11}.
 \end{split}                                                       \tag{9}
\]
They follow by solving the square root only through A^-4; every
coefficient needed in (9) is finite. The checker independently
rederives that truncation before forming the commutators.

For sigma nonzero, the decisive coefficients are
\[
 [V^{-21}A^2]\,[P_\sigma,B_6]
                   =-3828a_0^3\beta^3\eta^5\sigma\ne0,
 \qquad
 [P_\sigma,B_2]_{A^2}=a_0\beta\eta\sigma V^{-7}.             \tag{10}
\]
The first cannot be cancelled by B_2, so c_6=0 in (8), and
then its second coefficient forces c_2=0. If sigma=0, use instead
\[
 [V^{-3}A]\,[P_0,B_6]=-3a_0\beta\eta\mu^2\ne0,
 \qquad [P_0,B_2]_{A}=-2a_0\beta\eta\mu V^{-4}.             \tag{11}
\]
Again c_6=c_2=0. Every remaining term in (8) commutes. Thus
\[
 \{T:[P_\sigma,T]=0,\ \operatorname{ord}T\le8\}
                         =k\oplus kP_\sigma\oplus kP_\sigma^2,
                                                                  \tag{12}
\]
and no nonzero scalar mate of order at most eight exists. Descent
of (12) after a constant square-root extension uses highest
differential coefficients. Rational odd orders are impossible,
and trace of an algebraic-cover mate does not increase its order.
These exact first-order calculations are retained independently of
the universal exclusion that follows.

## A universal residue certificate closes the nonsymmetric case

Assume sigma is nonzero and put s=a_0*beta. Give V weight one and
A weight minus two; all Ore contractions preserve this canonical
weight. The lowest term of a_0*q is
\[
 L_-=s\eta^2V^{-7}C(\theta),\qquad
 C(x)=x(x-2),\qquad \theta=V^2A/\eta.                       \tag{13}
\]
For an arbitrary finite order the same leading subtraction as in (8)
gives T=q_0(P_sigma)+sum_(k=0)^N c_k B_(4k+2), with a finite sum
and scalar constants c_k. No bound on N is imposed. The formal
square root still has finite Laurent coefficients at each specified
negative A order, and each commutator with a positive odd root
power has order at most two.

Its other term a_0*A has weight five greater. Relative to the
weight-minus-fourteen leading part L_-^2 of P_sigma, the sigma
perturbation has gap ten and the mu perturbation gap fifteen.
At sigma=mu=0 the complete base operator is the actual square
g*q^2, so all its half-integer powers are the finite operators
a_0^(2k+1)*q^(2k+1). Their commutators vanish exactly, including
every contribution from the higher-weight term of q.

Consequently the uniquely most polar possible anomaly for
B_(4k+2)=(P_sigma^(k+1/2))_+ is linear in sigma and uses only
L_- as its base. Every base correction adds positive gap five,
a mu insertion has gap fifteen instead of ten, and a second
sigma insertion adds another gap ten. These statements follow
coefficientwise from the square-root recursion; its denominator
is only the leading monomial 2sV^-3.

The first sigma response of the square root at that lowest weight is
\[
             R_\sigma=\frac{\sigma}{s\eta}V^3F(\theta),
\qquad C(x+3)F(x)+C(x)F(x-7)=x,
\qquad F(x)=\tfrac12x^{-1}+O(x^{-2}).                       \tag{14}
\]
This is exactly L_-*R_sigma+R_sigma*L_-=sigma*V^-2*A.
Write mathcal L(B)=[x^-1]F(x)B(x), and let
\[
 p_n(x)=\prod_{j=1}^{n}C(x+7j).
\]
Multiplying (14) by an arbitrary polynomial and taking the formal
residue gives the annihilator
\[
 \mathcal L\{C(x+3)B(x)+C(x+7)B(x+7)\}=0.                 \tag{15}
\]
The polynomial right side of (14) has zero residue. Differentiating
the odd root power and translating x in each summand yields the
same moment reduction as in the
[square-face certificate](SQUARE_ORDERED_QUARTIC_FACE.md):
\[
 D_k=\mathcal L(p_{2k}),\qquad D_0=\tfrac12.                \tag{16}
\]
Explicitly its 2k+1 terms are mathcal L(b_i*p_j), i+j=2k,
where b_i=product_(l=0)^(i-1) C(x+3-7l). Equation (15) changes
each term to minus the next; their alternating odd sum is precisely
mathcal L(p_(2k)). Thus no infinite summation or unproved finite
pattern is involved.

For a symbolic n define
\[
 Q_n(x)=\frac{x^2}2+\frac{21n+23}2x+
                         \frac{245n^2+511n+246}4,
\qquad r_n=-\frac{(7n+2)(7n+4)(7n+6)(7n+11)}4.
\]
The exact polynomial certificate is
\[
 C(x+7(n+1))C(x+7(n+2))-r_n
       =C(x+3)Q_n(x)+C(x+7(n+1))Q_n(x+7).                 \tag{17}
\]
Multiply by p_n and apply (15). It proves
mathcal L(p_(n+2))=r_n*mathcal L(p_n) for every n>=0.
In particular
\[
 D_k=-2(7k-6)(7k-5)(7k-4)(14k-3)D_{k-1}\ne0,
                              \qquad k\ge1.               \tag{18}
\]
Every factor is nonzero in characteristic zero; D_1=-66.

The most polar A^-1 coefficient of the odd formal power is
D_k*sigma*s^(2k-1)*eta^(4k)*V^(1-14k). Only that negative
coefficient can produce A^2 in minus its commutator with P_sigma.
Since p_4=s^2V^-6, direct differentiation gives
\[
 [P_\sigma,B_{4k+2}]_{A^2}
 =(56k+2)D_k\sigma s^{2k+1}\eta^{4k+1}V^{-7-14k}
                    +\text{terms of greater V degree}.    \tag{19}
\]
For k=0 this is sigma*s*eta*V^-7; for k=1 it is exactly the
-3828 coefficient in (10). The nonzero term for the largest
index in any finite even-flow normal form is strictly more polar
than every lower-index flow. Thus all noncommuting flow
coefficients vanish. This proves
\[
 \operatorname{Cent}_{k(V)[A]}(P_\sigma)=k[P_\sigma]
                       \quad\hbox{if sigma}\ne0,            \tag{20}
\]
and excludes every finite scalar mate. Highest differential coefficients
descend the scalar coefficients of a polynomial in P_sigma from a
constant square-root extension back to k, so (20) holds over the
original rational coefficient field. Trace extends mate
nonexistence to finite algebraic coefficient covers. When sigma=0
the adjoint theorem already gives all-order mate exclusion; (20)
is only claimed for sigma nonzero.

## The full globally admitted mixed square

Use the actual global operators from
[the global admission note](GLOBAL_ADMISSION_AND_QUARTIC_CRITICAL_FACE.md),
with r=R, v=X^-1, z=a-2/(3r), W=X^3C^2 and lambda nonzero.
That note proves the exact identity
\[
 S_\lambda-\frac4{3r}
                  =-3r(z+\lambda hW)^2+2v.                 \tag{21}
\]
On the generic relative fiber, q_lambda=z+lambda*h*W is a finite
second-order rational differential operator with leading coefficient
\[
                            \lambda h r^2v^{-3}.            \tag{22}
\]
It is nonsquare up to constants, and g=-3r is a differential
constant. QAS1 therefore applies directly to the full S_lambda,
before any parameter rescaling. There is no finite generic
rational-position mate commuting with R, in any order, nor one
over a finite algebraic coefficient cover. In particular no
globally admitted finite conjugate completing this exact-R square
family exists. The twisted adjoint may have parameter poles;
that is harmless for a generic finite nonexistence proof.

For an h-independent scalar nu, adding nu*h*X^2*C is still globally
admitted and has critical term
sigma*V^-2*A with sigma=nu*r: under h=tau^5, v=tau^2*V,
z=tau*A one has tau^-2*h*X^2*C=r*V^-2*A-r*tau*V^-1.
For nu nonzero it escapes the adjoint
symmetry but is excluded by the universal residue proof. Indeed its
integral rescaling has critical reduction (7), whose full rational
centralizer is (20). At each finite parameter pole of a hypothetical
generic mate, subtract the corresponding polar polynomial in the
rescaled position. The pole strictly decreases; finitely many steps
give a forbidden critical mate. This is the finite pole descent of
[QCL1](CRITICAL_LAURENT_FACE_OBSTRUCTION.md), with (20) supplying
the exact centralizer hypothesis. Thus the nonsymmetric globally
admitted perturbation also has no finite scalar mate commuting with R,
in any partner order. This does not exclude arbitrary nonsymmetric
negative-Laurent corrections.

## The integer monomial-pole square family

The entire family
\[
 P_\ell=g(A+\beta V^{-\ell}A^2)^2+\mu V,
 \qquad \ell\ge1\text{ an integer},\qquad
 g\beta\eta\ne0,\qquad \mu\in k                            \tag{23}
\]
has no finite rational-coefficient scalar mate, nor a mate on a
finite algebraic coefficient cover. This includes even ell, when
the leading quadratic coefficient is a square up to constants and
the first parity proof alone does not apply. No full-centralizer
assertion for this integer family is made.

Use (3) to make q=A+beta*V^-ell*A^2 selfadjoint, irrespective of
whether its leading coefficient is square. The twist and the induced
action on H=V^2 A are
\[
 A^\dagger=-A-\ell\eta V^{-2}-\beta^{-1}V^\ell,\qquad
 H^\dagger=-H-(\ell+2)\eta-\beta^{-1}V^{\ell+2}.             \tag{24}
\]
A mate would therefore give a nonzero skew-adjoint finite mate T.
Expand its rational position coefficients at V=0, in the normal
form V^m D(H)+terms of greater V exponent. Here m is an integer
and D is a nonzero polynomial of degree e>=0. The last term of
H^dagger in (24) has strictly positive V exponent. Reversing the
order with V^m consequently sends the leading polynomial to
\[
                  D(-H-(m+\ell+2)\eta).                     \tag{25}
\]
The least exponent m is preserved. Skew-adjointness forces e to
be odd by comparing the highest H coefficients in (25).

Put n=ell+4>=5 and C_2(H)=H(H-2eta). The uniquely lowest
term of P_ell is
\[
 g\beta^2 V^{-2n} C_n(H),\qquad
                         C_n(H)=C_2(H-n\eta)C_2(H).         \tag{26}
\]
Its four roots on the eta-translation orbit are 0,2eta,n*eta,
(n+2)*eta. If m<2n, the negative V^(m-2n) coefficient in the
constant commutator must vanish, giving exactly
\[
             C_n(H+m\eta)D(H)=D(H-2n\eta)C_n(H).            \tag{27}
\]
Record divisors on that orbit by finite Laurent polynomials. For
C_n the character is c(t)=(1+t^2)(1+t^n); for D denote it by
d(t). Equation (27) gives the divisor identity
\[
                  (t^{-m}-1)c(t)=(t^{2n}-1)d(t).            \tag{28}
\]
This is a weaker character test than the primitive-2n criterion:
c vanishes at primitive 2n-th roots, so that stronger criterion
cannot be invoked here. Instead evaluate (28) at a primitive
n-th root xi. Then c(xi)=2(1+xi^2) is nonzero, since its
vanishing would force xi to have order four, whereas n>=5.
Thus n divides m. Comparing the coefficient of H^(e+3) in
(27) gives 4m+2ne=0. Hence e=-2m/n is even, contradicting
the skew-adjoint degree parity from (25).

If m>2n, the commutator has positive V valuation and cannot be
a nonzero constant. At the remaining boundary m=2n, its V^0
coefficient is g*beta^2 times
C_n(H+2n*eta)D(H)-D(H-2n*eta)C_n(H). It has degree e+3,
with nonzero leading coefficient
2n*eta*(e+4) times D's leading coefficient. It is not scalar.
This exhausts all m and proves (23). If roots of unity are
missing from k, extending constants for the character evaluation
is harmless. Finally apply normalized trace to exclude finite
algebraic coefficient-cover mates.

[The exact checker](check_twisted_adjoint_square.py), with
[its retained receipt](twisted_adjoint_square_controls.json), verifies
the rational anti-involution, the universal quadratic symmetry,
the mixed twist and its symmetry-breaking term, the complete
finite fractional-power controls and decisive commutator coefficients,
the symbolic-index certificate (17), and a positive square-leading
Airy control. It also verifies the monomial-pole bottom coefficient,
the weaker character test, adjoint reflection and boundary degree
controls for (23). All-order exclusion uses the written parity proof
or the separately proved universal recurrence, not an extrapolation
from the sixth-order calculation. No quotient-division-ring exclusion,
bounded full Long lift, DC1/DC2 solution, formal verification,
independent whole-proof replay, external review or novelty priority
is claimed.
