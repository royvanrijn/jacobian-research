# The primitive homogeneous resultant of three consecutive quadratic moments

Put \(Q=ax^2+bxy+cy^2\), use the factorial functional
\(\Lambda_2(x^iy^j)=i!j!\), and set
\[
 L_j=\frac{\Lambda_2(Q^j)}{(j!)^2},\qquad
 E_j=\binom{2j+1}{j}.
\]

**Theorem.** The forms \(L_j\) are primitive integer forms, and
for \(n\ge1\) their ternary homogeneous resultant is
\[
 \left|\operatorname{Res}(L_n,L_{n+1},L_{n+2})\right|
     =2^{4n^2+6n} E_n^{(n-1)(n+2)} E_{n+1}^{n^2}.          \tag{1}
\]
Thus the primes of the normalized deletion system are explicit
at every start. At odd primes their valuations are given by
the carries in \(n+(n+1)\) and \((n+1)+(n+2)\).
The formula applies equally to composite starting and terminal
orders.

Here is a local consequence for the affine coefficient scheme.
For a triple of consecutive orders \(I\) and a prime \(p\),
put
\[
 \omega_I(p)=\max_{j\in I,\ 1\le r\le j}
    \frac{v_p(j!)-v_p((j-r)!)+v_p(r!)}{r}.
\]
If \(p\) is outside the explicit support (1) for that triple,
a normalized phase \(Q+1\) with those three moments zero
has \(v_p(a),v_p(b),v_p(c)\ge-\omega_I(p)\)
at every place over \(p\).
Either of the two triples inside a four-window can give this
coefficient bound.

## Primitivity

The coefficient of \(a^i b^\ell c^s\) with \(i+\ell+s=j\)
in \(L_j\) is
\[
              \frac{(2i+\ell)!(2s+\ell)!}{j!\,i!\ell!s!}.
\]
It is an integer. For a prime-power term of the factorial
valuation formula, reduce \(i,\ell,s\) to their remainders
and divide those remainders by the prime power. The integer
parts cancel. For the resulting \(x,z,y\) in \([0,1)\),
the inequality needed is
\(\lfloor2x+z\rfloor+\lfloor2y+z\rfloor
 \ge\lfloor x+y+z\rfloor\).
The larger of \(2x+z,2y+z\) is at least their average,
and the other floor is nonnegative. The coefficient of
\(b^j\) is exactly one, proving primitivity.

## Factoring the resultant before evaluation

Set \(F_j=\Lambda_2(Q^j)/(2j)!\) and
\(C_j=\binom{2j}{j}\), so \(L_j=C_j F_j\).
Work temporarily in the coordinates \(a,t,c\) with \(b=2t\).
Write \(\alpha=c-t,\ \beta=a-t,\ D=a+c-2t\) and
\(\Delta=ac-t^2\). The homogeneous interval recurrence gives
\[
       D F_j=\frac{2j}{2j-1}\Delta F_{j-1}+T_j,\qquad
       T_j=\alpha c^j+\beta a^j .
\]
One can obtain it by integrating the derivative of
\((z-\alpha/D)(a z^2+2t z(1-z)+c(1-z)^2)^j\)
and then clearing denominators. This is a polynomial
identity; in particular it holds on \(D=0\).
[Homogeneous quadratic moment recurrence](BINARY_HOMOGENEOUS_QUADRATIC_STRONG_FACTORIAL.md#1-the-interval-identity)

Apply the usual resultant rules for addition of a homogeneous
multiple of another equation, multiplicativity in a slot,
and restriction to a plane. In absolute value,
\[
\begin{split}
 \operatorname{Res}(F_n,T_{n+1},T_{n+2})
   ={}&\operatorname{Res}(F_n,D,T_{n+2})\,
       \operatorname{Res}(F_n,F_{n+1},D)\,
       \operatorname{Res}(F_n,F_{n+1},F_{n+2}),\\
 \operatorname{Res}(F_n,T_{n+1},T_{n+2})
   ={}&\operatorname{Res}(F_n,T_{n+1},c-a)\,
       \operatorname{Res}(F_n,T_{n+1},\alpha)\,
       \operatorname{Res}(F_n,T_{n+1},c)^{n+1}.
\end{split}                                                \tag{2}
\]
For the first equality replace \(T_{n+1}\) by \(D F_{n+1}\),
split \(D\) from that slot and then replace \(T_{n+2}\)
in the factor with \(F_{n+1}\). For the second subtract
\(a T_{n+1}\) from \(T_{n+2}\), leaving
\(\alpha(c-a)c^{n+1}\). These operations can also be
seen from the multiplicity-weighted Poisson formula;
addition preserves each value on the complementary zero
scheme and products multiply the values. The polynomial
identities extend from finite generic intersections.
[Busé–D’Andrea–Sombra–Weimann, Propositions 2.1 and 2.2](https://www.maia.ub.edu/~sombra/publications/Flexes/Flexes.pdf)

Let \(\lambda_j=(2j+1)/(j+1)\). The five absolute
resultants needed in (2) are
\[
\begin{array}{c|c}
 \text{forms in }a,t,c &\text{absolute value}\\ \hline
 F_n,T_{n+1},c-a&2^{n(n+2)}(2n+1) C_n^{-n-1}\\
 F_n,T_{n+1},\alpha&2^{2n(n+1)}(2n+1) C_n^{-n-1}\\
 F_n,T_{n+1},c&2^{n(n+3)}C_n^{-n-2}\\
 F_n,F_{n+1},D&2^{n(n+1)}\lambda_n^{n+1}\lambda_{n+1}^n\\
 F_n,D,T_{n+2}&2^{n(n+3)-n}(n+1)^2\lambda_n^{n+3}.
\end{array}                                                \tag{3}
\]
For the first row set \(c=a\), giving
\(T_{n+1}=2(a-t)a^{n+1}\). At \(t=a\),
\(F_n=(2n+1)a^n\), and at \(a=0\),
\(F_n=2^n t^n/C_n\). For the second row set \(t=c\),
giving \(T_{n+1}=(a-c)a^{n+1}\). The specialization at
\(a=c\) is again \((2n+1)c^n\), and at \(a=0\) it is
\(4^n c^n/C_n\), from the beta integral of \(1-z^2\).
For the third row set \(c=0\), giving
\(T_{n+1}=(a-t)a^{n+1}\) and the two pure specializations
\(4^n a^n/C_n\) at \(t=a\) and \(2^n t^n/C_n\) at \(a=0\).

On \(D=0\), the quadratic interval polynomial is linear,
so \(F_j=\lambda_j G_j\) for
\(G_j=\sum_{i=0}^j a^i c^{j-i}\).
The binary resultant of \(G_j,G_{j+1}\) has absolute
value one from \(G_{j+1}-aG_j=c^{j+1}\).
Also \(T_{n+2}=(a-c)^2G_{n+1}/2\) on \(D=0\),
and \(G_n(a,a)=(n+1)a^n\).
Restriction along \(D=a+c-2t\) contributes
\(2^{de}\) to the resultant of the other two forms of
degrees \(d,e\), since the linear coefficient of \(t\)
is \(-2\). These identities give the remaining two rows.

Returning to \(b\) divides the absolute ternary resultant
by \(2^{n(n+1)(n+2)}\), from the determinant of the
coordinate change \(b=2t\). Multiplying the three forms
by \(C_n,C_{n+1},C_{n+2}\) contributes the powers
\((n+1)(n+2), n(n+2), n(n+1)\), respectively.
Putting (3) in (2) and using \(C_{j+1}=2E_j\) gives (1).

## Integral charts at the other primes

Choose a scalar \(w\) of valuation \(\omega_I(p)\) in
a finite extension of \(\mathbb Q_p\). In
\[
   \frac{\Lambda_2((Q+wW)^j)}{(j!)^2}
      =\sum_{r=0}^j
          \frac{(j-r)!}{j!\,r!}(wW)^r L_{j-r},
\]
every coefficient is integral and the leading forms at \(W=0\)
are exactly \(L_j\). In the chart \(W=1/w\),
the equations are the original normalized moments for \(Q+1\).
If any coefficient of \(Q\) had valuation less than
\(-\omega_I(p)\), scaling all projective coordinates to
have minimum valuation zero would give a point at \(W=0\)
in the reduced triple. This is impossible at a prime
outside (1), proving the stated bound. It is a weighted
version of the
[finite-free support deletion criterion](FACTORIAL_SUPPORT_DELETION_FINITE_FREE_MOMENTS.md).

The [exact replay](verify_binary_quadratic_primitive_homogeneous_resultant.py)
compares (1) with integer Macaulay determinants. The
unbounded proof is the factorization (2)–(3).
