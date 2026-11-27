# The binary triple-line scalar in a quintic Keller map

**Theorem QuinticKellerBinaryTripleLineScalarExclusion.** Over a field of
characteristic zero there is no three-variable Keller map \(F=(P,Q,R)\)
of ordinary degree at most five with
\[
                       R=x+x^2(z+xy)
\]
and independent binary fifth-degree parts \(P_5,Q_5\in\mathbb C[x,y]\)
after coefficient-field embedding. In particular this excludes the
[primitive degree-one residual](QUINTIC_KELLER_PURE_QUARTIC_AND_BINARY_TRIPLE_LINE.md#degree-four-exceptional-scalars).

Embed the finitely generated coefficient field in \(\mathbb C\).
Put \(v=z+xy\), write \(P,Q\) in the resulting coordinates, and
write \(r=x+x^2v\), \(B=\mathbb C[x,v]\). The Keller equation,
up to the constant choice of sign, is
\[
       [P,Q]=D(P)Q_y-P_yD(Q)=\kappa,\quad
       D=-x^2\partial_x+(1+2xv)\partial_v,\qquad\kappa\ne0 .
 \tag{1}
\]
The derivation fixes \(r\) and acts as \(-x^2\partial_x|_r\).
Thus its rational constants in \(\mathbb C(x,v)=\mathbb C(r,x)\)
are \(\mathbb C(r)\), and \(B[x^{-1}]=\mathbb C[r,x,x^{-1}]\).
For any nonzero univariate \(a(x)\),
\[
                   B[a^{-1}]\cap\mathbb C(r)=\mathbb C[r].
 \tag{2}
\]
A denominator root \(c\) in a reduced rational fraction in \(r\)
would force \(r-c\) to divide a power of \(a\), by substituting a
Bézout identity for its numerator and denominator. This is
impossible since \(r-c\) has positive \(v\)-degree.
There is no polynomial \(D\)-primitive of a nonzero constant:
such a primitive in \(\mathbb C[r,x,x^{-1}]\) would have the
form \(c/x+H(r)\), with \(c\ne0\) and \(H\) polynomial.

## Free-variable degrees

Retain all original ordinary degree bounds before the shear.
If a mate written in \((x,y,v)\) has \(y\)-degree \(m\), its
leading \(y\)-coefficient belongs to the space below, where
\(\mathbb C[x]_{\le j}\) denotes polynomials of degree at most \(j\):
\[
\begin{array}{c|l}
 m &\text{space}\\ \hline
 5&\mathbb C\\
 4&\mathbb C[x]_{\le4}\\
 3&\mathbb C[x]_{\le4}+\mathbb C v\\
 2&\mathbb C[x]_{\le4}+\mathbb C[x]_{\le2}v\\
 1&\mathbb C[x]_{\le4}+\mathbb C[x]_{\le2}v+\mathbb C v^2 .
\end{array}
 \tag{3}
\]
To check the bounds, group the original monomials by degree \(d\)
in \(y,z\). The pure \(y^d\) coefficient has \(x\)-degree at most
\(5-d\), and every coefficient with \(z\) has \(x\)-degree at
most \(4-d\). Set \(w=v/y\) in each group, so its expression
divided by \(y^d\) is \(A(x)+\sum_{k\ge1}a_k(x)(w-x)^k\).
If \(d>m\), the coefficients of \(w^j\) for \(j<d-m\)
are zero. For \(d=4\) the \(a_k\) are constant, so vanishing
at \(w=0\) makes \(a_k=0\) for \(k\ge2\), and vanishing also
of the coefficient of \(w\) makes the group zero.
For \(d=3\), vanishing of the constant coefficient forces
\(a_3\) constant by comparison at \(x^4\). If the \(w\)-coefficient
also vanishes, the \(x^3\)-coefficient of the constant term
and the \(x^2\)-coefficient of the \(w\)-term respectively give
\([x]a_2-a_3=0\) and \(-2[x]a_2+3a_3=0\). Thus \(a_3=0\)
and \(a_2\) is constant. For \(d=2\), vanishing of the constant
coefficient forces \(\deg a_2\le1\), so the coefficient of \(w\)
has degree at most two. The groups with \(d\le1\) and the pure
\(y^5\) term satisfy (3) directly. The bounds and the determinant
coefficient calculations below are also replayed by the
[exact checker](verify_quintic_binary_triple_line_scalar.py).

We use the degree-independent formal residue argument of
[Wang and Yan, Proposition 3.17](https://arxiv.org/pdf/2610.36495v1).
A mate independent of \(y\) would make (1) a product of a
\(D\)-derivative and a \(y\)-derivative equal to a unit, impossible.
A nonzero constant leading coefficient in \(y\) is impossible
as well. If \(P\) has degree \(m>0\) in \(y\), form
\(S=P^{-1/m}\) by expansion in \(y^{-1}\) whenever the chosen
\(m\)-th root of its leading coefficient is invertible in the
coefficient ring. The formal chain rule gives
\[
                 D(QP_y S)-\partial_y(QD(P)S)=-\kappa S .
 \tag{4}
\]
The second term has zero \(y^{-1}\)-coefficient, so the coefficient
of \(y^{-1}\) on the right must be a \(D\)-derivative in the
coefficient ring. A constant leading coefficient would give a
nonzero constant \(D\)-primitive in \(B\).
Two affine mates \(P=py+A,\ Q=qy+C\) are also impossible.
Their leading equation gives \(p/q\in\mathbb C(r)\).
The coefficients \(p,q\) cannot vanish simultaneously because
of (1), so are coprime. Substituting a reduced fraction for
\(p/q\) shows that both belong to \(\mathbb C[r]\). The constant
bracket would be \(D(qA-pC)=\kappa\).

Write \(m\ge n\ge1\) for the two \(y\)-degrees, with nonzero
leading coefficients \(p,q\). Thus \(m\le4\) by (3) and (4).
The leading equation of (1) gives
\[
                             p^n/q^m\in\mathbb C(r).
 \tag{5}
\]
For a reduced fraction in \(r\), a root of its numerator or
denominator respectively forces the corresponding coefficient
\(p\) or \(q\) to be divisible by \(r-c\).
Here every \(r-c\) is squarefree because \(r\) is a submersion.
In (3), no coefficient for \(m\ge3\) can be divisible by a full
\(r-c\), since this requires a nonzero \(v\)-coefficient
divisible by \(x^2\) in \(v\)-degree one. In either space
for \(m=1,2\), divisibility by \(r-c\) makes the coefficient
a scalar multiple of \(r-c\). A nonzero \(v^2\)-coefficient
would be constant and cannot be divisible by \(x^2\), and a
\(v\)-coefficient divisible by \(x^2\) can only be a scalar
multiple of \(x^2\). In these two spaces the \(D\)-invariant
coefficients are linear polynomials in \(r\), by (2) and (3).

If \(m=n\), then \(p/q\in\mathbb C(r)\) because
\(\mathbb C(r)\) is relatively algebraically closed in
\(\mathbb C(r,x)\). The affine pair was already excluded.
For \(m=n\ge3\) the ratio is constant by the absence of fibre
divisors. For \(m=n=2\) a nonconstant ratio would give a
full fibre factor in one coefficient and hence two linear
polynomial coefficients in \(r\). Their nonproportionality
would provide a constant linear combination with nonzero
constant leading coefficient, contradicting (4).
Thus in every case one can subtract a constant multiple of
one mate from the other to decrease one \(y\)-degree.
This preserves the original degree bound and preserves
independence and binary support of the quintic parts, since
it is an invertible constant target operation. Neither new
mate can be independent of \(y\). After reordering the degrees
are unequal, so the cancellation has terminated.

For \(m>n\) and \(m\ge3\), a full fibre factor in \(q\) would
make \(q\) a polynomial in \(r\). Equation (5) would then
force \(p\in\mathbb C[r]\), hence constant by (3).
Thus (5) is constant. Unique factorization gives
\(p=\lambda a^{m/d}\), \(q=\mu a^{n/d}\), with
\(d=\gcd(m,n)\) and nonzero scalars \(\lambda,\mu\).
Here \(a\) is univariate in \(x\): this follows from \(p\)
univariate when \(m=4\), and from its \(v\)-degree at most
one while \(m/d=3\) when \(m=3\).
In particular \(q(x)\) has degree at most two if \(n=2\),
and \(a\) is linear if \((m,n)=(4,3)\).

## The cases with smaller degree one

First take \(m=2,n=1\), with
\(P=py^2+hy+A,\ Q=qy+C\). If \(p/q^2\) is constant,
subtract a constant multiple of \(Q^2\) from \(P\) to get
the forbidden affine pair, for which the ordinary degree
bound is unnecessary. Otherwise (5) forces a full fibre
factor in \(q\) or in \(p\). In the second case \(p\) is a
multiple of that factor by (3), so (5) forces \(q\) invariant
and it is linear in \(r\). In both cases some \(r-c\) divides \(q\).
Reducing (1) modulo \(r-c\) gives
\(-\overline{P_y}\,\overline{DC}=\kappa\).
The factor \(\overline{DC}\) is a unit, by evaluating at \(y=0\).
Hence \(\overline p=0\) and \(p=\lambda q\), with \(q\)
a scalar multiple of \(r-c\).
The coefficient of \(y\) in (1) gives \(h=2\lambda C+L(r)\).
The constant coefficient then gives
\(D(q A-\lambda C^2-L C)=\kappa\), a contradiction.

For \(m\ge3,n=1\), \(Q=q(x)y+C\) with \(q\) nonconstant.
Use \(Q\) as variable over \(B[q^{-1}]\).
Equation (1) becomes \(qD(P)|_Q=\kappa\).
The positive \(Q\)-coefficients of \(P\) are invariant, thus
belong to \(\mathbb C[r]\) by (2). Subtract the polynomial
in \(Q,r\) given by those coefficients. The resulting
\(T\in B\) satisfies \(q DT=\kappa\), impossible.

## The cases with smaller degree two

Let \(m=3\) or \(4\), and write \(Q=q(x)y^2+By+C_0\) in this
paragraph, where the symbol \(B\) in the formula denotes its
linear coefficient. Put \(Y=y+B/(2q)\) and \(K=C_0-B^2/(4q)\).
Thus \(Q=qY^2+K\). Take the odd part \(fY+gY^3\) of \(P\) in
\(Y\) over \(\mathbb C[x,v,q^{-1}]\).
The even part commutes with \(Q\) because (1) is even and
its bracket with an even polynomial is odd. It is a
polynomial in \(Q\) over this localization, and commuting
with \(Q\) says that its coefficients are \(D\)-invariant.
By (2) they belong to \(\mathbb C[r]\).
In particular the odd part of \(P\) is polynomial in the
original free variable \(y\) over \(\mathbb C[x,v]\).
Also \(g\) is polynomial: for \(m=3\) it is \(p\), and
for \(m=4\) it is \(P_{[y^3]}-2\lambda q B\) since
\(p=\lambda q^2\) in that case.

The constant and quadratic \(Y\)-coefficients in (1) give
\[
             -f DK=\kappa,\qquad
             2q Df-fDq-3g DK=0 .
 \tag{6}
\]
In the localization, \(f\) and \(DK\) are units, so both
are univariate in \(x\). Substituting \(v=(r-x)/x^2\)
and comparing powers of \(r\) in \(DK\) gives
\(K=F(x)+G(r)\) with \(F\) rational and \(G\) polynomial.
If \(q(0)\ne0\), then \(F\) is regular at zero and
\(DK=-x^2F'\) vanishes there while \(f\) is a unit there.
This contradicts (6).

Let \(e=\operatorname{ord}_{x=0}q\in\{1,2\}\) and
\(k=\operatorname{ord}_{x=0}f\).
If \(F\) were regular there then \(k\le-2\) by (6).
But \(g=-f(2q Df-fDq)/(3\kappa)\) would have order
\(e+2k+1<0\), since \(2k-e\ne0\). Thus \(K\) has a pole
at \(x=0\). The expression \(K=C_0-B^2/(4q)\)
makes its pole order exactly \(e\) and makes
\(B_0=B(0,v)\) nonzero, since \(e\le2\).
It follows from \(DK=-x^2F'\) that \(k=e-1\).

If \(e=2\), then \(q=q_1x^2\) since \(\deg q\le2\).
The unit \(f\) is \(f_0 x\), so (6) gives \(g=0\).
But \(fY\) has the nonzero simple pole of
\(f B/(2q)\) at \(y=0\), impossible for the polynomial odd part.
If \(e=1\), write \(q=q_1 x+O(x^2)\) and
\(f=f_0+O(x)\). Comparing the pole term of \(K\) in (6)
gives \(B_0^2=4q_1\kappa/f_0\). In particular \(B_0\)
is a nonzero constant. The second equation of (6) gives
\(g=-f_0^2q_1x^2/(3\kappa)+O(x^3)\).
The coefficient of \(x^{-1}\) in the odd part at \(y=0\) is
\[
                  \frac{f_0B_0}{2q_1}
                     -\frac{f_0^2B_0^3}{24\kappa q_1^2}
                =\frac{f_0B_0}{3q_1}\ne0 .
\]
This is again impossible.

## The quartic and cubic free-variable degrees

Only \((m,n)=(4,3)\) remains. Its leading coefficients are
\(p=\lambda a^4,\ q=\mu a^3\) for linear \(a(x)\).
Use (4) with the roles exchanged and \(S=Q^{-1/3}\)
in \(\mathbb C[x,v,a^{-1}]((y^{-1}))\).
Its leading term is \(\mu^{-1/3}a^{-1}y^{-1}\).
Thus a rational function in \(\mathbb C(r,x)\) has
\(D\)-derivative equal to a nonzero scalar times \(1/a(x)\).
Its derivative in \(x\) at fixed \(r\) is a nonzero scalar
times \(1/(x^2a(x))\), which cannot have a simple pole
at a nonzero root of \(a\), since it is a rational derivative.
Consequently \(a\) is proportional to \(x\).

Rescale the two mates so the leading coefficients are \(x^4\)
and \(x^3\). If \(B_0\) denotes the coefficient of \(y^2\)
in \(Q\), set \(s=xy+B_0/(3x^2)\).
In \(\mathbb C[r,x,x^{-1}][s]\) write
\[
              Q=s^3+c s+d,\qquad
              P=s^4+\alpha s^3+\beta s^2+\gamma s+\delta .
\]
The derivative \(\partial_y s=x\) changes (1) to
\(D(P)Q_s-P_sD(Q)=\kappa/x\), with \(D\) at fixed \(s\)
and the nonzero constant \(\kappa\) rescaled as needed.
Equating powers five through two of \(s\) gives
\[
 \alpha\in\mathbb C[r],\quad
 \beta=4c/3+b,\quad
 \gamma=4d/3+\alpha c+e,\quad
 \delta=\alpha d+2c^2/9+2bc/3+f
\]
with \(b,e,f\in\mathbb C[r]\).
Put \(U=c+3b/2\) and \(V=d+3e/4\).
The last two coefficients give
\[
  UV\in\mathbb C[r],\qquad
  \frac4{27}U^3-\frac b3 U^2-\frac23 V^2
                       =\frac{\kappa}{2x^2}+H(r).
 \tag{7}
\]
Indeed the linear coefficient is \(-4D(UV)/3\) and
the constant is the \(D\)-derivative of the displayed
polynomial in \(U,V\).
Work in \(\mathbb C(r)[x,x^{-1}]\). If \(UV\ne0\), the two
factors are Laurent monomials of opposite exponents.
For \(U\) of positive exponent the highest power comes
from \(U^3\), and for \(U\) of negative exponent it comes
from \(V^2\). Neither is compatible with the right side.
Exponent zero in both also contradicts (7).
Hence \(UV=0\). If \(V=0\) and \(U\) is not constant
in \(x\), a positive highest exponent of \(U\) or a negative
lowest exponent of \(U\) gives the corresponding
uncancelled extreme exponent multiplied by three in (7).
Neither can match the exponents on the right. A constant
\(U\) is likewise impossible.
Thus \(U=0\), and (7) gives \(V=\eta/x\) for
a nonzero scalar \(\eta\), by comparing the extreme
powers and the coefficient of \(x^{-1}\).

It follows that
\[
                 Q=s^3-\tfrac32 b(r) s-\tfrac34 e(r)+\eta/x .
\]
Set the original free variable \(y=0\), so
\(s=B_0/(3x^2)\). If this is regular at \(x=0\) over
\(\mathbb C(v)\), the simple pole in \(Q\) persists.
If it has a pole, its cube has strictly higher pole order
than all other terms, since \(b(r),e(r)\) are regular.
Either case contradicts polynomiality of \(Q\).
This finishes the exclusion.
