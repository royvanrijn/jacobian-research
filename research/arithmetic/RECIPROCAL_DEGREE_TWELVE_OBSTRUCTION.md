# Reciprocal degree-twelve pencils: rational and real obstructions

There is no rational degree-twelve odd Laurent pencil with the good
reciprocal pairing of all its critical points described below. There is
also no real such pencil with every critical point real. These are
obstructions to the reciprocal ansatz; they do not exclude general
degree-twelve pencils with six paired simple critical values.

The rational obstruction is a projective modulo-31 certificate for one
explicit degree-twelve eliminant. It requires no integrality or
good-reduction hypothesis. The real obstruction follows from the same
eliminant and a uniform weighted-Rolle restriction.

## Setting and the even critical-degree operators

Let

\[
 f(z)=\frac{N(z^2)}z,\qquad N\in k[u]\text{ monic of degree }2k,
 \qquad N(0)\ne0,\qquad L=4k-1,
\]

over a field of characteristic zero. Put \(R_N=2uN'-N\). Assume:

* \(R_N\) has \(2k\) distinct nonzero roots.
* The values \(N(u)^2/u\) at those roots consist of exactly \(k\) distinct
  nonzero values, each taken twice.
* Their unique equal-value pairs are \(u\leftrightarrow\lambda/u\),
  for one \(\lambda\ne0\), and this involution has no fixed root.

Call this the good reciprocal locus. These hypotheses give \(2k\)
distinct finite branch values of \(f\), each with two simple critical
points; infinity has type \((4k-1,1)\).

For rational or real coefficients, the parameter \(\lambda\) is in the
same field: Galois preserves each unique equal-value pair, so it fixes
their common product. Choose \(\rho^2=\lambda\), set

\[
 F(v)=\rho^{-2k}N(\rho v),\qquad
 \Sigma(v)=\frac{R_N(\rho v)}{L\rho^{2k}},\qquad w=v+v^{-1}.
\]

Then \(\Sigma\) is monic reciprocal of degree \(2k\), with constant one,
and

\[
 2vF'-F=L\Sigma,\qquad F(0)=-L.                    \tag{1}
\]

At a paired root, \(F(v)=\pm vF(1/v)\). Define its plus factor by
\(F(v)+vF(1/v)=0\), and its minus factor by
\(F(v)-vF(1/v)=0\). Each factor is reciprocal with no fixed root, hence

\[
 \Sigma=v^k P(w)Q(w),\qquad \deg P=p,\quad\deg Q=q,\quad p+q=k,       \tag{2}
\]

with \(P,Q\) monic. The degree of each factor in \(v\) is even.

Define \(C_0=D_0=1\), \(C_1=w-1\), \(D_1=w+1\), with recurrence
\(Z_{j+1}=wZ_j-Z_{j-1}\). The Laurent-symmetric expressions

\[
 X(w)=\frac{F(v)+vF(1/v)}{v+1},\qquad
 Y(w)=\frac{F(v)-vF(1/v)}{v-1}
\]

are monic polynomials of degree \(2k-1\). They satisfy

\[
 \begin{split}
 X+2(w+2)X'&=L D_{k-1} P Q,\\
 Y+2(w-2)Y'&=L C_{k-1} P Q.                       \tag{3}
 \end{split}
\]

Indeed, differentiating their definitions gives
\(2R_F=(v-1)(X+2(w+2)X')+(v+1)(Y+2(w-2)Y')\).
The identity
\((v-1)D_{k-1}(w)+(v+1)C_{k-1}(w)=2v^k\)
then proves (3), or one may compare Laurent coefficients. The two
operators in (3) are polynomial isomorphisms: expand around \(w=-2\)
or \(w=2\), respectively, and divide the coefficient of degree \(j\)
by \(2j+1\).

At a root of \(P\), both \(X\) and its operator in (3) vanish.
The root is different from \(-2\), so \(X'\) vanishes there too.
Thus \(P^2\mid X\); similarly \(Q^2\mid Y\). Consequently \(p,q\ge1\):
a uniform-sign factor would have square degree \(2k>2k-1\).

If \(\lambda\) is nonsquare over the coefficient field, its quadratic
conjugation sends \(v\) to \(-v\) and interchanges the two sign factors.
They must have equal degrees in \(v\), namely \(k\); since each is even,
\(k\) must be even. In particular, when \(k\) is odd, \(\lambda\) is
necessarily a square. This gives a uniform obstruction in degrees
\(4k\equiv4\pmod8\).

For completeness, the polynomial conditions (3), \(P\mid X\),
\(Q\mid Y\) are also sufficient on the explicit open set: \(PQ\)
squarefree, \(PQ(2)PQ(-2)\ne0\), \(\gcd(P,Y)=\gcd(Q,X)=1\), and the
following \(k\) squared critical values distinct and nonzero. In
\(k[w]/(PQ)\), those values are

\[
 \frac{(w-2)Y^2}{4}\pmod P,\qquad
 \frac{(w+2)X^2}{4}\pmod Q.                       \tag{4}
\]

Equation (1) reconstructs \(F\) uniquely. Formula (4) follows directly from the definitions and
\(w-2=(v-1)^2/v\), \(w+2=(v+1)^2/v\), using the vanishing sign factor.
The conditions ensure exactly the required pairing and no additional
critical or branch degeneracies.

## Lemma 1. Uniform real balance in degree \(4k\)

If all the roots of \(P Q\) are real and larger than 2, then

\[
                  p\ge q,\qquad q\ge p-1.        \tag{5}
\]

Thus \(p=q\) when \(k\) is even and \(p=q+1\) when \(k\) is odd.
If all roots are smaller than \(-2\), interchange \(p,q\) in (5).

To prove this, on \(w>2\) put \(Z=\sqrt{w-2}\,Y\). Then \(Z(2)=0\),
each of the \(q\) roots of \(Q\) is a zero of \(Z\) of multiplicity at
least two, and

\[
 Z'=\frac{L C_{k-1}(w)P(w)Q(w)}{2\sqrt{w-2}}.
\]

The roots of \(C_{k-1}\) lie strictly inside \((-2,2)\); explicitly they
are \(2\cos((2j-1)\pi/(2k-1))\), \(1\le j\le k-1\).
Rolle's theorem requires a root of \(P\) before the first root of \(Q\)
and between every consecutive pair of roots of \(Q\). This proves
\(p\ge q\). Applying Rolle to \(\sqrt{w+2}\,X\) between consecutive
roots of \(P\) gives \(q\ge p-1\), since \(D_{k-1}\) also has all its
roots strictly inside \((-2,2)\). Replacing \(v\) by \(-v\) interchanges
the sign factors and sends \(w\) to \(-w\), proving the other assertion.
Only distinct roots on the good locus are used.

## Theorem 2. No rational good reciprocal model in degree twelve

Set \(k=3\), so \(L=11\). The preceding parity argument makes
\(\rho\) rational. Changing its sign interchanges \(p,q\), so we may
take \(p=1,q=2\). Write

\[
 P=w-c,\qquad Q=w^2+q_1w+q_0,\qquad Y=Q^2(w+a).
\]

The second equation in (3) becomes

\[
 4(w-2)(w+a)Q'+(3w+a-4)Q=11(w^2-w-1)(w-c).       \tag{6}
\]

Its two leading compatibility equations, with constant divisors 7 and 3,
give

\[
 Q=\frac{
 21w^2+(-27a-33c+27)w+
 45a^2+55ac-41a-55c+31}{21}.
\]

The remaining constant coefficient and \(X(c)=0\), after dividing its
automatic factor \(c+2\), give respectively

\[
 \begin{split}
 e_0={}&45a^3+55a^2c-5a^2-11ac-21a-11c-124,\\
 e_1={}&-1485a^2c^2+4455a^2c-2475a^2
 -1320ac^3+5313ac^2-4411ac-3289a\\
 &+360c^4+240c^3-3207c^2-107c-6241.
 \end{split}                                                   \tag{7}
\]

The removed factor \(c+2\) is nonzero on the good locus. No other
variable factor is divided out. Equations (7), with the opens above,
are an exact description of the reciprocal models in this normalization.
In particular, no condition \(F(r)=0\) occurs: zero is not a branch value
in this even critical-degree setting.

Exact elimination gives

\[
       \operatorname{Res}_a(e_0,e_1)=-57153600\,E(c),             \tag{8}
\]

where

\[
 \begin{split}
 E(c)={}&1125c^{12}-225c^{11}-16335c^{10}-43362c^9
 +207009c^8-100551c^7+2204367c^6\\
 &-15624026c^5+48834247c^4-83352841c^3
 +67437755c^2-3843328c+2647888.
 \end{split}                                                   \tag{9}
\]

Modulo 31, its values at \(0,1,\ldots,30\) are

\[
 (23,2,1,3,18,13,15,1,25,8,3,13,7,14,26,18,
  11,15,10,11,12,14,11,19,14,5,15,5,13,19,8).
\]

None is zero, and the homogeneous value at infinity is
\(1125\equiv9\pmod{31}\). Thus \(E\) has no point on
\(\mathbb P^1(\mathbb F_{31})\) and no rational root. Equations
(8)--(9) exclude every rational good model, with no assumption on
denominators or reduction of that model.

## Theorem 3. No all-real good reciprocal model in degree twelve

Suppose all twelve critical points of \(f\) are real. All six roots
of \(R_N\) are then positive, and hence \(\lambda>0\). Choose either
real square root \(\rho\), and normalize to \(p=1,q=2\) as above.

If \(\rho>0\), all roots of \(P Q\) exceed 2, contradicting Lemma 1.
If \(\rho<0\), all those roots are below \(-2\), so in particular
\(c<-2\). But \(E(c)>0\) for every \(c\le0\). To see this, set \(c=-x\),
\(x\ge0\). All coefficients of

\[
 E(-x)-x^8(1125x^4-16335x^2+207009)
\]

are nonnegative and its constant coefficient is positive. The remaining
quartic is positive because

\[
                 4(1125)(207009)-16335^2=664708275>0.
\]

This contradicts (8) and completes the proof over \(\mathbb R\).
The conclusion covers algebraic as well as transcendental real
coefficients.

As an independent exact control, a complete rational lexicographic
Groebner basis of (7) has one equation \(a=A(c)\) and one \(E(c)=0\).
Sturm isolation gives exactly two real candidate values,
\(c\approx2.924978061956\) and \(3.362643846155\). Exact interval
evaluation gives \(Q(2)\in(-10,-9)\) and \(Q(2)\in(-6,-5)\),
respectively. These controls are retained in the receipt; the preceding
short proof does not depend on numerical root approximations.

## Scope and replay

The all-real theorem restricts the reciprocal critical-point pairing,
rather than arbitrary real paired-critical Laurent covers. Its uniform
balance lemma also applies in higher degrees \(4k\), but does not prove
a general reciprocal obstruction there. It does not address the rational
degree-26 systems in the
[mixed coupled-ODE note](MIXED_RECIPROCAL_COUPLED_ODE.md).

The reduction is the even-critical-degree counterpart of that note and
the [reciprocal obstruction note](PAIRED_CRITICAL_RECIPROCAL_OBSTRUCTION.md).
All arguments needed here are included above. No novelty comparison or
external existence theorem is being asserted.

Replay from the repository root:

~~~text
.venv/bin/python research/arithmetic/verify_reciprocal_degree_twelve.py --write-receipt research/arithmetic/reciprocal_degree_twelve_verification.json
~~~

The checker binds (6)--(9) to the actual coefficients of \(F\), verifies
both symmetric operator identities, the exact resultant, the complete
modulo-31 table, the positive-polynomial certificate, and the independent
lexicographic/Sturm/interval controls. Its bounded general-operator
controls supplement the written uniform proof in Lemma 1.
