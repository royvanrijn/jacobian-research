# Univariate quadratic Strong Factorial: resonance and uniform prime-degree gates

Let \(K\) be a field of characteristic zero, and define \(\Lambda(x^j)=j!\). For a polynomial \(f\), let \(\mathcal N(f)\) be its number of nonzero monomials. Its Strong Factorial assertion asks for a nonzero moment in every positive block of length \(\mathcal N(f)\). In \(SFC(m)\), \(m\) counts variables, not this window length.

The phases here are dense univariate quadratics
\[
 f=ax^2+bx+c,\qquad abc\ne0. \tag{1}
\]
Their window has length three. They form a restricted family in \(K[x]\), also embedded in \(K[x,y]\). An arbitrary binary polynomial of total degree at most two may have six monomials and is a different target.

For an odd integer \(q\ge3\) with \(s_2(q-1)\ge4\), let \(r_4(q)\) be the fourth-lowest set-bit position of \(q-1\), positions numbered from zero, and put \(D(q)=2^{r_4(q)}\). An odd prime \(p\mid q\) is **admissible** if, with \(t=v_p(q)\), every base-\(p\) digit of \(q\) is at most \((p-1)/2\) and
\[
 L_p(q)<B_p(q),\qquad
 p^{L_p(q)}\le2q<p^{L_p(q)+1},\qquad
 B_p(q)=tp^t-\frac{p^t-1}{p-1}.
\]
These are integer tests; a floating-point logarithm is not involved. Let
\(J(q)=\prod_{p\mid q,\ p\text{ admissible}}p^{v_p(q)}\), with the empty product equal to 1, and define
\[
 \mathcal R=\{q\ge3:q\text{ odd},\ s_2(q-1)\ge4,
                  \ D(q)J(q)\le q-2\}. \tag{2}
\]

**Theorem 1 (partial Strong Factorial criterion).** If (1) has three consecutive zero factorial moments starting at \(n\ge1\), then
\[
 -c/b=n+2\in\mathcal R. \tag{3}
\]
Thus every dense phase with \(-c/b\notin\mathcal R\) satisfies its support-length-three window at every positive start. For a fixed \(q\in\mathcal R\), only \(n=q-2\) remains possible. The first integer not excluded by this criterion is \(q=91,n=89\). This is an uncovered index, not an example of vanishing.

Proposition 7 gives an exact reversible bridge between the remaining rational-polynomial coprimeness problem and the universal dense quadratic Strong Factorial assertion. Its two directions are proved with the local gates as contexts; neither open assertion is assumed as a proof input. The theorem is partial: it does not settle all \(q\in\mathcal R\), arbitrary inhomogeneous binary quadratic Strong Factorial, full \(SFC(2)\), or \(SFC(3)\).

## 1. The recurrence and the unique possible starting index

Put \(M_j=\Lambda(f^j)\), with \(M_0=1\). For every \(j\ge2\),
\[
 M_j=2aj(2j-1)M_{j-1}
   +j(j-1)(b^2-4ac)M_{j-2}
   +c^j+bj c^{j-1}. \tag{4}
\]
The boundary forcing has a plus sign. This is important: for \(f=x+1\), one has \(M_1=2,M_2=5\).

For completeness, over \(\mathbb C\) factorial integration gives
\(M_j=\int_0^\infty f(x)^j e^{-x}\,dx\). Integrate the derivative of \(f'f^r e^{-x}\), whose boundary value is \(-bc^r\), and use
\((f')^2=4af+(b^2-4ac)\). This gives
\[
 \int_0^\infty f'f^r e^{-x}\,dx
 =2a(2r+1)M_r+r(b^2-4ac)M_{r-1}+bc^r.
\]
Integration of the derivative of \(f^{r+1}e^{-x}\) also gives
\(M_{r+1}=(r+1)\int f'f^r e^{-x}\,dx+c^{r+1}\), proving (4). The resulting polynomial identity is over \(\mathbb Q\) and therefore holds over every characteristic-zero field.

If \(M_n=M_{n+1}=M_{n+2}=0\), equation (4) at \(j=n+2\) gives
\[
 c^{n+1}(c+(n+2)b)=0. \tag{5}
\]
For (1), this forces \(c=-(n+2)b\). Dividing the phase by \(b\), it remains to study
\[
 Ax^2+x-q,\qquad A=a/b,\qquad q=n+2\ge3. \tag{6}
\]
Thus every nonresonant coefficient ratio is already excluded at all starts. The next sections give the even-constant, uniform height-bounded prime and binary degree gates.

## 2. Even constants: a degree-preserved mod-2 resultant

**Lemma 2.** Fix an even integer \(C\), and put
\[
 U_m(T)=\frac{\Lambda((Tx^2/2+x+C)^m)}{m!}\in\mathbb Q[T].
\]
For every \(m\ge1\), the polynomials \(U_m,U_{m+1}\) are coprime over \(\mathbb Q\). Hence they have no common value-zero in any characteristic-zero field.

Let \(\mathbb Z_{(2)}\) denote the rationals with odd denominator. Formal exponential expansion, followed by coefficientwise factorial integration, gives
\[
 \sum_{m\ge0}U_m(T)z^m
 =e^{Cz}\sum_{i\ge0}(2i-1)!!\,T^i\frac{z^i}{(1-z)^{2i+1}}. \tag{7}
\]
There is no convergence assertion: this is an identity in \(\mathbb Q[T][[z]]\), with \((-1)!!=1\). The coefficient follows from
\((2i)!/(2^ii!)=(2i-1)!!\).

For \(r\ge1\), the evenness of \(C\) implies
\[
 v_2(C^r/r!)\ge r-v_2(r!)=s_2(r)\ge1,
\]
with zero terms interpreted as having infinite valuation. Thus \(e^{Cz}\) belongs to \(\mathbb Z_{(2)}[[z]]\) and reduces to 1 modulo 2. Every double factorial in (7) is odd, so \(U_m\in\mathbb Z_{(2)}[T]\), with odd leading coefficient \((2m-1)!!\). Its reduction \(p_m\) has generating series
\[
 \sum_{m\ge0}p_m(T)z^m
 =\frac{1+z}{1+Tz+z^2}\quad\text{over }\mathbb F_2. \tag{8}
\]
Therefore
\[
 p_0=1,\qquad p_1=1+T,\qquad p_m=Tp_{m-1}+p_{m-2}\ (m\ge2). \tag{9}
\]
The degree is exactly \(m\). A common root of \(p_m,p_{m+1}\) would propagate backward by (9) to \(p_0=0\), impossible. Their resultant is nonzero over \(\mathbb F_2\). Since both degrees are preserved, the rational resultant of \(U_m,U_{m+1}\) is a unit of \(\mathbb Z_{(2)}\), in particular nonzero. A rational polynomial Bezout identity persists over every characteristic-zero field. This proves Lemma 2.

In (6), if \(q\) is even, take \(C=-q\) and \(T=2A\). Even the first two zero equations are impossible. This closes every even starting index in the dense quadratic problem.

## 3. A uniform no-carry prime gate with a necessary height guard

For a fixed \(q\ge3\), write
\[
 R_h(A)=\frac{\Lambda((Ax^2+x-q)^h)}{h!},\qquad
 H_q(X)=\frac{X^qR_q(1/X)}{q!}
        =\sum_{j=0}^q a_j\frac{X^j}{j!}. \tag{10}
\]

**Lemma 3.** If \(p\mid q\) is admissible in the sense defined above, every rational factor of \(R_q\) has degree divisible by \(p^{v_p(q)}\).

Fix \(t=v_p(q)\ge1\), set \(L=L_p(q)\), and put \(B=tp^t-(p^t-1)/(p-1)\). Native factorial expansion gives
\[
 a_j=\sum_{k=0}^j\binom jk(-q)^k
                   \frac{(2q-j-k)!}{q!(q-j)!}. \tag{11}
\]
Let \(u_{j,k}\) be its \(k\)-th summand and \(N=2q-j\). The native term is the nonzero integer \(u_{j,0}=\binom Nq\); the exact ratio is
\[
 \frac{u_{j,k}}{u_{j,0}}
 =\binom jk\frac{(-q)^k}{(N)_k},\qquad
 (N)_k=N(N-1)\cdots(N-k+1). \tag{12}
\]
All factors are positive since \(k\le j\le q\le N\).

For \(1\le k\le N\), counting multiples of each power of \(p\) in an interval of length \(k\) gives
\[
 v_p((N)_k)\le v_p(k!)+\lfloor\log_pN\rfloor
             \le v_p(k!)+L. \tag{13}
\]
There are at most \(\lfloor k/p^e\rfloor+1\) multiples at the \(e\)-th level, and only \(\lfloor\log_pN\rfloor\) levels occur. For \(k\ge p^t\), Legendre's formula and \(s_p(k)\ge1\) imply
\[
 tk-v_p(k!)\ge k\left(t-\frac1{p-1}\right)+\frac1{p-1}
             \ge tp^t-\frac{p^t-1}{p-1}=B. \tag{14}
\]
Because \(L<B\), every large-shift ratio in (12) has valuation at least \(B-L\ge1\). Such terms are integral and strictly higher than the native term. No monotonicity of \(tk-v_p(k!)\) is being assumed.

For the small shifts \(1\le k<p^t\), retain every factorial level. For \(P=p^e\), define
\[
 \begin{split}
 \Phi_e={}&\lfloor j/P\rfloor-\lfloor k/P\rfloor
            -\lfloor(j-k)/P\rfloor\\
          &+\lfloor(2q-j-k)/P\rfloor-\lfloor q/P\rfloor
            -\lfloor(q-j)/P\rfloor .
 \end{split}
\]
Exactly,
\[
 v_p(u_{j,k})=tk+\sum_{e\ge1}\Phi_e. \tag{15}
\]
If \(e\le t\), then \(P\mid q\). With
\(\delta=\lfloor j/P\rfloor-\lfloor k/P\rfloor-\lfloor(j-k)/P\rfloor\in\{0,1\}\),
\(\Phi_e=\delta+\lceil j/P\rceil-\lceil(j+k)/P\rceil\ge-\lceil k/P\rceil\).
Allocate the term \(tk\) as \(+k\) at each of these \(t\) levels. Their contributions are nonnegative. If \(p^t\mid j\), the exact value is \(\Phi_e=-\lfloor k/P\rfloor\), and each allocated contribution \(k-\lfloor k/P\rfloor\) is strictly positive.

If \(e>t\), put \(u=q\bmod P\), \(v=j\bmod P\). The least nonzero digit of \(q\) is at position \(t\), and the digit bound gives
\(p^t\le u<P/2\). Thus \(k<p^t\le u\) and \(k<P\). Direct floor decomposition gives
\[
 \Phi_e=\mathbf1_{v<k}
      +\left\lfloor\frac{2u-v-k}{P}\right\rfloor
      +\mathbf1_{v>u}. \tag{16}
\]
For \(v\le u\), the floor is zero since \(0<u-k\le2u-v-k< P\). For \(v>u\), one has \(v>k\) and \(2u-v-k>-P\); the floor is at least \(-1\), canceled by the final indicator. Hence \(\Phi_e\ge0\) at every remaining level. Combining small and large shifts proves
\[
 a_j\in\mathbb Z_p\ (0\le j\le q),\qquad
 a_j\equiv\binom{2q-j}{q}\pmod p\quad\text{if }p^t\mid j. \tag{17}
\]

Take successive high-digit prefixes of \(q\), including 0 and \(q\). Each such index \(j\) is divisible by \(p^t\); the sum \(q+(q-j)=2q-j\) has no carry, since it doubles only the suffix digits. The factorial digit formula therefore makes the binomial coefficient in (17) a unit. In particular \(a_j\) is a unit at every such prefix, including the endpoints.

The complete lower polygon can be checked without an imported Coleman criterion. Let \(H\) be a high-digit prefix and let the next nonzero digit of \(q\) be \(d\) in position \(s\ge t\). The next prefix is \(H+dp^s\). The line between their factorial heights has slope
\[
 -\lambda_s=-\frac{p^s-1}{(p-1)p^s}. \tag{18}
\]
Every coefficient lies on or above this supporting line. Equivalently,
\[
 s_p(j)\ge s_p(H)+(j-H)/p^s,\qquad0\le j\le q. \tag{19}
\]
To prove it, write \(H=Qp^s\) with \(p\mid Q\), and \(j=vp^s+w\), \(0\le w<p^s\). Then \(v\le Q+d\). The integer function \(s_p(v)-v\) is nonincreasing, with increment \(-(p-1)v_p(v+1)\), and is constant from \(Q\) through \(Q+d\) since \(d<p\). Also \(s_p(w)-w/p^s\ge0\). These facts give (19). By Legendre's formula and the integrality of \(a_j\), this is precisely the supporting-line bound. Prefix endpoints attain their heights because their \(a_j\) are units. As the successive positions \(s\) decrease, the slopes in (18) increase, giving all lower edges.

A root valuation is the negative of one of these slopes; otherwise its vanishing sum has a unique term of least valuation. Thus each root of \(H_q\) has valuation
\(\lambda_s=(1+p+\cdots+p^{s-1})/p^s\), for some \(s\ge t\), with reduced denominator \(p^s\). A monic irreducible factor over \(\mathbb Q_p\) has conjugate roots of the same valuation. Its nonzero constant coefficient has integer valuation equal to its degree times \(\lambda_s\), so its degree is divisible by \(p^s\), hence by \(p^t\). The endpoint units ensure nonzero constant and full degree of \(H_q\); reciprocity preserves factor degrees. Summing local factor degrees, including multiplicities, proves Lemma 3.

The recurrence at index \(q\), divided by \(q!\), has zero boundary forcing. It yields
\[
 R_q=(1+4qA)R_{q-2}+2A(2q-1)R_{q-1}. \tag{20}
\]
Every rational common factor \(G\) of the two preceding moments divides \(R_q\). Applying Lemma 3 at every admissible prime simultaneously gives
\[
 J(q)\mid\deg G. \tag{21}
\]
Distinct prime powers in \(J(q)\) are coprime, so this is a product, not an assumption that one local factor is shared between different primes.

This subsumes the one-digit range: if \(q=up^t\) and \(2u<p\), then \(L=t\) and
\(B-t=(p^t-1)(t-1/(p-1))>0\). There is no carry, so the prime is admissible. In particular for every odd prime power \(q=p^t\), a common factor would have degree at least \(q\) but at most \(q-2\); hence
\[
 \gcd_{\mathbb Q[A]}(R_{q-2},R_{q-1})=1
             \qquad(q=p^t,\ p\text{ odd}). \tag{22}
\]
The height guard cannot simply be omitted from a no-carry deformation; Appendix A gives an explicit counterfamily to that stronger coefficient-integrality claim. It is not a counterexample to Strong Factorial.

## 4. Odd constants and at most three binary digits

**Lemma 4.** Fix any odd integer \(C\), and put
\(F_j(A)=\Lambda((Ax^2+x+C)^j)\).
If \(m\ge2\) is even and \(s_2(m)\le3\), then \(F_{m-1},F_m\) are coprime over \(\mathbb Q\), with the same all-characteristic-zero conclusion. More generally a common root would require
\[
 v_2(A)=-2+2^{-r},\qquad r\text{ a set-bit position of }m,
 \qquad s_2(m\bmod2^r)\ge3. \tag{23}
\]
Only the fourth or later set-bit cluster can therefore remain.

The literal coefficient expansion is
\[
 [A^i]F_j=\binom ji(2i)! H_{j,i}(C),
\quad H_{j,i}(C)=\sum_{t=0}^{j-i}\binom{j-i}t
       \frac{(2i+t)!}{(2i)!}C^{j-i-t}. \tag{24}
\]
For \(t\ge2\), the quotient contains the even factor \(2i+2\). The terms \(t=0,1\) give \(H_{j,i}\equiv1+j-i\pmod2\). Thus for even \(j\),
\[
 v_2([A^i]F_j)\ge2i+s_2(j-i)-s_2(j), \tag{25}
\]
with equality for every even \(i\), and at least one further unit of valuation for odd \(i\). This uses
\(v_2\binom ji=s_2(i)+s_2(j-i)-s_2(j)\) and
\(v_2((2i)!)=2i-s_2(i)\).
In particular \(F_j(0)\) is odd for even \(j\).

Let the set bits of \(m\) be \(r_1<\cdots<r_k\), all at least one, and let
\(I_h=\sum_{a=1}^h2^{r_a}\), \(I_0=0\). The exact lower polygon of the coefficient points of \(F_m\) has vertices
\[
 (I_h,2I_h-h),\qquad h=0,\ldots,k, \tag{26}
\]
with slopes \(2-2^{-r_h}\).
To verify all supporting lines, fix \(r=r_h\) and write
\(m=H+2^r+L\), where \(H\) is a multiple of \(2^{r+1}\) and \(L=I_{h-1}<2^r\). For \(t=m-i\), the lower bound (25) lies above the segment in (26) exactly when
\[
 s_2(t)-t/2^r\ge s_2(H)-H/2^r,\qquad0\le t\le m. \tag{27}
\]
Write \(t=q2^r+u\), \(0\le u<2^r\). Then the left side is
\(s_2(q)-q+s_2(u)-u/2^r\). Its last two terms are nonnegative. The integer function \(s_2(q)-q\) is nonincreasing, with increment \(-v_2(q+1)\). Also \(q\le H/2^r+1\), and its values at the last two indices coincide because \(H/2^r\) is even. This proves (27). All vertices in (26) attain (25) since their indices are even.

A root of \(F_m\) is algebraic and nonzero. Extend its valuation over \(\mathbb Q_2\). Unless its valuation is the negative of a polygon slope, a unique term in the vanishing sum has least valuation, impossible. Thus \(v_2(A)=-2+2^{-r}\) for a set bit \(r\) of \(m\).

Suppose now \(F_{m-1}(A)=F_m(A)=0\) and let \(s=v_2(m)\ge1\). Since \(v_2(A)>-2\), the factor \(1-4AC\) is a unit. The factors \(C+m\), \(C\), and \(m-1\) are also units. The recurrence (4) at \(j=m\) forces
\[
 v_2(F_{m-2}(A))=-s. \tag{28}
\]
For \(j=m-2\), however, (25) gives
\[
 v_2(F_j(A))\ge
 \min_i\{i/2^r+s_2(j-i)-s_2(j)\}
 =\frac{L'}{2^r}-s_2(L'),\qquad L'=j\bmod2^r. \tag{29}
\]
For the minimum formula, decompose \(j-i=q2^r+u\) as in (27). The same nonincreasing digit function makes the minimum occur at \(u=0\), \(q=\lfloor j/2^r\rfloor\); substituting \(j=Q2^r+L'\) proves (29). Cancellation can only increase the valuation, so this lower bound suffices.

If \(r=s\), then \(L'=2^s-2\), whose binary digit sum is \(s-1\). The bound is
\[
 v_2(F_{m-2}(A))\ge-s+2-2^{1-s}>-s. \tag{30}
\]
This includes \(m=2\), when \(F_0=1\).
If \(r>s\), set \(L=m\bmod2^r\) and \(t=s_2(L)\). Subtracting 2 from the lowest set bit gives \(s_2(L-2)=t+s-2\), so
\[
 v_2(F_{m-2}(A))\ge-s+(2-t)+(L-2)/2^r. \tag{31}
\]
For \(t=1\), the added quantity is positive. For \(t=2\), one has \(L\ge6\), again giving strict positivity. Both contradict (28). Hence a common root needs at least three lower set bits, proving (23). If \(m\) has at most three set bits, every cluster is excluded. Coprimality over \(\mathbb Q\) and then a rational Bezout identity prove Lemma 4 over every characteristic-zero field.

This is not a closure of the fourth-bit cluster. For \(m=30,r=4,s=1\), the bound (31) is \(-5/4\), below the forced value \(-1\); no contradiction follows. The bound does not assert that such a common root exists.

## 5. The general binary common-factor degree divisor

**Lemma 5.** Fix an odd integer \(C\) and even \(m\) with at least four set bits. Let \(r_4\) be its fourth-lowest set-bit position and \(D=2^{r_4}\). Every rational common factor of
\(F_{m-1}(A)=\Lambda((Ax^2+x+C)^{m-1})\) and \(F_m(A)\) has degree divisible by \(D\).

The common factor has no zero root since \(F_m(0)\) is odd. Factor it over \(\mathbb Q_2\) into monic irreducible components. By Lemma 4, each root has valuation
\[
 \rho=-2+2^{-r}=\frac{-2^{r+1}+1}{2^r},\qquad r\ge r_4. \tag{32}
\]
The numerator is odd. All conjugates of one local irreducible factor have the same valuation, because the valuation on an algebraic closure of the complete field is invariant under local Galois conjugacy. Its constant coefficient is nonzero and has integer valuation \((\deg h)\rho\), so \(2^r\mid\deg h\). In particular \(D\mid\deg h\). Summing the local factor degrees, with multiplicities, proves
\[
 D\mid\deg G. \tag{33}
\]
No identification of local factors at different primes is required when this is combined with (20): the total rational degree has both divisibilities.

## 6. Examples of the uniform criterion

The uniform criterion removes \(q=55,63,87\). Indeed:

| \(q,p\) | \(t\) | \(L_p(q)\) | \(B_p(q)\) | base-\(p\) digits of \(q\) | \(D(q)p^t\) |
|---|---:|---:|---:|---|---:|
| \(55,5\) | 1 | 2 | 4 | \(210_5\) | \(32\cdot5=160>53\) |
| \(63,7\) | 1 | 2 | 6 | \(120_7\) | \(16\cdot7=112>61\) |
| \(87,29\) | 1 | 1 | 28 | \(30_{29}\) | \(64\cdot29=1856>85\) |

Each row has the digit bound and strict height inequality. Lemmas 3 and 5 therefore give the contradiction directly, with no computational-certificate premise.

## 7. An exact reversible bridge for the remaining dense problem

**Proposition 7.** The following assertions are equivalent:

- Every dense univariate quadratic over every characteristic-zero field satisfies its support-length-three Strong Factorial window at every positive start.
- For every \(q\in\mathcal R\),
  \(\gcd_{\mathbb Q[A]}(R_{q-2}(A),R_{q-1}(A))=1\).

The first assertion can equivalently be imposed just over \(\mathbb C\). The converse from a single arbitrary field \(K\) is not asserted; a common algebraic root need not lie in that particular field.

To prove the forward direction, suppose a pair in the second assertion has a nonconstant rational common factor. It has a common root \(A\in\overline{\mathbb Q}\). This root is nonzero: at \(A=0\), its upper moment, of even order \(q-1\), is
\[
 \Lambda((x-q)^{q-1})
   =\int_0^\infty(x-q)^{q-1}e^{-x}\,dx>0. \tag{35}
\]
Thus \(f=Ax^2+x-q\) is dense. The first two moments at \(q-2,q-1\) vanish. At index \(q\), (4) has zero forcing, so its next moment also vanishes. This is a forbidden three-zero block for a dense phase over a number field, and also over \(\mathbb C\).

Conversely, a dense counterexample over any characteristic-zero field normalizes by (5)--(6) to a pair with integer \(q=n+2\). Theorem 1, proved independently in Section 8 from the local degree gates alone, places \(q\) in \(\mathcal R\). Rational coprimeness supplies a rational Bezout identity; it persists in the coefficient field and excludes that common value-zero. Therefore the second assertion implies the first. These are written implications with the local arithmetic results as contexts, not a bootstrap using the desired dense conjecture as a theorem.

## 8. The product criterion and its exact finite frontier

Suppose a dense quadratic has a three-zero block. Equation (5) gives \(q=-c/b=n+2\ge3\). If \(q\) is even, Lemma 2 excludes its first two zeros. If \(q\) is odd and \(s_2(q-1)\le3\), Lemma 4 excludes them. Otherwise \(D(q)\) is defined. A common value-zero is algebraic, since it is a root of a nonzero rational moment polynomial, and gives a nonconstant rational common factor \(G\). Lemma 5 and (21) imply
\[
 D(q)J(q)\mid\deg G\le q-2. \tag{34}
\]
The odd integer \(J(q)\) is coprime to the power of two \(D(q)\), justifying the product. Hence \(D(q)J(q)\le q-2\), exactly the remaining condition in (2). This proves Theorem 1 without Proposition 7 or the optional companion certificates as premises.

Odd prime powers are excluded by (22). Below 55, the only odd numbers with at least two distinct prime divisors are \(15,21,33,35,39,45,51\), whose decremented binary digit sums are \(3,2,1,2,3,3,3\). The complete odd composite list from 55 up to but not including 91 is
\[
 55,57,63,65,69,75,77,81,85,87. \tag{36}
\]
Their decremented digit sums are respectively \(4,3,5,1,2,3,3,2,3,4\). The uniform gate removes 55,63,87 as in Section 6; all other listed values are removed by the binary test, and 81 is also a prime power.

Finally \(91=7\cdot13\) and \(90=2+8+16+64\), so \(D(91)=64\). In base 7 one has \(91=160_7\), whose digit 6 exceeds 3. In base 13 one has \(91=70_{13}\), whose digit 7 exceeds 6. Neither prime is admissible. Thus \(J(91)=1\) and \(64\le89\), proving \(91\in\mathcal R\). The first start not excluded by the stated criterion is therefore \(n=89\). The finite trial-division/digit enumeration proves only that frontier statement; it asserts neither vanishing there nor the failure of every further method.

## 9. Other supports, references and logical boundaries

The statements above deliberately concern three nonzero coefficients. Some lower-support cases have elementary separate proofs: monomials never have a zero moment; an affine binomial satisfies \(U_j=jaU_{j-1}+c^j\), excluding adjacent zeros; and for \(ax^2+bx\) with \(ab\ne0\), the homogeneous recurrence (4) propagates a putative adjacent zero pair backward to \(M_0=0\), a contradiction. The constant-plus-square support has a previously proved Strong Factorial theorem; it is not a new consequence of the resonance lemmas here.

The standard meaning of the conjecture follows [Edo--van den Essen, *The Strong Factorial Conjecture*, arXiv:1304.3956v2, Definition 2.7 and Conjecture 2.8](https://arxiv.org/html/1304.3956v2). [Rocks, author preprint](https://people.math.sc.edu/rocks/BradyRocksSFC2015.pdf), particularly Theorems 3.10 and 3.24, provides the prime-power-variable-monomial and constant-plus-square context. Neither the general theorem nor any irreducibility theorem is a premise of the dense three-term proof above. All valuations, polygon supporting lines and degree gates needed here have been written explicitly. The earlier longer odd-prime-power initial-form argument is not a premise of the height-bounded uniform gate (Lemma 3).

Ordinary Factorial assertions exclude vanishing at every positive power or an infinite tail; they do not themselves establish the exact support-size windows. In particular [Liu--Sun, Theorem 2.6](https://doi.org/10.1017/S0004972719000546) concerns ordinary homogeneous binary Factorial, and the binary Factorial/EIT argument likewise has a different quantifier. No such implication is used above.

The accompanying finite exact checks validate the native factorial transcription and the explicitly finite first-index enumeration. The universal conclusions follow from the written arguments, not from a finite moment prefix. No explicit counterexample, full \(SFC(2)\), arbitrary inhomogeneous binary quadratic theorem, or \(SFC(3)\) conclusion is asserted. No novelty or external assurance claim is made.

The supplementary [portable standard-library checker](verify_univariate_quadratic_strong_factorial_resonance.py) and [retained finite verification receipt](univariate_quadratic_strong_factorial_resonance_verification.json) validate native transcriptions and the explicitly finite frontier arithmetic. From the directory containing this note, run a fresh replay without overwriting the retained receipt:

```text
python verify_univariate_quadratic_strong_factorial_resonance.py --output fresh_verification.json
```

These tests are supplementary; the universal assertions are proved in the written lemmas, not inferred from a checked moment prefix. The enormous counterfamily in Appendix A is handled by its all-range valuation proof, not by expanding a huge input.

## Appendix A. No carry alone does not preserve the actual coefficients

The height guard in Lemma 3 is a sufficient restriction on the actual deformed coefficients, not an inference from the undeformed polynomial. The following counterfamily shows why the guard matters. Fix an integer \(h\ge19\), put \(Q=5^h\), \(q=Q+10\), and \(j=Q-1=q-11\). The base-5 digits of \(q\) are 1 in position \(h\), 2 in position 1, and zero elsewhere. There is no doubling carry, but \(t=v_5(q)=1\), \(B_5(q)=4\), and \(L_5(q)=h\), outside the height bound.

For the actual coefficient \(a_j\) in (11), every \(\binom jk\), \(0\le k\le j\), is a 5-unit: the factorial complement identity for \(Q=5^h\) gives \(v_5\binom{Q-1}k=0\). In this case \(2q-j=Q+21\) and \(q-j=11\). For \(0\le k\le21\), the valuation of its \(k\)-th summand is exactly
\[
 k+v_5((21-k)!)-4\ge0. \tag{37}
\]
Indeed \(v_5((Q+r)!)=v_5(Q!)+v_5(r!)\) for \(0\le r\le21<Q\), and \(v_5(10!)=v_5(11!)=2\). The displayed bound follows directly for \(k=0,1,2,3\), and for \(k\ge4\) is immediate.

For \(k=22+r\), \(0\le r\le Q-23\), factorial complementation below \(Q\) gives the exact valuation
\[
 18-h+r-v_5(r!). \tag{38}
\]
Here \(v_5((Q-1-r)!)=v_5((Q-1)!)-v_5(r!)\) and \(v_5(q!)=v_5((Q-1)!)+h+2\). At \(r=0\), (38) equals \(18-h<0\). For \(r\ge1\), Legendre's formula gives \(r-v_5(r!)>0\), so every other tail term has strictly higher valuation. All terms in (37) have nonnegative valuation. A unique summand therefore has the least valuation, and
\[
 v_5(a_j)=18-h<0. \tag{39}
\]
This refutes the unrestricted no-carry coefficient-integrality assertion. It does not construct a common root of any moment pair or refute Strong Factorial. No enormous finite computation is needed for (37)--(39); the factorial identities and valuation bounds cover the entire indicated ranges.
