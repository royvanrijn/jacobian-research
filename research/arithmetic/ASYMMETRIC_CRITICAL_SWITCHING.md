# Critical-level switching for the modulus-12 binary problem

A half-dimensional sieve can be stopped just below the square-root scale.
Its remaining unsuccessful representations have exactly two large
forbidden prime factors. The estimate below controls their total weight
uniformly in the target, including targets with arbitrary prime-power
factors. It reduces a sufficient distribution hypothesis for two
modulus-12 summands to the square-root level, up to logarithmic losses.
That distribution hypothesis remains open.

Let
\[
 C=\{n\ge1:p\mid n\Rightarrow p\equiv1\pmod{12}\},\qquad
 A=\{n\ge1:p\mid n\Rightarrow p\equiv1\pmod3\},
\]
including 1, and put
\[
 U=\tau^2\mathbf1_C,\qquad V=\tau\mathbf1_A,\qquad
 X_N=N\mathfrak S_{\rm asym}(N).
\]
Here \(\mathfrak S_{\rm asym}\) is the exact positive singular series
in [Theorem AsymPointwise](ASYMMETRIC_GAUSSIAN_EISENSTEIN_SUMS.md).
Throughout, \(N\equiv2\pmod{12}\). Use the definitions from
[the half-dimensional sieve criterion](ASYMMETRIC_SEMILINEAR_SIEVE.md):
\[
 \begin{gathered}
 w_N(b)=U(N-b)V(b),\qquad
 A_d(N)=\sum_{1\le b<N,\ d\mid b}w_N(b),\\
 g_N(p)=\begin{cases}
 \dfrac{2p-1}{p^2-p+1},&p\equiv7\ (12),\ p\nmid N,\\
 0,&\text{otherwise},
 \end{cases}\qquad
 r_d(N)=A_d(N)-X_Ng_N(d),\\
 P_7(z)=\prod_{\substack{p<z\\p\equiv7\ (12)}}p,
 \qquad \mathcal V_N(z)=\prod_{p<z}(1-g_N(p)).
 \end{gathered}
 \tag{1}
\]
The density \(g_N\) is multiplicative on squarefree sieve divisors.

## Theorem AsymTwoLargePrimes

There is an absolute constant \(C_0\) such that, for every fixed
\(0<\epsilon<1/100\), every sufficiently large admissible target satisfies
\[
 B_z(N):=
 \sum_{\substack{1\le b<N\\(b,P_7(z))=1\\b\notin C}}w_N(b)
 \le C_0\epsilon^{3/2}X_N\mathcal V_N(z),
 \qquad z=N^{1/2-\epsilon}.
 \tag{2}
\]
The constant \(C_0\) is independent of \(\epsilon\). The onset may
depend on \(\epsilon\). There is no restriction on the prime factors
or target valuations of \(N\).

## Corollary AsymCriticalCriterion

Assume the following **open distribution hypothesis**: for every fixed
\(0<\epsilon<1/100\) and every \(A_0>0\), there are constants
\(B_0\ge0,C_{\epsilon,A_0},N_{\epsilon,A_0}\) such that
\[
 \sum_{\substack{d<D\\d\mid P_7(N^{1/2-\epsilon})}}|r_d(N)|
 \le C_{\epsilon,A_0}X_N(\log N)^{-A_0},
 \qquad D=N^{1/2}(\log N)^{-B_0},
 \tag{BVcrit}
\]
for every admissible \(N\ge N_{\epsilon,A_0}\).
Then every sufficiently large \(N\equiv2\pmod{12}\) belongs to
\(C+C\). For every \(\nu>0\), it has
\(\gg_\nu N^{1-\nu}\) distinct ordered value pairs in \(C^2\).

Only one sufficiently small fixed \(\epsilon\), chosen after the
absolute constants in the proof, is needed. The stated quantifiers in
(BVcrit) are a convenient sufficient formulation. An appropriate signed
lower-sieve remainder bound also suffices, as specified below.
Neither (BVcrit) nor this unconditional two-sided coverage is proved here.

## 1. Why only two forbidden primes remain

If \(w_N(b)>0\), then \(N-b\equiv1\pmod{12}\), so
\(b\equiv1\pmod{12}\). All prime factors of \(b\in A\) are 1 or
7 modulo 12. Hence
\[
 \Omega_7(b):=\sum_{p\equiv7\ (12)}v_p(b)\equiv0\pmod2.
 \tag{3}
\]
Every forbidden factor in (2) is at least \(z\), while
\(z^3>N\). A bad surviving value consequently has the unique form
\[
 b=cpq,\qquad c\in C,\quad c\le Y:=N^{2\epsilon},\quad
 z\le p\le q,\quad p,q\text{ primes }7\pmod{12}.
 \tag{4}
\]
Repeated factors \(p=q\) are allowed. Neither prime divides \(N\),
since it would then divide \(N-b\in C\). Since \((c,pq)=1\),
\[
 V(cpq)\le4V(c).
 \tag{5}
\]
Put \(w=N^{1/10}\); the chosen range of \(\epsilon\) gives \(Y<w<z\).
In particular every prime factor of \(c\) is below \(w\).

## 2. The corrected Henriot upper bound

We use K. Henriot's **New Theorem 5**, printed p.377 in the
[2014 erratum](https://doi.org/10.1017/S0305004114000280), with the
multiplicative growth-class definitions from
[the original paper, Section 2](https://arxiv.org/html/1102.1643v1#S2).
The corrected exact-valuation factors are essential here.

A useful observation concerns the source's local factor. For the two
primitive distinct linear forms used below, let \(P_\ell(\boldsymbol\nu)\)
be the exact-valuation probability in the corrected Euler sum at
\(\ell\ge5\). The \((0,0)\) term of that sum is 1, with no conditions.
The nonzero valuation events are disjoint and partition the complement
of the event where both polynomial values are units. Its measure is
\(1-\rho_Q(\ell)/\ell\). Thus, for the nonnegative local weights here,
\[
 \begin{split}
 &(1-\rho_Q(\ell)/\ell)
 \left(1+\sum_{\boldsymbol\nu\ne(0,0)}
 F(\ell^{\nu_1},\ell^{\nu_2})P_\ell(\boldsymbol\nu)\right)\\
 &\hspace{12mm}\le
 \int_{\mathbb Z_\ell}
 F\bigl(\ell^{v_\ell(Q_1(t))},
        \ell^{v_\ell(Q_2(t))}\bigr)\,dt.
 \end{split}
 \tag{6}
\]
Indeed, multiplying the nonzero, positive terms by
\(1-\rho_Q(\ell)/\ell\le1\) only decreases them. This compares the
source's full local factor with an actual integral; its Euler sum
alone is not that integral. A zero coordinate opposite a positive
coordinate imposes a unit, as in the corrected cross-coprimality rule.

Fix \(c\in C\), \(c\le Y\), and a prime \(p\ge z\),
\(p\equiv7\pmod{12}\), \(p\nmid N\). Let
\[
 h=(c,N)\in C,\quad F_h(n)=U(hn)/U(h),\quad
 f_w(n)=\mathbf1_{(n,\prod_{\ell\le w}\ell)=1}.
\]
The relevant prime \(q\)-sum is bounded by the rough sum for
\[
 Q_1(t)=12t+7,\qquad
 Q_2(t)=\frac Nh-\frac{cp}{h}(12t+7),\qquad
 F(u,v)=f_w(u)F_h(v),
 \tag{7}
\]
with external factor \(U(h)\). This is a multiplicative rough upper
bound for the prime condition, not a use of the prime indicator as a
multiplicative function.

Both forms in (7) are primitive. We have
\((N/h,cp/h)=1\), and \(Q_2(t)\equiv1\pmod{12}\).
Their product is primitive, has degree 2, and has coefficient norm
\(O(N)\). Uniformly in \(h\),
\[
 F_h(\ell^e)\le(e+1)^2\le4^e,\qquad
 F_h(m)\le\tau(m)^2\ll_\xi m^\xi.
 \tag{8}
\]
The corresponding multiplicative-growth ratios have the same bounds.
Thus the parameters of Henriot's \(\mathcal M_2\) class are fixed,
independent of \(h\); \(F\) itself is the coefficient majorant.

If \(Q=\lfloor N/(cp)\rfloor<z\), there is no relevant prime.
Otherwise cover \([z,Q]\) from the top by geometric intervals
\([2T/3,T]\), stopping once the lower endpoint lies below \(z\).
Their total lengths are \(O(Q)\), and their scales are at least a
constant times \(z\). After \(q=12t+7\), choose a half-open integer
interval with upper endpoint \(\lfloor(T-7)/12\rfloor\). Its parameters
satisfy \(x\asymp T\), \(y\asymp T/3\le x\), and all \(Q_2\)
values stay positive. No polynomial zero is present because \(p\nmid N\).
For sufficiently large \(N\), \(T\ge N^{0.48}\). The fixed choices
\[
 \alpha=1/2,\qquad\Delta=1/4,\qquad\xi=1/10000
 \tag{9}
\]
meet \(x\ge C\|Q_1Q_2\|^\Delta\), \(x^\alpha<y\le x\), and
\(\xi<\alpha/[50\cdot2(2+\Delta^{-1})]\).
All constants are independent of \(\epsilon,c,p,N,h\).

Drop the product cutoff in the source's nonnegative coefficient sum
and factor the larger sum over primes up to its cutoff. Inequality (6)
then permits actual local integrals as upper factors. At 2 and 3 the
two polynomial values are units; these finite factors are bounded.

For \(\ell>w\), no prime factor of \(c\) or \(h\) occurs.
Except at the single outer prime \(p\), the second polynomial has unit
slope. Since \(f_w\) has all local values 1 there, its local integral is
\[
 M_{U,\ell}=(1-\ell^{-1})\sum_{e\ge0}U(\ell^e)\ell^{-e},
 \tag{10}
\]
independently of \(v_\ell(N)\). At \(\ell=p\), the second polynomial
is a nonzero constant modulo \(p\), so the local integral is 1 rather
than \(M_{U,p}=1-1/p\). This costs one bounded factor. The product
of (10) over any range \(w<\ell\le x\) is bounded uniformly by the
fixed-progression Mertens estimates, or the fixed \(L\)-factorization
of \(U\). All primes dividing \(h\) are below \(Y<w\), and their
content factors will be included exactly in the next section.

## 3. Exact small-prime factors and cofactor restoration

For a good prime \(\ell\equiv1\pmod{12}\), \(\ell\le w\), put
\(u=1/\ell\), \(q_\ell=1-u\),
\(r=v_\ell(N)\), \(j=v_\ell(c)\), and \(A_i=(i+1)^2\).
Distribute the content factor \(U(h)\) locally. The actual prime-unit
\(q\)-variable integral, including that factor, is
\[
 K_{r,j}=\begin{cases}
 q_\ell A_{\min(r,j)},&j\ne r,\\
 (1-2u)A_r+q_\ell\displaystyle\sum_{k\ge1}A_{r+k}u^k,&j=r.
 \end{cases}
 \tag{11}
\]
For \(j<r\), the normalized target is divisible by \(\ell\), and a
unit prime variable makes the other polynomial a unit. For \(j>r\),
that polynomial is a nonzero constant modulo \(\ell\). For \(j=r\),
its constant and slope are units; splitting the unit-variable integral
at its root gives (11). In particular,
\[
 K_{0,0}=M_{U,\ell}-u>0,\qquad
 K_{r,0}=q_\ell\quad(r\ge1).
 \tag{12}
\]
The external content at this prime is \(A_{\min(r,j)}\), so (6)
indeed bounds the normalized source factor by (11), uniformly in both
valuations. At a bad prime \(\ell\equiv5,7,11\pmod{12}\),
\(\ell\le w\), both \(c,h\) are units and \(\ell\ne p\).
The local integral is \(1-1/\ell\) if \(\ell\mid N\), and
\(1-2/\ell\) otherwise. Their ratio is
\[
 j_\ell=\frac{\ell-1}{\ell-2}.
\]
Set
\[
 J_{\rm bad,w}(N)=
 \prod_{\substack{\ell\mid N,\ \ell\le w\\
                   \ell\equiv5,7,11\ (12)}}j_\ell,
 \qquad
 \Gamma_N(c)=
 \prod_{\substack{\ell\le w\\\ell\equiv1\ (12)}}
        \frac{K_{v_\ell(N),v_\ell(c)}}{K_{0,0}}.
 \tag{13}
\]
Only primes dividing \(cN\) change this finite product.

The generic small-prime baseline is \(K_{0,0}=M_{U,\ell}-1/\ell\)
at good primes and \(1-2/\ell\) at bad primes. Its first-order
coefficients are 2 and \(-2\), respectively, whose prime average is
\(-1\). Its product is \(\ll1/\log w\). Together with (10),
Henriot therefore gives
\[
 \sum_{\substack{z\le q\le N/(cp)\\q\text{ prime},\ q\equiv7\ (12)}}
       U(N-cpq)
 \ll\frac{N}{cp\log N}
          J_{\rm bad,w}(N)\Gamma_N(c).
 \tag{14}
\]
The same constant works for every target valuation and every \(c\le Y\).

Write \(I_{\ell,r}\) for the exact asymmetric local convolution factor.
At a good prime it is
\[
 \begin{split}
 I_{\ell,r}={}&q_\ell\sum_{i<r}A_i(i+1)u^i
   +(1-2u)A_r(r+1)u^r\\
 &+q_\ell\left((r+1)\sum_{i>r}A_i u^i
              +A_r\sum_{j>r}(j+1)u^j\right).
 \end{split}
\]
The following identity retains all valuations:
\[
 \sum_{j\ge0}(j+1)u^jK_{r,j}=I_{\ell,r}.
 \tag{15}
\]
Its terms with \(j<r\), \(j=r\), and \(j>r\) give exactly the
three corresponding parts of the preceding display.

By positivity, enlarge the harmonic cofactor sum to its Euler product.
Since every prime factor of \(c\le Y\) is below \(Y<w\), (15) gives
\[
 \begin{split}
 \sum_{\substack{c\le Y\\c\in C}}\frac{V(c)}c\Gamma_N(c)
 &\le\prod_{\substack{\ell\le Y\\\ell\equiv1\ (12)}}
            \frac{I_{\ell,v_\ell(N)}}{K_{0,0}}
   \prod_{\substack{Y<\ell\le w,\ \ell\mid N\\\ell\equiv1\ (12)}}
            \frac{q_\ell}{K_{0,0}}\\
 &\ll\sqrt{\log(2Y)}
   \prod_{\substack{\ell\le Y,\ \ell\mid N\\\ell\equiv1\ (12)}}
            \frac{I_{\ell,v_\ell(N)}}{I_{\ell,0}}.
 \end{split}
 \tag{16}
\]
The second product in the first line is at most 1 by (12).
For the generic first product,
\(I_{\ell,0}/K_{0,0}=1+2/\ell+O(\ell^{-2})\); fixed-progression
Mertens gives its \(\sqrt{\log(2Y)}\) upper bound.
The target factors in (16) are exactly those in AsymPointwise.
There is no excess logarithmic target factor hidden in this step.

## 4. The switched mass and the target comparison

Define
\[
 F_{N,Y,w}=J_{\rm bad,w}(N)
 \prod_{\substack{\ell\le Y,\ \ell\mid N\\\ell\equiv1\ (12)}}
           \frac{I_{\ell,v_\ell(N)}}{I_{\ell,0}}.
\]
Equations (4)–(5), (14), and (16) yield
\[
 \begin{split}
 B_z(N)
 &\ll\frac{NF_{N,Y,w}}{\log N}\sqrt{\log(2Y)}
      \sum_{\substack{z\le p\le\sqrt N\\p\equiv7\ (12)}}\frac1p\\
 &\ll\epsilon^{3/2}\frac{NF_{N,Y,w}}{\sqrt{\log N}}.
 \end{split}
 \tag{17}
\]
For each fixed \(\epsilon\), the prime sum is
\(\frac14\log((1/2)/(1/2-\epsilon))+o(1)\ll\epsilon\).
The implied constant is uniform for \(0<\epsilon<1/100\), after
allowing the onset to depend on \(\epsilon\).
Also \(\sqrt{\log(2Y)}\ll\sqrt{\epsilon\log N}\) for sufficiently
large \(N\). Squares \(p=q\) were included in (5).

The exact target product gives the uniform comparison
\[
 X_N\mathcal V_N(z)\gg
       \frac{NF_{N,Y,w}}{\sqrt{\log N}}.
 \tag{18}
\]
At target primes \(5,11\pmod{12}\), the asymmetric target ratio is
\(j_p\). At target primes \(7\pmod{12}\) below \(z\), the combined
singular-series and sieve ratio is also exactly \(j_p\), by equation
(4) of the half-dimensional sieve note. There are at most two target
primes \(7\pmod{12}\) at least \(z\), since \(z^3>N\). Their
uncancelled ratios
\[
 \frac{I_{p,r}}{I_{p,0}}=\frac{(p-1)^2}{p^2-p+1}
\]
have an absolute positive lower bound. Every good target ratio is at
least 1: the positive-increments identity (7) in that note gives
\[
 I_{p,r}-I_{p,0}
 =\sum_{\substack{i,j\ge1\\\min(i,j)\le r}}
   (U(p^i)-U(p^{i-1}))(V(p^j)-V(p^{j-1}))
                  p^{-\max(i,j)}\ge0.
\]
Thus retaining just the target factors in \(F_{N,Y,w}\) decreases the
product, apart from the bounded two uncancelled factors. The generic
singular product is a positive constant, and
\(\mathcal V_0(z)\asymp(\log z)^{-1/2}\). This proves (18).
Combining (17) and (18) proves Theorem AsymTwoLargePrimes.

## 5. The critical-level implication

Assume (BVcrit), take a fixed \(\epsilon\), and apply the lower
half-dimensional sieve from ASYMMETRIC_SEMILINEAR_SIEVE.md at
\[
 z=N^{1/2-\epsilon},\qquad
 D=N^{1/2}(\log N)^{-B_0},\qquad
 s=\frac{\log D}{\log z}\longrightarrow\frac1{1-2\epsilon}.
\]
Its upper Mertens constant is uniform in the target; its quantitative
error does not involve the target-dependent lower Mertens constant.
For \(1<s\le2\), the lower function is
\[
 f_{1/2}(s)=\sqrt{\frac{e^{\gamma_E}}{\pi s}}
       \int_1^s\frac{dt}{\sqrt{t(t-1)}}.
\]
Consequently there is an absolute \(c_1>0\) such that the sifted mass
is at least
\[
 c_1\sqrt\epsilon\,X_N\mathcal V_N(z)
 \tag{19}
\]
for sufficiently large \(N\), with fixed \(\epsilon\).
The remainder in (BVcrit) is negligible: take \(A_0>1\) and use
\(\mathcal V_N(z)\gg(\log N)^{-1/2}\).
The source's logarithmic main-term error also tends to zero.

Subtract (2) from (19). Choose one sufficiently small fixed
\(\epsilon\) so that \(C_0\epsilon^{3/2}<c_1\sqrt\epsilon/2\).
The remaining weighted \(C+C\) mass is positive for every sufficiently
large admissible target. Moreover
\(X_N\mathcal V_N(z)\gg N/\sqrt{\log N}\), by the same target
comparison as (18), or by (6) in the sieve note and monotonicity in \(z\).
Divisor bounds convert this to \(\gg_\nu N^{1-\nu}\) distinct pairs.
This proves the conditional corollary.

For the lower half-dimensional weights at these \(z,D\), it is enough
to prove the weaker, signed open condition
\[
 \sum_{d<D}\lambda^-_d r_d(N)
       \ge-o\bigl(X_N\mathcal V_N(z)\bigr).
 \tag{20}
\]
It need hold for only one sufficiently small fixed \(\epsilon\).
The little-o is uniform over the admissible targets; an absolute
little-o bound for the signed sum also suffices. The stronger absolute
sum in (BVcrit) avoids prescribing these weights in its statement.

## Dependencies and scope

The proof uses the exact local model and positivity from AsymPointwise,
the half-dimensional sieve and signed weights in
ASYMMETRIC_SEMILINEAR_SIEVE.md, Henriot's corrected New Theorem 5, and
prime Mertens estimates for the fixed modulus 12; for example,
[Loughran–Matthiesen, Lemma 2.4(2)](https://arxiv.org/html/1904.12845v2#S2)
applied to the four prime-class indicators. It uses no
varying-modulus prime number theorem and inserts no prime weight into
BGLM's signed modulus average.

The new unconditional result is the upper bound (2). The missing
square-root distribution estimate (BVcrit), or its signed substitute
(20), is essential to the binary conclusion. Neither the existing
pointwise asymmetric convolution nor its small-coefficient extension
supplies that estimate. The separate
[auxiliary-average theorem](ASYMMETRIC_AUXILIARY_AVERAGE.md) now proves a
positive \(L^1\) distribution level \(D=N^{\xi-\delta}\), for every
fixed \(0<\delta<\xi\), using a coefficient average through \(N^\xi\).
Its \(\xi\) is within the small auxiliary-modulus range, far below the
critical exponent \(1/2\). Expanding the companion divisor sum at
coefficients near the square-root level still creates combined moduli
near \(N^{3/4}\); no cancellation there is proved in this note.
No balance statement, effective cutoff, finite coverage completion,
or priority claim is made.
