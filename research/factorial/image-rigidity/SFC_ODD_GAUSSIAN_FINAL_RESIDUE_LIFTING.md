# Final Hermite residues at dominant initial prime powers

The quadratic Gaussian window is zero-free for every
\(q=(2k+1)p^e\) with \(k\ge0,e\ge1\) and odd prime
\(p>4k+2\). No quadratic-character or final-resultant exception
remains in this dominant initial-prime family. In particular both
\(q=5\cdot13^e\) and \(q=9\cdot89^e\) are covered at every
positive exponent. The unrestricted odd-start problem remains open.

At nonresidue primes, one secondary layer separates the two final
Hermite equations. When \(2\) is a quadratic residue, that layer
contributes a unit multiple of the original Wronskian. The next
layer has relative odd/even phase two and forces a Wronskian
valuation below one, contradicting its squared-gcd lower bound.

## 1. Statement and the exact reduction

Write \(G_m(A,B,C)=\mathbb E(A+BX_m+CX_m^2)^m\), where
\(X_m\sim\mathcal N(0,1/m)\). Put
\[
 S_k(Z)=\sum_{j=0}^k\frac{k_{\underline j}}{(2j+1)!}Z^j,
 \qquad P_k(Z)=\sum_{j=0}^k\frac{k_{\underline j}}{(2j)!}Z^j.
\]

**Final-residue theorem.** Let \(q=(2k+1)p^e\), with \(k\ge1\),
\(e\ge1\), and prime \(p>4k+2\). Then
\(G_q,G_{q+1},G_{q+2}\) have no common projective zero over
any characteristic-zero field, including when \((2/p)=1\)
and \(p\mid T_k\).

Together with the [internal shared-edge theorem](SFC_ODD_GAUSSIAN_SHARED_EDGE_LIFTING.md),
the dominant initial family now has the sole prime-size condition
\[
             q=(2k+1)p^e,\quad k\ge0,\ e\ge1,
                         \quad p>4k+2.                     \tag{1}
\]
The case \(k=0\) has no final shared slope and is already
covered by the internal shared-edge argument.
Every fixed homogeneous bivariate trinomial support consequently has
a computable factorial corner exclusion at all these prescribed
starts, by the exact Gaussian transfer. No new numerical cutoff is
evaluated here.

Use the [initial-order reduction](SFC_ODD_GAUSSIAN_LOCAL_WINDOW_TESTS.md#7-large-prime-powers-dividing-the-initial-order).
At a hypothetical common algebraic zero, in \(\overline{\mathbb Q}_p\),
\[
 B=1,\quad x=A/B,\quad h=C/B,\quad v(p)=1,\qquad
 v(h)\ge b=e/2-L/2,\quad L=\sum_{j=1}^e p^{-j}.
\]
The initial-order reduction and internal shared-edge theorem
exclude first slopes and internal slopes without any character
restriction. Thus the only remaining possibility for \(k>0\) is
\(v(x)=\alpha=L/2<v(h)\), on the final shared slope.

Set \(M=p^e\), \(\ell=ML=(M-1)/(p-1)\), \(N=(q+1)/2\), and
\[
 E(Y)=\sum E_iY^i,\quad E_i=\frac{(2(q+1))^iN_{\underline i}}{(2i)!},
 \qquad O(Y)=\sum O_iY^i,\quad O_i=\frac{(2(q+2))^iN_{\underline i}}{(2i+1)!}.
                                                               \tag{2}
\]
The normalized pure moments at the two unit orders are \(E(x^2)\)
and \(xO(x^2)\).

## 2. A common coordinate removes the first quadratic perturbation

This step is needed only when \(e=1\). For a unit order \(m\), let
\(F_m(x,h)=\mathbb E(x+X_m+hX_m^2)^m\) and \(Q_m(x)=F_m(x,0)\).
Gaussian integration by parts gives the exact polynomial identity
\[
 \partial_hF_m(x,0)=(1+x^2)Q'_m(x)-mxQ_m(x).
\]
For example, writing \(z=(x+z)-x\) in the differentiated moment
and using the Gaussian recurrence at variance \(1/m\) gives the
displayed formula. Hence the common change
\[
 u=x+(1+x^2)h,\qquad
 F_m(x,h)=(1-mxh)Q_m(u)+\mathcal R_m(x,h)             \tag{3}
\]
has remainder divisible by \(h^2\).

For completeness the valuation bound also follows coefficientwise.
If \(a_l\) is a pure coefficient, a nonzero coefficient of \(x^ih^j\)
in \(F_m\) is
\[
 a_{i+j}\binom{i+j}{j}
       \prod_{t=1}^j\frac{m-i-j+2t-1}{m}.
\]
The ratio is \(p\)-integral. At \(v(x)=\alpha<\beta=v(h)\), a
term with \(j\) quadratic factors lies at least \(j(\beta-\alpha)\)
above the pure minimum. Taylor expansion of \(Q_m(u)\) has the
same lower bound; multiplication by \(mxh\) only improves it.
After the exact cancellation of degrees zero and one in \(h\),
\(\mathcal R_m\) therefore lies at least \(2(\beta-\alpha)\)
above the pure minimum. The scalar in (3) is a unit and \(v(u)=\alpha\).
The change is the same for both unit orders; only the unit scalar
depends on \(m\).

For \(e=1\), use \(Y=u^2\). For \(e\ge2\), use \(Y=x^2\)
and the original first-order bound. In either case \(v(Y)=L\).

## 3. The next digit layer

Put \(r=(M+1)/2\) and
\(Z=(-1)^\ell2Y^M/p^\ell\), a unit. Give coefficients the weight
\(v(c_i)+iL\). The minimum in \(E\) and \(YO\) is
\(\gamma=-b\), at indices \(r+jM\), \(0\le j\le k\).

If \(e\ge2\), let
\[
 t=(M/p+1)/2,\qquad \delta=(p+1)/(2p),\qquad
 \epsilon=\frac{E_tY^t}{E_rY^r},\quad v(\epsilon)=\delta,
 \qquad\chi=(2/p).
\]
With \(o_v(\delta)\) denoting an element of valuation strictly
greater than \(\delta\), the exact coefficient expansions are
\[
 \frac{E(Y)}{E_rY^r}=S_k(Z)+\epsilon P_k(Z)+o_v(\delta),
\]
\[
 \frac{YO(Y)}{O_{r-1}Y^r(1+2Y)}
       =S_k(2Z)+\chi\epsilon P_k(2Z)+o_v(\delta).      \tag{4}
\]
When \(e=1\), the corresponding formulas have
\[
 \delta=b=(p-1)/(2p),\quad \epsilon=(E_rY^r)^{-1},
\]
\[
 \frac{E(Y)}{E_rY^r}=S_k(Z)+\epsilon P_k(Z)+o_v(\delta),
 \qquad
 \frac{YO(Y)}{O_{r-1}Y^r(1+2Y)}=S_k(2Z)+o_v(\delta). \tag{5}
\]

Here is the complete layer check behind (4)–(5). As in the previous
shared-edge proof,
\[
 v(E_i)=-v(i!)-\sum_{a=1}^e
          \mathbf1_{i\bmod p^a=(p^a+1)/2},\qquad
 v(O_i)=v(E_i)-v(2i+1).
\]
If the initial digit prefix has positive length \(a\), the even
weight is \((s_p(i)-i/M)/(p-1)-a\). At \(i=(p^a+1)/2\),
its gap above \(\gamma\) is
\[
                  (e-a)/2+\frac{1-p^{a-e}}{2(p-1)}.
\]
Writing \(i=(p^a+1)/2+p^a u\) adds exactly
\((s_p(u)-p^{a-e}u)/(p-1)\) to this base weight. All its
lower-digit contributions are positive, and its top-digit
contribution is zero. For \(a=e\), all indices \(r+jM\) are at
the minimum. For
\(a=e-1\), write the next digits as \(jp+v\); the extra gap is
\(v/p\), so exactly \(t+jM\) occur at the next layer \(\delta\).
Shorter positive prefixes have gap greater than one. Prefix zero
has nonnegative absolute weight, outside this layer for \(e\ge2\).
For \(e=1\), precisely \(jM\) have absolute weight zero; these
are the first correction at gap \(b\).

For \(YO\), write its index as \(i=j+1\). If \(v(2j+1)=0\),
the even weight at \(j\) is shifted by \(L\). If
\(j=(p^a-1)/2+p^av\), subtract \(a\) from that weight before
shifting. The resulting weight is
\[
 -a/2-\frac{p^a-1}{2M(p-1)}+L
             +\frac{s_p(v)-p^{a-e}v}{p-1}.
\]
Its base weight equals the even base weight above because
\(L=(1-1/M)/(p-1)\), and the same digit calculation applies.
No \(a>e\) occurs since \(2j+1<p^{e+1}\). Thus the complete
inventories
through the stated gap are

| Case | \(E\) indices | \(YO\) indices |
|---|---|---|
| \(e\ge2\) | \(r+jM\), \(t+jM\) | \(r+jM\), \(r+1+jM\), \(t+jM\) |
| \(e=1\) | \(r+jM\), \(jM\) | \(r+jM\), \(r+1+jM\) |

To compute the units, Wilson's block decomposition gives
\(\operatorname{unit}(n!)\equiv(-1)^{v(n!)}\prod n_a!\pmod p\)
for base-\(p\) digits \(n_a\). At \(r+jM\), the high factorial
digits give \(k_{\underline j}/(2j+1)!\); at \(t+jM\), or
\(jM\) in the first-exponent case, they give
\(k_{\underline j}/(2j)!\). The even and odd powers contribute
\(2^j\) and \(4^j\), respectively, with sign \((-1)^{j\ell}\).
This yields the shapes in (4)–(5). The child ratios are all
\(O_{r+jM}/O_{r+jM-1}\equiv2\pmod p\), giving \(1+2Y\).

For \(e\ge2\), the ratio of the normalized odd and even secondary
coefficients is
\[
 \frac{E_r/O_{r-1}}{E_t/O_{t-1}}
 \equiv2^{-(r-t)}=2^{-(p-1)p^{e-1}/2}=\chi\pmod p.
\]
This is the sign in (4). All coefficient congruence errors gain an
integer valuation step, which is beyond \(\delta<1\).

The original quadratic terms also lie beyond the layer. For
\(e\ge2\), their relative valuation is at least
\(b-\alpha=e/2-L>\delta\). For \(e=1\), (3) improves the
relative valuation to \(2(b-\alpha)=(p-2)/p>\delta\).
Consequently a common Gaussian zero satisfies both right-hand sides
of (4), or of (5), equal to zero up to \(o_v(\delta)\).

## 4. Compatibility at one common residue

The polynomial identities
\[
 P_k=S_k+2ZS_k',\qquad
 4ZS_k''+(6+Z)S_k'-kS_k=0                         \tag{6}
\]
follow directly from the coefficients. They identify the classical
Kummer differential equation, but no imported special-function
argument is needed. Since \(p>4k+2\), the reduction of \(S_k\)
has degree \(k<p\), constant term one, and is squarefree. Indeed
a common nonzero root of \(S_k,S_k'\) would, by differentiating
the second identity, annihilate every derivative through degree
\(k\), contradicting the unit leading coefficient and \(k!\ne0\).

Write \(A=S_k(Z)\), \(B=S_k(2Z)\), and let \(H\) be an integral
lift of their monic gcd modulo \(p\). If \(H=1\) there is no
common residue. At a common residue \(z\), squarefreeness implies
that the reduced quotient \(\bar A/\bar H\) is nonzero there.
Lift that quotient to an integral polynomial \(A_0\), so
\(A=HA_0+pA_1\) with \(A_0(Z)\) a unit. From either pair of
equations, \(v(A)\ge\delta\). Since \(\delta<1\), this implies
\(v(H(Z))\ge\delta\).

The Wronskian \(W=B'A-A'B\) is divisible by \(H^2\) modulo
\(p\): substitute \(\bar A=H A_0\), \(\bar B=H B_0\) and
cancel the two terms containing \(H'\). Therefore
\[
                       v(W(Z))\ge\min(2\delta,1)>\delta.    \tag{7}
\]
On the other hand (4) gives
\[
 W(Z)=-\epsilon\bigl(B'P_k(Z)-\chi A'P_k(2Z)\bigr)
                        +o_v(\delta).
\]
At a common residue, (6) reduces the coefficient in parentheses to
\[
                  4z(1-\chi)S_k'(z)S_k'(2z).
\]
For \(\chi=-1\) it is a unit. This forces \(v(W)=\delta\),
contradicting (7). For \(e=1\), equation (5) instead gives
\(W=-\epsilon B'P_k+o_v(\delta)\); its coefficient is
\(4zS_k'(z)S_k'(2z)\), again a unit, regardless of \(\chi\).
This proves the final-residue theorem.

## 5. The requested family and the next boundary

At \(k=2,p=13\), \(\chi=-1\) and
\[
 S_2=1+Z/3+Z^2/60,\qquad P_2=1+Z+Z^2/12.
\]
The reductions of \(S_2(Z),S_2(2Z)\) have gcd \(Z-2\).
At \(Z=2\), the compatibility coefficient for \(e\ge2\) is
\[
 \frac{14}{15}\frac{10}{3}+\frac25\frac{19}{3}
                     =\frac{254}{45}\equiv12\pmod{13}.     \tag{8}
\]
Equivalently the two equations require lift corrections \(9\) and
\(4\) modulo \(13\). The exact squared-residue identity is
\[
       90(B'A-A'B)=(Z-2)^2+13(Z+2).
\]
This excludes every exponent \(e\ge2\), with secondary layer
\(\delta=7/13\). The first exponent is excluded by (5), at
\(\delta=6/13\). Thus \(q=5\cdot13^e\) is settled for all
\(e\ge1\) without extrapolating a finite exponent check.

The [checker](verify_sfc_odd_gaussian_final_layer.py) audits ten
selected first-layer inventories and seven second-layer inventories,
their normalized shapes and signs, the differential identities,
Wronskian divisibility and 6948 exact straightening remainder
coefficients. It rejects an incorrect first-layer sign and an
incorrect second-layer ratio. The [receipt](sfc_odd_gaussian_final_layer_verification.json)
records an internal pass with Python 3.12.3. Run from the repository root:

```sh
python3 research/factorial/image-rigidity/verify_sfc_odd_gaussian_final_layer.py
```

These are finite exact controls of the written proof. They do not
constitute independent mathematical review or formal verification.

The factorial-unit calculation is the elementary modulo-prime case
of [Granville's prime-power treatment, Theorem 1](https://www.cecm.sfu.ca/organics/papers/granville/paper/binomial/html/node2.html).
The differential identity corresponds to
[DLMF 13.2.1](https://dlmf.nist.gov/13.2.E1) after
\(S_k={}_1F_1(-k;3/2;-Z/4)\). The local arguments above derive
the required specializations directly; no priority claim is made.

For \(e\ge2\) and \(\chi=1\), the first secondary compatibility
does vanish at a common residue. The profile \(q=9\cdot89^e\)
has \(\gcd(S_4(Z),S_4(2Z))=Z+25\) modulo \(89\), with residue
\(Z=64\). It remains a non-exclusion control for the first layer.
The following next-layer proof excludes its lift uniformly, and
also removes every other positive-character final-residue exception
under (1).

## 6. The second secondary layer at positive-character primes

Assume \(e\ge2\) and \(\chi=(2/p)=1\). Then \(p\ge7\).
Keep \(M,r,t,L,\gamma=-b\) from Section 3, and put
\[
 J=M/p,\qquad s=t+J,\qquad
 \delta_0=\frac{p+1}{2p},\qquad
 \delta_1=\frac{p+3}{2p}<1,\qquad
 \eta=\frac{E_sY^s}{E_rY^r}.                              \tag{9}
\]
The exact valuations are \(v(\epsilon)=\delta_0\) and
\(v(\eta)=\delta_1\). Through weight \(\gamma+\delta_1\),
the complete inventories, with \(0\le j\le k\), are

| Polynomial | Indices | Relative weight |
|---|---|---|
| \(E\), \(YO\) | \(r+jM\) | \(0\) |
| \(YO\) | \(r+1+jM\) | \(L\) |
| \(E\), \(YO\) | \(t+jM\) | \(\delta_0\) |
| \(E\), \(YO\) | \(s+jM\) | \(\delta_1\) |

Here the index in \(YO\) is its actual power of \(Y\).
In particular the child of the first secondary term is absent:
\(\delta_0+L>\delta_1\), because \(e\ge2\) gives \(L>1/p\).

To prove completeness, use the prefix calculation of Section 3.
A prefix of length \(a\le e-2\) has relative gap greater than
one. At prefix length \(e-1\), write the remaining digits as
\(v+pj\), \(0\le v<p\). Their exact additional gap is
\(v/p\); hence only \(v=0,1\) survive the cutoff in (9).
Neither equals \((p-1)/2\), so both have exactly this prefix
length. Prefix length \(e\) gives precisely \(r+jM\).
For prefix zero the absolute even weight is nonnegative:
writing \(i=jM+u\), \(0\le u<M\), gives
\((s_p(u)-u/M)/(p-1)\ge0\). It is outside the retained range,
since \(b>\delta_1\).

For \(YO\), write its coefficient as \(O_{i-1}Y^i\).
If \(v(2i-1)>0\), the parent-prefix calculation in Section 3
gives exactly the same base weights and higher-digit increments.
If \(v(2i-1)=0\), its weight is the even weight at \(i-1\)
plus \(L\). Only the leading parents then supply children;
every secondary child lies beyond the retained range. Also
\(2i-1<p^{e+1}\), so there is no longer prefix to inspect.
This proves the full inventory, including terms away from the hull.

### The same polynomial, with a doubled inter-row phase

At either secondary base \(u_v=t+vJ\), \(v=0,1\), the normalized
even coefficients at \(u_v+jM\) have shape \(P_k(Z)\).
Indeed \(N-u_v-jM=(k-j)M+((p-1)/2-v)J\), while
\(2(u_v+jM)=(2j)M+(2v+1)J+1\). Since \(2v+1<p\),
the top factorial digits give \(k_{\underline j}/(2j)!\)
without a carry. The lower digits do not depend on \(j\).
Wilson's block formula, including its factorial-valuation sign,
and the power \((2(q+1))^{jM}\) give precisely the scale
\(((-1)^\ell2Y^M/p^\ell)^j=Z^j\).

For the normalized odd/even ratio use the exact identity
\[
 \kappa_i:=\frac{E_i}{O_{i-1}}
 =\frac{(q+1)(N-i+1)}i
       \left(\frac{q+1}{q+2}\right)^{i-1}.                \tag{10}
\]
At all retained parents \(i\equiv1/2\pmod p\), so
\(\bar\kappa_i=2^{2-i}\). Thus
\[
 \frac{\kappa_r}{\kappa_{u_v}}
       \equiv2^{u_v-r}=\chi2^v=2^v\pmod p.              \tag{11}
\]
Changing the top index by \(jM\) adds the usual \(2^j\) in the
odd row, yielding \(P_k(2Z)\). The leading children still have
ratio two and factor as \(1+2Y\). Every coefficient congruence
error gains an integer valuation step, beyond \(\delta_1<1\).
Consequently the actual pure arrays satisfy
\[
 \begin{split}
 \frac{E(Y)}{E_rY^r}
     &=A+(\epsilon+\eta)P_k(Z)+o_v(\delta_1),\\
 \frac{YO(Y)}{O_{r-1}Y^r(1+2Y)}
     &=B+(\epsilon+2\eta)P_k(2Z)+o_v(\delta_1),
 \end{split}                                             \tag{12}
\]
where \(A=S_k(Z)\) and \(B=S_k(2Z)\).

The quadratic perturbations also lie beyond this layer. Uniformly,
\[
 \beta-\alpha\ge e/2-L
   >1-\frac1{p-1}>\frac{p+3}{2p}=\delta_1,
 \qquad \alpha=L/2,\quad\beta=v(h).                       \tag{13}
\]
The final strict inequality is equivalent to \(p^2-6p+3>0\),
valid for \(p\ge7\). The coefficientwise bound from Section 2
therefore applies to every term with a quadratic factor.
At a common quadratic Gaussian zero, both right sides of (12)
must vanish up to \(o_v(\delta_1)\).

### The first compatibility is exactly the Wronskian

Use primes for differentiation in \(Z\), and write
\(W=B'A-A'B\). Identity (6) gives, over the rationals,
\[
 P_k(Z)=A+2ZA',\quad P_k(2Z)=B+2ZB',\quad
 B'P_k(Z)-A'P_k(2Z)=W.                                  \tag{14}
\]
No derivative of \(\epsilon\) or \(\eta\) is taken: multiply
the two value equations (12) by \(B'\) and \(A'\) and subtract.
Their first secondary contribution is exactly \(\epsilon W\).
Thus
\[
 (1+\epsilon)W
 =-\eta\bigl(B'P_k(Z)-2A'P_k(2Z)\bigr)+o_v(\delta_1).
                                                               \tag{15}
\]
At a common residue \(z\), the coefficient in parentheses is
\(-2zA'(z)B'(z)\), a unit by squarefreeness and \(z\ne0\).
Since \(1+\epsilon\) is a unit, (15) forces \(v(W)=\delta_1<1\).
On the other hand (12) gives \(v(A)\ge\delta_0\); the same
lifted-gcd argument as in Section 4 gives
\[
                  v(W)\ge\min(2\delta_0,1)=1.
\]
This contradiction excludes every positive-character final residue.
Together with Sections 2--4 and the internal-edge theorem, it
proves (1) for all exponents and both quadratic characters.

For \(k=4,p=89\), the common residue is \(z=64\), with
\(A'(64)=14\), \(B'(64)=21\) modulo \(89\). The coefficient
in (15) is \(-2\cdot64\cdot14\cdot21\equiv15\ne0\).
The gaps are \(45/89\) and \(46/89\), and the inter-row ratios
are one and two, for every \(e\ge2\). Thus \(9\cdot89^e\)
has no common Gaussian coefficient direction at any positive exponent.

The checker exhaustively inspects the actual even and odd coefficient
inventories through \(\delta_1\) at seven selected inputs, including
\((k,p,e)=(4,89,2),(4,89,3)\), prime seven at exponent four,
and a negative-character phase control. Sparse coefficient units
computed by Wilson recursion are compared with independent
term-by-term factorial tables in the bounded controls. An incorrect
next-layer ratio of one at prime 89 is rejected. These finite
controls support the written uniform digit proof; they do not
replace it or constitute independent mathematical review.

The unrestricted odd-start problem, starts without a dominant prime
as in (1), and unrestricted trinomial resultant positivity remain
separate open questions.
