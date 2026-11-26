# Nearly linear kernel counts for the entire diagonal Padé sequence

The diagonal Padé fractions for a three-halves power give compatible pairs of square-full kernels. Increasing their degree might appear to produce many more pairs. Under the balanced modulus and witness-height constraints, however, every admissible degree is eventually at most six. The entire sequence has a nearly linear number of pairs.

The result concerns the specified diagonal sequence. A pair admitting a different rational witness is outside this count.

## 1. The sequence and the main theorem

For each integer d>=1 put m=2d+1 and define U_d,V_d in Z[t] by

\[
 (w-1)^m(w^2+mw+1)
                  =V_d(w^2)-w^3U_d(w^2).                \tag{1}
\]

The coefficients of w and w^{m+1} on the left are zero. Thus its even part defines V_d of degree at most d, while its odd part is divisible by w^3 and defines U_d of degree d. We use the same notation for their degree-d homogenizations. In particular,

\[
 U_d(1,0)=-1,\qquad V_d(0,1)=-1.
\]

The first two members are

\[
 \begin{aligned}
 U_1(x,y)&=5y-x,&V_1(x,y)&=5x-y,\\
 U_2(x,y)&=-x^2+14xy+35y^2,&
 V_2(x,y)&=35x^2+14xy-y^2.
 \end{aligned}
\]

For integers x,y with nonzero U_d(x,y),V_d(x,y), let h_d be their positive gcd and put

\[
             u_d=U_d(x,y)/h_d,\qquad v_d=V_d(x,y)/h_d.
\]

These coordinates are signed integers with gcd one. Replacing them by their absolute values does not change the congruences below. Their reduced height is max(|u_d|,|v_d|).

Fix positive constants c_-,c_+. For real T>=2 and an integer N with

\[
                         c_-T^5\le N\le c_+T^5,          \tag{2}
\]

let W_N(T) count ordered integer pairs satisfying

\[
 T\le x,y<2T,\quad x\ne y,\quad (xy,N)=1,                 \tag{3}
\]

for which at least one d>=1 has nonzero reduced coordinates satisfying

\[
 \max(|u_d|,|v_d|)\le4T^2,\qquad
                         x^3u_d^2\equiv y^3v_d^2\pmod N. \tag{4}
\]

The degree may be chosen separately for each pair.

**Theorem 1.** For every epsilon>0,

\[
                         W_N(T)\ll_{\epsilon,c_-,c_+}
                                           T^{1+\epsilon}. \tag{5}
\]

For all sufficiently large T, depending only on c_-, every admissible degree is at most six. The subfamilies of degrees d=3,4,5,6 have bounds with exponents

\[
                         \frac{20}{21},\quad
                         \frac{17}{18},\quad
                         \frac{52}{55},\quad
                         \frac{37}{39},                 \tag{6}
\]

respectively, with an arbitrary added epsilon. The degree-two bound is T^{1+epsilon}; it is the largest of the estimates used in (5).

The theorem does not require x,y to be squarefree, or u_d,v_d to be units modulo N. Imposing these extra restrictions gives the same upper bound.

## 2. The polynomial identity

Write

\[
 H_m(w)=(w-1)^m(w^2+mw+1),\qquad c_m=m^2-4,
\]
and
\[
 Q_m(x,y)=x^2-(m^2-2)xy+y^2=(x-y)^2-c_mxy.
\]

Since m is odd,

\[
 H_m(-w)=-(w+1)^m(w^2-mw+1).
\]

Multiplying this expression by (1) proves

\[
 V_d(w^2)^2-w^6U_d(w^2)^2
       =-(w^2-1)^m\{w^4-(m^2-2)w^2+1\}.
\]

After homogenization this is

\[
 x^3U_d(x,y)^2-y^3V_d(x,y)^2
                         =(x-y)^mQ_m(x,y).              \tag{7}
\]

If x,y belong to the same positive interval [B,2B), then

\[
 x^2+y^2\le\tfrac52xy<(m^2-2)xy
\]

because m>=3. Thus Q_m(x,y)<0 and

\[
                         |Q_m(x,y)|\le4m^2B^2.          \tag{8}
\]

In particular the right side of (7) is nonzero whenever x!=y.

## 3. Cancellation and height, uniformly in the degree

**Lemma 2.** If a,b are coprime positive integers and

\[
                   h=\gcd(U_d(a,b),V_d(a,b)),
\]

then

\[
                         h\le2^m m(m-2)^{m/2}.          \tag{9}
\]

**Proof.** A prime dividing h divides neither a nor b: reduction at a=0 or b=0 uses V_d(0,1)=-1 or U_d(1,0)=-1.

Fix p|h, write e=v_p(h), and extend the p-adic valuation, normalized by v(p)=1, to a field containing a root w of w^2=a/b. Then w is a unit. Homogeneity shows that U_d(w^2),V_d(w^2) both have valuation at least e. Consequently

\[
                   e\le\min\{v(H_m(w)),v(H_m(-w))\}.     \tag{10}
\]

Set Q_+=w^2+mw+1 and Q_-=w^2-mw+1. First suppose p is odd. If w-1 and w+1 are units, (10) gives

\[
 e\le\min(v(Q_+),v(Q_-))
                         \le v(Q_+-Q_-)=v_p(m).         \tag{11}
\]

Suppose instead k=v(w-1)>0. Then w+1 is a unit. Unless p divides m-2, Q_- is a unit and e=0. If A=v_p(m-2)>0, then p does not divide m+2, so Q_+ is a unit and

\[
 e\le\min\{mk,v(Q_-)\},\qquad
                         Q_-=(w-1)^2-(m-2)w.
\]

When 2k differs from A, the valuation of Q_- is min(2k,A), at most A. When 2k=A, the first term in the minimum is mA/2. Thus e<=mA/2 in every case. The argument near w=-1 is identical, interchanging Q_+ and Q_-.

The odd prime divisors of h therefore lie among those of m(m-2). These two integers are coprime, and their total contribution is at most m(m-2)^{m/2}.

For p=2 the residue of w is 1, because a/b is a unit and the equation w^2=1 has only the root 1 in characteristic 2. Both Q_+ and Q_- are units. Hence

\[
 e\le m\min\{v(w-1),v(w+1)\}\le m,
\]

since their difference is 2. This bounds the 2-part by 2^m and proves (9). If w=1 or -1, the same minimum is read using the finite valuation of the other factor. QED.

Put

\[
 \alpha=1+1/\sqrt2,\quad A_0=\alpha/2,\quad
 \beta=A_0^2=\frac{3+2\sqrt2}{8},\quad
 C_0=\frac{A_0}{5(1+2\sqrt2)}.
\]

**Lemma 3.** For coprime positive a,b in [B,2B), the reduced height is at least

\[
 \frac{C_0}{\sqrt{m-2}}
                         \left(\frac{\beta B}{m-2}\right)^d. \tag{12}
\]

**Proof.** The positive real w=sqrt(a/b) lies in [1/sqrt2,sqrt2]. On this interval

\[
                         mw-w^2-1\ge m/5.               \tag{13}
\]

Indeed the left side is concave in w. At the first endpoint, division by m gives at least 1/sqrt2-1/2>1/5; at the other it gives at least sqrt2-1>1/5.

It follows that

\[
 |H_m(-w)|=(1+w)^m(mw-w^2-1)\ge(m/5)\alpha^m.
\]

Using H_m(-w)=V_d(w^2)+w^3U_d(w^2), w^3<=2sqrt2, and homogeneity gives

\[
 \max(|U_d(a,b)|,|V_d(a,b)|)
                 \ge\frac{m}{5(1+2\sqrt2)}\alpha^mB^d.  \tag{14}
\]

Divide (14) by (9) and use m=2d+1. This is exactly (12). QED.

For distinct x,y in [T,2T), write g=(x,y), a=x/g, b=y/g and B=T/g. Homogeneity cancels g^d from the reduced fraction. Thus a reduced height at most 4T^2 implies

\[
 g\ge
 \frac{\beta T^{1-2/d}}
 {(2d-1)\{(4/C_0)\sqrt{2d-1}\}^{1/d}}
                         \gg T^{1-2/d}/d,               \tag{15}
\]

with an absolute implied constant. The uniformity follows because (2d-1)^{1/(2d)} is bounded for all positive integers d, as is the fixed constant (4/C_0)^{1/d}.

## 4. Every admissible degree is logarithmic

For the primitive a,b of the preceding section, (7) is a nonzero integer divisible by h^2, where h is the gcd in Lemma 2. Since 0<|a-b|<B, (8) gives

\[
                         h\le2mB^{m/2+1}.               \tag{16}
\]

Divide (14) by (16). Because m=2d+1 and B<=T, the reduced height is at least

\[
              \frac{\alpha^m}{10(1+2\sqrt2)B^{3/2}}
                 \ge\frac{\alpha^m}{10(1+2\sqrt2)T^{3/2}}.
\]

The height bound therefore forces

\[
 m\le
 \frac{\log\{40(1+2\sqrt2)\}+(7/2)\log T}{\log\alpha}.    \tag{17}
\]

In particular

\[
                              d\ll\log T                \tag{18}
\]

with an absolute constant. The assumption x!=y is essential: at x=y, (7) is zero and the reduced fraction is 1 for every degree.

## 5. A product-divisibility estimate at any interval scale

The counting estimate used below must apply at the smaller primitive scale B=T/g, rather than only at the balanced scale.

**Lemma 4.** Fix d>=1 and m=2d+1. For any real Y>=2 and positive integer N, the number E_d(N,Y) of ordered integer pairs

\[
 Y\le x,y<2Y,\quad x\ne y,\quad(xy,N)=1,\quad
                         N\mid(x-y)^mQ_m(x,y)
\]

satisfies

\[
                 E_d(N,Y)\ll_{d,\epsilon}
                            N^\epsilon Y^2N^{-1/m}.     \tag{19}
\]

**Proof.** Put c=c_m and

\[
                         N_*=\frac{N}{(N,c^m)}.
\]

For a prime p write e=v_p(N), A=v_p(c), l=v_p(L), where L=x-y, and let e_* = max(e-mA,0)=v_p(N_*).

Allocate p^{e_*} to S if ml>=e_*, and otherwise allocate it to M. In the latter case l<=A. Indeed if l>A then Q_m=L^2-cxy has valuation exactly A, by the unit hypothesis; consequently

\[
 e_*=e-mA\le ml+A-mA\le ml,
\]

a contradiction. For l<=A in the latter case, product divisibility gives

\[
                     v_p(Q_m)\ge e-ml\ge e-mA=e_*.
\]

Thus the allocation gives a coprime full-prime-power split

\[
                         N_*=SM,\quad S\mid L^m,\quad M\mid Q_m. \tag{20}
\]

The fixed loss satisfies N/c^m<=N_*<=N. It is important here not to allocate every bad prime power to L. When 2l=A the valuation of Q_m can be arbitrarily large; the allocation above permits it to enter M.

For a fixed split put

\[
                         A(S)=\prod_{p^e\parallel S}p^{\lceil e/m\rceil}.
\]

Then A(S)|L and A(S)>=S^{1/m}. There are at most 2Y/A(S) nonzero multiples L with |L|<Y, with no additional constant.

For fixed L set r=x/y modulo M. Since y is a unit, this is defined, and M|Q_m becomes

\[
                         r^2-(m^2-2)r+1\equiv0\pmod M.  \tag{21}
\]

For fixed m there are O_m(2^{omega(M)}) roots, where omega counts distinct prime divisors. To see every exceptional prime explicitly, there are no roots modulo 2: the polynomial there is r^2+r+1. At an odd p, completing the square reduces the count to roots of

\[
                         z^2\equiv\Delta\pmod{p^e},
             \qquad \Delta=m^2(m^2-4)\ne0 .
\]

Write A=v_p(Delta). If e<=A, the number is p^{floor(e/2)}. If e>A, it is zero when A is odd, and at most 2p^{A/2} when A is even. For the latter bound remove the fixed valuation p^{A/2} from z; unit square roots have at most two residues modulo p, each lifting uniquely, and the removed valuation gives only p^{A/2} unused final digits. This proves the claimed bound with a fixed factor depending only on m.

For any root in (21),

\[
                              (r-1,M)\mid c.             \tag{22}
\]

Indeed r is a unit. If min(v_p(r-1),e)>v_p(c), then (r-1)^2-cr has valuation v_p(c)<e, contradicting (21).

Recover y from

\[
                              (r-1)y\equiv L\pmod M.
\]

If soluble, this is one class modulo M/(r-1,M). Equation (22) therefore gives O_m(1+Y/M) choices of y in [Y,2Y), even when L is not a unit modulo M.

A fixed split contributes at most

\[
 O_m\left[
  2^{\omega(M)}\frac{Y}{S^{1/m}}\left(1+\frac{Y}{M}\right)
 \right].                                               \tag{23}
\]

For a contributing pair, Q_m is nonzero and |Q_m|<=4m^2Y^2. Hence M<=4m^2Y^2 and S>=N_*/(4m^2Y^2). The non-divisor factor in (23) is consequently at most

\[
 \frac{Y}{S^{1/m}}+\frac{Y^2S^{1-1/m}}{N_*}
       \ll_m N_*^{-1/m}(Y^{1+2/m}+Y^2)
       \ll_m Y^2N^{-1/m}.                               \tag{24}
\]

There are at most 2^{omega(N)} possible splits. Their total divisor factor is at most O_m(4^{omega(N)}), which is O_{m,epsilon}(N^epsilon): for all sufficiently large primes, 4<=p^epsilon, and the finitely many smaller primes contribute a constant. This proves (19). QED.

## 6. Reduced congruences eliminate every degree above six

Now assume (2)-(4) and write g=(x,y), a=x/g, b=y/g and B=T/g. Since (g,N)=1, the actual reduced congruence implies

\[
                              N\mid a^3u_d^2-b^3v_d^2.
\]

This integer is nonzero. In fact reducing the homogeneous forms cancels g^d, and (7) expresses it as the nonzero integer

\[
 \frac{(a-b)^mQ_m(a,b)}
      {\gcd(U_d(a,b),V_d(a,b))^2}.
\]

Its absolute value is less than 256B^3T^4, using a,b<2B and |u_d|,|v_d|<=4T^2. Thus

\[
 B\ge(c_-/256)^{1/3}T^{1/3},\qquad
                         g\le(256/c_-)^{1/3}T^{2/3}.    \tag{25}
\]

Combine this with the uniform lower bound (15). For d>=7 it would give

\[
                              T^{1/21}\ll_{c_-}d.
\]

But (18) bounds d by an absolute multiple of log T. These inequalities are incompatible for sufficiently large T. Therefore all admissible degrees eventually belong to {1,...,6}.

Two reductions play different roles: multiplying a reduced congruence by a gcd gives a necessary unreduced product congruence; cancelling the kernel gcd modulo N in this section is permitted because the kernels are units. No inverse of the numerator-denominator gcd has been used.

## 7. Counting at the forced primitive scale

Fix d in {3,4,5,6}. Equation (15) forces

\[
                              g\ge G_d=c_dT^{1-2/d}
\]

with c_d>0 fixed. For each actual g, the primitive pair a,b lies in [T/g,2T/g), is a unit pair modulo N, and satisfies the product congruence in Lemma 4. Dropping primitivity only enlarges the set.

Equation (25) ensures T/g>=2 for sufficiently large T. Apply (19) at Y=T/g, sum over the allowed gcds, and then enlarge this positive sum:

\[
 \begin{aligned}
 \#\{\hbox{pairs at degree }d\}
 &\ll_{d,\eta}N^{\eta-1/m}T^2
                              \sum_{g\ge G_d}g^{-2}\\
 &\ll_{\epsilon,c_-,c_+}
                  T^{\,1+2/d-5/(2d+1)+\epsilon}.
 \end{aligned}                                         \tag{26}
\]

This gives precisely the four exponents in (6), whose maximum is 20/21.

At d=2, applying Lemma 4 directly at Y=T gives O(T^{1+epsilon}). At d=1 the same bound is O(T^{1/3+epsilon}). All higher degrees have already been excluded. Their finite union proves (5) for large T. The bounded range of smaller T is absorbed into the implied constant using the trivial O(T^2) pair count. This completes the proof of Theorem 1.

## 8. Relation to the full compatible-kernel problem

The degree-two member is the fraction used in the [balanced square-full kernel construction](POWERFUL_SUM_GRAPHS.md#square-full-kernel-pairs-at-the-balanced-scale). That construction supplies infinitely many moduli with a nearly linear number of admissible squarefree unit pairs. The upper bound here has matching order, up to subpower factors, inside this diagonal Padé sequence.

The unrestricted count K_N(T) permits arbitrary fractions of height at most 4T^2. Theorem 1 neither bounds those additional witnesses nor gives a pointwise saving for K_N(T). Transformed or composed Padé constructions with extra parameters also require their own analysis.

The polynomial identity, cancellation bounds, root counts and counting argument have all been derived above. No external Padé approximation theorem is required for the upper bound.
