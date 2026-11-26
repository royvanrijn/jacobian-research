# The third open dyadic band for a bilinear three-term phase

Let \(f=\beta XY+dX+eY\) with all three coefficients nonzero, and
let \(P_n=\Lambda_2(f^n)/n!\), where
\(\Lambda_2(X^iY^j)=i!j!\). Over any finite extension of
\(\mathbb Q_2\), with \(v(2)=1\), the three-window holds at
every positive start if
\[
             \tfrac34< v(d)-v(\beta)<\tfrac78,\qquad
             v(d)=v(e).
\]
This theorem completes the third *open* band in
[dyadic residue symmetry](BINARY_BILINEAR_DYADIC_RESIDUE_SYMMETRY.md).
That proof already treats residue \(\overline{e/d}\ne1\).
Here the exact low-degree arithmetic is supplied by the
[finite coefficient checker](verify_binary_bilinear_third_dyadic_tie.py).
The reduction to it is uniform in the starting order and
the ramification degree.

## Remaining classes and the middle moment

Use \(b=\beta/d,\ q=e/d,\ s=-v(b)\) after dividing \(f\)
by \(d\), and retain \(P_n\) for the normalized moments.
A zero triple ending at \(k\) forces \(q^k=-1\). By the
parity layers of the preceding proof, it remains to
treat only \(\bar q=1\) and \(k\equiv0,\pm1\bmod16\).

For odd \(k\), one has \(q=-1\). Put \(k=2r+1\) and
\(X=b^2\). The preceding proof gives
\[
 E_r(X)=P_{2r}=\sum_{t=0}^r c_tX^t,\qquad
 c_t=(r+t)!/(r-t)!,\qquad
 P_{2r+1}/b=(r+1+X\partial_X)E_r .
\]
The residue classes of \(k\) imply \(r\equiv0,7\bmod8\).
Set \(a=r/2\) if \(r\) is even, or \(a=(r+1)/2\) if \(r\)
is odd, so \(a\) is a positive multiple of four.
For \(1\le t\le r\), the extra valuation
\(v(c_t/c_{t-1})-1\) is obtained from the sequence
\(a,a\pm1,a\mp1,a\pm2,a\mp2,\ldots\) by taking \(v_2\)
of each entry. This follows by factoring the unique even
integer of \((r+t)(r+1-t)\). The first seven extra valuations
are \(w,0,0,1,1,0,0\), where \(w=v_2(a)\ge2\).
Among \(1\le t\le7\), the slope \(v(c_t)/t\) is
uniquely least at \(t=7\), where it is \((9+w)/7\).
Beginning at \(t=8\), each block of eight extra valuations
is at least \((2,2,0,0,1,1,0,0)\), term by term where
the indices exist. Each prefix of that lower bound has
average at least \(3/4\). It follows that
\(v(c_t)-v(c_7)\ge(7/4)(t-7)\) for \(t\ge7\).
Since \(2s<7/4\), any least term of \(E_r(b^2)\)
with nonpositive valuation and positive index must
have index seven. Indeed, if \(v(c_7)/7>2s\) then
all terms with positive index have positive valuation.
Otherwise, for \(t<7\), their larger slope makes
\(v(c_t)-2st>v(c_7)-14s\), and for \(t>7\) the
suffix bound gives the same strict inequality.
A zero therefore needs the \(t=7\) term to cancel
the constant. Its coefficient in \(X E'_r(X)\) has
the same valuation, since \(7\) is odd, while every
other positive-index term has higher valuation.
This excludes a repeated zero and hence an odd terminal.

For even \(k\), put \(K=2^{v_2(k)}\ge16\). The element \(q\)
has order \(2K\) and \(\epsilon=1-q\) has valuation \(1/K\).
Write
\[
 P_n=\sum_{j=0}^n j!b^j p_{n,j},\qquad
 p_{n,j}=[T^n]F_j,\qquad
 F_j=\frac{T^j}{((1-T)(1-qT))^{j+1}},
\]
and use \(p_{n,j}=0\) when \(j>n\).
As shown by the residue at the middle index,
\(v(p_{k-1,j})=1-(j+1)/K\) for \(j<K\).
The cross-count cost is \(\phi(j)=v_2(j!)-js\), minimized
uniquely at \(j=7\). For \(j\ge K\) its lower bound is
\(\phi(7)+K(1-s)-1>\phi(7)+1\).
The middle comparison for \(j<K\) replaces \(s\) by
\(s+1/K\) and therefore has the unique minimizer \(j=7\)
or \(j=15\) except if
\[
                         s=\tfrac78-\tfrac1K .             \tag{1}
\]
There is enough room for both indices since \(K\ge16\).
The \(j=7\) contribution is below \(\phi(7)+1\), so the
tail \(j\ge K\) cannot cancel that unique minimum.
It remains to handle (1).

## A combination at the tie

At (1), consider
\[
        G(b)=P_k-(21-\epsilon)P_{k-1}/b .
\]
A zero triple forces \(G=0\). The coefficient of \(j!b^j\)
for \(j\ge0\) in \(G\) is
\[
                 g_j=p_{k,j}-(21-\epsilon)(j+1)p_{k-1,j+1}.
                                                               \tag{2}
\]
The coefficient of \(b^{-1}\) is the additional term
\(-(21-\epsilon)p_{k-1,0}b^{-1}\) and has positive
valuation by \(v(p_{k-1,0})=1-1/K\).
The constant \(g_0\) is a unit: \(p_{k,0}=(1+q)/(1-q)\)
is a unit and \(v(p_{k-1,1})=1-2/K>0\).
We show every positive-index summand in \(G\) has
positive valuation.

At (1) the cost is
\[
                  \phi(j)=j/8-s_2(j)+j/K .                \tag{3}
\]
Every \(j\ge32\) has cost at least
\(\phi(7)+32(1-s)-1>0\). Its shifted middle term in (2)
also has positive valuation, since the relevant cost is
\(\phi(j+1)+s>0\) and \(p_{k-1,j+1}\) is integral.
At \(K=16\) the same conclusion holds already for
\(j\ge16\), by \(\phi(7)+16(1-s)-1>0\).

## Orders \(K\ge64\)

Use the exact formula from the preceding proof
\[
 p_{n,j}=-\operatorname{Res}_{T=1}
           (T^{-n-1}-q^{n+1}T^{n+1})F_j\,dT .             \tag{4}
\]
For \(n=k\) the multiplier is \(T^{-k-1}+qT^{k+1}\);
replace it temporarily by \(T^{-1}+qT\). For \(n=k-1\)
replace \(T^{-k}+T^k\) by \(2\). With \(Z\) an indeterminate,
let \(N_j(Z)\) be defined by equation (8) of the
[preceding proof](BINARY_BILINEAR_DYADIC_RESIDUE_SYMMETRY.md#the-tie-at-s34-1k),
and put
\[
 S_\ell(Z)=\sum_{r=0}^{\ell}\binom{\ell}{r}\binom{\ell+r}{r}
                         (1-Z)^r Z^{\ell-r}.
\]
Taking residues in the local variable \(T-1\), the two
main parts of (4) are \((-1)^jN_j(\epsilon)/\epsilon^{2j+1}\)
and \((-1)^j2 S_j(\epsilon)/\epsilon^{2j+1}\), respectively.
Thus the main part of \(g_j\) is \((-1)^j L_j(\epsilon)/
\epsilon^{2j+3}\), where
\[
          L_j(Z)=Z^2N_j(Z)+2(21-Z)(j+1)S_{j+1}(Z).        \tag{5}
\]
For \(1\le j\le31\), each nonzero coefficient \(l_{j,i}\)
of \(Z^i\) in \(L_j(Z)\) satisfies
\[
    v_2(l_{j,i})+j/8-s_2(j)+(i-j-3)/K>0,\qquad K\ge64 .
\]
The checker verifies both endpoints \(1/K=0,1/64\)
for this affine inequality. It is exactly the weighted
lower bound for that term in (5), after using (3).

The exact residue corrections also have positive
weighted valuation. Write \(u=T-1\) and \(a=\log_2 K\).
The coefficient of \(u^h\) in \(T^{\pm k}-1\), at
\(1\le h\le j<K\), has valuation at least
\(a-\lfloor\log_2 j\rfloor\), by the binomial identity
\(\binom{x}{h}=(x/h)\binom{x-1}{h-1}\).
The correction in the terminal residue has valuation
at least \(a-\lfloor\log_2 j\rfloor-2j/K\), because
it starts at \(u^1\) and the local expansion of \(F_j\)
has denominator valuation at most \((2j+1-h)/K\) in
the term taken with multiplier degree \(h\).
The correction in the middle residue of index \(j+1\)
has multiplier \((T^k-1)^2/T^k\), which starts at \(u^2\).
In the relevant degrees \(h\le j+1\) it gains valuation
at least \(2(a-\lfloor\log_2(j+1)\rfloor)\). Its
denominator loss is at most \((2j+1)/K\).
The extra factors in (2) are integral. At \(1\le j\le15\)
the terminal and middle correction bounds are at least
\(3-30/K,\ 4-31/K\), respectively, and
\(\phi(j)\ge-17/8\).
For \(16\le j\le31\) the bounds are at least
\(2-62/K,\ 2-63/K\), with
\(\phi(j)\ge-9/8+23/K\) by the bit at place four.
In both ranges the totals are strictly positive for
\(K\ge64\). This proves the needed weighted bound on
\(g_j\) without evaluating any high-order moment.

## Orders \(K=16,32\)

Clear the common denominator \(\epsilon^{2j+3}\) in
(2) by using (4) in the local coordinate \(u=T-1\).
Its resulting numerator depends on \(k\) only through
integer binomials \(\binom{x\pm k}{h}\) with \(h\le j+1\),
where \(x\) is fixed once \(j\) is fixed. This follows by
expanding \(T^j\) times each multiplier in (4); the
factors \(q^{k}=-1\) have already been substituted.
The binomials are periodic modulo \(16\) under
\(k\mapsto k+512\) for \(h\le32\): every
\(\binom{512}{i}\) with \(1\le i\le32\) is divisible
by \(16\), so the binomial convolution proves the claim
for any integer upper argument. Negative shifts follow
from the same periodicity. The resulting error after
dividing by the denominator and adding the weight
is at least \(4-(2j+3)/K+\phi(j)\). This is positive
for \(K=16,\ 1\le j\le15\), and for
\(K=32,\ 1\le j\le31\), using
\(\phi(j)\ge-17/8+7/K\).

It suffices, therefore, to verify the residue classes
\(k\bmod512\) with \(v_2(k)=\log_2 K\).
The checker computes (2) directly from
\[
 p_{n,j}=\sum_{r=0}^{n-j}\binom{j+r}{j}\binom{n-r}{j}q^r
\]
modulo \(q^K+1\). After converting to the \(\epsilon\)-basis,
the valuations of the nonzero terms are distinct modulo
\(1\), since their exponents in \(\epsilon\) are below \(K\).
The least such valuation is therefore exact. The checker
verifies the needed strict weighted inequality for
every specified class and cross count. Together with
the tail and constant valuations this excludes \(G=0\)
at (1) and completes the third open band.

Run the arithmetic certificate from the repository root:

    python3 research/factorial/image-rigidity/verify_binary_bilinear_third_dyadic_tie.py

It uses only the Python standard library. It checks
589 polynomial coefficient bounds, 488 cyclotomic
controls, and one deliberately incorrect calibration.
The dyadic thresholds bounding this open band and
higher valuation bands are outside this theorem.
