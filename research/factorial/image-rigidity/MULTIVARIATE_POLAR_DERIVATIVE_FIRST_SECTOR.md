# A polar derivative criterion for the first sector of ordinary binomials

Let \(a,b\) be distinct nonnegative integer exponent columns and write
\[
 P_m(z)=\sum_{k=0}^m\binom mk F_{m,k}z^k,\qquad
 F_{m,k}=\prod_i((m-k)a_i+kb_i)! .
\]
For \(N=n+1\) define the positive function on \(0\le x\le n\)
\[
 g(x)=\prod_i\prod_{h=1}^{a_i}
              \bigl(na_i+(b_i-a_i)x+h\bigr).
\]
Empty products are one. For an upper-half-plane common root with
argument \(0<\theta<2\pi/N\), set \(j=\lfloor\pi/\theta\rfloor\).
Then \(1\le j<n\) and
\[
                              g(n)<g(j).                         \tag{1}
\]
The same condition holds with the columns exchanged. In particular
there is no common physical root in that sector or its conjugate
if \(g(n)\ge g(n-1)\) in either orientation. This applies to all
coordinatewise comparable columns. The criterion also covers some
mixed orientations by the degree-uniform bound below.

## The polar derivative comparison

Take \(A=P_n,\ B=P_N\), and \(C=B-zB'/N\). The binomial identity
\[
               (1-k/N)\binom Nk=\binom nk,\qquad 0\le k\le n,
\]
shows that \(C_k/A_k=F_{N,k}/F_{n,k}=g(k)\), where \(A_k,C_k\)
are the coefficients of the two polynomials. Every affine factor
in \(g(x)\) is positive on \([0,n]\), so \(\log g\) is concave.

Suppose \(A(z)=B(z)=0\) at \(z=r e^{i\theta}\) with
\(r>0\) and \(0<\theta<2\pi/N\). The identity \(\Im A(z)=0\)
forces a negative sine among \(k\theta\) for \(k\le n\), since
\(\sin\theta>0\) and the coefficients are positive. This proves
\(1\le j<n\). If \(g(n)\ge g(j)\), concavity gives
\(g(k)\le g(j)\) for \(0\le k\le j\) and \(g(k)\ge g(j)\)
for \(j\le k\le n\). The sines are nonnegative before the
cut at \(j\) and negative after it. Consequently
\[
 \Im C(z)
      =\sum_{k=0}^n[g(k)-g(j)] A_k r^k\sin(k\theta)\le0.
\]
On the other hand, writing \(B_k\) for the coefficients of \(B\),
the equation \(B(z)=0\) gives
\[
 \Im C(z)=-\frac{1}{N}\sum_{k=0}^{N}
          (k-\pi/\theta)B_k r^k\sin(k\theta)>0 .
\]
Each summand is nonpositive because \(N\theta<2\pi\), and the
summand for \(k=1\) is strictly negative. This proves (1).
Column exchange replaces \(z\) by \(1/z\), and conjugation returns
the angle to the upper half plane. If \(g(n)\ge g(n-1)\),
concavity makes the sequence \(g(0),\ldots,g(n)\)
nondecreasing, so (1) fails. When \(a_i\le b_i\) for all \(i\),
each factor in \(g\) is nondecreasing. Reversal supplies the
other comparable orientation.

## A degree-uniform sufficient bound for mixed orientations

For this bound, discard the rows with \(a_i=0\), which contribute
one to \(g\), and the equal rows, whose factors are constant in
\(x\). Assume that every retained row has \(b_i>0\). Set
\[
 d_i=b_i-a_i,\qquad
 \Delta=\sum_i d_i,\qquad
 V_b=\sum_i\frac{d_i^2}{b_i},\qquad
 W_a=\sum_{d_i>0}\frac{d_i}{a_i},
\]
where all these sums are over retained rows. If
\[
                  \Delta\ge\frac54 V_b+\frac5{24}W_a,           \tag{2}
\]
then for every \(n\ge2\) the physical moments are coprime in
the sector \(0<|\arg z|<2\pi/(n+1)\).

In fact for any retained row its contribution to
\((\log g)'(n)\) is
\[
 d_i\sum_{h=1}^{a_i}\frac1{nb_i+h}
 =d_i\sum_{h=1}^{a_i}\frac1{na_i+h}
  -n d_i^2 \sum_{h=1}^{a_i}\frac1{(nb_i+h)(na_i+h)} .
\]
The last subtracted term is at most \(d_i^2/(n b_i)\).
For \(a>0\), compare the decreasing convex function
\(t\mapsto(na+t)^{-1}\) on \([0,a]\) with its right endpoint
and trapezoidal unit sums. The integral is \(L=\log(N/n)\),
and the trapezoidal sum exceeds the right sum by
\((1/(na)-1/(Na))/2=1/(2nNa)\). Thus
\[
     L-\frac{1}{2nNa}\le\sum_{h=1}^{a}\frac1{na+h}\le L .
\]
The error is needed only when \(d_i>0\). Summing over the rows gives
\[
       (\log g)'(n)
          \ge\Delta L-\frac{V_b}{n}-\frac{W_a}{2nN}.
\]
Since \(L\ge 2/(2n+1)\) and \(\Delta\ge0\), a sufficient bound
is
\[
                  \Delta\ge
              \left(1+\frac1{2n}\right)V_b
                 +\frac{2n+1}{4nN}W_a .
\]
Both coefficients decrease on \(n\ge2\), and at \(n=2\) they
give (2). Concavity then makes \(g\) nondecreasing throughout
\([0,n]\), so the polar derivative criterion applies.

For example, the fully shared columns \(a=(10,10,10)\) and
\(b=(11,11,9)\) have \(\Delta=1\) and
\((5/4)V_b+(5/24)W_a=323/792<1\). Their mixed orientations
therefore satisfy (2). The test concerns a first sector at the
physical radius. Separate control of other sectors is required
for full consecutive coprimality.
