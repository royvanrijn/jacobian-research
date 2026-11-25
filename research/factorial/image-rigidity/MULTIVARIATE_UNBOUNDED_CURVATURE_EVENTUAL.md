# Uniform eventual coprimality at unbounded curvature

Let \(a,b\) be distinct ordinary exponent columns of arbitrary
finite dimension, and delete their equal rows. Put
\[
 S=\sum_i\frac{(b_i-a_i)^2}{\max(a_i,b_i)},\qquad
 V=\sum_{a_i b_i>0}\frac{(b_i-a_i)^2}{\min(a_i,b_i)}
       +\sum_{a_i b_i=0}\max(a_i,b_i).
\]
There exist absolute constants \(S_0,N_0\) such that, for
every terminal \(N\ge N_0\) with \(S\ge S_0\), the ordinary
factorial moment polynomials of orders \(N-1,N\) are coprime.
The profiles may have arbitrary private variables and total
degree difference. The constants here are existence bounds.

It suffices to prove that any sequence with \(N\to\infty\)
and \(S\to\infty\) has eventual coprimality. Indeed, failure
of the stated result would give counterexamples with both
\(N\) and \(S\) at least \(j\) for every integer \(j\).
The [explicit large-curvature criterion](MULTIVARIATE_HIGH_CURVATURE_ANGULAR_SLOPE.md)
already covers \(S\ge2\log N+5\). Work on the remaining
band along such a sequence. All estimates below are uniform
in the profile at a fixed term, and the profile is held fixed
when comparing its two orders. Write
\[
 n=N-1,\qquad D=N e^{-S/2},\qquad W=1+V,\qquad
 K=\lceil128(D+\log N)\rceil .
\]
The inequality \(V\le S(S+2)\) follows by writing a shared
row as \((l,l+d)\), where \(l\ge1\) and
\(d^2/(l+d)\ge d/l-1\). A private row contributes equally
to \(S,V\). Consequently
\[
          \frac{WK}{n}\to0,\qquad
          \frac{D+1}{n}\to0,\qquad K=o(n).                \tag{1}
\]
In the first limit the contribution from \(D\) is bounded
by a constant times \((1+S)^2 e^{-S/2}\), and the
contribution from \(\log N\) by \(O((\log N)^3/N)\).

## Coefficient halves and endpoint interpolation

For the two integer orders \(m=n,N\), write
\[
 C_{m,k}=\prod_i((m-k)a_i+kb_i)!,\qquad
 G_m=(C_{m,m}/C_{m,0})^{1/m},\qquad
 c_{m,k}=\binom mk\frac{C_{m,k}}{C_{m,0}G_m^k}.
\]
The normalized polynomial is \(T_m(z)=\sum_{k=0}^m c_{m,k}z^k\)
with \(z=G_m t\) for the original scalar variable \(t\).
Both endpoint coefficients equal one. The
[discrete row-curvature bound](MULTIVARIATE_HIGH_CURVATURE_ANGULAR_SLOPE.md#an-endpoint-coefficient-estimate)
gives
\[
                  c_{m,k}\le\binom mk
                   \exp\left(-\frac{S k(m-k)}{2m}\right). \tag{2}
\]
Put \(h=\lfloor m/2\rfloor\) and split
\[
 L_m(w)=\sum_{k=0}^{h} c_{m,k}w^k,\quad
 U_m(w)=\sum_{k=0}^{m-h-1} c_{m,m-k}w^k,\quad
 T_m(z)=L_m(z)+z^m U_m(1/z).
\]
Each half \(H_m\) has the expression of a left half in an
oriented profile, with a possible change of one in the upper
summation index. In either orientation and for \(K\le k\le m/2\),
\[
                         c_{m,k}\le8^{-k}                 \tag{3}
\]
for large \(N\). Indeed \(\binom mk\le(em/k)^k\) bounds
\(k^{-1}\log c_{m,k}\) by
\(1-\log w-(S/2)(1-w)\) at \(w=k/m\).
This is convex between \(K/m\) and \(1/2\).
At \(K/m\) the bound is at most
\(\log(eD/K)+SK/(2m)<-\log8\) for large \(N\),
using (1), and its value at \(1/2\) tends to \(-\infty\).

If a starting exponent is zero in the orientation for
a half, then \(H_m(w)=1+o(1)\) on \(|w|\le2\). To see
this, choose a row \((0,v)\), \(v\ge1\), and put
\(f_{m,k}=k!/(m!)^{k/m}\). Its normalized factorial is at
most \(f_{m,k}^v\), since \((vk)!/(k!)^v\) is log-convex
in \(k\) and starts at one. Log-convexity follows from the
increasing consecutive ratios
\(\prod_{j=1}^v(v-(v-j)/(k+1))\).
Now \(\binom mk f_{m,k}\le\exp(k)\) and
\(f_{m,k}\le\exp(-k(m-k)/(2m))\) by the curvature bound
in the row \((0,1)\). The remaining rows have curvature
\(S-v\). Thus the coefficient at \(k\le m/2\) in this
orientation is at most
\(\exp(k-(S-1)k(m-k)/(2m))\le\exp(-(S-5)k/4)\).
A geometric sum on \(|w|\le2\) gives the assertion.

Consider now an oriented half with all starting exponents
\(a_i\ge1\). Write \(d_i=b_i-a_i\), signed, and
\(B_m=m C_{m,1}/(C_{m,0}G_m)\). At integer orders this is
\(c_{m,1}\), and \(B_m\le2D\) for large \(N\) by (2).
Use the same formula with positive Gamma values for real
\(m\in[n,N]\). The exact interpolation at the positive endpoint is
\[
 \begin{aligned}
 \frac{c_{m,k}k!}{B_m^k}&=e^{f_m(k)},\\
 f_m(x)&=\log\Gamma(m+1)-\log\Gamma(m+1-x)-x\log m\\
 &\quad+\sum_i\bigl[\log\Gamma(ma_i+1+d_i x)-\log\Gamma(ma_i+1)
       -x(\log\Gamma(ma_i+1+d_i)-\log\Gamma(ma_i+1))\bigr].
 \end{aligned}                                             \tag{4}
\]
In particular \(f_m(0)=f_m(1)=0\).
Let \(L=\max(1,\max_i|d_i|/a_i)\), so \(L\le W\).
All Gamma arguments on \(0\le x\le K\) have real value at
least \(ma_i/2\) or \(m/2\) by (1).
The polygamma series gives, for \(r\ge2\),
\[
 |\psi_{r-1}(z)|\le (r-1)!\left(z^{-r}
                  +\frac{z^{1-r}}{r-1}\right)
                  \le2(r-2)!z^{1-r}
 \quad(z\ge r-1).
\]
For \(2\le r\le K\) its use in (4) is valid by \(K=o(n)\),
also after one \(m\) derivative, when \(r\) is replaced by
\(r+1\). Since
\(\sum_i|d_i|^r/a_i^{r-1}\le V L^{r-2}\), it follows that,
uniformly for these \(r,x,m\),
\[
 |\partial_x^r f_m(x)|\le H_r,\qquad
 |\partial_m\partial_x^r f_m(x)|\le(2r/n)H_r,\qquad
 H_r=(r-2)!\left(\frac{2L}{n}\right)^{r-2}\frac{4W}{n}.
                                                               \tag{5}
\]
Indeed an \(m\) derivative increases the polygamma index
by one and brings a factor \(a_i\) in its row, or a factor
one in the binomial term.

Define \(\nu_r(m)=\Delta_x^r f_m(0)\) by forward differences
of step one. The integral for an \(r\)-fold forward difference
over \([0,1]^r\) bounds \(\nu_r,\nu_r'\) by the two bounds
in (5), since the sum of the integration coordinates
lies in \([0,r]\). In particular \(\max_{2\le r\le K}H_r=o(1)\)
by \(2LK/n=o(1)\). Put \(w_r(m)=e^{\nu_r(m)}-1\).
Then
\[
 |w_r|\le M_r,\qquad |w_r'|\le(2r/n)M_r,\qquad
 M_r=2H_r.                                                 \tag{6}
\]
For integer \(k\le K\), Newton's forward-difference identity
and \(f_m(0)=f_m(1)=0\) give
\[
             e^{f_m(k)}
             =\prod_{\substack{A\subseteq\{1,\ldots,k\}\\
                               |A|\ge2}}(1+w_{|A|}).      \tag{7}
\]

## A connected-hypergraph estimate

For \(k\le K\), let \(\Sigma_k(m)\) be the sum of weights
of connected hypergraphs on the labeled vertices \(1,\ldots,k\),
with allowed edges \(A\) of size at least two and edge weight
\(w_{|A|}(m)\). The weight of a hypergraph is the product
over its edges, and one vertex is connected with weight one.
The right side of (7) is the sum over all such hypergraphs.
Partitioning the vertices into connected components gives
the following exact coefficients of the formal logarithm
of \(H_m\) through order \(K\):
\[
                    \ell_k(m)=\frac{B_m^k}{k!}\Sigma_k(m).
                                                               \tag{8}
\]
For clarity, the coefficient of \(w^k/k!\) in the exponential
of the connected exponential generating series counts
unordered collections of connected graphs on disjoint
nonempty label sets. Every hypergraph has exactly one such
decomposition, including the isolated vertices. This gives
(8) from (4) and (7).

We bound the absolute connected sums using edge weights
\((r-2)!\alpha^{r-2}\beta\), where
\[
             \alpha=2L/n,\qquad \beta=8eW/n.
\]
These weights equal \(eM_r\), and \(\alpha\le\beta\).
Their absolute sum over all possible edges on \(k\) labels
is at most \(2\beta k^2\), since for \(\alpha k\le1/2\)
\[
 \sum_{r=2}^k\binom kr(r-2)!\alpha^{r-2}\beta
    \le\beta k^2\sum_{r=2}^k
                \frac{(\alpha k)^{r-2}}{r(r-1)}
    \le2\beta k^2.
\]
For any fixed collection of distinct edges, summing over
hypergraphs containing that collection therefore costs at
most a factor \(e^{2\beta k^2}\) times the collection's
weight.

Every connected hypergraph contains a collection of edges
that can discover the vertices starting at vertex one.
At each step choose an edge joining an old vertex to new
vertices, assign one old vertex of that edge as parent,
and take all new vertices of the edge as one block of
children of that parent. After all vertices are discovered,
the parent assignments form a rooted labeled tree, the
children of each vertex are partitioned into blocks, and
each block and its parent specify one edge together with
additional vertices from the same edge. Dropping the
requirement that those additional vertices were already
discovered only increases the sum of the weights of the
collections. For a block of \(s\) children, summing over
its additional vertices gives at most
\[
 \begin{aligned}
 \sum_{t\ge0}\binom kt\beta(s+t-1)!\alpha^{s+t-1}
 &\le\beta(s-1)!\alpha^{s-1}(1-\alpha k)^{-s}\\
 &\le(s-1)!(2\beta)^s .
 \end{aligned}
\]
An original collection has distinct edges, so the previous
factor \(e^{2\beta k^2}\) applies. In estimating the sum of
collection weights by tree and block choices we may include
redundant encodings because all the majorants are positive.

For \(u\) children at a given parent, summing over the
set partitions into blocks gives at most \(u!(2\beta)^u\).
Indeed the product of \((s-1)!\) over blocks counts choices
of cycles on those blocks, and together the partitions
and cycles are the permutations of the \(u\) children.
The factor \(\prod u!\) over a rooted labeled tree counts
the ways of ordering all child lists. With the root fixed,
the number of ordered rooted labeled trees on \(k\) labels
is at most \((k-1)!4^{k-1}\). For instance the depth-first
traversal records a word of length \(2(k-1)\) for the
unlabeled shape, and the other labels can be assigned in
\((k-1)!\) ways. Thus the sum of the connected weights with
edge weights \(eM_r\) is at most
\[
                  (k-1)!(8\beta)^{k-1} e^{2\beta k^2}
                  \le(k-1)!\left(\frac{64e^2 W}{n}\right)^{k-1}.
                                                               \tag{9}
\]
The second inequality holds for \(2\le k\le K\) for
large \(N\), since \(2\beta K\le1/2\) by (1).
The bound for \(k=1\) is immediate.

The same majorant controls the order derivative of the
signed connected sum. Differentiating the product over the
edges of a hypergraph, (6) gives the absolute bound
\((2/n)(\sum_A |A|)\prod_A M_{|A|}\).
If \(X=\sum_A |A|\), then \(X\le k e^{X/k}\), and
\(e^{|A|/k}\le e\). Hence \(|\Sigma_k'(m)|\) is bounded
by \(2k/n\) times (9).
This argument handles edge weights that vanish at a particular
order, because it uses the absolute product bounds and does
not divide by \(w_r\).

## Comparison of the endpoint logarithms

The remaining scalar derivative satisfies, uniformly for real
\(m\in[n,N]\) at a positive starting endpoint,
\[
                          |B_m'/B_m|\le2/n.               \tag{10}
\]
For the raw ratio in \(B_m/m\), each factorial ratio has
finitely many factors differing from \(ma_i\) by at most
\(|d_i|\). Differentiating the logarithm after removing
\(\sum_i d_i\log(ma_i)\) gives a total error \(O(V/m^2)\):
each factor \(\log(ma_i+j)-\log(ma_i)\) has derivative
\(-j/(m(ma_i+j))\), where \(|j|\le|d_i|\).

For the center, put \(F(x)=x\log x-x\) with \(F(0)=0\),
and \(E(x)=\log\Gamma(x+1)-F(x)\) for \(x\ge0\).
Log-convexity of Gamma and its consecutive differences
give \(0\le E'(x)\le1/x\) for \(x>0\).
The trigamma series and integral comparison give
\(|E''(x)|\le1/x^2\). Writing \(\Delta=\sum_i d_i\), we have
\[
 \log G_m=\Delta\log m+\sum_i(F(b_i)-F(a_i))
                +\frac1m\sum_i(E(mb_i)-E(ma_i)).
\]
A shared row contributes at most \(|\log(b_i/a_i)|\) to
the sum of errors. A private row of degree \(v\) contributes
at most \(1+\log(mv)\), by \(E(1)=1\) and the bound on
\(E'\). For the derivative of the error sum, a shared row
contributes at most \(2|\log(b_i/a_i)|/m\), by differentiating
\(yE'(my)\) in \(y\). A private row contributes at most
\(1/m\). Therefore the derivative of the last term is
\(O((V+S\log m)/m^2)\). Combining with the raw ratio shows
\((\log B_m)'=1/m+O((V+S\log m)/m^2)\), proving (10).
In particular \(B_m\le3D\) for real \(m\in[n,N]\) by its
bound at the integer order \(n\).

Let \(C=64e^2\) and \(p_m=C W B_m/n\). By (1), \(p_m=o(1)\).
From (8)–(10),
\[
 |\ell_k(m)|\le\frac{B_m}{k}p_m^{\,k-1},\qquad
 |\ell_k'(m)|\le\frac{4 B_m}{n}p_m^{\,k-1}
 \qquad(1\le k\le K).                                    \tag{11}
\]
Put \(A_m(w)=\sum_{k=1}^K\ell_k(m)w^k\).
For large \(N\) this gives \(\sup_{|w|\le4}|A_m(w)|\le8B_m\),
\(\sup_{|w|\le2}|A_m(w)|\le4B_m\), and
\[
                  \sup_{|w|\le2}|A_N(w)-A_n(w)|
                       =O(D/n)=o(1).                    \tag{12}
\]
For either integer order the coefficients of \(e^{A_m(w)}\)
agree with \(H_m(w)\) through degree \(K\). On \(|w|\le2\),
the tail of \(e^{A_m(w)}\) after degree \(K\) has modulus
at most \(2^{-K}e^{8B_m}\) by the coefficient estimates on
the circle \(|w|=4\). After division by \(e^{A_m(w)}\)
this is at most \(2^{-K}e^{12B_m}=o(1)\). By (3), the
original polynomial's tail after degree \(K\) has modulus
\(O(4^{-K})\), which remains \(o(1)\) after multiplication
by \(e^{4B_m}\). Therefore
\[
        H_m(w)=e^{A_m(w)}(1+o(1))\qquad (|w|\le2),         \tag{13}
\]
with a uniform analytic relative error.
For an endpoint with a zero exponent, (12) and (13) hold
by taking \(A_m=0\), as shown in the first section.

Each half is thus zero-free on \(|w|\le2\) for large \(N\).
Take the branch \(\mathcal H_m=\log H_m\) with zero value
at \(w=0\). Equations (12) and (13) give on this disk
\[
          \sup|\mathcal H_m|=O(D+1),\qquad
          \sup|\mathcal H_N-\mathcal H_n|=o(1)             \tag{14}
\]
for each orientation. These branches are real on the
positive real radius, by the real coefficients and positivity
of the halves there.

## Root radii and angles

Set \(r_*=\exp(64(D+1)/m)\), eventually below \(3/2\).
On \(|z|=r_*\), the term \(z^m U_m(1/z)\) has modulus at
least \(\exp(64(D+1)-8D-o(1))\), by (13) and \(B_m\le2D\)
at the integer orders. This exceeds
\(|L_m(z)|\le\exp(8D+o(1))\). The dominating polynomial
has all its zeros in \(|z|<r_*\), since \(U_m\) is zero-free
on the disk of radius two. Rouché's theorem for \(T_m\) and
for its reciprocal yields
\[
                    |\log |z||=O((D+1)/n)=o(1)           \tag{15}
\]
for every root \(z\).

Write the two half logarithms as \(\mathcal L_m,\mathcal U_m\).
For real \(\lambda,\theta\) with \(|\lambda|\) small, put
\[
 J_m(\lambda,\theta)=\mathcal L_m(e^{\lambda+i\theta})
                         -\mathcal U_m(e^{-\lambda-i\theta}).
\]
By (14) and Cauchy's derivative estimates on a strictly
smaller disk, the first and second derivatives of either
half logarithm are \(O(D+1)\) on, for instance, \(|w|\le3/2\).
The first derivative of the two-order difference there is
\(o(1)\). Consequently \(\partial_\lambda J_m\) and
\(\partial_\lambda\partial_\theta J_m\) are \(O(D+1)\),
and \(J_N-J_n\) and \(\partial_\theta(J_N-J_n)\) are
\(o(1)\) uniformly in a fixed neighborhood of \(\lambda=0\)
and all real \(\theta\).

Suppose a common physical root exists. After conjugation,
its two normalized roots have the same angle
\(\theta\in(0,\pi]\), since \(G_m>0\) and positive real
roots are impossible. Write their log radii as \(\lambda_m\).
The exact root equations give
\[
       m(\lambda_m+i\theta)
            =J_m(\lambda_m,\theta)+(2j_m+1)i\pi,
            \qquad j_m\in\mathbb Z .
\]
Subtracting their real parts and using (15) yields
\[
 n|\lambda_N-\lambda_n|
 \le |\lambda_N|+o(1)+O(D+1)|\lambda_N-\lambda_n|.
\]
Since \(D+1=o(n)\), this implies
\(n|\lambda_N-\lambda_n|=o(1)\). For any \(\lambda\) under
consideration \(J_m(\lambda,0)\) is real. Integrating the
angular derivative of the order difference at fixed radius
gives an imaginary difference \(o(1)\theta\). Changing the
radius from \(\lambda_n\) to \(\lambda_N\) adds at most
\(O(D+1)|\lambda_N-\lambda_n|\theta=o(1)\theta\), by
integrating the mixed derivative starting at \(\theta=0\).
Subtracting the imaginary parts of the root equations
would therefore give \(2\pi(j_N-j_n)=\theta+o(1)\theta\),
which is strictly between zero and \(2\pi\) for large \(N\).
This is impossible.

Each original moment polynomial has integer coefficients.
Complex coprimality gives a rational Bézout identity, which
also holds over every characteristic-zero field.
