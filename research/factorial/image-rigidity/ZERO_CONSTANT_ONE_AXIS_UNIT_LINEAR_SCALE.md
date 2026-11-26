# The valuation-one boundary for one-axis quadratic windows

Let \(\Lambda_2(x^i y^j)=i!j!\) and
\[
 f=x^2/2+\beta xy+\gamma y^2/2+dx+ey,\qquad
 \bar\beta=1,\quad\bar\gamma=\bar d=\bar e=0 .
\]
Work in a finite extension of \(\mathbb Q_2\), with \(v(2)=1\);
bars denote reduction in the residue field.

**Theorem.** If \(v(e)=1\), the phase has no five consecutive zero
positive factorial moments. In particular, the exceptional
valuation-one case of the [symmetric approximation](ZERO_CONSTANT_ONE_AXIS_SYMMETRIC_APPROXIMATION.md)
is excluded over ramified as well as unramified fields.

## Precision forced by the three boundary equations

Suppose five zeros start at \(n\) and set \(j=n+2\).
The symmetric approximation gives \(n\) odd, \(v(\gamma)=v(d)=1\),
\(\overline{\gamma/2}=1\) and \(v(\beta-1)=v(d-e)=j\).
Put \(u=\beta-1,\ t=d-e\).
With the notation of the [weighted axis equations](ZERO_CONSTANT_ONE_AXIS_UNIT_POLE_REMOVAL.md#exact-interior-and-weighted-axis-equations),
write \(A_i,B_i\) for the normalized axis moments and set
\[
 r_i=A_i/A_{i-1},\qquad
 \rho=(\gamma d-\beta e)/(\gamma-\beta),\qquad
 z_i=(B_i+\rho B_{i-1})/A_{i-1}.
\]
Here \(A_i\) are units reducing to one, and \(v(B_i)\ge i\)
by the axis recurrence. Since \(v(\rho-e)=j+1\), we have
\(v(z_i)\ge i\). The three boundary equations read
\[
 u r_i+(t+ud)=(\gamma-\beta)z_i,\qquad i=j,j+1,j+2 .
\]
Writing \(q_i=r_{i+1}-r_i\), the exact recurrence
\(q_1=2-d/(1+d)\), \(q_i=2-d^2q_{i-1}/(r_{i-1}r_i)\) for \(i\ge2\)
shows \(v(q_i-2)\ge3\) for \(i\ge2\). Thus taking second
differences of the boundary equations gives
\[
 v(z_{j+2}-2z_{j+1}+z_j)\ge j+3,\qquad v(z_j)\ge j+2.
\]
Replacing \(\rho\) by \(e\) in the numerator of \(z_j\)
costs valuation at least \(2j\ge j+2\). Consequently
\[
                     v\bigl((B_j+eB_{j-1})/e^j\bigr)\ge2.    \tag{1}
\]

The equation at \(i=j+2\) also gives
\(v(u+t)\ge j+1\), since
\(u+t=(\gamma-\beta)z_i-u(r_i+d-1)\) and
\(v(r_i+d-1)\ge1\) by the axis recurrence. Let \(P_i^0\)
be the normalized moments of
\(f_0=x^2/2+xy+\gamma y^2/2+e(x+y)\).
The natural coefficient polynomials for the moments are integral.
Their derivatives in \(\beta,d\) at \(f_0\) are both one modulo
the integer two, for positive \(i\). Indeed \(\gamma,e\in2\mathcal O\),
and both derivatives of the universal matching series at
\((\beta,\gamma,d,e)=(1,0,0,0)\) equal \(T/(1+T)\) modulo two.
The linear Taylor part is therefore \(u+t\) to precision
\(j+1\), and the remaining terms have valuation at least \(2j\).
The exact zero of the original moment at \(j-1=n+1\) gives
\[
                              v(P_{j-1}^0/e^{j-1})\ge2.       \tag{2}
\]

## The axis condition fixes the parameter through integer-two precision

Use \(g=\gamma/e,\ a=2/e,\ b_i=B_i/e^i\),
\(c_i=b_i+b_{i-1}\) at \(i>0\), and \(p_i=P_i^0/e^i\).
Both \(g,a\) are units with \(\bar g=\bar a\). The exact recurrences are
\[
\begin{split}
 b_0&=1,\quad b_1=g+1,\quad b_i=(2i-1)g b_{i-1}+b_{i-2},\\
 p_0&=1,\quad p_1=a+g+2,\quad p_i=i a p_{i-1}+p_{i-2}+c_i
 \qquad(i\ge2).
\end{split}                                                   \tag{3}
\]
Choose a residue root \(\xi\) with \(\xi+\xi^{-1}=\bar g\).
The [axis torsion condition](ZERO_CONSTANT_ONE_AXIS_UNIT_POLE_REMOVAL.md)
gives \(\xi^j=1\): the ratio of the two roots is \(\xi^2\), and \(j\)
is odd. Also \(\xi\ne1\).
Choose the root of unity \(\lambda_0\) lifting \(\xi\) in a finite
unramified extension of \(\mathbb Q_2\), and work in its compositum
with the coefficient field. Its order is the order of \(\xi\).
Set \(g_0=\lambda_0+\lambda_0^{-1}\).
Let \(S_k(g)\) be the integral polynomial
\(\lambda^k+\lambda^{-k}\) when \(\lambda+\lambda^{-1}=g\).

Modulo four, the transition on \((b_{2m},b_{2m+1})\)
has matrix
\[
                  \begin{pmatrix}1&-g\\g&1-g^2\end{pmatrix}.
\]
Its trace is \(2-g^2\) and its determinant is one. Both parity
subsequences of \(c_i\) therefore satisfy the two-step recurrence
with this trace. For the even subsequence the first two values are
\(2-g^2\) and \(g^4-4g^2+2\) modulo four. Hence
\[
                       c_{2m}\equiv(-1)^m S_{2m}(g)\pmod4
                       \qquad(m\ge1).                         \tag{4}
\]
For the odd subsequence, put
\[
               (-1)^m c_{2m+1}\equiv S_{2m+1}(g)+2 H_m(g)
               \pmod4 .
\]
Here \(H_m\) is a polynomial modulo two. The odd initial values
\(c_1=g+2,\ c_3\equiv 2+g-2g^2-g^3\pmod4\)
give \(H_0=1,\ H_1=1+g+g^2\) modulo two.
The recurrence for \(H_m\) modulo two has trace \(g^2\)
and determinant one. At \(g_0\), its characteristic roots are
\(\xi^2,\xi^{-2}\), with difference \(\bar g^2\ne0\).
For \(M=(j-1)/2\), we have \(\xi^{2M}=\xi^{-1}\).
Solving the two-term recurrence at that value gives
\[
 H_M(\bar g_0)=
 \frac{\bar g(1+\bar g+\bar g^2)+\xi^3+\xi^{-3}}{\bar g^2}
 =1,
\]
where \(\xi^3+\xi^{-3}=\bar g^3+\bar g\).
At \(g_0\), the coefficients lie in the unramified extension,
so the residue equality implies \(H_M(g_0)\equiv1\pmod2\).
As \(S_j(g_0)=2\), this proves \(c_j(g_0)\equiv0\pmod4\).

The polynomial \(c_j(g)-(-1)^M S_j(g)\) has even integer
coefficients by the same congruence. If \(0<\delta=v(g-g_0)<1\),
choose the unit root \(\lambda\) with residue \(\xi\) of
\(Z^2-g Z+1\). The factorization of \(g-g_0\) as
\((\lambda-\lambda_0)(1-(\lambda\lambda_0)^{-1})\) shows
\(v(\lambda-\lambda_0)=\delta\).
Since \(j\) is odd, \(v(\lambda^j-1)=\delta\), and
\[
 S_j(g)-S_j(g_0)=(\lambda^j-1)^2/\lambda^j
\]
has valuation \(2\delta\). The change of the even-coefficient
correction has valuation at least \(1+\delta\), and \(c_j(g_0)\)
has valuation at least two. Thus \(v(c_j(g))=2\delta<2\),
contradicting (1). We conclude \(v(g-g_0)\ge1\).

## The even moment at its first nonzero digit

The even recurrence in (3) gives
\(\bar p_{2m}=\sum_{s=-m}^m\xi^{2s}\), using (4).
At each odd index the two terms \(\bar a\,\bar p_{2m}\) and
\(\bar c_{2m+1}=\xi^{2m+1}+\xi^{-2m-1}\) cancel, since
\(\bar a=\bar g=\xi+\xi^{-1}\). Thus all \(\bar p_{2m+1}=0\),
including \(\bar p_1=0\).
In the exact even recurrence the terms \(2m a p_{2m-1}\)
all have valuation strictly above one. Also
\(S_{2m}(g)=S_m(g)^2-2\), so \(g-g_0\in 2\mathcal O\) implies
\(S_{2m}(g)-S_{2m}(g_0)\in4\mathcal O\).
Equation (4) now gives
\[
 p_{j-1}\equiv
      \sum_{s=-M}^{M}(-1)^s\lambda_0^{2s}
 \quad\text{up to an error of valuation strictly above one.}
\]
The unsigned sum is zero because \(\lambda_0^{2j}=1\) and
\(\lambda_0^2\ne1\). The signed sum is therefore
\(-2\sum_{-M\le s\le M,\ s\ {\rm odd}}\lambda_0^{2s}\).
The last sum has unit residue: it is a geometric sum of ratio
\(\xi^4\ne1\), with \(M\) terms if \(M\) is even and \(M+1\)
if \(M\) is odd. Each count is positive and coprime to
\(j=2M+1\), so neither count is a multiple of the order
of \(\xi^4\). This gives \(v(p_{j-1})=1\), contradicting (2).
