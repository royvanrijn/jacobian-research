# Higher tanh moments conditioned on skew

Let \(0\le a_i\le A\) be a finite list with \(\sum_i a_i\le C\),
and assign signs \(\epsilon_i\in\{-1,1\}\). Assume \(A,C>0\) and
write \(P=\tfrac12\sum_i\epsilon_i a_i^2\). If either cap is zero,
all moments vanish. The list may be the normalized tanh factors
\(a_i=T_i/h\) in the Cayley model at a positive deficit \(h\).

For \(r\ge2\) define
\[
 \begin{aligned}
 m&=\lfloor C/A\rfloor,&
 D_r&=mA^r+(C-mA)^r,\\
 S_+(P)&=\frac{C^2+2P}{2C},&
 S_-(P)&=\frac{C^2-2P}{2C},\\
 H_r(S)&=\min\{A^{r-2}S^2,\ S^r\}.
 \end{aligned}
\]
Then the contributions of each orientation satisfy
\[
        \sum_{\epsilon_i=\pm1}a_i^r\le H_r(S_\pm(P)).
\]
In particular, for \(V_r=\sum_i a_i^r\) and
\(M_r=\sum_i\epsilon_i a_i^r\),
\[
 \begin{aligned}
 V_r &\le \min\{D_r,\ H_r(S_+(P))+H_r(S_-(P))\},\\
 -\min\{D_r,H_r(S_-(P))\}
     &\le M_r\le \min\{D_r,H_r(S_+(P))\}.
 \end{aligned}
\]
For \(P\in[P_-,P_+]\), one can use \(P_+\) in \(S_+\) and \(P_-\)
in \(S_-\). All upper bounds are rational at integer \(r\) when the
caps and skew interval endpoints are rational.

Let \(K_\pm=\sum_{\epsilon_i=\pm1}a_i^2\) and
\(s=\sqrt{K_+}+\sqrt{K_-}\). The triangle inequality within
each orientation gives \(s\le C\), while
\(K_+-K_-=2P\) gives \(s\ge\sqrt{2|P|}\).
For \(s>0\),
\[
                \sqrt{K_\pm}=\frac12\left(s\pm\frac{2P}{s}\right).
\]
Both functions on the right are nondecreasing for
\(s\ge\sqrt{2|P|}\). Hence \(\sqrt{K_\pm}\le S_\pm(P)\).
Within either orientation, termwise comparison bounds the \(r\)-th
moment by \(A^{r-2}K_\pm\), and comparison of the \(\ell^r\) and
\(\ell^2\) norms bounds it by \(K_\pm^{r/2}\).
The case \(s=0\) has only zero factors.

The independent bound \(D_r\) follows by concentrating the mass.
If two positive entries are below \(A\), replace the pair
by entries with the same sum and one at \(A\) or zero.
Convexity of \(t^r\) shows that this does not decrease the power sum.
Adding an entry can bring the sum to \(C\) if necessary.
The resulting list has \(m\) entries at \(A\) and at most
one fractional remainder.

At \(P=0\) this gives \(V_r\le2(C/2)^r\). For \(A\ge C/2\)
that bound is attained by a single equal pair with opposite signs.
The estimate reduces the odd-moment slack near zero skew,
where the usual bound \(C A^{r-1}\) can be much larger.
The signed and unsigned bounds can be used with the
[matched skew interval](MULTIVARIATE_TANH_EXCESS_LOSS_TRANSPORT.md)
at the next order. They also apply pointwise along the
[factorial-to-Gaussian annulus path](MULTIVARIATE_FACTORIAL_GAUSSIAN_ANNULUS_HOMOTOPY.md),
using its surviving Gauss factors and the total deficit at
each order.
