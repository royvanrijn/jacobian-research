# Wide upper mark windows at most maximal smaller exponents

Let \(A=\max_i\min(a_i,b_i)\) and \(B=\max_i\max(a_i,b_i)\)
for an arbitrary ordinary two-monomial profile after removing
equal rows. Use the protected upper window \((u_N,T_N]\) of the
[factorial envelope theorem](MULTIVARIATE_UNIVERSAL_ENVELOPE_MARKS.md#an-exact-mark-window-in-every-profile)
at terminal \(N=n+1\ge3\). Write \(G_N=T_N-u_N\).

**Prime criterion.** A prime \(p\in(u_N,T_N]\) with
\(p>\max(N,nB)\) proves coprimality of the adjacent physical
moments, uniformly over all other data of the profile.

In fact every lower argument is strictly below \(p\) and every
upper argument is at most \(NB<2p\), since \(N\le2n\).
Also \(p>NA\), so it is above every smaller-end argument
in the upper parent. After stripping row factorials, the
only possible upper contribution at \(p\) is the literal
\(p\), with valuation one each time. The lower coefficients
are all units, including their binomial factors. The upper
ones are units exactly on the nonempty interval
\(\{k:E_N(k)<p\}\) of index width at most one. A primitive
integral common factor has unit endpoints since both divide
the corresponding lower endpoints by Gauss's lemma. Its
reduction has full degree, has unit constant term, and divides
the upper reduction \(z^k(c+dz)\) with \(c\ne0\). Its degree
is therefore at most one. The even moment is strictly positive
on the real line, which excludes a linear factor. The rational
Bezout identity transfers the conclusion to any
characteristic-zero field.

## A uniform density-one result across orders

Let \(S=\sum_i(b_i-a_i)^2/\max(a_i,b_i)\) and fix \(S_0>0\).
Suppose the consecutive prime gaps \(g_j=p_{j+1}-p_j\) obey
\[
                   \sum_{p_j\le x} g_j^2
                               \ll_\epsilon x^{1+\sigma+\epsilon}
                                                                    \tag{1}
\]
for a constant \(\sigma\ge0\) and every \(\epsilon>0\).
Fix \((1/2+\sigma)/2<\beta\le1/2\) and an integer \(J\)
at least the [absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).

**Theorem.** Outside a set of
\(O_{S_0,J,\beta,\epsilon}(X^{3/2+\sigma-2\beta+\epsilon})\)
integers \(A\in[X,2X]\), every ordinary profile with this
maximal smaller exponent and \(S\le S_0\) has coprime
adjacent moments at each terminal \(3\le N<J\) for which
\(G_N\ge A^\beta\). The same exceptional set works for
all row counts, gaps, and orientations. In particular, at
those heights every such profile with least unequal gap
at least \(A^\beta\) satisfies SFC(2) in all degrees.

Indeed \(d^2\le S_0(A+d)\) for every absolute gap \(d\),
which gives \(B\le A+S_0+\sqrt{S_0 A}\). Thus for \(X\)
large relative to \(S_0,J\), one has \(NA>nB\) for all
\(3\le N<J\). Every upper mark window is contained in
\([NA, NA+C_0\sqrt X]\), where \(C_0\) depends only on
\(S_0,J\). A failure of the prime criterion in a window
of width at least \(A^\beta\) means the window is covered
by a gap \([p_j,p_{j+1}]\) with \(g_j\ge X^\beta\).
For such a gap the condition that it intersect the possible
window interval forces
\[
                      p_j-C_0\sqrt X\le NA\le p_j+g_j.
\]
It charges at most \(O_{S_0,J}(g_j+\sqrt X)\) integer \(A\),
and only gaps with \(p_j\ll_{S_0,J}X\) can be relevant.
Summing over those with \(g_j\ge X^\beta\) and using (1)
bounds the exceptional count for this \(N\) by
\[
 O_{S_0,J,\epsilon}\!\left(
    X^{1+\sigma+\epsilon-\beta}
   +X^{3/2+\sigma+\epsilon-2\beta}\right).
\]
The second exponent dominates for \(\beta\le1/2\).
Summing over the finitely many terminals proves the uniform
bound. The envelope theorem gives \(G_N\ge\min_i|b_i-a_i|\)
at every terminal \(N\ge3\). Terminal two follows from
strict positivity of the even moment at the real zero of
the linear moment. Terminals at and above \(J\) follow
from the absolute cutoff.

In (1) one can use \(\sigma=1/4\) from the mean-square
estimate of Peck and Maynard, yielding \(\beta>3/8\).
[Maynard, Theorem 3.1](https://arxiv.org/pdf/1201.1787)
The accepted version of Stadlmann's Theorem 1 gives
\(\sigma=23/100\), yielding \(\beta>73/200\).
[Stadlmann, Theorem 1](https://arxiv.org/pdf/2212.10867v2)
For either bound, taking \(\epsilon<2\beta-1/2-\sigma\)
gives a density-zero exceptional set of \(A\) by summing
over dyadic intervals. This is a one-parameter density
statement for maximal smaller exponents with uniform
quantification over profiles at each allowable height.
All-degree closure at arbitrary heights or very narrow
upper windows requires a further argument.
