# Fourth moments for groups of valuation layers

Use an anatomic box at height \(X\) with fixed layer
depth. For one letter write its complete positive
values as \(n=c\prod_i x_i^i\), with coefficient \(c\)
fixed and \(x_i\) of size at most \(O(X^{\alpha_i})\).
Put \(\alpha=\sum_i\alpha_i\) and \(A=X^\alpha\).
The usual divisor bound gives \(O(X^\varepsilon)\)
representations of each \(n\) at fixed depth, after
resizing \(\varepsilon\). Bounds remain valid on a
restricted collection of the allowed tuples.

## Sparse lifting on a grouped diagonal surface

Fix \(j\ge3\) and a set \(I\) of levels divisible by
\(j\). Write \(n=C(\mathbf z)v^j\), where
\(v=\prod_{i\in I}x_i^{i/j}\). Put
\[
        s=\sum_{i\in I}\alpha_i,\qquad
        h=\sum_{i\in I}(i/j)\alpha_i,\qquad V=X^h .
\]
The dictionary of possible roots \(v\) has cardinality
\(R\ll X^s\), and all roots are positive and \(O(V)\).
Each root has subpower multiplicity in the chosen
coordinates, again by the divisor bound. The number of
outside tuples \(\mathbf z\) is \(O(X^{\alpha-s})\).

Let \(S(\vartheta)\) be the additive exponential sum
on this letter with tuple multiplicity. Its fourth
moment counts \(n_1+n_4=n_2+n_3\). Terms with a cancelling
pair account for \(O_\varepsilon(X^\varepsilon A^2)\)
by the complete-value multiplicity bound, before the
outside tuples are fixed. After fixing the four
outside tuples, the other roots lie on a nonsingular
diagonal surface of degree \(j\) and no pairwise
subsum vanishes. Put \(\beta=3/\sqrt j\). If their
common divisor lies in \([G,2G)\), Salberger's bound
gives \(O_j((V/G)^\beta(1+\log V)^4)\) primitive
projective points. Each point has at most \(O(\min(G,R))\)
lifts: its first positive primitive coordinate fixes
the scaling as soon as its dictionary value is chosen.
For \(\beta\ge1\), the block costs at most
\(O_j(V^\beta(1+\log V)^4)\). For \(\beta\le1\), split
at \(G=R\); the bound for either block is at most
\(O_j(V^\beta R^{1-\beta}(1+\log V)^4)\).
This accounts for the gcd without replacing the
sparse dictionary size by its numerical height.
[Salberger, Corollary 0.7](https://doi.org/10.1112/plms.12508)

Summing over the outside tuples and the dyadic gcd
blocks gives
\[
 \int_0^1|S(\vartheta)|^4\,d\vartheta
        \ll_\varepsilon X^\varepsilon
         \left(X^{2\alpha}+X^{4\alpha-t_I}\right),
 \qquad t_I=(4-(1-\beta)_+)s-\beta h .                         \tag{GV1}
\]
Hence for fixed \(j\) the strongest bound obtainable
from these groups uses the positive coefficient levels:
\[
      t_j=\sum_{k\ge1}
            \big(4-(1-\beta)_+-\beta k\big)_+
                                      \alpha_{kj}.            \tag{GV2}
\]
A visible subsum is sufficient at finite depth. For
\(j\ge9\) these coefficients are
\((3-3(k-1)/\sqrt j)_+\), so several multiples of
a high degree can contribute. Each coefficient is
nonincreasing in \(\beta\) for \(k\ge1\), permitting
rational upper rounding of \(\beta\) in an exact
certificate.
[Zhang, Lemma 4.1 and equations (4.1)–(4.11)](https://doi.org/10.5281/zenodo.22151481)

## Permuted moment inequalities

For each letter \(r\), let \(\alpha_r\) be its total
mass, and let \(F_r\) be the largest available
divisible mass from the inherited BBLT selectors.
Their second moment costs \(X^{\alpha_r+\varepsilon}\)
and their fourth moment costs
\(X^{3\alpha_r-F_r+\varepsilon}\), uniformly in a fixed
box. Thus the best saving from \(4\alpha_r\) furnished
by those bounds and (GV2) is
\[
        Q_r=\max\big\{\alpha_r+F_r,\
                          \max_{j\ge3}\min(2\alpha_r,t_{r,j})\big\}.
\]
Orthogonality and Hölder with norms \(L^2,L^4,L^4\)
give, in every orientation,
\[
              D\le \alpha_A/2+\alpha_B+\alpha_C
                                      -(Q_B+Q_C)/4.          \tag{GV3}
\]
Averaging the three rows gives \(6D\le
5L-\sum_r Q_r\) for \(L=\sum_r\alpha_r\). The same
symmetric row follows by Cauchy–Schwarz between the
second and fourth moments on each letter and three
\(L^3\) norms. All the choices of indices and bounds
are fixed at the inner height limit.
[BBLT, Proposition 3.1 and (3.1)–(3.2)](https://arxiv.org/html/2410.12234v2#S3)

For a single exposed fifth layer the diagonal-size
bound is valid when
\(\alpha\le(4-3/\sqrt5)\alpha_5/2\), extending the
displayed \(157\alpha_5/120\) sufficient region in
Zhang. For a high-degree example, put three quarters
of a letter's mass in degree 11 and one quarter in
degree 13. The group \(j=11\) has
\(t_{11}=9\alpha/4\), yielding a fourth-moment saving
\(2\alpha\) instead of the inherited \(7\alpha/4\).
The improvement persists when such a group is used on
just one side of (GV3).

Even if every letter has the diagonal-size moment
\(Q_r=2\alpha_r\), (GV3) gives \(D\le L/2\).
On the equal-mass prime-degree model, any \(j\) meets at
most one occupied prime, with mass \(u/n\), and
\(t_j\le3u/n\). The inherited saving \(u+u/n\) is
therefore stronger for \(n\ge3\), and the mixed
Hölder rows permit \(D=2u-u/(2n)\).
[Precisely the formal profile used here](ABC_EXTERNAL_MOVING_ROWS_LIMIT.md)
