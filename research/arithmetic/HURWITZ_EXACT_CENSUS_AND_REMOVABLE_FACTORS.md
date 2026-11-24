# Exact sizes and removable factors in the odd polynomial Hurwitz domains

Fix \(r\ge3\), \(N=2r+1\), and an \(N\)-cycle \(P\). Let \(X_D\) consist
of the ordered \(r\)-tuples of double transpositions with product \(P\).
Let \(X_T\) consist of the corresponding tuples with a 3-cycle in the
first position and double transpositions elsewhere. Write \(d=(N-3)/2\).
The absolute domains have sizes
\[
 |X_D/\langle P\rangle|=N^{r-2}d^r,\qquad
 |X_T/\langle P\rangle|=N^{r-2}d^{r-1}.
 \tag{1}
\]
In the pure domain, the exact fraction of absolute tuples whose last
factor has one exclusive leaf on each transposition is
\[
 \left(1-\frac{2}{N}\right)^{r-2}
 \left(1-\frac{2}{N-3}\right)^{r-1}
 \ \longrightarrow\ e^{-2}.                                          \tag{2}
\]
By an exclusive leaf we mean a sheet moved by this transposition and
fixed by every earlier factor. For \(r=3\), one can define the same
condition on the degree-five residual domain without a monodromy
restriction. The proof below also gives the labeled counts in degree
five when \(r=2\).

Equation (1) agrees with the Goulden–Jackson cactus count. The derivation
below verifies the cardinalities directly on the fixed-product spaces,
including the free quotient when \(N\) is composite. Equation (2)
gives the size of the directly removable locus in
\(\texttt{GAP-ARITH-HURWITZ-UNIFORM-REMOVABLE-DT-REPRESENTATIVE}\).
Whether every braid orbit meets it remains open. [Goulden–Jackson, Theorems 3.1–3.2 and the last line of the proof](https://www-users.cse.umn.edu/~reiner/Classes/Math8680Fall2014Papers/GouldenJacksonOnCacti.pdf)

## The class product as a finite difference

Set \(L=N-1\) and write \(\chi_j\) for the character of the \(j\)-th
exterior power of the standard \(L\)-dimensional representation of
\(S_N\), for \(0\le j\le L\). The generating polynomial on a
permutation \(g\) is
\(\prod_{c}(1-(-s)^{|c|})/(1+s)\), with \(c\) ranging over its cycles.
These characters are distinct irreducibles. Here is an elementary
inner-product check valid in every degree. Counting permutations with
given cycle multiplicities gives the joint cycle index
\[
 \exp\left(\sum_{h\ge1}\frac{z^h}{h}
                  (1-(-s)^h)(1-(-t)^h)\right)
 =1+\frac{(1+s)(1+t)z}{(1-z)(1-stz)} .
\]
Indeed a cycle count vector of mass \(N\) has probability
\(\prod_h(h^{m_h}m_h!)^{-1}\), and expanding the exponential sums these
vectors. For \(N\ge1\), extract \(z^N\) and divide by
\((1+s)(1+t)\) to obtain \(\sum_{j=0}^{N-1}(st)^j\).
This says the exterior-power characters are orthonormal. At \(P\),
\(\chi_j(P)=(-1)^j\). The \(N\) squares of these values exhaust
the column norm \(|C_{S_N}(P)|=N\), so every other irreducible
character vanishes at \(P\).

For classes \(C_i\), the central class sum acts on an irreducible
of dimension \(d_j\) by \(|C_i|\chi_j(C_i)/d_j\). Take the trace of
\(P^{-1}\) times the product of these sums in the regular
representation and divide by \(N!\) to obtain the coefficient of \(P\).
As \(d_j=\binom Lj\), the exact total is
\[
 \frac{\prod_i|C_i|}{N!}\sum_{j=0}^L(-1)^j\binom Lj
                 \prod_{i=1}^r\frac{\chi_j(C_i)}{\binom Lj}.           \tag{3}
\]
For a double transposition \(D\) the generating polynomial is
\((1+s)^L-4s(1+s)^{L-2}\). For a 3-cycle \(T\) it is
\((1+s)^L-3s(1+s)^{L-2}\). Therefore
\[
 \frac{\chi_j(C_i)}{\binom Lj}
 =1-\xi_i\frac{j(L-j)}{L(L-1)},\qquad
 \xi_i=4\ \text{for }D,\qquad \xi_i=3\ \text{for }T.
\]
The product in (3) has degree \(2r=L\) in \(j\).
The alternating binomial sum is the \(L\)-th finite difference,
since \(L\) is even. It is \(L!\) times the leading coefficient.
Using \(|D|=3\binom N4\), \(|T|=2\binom N3\) gives
\[
 \#\{g_1\cdots g_r=P:g_i\in C_i\}
  =\frac1N\prod_i\frac{|C_i|\xi_i}{L(L-1)}
  =N^{r-1}d^{\,\#\{i:C_i=D\}}.                                       \tag{4}
\]
This argument works for any placement of the 3-cycles. It uses ordinary
finite-group character orthogonality, also used in the repository's
[degree-thirteen census](DEGREE_FOURTEEN_HURWITZ_ORBIT.md#an-independent-exact-total-bound).

Every nonidentity power of \(P\) has orbits all of one odd length greater
than one. It cannot preserve a four-element DT support, since that length
does not divide four. Consequently the simultaneous conjugation quotient
of each of the two stated domains is free. Dividing (4) by \(N\) gives (1).

## Counting the directly removable tuples

Compose permutations right to left. Suppose the last factor has
exclusive leaves \(a,b\), let \(R\) be the product of the earlier
factors, and write \(P=R g_r\). Since \(R\) fixes \(a,b\),
the two transpositions must pair \(a\) with \(P^{-1}(a)\) and \(b\)
with \(P^{-1}(b)\). Thus the two leaves are nonadjacent in the cyclic
order of \(P\). Conversely for each nonadjacent unordered pair, these
transpositions define \(g_r\), and \(P g_r\) fixes \(a,b\) and cycles
through all remaining labels in their inherited cyclic order. Indeed,
for a predecessor of either leaf the product skips that leaf and
continues at its successor, and otherwise it follows \(P\). Each
factorization of this residual long cycle by \(r-1\) double
transpositions lifts by fixing \(a,b\) in the earlier factors.

The leaves in each such tuple are determined uniquely. If both
endpoints of a transposition were fixed by every earlier factor,
\(R g_r\) would have that transposition as a two-cycle, contradicting
the long-cycle product. There are \(\binom N2-N=Nd\) choices of leaf
pair. The resulting labeled count is therefore
\[
 Nd (N-2)^{r-2}\left(\frac{N-5}{2}\right)^{r-1}.
\]
Dividing by (4) in the pure case gives (2). Dividing both counts by \(N\)
shows that the same fraction holds after quotienting. For \(N=15\)
the absolute totals are \(212\,576\,400\,000\) pure and
\(35\,429\,400\,000\) marked mixed. The pure count with directly
removable last factor is \(34\,808\,718\,750\). Direct enumeration\nof the labeled domains in degrees 5 and 7 gives pure totals 5 and 392,\nmarked totals 5 and 196, and directly removable counts 0 and 70,\nagreeing with the formulas.

## Scope of the conjugate-power mechanism

Féray's [pinned result](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/results/feray-conjugate-power-word/PROOF.md)
uses \(\sigma\tau\sigma^2\tau^2=X^2Y^{-1}XY\) for \(X=\sigma\tau\),
\(Y=\tau\), a bijective change of two independent uniform permutations.
Given the conjugacy class of \(X\), the power \(X^a\) is uniform on
its image class by equivariance. Conditional on any particular \(X\)
in that class, \(Y^{-1}X^bY\) is uniform on its image class, a law
independent of the choice of \(X\). This proves the conditional
independence used by the class Fourier argument. The same Schur
averaging mechanism yields the exact coefficient (3) for independent
class variables before conditioning on their product.

The mixing estimate in the external theorem requires few cycles,
including fixed points. A double transposition has \(N-2\) cycles.
The product of two DTs moves at most eight labels, whereas the uniform
law on \(A_N\) assigns probability at most
\(2\sum_{j=0}^8\binom Nj j!/N!\) to that event, tending to zero.
The variation distance for these two classes therefore tends to one.
The exact census and removable-locus density do not imply that a
fixed-product braid orbit intersects the removable locus. In
particular the pure and marked mixed uniform connectivity obligations
remain open.

For verification of the external rate, its two analytic inputs are the
uniform character bound in Teyssier–Thévenin Theorem 1.7 and the
fixed-rank zeta bound (6.3) with rank cutoff one. If both classes have
\(O(\log N)\) cycles, put \(q=(\log\log N+O(1))/\log N\). Schur
averaging and \(\ell^2\) Plancherel bound four times the squared TV
by \(\sum_{d_\chi>1}d_\chi^{-(2-4q)}\). The zeta bound at exponent one
and the \(\ell^s\) norm inequality give \(O(N^{-2}(\log N)^4)\).
For a uniform permutation \(X\), the cycle insertion identity gives
\(\mathbb E 2^{\#\mathrm{cycles}(X)}=N+1\). Powers of fixed nonzero
exponent multiply the cycle count by at most the absolute exponent.
Markov's inequality and parity averaging therefore give
\(O((\log N)^2/N)\) for the exact word \(X^2Y^{-1}XY\).
The two analytic inputs have no role in (1) or (2).
[Teyssier–Thévenin, v4](https://arxiv.org/html/2411.04347v4)

The leading-rate supplement retains the standard character \(F-1\).
For \(N\ge3\), counting marked disjoint cycles gives
\(\mathbb E[(F(X)-1)(F(X^2)-1)]=1\), since
\(F(X^2)=F(X)+2C_2(X)\), \(\operatorname{Var}F=1\), and
\(\operatorname{Cov}(F,C_2)=0\). Its sign twist vanishes when two
unused labels remain. Rank-two character polynomials in \(F,C_2\)
have bounded product moments at \(X,X^2\), so their Fourier
coefficients are \(O(N^{-2})\). At rank cutoff three the same
zeta estimate bounds the high-rank \(L^1\) remainder by
\(O(N^{-3}(\log N)^6)\) on the logarithmic-cycle event.
Using the threshold \(6\log_2 N\), the complement has probability
\(O(N^{-5})\) and discarded low ranks have dimensions \(O(N^2)\).
Thus the leading TV is the mean of \(|F-1|/(2(N-1))\) up to
\(O(N^{-2})\), or \(e^{-1}/(N-1)+O(N^{-2})\) by inclusion–exclusion.
[Pinned leading-rate proof](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/results/feray-conjugate-power-word/SHARP-RATE.md)
