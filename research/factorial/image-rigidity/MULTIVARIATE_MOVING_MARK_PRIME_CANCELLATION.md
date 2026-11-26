# Exact cancellation of moving factorial mark primes

The [factorial envelope marks](MULTIVARIATE_UNIVERSAL_ENVELOPE_MARKS.md)
give a necessary prime-overlap budget for a common factor. The
budget can be arbitrarily large even at bounded positive curvature.
Rows below the envelope can simultaneously make the two Newton
polygons horizontal at many moving mark primes. An additional
prime then gives all-degree separation for the same profiles.

## Simultaneous cancellation theorem

**Theorem.** Fix orders \(4\le n<N=n+1\), a finite ordinary core,
a positive target curvature \(\tau\), a target signed bias
\(\lambda\in[-\tau,\tau]\), an integer signed added degree gap
\(\Delta\), and \(0<\beta<1/2\). For each sufficiently large integer
\(A\), put \(D=\lfloor A^\beta\rfloor\). There is a finite ordinary
extension of the core with the following properties.

- The added rows have absolute curvature tending to \(\tau\),
  signed bias tending to \(\lambda\), largest relative gap tending
  to zero, and signed degree gap exactly \(\Delta\).
- At both orders the entire factorial-argument envelope is that
  of the single mirrored pair of rows \((A,A+D),(A+D,A)\). Each
  protected window has length \(D\).
- In each window there is a fixed-step subprogression of the
  permissible marks \(x\equiv1\bmod j\). At every prime \(p>CD\)
  dividing any of these selected marks in either parent, both
  stripped valuation arrays are constant in the coefficient index.
  The sets of these primes in the two windows are disjoint and
  each contains \(\gg D\) primes. Here \(C\) is fixed as \(A\)
  varies.
- For the full \(x\equiv1\bmod j\) progression of either window,
  the sum of off-mark prime loads in the necessary budget of the
  envelope mark theorem, divided by the logarithm of its mark
  product times its minimum outside multiplicity, tends to
  infinity. The two polynomials still have disjoint root
  arguments in every consecutive degree over \(\mathbb C\).
  Coprimality follows over every characteristic-zero field.

The signed bias convention here is \(\sum_i (b_i-a_i)/\max(a_i,b_i)\)
over the added rows. The midpoint convention has the same limit.
The theorem concerns the valuation tests. The reductions of the
units at these moving primes are not specified.

## The protected windows and low-height rows

Let \(j=n,N\), and put \(s_j=\lfloor j/2\rfloor\) and
\(q_j=\lceil j/2\rceil\). For the mirrored top pair the envelope at
index \(k\) is \(jA+D\max(k,j-k)\). Its minimum is \(jA+Dq_j\).
A triple of consecutive indices can have both endpoint ordinates
at most \(jA+D(q_j+1)\), and cannot have both at most the minimum.
Thus the protected window is
\[
             \mathcal W_j=(jA+Dq_j,\ jA+D(q_j+1)] .             \tag{1}
\]
At a mark \(x\) in this window the increasing row first contributes
\(x\) at index \(q_j+1\), and the decreasing row contributes it
through index \(s_j-1\).

Take \(T=4nN+1\) and a fixed prime \(L>N+1\). In (1) select
the marks \(x=T y\) with
\[
           y\equiv1\pmod{nN},\qquad
       (y-1)/j\equiv\epsilon_j\pmod L,\qquad
                      \epsilon_n=0,\quad\epsilon_N=1 .          \tag{2}
\]
This is a single residue class in \(y\bmod nNL\) and is
coprime to that modulus. Indeed \(L\) is prime to \(nN\)
and the \(L\)-residue is \(1+j\epsilon_j\) in \([1,L-1]\).
The selected marks satisfy \(x\equiv1\bmod j\).

For a mark \(x=T y\) at order \(j\) and an integer
\(1\le h<n\), adjoin a reflected unit-gap pair
\[
       (b,b+1),(b+1,b),\qquad b=h(y-1)/j.                        \tag{3}
\]
Its marked factorial literal at order \(j\) is \(hy=jb+h\).
All these bases are comparable to \(A\), and for large \(A\)
their larger exponents are less than \(A\), since \(T>4N\).
Any additional row with larger exponent less than \(A\) has
every factorial argument strictly below \(jA\). Such rows
preserve both full envelopes and every literal count in (1).

**Isolation lemma.** If \(p>CD\) is a prime factor of a selected
mark \(x'\) at order \(r\in\{n,N\}\), its only occurrences in
the original top ramps of both orders are the literal \(x'\)
in order \(r\). In all the ramps from (3), it occurs only in the
pairs assigned to this same order and mark, at that order's
literal \(h x'/T\). Its valuation there is \(v_p(x')\).

For the top ramps compare a factor \(z=j'A+\ell\) with the mark
using \(rz-j'x'\). Their normalized arguments both lie in
\([A,A+D]\), so the determinant has absolute value at most
\(N^2 D\). If \(p\mid z\), the determinant is zero for
\(C>N^2\). Since \(x'\) is coprime to \(r\) and the two orders are
coprime, equality requires \(j'=r\) and \(z=x'\).
For a factor \(z=j'b+\ell\) of a row in (3) compare instead
\[
 rTz-j'h x'
          = rj'h\bigl(x/j-x'/r\bigr)+rT(\ell-j'h/j).
\]
Its absolute value is at most \(N^3 D+TN^2\). At a common
prime \(p>CD\), for \(C>N^3+TN^2\), it must be zero. Reducing
the exact equality modulo \(r\) gives \(r\mid j'hx'\).
Because \(h<n\), \(\gcd(r,x')=1\), and the orders are
coprime, this requires \(j'=r\). It follows that \(z=h x'/T\).
The congruences \(x'/T\equiv1\bmod r\) and \(1\le\ell\le r\)
give \(\ell=h\). Hence \(b=h(x'/T-1)/r\), which with (2)
implies \(r=j\) and \(x'=x\). Taking \(A\) large enough
makes \(p\) coprime to \(T h\), giving the valuation assertion.
This proves the lemma, and also shows that no such prime
divides two selected marks. A fixed core has no such primes
in its stripped ramps for large \(A\).

## Exact cancellation and arbitrary off-mark load

Use \(R\ge1\) copies of the reflected top pair. For every mark
selected at \(j\) add \(R\) pairs of (3) with \(h=s_j\).
At a prime assigned to that mark, divide the stripped
valuations by \(e=v_p(x)\). The top pair contributes
\[
             \mathbf1_{k\ge q_j+1}+\mathbf1_{k\le s_j-1}
\]
per copy. The corresponding pair in (3) contributes
\(\mathbf1_{k\ge s_j}+\mathbf1_{k\le j-s_j}\).
Their sum is exactly two for every index. The unmarked
parent has zero contributions at this prime by isolation.
The binomial coefficients are units at these primes.

To add arbitrary extra load without changing any slopes,
put \(h_1=\lfloor(j+1)/2\rfloor\) and
\(h_2=\lceil(j+1)/2\rceil\). Both are strictly less than \(n\)
because \(n\ge4\). Add \(M R\) copies of each of the two
pairs in (3) with these offsets for each selected mark, taking
\(M=\lfloor\sqrt D\rfloor\). If the offsets agree the two
blocks are simply added with both multiplicities. Since
\(h_1+h_2=j+1\), the two valuation profiles sum to a constant
of \(2e\) at that prime. This is the exact complement of the
two indicator steps at each offset. The other parent is still
unaffected. Each of the four unit-gap rows contributes at
least \(e\) to the off-mark load of its corresponding selected
mark. Its literal \(h_i x/T\) is strictly smaller than \(x\),
so omitting \(x\) from that load leaves it present.

For clarity, the number \(Q_j\) of selected marks is
\(D/(T nN L)+O(1)\). Work first with their values \(y=x/T\)
and put \(P_j^{\rm sel}=\prod y\). At a prime \(p\le CD\)
dividing this progression, remove one term having maximum
\(p\)-valuation. The number of remaining terms divisible by
\(p^a\) is at most \(\lfloor(Q_j-1)/p^a\rfloor\), since their
indices are in a single residue class modulo \(p^a\). The
progression step is a unit by (2). Consequently the remaining
small-prime contribution divides \((Q_j-1)!\), and the removed
prime powers contribute at most \(\pi(CD)\log(O_{n,N}(A))\)
to the log product. The elementary Chebyshev estimate
\(\pi(t)=O(t/\log t)\), obtained from the binomial coefficient
bound in the [dense cluster theorem](MULTIVARIATE_DENSE_CLUSTERED_UNIT_GAP_CLOSURE.md#elementary-prime-estimates),
gives
\[
       \log (P_j^{\rm sel})_{>CD}
                    \ge (1-\beta-o(1)) Q_j\log A.              \tag{4}
\]
Every prime in this product belongs to just one marked number
of size \(O_{n,N}(A)\), so there are \(\gg D\) such primes.
There are also \(O_{n,N,\beta}(D)\), by bounding their logarithms
below by \(\log(CD)\) and the log product above by \(O(D\log A)\).

For the full progression of marks \(x\equiv1\bmod j\) in (1),
the minimum outside multiplicity is \(R\) at every mark.
The log of its whole product is \(O_{n,N}(D\log A)\).
The sum of off-mark loads at the selected primes is at least
\(4 M R\log(P_j^{\rm sel})_{>CD}\) by the preceding neutral
blocks. All these primes are above \(D\), as in the necessary
budget. Even for the whole progression, the ratio in the
theorem therefore tends to infinity. This uses literal factors
of the original ramps, including their actual prime valuations.

## Prescribing curvature and retaining a separating prime

Take
\[
 R=\max\left(1,\left\lfloor
                  \frac{(\tau-|\lambda|)A}{2D^2}
                           \right\rfloor\right).
\]
The top pair contributes curvature tending to \(\tau-|\lambda|\)
and signed bias zero. Its \(D^2/A\) tends to zero, including
when only one pair is present. The unit-gap pairs from (3)
have total curvature
\[
                 O\bigl((1+M)RD/A\bigr)=o(1),
\]
since their bases are comparable to \(A\),
\(R=O(A/D^2)\), and \(M=O(D^{1/2})\).
All these pairs have signed gap zero.

Here is a uniform way to choose further ordinary rows without
changing any valuations at the selected primes. Fix an integer
\(J>N\) at least the
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md),
and fix a prime \(p_0>J^2\). In an interval \([H,2H]\) of
possible baselines, restrict to any fixed residue class modulo
\(p_0\). A gap bounded by \(J\) has at most
\(F_J=J\sum_{2\le j<J}j\) forbidden residues modulo a
selected \(p>CD\) that could make any stripped ramp factor
divisible by \(p\) for an order \(j<J\). The counting bound
in the restricted progression is at most
\[
               F_J\left(\frac{H}{p_0p}+O(1)\right)
\]
for each \(p\), because \(p\) is coprime to \(p_0\) for large \(A\).
Writing \(K D\) for an upper bound on the total number of selected
primes from both windows, choose \(C\) also large enough that
\(K F_J/C<1/4\). The upper bound on their number just obtained
allows \(K\) independent of \(C\) once \(C\ge1\). If \(D=o(H)\),
the sum over all forbidden residues uses less than half the
available interval for large \(A\). Thus a baseline in the
prescribed residue class exists for every such choice of \(H\).
One can require avoidance for all gaps up to \(J\) by using
the longest ramp in this count.

When \(\lambda\ne0\), use a baseline \(B_1\) in such an interval
with \(D=o(H_1)\), \(H_1=o(A)\), and \(B_1\equiv0\bmod p_0\). Take
\(\lfloor |\lambda|(B_1+1)\rfloor\) unit-gap rows with sign
\(\operatorname{sgn}\lambda\). They give the missing curvature
and bias. For \(\lambda=0\) omit this block. Correct the resulting
signed gap to \(\Delta\) with unit-gap rows at a safe baseline
\(B_2\equiv0\bmod p_0\), \(H_1=o(H_2)\) and \(H_2=o(A)\).
This costs \(O((H_1+|\Delta|)/H_2)=o(1)\) curvature.
The two blocks are also prime-to-\(p_0\) in all stripped
ramps for \(j<J\).

Finally use gap \(J\) at a safe baseline \(B_3\equiv-1\bmod p_0\)
with
\[
          A D^{-1/2}\log A=o(H_3),\qquad H_3=o(A).
\]
Such \(H_3\) exists because \(D\) is a fixed positive power of \(A\).
Then \(D=o(H_3)\) as well. Use equal numbers \(R_3\) of the
two orientations. At \(p_0\) the only nonunit factor of either
ramp at order \(2\le j<J\) is \(j(B_3+1)\). It occurs once in
the combined endpoints and twice at every interior index, since
each gap is at least \(j\). Thus each combined pair adds at least
one to every interior valuation above the endpoint chord.
The stripped arrays of the fixed core, top pairs and rows from
(3) have at most \(O((1+M)RD+1)\) rows. The Legendre floor
formula bounds their total possible downward chord error by
\(O_J(((1+M)RD+1)\log A)\). Indeed a factorial layer's affine
part cancels from the chord and its fractional-part error
costs at most one per row. There are \(O_J(\log A)\) layers,
and the binomial coefficients are \(p_0\)-units. The prior
unit-gap fillers contribute no \(p_0\)-loss. Therefore an
integer \(R_3=O_J(AD^{-1/2}\log A)\) suffices to put every
interior valuation strictly above the chord simultaneously.
Its curvature cost is \(O(R_3/H_3)=o(1)\).

The [opposed-row angular argument](MULTIVARIATE_UNIVERSAL_OPPOSED_PADIC_ANGLE_MARKING.md#binomial-reductions-separate-root-quotients)
now gives all-degree separation. Rescale the argument over the
algebraic closure of the \(p_0\)-adic field to equalize the two
endpoint valuations, and divide by the constant coefficient.
The residual polynomials are \(1+z^j\).
They are squarefree because \(p_0>J\), and quotients of distinct
roots at adjacent orders reduce to disjoint nontrivial
root-of-unity sets. Equal nonreal arguments would give the
same conjugate quotient. Real shared arguments are excluded by
the strictly positive even moment. The first pair follows by
the same positivity and orders at or above \(J\) by the
absolute angular cutoff. All the witness and filler rows
are below the envelope, have vanishing relative gaps, and
avoid the selected moving primes by construction. The row
counts and the displayed curvature and bias estimates complete
the assertion.
