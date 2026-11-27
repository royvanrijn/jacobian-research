# A fixed quadratic field for fifth-twisted norm joins

The [moving-kernel descent](ABC_FIFTH_REFLECTION_MOVING_KERNELS.md)
requires a moment of \(5\)-torsion in \(\mathbf Q(\sqrt d)\) and
\(\mathbf Q(\sqrt{5d})\) varying together. To compare the witness
sets of the two fields, count primitive norm records that share a
prime, or that use the same pair of primes. Eliminating \(d\) in such
a record gives a factorization in the fixed field \(\mathbf Q(\sqrt5)\).
The prime or prime pair shared by both records can be removed from the
norm in that field. It survives as a generator in a two-dimensional
lattice with a modulus imposed by the two normalized coefficients.

## Uniform count of shared-prime and shared-edge records

**Proposition FJ1.** Fix an integer \(\ell\ge1\) and \(L>0\), and let
\(Z\ge2,\ V_1,V_2\ge1\). Count tuples \((d,p,q,r,u,v,s,w)\) with
\(d\ne1\) nonzero signed squarefree and prime to 5, with odd primes
\(p,q,r\in[Z,2Z)\) prime to \(10d\), and with \(p\ne q,r\).
Assume
\[
 \begin{split}
 u^2-dv^2 &=4(pq)^\ell,\qquad
 s^2-5dw^2=4(pr)^\ell,\\
 0<|v|&\le V_1,\qquad 0<|w|\le V_2,\qquad
 |u|,|s|\le LZ^\ell,\qquad (v,pq)=(w,pr)=1 .
 \end{split}                                                     \tag{1}
\]
For \(q\ne r\) their number is
\[
 \ll_{\ell,L,\epsilon} (ZV_1V_2)^\epsilon
                               \big(Z^2V_1V_2+Z^3\big),         \tag{2}
\]
and for \(q=r\) their number is
\[
 \ll_{\ell,L,\epsilon} (ZV_1V_2)^\epsilon
                               \big(V_1V_2+Z^2\big).            \tag{3}
\]
In both counts \(d\) varies over its full allowed range. In (2), for
fixed \(q,r\) and fixed normalized coefficient pair \((c,t)\) defined
below, there are \(O_{\ell,L,\epsilon}((ZV_1V_2)^\epsilon)\)
records when \(|ct|\gg Z\). In (3) both primes can vary with the same
per-coefficient bound when \(|ct|\gg Z^2\).

The primitive norm equations are the starting pair data in
[Heath-Brown–Pierce, Proposition 2.1 and §3](sources/heath-brown-pierce-1409.3177v2.pdf)
and [Chan–Koymans, §5 and Proposition 5.1](https://arxiv.org/pdf/2505.00611).
Their single-field counts use the trace in a square-root congruence
modulo the square of the corresponding coefficient. The join uses
two traces in the fixed-field factorization below.

### Removing the shared norm

Set \(\delta=\gcd(v,w)>0\), \(c=v/\delta\), \(t=w/\delta\), and
\(M=\mathbf Q(\sqrt5)\). In its ring of integers let
\[
 \beta=-s c+u t\sqrt5,\qquad
 N_{M/\mathbf Q}(\beta)=4p^\ell(r^\ell c^2-5q^\ell t^2).          \tag{4}
\]
The right side is nonzero by the parity of the 5-valuations of its
two terms. At \(p\), the coefficients \(u,s,c,t\) are units and the two
norm equations imply that \(d\) and \(5d\) are squares modulo \(p\).
Thus \(p\) splits in \(M\), and exactly one of the primes above \(p\)
divides \(\beta\): the norm is divisible by \(p\), and the two residues
\(-sc\pm ut\sqrt5\) cannot both vanish. All of the \(p\)-valuation of
(4) goes into that factor, at least \(\ell\). If \(q=r\), the same
argument applies at \(q\).

The class number of \(M\) is one, by the Minkowski bound \(\sqrt5/2\).
Let \(\nu\) generate the marked prime ideal above \(p\) for \(q\ne r\),
and their product above \(p,q\) for \(q=r\). Write \(J=p\) or \(pq\),
respectively. Multiplying by a unit makes the absolute values of the
two embeddings of \(\nu\) each \(\asymp \sqrt J\), with constant
depending only on \(M\). The vector of \(\nu\) in the integral basis
\(1,\omega\), \(\omega=(1+\sqrt5)/2\), is primitive since
\(|N(\nu)|=J\) is squarefree. Define \(\gamma=\beta/\nu^\ell\).
The possibilities for its norm are
\[
 N(\gamma)=
 \begin{cases}
 \pm4(r^\ell c^2-5q^\ell t^2), & q\ne r,\\
 \pm4(c^2-5t^2), & q=r .
 \end{cases}                                                     \tag{5}
\]
In (2), fix \(q,r,c,t\) before choosing \(\gamma\). In (3), fix only
\(c,t\). Both absolute embeddings of \(\gamma\) are
\(O_{\ell,L}(Z^{\ell/2}\max(V_1,V_2))\) in (2) and
\(O_{\ell,L}(\max(V_1,V_2))\) in (3). For each nonzero norm in (5)
there are divisor-many integral ideals of that absolute norm. Each
ideal has \(O_{\ell,L}(\log(2ZV_1V_2))\) generators meeting the embedding
bound, by the rank-one unit lattice and the nonzero integral norm.
This yields \(O_{\ell,L,\epsilon}((ZV_1V_2)^\epsilon)\) possibilities
for \(\gamma\), uniformly in the chosen coefficients and primes.

### Two coefficient moduli determine lattices for the generator

Write a bar for conjugation in \(M\). Since
\(\beta+\bar\beta=-2sc\) and \(\beta-\bar\beta=2ut\sqrt5\),
one has locally
\[
  \gamma(\nu/\bar\nu)^\ell\equiv-\bar\gamma\bmod c,\qquad
  \gamma(\nu/\bar\nu)^\ell\equiv\bar\gamma\bmod t.                 \tag{6}
\]
The ratio \(\nu/\bar\nu\) is a local unit at the moduli, since
\((ct,J)=1\). For an odd prime \(j\ne5\) dividing \(c\), the norm
in (5) is prime to \(j\) because \(j\nmid q,t\). At such a prime
dividing \(t\), it is prime to \(j\) because \(j\nmid r,c\).
Thus \(\gamma\) is invertible there. If \(2\mid ct\), exactly
one coefficient is even, so the norm in (5) has 2-valuation 2.
As 2 is inert in \(M\), division by \(\gamma\) loses at most one
power of 2 from the modulus. At \(5\mid t\) the norm is a
5-unit. At \(5\mid c\) its 5-valuation is one, so the valuation
of \(\gamma\) at the prime of \(M\) over 5 is one. Division loses
at most one rational power of 5 from the modulus. In particular
one can reduce the product of coefficient moduli by a factor at most 10
to obtain a prescribed \(\ell\)-th power of \(\nu/\bar\nu\)
modulo every remaining prime power.

Over \(\mathcal O_M/j^e\mathcal O_M\) there are
\(O_\ell(1)\) unit solutions of such a power congruence, uniformly
in \(j,e\). At primes away from \(2\cdot5\ell\), use the bound
for roots in the residue fields and unique lifting. At the remaining
finitely many primes the local unit group modulo \(1+j^2\mathcal O_{M,j}\)
is of fixed size, and on \(1+j^2\mathcal O_{M,j}\) the logarithm
turns the \(\ell\)-th power into multiplication by \(\ell\).
Its kernel modulo \(j^e\) has size at most
\(j^{2v_j(\ell)}\). Consequently the total number of choices
of root ratios across the moduli is \(O_{\ell,\epsilon}(|ct|^\epsilon)\).

For each fixed choice of ratios and one generator
\(\nu_0=a_0+b_0\omega\) with those ratios, every other
\(\nu=a+b\omega\) satisfies \(ab_0-a_0b\equiv0\) modulo the
remaining modulus. Indeed
\(\nu\bar\nu_0-\bar\nu\nu_0\) is
\(\pm\sqrt5(ab_0-a_0b)\). Divisibility by \(j^e\mathcal O_M\) of
this expression forces divisibility of the rational integer by \(j^e\),
also at the ramified prime 5. The integer pair \((a_0,b_0)\) is
primitive at every prime, so this is a lattice of index equal to
the remaining modulus, at least \(|ct|/10\).
A ratio with no generator contributes zero.

The vectors have size \(O(\sqrt Z)\) in (2) and \(O(Z)\) in (3).
For primitive lattice vectors in a ball of radius \(O(R)\) in rank
two, the usual successive-minima bound gives
\(O(1+R^2/[\mathbf Z^2:\Lambda])\). If the vectors in the ball
span rank one, a rational line carries at most two primitive vectors.
If they span rank two, the second minimum is \(O(R)\), and the product
of the minima is comparable to the covolume. This yields a bound of
\(O(1+Z/|ct|)\) or \(O(1+Z^2/|ct|)\) in the two cases, for each
choice of \(\gamma\) and the root ratios. In particular, two generators
with the same ratios cannot be linearly independent once the
modulus exceeds the bound on their determinant.

### Recovery and summation

The norm of \(\nu\) recovers \(p\), or recovers \(p,q\) up to order
when \(q=r\). From \(\nu^\ell\gamma=-sc+ut\sqrt5\) recover \(s,u\)
using nonzero \(c,t\). Then
\[
                   d\delta^2=\frac{u^2-4(pq)^\ell}{c^2}
\]
recovers both the signed squarefree radicand \(d\) and the positive
integer \(\delta\). Thus the removal of the common coefficient does
not leave an uncharged multiplicity over fields. Sum the bound
\(1+R^2/|ct|\) over \(1\le |c|\le V_1,\ 1\le |t|\le V_2\),
absorbing harmonic sums in the subpower. There are \(O(Z^2)\)
choices of separate endpoints in (2) and at most two orderings of
the recovered prime pair in (3). This proves both estimates.

## Effect at the fifth-power moment scale

At the diagonal \(5\)-torsion crossover \(H=X^{9/26}\), the
[coherent pair construction](ABC_COHERENT_PRIME_PATH_HIGHER_ODD_TORSION.md)
uses \(Z=X^{2/13+O(\eta)}\) and \(V_i\ll X^{7/26+O(\eta)}\).
By (3) the number of signed squarefree radicands with even one
such shared-edge record is \(O_\epsilon(X^{7/13+\epsilon})\),
with the auxiliary \(\eta\) sufficiently small. The marginal tail
there is \(O_\epsilon(X^{17/26+\epsilon})\). The shared-edge
subfamily hence has a power saving in a joint count, independent of
whether either field has large torsion. This applies also when the
first radicand is divisible by 5, by interchanging the two fields
and using the prime-to-5 member, which changes the field scale by
at most a fixed factor.

High torsion in both fields gives a lower bound on edges in each
prime graph separately. There is no shared-edge lower bound in that
argument. Even two equivalence relations with bounded block size
can each supply order as many edges as vertices and have empty
edge intersection, as with two distinct pair matchings on each
block of four vertices. Estimate (2) allows the two endpoints to
differ, but at the same scale its \(Z^2V_1V_2\) term has exponent
\(22/26\), before any witness multiplicity lower bound. A joint
argument that saves on those occupied generator lattices or forces
a sufficient shared-edge subfamily would be needed to improve the
twisted moment in (4) of the moving-kernel note.

## Exact controls

At \(\ell=5,d=-2,p=11,r=19\) use \((s,w)=(807014,307220)\).
The table gives an instance of each case with \(\beta=\nu^5\gamma\).

| \(q\) | \((u,v)\) | \((\delta,c,t)\) | \(\nu\) | \(\gamma\) |
| --- | --- | --- | --- | --- |
| \(3\) | \((10550,4756)\) | \((4,1189,76805)\) | \(4+\sqrt5\) | \(-81099754+36309330\sqrt5\) |
| \(19\) | \((1092486,448100)\) | \((20,22405,15361)\) | \(17+4\sqrt5\) | \(-818710+366878\sqrt5\) |

These test (1), normalization and both norms in (5) by integer
arithmetic in the basis \(1,\sqrt5\). Both rows have a factor 5 in a
normalized coefficient, so (6) in these examples includes the
ramified prime of the fixed field.
