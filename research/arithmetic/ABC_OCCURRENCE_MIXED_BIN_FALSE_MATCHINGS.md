# False induced matchings for joint real and modular tests

All three radical caps, comparable numerical heights and a least-repair
cutoff allow large induced matchings for a relaxed joint test. At the
critical combined scale of the two tests, the matching retains normalized
exponent \(1/3\). Every complete representative fails the exact sum.

[The joint real and modular occurrence test (IM3)](ABC_OCCURRENCE_INDUCED_ADDITIVE_MATCHINGS.md#a-short-real-filter-cuts-the-required-modulus-exponent-in-half)
· [False candidates for support-prime and auxiliary unit tests (AU1)](ABC_OCCURRENCE_AUXILIARY_UNIT_TESTS.md#au1-uniform-supply-of-false-candidates)

## An induced matching from three exact-power dictionaries

**MB1.** Fix \(e=3^k\) with \(k\ge1\), put \(L=3/e\), and take
\(n=X^L\). Fix \(c,\beta>0\) with \(c+\beta\le1/3\).
Write \(H=X^{1/e}\) and, for large \(X\), set
\[
 E(H)=\exp\!\left(20\frac{\log H}{\log\log H}\right),
 \qquad D=\left\lfloor Xn^{-c}E(H)\right\rfloor .
\]
Let \(Q\asymp n^\beta\) be any odd modulus with
\(3\nmid\varphi(Q)\). Such a modulus exists at every scale by taking
a power of \(5\). It may also be chosen as a prime congruent to \(2\)
modulo \(3\) by the prime number theorem in arithmetic progressions.

For all sufficiently large \(X\) at these fixed parameters there are
at least \(n^{c+\beta-o_X(1)}\) indexed rows \((A_i,B_i,C_i)\) with the
following properties.

- Every letter is an exact \(e\)-th power with root at most \(H\).
  All letters are comparable to \(X\), and \(X/2<C_i\le X\).
  Each radical is at most \(n^{1/3}\). Every mixed triple
  \((A_i,B_j,C_\ell)\) is pairwise coprime.
- With \(b(x)=\lfloor x/D\rfloor\), the joint test for a mixed triple
  is induced with carry zero
  \[
  b(A_i)+b(B_j)=b(C_\ell),\qquad
  A_i+B_j\equiv C_\ell\pmod Q
  \quad\Longleftrightarrow\quad i=j=\ell .
  \]
  The complete labels are distinct in each role.
- Every diagonal row has
  \(0<|A_i+B_i-C_i|<2D\). For every integer \(M\ge e\), its
  actual least \(M\)-powerful repair \(H_M(C_i)\) satisfies
  \(H_M(C_i)\le X^{M/e}\).

In particular \(D=Xn^{-c+o_X(1)}\), with the subpower factor providing
the room for uniform coprime choices at \(c+\beta=1/3\).
The budget \(L\) can be taken arbitrarily small through the indicated
sequence. Each height limit is taken at fixed \(L,c,\beta\).

### Address sets and complete representatives

Behrend's construction gives a set \(S\) of integers in
\([X/(3D),2X/(5D)]\) with no nonconstant three-term arithmetic
progression and with
\[
 |S|=(X/D)^{1-o_X(1)}.
\]
Let \(R_0\subseteq[1,\lfloor Q/20\rfloor]\) be a second such set of
size \(Q^{1-o_X(1)}\). On averaging the translates \(t+R_0\) over
\(t\bmod Q\), some translate contains at least
\(|R_0|\varphi(Q)/Q\) integers that are units modulo \(Q\).
Retain those integers as \(R\), using lifts in the original translated
interval. Thus \(|R|=Q^{1-o_X(1)}\), since
\(\varphi(Q)/Q\ge2^{-\omega(Q)}=Q^{-o(1)}\).
For \(r_1,r_2,r_3\in R\) the congruence
\(r_1+r_2=2r_3\bmod Q\) is an exact integer equality, since the
interval length is less than \(Q/20\). It forces \(r_1=r_2=r_3\).
[Behrend's progression-free construction](https://www.pnas.org/doi/10.1073/pnas.32.12.331)

For each \((s,r)\in S\times R\) prescribe the addresses
\[
 (b(A),A\bmod Q)=(s,r),\qquad
 (b(B),B\bmod Q)=(s,r),\qquad
 (b(C),C\bmod Q)=(2s,2r).
\]
Then \(b(A)\) and \(b(B)\) put their complete values inside
\([X/3,2X/5+D]\), while \(b(C)\) puts the complete target inside
\([2X/3,4X/5+D]\). These intervals lie below \(X\) for large \(X\).

For any bin \(j\) used here, its root interval
\([(jD)^{1/e},((j+1)D)^{1/e})\) has length
\(\asymp_e H D/X\) by the mean-value theorem. Raising to the \(e\)-th
power is a bijection on the unit group modulo \(Q\), since
\(\gcd(e,\varphi(Q))=1\). Each prescribed unit residue therefore
gives a single residue class for the root. The number of consecutive
terms of that class available in the root interval is
\[
 \gg_e HD/(XQ)\gg_{e,\beta} n^{1/3-c-\beta}E(H).
\]
Regard the implicit comparability constants for \(Q\) as fixed.

To choose the three roots \(x,y,z\) at an address, take any available
\(x\), then choose \(y\) coprime to \(x\), and then \(z\) coprime to
\(xy\). In a block of \(J\) consecutive terms of a progression of
step \(Q\), the number coprime to an integer \(d\) with
\(\gcd(d,Q)=1\) is
\[
 J\frac{\varphi(d)}d+O(2^{\omega(d)}).
\]
Indeed, for each squarefree divisor of \(d\), the terms divisible by
it form a single index class modulo that divisor. Inclusion-exclusion
gives the formula with error at most \(2^{\omega(d)}\). For
\(d\le H^2\) one has \(\omega(d)\le5\log H/\log\log H\) for large
\(H\), by comparing \(d\) with the product of the first
\(\omega(d)\) primes. Using \(\varphi(d)/d\ge2^{-\omega(d)}\),
the available block length is greater than \(4^{\omega(d)}\), even at
\(c+\beta=1/3\), by the choice of \(E(H)\). This gives both successive
coprime choices. All chosen roots are units modulo \(Q\).
Set \(A=x^e\), \(B=y^e\), and \(C=z^e\).

Each role address determines an index uniquely. For the target residues
this uses the invertibility of \(2\bmod Q\).
A mixed equality of bin labels forces the three \(s\)-coordinates
equal by construction of \(S\). A mixed congruence forces their
three \(r\)-coordinates equal by construction of \(R\).
Consequently the joint test is induced.

### Global role separation and exact repair

Assign all primes independently one of three role colours. A chosen
row survives if every prime of each letter has that letter's colour.
Its three complete representatives are pairwise coprime, so the
survival probability is \(3^{-\omega(xyz)}=n^{-o_X(1)}\) uniformly,
since \(xyz\le H^3\). Some colouring retains
\[
 |S||R|n^{-o_X(1)}=n^{c+\beta-o_X(1)}
\]
rows. Every mixed triple of their representatives is pairwise coprime,
and the mixed induced test survives restriction. If exact radical
labels must also be distinct within each role, one can thin further
by an inner subpower using the uniform fixed-support count.
[Uniform fixed-support counting](ABC_MINIMAL_REPAIR_SUPPORT_OCCURRENCE.md#exact-reduction-to-occurring-support-unions)

The bin equality implies \(|A_i+B_i-C_i|<2D\). Exact equality would
give a positive integer solution to the Fermat equation of exponent
\(3\), upon replacing the three roots by their \(e/3\)-th powers.
This is impossible. Finally, for \(C=z^e\) and \(M\ge e\), every
prime \(p\mid z\) satisfies \(\max(M,e v_p(z))\le Mv_p(z)\).
Hence
\[
 H_M(C)=\prod_{p\mid z}p^{\max(M,e v_p(z))}
                 \le z^M\le X^{M/e}.
\]
[The Fermat obstruction](https://annals.math.princeton.edu/1995/141-3/p01)

At \(c+\beta=1/3\), the original product radical
also satisfies the exact budget
\(\operatorname{rad}(ABC)\le xyz<z^3=C^L\).
The support counting bound shows that only a negligible
fraction of rows have any given role radical smaller than \(n^{1/3-\eta}\), for
any fixed \(\eta>0\). Indeed the distinct roots
in that role with the smaller radical have count
\(O(n^{1/3-\eta}X^\varepsilon)\). Use arbitrarily
small fixed \(\eta\) at each budget. This also
gives the lower target repair height
\(\operatorname{rad}(z)^M=H^{M-o_X(1)}\) in exponent
sense, since \(\operatorname{rad}(z^e)=\operatorname{rad}(z)\).

## Scope of the construction

The mixed group can thus support a full-exponent induced matching at
\(c+\beta=1/3\), even with three complete radical caps, global role
separation, and controlled exact target repair. The count uses
prescribed false representatives. The summand and target dictionaries
have all valuations divisible by the same \(e\), so the construction
lies on a valuation grid. The real width is
\(Xn^{-c+o_X(1)}\). In (IM3) the width is \(X^\theta\) for a fixed
small \(\theta\), corresponding to \(c=(1-\theta)/L\), which is much
larger than the range of (MB1) as \(L\) decreases. A bound at that
finer real scale, and a bound using the actual arithmetic
equalities, remain open.
