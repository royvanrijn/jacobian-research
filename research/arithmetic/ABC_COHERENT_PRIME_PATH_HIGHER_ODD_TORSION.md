# Coherent prime paths for higher odd torsion

The full-modulus prime-path count also applies to powers \(\ell\ge5\).
At the fifth power the rank-two normal produces a plane cubic in the
three moduli. A uniform count on its smooth fibers, with a separate
singular-fiber estimate, improves the signed quadratic \(5\)-torsion
moments. For every prime at least seven, a divisor count followed by a
prime-power congruence bypasses the lattice-rank decomposition. Real
witnesses use generators in one bounded interval of the unit logarithm.

## Signed tails and moments

For a sign \(\iota\in\{-1,1\}\), positive squarefree \(d\) with
\(\iota d\ne1\), and an odd integer \(\ell\ge5\), set
\[
 h_\ell(\iota d)=|\operatorname{Cl}(\mathbb Q(\sqrt{\iota d}))[\ell]|,
 \qquad
 N_{\ell,\iota}(H;X)=\#\{X\le d<2X:h_\ell(\iota d)>H\}.
\]
In the count \(d\) is always restricted as above. Write
\(\rho_5=7/4\) and \(\rho_\ell=2\) for \(\ell\ge7\).

**Theorem (higher odd prime-path tail).** For \(X\ge2,H\ge1\) and every
\(\varepsilon>0\), with constants allowed to depend on \(\ell,\varepsilon\),
\[
 \begin{split}
 N_{\ell,\iota}(H;X)\ll_{\ell,\varepsilon}X^\varepsilon
 \Big(&1+XH^{-2}+\min\big\{
       XH^{-1}+X^{\ell/2}H^{-\ell-1},\\
     & \mathbf 1_{\{\ell=5\}}X^{(3\ell-4)/2}H^{-(3\ell-1)}
       +X^{(1+(\ell-1)\rho_\ell)/2}H^{-(\ell\rho_\ell+1)}
                       \big\}\Big).                              \tag{1}
 \end{split}
\]
For every fixed \(1\le k\le\ell+1\), the sum of \(h_\ell(\iota d)^k\) for
positive squarefree \(d\le X\) satisfies
\[
 \ll_{\ell,k,\varepsilon}X^{\varepsilon}
 \left(X^{1+k a_\ell}+X^{1+(k-1)b_\ell}
             +X^{1+(k-2)c_\ell}+X^{k/2}\right),                    \tag{2}
\]
where
\[
   a_\ell=\frac{\ell-2}{2(\ell+1)},\qquad
   b_5=\frac9{26},\qquad b_\ell=\frac{2\ell-3}{4\ell}\quad(\ell\ge7),\qquad
   c_\ell=\frac{\ell-1}{2\ell}.
\]
In particular for \(\ell=5\) the exponents at orders \(3,4,5\) are
\(7/4,53/26,5/2\), for either sign. The corresponding exponents in the
imaginary-field bounds of [Heath-Brown–Pierce, Theorem 1.5](sources/heath-brown-pierce-1409.3177v2.pdf)
are \(9/5,11/5,13/5\). More generally (2) improves the higher-order branch
of that theorem for prime \(\ell\) and \((\ell^2-1)/(2\ell-1)<k<\ell+1\), since \(b_\ell<c_\ell\). For \(\ell\ge7\), the moment bound reaches \(X^{k/2+\varepsilon}\) for \(k\ge1+2\ell/3\), extending past \(\ell+1\) by the class-number bound. For example the signed seventh-torsion moment at order five has exponent \(18/7\), and the signed thirteenth-torsion moment at order ten has exponent \(5\).
For comparison, [Koymans–Thorner, Theorem 1.1](https://arxiv.org/pdf/2309.07204)
gives exponent \(1+k(1/2-1/(\ell+1))\) for \(1\le k\le\ell+1\) for
quadratic fields. These comparisons concern the cited bounds.

## Coherent witnesses with a bounded unit logarithm

Fix \(K=\mathbb Q(\sqrt{\iota d})\) and choose one ideal above each odd
unramified split prime in \([Z,2Z)\). Use the ordinary class group in
the imaginary case and the narrow class group \(C^+\) in the real case.
The kernel from narrow to ordinary classes is a \(2\)-group, so the map
identifies their \(\ell\)-torsion. Partition the chosen primes into cosets
of \(\ell\)-torsion. In one coset, choose an integral ideal \(J\) whose class
is inverse to the common class of the \(\mathfrak p_i^\ell\), and choose
generators \(\gamma_i\) of \(J\mathfrak p_i^\ell\). In the real case take
them totally positive. Write \(R^+\) for the logarithm of the first embedding
of a generator greater than \(1\) of the totally positive units. Normalize
\[
                  a_i=\tfrac12\log(\gamma_i/\bar\gamma_i)
\]
to \([0,R^+)\) by multiplying \(\gamma_i\) by a totally positive unit,
then subdivide that interval into intervals of length at most one. In
the imaginary case keep each class coset in one part. For each part set
\[
                    \alpha_{ij}=\bar\gamma_i\gamma_j/NJ .
\]
Its principal ideal is \((\overline{\mathfrak p}_i\mathfrak p_j)^\ell\).
Thus \(\alpha_{ij}\) is integral, \(\alpha_{ji}=\bar\alpha_{ij}\), and
\(\alpha_{ij}\alpha_{jr}=p_j^\ell\alpha_{ir}\). In the real case its
conjugate logarithms are \(\frac\ell2\log(p_i p_j)\pm(a_j-a_i)\), giving
\[
 \alpha_{ij}=(u_{ij}+v_{ij}\sqrt{\iota d})/2,\qquad
 |u_{ij}|\ll_\ell Z^\ell,\qquad |v_{ij}|\ll_\ell Z^\ell/\sqrt d .   \tag{3}
\]
In the imaginary case absolute value gives these bounds directly.
For \(i\ne j\), the ideal is not invariant under conjugation, so
both \(u_{ij}\) and \(v_{ij}\) are nonzero. Also \((p_i p_j,v_{ij})=1\)
and \(\alpha_{ij}\) is primitive in the integral basis, by its ideal
factorization. This chooses one common set of generators for all edges
of a part. The unit normalization adapts the interval construction of
[Chan–Koymans, proof of Theorem 2.1](https://arxiv.org/pdf/2505.00611)
from pairs to coherent triples.

There are \(O(\sqrt d\log d/h_\ell)\) parts in either sign, by the
class number formula and the bound \(hR\ll\sqrt d\log d\) in the real
case. Here \(h^+\le2h\), \(R\le R^+\le2R\), and \(R\) is bounded away
from zero. For \(Z\ge X^{1/(2\ell)}/4\) and \(Z\le X\),
[Heath-Brown–Pierce, equations (2.5)–(2.7)](sources/heath-brown-pierce-1409.3177v2.pdf)
give \(\gg Z/\log Z\) split primes outside \(O_{\ell,\eta}(X^\eta)\)
exceptional \(d\), for any fixed \(\eta>0\). The same proof counts
either signed fundamental discriminant up to \(8X\): the primitive
quadratic character has conductor \(|D|\le8X\), and the character large
sieve uses \(Z^{4\ell}\gg_\ell X^2\). Removing ramified or even primes
costs \(O(\log X)\).

Apply this at \(Z_0=X^{1/(2\ell)}/4\). Two primes in one part would
give \(0<|v_{ij}|\le 2\sinh(1)(2Z_0)^\ell/\sqrt d<1\) in the real
case, and the same inequality with \(2\) in place of \(2\sinh(1)\)
in the imaginary case. Hence each part has at most one prime and
\(h_\ell\ll_\ell X^{c_\ell}(\log X)^2\) outside that exceptional set.

Now take \(H\le X^{c_\ell+\eta}\) and \(Z=X^{1/2+4\eta}/H\), with
\(0<\eta\) sufficiently small in terms of \(\ell\).
For a nonexceptional \(d\) with \(h_\ell>H\), discarding parts of size
less than four costs at most \(O(\sqrt X\log X/H)=o(Z/\log Z)\).
The sum of the cubes of the remaining sizes is
\(\gg_\eta Z^3 H^2 X^{-1-\eta}\), by Hölder and the number of parts.
Each part of size \(s\ge4\) gives \(s(s-1)(s-2)\gg s^3\) ordered
triangles. Thus \(d\) gives
\[
                         \gg_\eta Z^3H^2X^{-1-\eta}               \tag{4}
\]
coherent ordered triangles satisfying (3). Distinct triples of
vertices give distinct ordered rational prime triples in this field.
The analogous square sum gives \(\gg_\eta Z^2 H X^{-1/2-\eta}\) pairs.
Their global count uses the same congruence as in
[Heath-Brown–Pierce, Proposition 2.3 and §3](sources/heath-brown-pierce-1409.3177v2.pdf)
and [Chan–Koymans, Proposition 5.1](https://arxiv.org/pdf/2505.00611).
For \(w=p_ip_j\), the norm equation puts \(u\) in \(v^2\)-residue
classes satisfying \(u^2\equiv4w^\ell\bmod v^2\); there are
\(O_\varepsilon(V^\varepsilon)\) classes by the unit condition
\((w,v)=1\). At an odd prime power there are at most two roots, and
at \(2\) dividing \(u\) by \(2\) reduces to a unit square-root
congruence with at most a bounded number of lifts. For \(u\ge0\),
\[
       |u-2w^{\ell/2}|=\frac{d v^2}{u+2w^{\ell/2}}
                             \ll_\ell Xv^2/Z^\ell
\]
in either sign. Thus each \(w,v\) yields
\(O_\varepsilon(V^\varepsilon(1+X/Z^\ell))\) traces, also allowing
negative \(u\). There are \(O(Z^2V)\) choices for \(w,v\), where
\(V=Z^\ell/\sqrt X\). Each record determines the field and the two
primes, with at most two orderings. Hence the global pair count is
\(O_{\ell,\varepsilon}(X^\varepsilon(Z^2\sqrt X+Z^{\ell+2}/\sqrt X))\).
The root counts and the unit logarithm intervals work for composite odd \(\ell\)
as well. Dividing by the pair lower bound gives
\(N_{\ell,\iota}\ll X^\varepsilon(1+XH^{-1}+X^{\ell/2}H^{-\ell-1})\).
The nonexceptional cap handles all larger \(H\), with error
\(O_{\ell,\eta}(X^\eta)\).

## Exact full moduli and the lattice count

For a triangle write \(P_i=p_i^\ell\) and
\((v_{12},v_{23},v_{13})=\delta(x,y,z)\), where \(\delta\) is the
positive common gcd. Put \(V=\max(1,Z^\ell/\sqrt X)\). The three signed
normalized integers are nonzero of size \(O_\ell(V)\) and are pairwise
coprime. For the last assertion, at each prime reanchor at a vertex
whose norm prime is different from it. In the integral basis the
edge coefficients, after removing the common integral-basis factor,
are determinants of primitive vertex vectors divided by the
anchor prime power. Divisibility at a prime power means equality of
the two projective directions of primitive vectors, an equivalence
relation. So the minimum of the three coefficient valuations is
attained at least twice. This gives pairwise coprimality after
division by \(\delta\), including at \(2\).

The pencil identity, obtained by taking the norm of the linear
relation on the anchor vectors, has the same expressions for both
signs. Put \(A=P_1y^2,B=P_2z^2,C=P_3x^2\). Specifically,
\[
 u_{12}=(A+B-C)/(yz),\qquad
 u_{23}=(B+C-A)/(xz),\qquad u_{13}=(B-A-C)/(xy).                   \tag{5}
\]
For example, the norm of \(v_{12}\alpha_{13}
=v_{13}\alpha_{12}-v_{23}P_1\), after dividing by \(P_1\delta^2\),
is \(C=B-u_{12}yz+A\). Reanchoring gives the other identities.
The quantities \(A,B,C\) are respectively divisible by \(y^2,z^2,x^2\).
The trace integrality therefore gives
\[
               x\mid A-B,\qquad y\mid C-B,\qquad z\mid A-C.         \tag{6}
\]
These are the entire moduli, with no removal of prime powers.
For an actual record, (5) determines \((4P_1P_2-u_{12}^2)/x^2
=-\iota d\delta^2\), so the data \((p_1,p_2,p_3,x,y,z)\) determine
the sign, squarefree radicand and scale uniquely.

## Counting by a coefficient divisor and the remaining prime

**Proposition (Rank-independent triangle count).** For a fixed odd prime \(\ell\) and either radicand sign, the number of coherent triangles with primes in \([Z,2Z)\) and \(|x|,|y|,|z|\ll_\ell V\) is
\[
                         \ll_{\ell,\varepsilon}
                                 (ZV)^\varepsilon(Z^2V^2+Z^3)      \tag{7a}
\]
when \(Z\ge2,V\ge1\). The count is uniform over all fields and retains the exact field recovery in (5).

Fix distinct \(p_1,p_2\) and nonzero signed \(y,z\). Then \(A-B=P_1y^2-P_2z^2\ne0\): the valuations of its two terms at \(p_1\) have different parity because \(\ell\) is odd and \(p_1\ne p_2\). Hence \(x\mid A-B\) leaves divisor-many signed \(x\), with subpower loss in \(ZV\). For each one, the full moduli \(y,z\) combine to give
\[
             P_3 x^2\equiv P_2z^2+P_1y^2\pmod{|yz|}.               \tag{7b}
\]
The right side and \(x^2\) are units modulo \(|yz|\). Indeed \(p_2\) is prime to \(y\), \(p_1\) is prime to \(z\), and the normalized coefficients are pairwise coprime. Thus \(p_3\) is in a subpower number of unit classes modulo \(|yz|\): the number of \(\ell\)-th roots is at most \(\ell\) at each odd prime power and at most one at powers of \(2\). In the interval for \(p_3\) this contributes a subpower multiple of \(1+Z/|yz|\). Summing over \(y,z\) and over the first two primes gives (7a). For each resulting tuple the recovery in (5) recovers the traces, sign, squarefree radicand and scale. A fixed ordered prime tuple and field can therefore be counted in the selected witnesses by this bound.

This uses the same full moduli as the lattice method. It counts the coefficient on edge \(12\) by divisibility before choosing the third prime via the other two moduli. In particular the whole \(V^3\) rank-one term can be dropped when the final estimate allows the cost \(Z^2V^2\).

## The lattice count at the fifth power

Write \((q_1,q_2,q_3)=(|x|,|y|,|z|)\), \(Q=q_1q_2q_3\).
Coprimality of the coefficients and the incident prime units turn (6)
into at most \(O_{\ell,\varepsilon}(V^\varepsilon)\) homogeneous
lattices of index \(Q\) for a fixed modulus triple, each with
\[
              p_1=c_1p_2\bmod q_1,\quad
              p_3=c_2p_2\bmod q_2,\quad
              p_1=c_3p_3\bmod q_3.                                \tag{7}
\]
The coefficients \(c_i\) are unit roots, with
\(y^2c_1^\ell=z^2\bmod q_1\), \(x^2c_2^\ell=z^2\bmod q_2\),
and \(y^2c_3^\ell=x^2\bmod q_3\). Each unit power equation has at most
\(\ell\) roots at an odd prime power, and at most one at powers of \(2\).
A modulus \(1\) imposes no restriction. The index follows by CRT,
since at each modulus one coordinate is fixed by the other two.

Count primitive prime vectors in a ball of radius \(O(Z)\).
Rank one contributes at most one positive primitive vector for each
lattice and hence \(O_{\ell,\varepsilon}(V^{3+\varepsilon})\) records.
Rank three gives \(O(Z^3/Q)\) points per lattice and
\(O_{\ell,\varepsilon}(Z^3 V^\varepsilon)\) in total by summing the three
harmonic sums. In rank two use the full plane intersection lattice and
its integer normal \(n\), the cross product of a basis. Its covolume is
\(P=\|n\|\), and the count per lattice is \(O(Z^2/P)\). These lattice
bounds are the successive-minima estimates for the rank of the span of
vectors in the ball. They follow by using the first \(r\) minima in that
span and its rank-\(r\) lattice covolume, since each of those minima is
at most the ball radius.

The normal has \(q_1\mid n_3,q_2\mid n_1,q_3\mid n_2\), by substituting
(7) into the integer minors of its basis. If one coordinate is zero,
the two moduli whose defining relations involve that coordinate
both divide the other two coordinates, because their coefficients
in (7) are units. Two zero coordinates cannot support positive primes
in the plane. Summing the one-zero case over each dyadic box of
moduli costs \(O_{\ell,\varepsilon}(Z^2 V^{1+\varepsilon})\), since
\(P\) is at least the product of those two coprime moduli. For all
coordinates nonzero, \(P\gg w=\max_i q_i\). The next subsection uses an additional relation among these nonzero coordinates.

### The cubic in the moduli at \(\ell=5\)

Write \((n_3,n_1,n_2)=(a_1q_1,a_2q_2,a_3q_3)\). The three normal
congruences in (7) give
\[
                         Q\mid a_1^5q_1^3+a_2^5q_2^3+a_3^5q_3^3. \tag{8}
\]
For example modulo \(q_1\) the normal satisfies \(n_2=-c_1 n_1\).
Raise to the fifth power, use \(y^2c_1^5=z^2\), substitute the
signed normal, and cancel \(q_2^2q_3^2\), which is a unit.
This gives \(a_2^5 q_2^3+a_3^5 q_3^3\equiv0\bmod q_1\).
The analogous normal relations modulo \(q_2,q_3\) give the other
congruences. The argument includes \(2\).

Fix dyadic sizes \(Q_i\) for \(q_i\), \(Q_*=\prod Q_i\),
\(u=\min Q_i\), and a dyadic size for \(P\). If \(P\gg Q_*/V^{7/4}\),
the unrestricted number \(O(Q_*)\) of modulus triples already gives
\(O(Z^2 V^{7/4+\varepsilon})\) after the lattice weights.
Otherwise use (8) to write
\[
                 a_1^5 q_1^3+a_2^5 q_2^3+a_3^5 q_3^3=mq_1q_2q_3.
                                                                    \tag{9}
\]
There are \(O(P^3/Q_*)\) nonzero signed triples \(a\), with no
primitive restriction on \(a\), and for each triple at most
\(O(1+P^5/(Q_* u^2))\) integers \(m\). The bound for \(m\) follows
from \(|a_i|q_i\ll P\), so each term on the left has size
\(O(P^5/u^2)\).

For fixed \(a,m\) this is a projective plane cubic. A singular point
has no zero coordinate, by the three partial derivative equations.
Multiplying those equations gives \(m^3=27(a_1a_2a_3)^5\).
There is therefore at most one integer \(m\) with a singular cubic for
each \(a\). For this value \(O(V^2)\) modulus triples suffice, by fixing
two coordinates and solving a nonzero cubic in the third. All other
cubics are smooth and geometrically irreducible. By
[Walsh, Theorem 1.1](sources/walsh-1308.0574v2.pdf) each has
\(O(V^{2/3})\) rational points of height \(O(V)\), with an absolute
constant. The positive coprime moduli have only one representative per
rational point. The possible root classes for (7) contribute only
\(O_\varepsilon(V^\varepsilon)\) per modulus triple.

After division by \(P\), the respective cost of the singular and smooth
curves, before the common factor \(Z^2 V^\varepsilon\), is at most
\[
           V^2 P^2/Q_*,\qquad
           V^{2/3}\big(P^2/Q_*+P^7/(Q_*^2u^2)\big).
\]
In the low-\(P\) range, \(P^2/Q_*\ll Q_*/V^{7/2}\ll V^{-1/2}\) and
\(P^7/(Q_*^2u^2)\ll Q_*^5/(V^{49/4}u^2)\ll V^{3/4}\).
The last bound uses \(Q_*^5/u^2\ll V^{13}\).
Thus both costs are \(O(V^{7/4})\). Summing the dyadic boxes gives
\(O_\varepsilon(Z^2V^{7/4+\varepsilon})\) for rank two at \(\ell=5\).
The singular case is counted with its own bound even when it contains
rational lines.

The [refined cofactor count](ABC_FIFTH_COFACTOR_CUBIC_COUNTS.md#the-rank-two-tail-term) identifies the rational lines in the singular case and uses the primitive coefficients in the smooth count. Its [rank-one estimate](ABC_FIFTH_COFACTOR_CUBIC_COUNTS.md#small-cofactor-quotients-in-rank-one) covers singular and small-quotient records at the crossover.

## Summation and comparison

The triangle total is
\(O_{\ell,\varepsilon}(X^\varepsilon(Z^3+\mathbf 1_{\{\ell=5\}}V^3+Z^2V^{\rho_\ell}))\).
Recovery from (5) permits summing this triple count over all signed
squarefree fields in the window. Divide by (4), set
\(Z=X^{1/2+4\eta}/H\) and \(V=Z^\ell/\sqrt X\), absorbing the constant part of \(V\) into \(Z^3\), and take the
auxiliary allowances sufficiently small. This gives the second tail
inside the minimum in (1), together with \(XH^{-2}\) and the exceptional
term. At larger \(H\) the cap again leaves only \(O(X^\eta)\) exceptional
fields. Combining with the pair bound proves (1).

For the moments, discard the small exceptional set at \(Z_0\) using
the class-number bound \(O(X^{1/2+\eta})\). There are only \(O(\log X)\)
dyadic \(H\)-boxes for nonexceptional fields up to
\(H\le X^{c_\ell+\eta}\). Their additional split-prime exceptions likewise
cost \(O_{\ell,\eta}(X^{k/2+O_k(\eta)})\) in total.
Write \(H=X^h\), suppressing the allowed subpowers. For \(h\le a_\ell\)
the trivial tail gives exponent at most \(1+k a_\ell\). Between \(a_\ell\)
and \((\ell-2)/(2\ell)\), the decreasing pair term contributes at most
\(1+k a_\ell\) for \(k\le\ell+1\), since it has decay degree \(\ell+1\)
and equals \(X\) at \(h=a_\ell\). The pair term \(XH^{-1}\) is sufficient
from there up to \(h=b_\ell\), contributing at most \(1+(k-1)b_\ell\).
For \(\ell=5\), the rank-one term crosses \(XH^{-1}\) at \(h=b_5\) and the rank-two term crosses earlier. For \(\ell\ge7\), the triangle term from (7a) crosses \(XH^{-1}\) at \(h=b_\ell\). These terms have decay degree at least \(\ell+1\), so above \(b_\ell\) their moment contributions are bounded by their values there.
The \(XH^{-2}\) term for \(k\ge2\) has its maximum at the cap,
giving \(1+(k-2)c_\ell\); for \(1\le k<2\) its contribution is at most
\(1\). This proves (2) in a dyadic \(d\)-window, and the sum over the
windows costs a subpower.

At the old rank-one crossover \(V=Z^{(\ell+2)/4}\), (7a) saves a factor \(V^{(\ell-6)/(\ell+2)}\) over \(V^3\) for \(\ell\ge7\), up to the permitted subpower factors. For \(\ell=3,5\) the same bound is too large there. The occupied-lattice problem still governs the rank-one terms for those exponents.

## Exact controls

For \(\ell=5\), take generators \(\pi_i\) of primes of norms \(3,11,19\)
in \(\mathbb Q(\sqrt{-2})\) given by \(1+\sqrt{-2},3+\sqrt{-2},1+3\sqrt{-2}\).
Using \(\alpha_{ij}=(\bar\pi_i\pi_j)^5\) gives, in order \(12,23,13\),
coefficients \((-4756,-871856,32596)\) and traces
\((-10550,-273582,-16786)\), with \(\delta=116\).
The same construction in \(\mathbb Q(\sqrt2)\), using
\(11+7\sqrt2,13+7\sqrt2,15+7\sqrt2\) of norms \(23,71,127\), gives
coefficients \((-800654092,-13431133772,-9721697944)\) and traces
\((1152625050,24479090914,13779455974)\), with \(\delta=28\).
These test (5)–(6) and recovery of \(-\iota d\delta^2\), with opposite
signs. The fixtures only test the algebraic identities.

Taking seventh powers of the same prime products gives common scales \(4\) and \(196\) in the imaginary and real fixtures, respectively. Both give unit power congruences in (7b), with the exponent \(7\) on the unknown prime, and recover the same signed fields. In either exponent the congruence follows numerically by substituting the displayed prime generators.

## Prime occupancy with the field window retained

The estimate (7a) charges one possible small prime to every compatible
class when \(|yz|>Z\). At \(\ell=5,V=Z^{7/4}\), its first term is
\(Z^{11/2}=Z^{1/4}V^3\). A saving strictly larger than \(Z^{1/4}\)
in that term would therefore improve the old crossover. The saving
below is an **open estimate**. The exact reduction, removal of small
coefficients, and conditional moment implication are proved here.

### An explicit occupancy weight

Fix a constant \(B\ge128\), absorbing the fixed witness constants.
Let \(a=p_1^5,b=p_2^5\), choose nonzero signed \(y,z\) of absolute
value at most \(BV\), and choose a signed divisor
\[
             x\mid ay^2-bz^2,\qquad 0<|x|\le BV.                 \tag{PO1}
\]
Require pairwise coprimality of \(x,y,z\), and
\((x,p_1p_2)=(y,p_2)=(z,p_1)=1\). Put \(m=|yz|\), and let
\[
 R(p_1,p_2;x,y,z)=
 \{r\in(\mathbb Z/m\mathbb Z)^\times:
             r^5x^2\equiv ay^2+bz^2\pmod m\}.                  \tag{PO2}
\]
For \(m=1\) use its one residue class. There are
\(O_\varepsilon((ZV)^\varepsilon)\) such classes. With an integer
representative \(r\), enumerate its possible occupants as
\(t=r+km\in[Z,2Z)\), and retain only distinct odd primes
\(t\ne p_1,p_2\). In particular \((t,yz)=1\).

**Lemma PO0 (Exact field test).** For each such occupant the traces in
(5), with \(P_3=t^5\), are integers, and
\[
 E(t)=\frac{4ab-u_{12}(t)^2}{x^2}
 =\frac{4ab\,y^2z^2-(ay^2+bz^2-t^5x^2)^2}{x^2y^2z^2}
                  \in\mathbb Z\setminus\{0\}.                  \tag{PO3}
\]
Write uniquely \(E(t)=-\iota d\delta^2\), where \(d,\delta>0\)
and \(d\) is squarefree. Thus the field-window condition is exactly
\(X\le d<2X\), with \(d\ne1\) for \(\iota=1\).

The two reductions of (PO2), together with (PO1), give (6), hence
trace integrality. Set \(A=ay^2,B'=bz^2,C=t^5x^2\). The numerator
in (PO3) is
\(-(A-B')^2+2C(A+B')-C^2\), divisible by \(x^2\).
Since \((x,yz)=1\) and the trace is integral, this proves \(E\) is
integral. It is nonzero because \(4(p_1p_2)^5\) is not a square.
The other two norm defects are the same \(E\), by expansion of (5).
The traces and coefficients \(\delta(x,y,z)\) consequently give
integral algebraic numbers of the required norms and coherent product.
If also \((\delta,p_1p_2t)=1\), their incident coefficient units give
the one-sided prime ideal factorizations: a supporting prime cannot
divide both the element and its conjugate. Its ramification would
force its norm valuation to be one, contrary to the fifth-power norm.
The coherent product forces opposite orientations at the middle prime.
This proves the field test and the stated norm interpretation.

Define \(W_{\iota,X}(t)\) to be one when this squarefree-part test and
\((\delta,p_1p_2t)=1\) hold, and when
\[
 \max_{ij}|u_{ij}|\le BZ^5,\qquad
       |E(t)|\max(x^2,y^2,z^2)\le B^2Z^{10}.                    \tag{PO4}
\]
Otherwise put it zero. These are explicit necessary conditions for
the selected witnesses, once \(B\) is fixed large enough. Their sum
is a majorant; passing them does not assert membership in the particular
real unit interval or the previously chosen path selection. Those
additional restrictions can be retained if needed.

This retains the squarefree part of \(E\), rather than replacing it
by \(|E|\). In a balanced coefficient box, rescaling
\(p_i=Zr_i\) and \((x,y,z)=V(s_x,s_y,s_z)\) gives
\(E=(Z^{10}/V^2)\mathcal E(r,s)\). Even with \(\delta=1\), the
relaxed condition \(1<|\mathcal E|<2\) contains open real regions
of both signs. For example take all \(r_i=r\in(1,2)\),
\(s_x=1,s_y=s_z=s\), so
\(\mathcal E=r^{10}(4-s^{-2})\), and then perturb inside the box.
This only checks the archimedean region; it makes no assertion about
the density of integer, prime or squarefree records. Magnitude alone
does not furnish the desired power saving by an interval-length argument.

### Small normalized coefficients already meet the target

**Proposition PO1 (Small-coefficient removal).** For a fixed odd prime
exponent, \(1\le U\le BV\), and either sign, the coherent triangles
with \(\min(|x|,|y|,|z|)\le U\) number
\[
        O_{B,\ell,\varepsilon}((ZV)^\varepsilon
                                  (Z^2UV+Z^3)).                 \tag{PO5}
\]
The bound is uniform under additional field-window and path restrictions.

Permute the vertices so a small coefficient is \(y\), leaving it
free in the divisor argument. For fixed \(p_1,p_2,y,z\), (PO1)
gives divisor-many \(x\). The remaining prime costs at most a
subpower multiple of \(1+Z/|yz|\). Summing \(0<|y|\le U\),
\(0<|z|\le BV\) gives \(UV+O(Z\log(2U)\log(2BV))\).
The first two primes cost \(O(Z^2)\). The finite choices of vertex
permutation and signs complete the proof. If \(U<1\) the sector is empty.

In particular \(U=V/Z^{1/3}\) removes this sector at cost
\(O_\varepsilon((ZV)^\varepsilon(Z^{5/3}V^2+Z^3))\).
The residual has **all three** coefficients greater than this threshold;
they need not lie in one dyadic box.

### A uniform open estimate and its fourth-moment consequence

Let \(\mathcal O_{5,\iota}^{\rm res}(Z,V;X)\) be the sum, over
the divisor records (PO1) with
\(\min(|x|,|y|,|z|)>V/Z^{1/3}\), of
\[
 \sum_{r\in R(p_1,p_2;x,y,z)}
 \ \sum_{\substack{k\in\mathbb Z:\ t=r+k|yz|\in[Z,2Z)\\
                  t\text{ odd prime},\ t\ne p_1,p_2}}
                                      W_{\iota,X}(t).           \tag{PO6}
\]
The first two primes are distinct odd primes in \([Z,2Z)\).
This is an explicit sufficient majorant for the selected residual paths.

**Open estimate PO2 (Uniform residual prime occupancy).** For each
fixed \(B\) and \(\varepsilon>0\), prove
\[
 \mathcal O_{5,\iota}^{\rm res}(Z,V;X)
       \ll_{B,\varepsilon}(ZV)^\varepsilon
                             (Z^{5/3}V^2+Z^3)                  \tag{PO7}
\]
for \(X=Z^{10}/V^2\ge2\), both signs, and the entire band
\[
                         Z^{17/10}\le V\le Z^{9/5}.             \tag{PO8}
\]
It would also suffice to prove this with the original selection
restrictions retained. An estimate at \(V=Z^{7/4}\) alone does not
supply (PO7). Throughout the residual band each coefficient is
\(>Z^{41/30}\), and \(|yz|>Z^{41/15}\). For sufficiently large
\(Z\), every class therefore has at most one integer occupant in
the prime interval. The missing bound concerns how often that occupant
exists and passes the prime and norm tests, averaged over these
correlated divisor records.

**Conditional corollary PO3 (Fourth moment from occupancy).** If
(PO7) holds throughout (PO8) for a fixed sign, then
\[
 \sum_{\substack{0<d\le X\\ d\text{ squarefree},\ \iota d\ne1}}
       h_5(\iota d)^4\ll_\varepsilon X^{59/29+\varepsilon}.      \tag{PO9}
\]
This implication does not assert (PO7) or improve the current
unconditional exponent \(53/26\).

Combine (PO5) and (PO7), divide by the path lower bound (4), and
substitute \(Z=X^{1/2+4\eta}/H,V=Z^5/\sqrt X\). In (PO8) this gives
\[
 N_{5,\iota}(H;X)\ll_\varepsilon X^\varepsilon
            (1+XH^{-2}+X^{13/3}H^{-32/3}).                     \tag{PO10}
\]
Choose the auxiliary \(\eta>0\) sufficiently small for the final
\(\varepsilon\), as in (1). Ignoring these absorbable subpowers,
the band corresponds to
\(11/32\le h=\log_XH\le23/66\). The new tail crosses \(X/H\)
at \(h=10/29\), strictly inside the band, where
\(V=Z^{16/9}\). Multiplying the tail by \(H^4\) gives the exponent
\(1+3(10/29)=59/29\) at this crossing and a decreasing main
contribution above it.

Below the crossing the pair tail, together with the trivial count
below \(H=X^{1/4}\), gives at most
\(\max(2,1+3(10/29))\). The decreasing \(X^{5/2}H^{-6}\)
pair term contributes at most \(X^2\). Above the occupancy band,
use the original fifth-power triangle tail (1). Its rank-one and
rank-two fourth-moment contributions at \(h=23/66\) have exponents
\[
 \frac{11}{2}-10\frac{23}{66}=\frac{133}{66},\qquad
 4-\frac{23}{4}\frac{23}{66}=\frac{527}{264},
\]
both strictly below \(59/29\), and both decrease with \(h\).
The term \(XH^{-2}\) contributes at most \(X^{9/5+\varepsilon}\)
up to the nonexceptional cap \(H\le X^{2/5+\eta}\).
Exceptional fields cost \(X^{2+\varepsilon}\), as before.
The strict margins at the band endpoints absorb the \(\eta\)
allowances. Dyadic summation proves (PO9). The improvement would be
\(53/26-59/29=3/754\). Frey's comparison in the
[five-Selmer application](X011_FIVE_SELMER_MOMENTS_FROM_QUADRATIC_TORSION.md#transferred-moments)
would then transfer the negative-sign version to its second moment.

More generally, for \(0<\theta\le1\), a bound
\(Z^{2-\theta}V^2+Z^3\) gives the tail term
\(X^{(9-\theta)/2}H^{-(11-\theta)}\), with crossing
\(b_\theta=(7-\theta)/(20-2\theta)\). This is below \(9/26\)
exactly when \(\theta>1/4\), and its crossover is
\(V=Z^{(5+\theta)/3}\). This algebra presupposes a bound uniform
around that new scale and enough surrounding parameter coverage;
it does not establish another moment theorem by itself.

### Primary-source comparison and exact controls

The closest retained pair argument is Heath-Brown--Pierce,
Proposition 2.3 and Section 3, in
[arXiv:1409.3177v2](sources/heath-brown-pierce-1409.3177v2.pdf).
It uses a unit square-root congruence and the trace interval to count
two-prime norm records. Chan--Koymans,
[arXiv:2505.00611v1, Proposition 5.1](https://arxiv.org/pdf/2505.00611v1),
supplies the analogous real quadratic pair bound with bounded witnesses.
Neither statement counts the third-prime occupancy in (PO6).

The modular-root comparison inspected here is Shkredov--Shparlinski--Zaharescu,
[Theorem 1.1, published pp. 2--3](https://doi.org/10.1007/s00209-024-03436-5).
It bounds average discrepancy for \(s^2\equiv\lambda p\pmod q\),
over prime moduli \(q\le Q\), by
\(P(P^{-1/12}+P^{-1/5}Q^{1/10})Q^{o(1)}\), for \(1\le P\le Q\).
This concerns quadratic roots with an unweighted modulus average.
Our classes are fifth roots for composite \(|yz|\), with divisor and
field weights coupled to the same variables. No reduction to its
hypotheses is established. With \(Q\asymp V^2,P\asymp Z\) at the
critical scale, even its displayed error is not a saving over \(P\).

For an arithmetic-progression comparison, Li,
[arXiv:2611.29075v1, Theorem 1.1](https://arxiv.org/pdf/2611.29075v1),
assumes a fixed nonzero residue \(a\), divisor-bounded product weights,
and in particular \(Q_1Q_2\le x^{8/15-\varepsilon}\), with three
additional restrictions on the factors. Here the prime scale is
\(x\asymp Z\), the classes vary with the divisor record, and already
\(|yz|>Z^{41/15}\) in the residual band. That theorem does not
directly apply. These are scoped comparisons, not a complete priority
survey or proof that other methods cannot establish (PO7).

[The occupancy checker](verify_fifth_prime_occupancy.py) uses only the
Python standard library. It checks the two signed norm fixtures,
rejects corrupted prime tuples, compares CRT enumeration with direct
enumeration at \(Z=10,V=7\), and checks the rational crossing and
endpoint exponents. Its census counts compatible
classes, integer occupants, distinct prime occupants and the exact
field-window survivors separately, using positive coefficients only.
The canonical command is
`python3 research/arithmetic/verify_fifth_prime_occupancy.py`.
The direct control matches 402 tuples. At \(Z=10,V=56\), the census
has 56,554 divisor records, 57,262 compatible classes, 6,776 integer
occupants and 1,320 distinct-prime occupants; none reaches the residual
coefficient sector or the prescribed field window. At the exact critical
scale \(Z=16,V=128,X=2^{26}\), the command with `--Z 16 --V 128`
gives 616,886 divisor records, 625,852 compatible classes, 29,212 integer
occupants, 5,466 distinct-prime occupants and twelve residual prime
occupants. None passes the field window. One residual example is
\((p_1,p_2,p_3)=(17,23,29)\), \((x,y,z)=(57,113,53)\), whose
squarefree defect is \(3,304,670,587>2X\). A separate norm fixture
over \(\mathbb Q(\sqrt{-3})\) checks acceptance of the field test
and rejection after changing the window; it is not a dyadic prime-box
fixture. These finite controls test reconstruction, coverage and rejection,
and supply no asymptotic occupancy exponent or new moment bound.
