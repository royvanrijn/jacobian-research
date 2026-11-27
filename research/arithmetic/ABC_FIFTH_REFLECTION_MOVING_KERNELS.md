# Fifth reflection and moving quadratic kernels

Let \(h_5(L)=|\operatorname{Cl}(L)[5]|\), with the ordinary ideal class
group. For nonzero signed squarefree \(\kappa\), set
\(\kappa^*=\operatorname{sf}(5\kappa)\), where \(\operatorname{sf}\)
retains the sign, and put
\[
 k_\kappa=\mathbf Q(\sqrt\kappa),\qquad
 B_\kappa=\mathbf Q(\zeta_5,\sqrt\kappa),\qquad
 F_{a,\kappa}=\mathbf Q(\zeta_{5^a},\sqrt\kappa).
\]
The two radicands \(1,5\) give fixed fields and can be separated in every
sum. This is the family used by the
[degree-five cyclic isogeny descent](ABC_ROOT_RATIO_MIXED_COVERS.md#cyclic-isogeny-descent-for-2p2pa)
for signature \((10,2,5^a)\).

## Class ranks at the first layer

**Proposition FR5.** For \(\kappa\ne1,5\), the first layer satisfies
\[
 h_5(B_\kappa)
   =5^{\,2{\bf1}_{\kappa>0}-\delta_1-\delta_2}
      \big(h_5(k_\kappa)h_5(k_{\kappa^*})\big)^2,\qquad
 \delta_1,\delta_2\in\{0,1\}.
 \tag{1}
\]
Here the two \(\delta_i\) are the local terms in Kishi's reflection
theorem. For fixed \(a\ge1\),
\[
 h_5(F_{a,\kappa})\ll_a
       \big(h_5(k_\kappa)h_5(k_{\kappa^*})\big)^{2\cdot5^{a-1}} .
 \tag{2}
\]

Take \(d\) to be the member of \(\{\kappa,\kappa^*\}\) prime to 5.
Then \(d\ne1\), and \(k_1=\mathbf Q(\sqrt d)\) and
\(k_2=\mathbf Q(\sqrt{5d})\) are the two moving quadratic subfields.
The cyclic quartic mirror is
\(M=\mathbf Q((\zeta_5-\zeta_5^{-1})\sqrt d)\).
Write \(\tau\) for the generator of the cyclotomic group acting by the
second power on \(\zeta_5\), and \(\sigma\) for the involution acting
nontrivially on \(k_1\). The respective fixing subgroups of
\(k_1,k_2,M\) in \(G=\langle\tau\rangle\times\langle\sigma\rangle\)
are \(\langle\tau\rangle\), \(\langle\tau\sigma\rangle\), and
\(\langle\tau^2\sigma\rangle\).

For any subgroup \(H\le G\) the extension and ideal norm on the
5-primary class groups satisfy \(Nj=[B_\kappa:B_\kappa^H]\) and
\(jN=\sum_{h\in H}h\). The degree is invertible on these groups, so the
extension identifies \(\operatorname{Cl}(B_\kappa^H)[5]\) with the
\(H\)-invariants. The action splits over \(\mathbf F_5\). Index
characters by \((i,e)\) with \(\tau\mapsto2^i\) and
\(\sigma\mapsto(-1)^e\). All \(e=0\) components vanish because
\(\mathbf Q(\zeta_5)\) has class number one (the Minkowski bound is
less than two). The \((0,1)\) and \((2,1)\) components give the two
quadratic ranks. The \((1,1)\) and \((3,1)\) components give the mirror
rank. In particular
\[
 r_5(B_\kappa)=r_5(k_1)+r_5(k_2)+r_5(M).
\]
Kishi's Main Theorem for \(d\ne1\) squarefree and prime to 5 gives
\(r_5(M)=r_5(k_1)+r_5(k_2)+2{\bf1}_{d>0}-\delta_1-\delta_2\),
which proves (1) for both signs, including when \(5\mid\kappa\) by
interchanging the quadratics.
[Kishi, Main Theorem, pp. 1–2](https://www.math.okayama-u.ac.jp/mjou/mjou47/_01_Kishi.pdf)

Let \(G_{B,S}\) be the group of the maximal pro-5 extension of
\(B_\kappa\) unramified away from the primes \(S\) above 5.
Since \(B_\kappa\) contains \(\mu_5\), Kummer theory bounds
its generator rank by \(r_5(B_\kappa)+O(1)\). Indeed, the
classes of radicands have valuations divisible by 5 off \(S\), and
their dimension is
\(\dim(\mathcal O_{B_\kappa,S}^\times/
(\mathcal O_{B_\kappa,S}^\times)^5)+
\dim\operatorname{Cl}(\mathcal O_{B_\kappa,S})[5]\).
The first term is bounded by degree and \(\#S\), and the second
by the ordinary 5-rank. The index of the subgroup over
\(F_{a,\kappa}\) is \(5^{a-1}\), because \([B_\kappa:\mathbf Q]\)
is prime to 5 and the cyclotomic increment is cyclic of this degree.
This subgroup is the full maximal pro-5 group unramified away from 5
over \(F_{a,\kappa}\). The conjugate compositum of any finite Galois
5-extension of \(F_{a,\kappa}\) unramified away from 5 is still a
5-extension over \(B_\kappa\) with the same ramification restriction.
Schreier bounds the upper generator rank by
\(5^{a-1}r_5(B_\kappa)+O_a(1)\). The class group rank at the
upper layer is at most this generator rank. Equation (2) follows.

## Transfer of moments and the support budget

**Corollary FR6.** Suppose that for the sign under consideration the
quadratic \(5\)-torsion moment of order \(4\cdot5^{a-1}\nu\)
is \(O_\epsilon(X^{\Theta+\epsilon})\) when ordered by absolute
discriminant, where \(\nu>0\) and \(\Theta\ge1\). Then
\[
 \sum_{\substack{1\le|\kappa|\le X\\
                 \kappa\ \mathrm{squarefree}}}
        h_5(F_{a,\kappa})^\nu\ll_{a,\nu,\epsilon} X^{\Theta+\epsilon}
 \tag{3}
\]
over that sign. In particular, for \(a=1\) the signed first-layer
moments at orders \(3/4,1,5/4\) are bounded with exponents
\(7/4,53/26,5/2\), respectively, up to \(\epsilon\).

To see (3), use (2) and Cauchy–Schwarz in \(\kappa\). The involution
\(\kappa\mapsto\kappa^*\) preserves the sign, is injective, and
satisfies \(|\kappa^*|\le5|\kappa|\). Both quadratic discriminants
therefore have absolute value at most \(20X\) throughout the sum.
For the three examples the needed moments are of orders \(3,4,5\)
and supplied by the
[signed coherent prime-path bounds](ABC_COHERENT_PRIME_PATH_HIGHER_ODD_TORSION.md).

One may also order the nonexceptional first layer by
\(|\operatorname{Disc}(B_\kappa)|\). The four characters of \(G\)
nontrivial on \(\sigma\) all have conductor divisible by the
prime-to-10 part of \(|\kappa|\), exactly to the first power at
each such prime. The conductors of the cyclotomic characters and the
conductor factors at 2 and 5 are bounded, since the local characters
come from \(\mathbf Q(\zeta_5)\) and a quadratic character. The
conductor-discriminant formula gives
\(|\operatorname{Disc}(B_\kappa)|\asymp|\kappa|^4\). The two possible
radicands for each field interchange under \(\kappa\mapsto\kappa^*\).
Consequently the first moment on the distinct nonexceptional fields
with absolute discriminant at most \(Y\) is
\(O_\epsilon(Y^{53/104+\epsilon})\).

For (RR6), write \(r=\operatorname{rad}(AB)=AB\) in this signature.
Each \(r\) has \(2^{\omega(r)}\) allocations between the two squarefree
kernels and \(\kappa\) has absolute value \(r\). For a fixed sign,
dyadic summation and (3) bound the weighted moment with order \(\nu\)
by \(O_{a,\nu,\epsilon}(Z^{\Theta-1+\epsilon})\); the allocation
factor is absorbed in \(Z^\epsilon\). Thus a near-linear quadratic
moment at order \(4\cdot5^{a-1}\nu\), for the fixed order \(\nu\)
required by (RR5), suffices for the moving-kernel defect-budget bound
in signature \((10,2,5^a)\). The new unconditional first-layer
first moment only bounds its weighted cost with exponent \(27/26\).

At \(a=1\), (1) in fact identifies the exact arithmetic moment behind
(RR6) for every fixed \(\nu>0\). The weighted estimate with exponent
\(\epsilon\) for every \(\epsilon>0\) holds for a prescribed sign if
and only if
\[
 \sum_{\substack{1\le|\kappa|\le X\\
                 \kappa\ \mathrm{squarefree}}}
     \big(h_5(k_\kappa)h_5(k_{\kappa^*})\big)^{2\nu}
       \ll_{\nu,\epsilon} X^{1+\epsilon}
 \tag{4}
\]
over that sign. Sufficiency follows by the allocation bound and
dyadic summation. Necessity follows by retaining one allocation for
each \(r\), using the lower bound in (1) and multiplying the weighted
sum by \(Z\). The two quadratics vary with the same radicand up to a
fixed twist. A bound for one field or a bounded catalogue of fields
would leave (4) open, as would the unconditional bound (3).

## Exact illustrations

Kishi's Examples 7.1 and 7.2 give
\((r_5(k_1),r_5(k_2),r_5(M))=(0,1,3)\) at \(d=2723\) and
\((2,0,2)\) at \(d=-14606\). Thus the constant \(25\) in the
positive upper bound (1) is attained, and in the negative example
the upper bound with constant \(1\) is attained. These ranks can be
replayed with the defining polynomials \(x^2-d,\ x^2-5d,\
x^4+5dx^2+5d^2\). A defining polynomial for \(B_\kappa\) is given
by composing \(x^2-d\) with \(x^4+x^3+x^2+x+1\).
[Kishi, Examples 7.1–7.2](https://www.math.okayama-u.ac.jp/mjou/mjou47/_01_Kishi.pdf)
