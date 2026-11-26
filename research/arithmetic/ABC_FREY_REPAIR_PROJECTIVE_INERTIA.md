# Frey inertia under minimal target repair

A common multiplier preserves the finite-level projective Galois representation of the Frey curve. The inertia orders at primes away from the level and \(2\) measure the actual primitive valuation defects. Making the repaired target a high power preserves these inertia orders.

## Projective representations and inertia orders

For a positive primitive triple \(a+b=c\), use
\[
                  E_{a,b}\colon y^2=x(x-a)(x+b).
\]
For \(m\ge2\), let \(\mathbb P\rho_{E,m}\) denote the action on \(E[m]\) modulo scalar matrices in \(\mathrm{GL}_2(\mathbb Z/m\mathbb Z)\). For any positive integer \(Z\),
\[
             \mathbb P\rho_{E_{Za,Zb},m}
                    \simeq \mathbb P\rho_{E_{a,b},m}.          \tag{FR1}
\]
Indeed \((x,y)\mapsto (Zx,Z\sqrt Z\,y)\) is an isomorphism over \(\mathbb Q(\sqrt Z)\). Its Galois conjugates differ by the involution \([-1]\) according to the quadratic character of \(Z\), so the torsion action is twisted by that scalar character. Quotienting by the scalars proves (FR1), also at even levels. The projective torsion fields are equal in a fixed algebraic closure. [Siksek, §§7–9](https://samirsiksek.github.io/siksek.github.io/modcurves2.pdf#page=10)

For a rational prime \(\ell\nmid 2m\), put \(n_\ell=v_\ell(abc)\). For either curve in (FR1) the projective inertia image has order
\[
                      \frac{m}{\gcd(m,2n_\ell)} .           \tag{FR2}
\]
At \(n_\ell=0\) the right side is one by the convention \(\gcd(m,0)=m\). The primitive model has discriminant \(16(abc)^2\) and \(c_4=16(a^2+ab+b^2)\). At \(n_\ell>0\), exactly one of \(a,b,c\) is divisible by \(\ell\), so \(c_4\) is a unit. This gives a minimal model with multiplicative reduction and minimal discriminant valuation \(2n_\ell\). Off the three supports it has good reduction.

In the multiplicative case pass to the unramified extension of degree at most two that makes the reduction split. The Tate parameter \(q_{\rm Tate}\) has valuation \(2n_\ell\). Over the maximal unramified extension, \(\mu_m\) and the \(m\)-th roots of units are available because \(\ell\nmid m\). Adjoining \(q_{\rm Tate}^{1/m}\) therefore has ramification degree \(m/\gcd(m,2n_\ell)\). The Tate basis of torsion given by a primitive root of unity and \(q_{\rm Tate}^{1/m}\) makes this inertia action unipotent, with a scalar matrix only for the identity. This proves (FR2) for the primitive model. Use (FR1) to get the same order on the multiplied triple. [Tate uniformization, Silverman Appendix C, §14](https://www.math.ens.psl.eu/~obenoist/refs/Silverman.pdf#page=445) · [Siksek, Corollary 9.1](https://samirsiksek.github.io/siksek.github.io/modcurves2.pdf#page=12)

## Relation to a fixed power order

For a fixed integer \(q\ge2\), write every primitive letter exactly as \(a_i=A_i x_i^q\) with \(0\le v_\ell(A_i)<q\). The coefficient radical \(\operatorname{rad}(A_1A_2A_3)\) records precisely the primes where the actual valuation on the unique supported letter fails divisibility by \(q\). At torsion level \(2q\), (FR2) says that the ramification support away from \(2q\) is exactly this coefficient radical away from \(2q\). At torsion level \(q\), the same assertion holds for odd \(q\). At even level \(q\) the support in this formula tests divisibility by \(q/2\). These are pointwise identities at fixed level and require primitivity of the original complete letters.

For an odd prime \(p\), there is also a linear conductor identity. The conjugation representation \(\operatorname{Ad}^0\rho_{E,p}\) on the three-dimensional trace-zero matrix space is unchanged by the twist. Away from \(2p\) its conductor is
\[
                 \prod_{\substack{\ell\mid abc,\ \ell\nmid2p\\
                                    p\nmid v_\ell(abc)}}\ell^2 .
                                                               \tag{FR3}
\]
At any prime in the product, inertia acts through a nontrivial unipotent group of order \(p\). The trace-zero subspace that centralizes it has dimension one, so tame conductor exponent is two. The other primes away from \(2p\) have trivial inertia. This calculation uses the torsion action itself.

In the least target repair use \(Z=Z_M(c)=\prod_{\ell\mid c}\ell^{(M-v_\ell(c))_+}\). If \(q\mid M\), every repaired target valuation below the cutoff becomes divisible by \(q\). For any such prime \(\ell\nmid 2q\) with \(q\nmid v_\ell(c)\), (FR2) gives nontrivial inertia at level \(2q\) on the repaired Frey curve. For example the primitive triple \(1+4=5\) has projective inertia of order \(3\) at \(5\) on its \(3\)-torsion. The least order-\(3\) repair is \(25+100=125\) and has the same order-\(3\) inertia there. Ordinary residual level lowering on a primitive power equation uses the original discriminant valuations at its multiplicative primes. The repaired exponents cannot supply those divisibility hypotheses. [Least target repair](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target) · [Dieulefait–Soto, §3](https://arxiv.org/pdf/1811.05957v2)
