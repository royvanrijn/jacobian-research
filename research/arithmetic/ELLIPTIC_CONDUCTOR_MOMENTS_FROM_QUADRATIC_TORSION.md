# Elliptic-curve conductor moments from quadratic 3-torsion

Write \(E(q)\) for the number of \(\mathbb Q\)-isogeny classes of elliptic curves of conductor \(q\). Let \(\gamma=\log(3)/(2\beta(0))=1.9745\ldots\) with \(\beta(0)\) as in Helfgott–Venkatesh, Theorem 3.8, and [Pierce–Turnage-Butterbaugh–Wood, Theorem 1.12](sources/pierce-turnage-butterbaugh-wood-1902.02008v3.pdf).

## Discriminant conversion and the transferred moments

Theorem 1.12 requires a bound \(\ll_{k,\varepsilon} X^{\Theta_k+\varepsilon}\) for the \(k\)-th moment of quadratic class-group \(3\)-torsion ordered by absolute **field discriminant**, for a fixed real \(k\ge 1/\gamma\) and \(\Theta_k\ge1\). It gives the same exponent in \(Q\) for \(\sum_{q\le Q}E(q)^{\gamma k}\). Its proof in §6.3 actually bounds the larger number of isomorphism classes using the same exponent.

Every quadratic field has a unique signed squarefree radicand \(d\ne1\). Its fundamental discriminant \(D\) is \(d\) for \(d\equiv1\bmod4\), and \(4d\) otherwise. Thus \(|D|\le X\) gives \(|d|\le X\). Positivity lets us apply the [signed radicand moment bound](ABC_QUADRATIC_THREE_TORSION_FERMAT_TAIL.md#moment-consequences-and-remaining-terms) without a lower cutoff or congruence restriction. The finitely small \(X\) are absorbed in the constant. In particular \(k=2,4\) are both admissible. The resulting bounds, for \(Q\ge2\) and every \(\varepsilon>0\), are
\[
 \sum_{q\le Q}E(q)^{2\gamma}\ll_\varepsilon Q^{17/14+\varepsilon},
 \qquad
 \sum_{q\le Q}E(q)^{4\gamma}\ll_\varepsilon Q^{23/14+\varepsilon}.
\]
More generally, for each fixed \(k\ge1\), the exponent supplied by the same argument is \(\max(1+3(k-1)/14,\,1+(1597/5000)(k-2))\). These are proof consequences of the retained prime-path argument. The source review and publication status of that argument are inherited only as dependencies. In particular, fixed moments supply bounds on the number of unusually large conductors by Markov's inequality. The Brumer–Silverman pointwise conjecture \(E(q)\ll_\varepsilon q^\varepsilon\) remains open.

## Comparison at the same moments

Corollary 1.13 of the transfer paper gives \(23/18\) at \(2\gamma\) and \(11/6\) at \(4\gamma\), using [Heath-Brown–Pierce, Corollary 1.4](https://arxiv.org/pdf/1409.3177v2). The decreases from those two exponents are \(4/63\) and \(4/21\), respectively.

Theorem 1.1 of [Koymans–Thorner](https://arxiv.org/pdf/2309.07204) specializes to exponents \(3/2\) and \(2\) at quadratic \(3\)-torsion orders \(2\) and \(4\), respectively, using the quadratic-field count \(O(X)\). The improved pointwise bound \(\kappa=1597/5000\) of [Chan–Koymans, Theorem 1.1](https://arxiv.org/pdf/2505.00611) can be interpolated with the first moment and with Heath-Brown–Pierce. For example, interpolation from the second quadratic moment would bound the fourth moment with exponent \(23/18+2\kappa\), above both the old \(11/6\) and the new \(23/14\). Chan–Koymans's real-field average in Theorem 1.2 concerns primes \(\ell\ge5\). [Shankar–Shankar–Wang, Theorem 1.1](https://www.math.uwaterloo.ca/~x46wang/Papers/Conductor.pdf) counts restricted families by bounded conductor, so its exponents cannot be compared as upper bounds for these unrestricted conductor moments. The comparison here establishes the improvement over the cited corollary. Historical priority beyond the inspected statements remains to be checked.

## Internal checks on the prime-path input

In the witness construction the ideal of \(\bar\beta_i\beta_j/p_0^3\) is \((\bar{\mathfrak p}_i\mathfrak p_j)^3\), and \(\alpha_{12}\alpha_{23}=p_2^3\alpha_{13}\) follows from \(\beta_2\bar\beta_2=p_0^3p_2^3\) including when the middle vertex is the anchor. The signed coefficient minors satisfy the Plücker relation. Three products of equal edge colors would share a common scalar \(c^2\), where \(c\) is the common positive cube-free part, and canceling \(c^2\) gives the signed Fermat equation. Counting nonmonochromatic four-subsets then counts adjacent edge pairs with multiplicity at most \(s-3\) in a coset of size \(s\). The split-prime estimate in Heath-Brown–Pierce (2.5)–(2.7) is used at a common \(Z\) in each \((X,H)\) box, with \(Z\gg X^{1/6}\). Summing a moment uses only logarithmically many \(H\)-boxes, or the pointwise bound on the exceptional portion.

For the full-modulus count, the even part of a \(q_i\) requires \(v_2(y_i)\ge3\). In the corresponding norm congruence the left side has valuation two, forcing \(v_2(j)=1\). Dividing by four at 2 yields a unit cube congruence. An even \(r\) instead forces both \(y_i\) odd. Thus no two retained moduli share a 2-part. The cofactor constraints in the rank-two branch come from substituting the congruences into integer minors. They do not require reducing a basis of the plane to a linearly independent pair modulo a composite modulus. After division by the common gcd of the cofactor parameters, scaling the primitive coefficient-plane lattice by the reciprocal box widths gives the covolume in the proof. The primitive modulus triple bound in a line is needed to handle the first minimum separately from the rank-two bound. These checks are internal. Independent review of the lattice count, the witness construction and the complete import chain remains open.
