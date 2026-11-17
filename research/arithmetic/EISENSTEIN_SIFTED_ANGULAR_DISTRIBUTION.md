# Sifted Eisenstein norm shells and their angular distribution

Let \(R=\mathbf Z[\omega]\), where \(\omega^2+\omega+1=0\), and put \(N(z)=z\bar z\). Let \(A\) consist of the positive integers whose prime divisors are all \(1\bmod3\), including \(1\). The notation \(3\cdot A\) means \(\{3a:a\in A\}\). For \(n\equiv2\bmod6\), let \(G_n\) consist of the ordered pairs \((z,w)\in R^2\) with \(N(z)+N(w)=n\) and \(N(z),N(w)\in A\). For \(n\equiv4\bmod6\), require \(N(z)\in A,\ N(w)\in3\cdot A\), fixing which component carries the factor of 3.

**Theorem EisAngle.** As \(n\) tends to infinity along either class, the unit vectors \((z,w)/\sqrt n\) for \((z,w)\in G_n\), each counted once, become equidistributed in the normalized surface measure on \(S^3\subset\mathbf C^2\). For each fixed interval \(I\subset[0,1]\) of positive length, the proportion with \(N(z)/n\in I\) tends to \(|I|\). For every \(\epsilon>0\), such witnesses give \(\gg_{I,\epsilon} n^{1-\epsilon}\) distinct ordered pairs of positive summand values. In particular, each sufficiently large target in either class admits a decomposition with both summands between \(n/3\) and \(2n/3\).

The proof uses the all-ideal Hecke norm estimate and torus-normalizer averaging for odd inert primes from the unpublished binary classification. It imposes a fixed finite set of local conditions before applying harmonic theta series. A uniform estimate for the primes above that finite set then gives angular equidistribution after the full sieve. [Norm count](../../unpublished/eisenstein-restricted-prime-sums/sections/norm-count.tex), [Hecke lemma and orbit averaging](../../unpublished/eisenstein-restricted-prime-sums/sections/prime-bounds.tex), [odd-component filter](../../unpublished/eisenstein-restricted-prime-sums/sections/mixed-local.tex), [large-prime and divisor bounds](../../unpublished/eisenstein-restricted-prime-sums/sections/elementary-estimates.tex).

## Proposition EisEta. The obstruction at 3

Write \(a_j=\sum_{z\in R}q^{jN(z)}\), \(q=e^{2\pi i\tau}\). Writing \(x-y\omega\) instead of \(x+y\omega\) interchanges the two signs of the mixed term in the quadratic form. Multiplication by \(1-\omega\) maps \(R\) bijectively to the norm-divisible-by-3 sublattice, multiplying norms by 3. Multiplication by 3 maps \(R\) bijectively to the norm-divisible-by-9 sublattice, multiplying norms by 9. Hence \(a_1-a_3\) counts precisely the norm-1-modulo-3 components, including no zero. Write \(r(n)\) for the total ordered, signed quaternion shell count, \(B_3(n)\) for the number with some component norm divisible by 9, and \(F(\tau)=\eta(3\tau)^2\eta(9\tau)^2=\sum_{n\ge1}c_nq^n\), normalized with \(c_1=1\).

If \(n\equiv1\pmod3\), exactly one component norm is 1 modulo 3. The other is divisible by 3. In particular the two possible orientations for divisibility by 9 are disjoint, including representations with a zero component. Thus, *as series supported on* \(n\equiv1\pmod3\),
\[
 \sum r(n)q^n=2(a_1-a_3)a_3,\qquad
 \sum B_3(n)q^n=2(a_1-a_3)a_9 .
\]
In Borwein, Borwein and Garvan's cubic notation, Lemma 2.1(ii),(iii) gives \(a_1-a_3=2c(q^3)\), \(3a_9-a_3=2b(q^3)\). Proposition 2.2(i),(ii) gives \(b(q)=\eta(\tau)^3/\eta(3\tau)\), \(c(q)=3\eta(3\tau)^3/\eta(\tau)\). Consequently
\[
 \sum_{n\equiv1(3)}(B_3(n)-r(n)/3)q^n
   ={2\over3}(a_1-a_3)(3a_9-a_3)
   ={8\over3}b(q^3)c(q^3)=8F(\tau).
\]
Since \(r(n)=12\sigma(n)\) for \(3\nmid n\), independently obtained by counting ideals of the maximal class-number-one order, this proves, for **every positive \(n\equiv1\bmod3\)**,
\[
  B_3(n)=4\sigma(n)+8c_n . \tag{1}
\]

Martin and Ono, Theorem 2 (level 27 row) and Theorem 3, identify \(F\) as a normalized weight-two, level-27 newform with trivial character and \(c_p=0\) for all \(p\equiv2\bmod3\), including 2. They also give the split coefficients as a CM trace \(2m+r\) with \(p=m^2+mr+r^2\), hence \(|c_p|\le2\sqrt p\). For \(p\ne3\) the Euler denominator is \(1-c_pX+pX^2=(1-\alpha_pX)(1-\beta_pX)\). Its parameters \(\alpha_p,\beta_p\) have absolute value \(\sqrt p\) by this elementary trace bound (including the equal-root case). Multiplicativity gives \(|c_n|\le\tau(n)\sqrt n\) for every \(3\nmid n\). For an all-inert integer \(n\), the local reciprocal factor is \((1+pX^2)^{-1}\). Its coefficients are zero at odd exponents, and \((-p)^{a/2}\) at even exponents \(a\). This proves exactly the binary classification's reduced-target formula \(c_n=\lambda_n\), independently of its six-point and Jacquet–Langlands arguments. In particular \(c_7=-1\), giving \(B_3(7)=24\) by (1). Using the all-inert nonsquare formula at 7 would incorrectly give 32.

For even \(n\equiv1\bmod3\), the two component norms have the same parity. Both are even precisely when the quaternion lies in \(2(R\oplus R)\). Let \(B_3^o\) count the ramified obstruction among the odd-component shell points. Scaling by 2 preserves divisibility by 9. The number of odd-component shell points is \(R_o(n)=12S(n)\), where \(S(n)=\sigma(n)-\sigma(n/4)\), with every quantity at \(n/4\) interpreted as zero if \(4\nmid n\). Then
\[
 B_3^o=4S+8(c_n-c_{n/4}),\qquad
 B_3^o/R_o=\tfrac13+O_\epsilon(n^{-1/2+\epsilon}).
\]
Here \(R_o\ge18n\) since \(S=3\cdot2^{e-1}\sigma(m)\) when \(n=2^e m\) and \(m\) is odd. This establishes an unrestricted mixed main term. For all-inert square \(n\), \(c_{n/4}=-c_n/2\) recovers exactly \(4S+12\lambda_n\).

Borwein, Borwein and Garvan, *Some cubic modular identities of Ramanujan*, Lemma 2.1(ii),(iii) and Proposition 2.2(i),(ii), [author PDF](https://qseries.org/fgarvan/papers/cubic.pdf). Martin and Ono, *Eta-quotients and elliptic curves*, Theorems 2 and 3, [author PDF](https://uva.theopenscholar.com/files/ken-ono/files/026_8.pdf).

## Lemma EisTail. A uniform tail for arbitrary targets

Fix \(n\equiv2,4\bmod6\), with no hypothesis on its split prime factors. For \(p\ge5\), \(p\equiv2\bmod3\), let \(b_p^o\) count ordered odd-component norm-\(n\) shell points with \(p\mid N(z)N(w)\). Write \(\mu_p=2/(p(p-1))\). If \(p\nmid n\), apply the binary classification's Lemma on the operator \(T_m\), after projection off constants and the determinant character on \(\Gamma_p\backslash{\rm PGL}_2(\mathbf F_p)\). Its hypotheses hold at \(m=n\) and \(m=n/4\), whenever integral. The scalar 2 reduces projectively to one. Thus \(D_n=T_n-T_{n/4}\) counts the odd shell, has constant eigenvalue \(S\), and has complementary norm at most \(2\tau(n)\sqrt n\) by the triangle inequality. The torus-normalizer orbit average gives
\[
 b^o_p/R_o\le \mu_p+O_\epsilon(n^{-1/2+\epsilon}).
\]
Here the divisor estimate \(\tau(n)\ll_\epsilon n^\epsilon\) applies to every target. If \(p\mid n\) and spoils a component, the other is also zero modulo \(pR\), because the norm form is anisotropic modulo \(p\). The spoiled count vanishes if \(p^2\nmid n\). Otherwise the ratio is exactly \(S(n/p^2)/S(n)=\sigma(p^{a-2})/\sigma(p^a)<p^{-2}\) for \(p^a\Vert n\). This is at most \(\mu_p\).

Let \(M(n)=\max_{1\le u\le n}\prod_{q^a\Vert u,\ q\equiv1(3)}(a+1)\), the largest split-prime divisor product through \(n\). For every odd inert \(p\), the lattice bound in the binary classification gives \(b_p^o/R_o\le 4M(n)/p^2\), \(M(n)\le\max_{t\le n}\tau(t)\ll_\epsilon n^\epsilon\). Summing up to \(P=n^{1/4}\), using at most \(P\) integers, and summing \(p>P\) by \(\sum_{k>P}k^{-2}\ll P^{-1}\), gives for any fixed \(Y\ge 2\), and any \(0<\epsilon<1/4\),
\[
 \limsup_{n\to\infty}\sum_{p>Y} b_p^o(n)/R_o(n)
 \le \sum_{\substack{p>Y\\p\equiv2(3)}}\mu_p .
 \tag{2}
\]
Both sums here and below run over odd inert primes. The two discarded errors are \(O_\epsilon(n^{-1/4+\epsilon})\), uniformly in \(Y\le P\). In particular the right side tends to zero with \(Y\), and the whole odd-prime loss is at most \(3/20+o(1)\), using the binary classification's elementary wheel bound.

The pure case satisfies \(|G_n|/R_o\ge17/20-o(1)\). For mixed targets, (1) and the preceding odd subtraction show that the fraction surviving the condition at 3 and all odd inert primes, in either orientation, is at least \(1-1/3-3/20-o(1)=31/60-o(1)\). Swapping components is a bijection between the two disjoint orientations. Hence \(|G_n|/R_o\ge31/120-o(1)\) for the fixed orientation. In both cases \(|G_n|\gg R_o\ge18n\) as \(n\to\infty\).

## Fixed sublattices and harmonics

Give \(R^2\subset\mathbf C^2\simeq\mathbf R^4\) its Euclidean quadratic form \(|z|^2+|w|^2\). For fixed \(Y\), let \(G_{n,Y}\) impose all of the preceding conditions except for primes \(p>Y\). On even targets, oddness is the complement of the sublattice \(2R\oplus2R\). Avoiding an inert \(p\) in both coordinates is the intersection of the complements of \(pR\oplus R\) and \(R\oplus pR\). In the mixed case the orientation is \(w\in(1-\omega)R\) and \(w\notin3R\). Because \(n\equiv1\bmod3\), this already forces \(z\notin(1-\omega)R\). All intersections of the finitely many sublattices so obtained are full rank. Inclusion-exclusion therefore writes \(1_{G_{n,Y}}\), on the shell, as a fixed finite integral combination of ordinary lattice indicators.

For any one of these fixed lattices \(L\) and any homogeneous harmonic polynomial \(H\) of degree \(d>0\), the spherical theta series
\[
 \Theta_{L,H}(\tau)=\sum_{v\in L} H(v)q^{|v|^2}
\]
is a holomorphic cusp form of weight \(k=2+d\) on a congruence subgroup. For odd \(d\) the series vanishes under \(v\mapsto-v\). For even \(d\) it has trivial character on \(\Gamma_0(N_L)\), for a fixed level. Its integral even Gram matrix uses the bilinear form \(2\langle u,v\rangle\), and has determinant \(9[R^2:L]^2\), a square. The even-rank quadratic character is the Kronecker symbol of \((-1)^{4/2}\det L\). Equivalently each lattice is a direct sum of two scaled hexagonal lattices whose weight-one quadratic characters multiply to one. Walling's Theorems 1 and 2 give respectively the modular transformation and cuspidality. The half-norm convention \(e^{\pi i\tau(v,v)}\) is \(q^{|v|^2}\) for this Gram matrix.

At fixed level, \([q^n]\Theta_{L,H}\ll_{L,H,\epsilon}n^{(k-1)/2+\epsilon}\) for **all** integers \(n\). Decompose into finitely many oldforms \(f(t\tau)\) from normalized eigennewforms. Their coefficients are zero if \(t\nmid n\) and \(a_f(n/t)\) otherwise. At good primes use Deligne. For a prime dividing the newform level once, the trivial-character newform eigenvalue has magnitude \(p^{k/2-1}\). If the square of the prime divides the level, it vanishes. The Euler denominator at each bad prime has degree at most one. The bound therefore extends to all prime powers by multiplicativity, with \(\tau(n)\) loss. These normalizations for arbitrary levels are recorded together in Barrett, Burkhardt, DeWitt, Dorward and Miller, equations (2.3), (2.7), (2.9)–(2.11). Thus
\[
 \sum_{\substack{v\in L\\ |v|^2=n}}H(v/\sqrt n)
  =n^{-d/2}[q^n]\Theta_{L,H}
  =O_{L,H,\epsilon}(n^{1/2+\epsilon}).
\]
The fixed finite inclusion-exclusion proves the same bound (with constant depending on \(Y,H,\epsilon\)) on \(G_{n,Y}\).

L. Walling, *The Eichler Commutation Relation for theta series with spherical harmonics*, Theorems 1 and 2, [article PDF](https://bibliotekanauki.pl/articles/1391780.pdf). Her correction concerns the corollary to Theorem 3 and Section 4. [Walling, correction](https://doi.org/10.4064/aa-66-4-391). A. O. L. Atkin and J. Lehner, *Hecke operators on \(\Gamma_0(m)\)*, Theorem 3, [DOI](https://doi.org/10.1007/BF01359701). For the arbitrary-level newform and Fourier normalizations [Barrett et al., §2](https://web.williams.edu/Mathematics/sjmiller/public_html/math/papers/GeneralN_21.pdf). The unramified estimate uses [Deligne, modular forms, Theorems 4.9 and 5.1](https://www.numdam.org/item/SB_1968-1969__11__139_0.pdf) and [Deligne, Weil I, Theorem 1.6](https://www.numdam.org/item/PMIHES_1974__43__273_0.pdf).

## Proof of Theorem EisAngle

The set difference \(G_{n,Y}\setminus G_n\) is bounded in size by \(\sum_{p>Y}b^o_p(n)\). For fixed \(H\), divide the fixed-lattice estimate by \(R_o\ge18n\), and bound the effect of the set difference by \(\|H\|_{\infty,S^3}\) times (2). First let \(n\to\infty\) within its admissible class, then \(Y\to\infty\). This yields \(\sum_{v\in G_n}H(v/\sqrt n)=o(R_o)=o(|G_n|)\). The constant function 1 has average one, and positive-degree spherical harmonics have spherical mean zero. Finite sums of their restrictions and constants are dense in \(C(S^3)\) by polynomial approximation and the harmonic decomposition of homogeneous polynomials. The probability measures giving equal mass to each signed, ordered witness in \(G_n\), rescaled by \(\sqrt n\), thus converge weakly to normalized surface measure on \(S^3\).

Writing \(z=e^{i\theta}\cos\phi,\ w=e^{i\psi}\sin\phi\) on \(S^3\), surface area is proportional to \(\sin\phi\cos\phi\,d\phi\,d\theta\,d\psi\). Hence \(|z|^2\) is uniform on \([0,1]\). For any fixed interval \(I\subset[0,1]\) of positive length (any endpoint convention),
\[
 \#\{(z,w)\in G_n:N(z)/n\in I\}
       = (|I|+o_I(1))|G_n|.
\]
Its two boundary slices have zero surface measure. Thus each sufficiently large target in either residue class has a decomposition with both actual summand values between \(n/3\) and \(2n/3\). In the mixed problem the two balanced values are \(a\) and \(3b\). Any fixed positive interval for \(a/n\) is covered. Each ordered pair of positive norm values accounts for at most \(36M(n)^2\ll_\epsilon n^\epsilon\) witnesses (use the divisor estimate with \(\epsilon/2\)). The interval therefore contains \(\gg_{I,\epsilon}n^{1-\epsilon}\) distinct ordered pairs. The limiting proportions count each value pair with its number of coordinate witnesses.

The threshold in Theorem EisAngle depends on the interval. Its proof gives no numerical onset. For ordinary balanced sums the recorded finite certificate finds precisely two exceptions, \(8\) and \(44\), below \(10^8\). Determining every exception above that range remains open. The angular result has internal proof support. External review and priority remain open. [Finite balanced calculation](../../unpublished/eisenstein-restricted-prime-sums/sections/balanced-summands.tex).

The Gaussian sieve strategy has earlier primary precedents: [Lucia's answer,
3 September 2013](https://mathoverflow.net/a/141096), works with local forbidden-prime
fractions and a small/large-prime split; [Dietmann and Elsholtz's answer,
5 September 2013](https://mathoverflow.net/a/141306), describes a modular-form
route with a possible finite completion. These are method precedents, not
identified proofs of the Eisenstein statements here. The 2026 Gaussian release
remains the immediate starting point of the manuscript's exposition.
