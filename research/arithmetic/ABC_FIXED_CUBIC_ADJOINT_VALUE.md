# A fixed cubic adjoint value for abc over the rationals

Let \(a+b=c\) be a primitive positive triple, \(R=\operatorname{rad}(abc)\),
\(A=\log R\), and \(h=\log c\). Use the Frey curve
\[
E:y^2=x(x-a)(x+b).
\]
Fix either primitive cubic Dirichlet character \(\chi\) modulo \(7\) once
for the entire family. For \((a,b,c)\ne(1,1,2)\), let \(f_E\) be the
rational weight-two newform of \(E\), and put
\(\sigma_E=\operatorname{Ad}\pi_{f_E}\otimes\chi\) in unitary normalization.
Write \(P(E)=\frac{i}{2}\int_{E(\mathbb C)}\omega_E\wedge\bar\omega_E\)
for the area of the period lattice of a minimal Néron differential.
The quantity
\[
X_E=\frac{G(\bar\chi)^2 L(1,\sigma_E)}{\pi P(E)}
\]
belongs to the fixed field \(\mathbb Q(\zeta_3)\), is nonzero, and satisfies
\[
\log|X_E|=h+O(\log(h+2)).
\tag{CA1}
\]
The completed value satisfies
\[
\log|\Lambda(1,\sigma_E)|=A+O(\log(A+2)).
\tag{CA2}
\]
Both estimates are uniform over the triples, including when \(7\mid abc\).
Here \(G(\bar\chi)=\sum_{r\bmod 7}\bar\chi(r)e^{2\pi i r/7}\).

## Conductor and cuspidality

At every odd prime dividing \(R\), the displayed model has \(c_4\) a unit
and positive discriminant valuation, so \(E\) is multiplicative. At all
other odd primes it has good reduction. The two invariants used here are
\(c_4=16(a^2+ab+b^2)\) and \(\Delta=16a^2b^2c^2\).
The good and multiplicative alternatives each give trivial inertia in
the Weil part of the adjoint Weil–Deligne representation. The
multiplicative case has one monodromy block of size three, so its
adjoint conductor exponent is two. Twisting by \(\chi\) preserves
these exponents away from \(2,7\). At \(7\) the twist is tame with
nontrivial scalar inertia of order three, so its exponent is exactly
three, whether the reduction is good or multiplicative. Inertia on the
dual is the inverse scalar. This proves that \(\sigma_E\) is non-self-dual
even though its central character is trivial.
[Weil–Deligne conductor formula for special blocks and twists](https://arxiv.org/abs/2309.02782).

There is an odd multiplicative prime unless \(R=2\).
Pairwise coprimality makes \((1,1,2)\) the only positive triple with
\(R=2\). At a multiplicative prime, \(\pi_{f_E}\) is an unramified twist
of Steinberg and has nonzero monodromy. Any nontrivial self-twist of
\(\pi_{f_E}\) would be quadratic, by the central character, and would
make it an induction of a character of a quadratic field. Such an
induction has zero monodromy at every finite place. Thus the
Gelbart–Jacquet theorem gives a cuspidal adjoint lift, and a character
twist preserves cuspidality.
[Gelbart–Jacquet, Theorem 9.3 and Remark 9.9](https://www.numdam.org/articles/10.24033/asens.1355/).

At \(2\), the adjoint conductor exponent is absolutely bounded. Indeed
the curve acquires semistable reduction over \(\mathbb Q_2(E[5])\), of
degree at most \(|\operatorname{GL}_2(\mathbb F_5)|\). On inertia over
this field the \(5\)-adic matrices are congruent to the identity
modulo \(5\). Their eigenvalues are roots of unity by quasi-unipotence,
and the congruence forces each eigenvalue \(\lambda\) to satisfy
\(v_5(\lambda-1)\ge1\). A nontrivial root of unity has
\(v_5(\lambda-1)\le1/4\). The inertia action is therefore unipotent,
and the semistable reduction criterion applies. The finite Weil part
of the adjoint then has trivial inertia over this extension. There
are finitely many extensions of \(\mathbb Q_2\) of the bounded degree,
so their upper ramification breaks bound the Swan conductor uniformly.
The tame and monodromy terms are bounded by the rank. It follows that
the arithmetic conductor is
\[
C_E=7^3\,2^{O(1)}\!
\left(\prod_{\substack{p\mid R\\p\ne2,7}}p\right)^2,
\qquad \tfrac12\log C_E=A+O(1).
\tag{CA3}
\]

## A two-sided bound at the edge

The finite Euler parameters of the unitary cuspidal twist have modulus
at most one. At potentially good places this follows from pure
weight zero for the adjoint after restriction to a field of good
reduction. At potentially multiplicative places the adjoint has a
single special block with unitary twisting character, so any surviving
Euler parameter on inertia invariants in the monodromy kernel has
modulus at most one. The archimedean parameters are fixed at each real
place. For the weight-two adjoint and an even character the factor is
\(\Gamma_{\mathbb R}(s+1)\Gamma_{\mathbb C}(s+1)\).

Let \(Q\) be the analytic conductor at \(s=1\). Brumley's theorem,
applied with the trivial \(\mathrm{GL}_1\) representation as the
self-dual factor and the given non-self-dual cuspidal representation as
the other factor, yields a zero-free region of width
\(\gg1/\log(3+Q)\) for \(|\operatorname{Im}s|\le1\). The central
characters are unitary and trivial on the positive real scalars.
[Brumley, Theorem A.1 and Remark A.2 in Humphries](https://arxiv.org/abs/1703.05450).

For a sufficiently small fixed \(c_0>0\), put
\(\delta=c_0/\log(3+Q)\). The disk centered at \(1+\delta\) of radius
\(3\delta\) is zero-free. The Euler product bounds the canonical
\(\log L(1+\delta)\) by \(O(\log(2/\delta))\), using
\(\zeta(1+\delta)^3\) as a majorant for the Euler product. On the line
\(\operatorname{Re}s=-\delta\) the functional equation and the same
Euler bound give
\[
|L(s)|\ll \delta^{-O(1)}Q^{1/2+\delta}(2+|\operatorname{Im}s|)^{O(1)}.
\]
The powers and constants are uniform by Stirling's formula and the
fixed archimedean types. Apply the three-lines bound in the strip
\([-\delta,1+\delta]\), dividing by a sufficiently large fixed power
of \(s+2\) to account for polynomial growth. On the part of the disk
inside this strip the exponent of \(Q\) in the interpolated bound is
\(O(\delta)\). The rest of the disk uses the Euler product. Thus the
real part of the analytic logarithm is \(O(\log(2/\delta))\) from above
on the disk. Borel–Carathéodory, with the logarithm fixed by its Euler
value at the center, bounds the absolute value of the logarithm at
\(1\) by the same order. In particular
\[
\big|\log|L(1,\sigma_E)|\big|\ll\log\log(3+Q)\ll\log(A+2).
\tag{CA4}
\]
Combining (CA3) and (CA4), with the fixed archimedean factor at \(1\),
proves (CA2).

## A value in one fixed quadratic field

The even character makes \(s=2\) critical for the weight-two symmetric
square in classical normalization. Its Euler roots have weight two,
so this value is \(L(1,\sigma_E)\) in unitary normalization. Apply
Schmidt's rationality theorem at classical \(s=k=2\), with
\(\delta=1\) and trivial nebentypus. The resulting algebraic number is
\(-G(\bar\chi)^2 L(1,\sigma_E)/(4\pi^3\langle f_E,f_E\rangle_1)\).
It lies in the field generated by the coefficients of \(f_E\) and the
values of \(\chi\), with Galois equivariance in both. Here
\(\langle\ ,\,\rangle_1\) is the Petersson integral for
\(\Gamma_1(N_E)\), and the value is for the primitive finite
\(L\)-function.
[Schmidt's rationality theorem in the normalization of Loeffler–Zerbes, Proposition 2.2.1 and Theorem 2.2.3](https://wrap.warwick.ac.uk/id/eprint/83189/7/WRAP-Iwasawa-theory-symmetric-form-Loeffler-2017.pdf).

For a modular parametrization \(\phi:X_0(N_E)\to E\) over
\(\mathbb Q\), \(\phi^*\omega_E=m_\phi\,2\pi i f_E(z)\,dz\) has
\(m_\phi\in\mathbb Q^\times\) by the rational differential structures.
Integrating the two-form gives
\(\deg(\phi)P(E)=4\pi^2m_\phi^2\langle f_E,f_E\rangle_0\),
where the index from \(\Gamma_0\) to \(\Gamma_1\) is a rational
integer. Therefore \(\pi^2\langle f_E,f_E\rangle_1/P(E)\) is rational
and \(X_E\in\mathbb Q(\zeta_3)^\times\). This uses the period
normalization even for a nonoptimal curve in an isogeny class.

The displayed Weierstrass equation is minimal at odd primes
and its differential differs from the minimal Néron differential by
an absolute factor at most two. At \(2\), the valuation of \(c_4\) is
exactly four, since \(a,b\) are coprime. If the scale change to a
minimal equation is \(u\), its exponent at \(2\) is at most one
because \(c_4\) changes by \(u^{-4}\). The area computed with
\(dx/(2y)\) is
\[
\frac{4}{c}
 K\big(\sqrt{a/c}\big)\,K\big(\sqrt{b/c}\big),
\quad
K(k)=\int_0^{\pi/2}(1-k^2\sin^2\theta)^{-1/2}\,d\theta.
\]
The real and imaginary cycles are given by doubling the integrals
on \([a,\infty)\) and \([0,a]\), respectively. The two factors are
bounded below by absolute constants. One factor has squared modulus
at most \(1/2\) and is bounded above. The other is
\(O(1+\log c)\) by \(1-k^2\ge1/c\) and comparison of the integrand
near \(\pi/2\) with \(1/\sqrt{1-k^2+(\pi/2-\theta)^2}\).
Consequently \(\log P(E)=-h+O(\log(h+2))\).
Since \(A\le\log(abc)<3h\), equation (CA4) proves (CA1).

## An equivalent valuation bound

Conjugation sends \(X_E\) to the value for \(\bar\chi\). The signed
rational valuations of its norm satisfy the exact identity
\[
\tfrac12\sum_p v_p\big(\mathrm N_{\mathbb Q(\zeta_3)/\mathbb Q}X_E\big)
     \log p=\log|X_E|.
\tag{CA5}
\]
Full abc is equivalent to proving, for each \(\varepsilon>0\),
\[
\log|X_E|\le (1+\varepsilon)\log R+O_\varepsilon(1)
\tag{CA6}
\]
uniformly over the triples. In either direction, (CA1) and
\(\log(h+2)=o(h)\) allow the error to be absorbed by choosing a
slightly smaller \(\varepsilon\) in the assumed inequality. The
single exceptional triple only changes the constant.

The missing estimate is an upper bound for the signed sum in (CA5)
with the coefficient in (CA6), or a sequence for which that sum
exceeds the coefficient by a fixed positive amount with unbounded
excess. This formulation uses a fixed rational base, a fixed cubic
character and a fixed coefficient field. Rationality and the
two-sided analytic estimate supply the reduction. The size of the
signed sum of finite valuations remains open as a function of the
radical.

As a comparison, the Euler-system bound of Loeffler and Zerbes
applies to the untwisted adjoint at a good ordinary prime of big
residual image. Theorem A controls a cyclotomic characteristic ideal
by the adjoint \(p\)-adic \(L\)-function, with an undetermined power of
that prime. The target here requires an aggregate upper bound across
rational primes for the fixed cubic twist in (CA5).
[An Euler system for the adjoint of a modular form, Theorem A](https://arxiv.org/abs/2312.04665).

The fixed cubic adjoint has a uniformly bounded product of
global torsion orders, and the product of its local Tamagawa
factors for \(p\ne\ell\) has logarithm
\(O(h/\log(h+2))\). These factors can be divided out of the
norm valuation problem without changing the leading coefficient.
The integral comparison at \(p=\ell\) and the global Selmer or
congruence estimate remain open.
[Local Tamagawa factors and global torsion](ABC_CUBIC_ADJOINT_TAMAGAWA_FACTORS.md).
