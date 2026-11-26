# Pairwise Sato–Tate bounds for small Frobenius cofactors

Let \(A/K\) be an abelian variety of dimension \(g\ge3\) over a totally
real number field. Suppose that its full geometric rational endomorphism
algebra is a totally real field \(F\) of degree \(g\), with all its
endomorphisms defined over \(K\). At a good prime \(v\) put \(q_v=Nv\) and
write \(X^2-a_vX+q_v\) for its relative Frobenius polynomial. Label the
real embeddings of \(F\) and write \(x_i(v)=\sigma_i(a_v)/(2\sqrt{q_v})\).
The normalized absolute trace and discriminant norms are
\[
 T_v=\prod_{i=1}^g|x_i(v)|=\frac{|N_F(a_v)|}{2^g q_v^{g/2}},\qquad
 D_v=\prod_{i=1}^g(1-x_i(v)^2)=\frac{(-1)^gN_F(a_v^2-4q_v)}{4^g q_v^g}.
\]

**Theorem.** Let \(M/K\) be a fixed finite Galois extension with \(M\)
totally real or CM, and let \(\mathcal C\) be a conjugacy class of
its Galois group. For any distinct \(i,j\), the rectangle discrepancy
of \((x_i(v),x_j(v))\) at good absolute
degree-one primes \(v\) of norm at most \(X\) in the class \(\mathcal C\) is
\(O_{A,K,M}(\delta_X)\) from product Sato–Tate measure on \([-1,1]^2\), where
\[
                 \delta_X=\frac{\sqrt{\log\log X}}{(\log X)^{1/4}}.
\]
Write \(P_{\mathcal C}(X)\) for proportion in this class and range. For
\(0<u\le1\), uniformly as \(X\to\infty\),
\[
\begin{split}
 P_{\mathcal C}(X)(T_v\le u)&\ll_{A,K,M}u^{2/g}+\delta_X\log(2/\delta_X),\\
 P_{\mathcal C}(X)(D_v\le u)&\ll_{A,K,M}u^{3/g}+\delta_X\log(2/\delta_X).
\end{split}
\]
Both exponents of \(u\) are sharp for inference from pairwise
independence and the Sato–Tate marginals alone.
In dimension three, for \(0<u,w<1/8\), there is also the joint bound
\[
 P_{\mathcal C}(X)(T_v\le u,\,D_v\le w)
       \ll_{A,K,M}u w^{3/4}+u^{1/2}w^{3/2}+\delta_X.
\]

For the split one-failure real-cubic problem, take \(M\) to be the
compositum with the normal closure of \(F\) and \(\mathcal C=\{1\}\).
At absolute degree-one primes this is the completely split coefficient
condition. With \(b_v=N_F(a_v)/p\) and \(d_v=-N_F(a_v^2-4p)/p\), one has
\(T_v=|b_v|/(8\sqrt p)\) and \(D_v=d_v/(64p^2)\). Among split primes,
the one-failure primes with \(|b_v|\le\epsilon(p)\sqrt p\) or
\(d_v\le\epsilon(p)p^2\) have relative density zero for every function
\(\epsilon(p)\to0\). This concerns a fixed characteristic-zero abelian
variety. The possible failures with cofactors of full Weil order remain.
[Frobenius cofactors and the split residual case](CM_NEWTON_SPLIT_COFACTOR_NORM_AND_TRACE_COUNT.md).

## Pairwise distribution in a prescribed Frobenius class

The \(g\) rank-two compatible systems from the \(F\)-action have weights
\(0,1\), cyclotomic determinant, and connected derived monodromy
\(\mathrm{SL}_2^g\) over an algebraic closure of the auxiliary coefficient
field. The Hodge–Tate weights follow from the real multiplication Hodge
decomposition, and the monodromy follows from Lombardo's common-determinant
image theorem at sufficiently large split auxiliary primes. Each
system remains irreducible after every finite extension. A scalar twist
acts trivially on the projective representation, so independence of the
derived factors rules out twist equivalence between different labels
after any finite extension. Each system is regular and totally odd,
the latter because its rank-two alternating pairing has the cyclotomic
similitude and every real conjugation has determinant minus one.
[Lombardo, Theorem 1.4](https://arxiv.org/pdf/1411.1703).

Apply simultaneous potential automorphy over \(M\) if it is totally
real, or over its maximal real subfield \(M^+\) otherwise, taking \(K\)
for the ground field in the Galois clause. This gives one totally real
finite Galois \(L/K\) over which the labeled systems are cuspidal Hilbert
modular. For any intermediate field \(k\) with \(L/k\) Galois and solvable,
the same holds over \(k\): descend
one prime-cyclic step at a time. At each step the representation over
the larger field is invariant by strong multiplicity one because its
Galois system is a restriction from the smaller field. A cuspidal
cyclic descent is still regular algebraic by the archimedean local
parameters. The associated Galois system agrees with the required one
on the larger field, by base change and Chebotarev. Their restrictions
there are irreducible, so the two systems on the smaller field differ
by a character of the cyclic group. Twist the descent by that character
via class field theory. This gives the original system on each fixed
field. The unitary normalizations have trivial central character by the
cyclotomic determinant, and are nondihedral and pairwise
twist-inequivalent by monodromy.
[Barnet-Lamb–Gee–Geraghty–Taylor, Theorem 5.4.1](https://www.ma.imperial.ac.uk/~tsg/Index_files/potentialautomorphy.pdf),
[Arthur–Clozel, Chapter 3, Theorem 4.2(d)](https://www.claymath.org/library/cw/arthur/pdf/30.pdf),
[Barnet-Lamb–Gee–Geraghty, §7.1](https://arxiv.org/pdf/0912.1054).

Let \(U_m(\cos\theta)=\sin((m+1)\theta)/\sin\theta\), and put
\(\widetilde L=LM\), \(G=\operatorname{Gal}(\widetilde L/K)\). Inflate the
indicator of \(\mathcal C\) to \(G\) and expand it in irreducible characters.
Brauer induction writes each character as a combination of
\(\operatorname{Ind}_H^G\eta\) with \(H\) solvable and \(\eta\) one-dimensional.
When \(M\) is CM, \(\widetilde L\) is CM with maximal real subfield \(L\).
Its complex conjugation \(c\) is central in \(G\). For a given irreducible
character, project its Brauer decomposition onto the central eigenspace
with the same sign for \(c\). This keeps or kills a term with \(c\in H\),
according to \(\eta(c)\). If \(c\notin H\), then the character induced to
\(H\langle c\rangle=H\times\langle c\rangle\) splits into the two linear
extensions of \(\eta\). The projection keeps the one with the required
sign. All resulting solvable subgroups contain \(c\). Thus in both cases
the intermediate field \(k=\widetilde L^H\) is totally real and \(L/k\)
is solvable and Galois. The character \(\eta\) gives a finite-order ray
class character of \(k\) by reciprocity, including the possible sign
at real places.

For each induced term, the fixed cosets of Frobenius at \(v\) give its
degree-one relative primes \(w\) in \(k\), with \(\eta\) evaluated on
their Frobenius classes in \(\widetilde L/k\). Hence the prime sum of
\(U_m(x_i(v))U_n(x_j(v))\), weighted by
the induced character, is the prime sum of the corresponding symmetric
powers over \(L^H\) twisted by the finite-order ray character \(\eta\),
up to higher residue degree and finitely many omitted primes. The total
omission for fixed fields is \(O_{A,K,L}((m+1)(n+1)\sqrt X)\), since
\(|U_m U_n|\le(m+1)(n+1)\) and the prime below an omitted higher-degree
prime has norm at most \(\sqrt X\). This includes the restriction to
absolute degree one in \(K\).
[Brauer, character induction theorem](https://doi.org/10.2307/1969183),
[Wong, Theorem 1.1 and §3](https://www-math.nsysu.edu.tw/~pjwong/stuff/CST.pdf).

Newton and Thorne supply the cuspidal symmetric powers over the totally
real fields \(L^H\). Thorner's Proposition 6.2 gives a prime number
estimate for their Rankin–Selberg convolution, including the ray
character twists in the preceding paragraph and the case when one
symmetric power is zero. With \(\mathcal M=\max(m,n)\ge1\) its error
for each fixed field, representations and character is bounded by
\[
 X\exp\!\left(-\frac{c\log X}{\mathcal M^4\log(C\mathcal M)
                         +\mathcal M\sqrt{\log X}}\right)
                (\mathcal M^4\log(C\mathcal M X))^6
\]
up to fixed constants, taking the worse of the allowed conductor
exponents. In particular, for
\(1\le\max(m,n)\le c_0(\log X)^{1/4}/\sqrt{\log\log X}\), with \(c_0>0\)
sufficiently small, the bound is \(O_{A,K,M}(X/(\log X)^{10})\). Indeed
the first term of the denominator is \(O(c_0^4\log X/\log\log X)\),
and the second is smaller for large \(X\). All the fields and characters
used in the finite Brauer expansion were fixed before choosing \(c_0\).
[Newton–Thorne, Theorem A](https://arxiv.org/pdf/2212.03595),
[Thorner, Proposition 6.2](https://arxiv.org/html/2404.06482v2#S6.SS1).

The zero Chebyshev mode is the count of the selected primes. Finite
Chebotarev gives its asymptotic \((|\mathcal C|/[M:K])\operatorname{Li}(X)\).
The two-dimensional interval approximations in the proof of Thorner's
Theorem 2.1 use Chebyshev modes of individual degree at most \(N\), with
bounded coefficients and a constant-term error \(O(1/N)\), uniformly in
the endpoints of a closed rectangle. Apply them with
\(N=\lfloor c_0/\delta_X\rfloor\). The \(O(N^2)\) nonzero modes above
contribute \(o(\delta_X\operatorname{Li}(X))\). Dividing by the selected
count proves the discrepancy. The rate is uniform in rectangles for
fixed \(A,K,M\); potential automorphy is used for existence of the finite
fields entering its constants.
[Thorner, §6.2](https://arxiv.org/html/2404.06482v2#S6.SS2).

## From two coordinates to small cofactor norms

Under Sato–Tate measure \(2\pi^{-1}\sqrt{1-x^2}\,dx\) on \([-1,1]\),
both \(Z=|x|\) and \(Z=1-x^2\) take values in \([0,1]\), with
\(\Pr(Z\le z)\asymp z^\lambda\) for small \(z>0\), where \(\lambda=1\)
and \(\lambda=3/2\), respectively. An interval condition on \(Z\) is a
union of at most two interval conditions on \(x\), so the same
discrepancy bound holds for pairs of such variables.

Sort the values of the \(g\) transformed coordinates and call the two
smallest \(z_1,z_2\). A product at most \(u\) entails
\(z_1z_2^{g-1}\le u\). Put \(u'=\max(u,\delta_X^{g/(2\lambda)})\) and
\(r=(u')^{1/g}\). For \(z_2\le r\), the pair lies in \([0,r]^2\),
whose product measure is \(O(r^{2\lambda})\). On the shell
\(2^kr<z_2\le2^{k+1}r\), it lies in a rectangle with upper first
coordinate \(u'/(2^kr)^{g-1}\). That rectangle has product measure
\(O((u')^{2\lambda/g}2^{-k(g-2)\lambda})\). Sum over the finitely many
coordinate pairs and over the shells. The measures sum to
\(O((u')^{2\lambda/g})\) because \(g>2\), and there are
\(O_g(\log(2/\delta_X))\) rectangles. Truncation of \(u\) and the
discrepancy errors give both claimed bounds. For a varying threshold
\(\epsilon(p)\to0\), first discard finitely many primes to bound it by
any fixed \(u>0\), then let \(u\to0\).

For the joint bound in dimension three, classify a coordinate as
small when \(|x_i|\le1/\sqrt2\). If all are small then \(D_v\ge1/8\);
if none is small then \(T_v>1/8\). With exactly one small coordinate
it must have \(|x_i|\ll u\); one of the two other coordinates must
have \(1-x_j^2\ll w^{1/2}\). The product Sato–Tate measure in this
pair is \(O(uw^{3/4})\). With exactly two small coordinates one has
\(|x_i|\ll u^{1/2}\) on one of them, and the remaining coordinate
has \(1-x_j^2\ll w\). The pair measure is \(O(u^{1/2}w^{3/2})\).
In both cases only a bounded number of pair rectangles is used.

For sharpness let \(U,V\) be independent uniform variables on
\(\mathbb R/\mathbb Z\). For \(i=0,\ldots,g-1\), put
\(W_i=(U+iV)\bmod1\). Any pair of distinct coordinates is uniform
on the two-dimensional torus because the associated integer matrix
has nonzero determinant and is therefore surjective on that torus.
Take \(Z_i\) to be the Sato–Tate quantile of the indicated transform
at \(W_i\). This gives pairwise independence with the required
marginals. The original real traces can be recovered by taking
independent uniform signs times \(Z_i\) or times \(\sqrt{1-Z_i}\),
respectively. On \(0<U,V<c r^\lambda\) for a sufficiently small constant
\(c>0\), all \(W_i\) lie below the transform's distribution function
at \(r\), so all \(Z_i\le r\). This event has probability
\(c^2r^{2\lambda}\). Setting \(r=u^{1/g}\) proves the sharpness assertion.
The example is a probability law, with no assertion about a
compatible Frobenius system realizing it.

Full joint weak equidistribution together with every fixed auxiliary
Frobenius distribution is already compatible with varying one-failure
finite-field examples. The present result applies to fixed systems and
discards the boundary-cofactor ranges. It does not address the diagonal
residue-characteristic divisibility on the remaining interior.
[Local simultaneous realizability](CM_NEWTON_THREE_BALANCED_GLOBAL.md#simultaneous-archimedean-and-auxiliary-frobenius-data).
