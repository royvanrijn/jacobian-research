# An unconditional Frobenius cofactor bound for real multiplication

Let \(B/K\) be an abelian \(g\)-fold whose geometric rational
endomorphism algebra is a totally real field \(F\) of degree \(g\),
with its endomorphisms defined over \(K\). At a good degree-one prime
\(v\) of characteristic \(p\) unramified in \(F\) and prime to the
endomorphism-order index, write \(X^2-a_vX+p\) for the
[relative Frobenius polynomial](CM_NEWTON_FROBENIUS_RAMIFICATION_AND_DENOMINATORS.md),
\(D_v=a_v^2-4p\), and \(r_v=g-f_p(B_v)\) for the loss of \(p\)-rank.
For odd \(p\) the discriminant has local valuation zero at a coefficient
prime where the trace is a unit, and valuation one where the trace is
in the maximal ideal. In the first case the slopes are \((0,1)\);
in the second they are \((1/2,1/2)\). Each pair contributes with
its coefficient residue degree to the absolute characteristic
polynomial. Hence \(v_p(N_{F/\mathbb Q}(D_v))=r_v\).
The Weil bound makes \(D_v\) totally negative, since equality
\(a_v^2=4p\) is impossible at an unramified coefficient prime over
\(p\). Let \(c_v\) be the positive squarefree part of the positive
integer \((-1)^g N_{F/\mathbb Q}(D_v)/p^{r_v}\).

**Theorem.** For \(0<\delta<\theta<1/(6g)\), the number of these
good degree-one primes with \(p\le X\) and \(c_v\le(\log p)^\delta\) is
\[
 O_{B,K,\theta,\delta}\!\left(
     \frac{X\log\log X}{(\log X)^{1+\theta-\delta}}\right).
\]
In particular, for a fixed real-cubic threefold the locus of split
one-failed-pair primes with \(c_v\le(\log p)^\delta\) has relative
density zero among all degree-one split primes for every
\(\delta<1/18\). Here \(c_v\) is precisely the cofactor used in the
[real-cubic square sieve](CM_NEWTON_FROBENIUS_DISCRIMINANT_SIEVE.md).

## Uniform Chebotarev averaging

Let \(L\) be the normal closure of \(F\). For a sufficiently large
rational prime \(\ell\) split in \(F\), the projective auxiliary group is
\[
 P_\ell=\left\{(g_1,\ldots,g_g)\in\operatorname{PGL}_2(\mathbb F_\ell)^g:
       \det(g_1)=\cdots=\det(g_g)
       \text{ in }\mathbb F_\ell^\times/\mathbb F_\ell^{\times2}\right\}.
\]
For distinct such primes \(\ell,\ell'\in[Y,2Y]\), let \(M\) be the
compositum over \(K\) of \(KL\) and the two projective division fields.
Lombardo's common-determinant theorem and Serre's almost-independence
theorem apply in this dimension. The product argument and local
discriminant estimate in the [real-cubic sieve](CM_NEWTON_FROBENIUS_DISCRIMINANT_SIEVE.md#effective-auxiliary-distributions-on-the-identity-class)
apply with \(|P_\ell|\le\ell^{3g}\), giving
\[
 G:=\operatorname{Gal}(M/K)
    =\operatorname{Gal}(KL/K)\times P_\ell\times P_{\ell'},\qquad
 n_M\ll_{B,K}Y^{6g},\quad
 \log|D_M|\ll_{B,K}Y^{6g}(\log Y)^2.
\]
For the product assertion, the kernel of the projection to
\(\operatorname{Gal}(KL/K)\) surjects onto the unrestricted product of
the two projective groups, and the projection onto
\(\operatorname{Gal}(KL/K)\) is surjective.
[Lombardo, Theorem 1.4](https://arxiv.org/pdf/1411.1703) ·
[Serre, Theorem 1 and §3.1](https://www.college-de-france.fr/media/jean-pierre-serre/UPL8570957518086251784_independance.pdf).

Take \(Y=(\log X)^\theta\). For a class function \(W\) on \(G\) with
\(|W|\le1\), Thorner–Zaman's estimate on \(X/2<Nv\le X\) gives
\[
 \sum_{X/2<Nv\le X} W(\operatorname{Frob}_v)
  = I(X)\,\mathbb E_G W
    - I_\beta(X)\,\mathbb E_G(W\eta)
    + O_{B,K,\theta,A}\bigl(I(X)(\log X)^{-A}\bigr).
\]
Here the sum is over unramified degree-one primes,
\(I(X)=\operatorname{Li}(X)-\operatorname{Li}(X/2)\), and
\(I_\beta(X)=\operatorname{Li}(X^\beta)-\operatorname{Li}((X/2)^\beta)\)
when there is a sufficiently close real zero. In that case \(\eta\)
is a single quadratic or trivial character of \(G\), and
\(I_\beta(X)\ll I(X)\). In the absence of that zero set \(I_\beta(X)=0\).

To see the common character explicitly, set
\(U=\log(|D_M|n_M^{n_M})\) and apply Thorner–Zaman's Theorem 1.4
to each class using \(H=\langle g\rangle\) on a representative \(g\).
Conductor–discriminant for \(M/M^H\) gives
\(\log(|D_{M^H}|\mathcal Q n_{M^H}^{n_{M^H}})\le U\).
Since \(U\log\log X=o(\log X)\) and
\(n_M(\log\log X)^2=o(\log X)\), both error terms save any fixed
power of \(\log X\), uniformly in the class. Any term from a real zero
with \(1-\beta\ge1/(8U)\) has the same saving. A closer zero of
\(\zeta_M\) is real and simple by Stark's zero-free region. Stark's
Theorem 3 assigns it an intermediate field \(T/K\), either trivial
or quadratic, with \(\zeta_{M^H}(\beta)=0\) exactly when
\(T\subseteq M^H\). If \(T\subseteq M^H\), the Hecke character in
Theorem 1.4 is trivial on \(H\). Otherwise the quadratic character
of \(TM^H/M^H\) contributes the zero, by factoring \(\zeta_{TM^H}\).
The character on \(H\) is thus in both cases the restriction of the
character \(\eta\) of \(T/K\). Its zero is in the specified
Hecke exceptional range since \(1-\beta<1/(8U)\). This proves the
averaging formula by summing the errors weighted by class size.
Removing primes of higher absolute degree costs \(O_K(X^{1/2})\).
[Thorner–Zaman, Theorem 1.4](https://msp.org/ant/2019/13-5/ant-v13-n5-p02-s.pdf) ·
[Stark, Lemma 3 and Theorem 3](https://doi.org/10.1007/BF01405166).

## A squareclass test that includes zero discriminants

Fix an integer \(r\) between zero and \(g\). For \(\chi_\ell\) the
quadratic character extended by zero, use the projective test
\(H_\ell=\chi_\ell((-1)^gd^{-r}\prod_i(t_i^2-4d))\), as in
the [real-cubic square sieve](CM_NEWTON_FROBENIUS_DISCRIMINANT_SIEVE.md#a-square-sieve-for-the-cofactor),
where \(d\) is the common determinant of linear representatives and
\(t_i\) are their traces. Define
\[
   U_\ell(c)=\chi_\ell(c)H_\ell+\mathbf 1_{H_\ell=0}.
\]
For \(\ell\ne p\), if \((-1)^gN(D_v)/p^r\) has squareclass \(c\), the test
is one on \(\operatorname{Frob}_v\): when \(H_\ell\ne0\) the
squareclass can be read modulo \(\ell\), and when it is zero the
indicator supplies one. Here the normalized norm is \(\ell\)-integral, so
if \(H_\ell\ne0\) its valuation is zero and the rational square
relating the two representatives is an \(\ell\)-adic unit.

Every quadratic character of \(P_\ell\) factors through the common
determinant squareclass because the kernel is a product of \(g\)
perfect groups \(\operatorname{PSL}_2(\mathbb F_\ell)\) for \(\ell\ge5\).
After conditioning on the common determinant in the full linear group,
the discriminant calculation gives
\(\mathbb E(H_\ell\mid d)=O(\ell^{-2g})\), and the probability of
\(H_\ell=0\) is \(O(\ell^{-1})\). The latter follows from the at
most two repeated-root traces in each factor, each represented by
\(\ell^2\) matrices out of \(\ell(\ell^2-1)\).
The full linear group pushes forward uniformly to the projective
group. For every quadratic or trivial character \(\eta_\ell\) this gives
\[
          \big|\mathbb E(U_\ell(c)\eta_\ell)\big|\ll\ell^{-1},
\]
uniformly in \(c\) and \(r\). Thus for \(W=U_\ell(c)U_{\ell'}(c)\)
both \(\mathbb E_G W\) and \(\mathbb E_G(W\eta)\) are \(O(Y^{-2})\)
by the product decomposition of \(G\). Multiplication of \(W\) by
the indicator of the identity in \(\operatorname{Gal}(KL/K)\) preserves
the estimate if one counts only split primes.

Write \(\mathcal P\) for the split auxiliary primes in \([Y,2Y]\),
so \(J=|\mathcal P|\asymp_F Y/\log Y\). Among degree-one primes
in \(X/2<p\le X\), the unnormalized second moment of
\(S_v(c)=\sum_{\ell\in\mathcal P}U_\ell(c)(\operatorname{Frob}_v)\)
is at most
\[
   O_{B,K,\theta,A}\bigl(I(X)(J+J^2Y^{-2}+J^2(\log X)^{-A})\bigr)
       = O_{B,K,\theta}\bigl(I(X)J\bigr)
\]
upon choosing \(A>1\). Diagonal terms use \(|U_\ell|\le1\).
The prime characteristic in this interval exceeds \(2Y\); all the
division fields in the sum are unramified there after discarding
the fixed bad primes. For a prescribed squareclass, \(S_v(c)=J\),
so its count in the interval is \(O(I(X)/J)\). There are at most
\((\log X)^\delta\) candidate values of \(c\) and \(g+1\) possible
values of \(r_v\). Sum over these choices and then the dyadic
intervals to prove the theorem.

The unconditional argument controls cofactors growing as a small power
of \(\log p\). The Weil bound gives \(c_v\le4^g p^{g-r_v}\), which
is \(64p^2\) for a real-cubic threefold with one failed pair. The density
of the entire split one-failed-pair locus remains open.
