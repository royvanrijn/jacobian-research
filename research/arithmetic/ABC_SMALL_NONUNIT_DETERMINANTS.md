# Coprime pairs at small nonzero distance and small-determinant abc triples

## A fixed-parameter pair estimate

Fix \(\theta\ge3\) and \(\epsilon>0\). Let \(W\) be sufficiently large in terms of \(\theta,\epsilon\). For every integer \(h\) with
\[
0<|h|\le W^{1/(4000\theta^2)},
\]
the number of positive coprime pairs \(m_1,m_2\le W^\theta\) with \(m_2-m_1=h\) and
\[
m_1=d^2u,\qquad m_2=e^2v,\qquad u,v\ \text{square-free},\qquad W\le uv<2W
\]
is
\[
\ll_{\theta,\epsilon} W^{1-1/(400\theta^2)+\epsilon}.
\tag{P}
\]
The implied constant is uniform in \(h\). The conclusion also holds with a fixed constant multiplying the height cap or the shift cap. The implicit constant may then depend on those constants. The proof uses the source lemmas of Diao [1], specifically the balanced approximate determinant and square-sieve argument in §5 and the quadratic-field argument in §§6–7. We give the changes required for arbitrary \(h\).

Put \(P=UV\) for the dyadic sizes \(u\asymp U,\ v\asymp V\), with dyadic sizes \(d\asymp D,\ e\asymp E\). Then \(P\asymp W\). There are \(O_\theta((1+\log W)^4)\) boxes. For a surviving box, \(e^2v\) and \(d^2u\) are comparable, since their product is at least \(W\) and their difference is \(O(W^{1/(4000\theta^2)})\). Thus
\[
Y_0=DE\sqrt P\asymp D^2U\asymp E^2V.
\tag{1}
\]
In particular \(Y_0\ll W^\theta\). Coprimality gives \(\gcd(d,e)=\gcd(u,v)=1\) and \(\gcd(h,uvde)=1\).

### Balanced kernel sizes

Assume \(U,V\ge W^{1/(200\theta^2)}\) up to fixed factors. If \(DE\le P^{1-1/(200\theta^2)}\), freeze \(d,e\). Their coprimality makes the solutions of \(e^2v-d^2u=h\) a single arithmetic progression in \(u,v\), with steps \(e^2,d^2\). By (1) this gives
\[
\ll DE(1+P/Y_0)=DE+\sqrt P
\]
points in the box, which suffices. Otherwise use
\[
M=\lceil P^{1-1/(100\theta^2)}\rceil,\qquad B=\lceil\theta^2\rceil .
\]
The proof of [1, Proposition 5.2] uses ratio intervals of lengths \(D/(EM)\) and \(V/(UM)\) in \(d/e\) and \(v/u\). The normalized ratios satisfy \(t=s^2+O(|h|/Y_0)\), since the exact error is \(h/(e^2u)\) before normalization. Thus each interval for either ratio pairs with \(O(1)\) intervals for the other, because \(M|h|\ll Y_0\). In fact \(Y_0>P^{3/2-1/(200\theta^2)}\).

The logarithmic determinant coordinates in [1, Lemma 4.1] now have sizes \(O(1/M)\) and \(O(|h|/Y_0)\). Apply [1, Lemma 4.2] with its error parameter \(T=|h|\). In the proof of [1, Proposition 5.2], this increases the logarithmic bound for the integer determinant by at most \(B\log|h|\) and otherwise gives the same expression (5.7)–(5.8). In detail, put \(p=\log M/\log P,\ q=\log Y_0/\log P\); for sufficiently large \(W\), \(p\ge9/10,\ 1-p\le1/(50\theta^2),\ q\le7\theta/6\). The original margin
\[
\Phi=\frac{B+1}{2}-\frac{q^2}{4p}-(1-p)B(B+1)\ge \frac{\theta^2}{12}
\]
dominates \(B\log|h|/\log P\le B/(4000\theta^2)+o(1)\). All required integer determinants therefore vanish, giving the same auxiliary forms of bidegree at most \((1,B)\) with the coefficient bounds of [1, Lemma 4.3].

The proof of [1, Proposition 5.1] now applies with the following changes. A factor independent of \(d,e\) contributes only \(O(B)\) primitive directions for \((u,v)\) per interval pair. For each such direction the binary quadratic equation with right-hand side \(h\) contributes \(O_\theta(\tau(|h|))\) points by [1, Lemma 3.4], since \(uv\) is square-free, comparable to \(W\), and the coordinates are polynomially bounded in \(W\). This costs \(W^\epsilon\) after shrinking the auxiliary epsilon. For every remaining factor, [1, Lemma 4.4] gives a nonzero square value of an odd-degree binary form; that lemma has arbitrary nonzero right-hand side. The additional \(\log |h|\) in [1, (4.16)] preserves the polylogarithmic square-sieve hypothesis. The interval sizes, lattice argument in [1, Lemma 3.3], and square sieve [1, Corollary 3.2] give the bound (5.1)–(5.5) of that proof in the form
\[
\ll_{\theta,\epsilon}
 W^\epsilon(P/M)^{1/4}
 \big((PM)^{1/2}+U+V\big).
\]
Here the term \(O_\theta(\tau(|h|)M)\) is absorbed because \(M\le P\). Since \(U+V\ll P^{1-1/(200\theta^2)}\) and \(M\asymp P^{1-1/(100\theta^2)}\), the resulting power is at most \(P^{1-1/(400\theta^2)+\epsilon}\).

### Unbalanced kernel sizes and a nonsplit quadratic field

Exchange the two letters, and change the sign of \(h\) accordingly, so the small kernel is \(v\ll W^{1/(200\theta^2)}\). For fixed \(v,h\) set
\[
n=ev,\qquad w=uv,\qquad A=-hv,\qquad d^2w=n^2+A.
\tag{2}
\]
We have \(\gcd(n,d)=1\) and \(0<|A|\le 2W^{1/(100\theta^2)}\) for large \(W\). For fixed \(v\), the original tuple is recovered from \(d,w,n\). The dyadic \(d\)-scale satisfies \(D\ll W^{\theta/2}\le W^{\theta-1/2}\) at large \(W\).

For \(D\le W^{1-1/(32\theta^2)}\), the proof of [1, Lemma 7.1] applies with arbitrary \(A\), because \(\gcd(n,d)=1\) forces \(\gcd(A,d)=1\) in any contributing term. The congruence \(n^2\equiv-A\pmod{d^2}\) then has at most \(4\cdot2^{\omega(d)}\) roots by Hensel lifting and the Chinese remainder theorem. For a fixed \(d\), the interval of positive \(n\) has length \(O(d\sqrt W)\) by comparing the squares \(d^2w-A\). The result is \(O((W^{1/2}+D)\log(2D))\).

For the larger range write \(A=A_0 q_0^2\), with \(A_0\) signed square-free and \(q_0>0\). Suppose first that \(A_0\ne-1\). We use the following variant of [1, Proposition 7.2], restricted to \(\gcd(n,d)=1\):
\[
 N_A^{\rm prim}(D,W)
 \ll_{\theta,\epsilon}
 |A_0|^{7/8} W^{1-1/(32\theta^2)+\epsilon}.
\tag{3}
\]
The source argument, with the following explicit changes, proves (3).

1. In [1, §6.1] use \(K=\mathbb Q(\sqrt{-A_0})\) with absolute discriminant \(\Delta\asymp|A_0|\). Put \(d_2=\gcd(d,2)\), \(d'=d/d_2\), and \(\zeta=(n+q_0\sqrt{-A_0})/d_2\). This element is integral. If \(d\) is even, (2) and primitivity force \(n,q_0\) odd and \(-A_0\equiv1\pmod4\), which is the needed integrality condition. The coefficient \(z_2\) of the nonconstant integral basis element is nonzero and has size at most \(2q_0\). For fixed \(A\) it takes at most four values, from the basis choice and \(d_2\). The elements \(\zeta,\bar\zeta\) have relatively prime ideals at primes above \(d'\), because their sum \(2n/d_2\) is coprime to \(d'\). Thus the ideal \(\mathfrak d\) with norm \(d'\) and \(\mathfrak d^2\mid(\zeta)\) used in [1, Lemma 6.1] is obtained by the same prime-ideal factorization. Minkowski gives \(\mu\in\mathfrak d\setminus\{0\}\) with both embedding sizes \(O((d')^{1/2}\Delta^{1/4})\). With \(k=|\operatorname{Nm}\mu|/d'\) and \(\nu=k^2\zeta/\mu^2\), we have
\[
 k\ll\Delta^{1/2},\qquad
 \mu^2\nu=k^2\zeta,\qquad
 |\operatorname{Nm}\mu|=d'k,\qquad
 \operatorname{Nm}\nu=wk^2 .
\]
These are the identities (6.2). With \(\mu=r\alpha,\ \nu=\ell\beta\) and primitive integral coordinate pairs for \(\alpha,\beta\), the coefficient comparison is (6.3) with \(\xi=k^2z_2/(r^2\ell)\ne0\). The proof of (6.10) gives \(\ell\mid k\) from the squarefreeness of \(w\). Fixing \(k,r,\ell,z_2\) and the two coordinate pairs determines \(\zeta\), hence \(n,d,w\), up to at most the two choices for \(d_2\).
2. In [1, §6.2], the bound on the archimedean error \(t_2-2t_1\) comes from \(\zeta/\bar\zeta\) and is still \(O(|A|^{1/2}/n)\). The embedding rescalings, slope parametrization, and Jacobian bounds of [1, Lemmas 6.2–6.5] are unchanged. Where a lattice determinant occurs, use \(\Delta\), which is the actual discriminant of \(K\). In the size bounds one may use \(\Delta\ll |A|\). The nondegenerate fixed-direction estimate [1, Lemma 6.6] only requires the identity \(\operatorname{Nm}\beta=w(k/\ell)^2\) with positive square-free \(w>1\), and still gives \(O_\theta(\tau(|\xi|))\).
3. For [1, Lemma 7.3], \(M=\lceil W^{1-1/(8\theta^2)}\rceil\) and \(B=\lceil\theta^2\rceil\) are unchanged. The error penalty is \(B\log|A|/2\) from [1, Lemma 4.2]. The other discriminant terms in the integer determinant have upper bound
\[
 (B+1)(2B+1)\log k-\tfrac12(B+1)^2\log\Delta
 \ \le\ \tfrac12B(B+1)\log\Delta+O(B^2),
\]
by \(k\ll\Delta^{1/2}\). Since \(\Delta\ll|A|\), the original upper bound in Lemma 7.3 and its negative margin remain valid. The same parameter \(M\) obeys \(M|A|^{1/2}\ll D\) for the stated large-\(D\) range. In [1, Lemma 7.4], the larger \(z_2\) still gives \(|\xi|\le k^2(2q_0)\ll |A|^{3/2}\), sufficient for the unchanged logarithmic bounds for the auxiliary coefficients and squared-value form.
4. The local counts in [1, Lemmas 7.6–7.7] therefore apply before summation in \(k,r,\ell,z_2\). For each fixed \(A\), apply the divisor inequality \(\tau(k^2z_2)^j\le \tau(k^2)^j\tau(|z_2|)^j\). The second factor costs \(W^\epsilon\) and the first is summed by [1, Lemma 7.5] with \(z_2=1\). For \(S_1\) and the sparse interval pairs, this gives \(\Delta^{1/2}M W^\epsilon\) up to \(\theta\)-dependent constants. For \(S_2\), the local estimate (7.6) is unchanged, namely \(\ll \log^4(2|A|W)R^{3/4}M^{1/4}\Delta^{-3/8}\) with \(R=W(k/\ell)^2\). The sum of \((k/\ell)^{3/2}\) over \(k\ll\Delta^{1/2}\) and \(r^2\ell\mid k^2z_2\) is \(\ll \Delta^{5/4}W^\epsilon\), using \(\sum \ell^{-3/2}<\infty\) and the divisor moment above. Thus \(S_2\ll_{\theta,\epsilon}\Delta^{7/8}W^{3/4+\epsilon}M^{1/4}\). The sparse and fixed-direction contributions are smaller than this bound. Since \(\Delta\asymp |A_0|\), substitution of \(M\) proves (3).

Sum the elementary bound and (3) over \(O(W^{1/(200\theta^2)})\) possible \(v\)'s at fixed \(h\). The exponent in (3) after summing is at most
\[
 1-\frac1{32\theta^2}+\frac7{800\theta^2}+\frac1{200\theta^2}
 =1-\frac7{400\theta^2}.
\]
The small-\(D\) bound gives at least as much saving. This leaves room for the requested exponent of (P) and polylogarithmic boxes.

### The split quadratic case

If \(A_0=-1\), then \(hv>0\) is a square. Since \(h\) and the square-free \(v\) are coprime, \(v=1\) and \(h=\kappa^2\) with \(\kappa>0\). In particular the second term is \(e^2\). Now
\[
 d^2u=e^2-\kappa^2=(e-\kappa)(e+\kappa).
\]
Here \(e>\kappa\), and \(\gcd(e,\kappa)=1\) by primitivity. Thus \(g=\gcd(e-\kappa,e+\kappa)\) is either \(1\) or \(2\). The new positive coprime pair
\[
 (m'_1,m'_2)=\left((e-\kappa)/g,(e+\kappa)/g\right)
\]
has difference \(h'=2\kappa/g\), the same square-free kernel product \(u=w\), and height at most \(2\sqrt{\max(m_1,m_2)}\). The kernel product remains \(w\) because \(m'_1m'_2=d^2u/g^2\) and \(g\mid d\). At fixed \(h\) and \(g\), \(e\) and the original pair are determined by the new pair. The new shift satisfies \(|h'|\le\max(|h|,2)\); this remains within the modified shift cap with an absolute constant. Orientation changes and the two possible gcds create only a bounded factor per step. Repeat if the new coprime pair is again in the split case. After \(K=K(\theta)\) steps the height would be at most \(4W^{\theta/2^K}<\sqrt W\), which contradicts the unchanged kernel product at least \(W\), since the product of the two complete values would be smaller than \(W\). Thus every branch reaches a previously counted region in at most \(K\) steps. These bounds are uniform in all encountered shifts and apply at the same \(\theta\) with a fixed constant in the height cap. The bounded depth and bounded multiplicities complete (P).

## Consequence for small determinants

Let
\[
 o(m)=\prod_{\substack{p\mid m\\ v_p(m)\ {\rm odd},\ v_p(m)\ge3}}p.
\]
There is an absolute \(\eta_0>0\) such that for every fixed \(0<\lambda\le1\) the number of ordered primitive positive triples \(a+b=c\le X\) satisfying
\[
 \operatorname{rad}(abc)<c^\lambda,\qquad
 o(a)o(b)o(c)\le c^{\eta_0\lambda^3},
\]
and admitting positive multipliers \(u,v\le c^{\eta_0\lambda^3}\) with
\[
 0<|va-ub|\le c^{\eta_0\lambda^3}
\]
is \(O_{\lambda,\epsilon}(X^{\lambda/3-\eta_0\lambda^3+\epsilon})\). The determinant can have any prime factorization.

Put \(\rho=\lambda/3,\ \delta=10^{-9}\rho^3,\ H=Y^\delta\) in a height shell \(Y/2<c\le Y\). A letter with radical at most \(Y^{\rho-6\delta}\) has \(O_\epsilon(Y^{\rho-6\delta+\epsilon})\) possible complete values by the usual radical box count. Selecting the two small multipliers at most \(2H\) and a nonzero shift at most \(H\) counts any opposite value with at most \(O(H^3)\) choices for the data. The relations \((u+v)a-uc=va-ub\) and \((u+v)b-vc=ub-va\) ensure that this argument applies regardless of the specified letter. The shell contribution here is \(O_\epsilon(Y^{\rho-3\delta+\epsilon})\).

Otherwise the two smallest radicals are both at most \(Y^{\rho+3\delta}\), since their product with the third is at most \(Y^{3\rho}\) and the smallest exceeds \(Y^{\rho-6\delta}\). Pick those letters, scale by the corresponding multipliers \(j_1,j_2\le2H\) in the difference identities, and divide the two values by their gcd \(g\). Since \(g\) divides the nonzero determinant, \(g\le H\). The resulting positive coprime values \(m_1,m_2\) differ by a nonzero integer of size at most \(H\). Their common dyadic height \(M\) lies in \([Y/(8H),2HY]\) up to absolute factors. Indeed the selected pair contains a letter at least \(Y/4\), and the two scaled values are within \(H\) of each other. Each resulting pair has \(O_\epsilon(HY^\epsilon)\) original triples, by choosing \(g\) and the divisors \(j_i\) of \(gm_i\). The radicals of the normalized values are at most \(2Y^{\rho+4\delta}\), and
\[
 o(m_1)o(m_2)\le 4Y^{5\delta}.
\]
For the latter bound, passing from a letter \(n\) to \(jn/g\) can introduce new primes in \(o\) only from \(jg\), and the original odd-overlap product is at most \(H\).

Write the normalized values as \(m_i=d_i x_i^2\) with square-free kernels \(d_i\) and localize \(d_1d_2\asymp W=Y^s\). For each fixed nonzero shift, their fixed \((d_1,d_2)\) fibres are \(Y^\epsilon\), by the uniform quadratic norm count [1, Lemma 3.4] when \(d_1d_2>1\) and elementary difference of squares when \(d_1d_2=1\). Thus for \(s\le\rho-3\delta\), choose the kernel product, factor it, and sum over the at most \(2H\) shifts and \(O_\epsilon(HY^\epsilon)\) preimages. This costs \(O_\epsilon(Y^{s+2\delta+\epsilon})\).

For large \(W\) use the support-corrected elementary cover. The identity
\[
 \operatorname{rad}(x_i)=\operatorname{rad}(m_i)o(m_i)/d_i
\]
bounds the possible \(x_1,x_2\), after dyadic localization of \(d_i,o(m_i)\), by
\[
 \ll_\epsilon Y^\epsilon
 \frac{(Y^{\rho+4\delta})^2Y^{5\delta}}{W}.
\]
For each fixed \(x_1,x_2\) and shift, coefficient solutions of \(d_2x_2^2-d_1x_1^2=h\) form one progression with steps \((x_2^2,x_1^2)\), since \(\gcd(x_1,x_2)=1\). It contributes \(O(1+W/M)\) points in the two coefficient intervals. This factor is bounded: \(W\ll Y^{2\rho+8\delta}\) and \(M\gg Y^{1-\delta}\). Including shifts and preimages now gives \(O_\epsilon(Y^{2\rho+15\delta-s+\epsilon})\), which suffices for \(s\ge\rho+16\delta\).

In the remaining range \(\rho/2\le s\le\rho+16\delta\), use the single fixed parameter \(\theta=8/\rho\) in (P). The height cap holds since \(s\theta\ge4\). The shift cap holds as
\[
 \frac{s}{4000\theta^2}\ge\frac{\rho^3}{512000}>\delta.
\]
The saving in the \(Y\)-exponent at a fixed shift is \(s/(400\theta^2)\ge\rho^3/51200\). This covers the distance and preimage costs \(2\delta\) and the possible excess \(s-\rho\le16\delta\), with at least \(\delta\) to spare. Summing the dyadic boxes and height shells proves the claim with \(\eta_0=10^{-9}/27\).

## Sources and scope

[1] Yijie Diao, *On the average of class numbers of indefinite binary quadratic forms*, arXiv:2610.28001v1, 20 October 2026. [Theorem 1.3 and §§2–8](https://arxiv.org/html/2610.28001v1). The altered proof above uses Lemmas 3.3–3.4, 4.2–4.4, Corollary 3.2, and the arguments of §§5–7. The widening to arbitrary shifts uses these arguments, rather than an application of Diao's theorem about the norm-one equation to a variable norm.

[2] Chunping Zeng, *Power-saving bounds for balanced radical triples via small Pell discriminants*. [Hexagon 2611.00611v1, especially §§4 and 6](https://hexagonmath.org/pdf/2611.00611v1). The support-corrected large-kernel cover above uses the existing project identity for \(\operatorname{rad}(x_i)\), which pays for varying square-free kernels. Zeng's fixed-support lemma alone leaves the count over varying supports in Lemma 4.1 unproved.

