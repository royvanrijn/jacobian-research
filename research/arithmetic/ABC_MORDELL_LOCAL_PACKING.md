# Local packing for the \((2,3,6)\) residual radical

In an anatomic box for primitive \(a+b=c\le X\), assign the square, cube and sixth-power roles to an ordered permutation \((S,T,F)\) of the letters. Write \(\alpha_{J,i}\) for the mass of the exact valuation class \(i\), \(L_J\) for the total mass on \(J\), and \(D\) for the count exponent. From the first six classes define
\[
\begin{aligned}
 E_J&=\alpha_{J,2}+\alpha_{J,4}+\alpha_{J,6},
 & C_J&=\alpha_{J,3}+\alpha_{J,6},\\
 \rho_+&=L_S+L_T+L_F-E_S-C_T-\alpha_{F,6},
 & \sigma_+&=L_S+L_F-E_S-E_F,\\
 w_-&=\sum_{J\in\{S,T,F\}}
       \sum_{\substack{1\le i\le6\\ e_J\nmid i}}
                    \frac{i}{e_J}\alpha_{J,i},
 & (e_S,e_T,e_F)&=(2,3,6).
\end{aligned}
\]
For rational \(k=a/b\ge1\) and rational \(w_0>0\) with
\(3^a(3w_0)^b\ge1\), the domain \(w_-\ge w_0\) has the affine bound
\[
 D\le \rho_++(f(k)-1)\sigma_+,\qquad
 f(k)=
 \begin{cases}
 (5k+13)/18 & 1\le k\le4,\\
 11/6+(1597/5000)(k-4) & k\ge4 .
 \end{cases}
\]
The two branches are \(w_-\le w_0\), and \(w_-\ge w_0\) with the displayed bound. In particular \(k=1,\ w_0=1/9\) gives \(D\le\rho_+\). The test uses retained masses and has no short-height hypothesis.

**Local separation.** Write the signed canonical equation \(Ax^2+By^3=Cz^6\), with square-free, cube-free and sixth-free coefficient magnitudes. Put \(\Delta=A^3 B^2 C\), \(U=-AB y/z^2\) and \(V=A^2B x/z^3\). On \(E_\Delta:V^2=U^3+\Delta\), the ratio \(t=U^3/\Delta=-By^3/(Cz^6)\) determines the primitive triple and \(t+1=Ax^2/(Cz^6)\). Scaling to \(E_1\) by a sixth root of \(\Delta\) gives \(u^3=t,\ v^2=t+1\) and \(\widehat h=h(t)/6+O(1)\). Counted points have height at most \(\log X/6+O(1)\). Rational torsion is uniformly bounded by Northcott on \(E_1\), since the scaled points have degree at most six and bounded height. [Silverman, VIII.9](https://www.math.ens.psl.eu/~obenoist/refs/Silverman.pdf#page=264)

Fix the coefficients and the retained valuation of each residual prime used by \(w_-\). The residual radical \(R=\operatorname{rad}|\Delta|\) determines at most \(8^{\omega(R)}\) coefficient magnitudes, and the retained classes or omitted tail give at most \(7^{\omega(R)}\) assignments. At \(p\nmid6\) with locked valuation \(i\), primitivity gives \(v_p(t+1)=i\) on \(S\), \(v_p(t)=i\) on \(T\), and \(v_p(t)=-i\) on \(F\). On \(E_1\) these are residue discs about \((u,0)\) with \(u^3=-1\), about \((0,\pm1)\), and about infinity, respectively. Their local parameters \(v,u,u/v\) have valuations \(i/2,i/3,i/6\). Partition by the disc at each locked prime. Translation by the corresponding torsion point and the integral formal group law put a difference \(P-Q\) in the infinity disc with parameter valuation at least \(i/e_J\). The function \(t\) has pole order six there. Thus the rational \(t(P-Q)\) has denominator divisible by \(p^{6i/e_J}\). Omitting \(2,3\) loses at most \(O(1)\) in the log weight of the retained classes. There are \(X^{o(1)}\) parts, and distinct Mordell–Weil vectors in each part satisfy \(\widehat h(P-Q)\ge (w_--o(1))\log X\). [Silverman, IV.1–3 and VII.2](https://www.math.ens.psl.eu/~obenoist/refs/Silverman.pdf#page=134)

**Packing and averaging.** Scale the Mordell–Weil vectors into the unit ball and lift \(v\) to \((v,\sqrt{1-\|v\|^2})\) on the unit sphere. Distances cannot decrease. Their pairwise dot products are at most \(\eta+o(1)\), where \(\eta=\max(0,1-3w_0)\). A random cap of threshold \(\sqrt{\eta+\delta}\) contains each point with probability \(\gg_\delta (r+1)^{-1}(1-\eta-\delta)^{(r+1)/2}\), where \(r\) is the rank. Projection of the cap to its axis complement gives pairwise dot products at most \(-\delta/2\) for large \(X\). Squaring the sum of the projected vectors bounds the cap occupancy by \(1+2/\delta\).

The [CM rank bound](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#the-3-class-group-in-the-rank-bound) gives \(3^r\ll C^{1+\omega(R)}h_3(d)^2\), with \(d=\operatorname{sf}(\Delta)\) the signed squarefree kernel. In particular \(r=O(\log X)\). Choosing fixed \(\delta\) small in terms of \(\varepsilon,w_0\) absorbs the excess in the cap exponent into \(X^\varepsilon\). Since \(3^{-k}\le\min(1,3w_0)=1-\eta\), the resulting count per coefficient slice is \(O_\varepsilon(X^\varepsilon h_3(d)^k)\). This also handles equality at the affine threshold.

The signed first moments and fourth moments of \(h_3\), interpolated for \(1\le k\le4\) and combined with the pointwise bound for \(k>4\), give \(\sum_{|d|\le Z}h_3(d)^k\ll_{k,\varepsilon}Z^{f(k)+\varepsilon}\) over signed squarefree kernels. The split kernel contributes one. Since \(|d|\mid R\), there are at most \(Y/|d|\) possible radicals \(R\le Y\) per kernel. Dyadic summation gives \(O_\varepsilon(X^\varepsilon Y Z^{f(k)-1})\). Apply \(Y=X^{\rho_++o(1)}\) and \(Z=X^{\sigma_++o(1)}\). The parity kernel comes from odd valuations on \(S,F\), and discarding tail removals gives the stated upper masses. [Heath-Brown–Pierce, Theorem 1.2 and Corollary 1.4](https://arxiv.org/pdf/1409.3177) · [Chan–Koymans, Theorem 1.1](https://arxiv.org/pdf/2505.00611)
