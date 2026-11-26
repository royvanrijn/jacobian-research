# Even-start collinear ternary windows under large common translation

Let \(a\) be a strictly positive integer exponent vector. Let \(b_1,b_2,b_3\) be distinct collinear integer exponent offsets of the same total degree. For each fixed even \(m\ge2\), there is a threshold such that the ordinary support \(Ta+b_1,Ta+b_2,Ta+b_3\) has no three consecutive zero factorial moments at \(m,m+1,m+2\) for any integer \(T\) above that threshold. The coefficients are arbitrary complex numbers, not all zero. The threshold may depend on the direction, offsets, and prescribed start.

The [quadratic Gaussian reduction](SFC_HOMOGENEOUS_TERNARY_SHIFTED_GATES.md#collinear-offsets) reduces the assertion to excluding a common solution of the three equations
\[
 G_j(A,B):=\mathbb E[(A+B X_j+X_j^2)^j]=0,\qquad
                  X_j\sim\mathcal N(0,1/j).
\]
We prove exclusion for any three positive orders consisting of an odd order and two even orders with distinct two-adic valuations. If a common zero exists at an odd start \(m\), the same argument and a centered expansion impose a single pair of valuations at each place over two. The final comparison refines this pair by explicit quarter-valuations on the odd and even orders.

## Even orders and distinct valuations

The three \(G_j\) are rational polynomials. A complex solution would give an algebraic solution by the weak Nullstellensatz over \(\overline{\mathbb Q}\). Embed the field of its coefficients into \(\overline{\mathbb Q}_2\) and put \(a_0=v_2(A)\), \(b_0=v_2(B)\), allowing the value \(+\infty\) for a zero parameter. Write \(e=v_2(j)\). The even Gaussian moments are \(\mathbb E X_j^{2r}=(2r-1)!!/j^r\), with odd numerator and \(\mathbb E X_j^0=1\). In \(j^jG_j\), the term taking \(p\) factors \(A\), \(q\) factors \(BX_j\), and \(j-p-q\) factors \(X_j^2\) vanishes for \(q\) odd. Otherwise its valuation is
\[
 p(a_0+e)+q(b_0+e/2)+v_2\binom{j}{p,q,j-p-q}.
\]
When \(j\) is even, all three pure terms occur, with multinomial coefficient one. If a single line among \(a_0+e,\ b_0+e/2,\ 0\) is strictly lowest, its pure term is the unique term of least valuation. Thus an even zero requires a tie on the lower envelope.

Two ties at distinct positive \(e\) require all three lines, with the slopes decreasing from \(1\) to \(1/2\) to \(0\) as \(e\) increases. In particular \(a_0<\min(b_0,0)\) at \(e=0\). To see the slope conclusion, a line active at the first tie with higher slope than a second active line lies strictly above that second line at every later point. A second tie therefore requires a third line of still smaller slope. Similarly the higher slope at the first tie is uniquely active below that tie. For the odd order, the pure \(A^j\) term exists and is the unique term of least valuation, a contradiction. At an even start the two even orders are two apart, so their valuations are distinct. This proves the stated translation theorem.

## Odd starts and the remaining valuation pair

Let \(m\) be odd and put \(e=v_2(m+1)>0\). At an odd order, a unique least weight among \(a_0,b_0,0\) already excludes a zero when it is the weight of \(A\) or \(X_j^2\). If \(b_0<\min(a_0,0)\), the smallest possible weights come from \(j-1\) factors \(BX_j\) and a single remaining factor \(A\) or \(X_j^2\). Both multinomial coefficients are \(j\), a unit. All terms using fewer \(B\)'s have strictly larger valuation than the lower of those two. Thus \(a_0=0\) is necessary in this case. The same comparison works at \(j=1\), where only the two indicated terms occur. The only possibilities for an odd zero are therefore \(a_0=0>b_0\), \(a_0=b_0<0\), or \(a_0,b_0\ge0\). The last possibility makes the line \(0\) uniquely lowest at the even order \(m+1\). In either of the first two possibilities \(a_0+e>b_0+e/2\). The even-order tie consequently requires \(b_0=-e/2\).

We now exclude \(a_0=0\) for every odd start. Put \(u=B^{-2}\), so \(v_2(u)=e\). At an odd order \(j\), divide \(G_j\) by the nonzero \(jB^{j-1}\mathbb E X_j^{j-1}\), obtaining
\[
 F_j(A,u)= A+1+\sum_{\ell=1}^{(j-1)/2}
       u^\ell c_{j,\ell}(A),\qquad
 c_{j,\ell}(A)=
 \frac{\binom{j}{2\ell+1}}{j}\,
 \frac{\mathbb E[X_j^{j-2\ell-1}(A+X_j^2)^{2\ell+1}]}
      {\mathbb E X_j^{j-1}} .
\]
The coefficient \(A+1\) uses the exact Gaussian ratio \(\mathbb E X_j^{j+1}/\mathbb E X_j^{j-1}=1\). All \(c_{j,\ell}\) belong to \(\mathbb Z_2[A]\).

In fact \(c_{j,\ell}(-1)\in4\mathbb Z_2\) for every \(\ell\ge1\). Under the Gaussian tilt by \(X_j^{j-2\ell-1}\), the \(r\)-th moment of \(Y=X_j^2\) is
\[
 \mu_r=\prod_{h=0}^{r-1}\frac{j-2\ell+2h}{j}
 \equiv 1+\frac2j\left(\binom r2-\ell r\right)\pmod4.
\]
The centered moment \(\mathbb E[(Y-1)^{2\ell+1}]\) under this tilt is the finite difference of order \(2\ell+1\) of \(\mu_r\) at zero. Its reduction modulo four vanishes, since the displayed reduction is quadratic in \(r\).

Set \(A=-1+x\). The equation is \(x+u(4 f(u)+xg(x,u))=0\) with \(f,g\) integral over \(\mathbb Z_2\). Formal coefficient recursion gives a solution \(x_j(u)\in4u\mathbb Z_2[[u]]\), convergent for \(v_2(u)>0\). It is the unique solution with \(x\) integral, since the difference of the left sides at two integral arguments is their difference times a unit. For \(j\ge3\), direct calculation at \(\ell=1\) gives \(c_{j,1}(-1)=-2(j-1)/(3j)\). Indeed the first three tilted moments are \((j-2)/j,(j-2)/j,(j^2-4)/j^2\), so the centered third moment is \(-4/j^2\). The ratio of the original Gaussian moments is \(j/(j-2)\). At \(j=1\) the same formula gives zero, consistent with the empty sum.

The coefficient of \(u\) in \(x_{j+2}-x_j\) is \(4/(3j(j+2))\), of valuation two. Every higher coefficient is divisible by four. At \(v_2(u)>0\) the linear term is therefore strictly minimal and \(x_j(u)\ne x_{j+2}(u)\). This excludes the unit-\(A\) branch simultaneously at both odd orders.

In summary, an odd-start solution can only have
\[
 v_2(A)=v_2(B)=-\tfrac12 v_2(m+1).
\]
In the original gate variable \(U=B^2\), the second equality becomes \(v_2(U)=-v_2(m+1)\). The following additional conditions are necessary in this remaining pair. The even-start theorem is for each fixed window, with no common translation threshold established for the collection of all windows.

## Quarter-valuations in the remaining odd-start case

Let \(n=m+1\), \(e=v_2(n)\), \(x=A/B\), and \(h=1/B\), so \(x\) is a unit and \(v_2(h)=e/2\). Choose \(\pi^4=2\) and \(\delta^2=n\) in the algebraic closure. Put \(c=B\delta\), \(d=nA\), and \(D=1-c+d\). At any remaining common zero,
\[
 v_2(x+h-1)=\tfrac14,\qquad \overline{(x+h-1)/\pi}=1 .
\]
For the even order there are the following additional possibilities:
\[
 \begin{array}{ll}
 4\mid n:&v_2(c-1)=\tfrac14,\quad \overline{(c-1)/\pi}=1,\\
 n\equiv2\bmod4:&
 v_2(c-1)=\tfrac14,\quad \overline{(c-1)/\pi}=1,
 \quad\text{or}\quad v_2(c-1)=\tfrac12,\quad v_2(D)\ge1.
 \end{array}
\]
These conditions are necessary and do not establish exclusion in the indicated cases.

For odd \(j\), define \(M_k=\mathbb E(1+X_j)^k\) with \(X_j\sim\mathcal N(0,1/j)\). Gaussian integration by parts, or the formal exponential generating function \(e^{t+t^2/(2j)}\), gives \(M_{k+2}=M_{k+1}+(k+1)M_k/j\). Write
\[
 I_r=M_{4r}/2^r,\quad J_r=M_{4r+1}/2^r,\quad
 K_r=M_{4r+2}/2^{r+1},\quad L_r=M_{4r+3}/2^{r+1}.
\]
The same recurrence reads
\[
 \begin{split}
 K_r&=\tfrac12(J_r+(4r+1)I_r/j),&
 L_r&=K_r+(2r+1)J_r/j,\\
 I_{r+1}&=(1+(4r+3)/j)K_r+(2r+1)J_r/j,&
 J_{r+1}&=I_{r+1}+(4r+4)L_r/j .
 \end{split}
\]
Starting with \(I_0=J_0=1\), this proves \(I_r,J_r\) are odd and congruent modulo four, while \(K_r,L_r\) are integral. For \(j=1\), \(K_r\) is odd and \(L_r\) is even, by reducing the first two equations modulo two using \(J_r\equiv I_r\pmod4\). In particular, for every odd \(j\),
\[
 v_2(M_k)\ge (k-1)/4,\quad
 \text{with equality exactly when } k\equiv1\pmod4 .
\]

For the odd order \(j\equiv1\pmod4\), set \(t=x+h-1\) and \(U=1+X_j\). Reduction of \(G_j/B^j\) gives \(\bar x(\bar x+1)^{j-1}\) in the residue field, since only even Gaussian powers contribute. This excludes the unit \(x\) when \(j=1\). Otherwise it gives \(v_2(t)>0\). The unreduced expression is
\[
 \mathbb E[(t+(1-2h)U+hU^2)^j].
\]
Terms with \(k\) nonconstant factors, of which \(q\) use \(hU^2\), have valuation at least
\[
 (j-k)v_2(t)+(k-1)/4+q(v_2(h)+1/4).
\]
If \(0<v_2(t)<1/4\), the term with \(k=1,q=0\) has strictly least valuation, including comparison with the constant-only term \(t^j\). If \(v_2(t)>1/4\), the term with \(k=j,q=0\) is strictly least, since \(j\equiv1\bmod4\). This comparison applies also at \(t=0\). At \(v_2(t)=1/4\), division by \(\pi^{j-1}\) retains only \(q=0,\ k\equiv1\bmod4\) in the residue. For \(j=4r+1\) the result is \((\overline{t/\pi}^{\,4}+1)^r\), by \((1+z)^4=1+z^4\) in characteristic two. A zero forces \(\overline{t/\pi}=1\).

For the even order set \(Y\sim\mathcal N(0,1)\), \(U=1+Y\), and \(C=c-2\). The divided middle equation is
\[
 \mathbb E[(Y^2+cY+d)^n]=\sum_{k=0}^n\binom nk D^{n-k} V_k=0,
 \qquad V_k=\mathbb E[(CU+U^2)^k].
\]
Here \(d=c\delta x\) has valuation \(e/2\). Reduction gives \((c+1)^n\), so \(\bar c=1\) and \(C\) is a unit. Expanding \(V_k\) using the moments just proved for \(j=1\) yields \(v_2(V_k)\ge(k-1)/4\) for all \(k\). If \(4\mid k\), the term with \(M_k\) uniquely has least valuation \(k/4\). The only other candidate at \(k/4\) has \(M_{k+1}\), and its coefficient \(k C^{k-1}\) is even. For \(k=4r+2\), the only terms that can have valuation \(r+1\) use \(M_k\) or \(M_{k+2}\). The term using \(M_{k+3}\) has even binomial coefficient, and the one using \(M_{k+1}\) has higher moment valuation. Therefore
\[
 V_k/2^{r+1}=C^{k-2}(uC^2+v)+2W(C),
 \qquad u,v\in\mathbb Z_2^\times,\quad W\in\mathbb Z_2[C].
\]
For \(w=v_2(c-1)<1/2\), \(C^2+1\) has valuation \(2w\), since \(C-1\) also has valuation \(w\). Thus \(v_2(V_k)=r+1+2w\) if \(w<1/2\), and \(v_2(V_k)\ge r+2\) otherwise.

If \(4\mid n\), then \(\binom nk\) is even for \(4\nmid k\). Each coefficient \(\binom nk V_k\) has valuation at least \(k/4\), with strict inequality for \(4\nmid k\), and equality at \(k=0,n\). Consequently \(v_2(D)=1/4\) is necessary, and reduction after division by \(\pi^n\) gives \((\overline{D/\pi}^{\,4}+1)^{n/4}\). Since \(v_2(d)=e/2\ge1\), this gives the stated condition on \(c-1\).

If \(n\equiv2\bmod4\), \(\binom nk\) is even for odd \(k\), and \(\binom n{n-2}\) is odd. Write \(s=v_2(D)\). When \(w<1/2\), \(s=w\). For \(s<1/4\), the \(k=0\) term is strictly least by the coefficient bounds. For \(1/4<s<1/2\), the \(k=n-2\) term uniquely has least valuation \(2s+(n-2)/4\). Indeed, smaller \(4\)-divisible indices give strictly larger values, while all \(k\equiv2\bmod4\) give at least \(n/4+1/2+2s\), and the odd indices give at least \(s+1+(n-2)/4\). At \(s=1/4\) only \(4\)-divisible \(k\) contribute after division by \(\pi^n\), giving \(\overline{D/\pi}^{\,2}(\overline{D/\pi}^{\,4}+1)^{(n-2)/4}\). A zero forces the displayed quarter residue, and is impossible there when \(n=2\). Finally, if \(w\ge1/2\) and \(s<1\), then \(1/2\le s<1\). Again \(k=n-2\) is strictly least. The least bound for \(k\equiv2\bmod4\) is \((n-2)/4+2\), and the bound for odd \(k\) stays \(s+1+(n-2)/4\). Hence \(s\ge1\) is necessary. Since \(v_2(d)=1/2\), this requires \(w=1/2\).
