# Signed norm factors with variable squarefree kernel

Fix a dyadic integer \(q=2k\ge4\). Count primitive positive triples \(a+b=c\le X\) admitting
\[
 a=D x^2,\qquad b=B R^{2k},\qquad c=z^2,
\]
where \(D>1\) is squarefree and \(B\) is an actual square. Assume \(D\le X^{f\lambda}\), \(\operatorname{rad}(B)\le X^{g\lambda}\), and \(\operatorname{rad}(a)\le X^{\ell\lambda}\). The same method allows any additional radical caps on \(b,c\). The original count satisfies
\[
                  \ll_{k,\varepsilon} X^{((\ell+f)/2+g)\lambda+\varepsilon}.
 \tag{RN1}
\]
No short-term height assumption is used. In particular, in the repaired target band with \(\operatorname{rad}(c)\ge c^{(1/3-\xi)\lambda}\), the exponent is at most \(((\ell+f)/2+g)/((1/3-\xi)M)\) up to arbitrary small exponent allowances. Here the displayed original caps can be given relative to \(c\) rather than \(X\) by taking dyadic shells.

## Primitive ideals and an actual anchor

Fix \(D\) and the coefficient \(B\). Let \(K=\mathbb Q(\sqrt D)\) and \(\alpha=z+x\sqrt D\). Then \(N\alpha=b>0\); both embeddings are positive since \(z^2>D x^2\). Away from primes above \(2D\), the ideals of \(\alpha,\bar\alpha\) are coprime by \(\gcd(z,x)=1\). No prime of \(D\) divides \(b\), by primitivity of the original triple. Thus
\[
                  (\alpha)=\mathfrak d\mathfrak b^{2k},
\]
where \(\mathfrak d\) is supported on \(2B\), with each prime-ideal valuation below \(2k\). Indeed at good primes all of the rational norm exponent comes from one ideal and is a multiple of \(2k\) outside \(B\). We may prescribe all residual valuations at bad primes as well. There are \(O_k(1)^{\omega(2B)}=X^{o(1)}\) choices of \(\mathfrak d\) for fixed \(D,B\).

The classes of \(\mathfrak b\) belong to a coset of \(\mathrm{Cl}(K)[2k]\) when nonempty. For real quadratic \(K\), \(|\mathrm{Cl}(K)[2]|\le 2^{1+\omega(D_K)}\). To see this, a 2-torsion class has a representative \(\mathfrak e=(v)\bar{\mathfrak e}\) with \(N(v)=\pm1\). If a norm-minus-one unit exists then we can arrange \(N(v)=1\). Otherwise the sign is well-defined on 2-torsion classes and multiplicative, so at most two cosets arise. In the norm-one case write \(v=w/\bar w\), by \(w=1+v\) or \(w=\sqrt D\) if \(v=-1\). Then \(\mathfrak e/(w)\) is fixed by conjugation. Rational ideals remove inert primes and pairs of split primes; what remains at ramified primes is determined by parity of their exponents. This proves the bound. Iterating shows \(|\mathrm{Cl}(K)[2k]|=X^{o(1)}\), since \(2k\) is fixed dyadic.

For a fixed class choose any actual reference solution \(\alpha_*,\mathfrak b_*\). The quotient ideal \(\mathfrak b/\mathfrak b_*\) is principal; choosing a generator makes \(\alpha/\alpha_*\) equal its \(2k\)-th power times a totally positive unit. The unit modulo \(2k\)-th powers of totally positive units gives at most \(2k\) classes. One can use the elementary discreteness of the first embedding of norm-one units, as in the real quadratic norm count for Erdős #943. For each nonempty class of that partition choose an actual solution \(\eta=\alpha_0\) and its associated ideal \(\mathfrak b_0\). Then all solutions in this class have
\[
       \alpha=\eta\gamma^{2k},\qquad
       \gamma\in\mathfrak b_0^{-1},\qquad
       \sqrt{N\eta}=s\in\mathbb Z_{>0}.
\]
Here \((\eta)=\mathfrak d\mathfrak b_0^{2k}\), and the norm is square since \(b=B R^{2k}\) and \(B\) is square. The anchor is an actual solution of polynomial height, so no bound on a fundamental unit or an arbitrarily generated principal ideal is needed.

Use a \(\mathbb Z\)-basis \(\omega_1,\omega_2\) of \(\mathfrak b_0^{-1}\), of polynomial rational-coordinate height. Such a basis follows by Hermite normal form on \(\overline{\mathfrak b_0}=N\mathfrak b_0\,\mathfrak b_0^{-1}\) in the basis of \(\mathcal O_K\), since \(N\mathfrak b_0\le X^{1/(2k)}\). Write \(\gamma=u\omega_1+v\omega_2\). Both embeddings of \(\gamma\) have polynomial height bounds above in absolute value because the ratio \(\alpha/\eta\) does, and the embeddings of \(\eta\) are bounded below by inverse polynomial size from its nonzero integral norm. Thus \(u,v\) are polynomially bounded too.

## Rational factor split and cancellation of the anchor ideal

As polynomials in the two integer variables \(u,v\), put
\[
          P_\pm=\eta\gamma^k\pm s\bar\gamma^k .
\]
Their coefficients lie in \(\mathfrak b_0^k\). For the first term use \((\eta)=\mathfrak d\mathfrak b_0^{2k}\) and \(\gamma\) with coefficients in \(\mathfrak b_0^{-1}\). For the second use that the rational ideal \((s)\) is divisible by \((N\mathfrak b_0)^k=(\mathfrak b_0\bar{\mathfrak b_0})^k\); indeed \(s=(N\mathfrak d)^{1/2}(N\mathfrak b_0)^k\), with integral square root since \(N\eta\) is square. Conjugation acts by
\[
                   \bar P_\pm=\pm(s/\eta)P_\pm .
\]
Consequently the ratio of any two nonzero coefficients is rational. Write \(P_\pm=\sigma_\pm F_\pm\) with \(F_\pm\) primitive integer homogeneous forms of degree \(k\). By Bezout on the integer coefficients of \(F_\pm\), each \(\sigma_\pm\) itself lies in \(\mathfrak b_0^k\). Their rational form coefficients have polynomial height uniformly, since the original algebraic coefficients and their ratios have polynomial height and the degree is fixed.

The exact product is
\[
                   P_+P_-=2\eta\sqrt D\,x .
\]
Outside primes of \(2DB\), \(\sigma_+\sigma_-/\eta\) is integral locally, since \(\mathfrak d\) is supported on \(2B\). Therefore any rational prime outside \(2DB\) dividing either \(F_\pm(u,v)\) divides \(x\). If such a prime \(p\) divides both, then \(P_+\) and \(P_-\) both lie in \(p\mathfrak b_0^k\) locally, and hence so does \(\eta\gamma^k\). At each prime ideal above this unramified \(p\), this forces \(v_{\mathfrak p}(\gamma)\ge -v_{\mathfrak p}(\mathfrak b_0)+1\): the extra \(1/k\) in the valuation inequality rounds up. Thus \(\eta\gamma^{2k}=\alpha\) has valuations at least \(2k\) at all primes above \(p\), forcing \(p\) to divide both \(z,x\) since \(p\nmid2D\). This contradicts primitivity. In particular there is no loss from the possibly large norm or radical of the chosen anchor ideal. The outside radicals of the two form values are disjoint and their product divides \(\operatorname{rad}(x)\) outside the fixed support \(2DB\), which is at most \(\operatorname{rad}(a)/D\).

Both forms are squarefree as binary forms, since the projective ratio \(\gamma/\bar\gamma\) is a nonconstant fractional linear function and the equations \((\gamma/\bar\gamma)^k=\pm s/\eta\) have distinct roots. Each form has degree \(k\ge2\). A specified nonzero form value gives subpower many integer pairs of polynomial height uniformly. If the form has a nonlinear irreducible factor, fix that factor's value by its divisibility into the specified value, then use the uniform fixed-value norm count for binary forms. Otherwise fix values of two distinct linear factors and solve their independent linear system. The form values here are nonzero because their product is proportional to \(x>0\).

It follows that, for fixed \(D,B\), a smallest outside radical among the form values costs at most \((X^{\ell\lambda}/D)^{1/2}\) choices. Its parts supported on \(2DB\) cost subpower by the fixed-support estimate, and give its exact value. Sum over the subpower many choices of ideal data, then over \(D\le X^{f\lambda}\) using \(\sum_{D\le Z}D^{-1/2}\ll\sqrt Z\). Enumerating \(B\) by its radical gives \(X^{g\lambda+o(1)}\) possibilities. This proves (RN1).

For arbitrary fixed even \(2k\), replacing the dyadic torsion bound by the elementary real quadratic class-number bound \(D^{1/2+o(1)}\) gives exponent \((\ell/2+f+g)\lambda+\varepsilon\). The dyadic case is sharper. For example \(\ell=1/3+1/100\), \(f=1/12\), \(g=1/48\), and target lower band \(1/3-1/100\) give repaired coefficient less than \(0.725/M\). All three original caps may hold. The condition is on genuine primitive parity and target/summand order structure, and cannot be forced by an artificial common repair.

## Positive-definite target variant

The same bound holds with the near-power condition assigned to the actual target:
\[
             a=x^2,\qquad b=D y^2,\qquad c=B R^{2k},
\]
where \(B\) is a square and \(D\) is positive squarefree. Here use \(\operatorname{rad}(b)\le X^{\ell\lambda}\), and \(D,\operatorname{rad}(B)\) have the same caps \(f,g\). Use \(\alpha=x+y\sqrt{-D}\) in the imaginary quadratic field, with \(N\alpha=c\). Its ideal and its conjugate are coprime away from \(2D\), and primes of \(D\) do not divide \(c\) by primitivity. The residual ideal is once again supported on \(2B\), and each chosen solution's norm is square. In the imaginary field the 2-torsion proof has no possible negative norm sign, and the unit quotient is automatically bounded because the whole unit group is finite. The same actual-anchor argument gives
\[
 P_+P_-=2\eta\sqrt{-D}\,y.
\]
Primitive integral forms descend and their anchor-ideal contents cancel exactly as in the signed case. In particular their outside-support radicals have product at most \(\operatorname{rad}(b)/D\), with no angular or shortness hypothesis. This proves the same counting and repaired exponents for the genuine even-power target with a square defect coefficient, one square summand, and the partial parity condition on the other. The hypotheses apply to primitive valuations, not to a repaired representative.
