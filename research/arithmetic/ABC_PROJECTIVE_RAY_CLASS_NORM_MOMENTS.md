# Projective ray class moments for fixed-degree norm forms

A full divisor of a primitive binary-form input constrains the ideal of the associated norm element modulo a rational scalar. In degree at least three, a divisor of the form value also gives such a constraint by retaining the prime ideals away from the unique prime carrying that value. A compact group keeps the residue constraint together with the centered logarithms of the embeddings. These conditions make high moments of the norm-factor dictionaries diagonal up to class-group torsion. For weighted \(q\)-th powers this gives a strict saving at the balanced abc endpoint when the two \(q\)-free coefficient radicals have sufficiently small combined mass.

## Counting primitive inputs of an irreducible binary form

Fix \(d\ge2\), \(B_0>0\), \(\nu>0\) and nonnegative \(g,\ell,\mu\). For every \(\eta>0\) and sufficiently small fixed \(\lambda>0\), let \(F\in\mathbb Z[U,V]\) be irreducible, homogeneous of degree \(d\), with coefficients of magnitude at most \(X^{B_0}\). Put \(a_0=F(1,0)\) and \(a_d=F(0,1)\). Both are nonzero. Let \(\theta\) be a root of the monic integer polynomial \(a_0^{d-1}F(Z/a_0,1)\), and let \(K=\mathbb Q(\theta)\). Suppose
\[
          |\operatorname{Disc}(K)|\ll_d X^{g\lambda}.
\]
The discriminant \(\operatorname{Disc}(F)\) below means the discriminant of \(F(Z,1)\). Count primitive positive \(u,v\le X^{B_0}\) with
\[
 \prod_{\substack{p\mid F(u,v)\\p\nmid P}}p\le X^{\ell\lambda}
 \quad\text{and}\quad
 M_1\mid v,\quad M_2\mid u,\quad M=M_1M_2\ge X^\nu,\quad
 \operatorname{rad}(M)\le X^{\mu\lambda},\quad
 \gcd(M,2a_0a_d\operatorname{Disc}(F))=1.                  \tag{PR1}
\]
Alternatively, when \(d\ge3\), put the divisor on the form value:
\[
 \prod_{\substack{p\mid F(u,v)\\p\nmid PM}}p\le X^{\ell\lambda},
 \quad M\mid F(u,v),\quad M\ge X^\nu,\quad
 \operatorname{rad}(M)\le X^{\mu\lambda},\quad
 \gcd(M,2a_0a_d\operatorname{Disc}(F))=1.                 \tag{PR1'}
\]
In this version the cap counts the radical outside the divisor. The divisor record may vary in either version. The optional support \(P\) is a positive integer at most \(X^{B_0}\) prescribed after the form and divisor record are fixed. Uniformly in the field and the coefficients, the count for fixed \(F\) is
\[
       \ll_{d,B_0,g,\ell,\mu,\nu,\eta,\lambda,\varepsilon}
         X^{(\mu+\ell/2+\eta)\lambda+\varepsilon}.      \tag{PR2}
\]
Every sufficiently small bound on \(\lambda\) used for this estimate depends only on the fixed parameters. A fixed divisor record can be counted without the cost \(\mu\). The proof retains the separate \(e\)-torsion factors of the class group, given in (PR5) below.

Write \(\alpha=a_0u-\theta v\), so \(N_{K/\mathbb Q}\alpha=a_0^{d-1}F(u,v)\). The element \(\alpha\) is nonzero, since \(F\) is irreducible of degree at least two and the pair is primitive. Both the coefficients of its defining data and all conjugates of \(\alpha\) are of polynomial height in \(X\), uniformly at fixed \(d,B_0\). The nonzero integral norm makes every conjugate's absolute value at least an inverse power of \(X\). At a prime dividing \(F(u,v)\) and outside \(2a_0a_d\operatorname{Disc}(F)\), both \(u,v\) are units. The displayed monic polynomial for \(\theta\) has squarefree reduction at this prime because its discriminant has no prime divisors outside \(a_0\operatorname{Disc}(F)\). The ideal of \(\alpha\) can therefore use only the degree-one prime specified by \(\theta\equiv a_0u/v\). This specifies one prime even when there are several degree-one primes above the same rational prime. There is at least one other prime above it since the prime is unramified and \(d\ge2\).

For the initial enumeration, a positive integer bounded by a fixed power of \(X\), on prescribed prime support of radical at most a fixed power of \(X\), has \(X^{o_X(1)}\) possibilities. For each fixed \(s>0\), apply the Rankin bound with Euler product \(\prod_p(1-p^{-s})^{-1}\) over this support. Its logarithm is \(o_X(\log X)\) uniformly by splitting the primes at a sufficiently large constant and using \(\omega(P)=O(\log X)\) for the large primes. Taking \(s\) arbitrarily small proves the bound. Ideals of polynomial norm supported above this same support also have \(X^{o_X(1)}\) possibilities, since for each norm \(n\) the number of exponent allocations upstairs is at most \(\tau(n)^d\). Thus fixing the divisor record costs \(X^{\mu\lambda+o_X(1)}\). For each fixed record include all primes over \(2a_0a_d\operatorname{Disc}(F)PM\) in a fixed ideal factor of \((\alpha)\), at subpower cost. In (PR1) no prime over \(M\) divides \((\alpha)\), since either \(\alpha\) or \(\alpha/\theta\) reduces to a nonzero rational scalar. In (PR1') the fixed factor records the unique degree-one prime \(\mathfrak p_p\) over each \(p\mid M\) carrying the valuation of \(F(u,v)\).

At every rational prime outside the fixed support that has a degree-one prime, choose one such ideal as the allowed ideal. Independent uniform choices give a choice retaining at least an \(X^{-o_X(1)}\) proportion of the rows, since each row uses at most \(O(\log X/\log\log X)\) rational primes. These choices can be restricted to the finitely many primes up to a polynomial bound for the norms. Work in one such choice, bounding the original count by the loss factor times the bound for its retained rows. Fix a degree cutoff \(J\) large enough that the full part of \(|F(u,v)|\) supported on valuations \(>J\) has radical at most \(X^{\eta\lambda/8}\), using the polynomial height bound for \(F(u,v)\). Fix that whole part and its corresponding ideal factor at cost \(X^{\eta\lambda/8+o_X(1)}\). The remaining equality is
\[
            (\alpha)=\mathfrak d\prod_{1\le e\le J}\mathfrak a_e^e ,
\]
with \(\mathfrak d\) fixed, and each \(\mathfrak a_e\) squarefree and supported on the allowed degree-one prime ideals outside the fixed support. Bin the norms with \(N\mathfrak a_e\le X^{\lambda\vartheta_e}\) in a finite cap mesh such that \(H_0=\sum_e\vartheta_e\le \ell+\eta/4\). The cutoff and the mesh are fixed in the height limit.

## The compact group and the diagonal moment

For (PR1) use the modulus \(\mathfrak m=M\mathcal O_K\). For (PR1') use
\[
            \mathfrak m=\prod_{p^k\parallel M}
                     \prod_{\substack{\mathfrak p\mid p\\\mathfrak p\ne\mathfrak p_p}}
                              \mathfrak p^{k}.
\]
Every retained ideal factor, including \(\mathfrak d\), is prime to \(\mathfrak m\). Let \(I(\mathfrak m)\) be the discrete group of all fractional ideals prime to \(\mathfrak m\), and put
\[
 V_\infty=\{(t_\sigma)_\sigma:\ t_{\bar\sigma}=t_\sigma,\
                            \textstyle\sum_\sigma t_\sigma=0\},
 \qquad
 L(\gamma)_\sigma=\log|\sigma(\gamma)|-\frac1d\log|N_{K/\mathbb Q}\gamma|.
\]
Here the index runs over all embeddings, counting a complex-conjugate pair twice. For \(\gamma\in K^\times\), the condition that \(\gamma\) reduce to a rational unit modulo \(\mathfrak m\) means \(\gamma\equiv a\bmod\mathfrak m\mathcal O_{K,\mathfrak p}\) for every \(\mathfrak p\mid\mathfrak m\), for a single \(a\in\mathbb Z\) coprime to \(M\). Define
\[
 G=\bigl(I(\mathfrak m)\times V_\infty\bigr)/
 \{((\gamma),-L(\gamma)):\gamma\text{ reduces to a rational unit modulo }\mathfrak m\}.
                                                               \tag{PR3}
\]
There is no real positivity condition. In each ideal component the subgroup in the denominator, when nonempty, is a translate of the log lattice of the units reducing to rational units. Dirichlet's theorem gives a full lattice since these units contain the finite-index subgroup of units congruent to one. The denominator is therefore closed. Put \(g(\mathfrak i)=[(\mathfrak i,0)]\) and \(h(t)=[(1,t)]\). The image \(H=h(V_\infty)\) is a torus of dimension \(r=\dim V_\infty\). The quotient \(G/H\) is a quotient of the ordinary ray class group for modulus \(\mathfrak m\), so \(H\) is open of finite index and \(G\) is compact. The same assertions hold with \(H=\{0\}\) when \(r=0\). [Degree-zero Arakelov groups, Schoof, Section 2](https://arxiv.org/pdf/0801.3835) · [Ray class exact sequence, Stevenhagen, Sections 3–4](https://websites.math.leidenuniv.nl/algebra/Stevenhagen-CFT.pdf)

For each fixed \(e\ge1\) its \(e\)-torsion satisfies
\[
     |G[e]|\le e^r\,|\mathrm{Cl}(K)[e]|\,
                         |(\mathcal O_K/\mathfrak m)^\times[e]|
             \le |\mathrm{Cl}(K)[e]|X^{o_X(1)}.              \tag{PR4}
\]
Indeed the order of \(e\)-torsion of a finite abelian group equals the order of the quotient by \(e\)-multiples, so passing to a quotient cannot increase it. Apply this to \(G/H\) and then the kernel in the ray exact sequence. Multiplication by \(e\) on \(H\) has kernel of order \(e^r\). For a rational prime \(p\mid M\), all local factors are odd and unramified, since \(p\nmid 2a_0\operatorname{Disc}(F)\). For each factor of residue degree \(f\) and modulus \(p^k\), the residue multiplicative group contributes at most \(e\) roots. On principal units the usual binomial expansion gives \(v(u^e-1)=v(u-1)+v_p(e)\) in the unramified local valuation. Apply the expansion first to the prime-to-\(p\) part of \(e\) and then successively to multiplication by \(p\). Hence the principal unit kernel has size at most \(p^{f v_p(e)}\), regardless of \(k\). The product over factors is at most \(e^{2d}\) per rational prime. Since \(M\) has polynomial height in \(X\), the result is subpower uniformly, including nonsquarefree \(M\).

In (PR1), choose \(\beta\in\mathcal O_K\) with \(\beta\equiv1\bmod M_1\mathcal O_K\) and \(\beta\equiv\theta\bmod M_2\mathcal O_K\). The two residue values are units at the relevant primes, and the Chinese remainder theorem shows that \(\alpha/\beta\) reduces to a rational unit modulo \(\mathfrak m\). In (PR1'), the completion of \(K\) at each \(\mathfrak p_p\) is \(\mathbb Q_p\), since it is unramified of degree one. For \(p^k\parallel M\), let \(r_p\) be the residue of \(\theta\) in \(\mathcal O_{K,\mathfrak p_p}/\mathfrak p_p^k\cong\mathbb Z/p^k\mathbb Z\). Then \(a_0u/v\equiv r_p\bmod p^k\). With \(r_*\in\mathbb Z\) given by the Chinese remainder theorem, \(\alpha\equiv v(r_*-\theta)\bmod M\mathcal O_K\). Put \(\beta=r_*-\theta\), which is a unit at each prime of \(\mathfrak m\) by squarefree reduction and the exclusion of the active prime. Hence \(\alpha/\beta\) reduces to a rational unit modulo \(\mathfrak m\) here as well.

In both cases \(g((\alpha))=g((\beta))+h(L(\alpha)-L(\beta))\). The vector \(L(\alpha)\) is \(O_{d,B_0}(\log X)\). With \(W=\nu\log X/32\), write \(C_W=\{t\in V_\infty:\|t\|_\infty\le W\}\). Cover the range of \(L(\alpha)\) by \(O_{d,B_0,\nu}(1)\) translates of \(C_{W/2}\). In each translate, the sum \(\sum_e e\,g(\mathfrak a_e)\) is constrained to a fixed translate of \(h(C_{W/2})\).

For a positive Fourier test on this set use the autocorrelation of \(\mathbf 1_{h(C_W)}\), divided by the probability Haar measure of \(h(C_W)\) in \(G\). Write \(K_0\) for this function. It is supported on \(h(C_{2W})\), bounded above by one, and at least \(5^{-r}\) on \(h(C_{W/2})\). To see the last assertion, \(C_W\) is covered by at most \(5^r\) translates of \(C_{W/2}\) by packing \(C_{W/4}\) about a maximal \(W/2\)-separated set in \(C_W\). For \(t\in C_{W/2}\) the intersection \(h(C_W)\cap(h(C_W)+h(t))\) contains \(h(C_{W/2})\). This argument applies regardless of the size and shape of the log lattice. Parseval gives \(\widehat K_0(\chi)\ge0\) and \(\sum_\chi\widehat K_0(\chi)=1\). Its Fourier series is absolutely convergent, and translates of this test change the coefficients only by phases.

For each nontrivial slot \(e\), let \(\mathcal A_e\) be the dictionary of possible ideals of norm at most \(X^{\lambda\vartheta_e}\) and put \(F_e(\chi)=\sum_{\mathfrak a\in\mathcal A_e}\chi(e\,g(\mathfrak a))\). The orientation fixed above makes each ideal determined by its norm. Fix \(0<\xi\le1/2\) and write \(m_e=\lceil H_0/(2\xi\vartheta_e)\rceil\), so \(m_e\vartheta_e\le(1+1/(2\xi))H_0\) and \(\sum_e(2m_e)^{-1}\le\xi\). Slots of zero cap are omitted, and if all are omitted there is only one tuple to count. A \(2m_e\)-th moment against \(\widehat K_0\) tests pairs of products of \(m_e\) ideals for
\[
                        e\,g(\mathfrak A/\mathfrak B)\in h(C_{2W}).
\]
For each fixed \(\mathfrak A\), the possible \(g(\mathfrak B)\) lie in at most \(|G[e]|\) translates of \(h(C_{2W/e})\), since subtracting \(h(t/e)\) from an \(e\)-preimage of \(h(t)\) gives an element killed by \(e\). For two products in the same translate, \(\mathfrak B_1/\mathfrak B_2=(\gamma)\), where \(\gamma\) reduces to a rational unit modulo \(\mathfrak m\) and \(\|L(\gamma)\|_\infty\le4W/e\). Both assertions are direct consequences of the kernel in (PR3).

The integer \(D_2=N\mathfrak B_2\) lies in \(\mathfrak B_2\), so \(D_2\gamma\) is integral. With \(N\mathfrak B_i\le X^{(1+1/(2\xi))H_0\lambda}\), all embeddings satisfy
\[
 |\sigma(D_2\gamma)|
     \le (N\mathfrak B_1)^{1/d}(N\mathfrak B_2)^{1-1/d}\exp(4W/e)
     \le X^{(1+1/(2\xi))H_0\lambda+\nu/8}.
\]
Write \(\widetilde K\) for a normal closure of degree \(D\) over \(\mathbb Q\), and fix two embeddings \(\sigma,\tau:K\hookrightarrow\widetilde K\). For (PR1), the difference \(\Delta=\sigma(D_2\gamma)-\tau(D_2\gamma)\) lies in \(M\mathcal O_{\widetilde K}\), since \(D_2\gamma\) reduces modulo \(M\mathcal O_K\) to a rational integer. If \(\Delta\ne0\), this gives \(|N_{\widetilde K/\mathbb Q}\Delta|\ge M^D\).

For (PR1'), take \(p^k\parallel M\). The field \(\widetilde K\) is unramified at \(p\) since each completion of \(K\) at \(p\) is unramified and the compositum of their local conjugates is unramified. Among the primes of \(\widetilde K\) above \(p\), the fraction restricting under \(\sigma\) to the distinguished prime \(\mathfrak p_p\) is \(1/d\). Indeed if the common residue degree in \(\widetilde K\) is \(f_L\), the total number of primes is \(D/f_L\) and the number above that degree-one prime of \(\sigma(K)\) is \((D/d)/f_L\). The same holds for \(\tau\). At a fraction at least \(1-2/d\) of primes above \(p\), both restrictions are retained in the modulus. The two conjugates of \(D_2\gamma\) reduce there to the same rational integer modulo the \(k\)-th prime power. Thus if \(\Delta\ne0\) then
\[
    |N_{\widetilde K/\mathbb Q}\Delta|\ge M^{D(1-2/d)}.
\]
For the form-value case \(d\ge3\) gives \(1-2/d\ge1/3\). In both cases the lower bound contradicts \(|\sigma'(\Delta)|\le2X^{(1+1/(2\xi))H_0\lambda+\nu/8}\) for every embedding \(\sigma'\) of \(\widetilde K\), provided \(\lambda\) is sufficiently small in terms of \(\ell,\eta,\nu,\xi\). Every such difference vanishes, so \(\gamma\in\mathbb Q^\times\). At each allowed rational prime the ideal quotient \(\mathfrak B_1/\mathfrak B_2\) is supported on at most one prime upstairs, while a nonzero valuation of a rational \(\gamma\) would use every prime upstairs. Therefore \(\mathfrak B_1=\mathfrak B_2\). This is the diagonal step.

A fixed ideal product has subpower many ordered decompositions since \(m_e\) is fixed at the height limit and the product norm is polynomially bounded. There are at most \(X^{m_e\lambda\vartheta_e}\) choices of the first tuple. Fourier inversion now gives
\[
 \left(\sum_\chi\widehat K_0(\chi)|F_e(\chi)|^{2m_e}\right)^{1/(2m_e)}
       \le X^{\lambda\vartheta_e/2+o_X(1)}
                            |\mathrm{Cl}(K)[e]|^{1/(2m_e)}.
\]
The absolute Fourier expansion with Hölder bounds the number of tuples in the fixed translate by a constant times the product of these norms. The unused part of the Hölder exponents is filled by the constant function one in the probability measure. Thus the resulting class-group factor is
\[
                  \prod_e|\mathrm{Cl}(K)[e]|^{1/(2m_e)}.
                                                               \tag{PR5}
\]
This factor is at most \(|\mathrm{Cl}(K)|^\xi\),
which is \(O_{d,\varepsilon}(X^{g\lambda\xi/2+\varepsilon})\)
by Minkowski's theorem and the ideal divisor bound.
Every class has an integral representative of norm
\(O_d(\sqrt{|\operatorname{Disc}(K)|})\), and each integer
\(n\) occurs as norm at most \(\tau(n)^d\) times. Fix
\(\xi\) small enough that \(g\xi/2<\eta/4\) before taking
\(\lambda\) small for the amplified norm cutoff. This
absorbs the entire class-group factor in the allowance.
[Amplified moments in full ray groups](ABC_BINARY_FORM_RAY_AMPLIFICATION.md#amplified-moments-and-the-class-number-price)

With fixed ideal factors the possible \(\alpha\)'s differ by units. Unit log vectors within a bounded distance are uniformly separated at fixed degree, because a unit with bounded conjugate moduli belongs to a fixed finite set by the coefficient bounds for its minimal polynomial. The zero-log units are uniformly finite for the same reason. The allowed log range \(O_{d,B_0}(\log X)\) therefore gives \(O_{d,B_0}((1+\log X)^d)\) generators. Each \(\alpha=a_0u-\theta v\) determines \(u,v\) since \(\theta\) is irrational. Adding the costs of the records, the tail, the cap mesh and (PR5) proves (PR2).

### Combining divisors on the inputs and the form value

One may use an input divisor \(M_{\rm in}=M_1M_2\) and a form-value divisor \(M_{\rm out}\) together. Assume \(d\ge3\), the same coprimality with \(2a_0a_d\operatorname{Disc}(F)\), the radical cap \(\operatorname{rad}(M_{\rm in}M_{\rm out})\le X^{\mu\lambda}\), and height condition
\[
                    M_{\rm in}M_{\rm out}^{\,1-2/d}\ge X^\nu .
\]
Use the radical cap on \(F(u,v)\) outside \(PM_{\rm out}\). Then (PR2) holds with these data. The two moduli are coprime because \(F(u,v)\) is prime to the input divisor. Combine the two constructions of \(\beta\) by the Chinese remainder theorem. In the diagonal step the prime powers from the input divisor contribute \(M_{\rm in}^D\) to the norm of the conjugate difference, and those from the form-value divisor contribute \(M_{\rm out}^{D(1-2/d)}\). Their product gives the required contradiction. The record enumeration and torsion bounds apply to the union of the prescribed supports.

### The target-divisor statement fails in degree two

Take \(F(U,V)=U^2+V^2\) and \(M=5^k\) with \(X^\nu\le5^k<5X^\nu\), where \(0<\nu<1\) is fixed. Choose primitive positive \(a,b\le \lfloor X^{\ell\lambda/2}/4\rfloor\) with \(a\equiv1\bmod5\) and \(b\equiv0\bmod5\). Möbius inversion in these two progressions gives \(\gg X^{\ell\lambda}\) choices. For each divisor \(d\) prime to \(5\), the unsigned main term for \(d\mid a,b\) is \(X^{\ell\lambda}/(400d^2)\) up to the rounding of the box side, and the summed error is \(O(X^{\ell\lambda/2}\log X)\). The sum of the main-term coefficients is positive.

Multiply \(a+ib\) by \((2+i)^k\) and rotate by a Gaussian unit to obtain positive real and imaginary parts \(u,v\). They are coprime. A rational prime other than \(5\) is excluded by invertibility of \(2+i\) modulo that prime. At \(5\) the product uses only one of the two conjugate primes, since its original norm was prime to \(5\). A coordinate cannot vanish because the pair is primitive and the norm exceeds one. Each rotated pair has at most four preimages. Its norm is divisible by \(M\) and the radical outside \(5\) is at most \(a^2+b^2\le X^{\ell\lambda}\). Also \(u,v\le X\) for sufficiently small fixed \(\lambda\). Thus the \(\ell\) coefficient is attained with a fixed divisor record and constant field discriminant, which rules out the \(\ell/2\) coefficient of (PR2) under (PR1') in degree two.

## Weighted powers and the balanced endpoint

Let two complete letters of a positive primitive triple \(a+b=c\le X\) be
\[
                      U=A x^q,\qquad V=B y^q,\qquad q\ge2,
\]
with \(A,B\) their canonical \(q\)-free coefficients, and call the remaining letter \(T\). Use \(T=V-U\) with \(V>U\) for opposite letters and \(T=V+U\) for same-side letters. Fix \(q\) and suppose \(\operatorname{rad}(AB)\le X^{f\lambda}\), \(\operatorname{rad}(T)\le X^{\ell\lambda}\), and a divisor record as in (PR1) exists on the roots \(x,y\) with \(M\) prime to \(2qAB\). If the binomial \(B u^q\pm A v^q\) has a nonlinear irreducible factor of degree \(d\), then
\[
                   \#\{(a,b,c)\}
             \ll_{q,f,\ell,\mu,\nu,\eta,\lambda,\varepsilon}
                X^{(\mu+\ell/2+f+\eta)\lambda+\varepsilon}.
                                                               \tag{PR6}
\]
Indeed \(A,B\) cost \(X^{f\lambda+o_X(1)}\). A primitive integral factor has polynomial-height coefficients by compactness of the coefficient unit spheres in fixed degree and the nonvanishing of a product of two nonzero polynomials. Its end coefficients divide the binomial's end coefficients, its discriminant divides the binomial discriminant by the discriminant product formula, and its value divides \(T\). Thus \(M\) is prime to all its excluded primes. The field of the factor is generated by a root of \(Z^q\pm A/B\), or equivalently by its \(B\)-multiple, which is a root of \(Z^q\pm A B^{q-1}\). The field discriminant satisfies
\[
                  |\operatorname{Disc}(K)|
                       \ll_q\operatorname{rad}(AB)^{d-1}.
\]
At \(p\nmid q\) adjoining the \(q\)-th roots of unity is unramified and the subsequent Kummer extensions have degree dividing \(q\). Hence the factor field is tame there, with discriminant exponents \(\sum_{\mathfrak p\mid p} f_{\mathfrak p}(e_{\mathfrak p}-1)\le d-1\), and is unramified outside \(qAB\). At \(p\mid q\), the discriminant of the full monic polynomial has bounded valuation because \(A,B\) are \(q\)-free. The discriminant of the monic minimal polynomial of the chosen root divides it. These bounds verify (PR2) with \(g=(d-1)f\).

A nonlinear factor always exists for \(q\ge3\), because the roots of \(Z^q\pm A/B\) cannot all be rational. For \(q=2\) it exists with the same-side sign, and with the opposite sign it fails precisely when the two primitive letters are exact squares. Mixed fixed powers \(U=A x^j,V=B y^k\) also enter this result when \(j,k\) have a common divisor \(q\ge2\). Use the canonical coefficients and roots of the resulting \(q\)-th powers at order \(q\). The hypothesis is the actual kernel cap at that order.

There is a target-divisor version of (PR6) when the chosen factor has \(d\ge3\) and \(M\) divides its value. In that version \(\ell\) bounds the radical of the factor value outside \(M\), and the divisor may come entirely from the third letter. If the full cap on \(T\) is \(X^{\ell\lambda}\) and the modulus radical is in a narrow bin of mass \(s\), the modulus and factor together cost \(\ell/2+s/2\), up to the bin width and \(\eta\). The corresponding bound follows by applying (PR1') with the remaining cap \(\ell-s\).

Consequently there is a strict improvement below the small-budget coefficient \(1/3\) on the unrestricted-height sector
\[
        \prod_{\substack{p\mid UV\\q\nmid v_p(UV)}} p\le c^{f\lambda},
        \qquad f<\frac16,
        \qquad \operatorname{rad}(abc)<c^\lambda .             \tag{PR7}
\]
More precisely, for each fixed \(q,f\) as above there is \(\delta(q,f)>0\) such that the sector is \(O_{q,f,\lambda,\varepsilon}(X^{(1/3-\delta(q,f))\lambda+\varepsilon})\) for all sufficiently small fixed \(\lambda\). To prove it when a nonlinear factor exists, take a violating sequence of fixed-budget inner height limits followed by \(\lambda\to0\). Work in a near-saturating anatomic box. Its actual total radical
mass divided by the ambient budget tends to one, by the limiting
one-third bound applied to its smaller budget if the ratio were
bounded away from one. The one-letter quarter inequality gives balanced caps \(1/3+o(1)\) on all three letters. The nonatomic mass on the chosen pair is at least \(1/3\). Removing \(q\)-defect primes costs at most \(f<1/3\) of this mass, and removing primes dividing \(2q\) has vanishing cost. To formalize this, split the scaled degree measures by divisibility of each exact index by \(q\) and pass to weak limits of the submeasures. The divisible submeasures on the pair still have positive nonatomic mass. A compact positive-degree continuity band of arbitrarily small positive mass on one of the chosen letters gives a full divisor of its root. Its radical cap \(\mu>0\) can be made arbitrarily small, and its fixed positive lower numerical weight gives \(\nu>0\), after division by \(q\) for the root. The band lies below the packed cutoff, whose scaled degree tends to
infinity, so it gives genuine exact valuation slices. The band and the
constants are fixed before applying the small-budget threshold and
taking the ensuing height limits. Apply (PR6) with the third-letter cap approaching \(1/3\). Its coefficient can be made arbitrarily close to \(1/6+f<1/3\), a contradiction. [One-letter quarter inequality and pair nonatomic bound](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)

The only remaining case is a split opposite quadratic. Primitivity forces the complete values to be exact opposite squares, whose cyclotomic count is strict below \(1/3\). [Cyclotomic count for opposite square letters](ABC_MINIMAL_REPAIR_CYCLOTOMIC_SPLIT.md)

In particular at every balanced saturating sequence and for each fixed \(q\), the combined \(q\)-defect radical mass is at least \(1/6\) on every pair. The argument fixes \(q\) before both limits. The binary-form extension also handles reducible forms directly, arbitrary actual projective residue lists and ramified root locks. [Binary-form ray amplification](ABC_BINARY_FORM_RAY_AMPLIFICATION.md)

## Minimal target repair

For the least multiplier making \(c\) \(R\)-full, write \(H_R(c)=\prod_{p\mid c}p^{\max(v_p(c),R)}\). The height condition \(H_R(c)\le Y\) implies \(c\le Y\) and \(\operatorname{rad}(c)\le Y^{1/R}\). When the two selected letters are the summands, a nonlinear factor always exists. Apply (PR6) with \(X=Y,\lambda=1/R,\ell=1\). For every sufficiently large fixed \(R\), if the coefficients have \(\operatorname{rad}(AB)\le Y^{f/R}\) and the roots admit a divisor record with \(M\ge Y^\nu,\ \operatorname{rad}(M)\le Y^{\mu/R},\ \gcd(M,2qAB)=1\), the minimally repaired primitive quotient count is
\[
                    \ll Y^{(\mu+1/2+f+\eta)/R+\varepsilon}.
                                                               \tag{PR8}
\]
For \(\mu+f<1/2\) this saves a fixed fraction of \(1/R\). The least multiplier is determined by \(c\), so the count contains no free common-multiplier factor. This direct target bound uses the target repair height and has no separate lower bound on \(\operatorname{rad}(c)\). If a nonlinear factor of degree at least three has a good divisor \(M\ge Y^\nu\) supported on its value, the target-divisor version gives the same bound with \(M\) on that factor instead of on the input roots. In a bin where its radical has mass \(s\) relative to \(Y^{1/R}\), the exponent from the modulus and the target is \((1+s)/(2R)\), up to the bin width and the allowance \(\eta/R\).

For the balanced reverse-transfer sector of primitive total budget at most \((1+2\xi)\lambda\) and lower target cap \(\operatorname{rad}(c)\ge c^{(1/3-\xi)\lambda}\), use (PR7) with enlarged budget \((1+2\xi+\zeta)\lambda\) for fixed \(\zeta>0\). The same \(f\) still bounds the defect mass relative to the enlarged budget. The lower cap gives \(c\le Y^{1/[R(1/3-\xi)\lambda]}\), so the bound has exponent
\[
             \frac{(1/3-\delta(q,f))(1+2\xi+\zeta)}
                  {R(1/3-\xi)}+\varepsilon .
\]
For sufficiently small \(\xi,\zeta\) this is strictly below \(1/R\) up to \(\varepsilon\). It applies to any chosen pair of letters and to the reducible cases above. [Minimal-target-repair equivalence](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)
