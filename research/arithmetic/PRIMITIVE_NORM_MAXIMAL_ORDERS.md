# Primitive norm sums in every maximal imaginary quadratic order

Fix an imaginary quadratic field \(K\), fundamental discriminant \(D<0\),
maximal order \(\mathcal O_K\), ideal class group \(G\) and class number \(h\).
All constants below may depend on the fixed field. Write
\[
 P=\{\operatorname N_{K/\mathbb Q}(x+y\omega)>0:\gcd(x,y)=1\},
 \qquad R_P(N)=\#\{(a,b)\in P^2:a+b=N\}.
\]
The integral basis \(1,\omega\) is fixed; changing it preserves primitivity.
These are ordered pairs of numerical values, with no coordinate or ideal
multiplicity. Let \(\Lambda\) be the positive integers represented
primitively by the principal norm form over every \(\mathbb Z_p\).
Call \(N\ge2\) admissible if the sum of two principal norm forms has a
pair-primitive solution at every finite prime; denote this set by
\(\mathcal A_K\). Put
\[
 J_K(N)=
 \prod_{\substack{q\mid N\\q\nmid2D,\ (D/q)=-1}}\frac{q-1}{q-2},
 \qquad R_\Lambda(N)=\#\{(a,b)\in\Lambda^2:a+b=N\}.
 \tag{M.1}
\]

**Theorem MaximalOrderValueCounts.** There are \(c_K,C_K>0,N_K\) such
that every \(N\ge N_K\) in \(\mathcal A_K\) satisfies
\[
             c_K\frac{N J_K(N)}{\log N}
       \le R_P(N)\le C_K\frac{N J_K(N)}{\log N}.
 \tag{M.2}
\]
There is no restriction on \(G\), including its even part or cyclicity.
Every target outside \(\mathcal A_K\) has count zero. The constants and
threshold are uniform in all permitted valuations of the target, but
are not asserted uniform in varying \(K\).

**Theorem ClassErrorComparison.** If \(2G=0\), then \(P=\Lambda\).
Otherwise, for every fixed \(0<\delta<1/h\), uniformly for \(N\ge3\),
\[
 0\le R_\Lambda(N)-R_P(N)
       \ll_{K,\delta}\frac{N J_K(N)}{(\log N)^{1+\delta}}.
 \tag{M.3}
\]
In particular \(\delta=1/(2h)\) is permissible. On sufficiently large
admissible targets,
\[
 R_P(N)=R_\Lambda(N)\left(1+O_{K,\delta}((\log N)^{-\delta})\right).
 \tag{M.4}
\]
These are unweighted value-count statements, stronger than converting a
coordinate asymptotic using a maximum representation multiplicity.
They do not evaluate an unweighted leading constant or numerical onset.

The proof has two independent parts. Generic prime-factor saturation
makes every hidden class failure deficient in one positive-density
ideal-class prime family. An exact prime-square sieve bounds those
failures on every additive line. Separately, the genus shadow is a
finite collection of quadratic-splitting cores with periodic character
conditions. The existing fixed-period theorem supplies a sharp
unweighted lower bound on a positive core for every admissible target.
No field-specific finite absorber pool is constructed.

## 1. Primitive ideals and the genus shadow

A primitive integral ideal is one not contained in \(p\mathcal O_K\)
for any rational prime \(p\). Unique ideal factorization gives its
possible numerical norms: no inert prime factor, ramified exponents
at most one, and every entire split prime power placed on one of its
two conjugate prime ideals. Conversely these conditions produce
primitive ideals. Denote this numerical norm set by \(\mathcal I\).

For \(n\in\mathcal I\), orient every entire split prime power. Its
possible ideal classes form
\[
 \mathcal S(n)=
 \sum_{p\ {\rm ram}}v_p(n)[\mathfrak p_p]
    +\sum_{\substack{q\mid n\\q\ {\rm split}}}
             \{v_q(n)g_q,-v_q(n)g_q\}\subset G.
 \tag{M.5}
\]
In particular, \(n\in P\) exactly when \(0\in\mathcal S(n)\).
A principal ideal in this set has an integral generator of norm \(n\);
a rational prime cannot divide that generator without dividing the
ideal. Orientations of separate copies of a prime power are not
permitted. Every member of \(\mathcal S(n)\) has the same class modulo
\(2G\), because changing an orientation changes its class by \(2v_q(n)g_q\).

We give the principal-genus criterion directly from norm theorems, so
that its use at ramified primes and two is explicit. Write
\(H_p(t)=(t,D)_p\) for the quadratic Hilbert symbol.

**Lemma HilbertGenusKernel.** For a fractional ideal \(\mathfrak a\), the
characters
\[
 \Theta_p([\mathfrak a])=H_p(\operatorname N\mathfrak a),
       \qquad p\mid2D,
 \tag{M.6}
\]
are well defined on \(G\), and their joint kernel is exactly \(2G\).
Moreover
\[
 \Lambda=\{n\in\mathcal I:
       [I(n)]\in2G\text{ for any orientation }I(n)\}.
 \tag{M.7}
\]

**Proof.** Principal ideal multiplication changes the positive rational
ideal norm by a global element norm, hence by a local norm at every
place. The Hilbert norm criterion and multiplicativity make (M.6)
well defined, quadratic, and trivial on \(2G\).

Conversely let every \(\Theta_p([\mathfrak a])\) be one, take an integral
representative and put \(n=\operatorname N\mathfrak a>0\).
At an unramified inert prime its norm valuation is even, the local
norm condition; split places impose no condition. The hypothesis
handles \(p\mid2D\), and positivity handles infinity.
The Hasse norm theorem for the cyclic extension \(K/\mathbb Q\)
gives \(\alpha\in K^*\) with \(\operatorname N\alpha=n\).
The fractional ideal \(\mathfrak b=\mathfrak a/(\alpha)\) has norm one.
Its unique prime-ideal exponent over an inert or ramified rational
prime is zero; at a split pair the exponents are \(t,-t\).
Thus \(\mathfrak b=\mathfrak j/\overline{\mathfrak j}\) for a
fractional ideal \(\mathfrak j\). Conjugation inverts ideal classes,
so \([\mathfrak a]=[\mathfrak b]=2[\mathfrak j]\).
There is no assertion that \(\alpha\) is integral.

For (M.7), primitive local norms have no inert prime factors and
ramified norm valuation zero or one, and are local norms everywhere.
Apply the same norm-one ideal argument. Conversely, membership in
\(2G\) supplies \(H_p(n)=1\). At an inert prime \(n\) is a unit,
so a local norm preimage is a unit. At a ramified prime its prime-ideal
valuation is \(v_p(n)\in\{0,1\}\), so it is integral and outside
\(p\mathcal O_K=\mathfrak p_p^2\), including at two.
At a split prime use \((1,n)\) in
\(\mathcal O_K\otimes\mathbb Z_p\simeq\mathbb Z_p^2\).
This gives primitive local representatives everywhere. \(\square\)

The exact imported norm statements are Milne, Class Field Theory
v4.03, III Theorem 4.4(d), printed p.113, and VIII Theorem 3.1,
printed p.234. The proof has not replaced integral representation
by a rational norm assertion: the latter identifies only the class
modulo \(2G\), and the separate valuation argument proves local
integrality and primitivity.

## 2. A generic saturation criterion

Choose one representative from every unordered class pair
\(\{g,-g\}\) with \(2g\ne0\), and call the resulting finite set \(T_G\).
Set
\[
 \Omega_g=\{q\text{ odd unramified split}:
             [\mathfrak q],[\overline{\mathfrak q}]=g,-g\},
 \quad m_g=\operatorname{ord}(2g)-1,\quad
 k_g(n)=\#\{q\in\Omega_g:v_q(n)=1\}.
 \tag{M.8}
\]

**Lemma ClassSaturation.** If \(n\in\mathcal I\) and
\(k_g(n)\ge m_g\) for every \(g\in T_G\), then
\(\mathcal S(n)\) is the entire coset of \(2G\) determined by \(n\).
Consequently
\[
 \Lambda\setminus P\subset
       \bigcup_{g\in T_G}\{n\in\mathcal I:k_g(n)\le m_g-1\}.
 \tag{M.9}
\]
If \(2G=0\), the union is empty and \(\Lambda=P\).

**Proof.** Reserve \(m_g\) distinct exponent-one primes from each
\(\Omega_g\) and orient all other entire powers arbitrarily.
The \(\Omega_g\) are disjoint. For the reserved block the available sums
\[
              (m_g-2j)g,\qquad0\le j\le m_g,
\]
are precisely \(m_gg+\langle2g\rangle\).
The subgroups \(\langle2g\rangle\), over all nonambiguous inverse pairs,
generate \(2G\). Their independent sums therefore fill a translate of
\(2G\), which is the already determined genus coset. In that coset the
principal class is attainable exactly when \(n\in\Lambda\).
Ramified classes and every power of split two are part of the residual
class and are corrected in the same argument. \(\square\)

For each such inverse pair the rational-prime density is \(1/h\), not
\(1/(2h)\). The fixed-class prime theorem in Bourgain--Fuchs,
arXiv:1105.3698v2, Section 5 Lemma 2, (5.2), gives by partial summation
\[
 \sum_{\substack{q\le y\\q\in\Omega_g}}\frac1q
                  =\frac1h\log\log y+O_K(1).
 \tag{M.10}
\]
The source counts primes represented by an unoriented form class;
its ambiguous-class factor is two, and the classes used here are
nonambiguous. Its fixed-discriminant specialization applies for all
sufficiently large \(y\); smaller primes change only \(K\)-constants.
This Mertens conclusion uses the actual prime-count estimate, not
merely prime-class existence or Dirichlet density.

## 3. The class-error sieve uses exact valuations

Fix \(g\in T_G\), put \(m=m_g,\Omega=\Omega_g\), and count integers
\(1\le a<N\) for which neither \(a\) nor \(N-a\) has an odd inert
prime factor and \(k_g(a)\le m-1\).
Every pair in \(\Lambda^2\) with its first summand deficient in this
class is included. All ramified, dyadic and genus restrictions may
be dropped for this upper bound; no variable ramified coefficient
or squarefull decomposition is necessary.

For odd inert \(q\), the forbidden residues of \(a\) modulo \(q\) are
\(0,N\), with
\[
                       \nu_N(q)=2-\mathbf1_{q\mid N}.
 \tag{M.11}
\]
For a marked prime \(q\in\Omega\), forbid \(v_q(a)=1\), exactly the
\(q-1\) residues \(q,2q,\ldots,(q-1)q\bmod q^2\).
Its forbidden density is
\[
                         f_\Omega(q)=\frac{q-1}{q^2}.
 \tag{M.12}
\]
Thus \(q^2\mid a\) is allowed, including all higher powers.
At a marked split prime this density does not change when \(q\mid N\)
or when its target valuation increases.

**Lemma MixedPrimeSquareSieve.** For a finite mark set
\(\mathcal M\subset\Omega\cap[1,y]\), the number \(S_{\mathcal M}(X;y)\)
of surviving integers in an interval of length \(X\) satisfies,
with \(y=z^{1/16}\),
\[
 S_{\mathcal M}(X;y)
 \ll_K \frac{XJ_K(N)}{\log y}
            \prod_{q\in\mathcal M}(1-f_\Omega(q))
                  +O_\epsilon(z^{4+\epsilon}).
 \tag{M.13}
\]
The constants are uniform in the mark set, root positions and all
target valuations. With no marked primes the same assertion holds
without choosing any ideal-class prime family.

**Proof.** At inert primes use modulus \(\ell_q=q\) and at marked
primes use \(\ell_q=q^2\). For a squarefree prime index \(d\), let
\(\ell(d)=\prod_{q\mid d}\ell_q\),
\(\nu(d)=\prod_{q\mid d}\nu(q)\) and
\(f(d)=\nu(d)/\ell(d)\). CRT gives
\[
 \#\{\text{all forbidden conditions for }q\mid d\}
                      =Xf(d)+O(\nu(d)).
\]
At inert primes \(\nu(q)\le2\); at marked primes \(\nu(q)=q-1\).
Thus \(\nu(d)\le d\). All these primes are odd and
\(0<f(q)\le2/q<1\).

For squarefree indices \(d\le z\), Selberg's quadratic majorant has
\[
 S_{\mathcal M}\le XQ(\lambda)
 +O\!\left(\sum_{d,e\le z}|\lambda_d\lambda_e|\nu([d,e])\right),
 \quad
 Q(\lambda)=\sum_{d,e\le z}\lambda_d\lambda_e f([d,e]).
 \tag{M.14}
\]
Put \(w(d)=\prod_{q\mid d}f(q)/(1-f(q))\) and
\(G(z)=\sum_{d\le z}\mu^2(d)w(d)\).
For completeness, with \(x_d=\lambda_df(d)\) and
\(y_u=\sum_{u\mid d}x_d\), multiplicativity gives
\[
 Q=\sum_{u\le z}\frac{y_u^2}{w(u)},\qquad
 1=\lambda_1=\sum_{u\le z}\mu(u)y_u.
\]
Cauchy is attained by \(y_u=\mu(u)w(u)/G(z)\).
Finite Möbius inversion supplies the actual admissible weights
\[
 \lambda_d=\mu(d)\prod_{q\mid d}(1-f(q))^{-1}
                    \frac{G_d(z/d)}{G(z)},\qquad
 G_d(t)=\sum_{\substack{e\le t\\(e,d)=1}}\mu^2(e)w(e).
\]
They have \(Q=1/G(z)\) and
\(|\lambda_d|\ll_\epsilon d^\epsilon\), uniformly in \(N,\mathcal M\).
Since \(\nu([d,e])\le\nu(d)\nu(e)\le de\), the endpoint error is
\(O_\epsilon(z^{4+\epsilon})\). This explicitly pays for the
prime-square moduli.

Let
\[
 E_{\mathcal M,N}(y)=
 \prod_{\substack{q\le y\\q\ {\rm odd\ inert}}}
                    (1-\nu_N(q)/q)^{-1}
 \prod_{q\in\mathcal M}(1-f_\Omega(q))^{-1}.
\]
Normalize its positive expansion on squarefree prime indices.
A prime is included with probability \(f(q)\), so
\[
 \mathbb E(\log d)=\sum_{q\le y}f(q)\log q
                   \le2\sum_{q\le y}\frac{\log q}{q}\le4\log y.
\]
The last inequality follows from the prime contributions to
\(\log(\lfloor y\rfloor!)\) and
\(\lfloor n/q\rfloor\ge n/(2q)\).
With \(y=z^{1/16}\), Markov retains at least three quarters of
the Euler mass on prime indices \(d\le z\). Hence
\(G(z)\ge(3/4)E_{\mathcal M,N}(y)\).

The inert target correction is exact:
\[
 \prod_{\substack{q\le y\\q\ {\rm odd\ inert}}}(1-\nu_N(q)/q)^{-1}
 =
 \frac{\prod_{q\le y,\ q\ {\rm odd\ inert}}(1-2/q)^{-1}}
      {\prod_{q\le y,\ q\ {\rm odd\ inert},\ q\mid N}(q-1)/(q-2)}.
 \tag{M.15}
\]
Fixed quadratic-character Mertens gives the numerator
\(\asymp_K\log y\); the denominator is at most \(J_K(N)\).
Combine this with (M.14) and the lower bound for \(G\) to prove
(M.13). No arithmetic-progression distribution estimate for the
moving mixed modulus is used. \(\square\)

Independently mark each \(q\in\Omega\cap[1,y]\) with probability
\(t\in(0,1)\). This is a finite probability space, not a hypothesis
about random primes. A deficient integer has at most \(m-1\) exact-one
\(\Omega\)-primes in total; it survives with probability at least
\((1-t)^{m-1}\). Thus its count is at most
\((1-t)^{-(m-1)}\mathbb E S_{\mathcal M}\).
Independence and (M.10) give
\[
 \mathbb E\prod_{q\in\mathcal M}(1-f_\Omega(q))
 =\prod_{\substack{q\le y\\q\in\Omega}}(1-tf_\Omega(q))
                   \asymp_{K,t}(\log y)^{-t/h}.
 \tag{M.16}
\]
The correction \(f_\Omega(q)-1/q=-1/q^2\) is summable.
Choose \(z=N^{1/10}\), \(y=N^{1/160}\), and endpoint
\(\epsilon=1/2\). Every deficient family contributes
\[
 \ll_{K,t}\frac{N J_K(N)}{(\log N)^{1+t/h}}+N^{9/20}.
 \tag{M.17}
\]
Charge a pair in \(\Lambda^2\setminus P^2\) to one of its bad
summands and then to an inverse class in (M.9).
There are finitely many such classes; \(J_K(N)\ge1\).
The power error is absorbed for large \(N\), and increasing the
constant covers all \(N\ge3\). Taking \(t=h\delta<1\)
proves (M.3). When \(2G=0\), (M.7) proves the zero-error assertion.

## 4. Explicit finite characters and sufficient core strata

Let \(R\) be the product of the odd ramified primes, and put
\[
                         L=8R,\qquad S=2R.
 \tag{M.18}
\]
Write \(D=2^\beta v\), \(v\) odd, with \(\beta\in\{0,2,3\}\).
For odd \(a\) coprime to \(D\), define
\[
 \psi_p(a)=\left(\frac ap\right)\quad(p\mid D,\ p\ {\rm odd}),\qquad
 \psi_2(a)=(-1)^{\epsilon(a)\epsilon(v)+\beta\omega(a)},
 \tag{M.19}
\]
where \(\epsilon(u)=(u-1)/2\bmod2\) and
\(\omega(u)=(u^2-1)/8\bmod2\).
These are the unit characters \(H_p(a)\), periodic modulo \(L\).
Milne's dyadic Hilbert formula is in VIII, printed p.248.

For a squarefree prefix \(r\mid S\), set
\(\alpha_r=(H_p(r))_{p\mid2D}\).
At an odd ramified prime, for \(e=v_p(r)\), its entry is
\[
 H_p(r)=
       \left(\frac{(r/p^e)(-D/p)^e}{p}\right).
\]
At two, writing \(r=2^e u\) with \(u\) odd gives
\[
 H_2(r)=(-1)^{\epsilon(u)\epsilon(v)+e\omega(v)+\beta\omega(u)}.
 \tag{M.20}
\]
Let \(\boldsymbol\psi\) be the full vector in (M.19),
\(\chi=\chi_D\bmod L\), \(H=\ker\chi\), and \(\mathcal C\)
the semigroup supported on odd unramified split primes.
Hilbert reciprocity gives
\[
 \chi(a)=\prod_{p\mid2D}\psi_p(a)\quad((a,L)=1,\ a>0),
 \qquad \prod_{p\mid2D}H_p(r)=1\quad(r\mid S).
 \tag{M.21}
\]
Outside \(2D\) the first product's omitted factors are
\(\chi(q)^{v_q(a)}\), and the infinite factor is one.
For \(r\) supported on \(S\), every place outside \(2D\)
already has symbol one. These identities supply compatibility
with \(H\), rather than assuming it.

Let \(\mathcal T\) consist of squarefree \(r\mid S\), excluding
even \(r\) when two is inert. For each such prefix put
\[
 A_r=\{a\in(\mathbb Z/L)^*:\boldsymbol\psi(a)=\alpha_r\}
        \subset H,\qquad
 \mathcal C_r=\{c\in\mathcal C:\boldsymbol\psi(c)=\alpha_r\}.
\]
Then
\[
                              r\mathcal C_r\subset\Lambda.
 \tag{M.22}
\]
Indeed \(rc\) satisfies the prime-valuation conditions and every
Hilbert norm symbol is one. Every \(A_r\) is nonempty: an ideal
of norm \(r\) exists, and an odd degree-one class prime away from
\(2D\) in its class supplies that character vector.
If two splits, (M.22) is only a sufficient family; the full
shadow may have arbitrarily high powers of two.

## 5. Finite-vector fixed-period transfer and positive unweighting

We state the analytic interface with all states retained.
It uses the quadratic case of
[Asymmetric subgroup splitting](ASYMMETRIC_SUBGROUP_SPLITTING.md),
its fixed-period box input, and the explicit genus-sign transfer
in the [class-number-two](PRIMITIVE_NORM_CLASS_NUMBER_TWO.md) and
[class-number-four](PRIMITIVE_NORM_CLASS_NUMBER_FOUR.md) proofs.

**Lemma FiniteGenusFilter.** Fix coprime positive integers \(A,B\)
supported on primes dividing \(L\). Let
\(\mathcal C_\alpha=\{c\in\mathcal C:\boldsymbol\psi(c)=\alpha\}\).
For every sufficiently large \(M\) with
\[
 \Gamma_{A,B;\alpha,\beta}(M)=
 \frac{L}{|H|^2}\#\{(a,b)\in H^2:
        Aa+Bb=M\bmod L,\
        \boldsymbol\psi(a)=\alpha,\boldsymbol\psi(b)=\beta\}>0
 \tag{M.23}
\]
one has
\[
 \#\{Ax+By=M:x\in\mathcal C_\alpha,\ y\in\mathcal C_\beta\}
               \asymp_{K,A,B}\frac{M J_K(M)}{\log M}.
 \tag{M.24}
\]
This is uniform in target valuations. Zero finite density gives
an identically zero count.

**Proof of the weighted interface.** Set \(V=\tau\mathbf1_\mathcal C\)
and \(F_T(x)=V(Tx)/V(T)\). Retain every feasible vector state
throughout the family. For good shifts \(T_1,T_2\) and good
coefficients \(u,v\in\mathcal C\), the required weighted family is
\[
 \begin{split}
 &\sum_{\substack{Au x+Bv y=M\\
       \boldsymbol\psi(x)=\alpha,\ \boldsymbol\psi(y)=\beta}}
       F_{T_1}(x)F_{T_2}(y)\\
 &\quad=\frac{M}{ABuv}\,
       \Gamma_{A,B;u,v;\alpha,\beta}(M)\rho_{T_1}\rho_{T_2}
       \prod_{q\nmid L}\frac{D_q}{M_{1,q}M_{2,q}}
            +O(M^{1-\eta}).
 \end{split}                                                    \tag{M.25}
\]
Here \(\Gamma\) is (M.23) with coefficients \(Au,Bv\);
the good local radial line is \(AuX+BvY=M\).
The residue means and good factors are exactly (H3)--(H4)
of the subgroup theorem. Shifts, coefficients and all tested
divisor products are within a sufficiently small power \(M^\theta\),
chosen after \(\eta>0\), as in (H17)--(H18).

Insert on the original first unit variable the filter
\[
 \mathbf1_{\boldsymbol\psi(x)=\alpha}
 \mathbf1_{B\mid M-Aux}
 \mathbf1_{(z,L)=1,\ \chi(z)=1,\
              \boldsymbol\psi(z)=\boldsymbol\psi(v)\beta},
                  \qquad z=(M-Aux)/B.                       \tag{M.26}
\]
It has fixed period \(\ell=BL\). Expand it in the finitely many
characters \(\xi\bmod\ell\). For \(\xi\in\{1,\chi\}\),
\(F_T\xi=F_T\) pointwise. For every other \(\xi\), both
\(\xi,\xi\chi\) are nonprincipal, and the source bound (H9)
gives zero principal mean and a power-saving bound.
Kernel indices need only be \(\ell\)-units of \(\chi\)-value one.
They need not have trivial genus vector, and all kernel phases
\(\xi(jc_0)\) are retained.

The filter guarantees the second paired-divisor quotient has
\(\chi=1\), as required by (H11)--(H18).
The congruence \(Aux=M\bmod Bvka\) has its entire \(B\)-part
in (M.26); the remaining varying modulus \(vka\) is an
\(\ell\)-unit. The true hyperbola boundary is
\(Aux+Bvka^2=M\). Fixed-period boxes apply at period \(\ell\).

There is one additional coefficient check beyond the earlier
\((1,t)\) proofs. In the medium-variable BGLM model the parameters
are \(a_1=M,\ a_2=Aujc_0,\ a_3=1,\ d=vka_{\rm aux},\ c=\ell\).
The actual arXiv:2111.01601v5 Proposition 5.1, printed p.15,
does not require \((a_2,c)=1\). Its varying \(q\) and short
variable \(n\) are \(\ell\)-units, so the extra factor \(A\)
preserves \((q,a_1a_2d)=1\) and \((n,a_2)=1\).
The remaining coprimalities are the source's original extractions.
Since \(A,B,\ell\) are fixed, decreasing \(\theta\) and taking
the threshold larger preserves the size hypotheses.
All old small-summand, smooth, kernel, gcd and bilinear errors
therefore fit the same (H17)--(H18) hierarchy. The finite vector
adds only finitely many bounded character coefficients.

For normalization let \(H_\ell\) be the lift of \(H\);
\(|H_\ell|=B|H|\). If \(\mu\) is the average of (M.26) on
\(H_\ell\), each surviving first residue gives \(B\) second lifts.
Thus
\[
 \frac{\ell}{|H_\ell|^2}\#\{\text{finite pairs modulo }\ell\}
             =\frac{2BL}{\varphi(L)}\mu.
\]
Each root modulo \(L\) has exactly \(B\) pair lifts modulo \(BL\):
at a prime of \(B\), \(Au\) is a unit. This equals the
\(\Gamma\) in (M.23). The first real interval contributes \(1/A\);
the raw fixed factor is \(2(\varphi(L)/L)\mu\).
Using (H10) gives precisely
\(M\Gamma\rho_{T_1}\rho_{T_2}/(ABuv)\), with no missing \(B\).

Positive \(\Gamma\) forces \((M,AB)=1\). Every source
coefficient/gcd reduction now carries the vector state:
\(x=e x'\) changes \(\alpha\) to
\(\alpha\boldsymbol\psi(e)\), and similarly on the other side.
Common good factors divide both coefficients and the target;
cross-forced factors and target gcds use the same exact
substitutions as (H15)--(H16). All removed factors are units
modulo \(\ell\). Multiplication by them bijects the finite
solution pairs, so their correctly transported \(\Gamma\)
is unchanged. Holding the original states fixed after removal
would be an invalid substitution.

**Positive unweighting.** Apply the positive section-6 proof of
the subgroup theorem to the full family (M.25), not just its
base weighted sum. For actual restored values
\(x=d u x',y=d v y'\), retain the residual states
\[
 \alpha\boldsymbol\psi(du),\quad
 \beta\boldsymbol\psi(dv),\qquad
               Au x'+Bv y'=M/d.
\]
The unit map \((x',y')\mapsto(du x',dv y')\) bijects the
finite pairs with the original ones. Hence the same positive
\(\Gamma\) factors out of every lower-model term.
Primes with nontrivial genus vector are permitted; discarding
them would lose the required restoration exponent.

At split good primes the rough weighted bad density is
\(4/q+O(q^{-2})\); their density is \(1/2\), giving rough
dimension two. The two directed valuation-one restoration
options each have ratio \((q-1)/(q(q-2))\), so together they
restore exponent one over all split primes. The target inert
factor left over is exactly \((q-1)/(q-2)\).
Fixed coefficients are units at every such prime.
The source's recoverable patterns, injection into actual
values, residual weight bound, positive cofactor model and
Markov cap are unchanged; filtering only decreases their
upper bounds. The uniformly shifted estimates (M.25)
supply the same summed error budget.
The finitely many positive \(\Gamma\)'s have a positive
minimum. This proves the lower bound (M.24).
Dropping the vector states and sieving only good inert
primes proves the matching dimension-one upper bound.
This is a finite-state extension of the existing proof;
it does not infer unweighted counts by dividing a weighted
asymptotic by divisor multiplicities. \(\square\)

## 6. Finite profiles cover every local target

For \(r,s\in\mathcal T\), put
\[
 c=(r,s),\quad A=r/c,\quad B=s/c,\quad M=N/c.
\]
If \(c\nmid N\), this profile is unavailable. Otherwise define
\[
 \Gamma_{r,s}(N)=\frac{L}{|H|^2}
       \#\{(a,b)\in A_r\times A_s:
                            Aa+Bb=M\bmod L\}.             \tag{M.27}
\]
The labels are those of the original prefixes \(r,s\),
not the reduced \(A,B\).

**Lemma ShadowProfileCoverage.** For every \(N\ge2\),
\[
 N\in\mathcal A_K
 \quad\Longleftrightarrow\quad
 \text{some }r,s\in\mathcal T\text{ have }c\mid N,\
                                  \Gamma_{r,s}(N)>0.
 \tag{M.28}
\]
In particular this is a finite local test modulo \(16R^2\).

**Proof.** Choose actual primitive local norm summands at all
primes dividing \(2D\). At ramified primes their valuations
are zero or one. At inert two both are odd.
At split two we may choose both odd if \(N\) is even
(use norm values \(1,N-1\)), and valuations zero and one
if \(N\) is odd (use \(N-2,2\)). Every value in the split
algebra is a primitive norm. These choices select
\(r,s\in\mathcal T\).

The quotients of the local norm summands by their full
prefixes are units at each bad prime, and have the respective
Hilbert characters \(\alpha_r,\alpha_s\).
Their equations divided by \(c\) give \(Aa+Bb=N/c\).
CRT on the unit residues modulo \(L\) supplies the pair in
(M.27). Compatibility with \(H\) follows explicitly from
(M.21); no additional splitting-residue condition is assumed.

Conversely, at each bad prime a residue root in (M.27)
lifts to an exact linear solution in the same unit residue
classes, since \((A,B)=1\). Its Hilbert characters and
primitive valuations remain fixed, so multiplying back
by \(r,s\) gives primitive local norm values summing to \(N\).
At every odd prime outside \(2D\), every target is a sum
of two primitive binary norms: inert unit norms cover
all targets by two units, and the split algebra imposes
no restriction. Two was already treated. This proves
(M.28). Finally \(c\mid2R\), so \(N\bmod16R^2=(2R)L\)
determines both \(c\mid N\) and \(N/c\bmod L\).
Every positive finite density is at least \(L/|H|^2\).
Arbitrarily deep allowed target cancellation is retained.
\(\square\)

**Theorem GenusShadowCounts.** For every sufficiently large
\(N\in\mathcal A_K\),
\[
                    R_\Lambda(N)\asymp_K
                              \frac{N J_K(N)}{\log N}.     \tag{M.29}
\]

**Proof.** Select a positive profile in (M.28). Apply (M.24)
to its \(A,B,\alpha_r,\alpha_s\). The injective numerical map
\((x,y)\mapsto(rx,sy)\) gives shadow values summing to \(N\).
All prefixes belong to a fixed finite set,
\(M=N/c\asymp_K N\), and \(J_K(M)=J_K(N)\).
Take the minimum lower constant and maximum threshold over
the positive finite family.

For the upper bound do not decompose all possible powers
of split two. Drop genus, ramified and dyadic restrictions,
and sieve the two original values only at good odd inert
primes. Their two roots merge exactly at target primes,
giving \(O_K(NJ_K(N)/\log N)\).
Outside \(\mathcal A_K\) the shadow convolution is exactly
zero because its values themselves have primitive local
representatives. \(\square\)

Combining (M.3) with (M.29) proves (M.4) and the lower
bound (M.2); the same upper sieve proves its upper bound.
This completes the two theorems stated at the start.

## 7. Sources, finite evidence and precise limits

The arithmetic strategy of many class primes and independent
conjugation is classical. Bredikhin--Linnik 1966, introduction
pp.145--147 and Lemma 5 pp.157--160, explicitly uses this
class-mixing mechanism. Its formal A--C results concern
prime plus an ideal norm, with ideal/coordinate weights.
Bourgain--Fuchs v2 Theorems 2'/3' concern squarefree
within-genus exceptions and fixed translated primes.
Those inspected statements are not the all-target
unweighted two-norm comparison (M.3). This note supplies
the required first-prime-power sieve and finite-vector
transfer explicitly; no historical-priority claim follows.

Primary editions and exact application scopes are pinned in
[the source ledger](primitive_norm_maximal_order_source_ledger.json).
The genus-filter proof uses the actual BGLM v5 Proposition 5.1
hypotheses, together with the existing fixed-period
quadratic splitting and positive restoration proofs.
The class-prime Mertens input is genuinely a fixed ideal-class
prime-distribution statement, stronger than choosing finitely
many class primes. It can involve nonabelian splitting
data over \(\mathbb Q\); the genus-shadow additive estimate
itself uses quadratic Dirichlet characters.

[The finite checker](verify_primitive_norm_maximal_orders.py)
audits generic sign saturation in its stated finite groups,
mixed prime/prime-square CRT counts, independent rational
Selberg optimization, finite probabilistic majorants and
purposeful corruptions. Its
[retained receipt](primitive_norm_maximal_order_checks.json)
specifies the exact finite domains. These checks do not
prove the imported prime-class, norm or infinite analytic
theorems, and are not external or formal verification.

Certificate verification requires active Python assertions. The checker rejects optimized execution (-O, -OO, or optimization enabled through PYTHONOPTIMIZE) before generating a receipt.

Run from the repository root:

    python3 research/arithmetic/verify_primitive_norm_maximal_orders.py --output research/arithmetic/primitive_norm_maximal_order_checks.json

This is an internally written and source-audited theorem for
every fixed maximal order. No conductor extension, mixed-field
theorem, exact global exception list, leading constant,
numerical cutoff or uniformity in varying discriminant is
asserted. At conductor primes primitive order norms can
have unbounded valuations and require a different
ideal/local interface; none is silently imported here.
