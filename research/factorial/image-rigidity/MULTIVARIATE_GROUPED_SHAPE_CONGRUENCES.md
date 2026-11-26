# Coefficientwise shape congruences from several Newton slopes

The projective polynomial of a possible ordinary factorial common factor admits a separate denominator bound for each coefficient. At primes dividing these bounds, the scaled factors in the individual Newton slope groups determine additional residue classes. This extends the [split quadratic formula](MULTIVARIATE_DENOMINATOR_PRIME_SHAPE_CONGRUENCES.md) to higher factor degrees and to nonintegral root valuations. It preserves the condition that the same local allocation must account for every coefficient.

## The local formula

Let \(P_n,P_N\) be adjacent ordinary factorial parents, and suppose a possible monic rational common factor \(H\) has degree \(d\). Its roots \(\alpha_1,\ldots,\alpha_d\) are counted with multiplicity, and their product \(c\) is nonzero. Write
\[
          \prod_{i=1}^d(1+\gamma_iY)=\sum_{k=0}^d E_kY^k,
          \qquad \gamma_i=\alpha_i^d/c .
\]
The \(E_k\) are rational and \(E_0=E_d=1\). This is the same projective shape as the monic \(\prod_i(Y-\gamma_i)\), with the coefficients read in the opposite order and alternating signs. [Projective denominator rule](MULTIVARIATE_DIFFUSE_PROJECTIVE_FACTOR_SIEVE.md#exact-projective-denominators)

Fix a prime \(p\) and a possible common-root valuation allocation. Let its distinct root valuations be \(s_g=a_g/e\) in a common denominator \(e>0\), with multiplicities \(m_g\). Use the Eisenstein field \(K=\mathbb Q_p(\pi)\), where \(\pi^e=p\), normalize \(w(\pi)=1\), and let \(\mathcal O\) be its valuation ring. Define the monic unit-root polynomial
\[
       H_g(X)=\prod_{v_p(\alpha_i)=s_g}(X-\alpha_i/\pi^{a_g})
                         \in\mathcal O[X].
\]
The root group is invariant over \(K\), and its normalized roots have valuation zero. These facts give the assertion about the coefficients and the unit constant term. Set
\[
 c_g=(-1)^{m_g}H_g(0),\quad U=\prod_gc_g,\qquad
 b_g=d a_g-\sum_hm_h a_h .
\]
For each group put
\[
 R_g(Y)=\operatorname{Res}_X(H_g(X),1+YX^d)
              =\sum_{\ell=0}^{m_g} C_{g,\ell}Y^\ell .
\]
These coefficients are polynomials over the integers in the coefficients of \(H_g\). The monic resultant formula gives
\[
  E_k=U^{-k}
       \sum_{\substack{0\le\ell_g\le m_g\\\sum_g\ell_g=k}}
         \pi^{\sum_g b_g\ell_g}\prod_g C_{g,\ell_g}.             \tag{1}
\]
Indeed \(c=\pi^{\sum_h m_h a_h}U\), and each normalized root in group \(g\) contributes \(\pi^{b_g}\) times its \(d\)-th power divided by \(U\).

For this allocation define \(S_k\) to be the sum of the \(k\) smallest entries in the list of \(b_g\)'s repeated \(m_g\) times. Equation (1) gives \(w(E_k)\ge S_k\). Its valuation is a multiple of \(e\) because \(E_k\) is rational. Thus the denominator exponent at \(p\) satisfies
\[
              v_p(\operatorname{den} E_k)
                      \le \max(0,\lfloor-S_k/e\rfloor).         \tag{2}
\]
Maximize this bound over all permitted allocations to get \(d_{k,p}\), and put \(D_k=\prod_p p^{d_{k,p}}\). The product is finite since outside the primes of the parent endpoints all root valuations are zero. The same allocation is used for all \(k\) when evaluating (1), even if the separate maxima in the bounds came from different allocations.

## Passing truncated face data to rational residues

Suppose that at the given allocation each \(H_g\) is specified modulo \(\pi^L\), where \(L\ge1\), and let \(T_k=D_k E_k\in\mathbb Z\). Use any unit-constant monic representatives of the local group polynomials in (1). Denote the resulting expression for \(D_kE_k\) by \(\widetilde T_k\in K\). Every \(C_{g,\ell}\) and \(U^{-1}\) changes by an element of \(\pi^L\mathcal O\) when the representatives change. Consequently
\[
                w(T_k-\widetilde T_k)\ge L+e d_{k,p}+S_k.      \tag{3}
\]
When the right side is positive, an integer congruent to \(\widetilde T_k\) at that precision has a uniquely specified residue modulo
\[
                    p^{\left\lceil(L+e d_{k,p}+S_k)/e\right\rceil}.
                                                                 \tag{4}
\]
If there is no such integer, the candidate local data are impossible. This is a rationality condition on the residue over \(K\). The projection is exact: \(\pi^b\mathcal O\cap\mathbb Q_p=p^{\lceil b/e\rceil}\mathbb Z_p\) for integer \(b\). Formula (3) still applies if some summands have negative \(w\)-valuation. In that case the sum with its prescribed precision is formed in \(K\) before projecting to rational integers.

There is an exact finite source for the candidate \(H_g\bmod\pi^L\). In each parent substitute \(z=\pi^{a_g}X\), divide by the minimum coefficient power of \(\pi\), and truncate at precision \(L\). The monic \(H_g\) must divide both scaled parents in \((\mathcal O/\pi^L)[X]\). This gives finite sets of unit-constant monic polynomial candidates at prescribed degrees. The [truncated Newton comparisons](MULTIVARIATE_PHYSICAL_TRUNCATED_NEWTON_FACES.md) can first restrict their admissible degrees. Enumerating the candidate slope allocations and factors gives a union of compatible residue vectors for \((T_1,\ldots,T_{d-1})\) at \(p\). Different coordinates from the same choice are kept together. The construction only needs necessary local conditions and permits further restrictions from the exact local branches.

### Extra agreement at a repeated local block

There is sometimes more precision on the shape than on the individual
group factors. Suppose \(p\mid d\) and \(p\mid k\), and put
\(r=v_p(\gcd(d,k))\). If two surviving choices of the same valuation
allocation agree on all \(H_g\bmod\pi^L\), then their expressions
for \(T_k\) in (1) agree at precision
\[
               \min(L+e,p^rL)+e d_{k,p}+S_k ,                 \tag{5}
\]
in place of the precision in (3).

In characteristic \(p\), the polynomial in the coefficients of
\(H_g\) giving the elementary symmetric function of the \(d\)-th
powers of its roots is a \(p^r\)-th power. Indeed it is the
\(p^r\)-th power of the corresponding polynomial for the
\((d/p^r)\)-th powers, by the formal symmetric-function identity.
Their difference over the integers is \(p\) times an integer
polynomial. When the inputs vary by \(\pi^L\), this last term varies
by at least \(L+e\). The binomial expansion for the difference of
two \(p^r\)-th powers gives the bound \(\min(L+e,p^rL)\) on the
other term. The same bound holds for the change in \(U^{-k}\),
because it is the \(p^r\)-th power of a unit. Multiplying by the
fixed powers of \(\pi\) in (1) gives (5). This can determine a
further residue of a shape coefficient without distinguishing the
underlying local factors.


The same divisibility gives a separate improvement to (2). Write
\(q=p^r\). For \(\sum_g\ell_g=k\), the exponent
\[
                  \sum_g b_g\ell_g
                     =d\sum_g a_g\ell_g-k\sum_h m_h a_h
\]
is divisible by \(q\). Use that exponent divided by \(q\),
replace the \(d\)-th powers by \((d/q)\)-th powers, and replace
\(U^{-k}\) by \(U^{-k/q}\) in (1), forming \(A_k\).
The Frobenius identities give \(E_k=A_k^q+pB_k\) with
\(w(B_k)\ge S_k\). Indeed every pure power and every intervening
multinomial term carries valuation at least \(S_k\), and all the
intervening multinomial coefficients are divisible by \(p\).
If \(w(E_k)<S_k+e\), its leading valuation is therefore a
multiple of \(q\). The only rational valuation in this first
interval is \(v_p(E_k)=\lceil S_k/e\rceil\). If this integer
is not a multiple of \(q/\gcd(q,e)\), the true valuation is
at least \(\lceil S_k/e\rceil+1\), saving one denominator digit
whenever (2) was positive. This refinement requires only the
candidate slope allocation.

For example \(z^4+z^3+4z^2+4z+4\) has, at \(2\), three roots of
valuation \(2/3\) and one unit root. For \(d=4,k=2,e=3\) one has
\(S_2=-4\) and \(q=2\). The initial denominator cap has exponent
one, which drops to zero by the preceding criterion.

## Removing the denominator loss in an interval test

Coefficients at good native moduli give congruences on the same \(T_k\), using [joint shape reconstruction](MULTIVARIATE_JOINT_NATIVE_SHAPE_RECONSTRUCTION.md). At a prime \(p\mid D_k\), the sufficient precision \(L\ge -S_k\) makes (4) at least \(p^{d_{k,p}}\). Suppose the union of permitted choices gives one such residue class at each denominator prime, and the coprime good modulus \(M\) determines another class. Then \(T_k\) lies in one progression of modulus at least \(D_k M\). For any certified real interval \(I_k\) containing \(E_k\), exclusion follows if the progression misses \(D_k I_k\). If \(|I_k|<M\), at most one value of \(E_k\) survives, independently of the size of the coefficient denominator budget. Without uniqueness of the local residue class, the corresponding finite union of progressions still gives an exact necessary set.

An interval can already be obtained from bounds \(0<r_{\min}\le|\alpha_i|\le r_{\max}\). Each of \(\gamma_i\) and its inverse has absolute value at most \((r_{\max}/r_{\min})^d\), and \(\prod_i\gamma_i=1\). Hence
\[
                  |E_k|\le \binom dk
                         (r_{\max}/r_{\min})^{d\min(k,d-k)}.
\]
Tighter intervals may use simultaneous real root boxes. The argument gives no general guarantee that the local groups are determined to the required precision or that the joint residue sets give an empty real intersection.

## Bounding several compatible residue classes

Use a common scale \(D=\operatorname{lcm}_kD_k\) and the integer vector
\(Z=(D E_1,\ldots,D E_{d-1})\). Data on \(T_k=D_k E_k\) modulo
\(p^{e_{k,p}}\) give data on the corresponding \(Z_k\) modulo
\(p^{e_{k,p}+v_p(D/D_k)}\). Thus the correlated joint classes of \(Z\)
are available modulo \(p^{e_p}\), with
\(e_p=\min_k(e_{k,p}+v_p(D/D_k))\). In particular each coordinate
calibrated through its own denominator gives \(e_p\ge v_p(D)\).

Let \(r_{p,j}\ge1\) be the number of jointly allowed classes of \(Z\)
modulo \(p^j\), for \(1\le j\le e_p\). Suppose the real coefficient
intervals have maximum width \(\Delta\), and put \(W=D\Delta\).
If \(W<1\), at most one integer vector remains. For \(W\ge1\) set
\[
 A=\sum_p\sum_{j=1}^{e_p}\frac{\log p}{r_{p,j}},
          \qquad B=\sum_p e_p\log p .
\]
If \(A>\log W\), the number of candidate integer vectors is at most
\((B-\log W)/(A-\log W)\). This is a correlated vector form of
[Gallagher's larger sieve, Theorem 1](https://www.impan.pl/en/publishing-house/journals-and-series/acta-arithmetica/all/18/0/97714/a-larger-sieve).

Indeed, among \(q\) candidate vectors, a level admitting \(r\) classes
contains at least \(q(q/r-1)/2\) pairs in the same class, by
Cauchy–Schwarz. This lower bound may be replaced by the exact minimum
\(r\binom{\lfloor q/r\rfloor}{2}+(q\bmod r)\lfloor q/r\rfloor\).
For each pair of distinct vectors the positive gcd of their
coordinate differences is at most \(W\). Each equal residue
level contributes one factor of \(p\) to that gcd. Multiplying over
the pairs yields \(q(qA-B)\le q(q-1)\log W\), as required.
The exact minimum gives a comparison of integer products and a
rational power of \(W\), without logarithms.

If all denominator prime powers are covered and
\(M_{\rm extra}=\prod_p p^{e_p-v_p(D)}\), the gap in the bound is
\[
   A-\log W = \log M_{\rm extra}-\log\Delta
           -\sum_p\sum_{j=1}^{e_p}(1-1/r_{p,j})\log p .
\]
Thus uniquely calibrated levels contribute their entire depth.
The remaining loss comes only from levels with several allowed
joint classes.

## An ordinary quartic calculation

Take three ordinary rows \((0,1)\) and three \((1,0)\). Its primitive order-four parent, in increasing coefficient order, is \((144,9,4,9,144)\). At \(p=3\) it has two roots of valuation \(1\) and two of valuation \(-1\). In its degree-four projective shape, the list of normalized exponents is \((-4,-4,4,4)\). The separate denominator exponents for \(E_1,E_2,E_3\) are exactly \(4,8,4\), respectively. Clearing all coefficients by the full shape denominator would spend depth eight even on the first and third coordinates.

The [exact checker](verify_multivariate_grouped_shape_congruences.py) recomputes the ordinary quartic from factorials and its projective coefficients by Newton identities. It lifts its degree-two unit-root factor after \(z=3X\), obtains the opposite-slope factor by reciprocity, and recomputes the scaled coefficient residues by (1) through precision eight. It checks extra Frobenius precision on ambiguous group factors and a correlated two-coordinate counting bound with different original coordinate precisions. A separate control uses the real-positive polynomial \((z^4+3)(z^2+1)\) to check nonintegral Newton valuations and the rounding in (2). The quartic is used as a parent calculation and does not assert a common factor of ordinary moments.
