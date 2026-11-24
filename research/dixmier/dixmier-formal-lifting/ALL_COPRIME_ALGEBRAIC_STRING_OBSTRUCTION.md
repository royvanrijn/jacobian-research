# Coprime algebraic string operators and every prime-order scalar repair

Let \(F\) be a one-variable algebraic function field over a characteristic-zero
field of constants \(C=\ker\delta\), with nonzero derivation \(\delta\). Write
\[
 F[D;\delta],\qquad [D,f]=\delta f
\]
for the ring of finite scalar differential operators. Constants, coefficient
degrees, algebraic covers and differential orders are unrestricted unless a
condition is stated explicitly.

**Theorem CoprimeAlgebraicStringExclusion.** If \(P,Q\in F[D;\delta]\) satisfy
\[
 [P,Q]=\kappa\in C^\times,\qquad n=\operatorname{ord}P\ge2,
 \qquad m=\operatorname{ord}Q\ge2,                         \tag{1}
\]
then \(\gcd(n,m)>1\).

**Corollary PrimeOrderMateNormalForm.** If \(\operatorname{ord}P=p\) is prime
and \([P,Q]=\kappa\ne0\), every finite mate has the form
\[
 Q=A(P)+Q_1,\qquad A\in C[z],\qquad \operatorname{ord}Q_1=1. \tag{2}
\]
After a finite algebraic coefficient extension and an affine change of the
differential generator, one can write
\[
 \begin{split}
 P&=E^p+\sum_{i=2}^{p-1}c_iE^{p-i}+w,\\
 Q&=A(P)+aE+b,\qquad a,b,c_i\text{ constant},\quad a\ne0,
 \qquad w'=-\kappa/a .
 \end{split}                                             \tag{3}
\]
The sum is empty when \(p=2\). The constants in this normalized description may belong to the extended
constant field; (2) holds over the original field. Every displayed possibility
occurs. For example \(P=E^p+x\), \(Q=-E+A(P)\), \([E,x]=1\), has commutator one.

These results replace the separate quadratic and cubic order arguments by an
all-prime-order statement. The algebraic input missing from a direct infinity
argument is rationality in the correct time coordinate. A finite residue
identity constructs that coordinate. Scalar Stokes injectivity then removes
algebraic branching. A single weighted leading bracket finishes the proof;
no order-dependent list of Painleve equations or degree cutoff is used.

## 1. Normalization and finite fractional-power subtraction

First extend the constants to an algebraic closure and choose a coefficient
field component. The descent to arbitrary characteristic-zero constants is
justified in Section 5. Adjoin an \(n\)-th root \(\rho\) of the leading
coefficient of \(P\). The change \(E=\rho D\), followed by a scalar shift of
\(E\), puts \(P\) in the form
\[
 P=E^n+u_2E^{n-2}+\cdots+u_n,\qquad [E,f]=f'.              \tag{4}
\]
The new derivation is a nonzero algebraic multiple of the original one. These
changes preserve differential orders and the scalar commutator. They are used
on the generic coefficient field, with no assertion of regularity in an
additional deformation parameter.

The formal pseudodifferential root
\[
 W=P^{1/n}=E+\sum_{j\ge1}w_jE^{-j}
\]
exists uniquely. Its coefficients are differential polynomials over
\(\mathbb Q\) in the \(u_i\): the next coefficient in \(W^n=P\) occurs with
nonzero coefficient \(n\). Formal multiplication uses
\[
 E^i f=\sum_{k\ge0}\binom{i}{k}f^{(k)}E^{i-k},
 \qquad i\in\mathbb Z.                                   \tag{5}
\]
Each coefficient uses only finitely many terms. Define \(B_j=(W^j)_+\).
Since \(W^j\) commutes with \(P\),
\[
 [P,B_j]=[(W^j)_-,P],\qquad
 \operatorname{ord}[P,B_j]\le n-2.                        \tag{6}
\]
If an operator with commutator of order at most \(n-2\) has leading term
\(bE^j\), \(j\ge0\), the coefficient of \(E^{n+j-1}\) in that commutator is
\(nb'\). Hence \(b\) is constant. Subtracting \(bB_j\) lowers the order and
preserves the commutator bound. Repeating finitely many times gives
\[
 Q=\sum_{j=0}^{m}c_jB_j,\qquad c_j\in C,\qquad c_m\ne0.    \tag{7}
\]
Thus every constant lower flow is retained. Terms with \(n\mid j\) are
polynomials in \(P\) and commute with it.

## 2. The commutator itself supplies an algebraic time

An arbitrary derivation on an algebraic curve need not admit an algebraic
primitive of one. That property must be proved here.

Put \(B=QW^{-1}=\sum_j b_jE^j\), and write \(P=\sum_{i=0}^n p_iE^i\).
Since \(W^{-1}\) commutes with \(P\),
\[
 [P,B]=\kappa W^{-1},\qquad
 \operatorname{res}_E[P,B]=\kappa.                        \tag{8}
\]
The residue on the left is an exact derivative in the same coefficient field.
A finite primitive is
\[
 \mathcal T(P,B)=
 \sum_{i=1}^{n}\sum_{j=-i}^{-1}\binom{i}{i+j+1}
 \sum_{\ell=0}^{i+j}(-1)^\ell
 p_i^{(\ell)}b_j^{(i+j-\ell)}.                            \tag{9}
\]
Indeed, for \(k=i+j+1\ge1\),
\[
 \binom{j}{k}=(-1)^k\binom{i}{k},\qquad
 ab^{(k)}-(-1)^kba^{(k)}
 =\left(\sum_{\ell=0}^{k-1}(-1)^\ell
 a^{(\ell)}b^{(k-1-\ell)}\right)'.
\]
The residue terms with \(k=0\) cancel, and all other index pairs contribute
zero. This proves \(\mathcal T'=\operatorname{res}[P,B]\). Consequently
\[
 x=\mathcal T(P,QW^{-1})/\kappa\in F,\qquad x'=1.          \tag{10}
\]
Although \(W^{-1}\) is formal, (9) uses only its finitely many coefficients
needed for \(b_{-1},\ldots,b_{-n}\). No formal operator is being declared
finite. The nonconstant \(x\) makes \(F/C(x)\) finite, and the derivation is
the unique extension of \(d/dx\).

## 3. Scalar Stokes injectivity removes algebraic branching

Assume \(\gcd(n,m)=1\). The primary analytic input is
[Beals--Sattinger, *Integrable Systems and Isomonodromy Deformations*,
arXiv:solv-int/9801010v1](https://arxiv.org/abs/solv-int/9801010),
Section 9, printed pages 43--47. The page-43 remark permits precisely the
finite constant linear combinations in (7). Theorem 9.1 and Lemma 9.3 give
the forward sectorial solutions, their phase and constant Stokes matrices.
Theorem 9.4 is scalar injectivity. No inverse-data existence or genericity
assertion from Theorem 9.5 is used.

The compatible scalar equations, with the sign convention (1), are
\[
 P\psi=\lambda\psi,\qquad
 \partial_\lambda\psi=Q\psi/\kappa.                       \tag{11}
\]
Reducing derivatives modulo \(P-\lambda\) gives a rank-\(n\) companion system
whose spectral coefficient matrix is polynomial in \(\lambda\). After
\(\lambda=z^n\), the original companion spectral system is ordinary at every
finite \(z\), including zero.

Put \(J=\operatorname{diag}(1,\omega,\ldots,\omega^{n-1})\), where \(\omega\)
is a primitive \(n\)-th root of unity, and let \(\Lambda\) be the associated
Vandermonde matrix. The normalized sectorial fundamental matrices have the form
\[
 \begin{split}
 \Psi_\nu&=d(z)\Lambda(I+O(z^{-1}))e^{\Phi(x,z)},\\
 d(z)&=\operatorname{diag}(1,z,\ldots,z^{n-1}),\qquad
 \Phi(x,z)=\Phi_0(z)+xzJ .
 \end{split}                                             \tag{12}
\]
The leading spectral eigenvalue coefficient is a nonzero constant multiple of
\(J^m\). Its eigenvalues are distinct because \(\gcd(n,m)=1\); lower flows
do not change them. The finite polynomial and logarithmic coefficients of
\(\Phi_0\) are constants. One can also see their constancy by simultaneously
comparing the diagonal parts of the spectral equation and the \(x\)-equation,
whose sheared form is \(zJ+O(z^{-1})\). Formal diagonalization divides by
differences of the fixed distinct leading eigenvalues, so the finite phase
coefficients belong to the algebraic coefficient field. Their \(x\)-derivatives
vanish, and hence they belong to \(C\).

The right jumps \(\Psi_{\nu+1}=\Psi_\nu S_\nu\) are constant in \(x\).
For two coefficient branches with identical phase and jumps, take the left
ratio \(G=\widetilde\Psi_\nu\Psi_\nu^{-1}\). Equal right jumps cancel. The
ordinary finite spectral system makes \(G\) entire, and (12) gives
\[
 G_{ij}=\delta_{ij}+O(z^{i-j-1}),\qquad 0\le i,j<n.        \tag{13}
\]
The sectorial bounds cover infinity. Therefore \(G\) is a unipotent lower
triangular polynomial matrix and its first row is exactly \((1,0,\ldots,0)\).
The scalar first rows of the two fundamental matrices agree. Their successive
\(x\)-derivatives give all companion rows, so all coefficients of \(P\) agree.
This is scalar injectivity, with the ratio orientation fixed by the right-jump
convention. Matching the phase is essential: translating \(x\) can change a
phase constant while retaining the Stokes matrices.

Continue the algebraic coefficients \(u_i\) around any finite \(x\)-value,
along a loop avoiding their finitely many singular values. The formal phase
constants and all Stokes matrices remain unchanged. The continued and original
branches have a common local domain; the preceding argument identifies them.
This also applies to loops around pole branches. Thus
\(C(x,u_2,\ldots,u_n)/C(x)\) is unramified at every finite point. If its degree
is \(d\), genus is \(g\), and it has \(b\) points over infinity,
[Riemann--Hurwitz](https://stacks.math.columbia.edu/tag/0C1B) gives
\(2g-2=-d-b\). Hence \(d=b=1\), and
\[
 u_2,\ldots,u_n\in C(x).                                 \tag{14}
\]
Only existing algebraic branches were continued. There is no need to reconstruct
a solution from arbitrary Stokes data or control a degenerate inverse problem.

## 4. One weighted bracket excludes every coprime pair

If every \(u_i\) is bounded at rational infinity, all coefficients of every
\(B_j\) in (7) are bounded there. The coefficients of \(Q\) are therefore
bounded too. In the commutator every surviving term differentiates at least
one rational coefficient, so every coefficient tends to zero at infinity.
This contradicts the nonzero constant in (1).

Otherwise set
\[
 R=\max_{2\le i\le n}\frac{\deg_\infty u_i}{i}>0.
\]
Rationality implies \(R\ge1/n\). Assign weight \(R\) to the differential
symbol \(\xi\) and use rational degree at infinity for coefficients. The
leading symbols are
\[
 P_*=x^{Rn}f(z),\qquad Q_*=c_mx^{Rm}g(z),\qquad z=\xi x^{-R},
 \qquad g=(f^{m/n})_+,                                   \tag{15}
\]
where \(f=z^n+\sum_i u_{i,*}z^{n-i}\) is monic, has zero \(z^{n-1}\)
coefficient, and is not \(z^n\).

To justify (15), give \(u_i^{(k)}\) formal weight \(i+k\). Fractional-power
coefficients are homogeneous in these weights. A monomial containing a total
of \(k\) derivatives has rational degree at most its derivative-free weighted
degree minus \((R+1)k\). Thus derivative corrections are strictly lower when
\(R>0\). The lower flows \(j<m\) have strictly smaller weight than the top flow.
The leading fractional power is consequently the ordinary commutative power
in (15). Rational powers of \(x\) only express a weighted symbol; they can
equally be interpreted after a finite Puiseux extension at infinity.

The first-order term of the Ore commutator is the leading Poisson bracket.
Every further derivative lowers weight by another \(R+1\). Direct calculation
gives
\[
 \{P_*,Q_*\}=
 c_mR\,x^{R(n+m-1)-1}(mf'g-nfg').                         \tag{16}
\]
Primes inside the parentheses mean differentiation in \(z\). The weight
\(R(n+m-1)-1\ge(m-1)/n\) is positive. Equality (1) forces the polynomial
factor in (16) to vanish. It follows that \(g^n/f^m\) is constant; monicity
makes that constant one. Coprimeness then gives
\[
 f=h^n,\qquad g=h^m,\qquad \deg h=1.
\]
The missing \(z^{n-1}\) coefficient forces \(h=z\), contradicting \(f\ne z^n\).
This proves CoprimeAlgebraicStringExclusion.

The same algebra records the exact boundary without coprimeness: if
\(d=\gcd(n,m)\), a vanishing factor in (16) is equivalent to
\(f=h^{n/d}\), \(g=h^{m/d}\), with monic \(h\) of degree \(d\). This is a
leading-symbol condition, not a factorization of the differential operators.
It does not settle the composite common-divisor cases.

## 5. Constant-field descent and the prime-order conclusion

For arbitrary characteristic-zero constants, first use the finite algebraic
normalization above and (10). Choose a finitely generated constant subfield
\(C_0\) supporting the finitely many algebraic equations and rational
expressions for the coefficients, normalization data and nonzero \(\kappa\).
Adjoining those coefficients to \(C_0(x)\) gives a finite extension \(F_0\).
Differentiation of their separable minimal equations makes it stable under
\(d/dx\). Embed \(C_0\) in \(\mathbb C\), keeping \(x\) transcendental, and
choose a field component of \(F_0\otimes_{C_0(x)}\mathbb C(x)\). The map from
\(F_0\) is injective and the derivations agree. The complex argument applies,
so a hypothetical coprime pair over the original constants is impossible.

Now let \(\operatorname{ord}P=p\) be prime, before normalization. If
\(\operatorname{ord}Q=kp\), comparison of the leading commutator coefficient
shows that
\[
 c=\operatorname{lc}(Q)/\operatorname{lc}(P)^k\in C.
\]
Subtract \(cP^k\). Repetition terminates with a remainder of order not
divisible by \(p\). An order at least two is coprime to \(p\), and has just
been excluded. Order zero is impossible: after normalization its coefficient
has zero derivative by the coefficient of \(E^{p-1}\), so it commutes with
\(P\). The remainder has order one. This proves (2) over the original \(C\).
In normalized coordinates (7) makes that remainder \(aE+b\) with \(a,b\)
constant. Comparing coefficients of \([P,aE+b]\) gives (3).

## 6. A specialization obstruction without a unit leading coefficient

**Theorem PrimeOrderSpecializationObstruction.** Let \(\mathcal O\) be a
discrete valuation ring in such an algebraic differential field, with
\(\delta\mathcal O\subseteq\mathcal O\). Suppose \(P,Q\in\mathcal O[D;\delta]\)
have nonzero constant commutator, prime generic order for \(P\), and reductions
of differential orders \(s\ge2\) and \(t\ge2\), respectively. Assume the
valuation has a constant uniformizer, as in an \(h\)-adic coefficient ring.
Then \(s\mid t\). The leading coefficient of \(P\) may vanish on reduction.

By (2), absorb the constant term of \(A\) into \(Q_1\) and write
\(Q=\sum_{j\ge1}a_jP^j+Q_1\), with \(\operatorname{ord}Q_1=1\).
If some \(a_j\) had negative valuation, multiply by a constant power of the
uniformizer that makes their minimum valuation zero. Reduction produces a
nonzero polynomial of positive degree in \(\overline P\). Its differential
order is at least \(s\ge2\); it cannot be cancelled by the reduction of an
order-at-most-one operator. The latter reduction is legitimate because
\(Q_1=Q-\sum a_jP^j\) has no coefficient with worse valuation than the minimum
just chosen. This contradicts the vanishing reduction of the scaled \(Q\).
Thus every \(a_j\) is regular, and so is \(Q_1\).

After reduction, either a nonconstant polynomial in \(\overline P\) survives,
giving \(t=js\) for its largest nonzero index \(j\), or only the order-at-most-one
remainder survives. The latter contradicts \(t\ge2\). This proves the theorem.

**Corollary CompositeLongScalarOrders.** In the
[retained exact-first-image scalar Long chart](QUINTIC_PHASE_AND_BRANCH_OBSTRUCTION.md),
\[
 [a,v]=h/v,\qquad
 S_0=4a-3ra^2+2v,\qquad
 T_0=(ra^3-a^2-va)/2 ,
\]
the parameters \(r,h\) are differential constants. Work generically in \(r\),
so these classical differential orders are two and three. Any regular finite
scalar repair \(S,T\), with these reductions and \([S,T]=h\), has **composite
generic differential orders with a common factor**. Indeed apply the
specialization theorem in both directions, and then apply
CoprimeAlgebraicStringExclusion. The generic nonzero commutator is permitted
to vanish under \(h=0\) specialization.

Consequently every prime-order repair, of arbitrary coefficient degree and
arbitrary finite partner order, is excluded. Special parameter values such
as \(r=0\) are not used to compute the generic classical orders.

This conclusion uses the one-variable algebraic scalar chart and its regular
classical specialization. It does not normalize a changed first image into
that chart, treat matrix operators, exclude all composite-order repairs, or
resolve either rank-one or rank-two Dixmier.

## Evidence and scope

The [source ledger](all_coprime_algebraic_string_source_ledger.json) records
the actual primary pages and hypotheses, the pinned PDF and the classical
Long chart. The [checker](verify_all_coprime_algebraic_string.py) controls the
finite residue primitive, fractional-power subtraction, spectral companion
polynomiality, weighted leading bracket, common-power boundary and positive
prime-order examples. Its receipt distinguishes those finite identities from
the universal written proof and the analytic input.

No finite sample proves the theorem for all orders. No historical-priority,
external-review or formal-verification claim is made.

