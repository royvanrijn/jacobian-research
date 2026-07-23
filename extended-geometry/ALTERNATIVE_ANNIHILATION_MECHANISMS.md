# Beyond prime reduction: alternative lower-face annihilation mechanisms

## 1. Executive conclusion

The prime in the lower-face proof is not essential as a prime, but its two
jobs must be distinguished.

Starting from
\[
F_r(U)=\operatorname{CT}_T(P^r)
      =cU^n+\sum_{j>n}c_jU^j,\qquad c\ne0,                 \tag{1.1}
\]
the characteristic-\(p\) argument uses
\[
F_{rp}=F_r^p\pmod p                                      \tag{1.2}
\]
to remove cross terms, and then uses
\[
p\mid \frac{(jp)!}{(np)!}\qquad(j>n)                     \tag{1.3}
\]
to kill every higher radial layer after factorial normalization.

Thus a replacement needs an accessible family of operations with:

1. a nonzero image of \(cU^n\);
2. suppression of cross terms created by dilation or powering; and
3. a vanishing ideal, valuation, or asymptotic scale that contains every
   higher layer but not the normalized lowest layer.

There are two genuinely promising replacements.

* **Cyclotomic/q replacement.** A \(q\)-binomial or \(q\)-shuffle
  deformation, reduced modulo \(\Phi_N(q)\), reproduces both (1.2) and
  (1.3), for arbitrary \(N>1\), not only primes.
* **Scale-family replacement.** If moment vanishing is available in a
  one-parameter variance family, then \(R\to\infty\) isolates the lowest
  layer directly. This gives a clean characteristic-zero analytic theorem.

Ordinary roots-of-unity filters and finite-difference operators are useful
projectors, but by themselves they do not follow from the original
fixed-variance moment hypothesis. Function-field valuations are the right
abstract language for a deformation, but a valuation without an accessible
deformed vanishing identity does not prove anything new.

## 2. What the prime proof actually suppresses

The mixed terms inside the exposed face should **not** be killed. They are
what produce
\[
c=\operatorname{CT}_T(P_0^r)\ne0
\]
through the Duistermaat--van der Kallen theorem. The mixed terms killed by
Frobenius are the later cross terms among the radial coefficients
\(c,c_{n+1},\ldots\) when \(F_r\) is raised to the \(p\)-th power.

This distinction rules out several superficially plausible replacements.
An operator that projects angular weight zero merely repeats
\(\operatorname{CT}_T\); it does not isolate the radial coefficient \(c\).
An operator that kills all mixed face monomials may also kill the only
reason \(c\) is nonzero.

The exposed-face step itself remains characteristic-free. See the original
[constant-term theorem](https://webspace.science.uu.nl/~kalle101/powers.pdf)
and the repository's
[lower-face theorem](TWO_REAL_GMC_LOWER_FACE_THEOREM.md).

## 3. An abstract annihilation datum

The arithmetic second half of the proof can be packaged as follows. Let
\(\mathfrak a_N\) be an ideal, or let \(v_N\) be a valuation, on a
coefficient ring \(A_N\). Suppose there are operations \(D_N\) and
normalized functionals \(\widetilde{\mathcal L}_N\) such that
\[
D_N(F_r)
=c^{[N]}U^{nN}+E_N,                              \tag{3.1}
\]
where
\[
\widetilde{\mathcal L}_N(c^{[N]}U^{nN})
\notin\mathfrak a_N,\qquad
\widetilde{\mathcal L}_N(E_N)\in\mathfrak a_N.   \tag{3.2}
\]
If moment vanishing implies
\(\widetilde{\mathcal L}_N(D_N(F_r))=0\), reduction modulo
\(\mathfrak a_N\) gives a contradiction.

The valuation version replaces (3.2) by a unique strict minimum:
\[
v_N\!\left(\widetilde{\mathcal L}_N(c^{[N]}U^{nN})\right)
<
v_N\!\left(\widetilde{\mathcal L}_N(E_N)\right). \tag{3.3}
\]
The analytic version replaces it by
\[
\widetilde{\mathcal L}_R(D_R(F_r))
=\kappa c+o(1),\qquad \kappa\ne0.                \tag{3.4}
\]

This formulation makes the real issue visible: not the existence of a
projector, but the implication from the assumed moment identities to the
vanishing of its output.

## 4. Cyclotomic divisibility: the closest algebraic replacement

Put
\[
[m]_q=\frac{1-q^m}{1-q},\qquad
[m]_q!=\prod_{a=1}^m[a]_q.
\]
For every \(N>1\),
\[
v_{\Phi_N}([m]_q!)=\left\lfloor\frac mN\right\rfloor.       \tag{4.1}
\]
Indeed, \([a]_q\) contains \(\Phi_N(q)\) exactly when \(N\mid a\), and
then with multiplicity one. Consequently,
\[
v_{\Phi_N}\!\left(\frac{[jN]_q!}{[nN]_q!}\right)
=j-n>0\qquad(j>n).                               \tag{4.2}
\]
This is the exact cyclotomic analogue of (1.3). More general cyclotomic
valuations of \(q\)-Pochhammer symbols are studied by
[Adamczewski--Bell--Delaygue--Jouhet](https://arxiv.org/abs/2209.11075).

The mixed-term analogue comes from Gaussian binomial coefficients:
\[
\binom Nk_q\equiv0\pmod{\Phi_N(q)}
\qquad(0<k<N).                                   \tag{4.3}
\]
Thus, for \(q\)-commuting variables \(YX=qXY\), the quantum binomial
theorem gives
\[
(X+Y)^N=X^N+Y^N\pmod{\Phi_N(q)}.                 \tag{4.4}
\]
The same statement holds for the appropriate \(q\)-multinomial or shuffle
product. Congruences at primitive roots of unity are treated in
[Straub, *Congruences for q-binomial coefficients*](https://link.springer.com/article/10.1007/s00026-019-00461-8);
the broader q-Lucas context is surveyed
[here](https://arxiv.org/abs/1409.3820).

Equations (4.2) and (4.4) give a literal blueprint:

1. replace ordinary multiplication of the radial summands by a
   \(q\)-shuffle or \(q\)-commuting product;
2. replace \(j!\) by \([j]_q!\);
3. use the \(N\)-fold quantum power;
4. divide by \([nN]_q!\) before reducing modulo \(\Phi_N(q)\).

The lowest contribution becomes a unit multiple of \(c^N\); every mixed
quantum multinomial is divisible by \(\Phi_N\); every higher pure layer
acquires at least one additional factor of \(\Phi_N\).

### Limitation

This does not reprove the original fixed Gaussian theorem. Vanishing at
\(q=1\) does not imply a polynomial identity in \(q\), and ordinary
commuting powers have ordinary binomial coefficients, not Gaussian ones.
The construction therefore belongs naturally to a **q-Gaussian or
q-shuffle moment problem**. The main research task is to choose a
q-deformation whose moment identities are intrinsic rather than inserted
solely to make the proof work.

## 5. Quantum Frobenius at a root of unity

Quantum Frobenius is the structural version of Section 4. In Lusztig's
divided-power setting, specialization at a root of unity separates a
finite "small quantum" kernel from a classical enveloping-algebra image.
Unlike ordinary Frobenius, it is available at many composite orders.
For arbitrary roots of unity and the necessary qualifications, see
[Lentner](https://arxiv.org/abs/1406.0865). Its use as a characteristic-zero
replacement for geometric Frobenius is already established in a different
setting by
[Kumar--Littelmann](https://annals.math.princeton.edu/2002/155-2/p05).

For the present problem this machinery is probably unnecessary unless the
desired analogue is genuinely noncommutative. The scalar calculation
(4.2)--(4.4) already contains the needed annihilation. Quantum Frobenius
becomes valuable if:

* the angular/radial algebra has a natural braided multiplication;
* divided powers are part of the moment functional; or
* one wants functoriality across representations rather than a single
  polynomial identity.

The danger is categorical overhead: a quantum Frobenius map on a quantum
group does not automatically define the required operation on the
commutative Gaussian coefficient algebra.

## 6. Roots-of-unity filters

For a primitive \(N\)-th root \(\zeta\),
\[
\Pi_{a,N}f(U)
=\frac1N\sum_{s=0}^{N-1}\zeta^{-as}f(\zeta^sU)
=\sum_{j\equiv a\pmod N}f_jU^j.                 \tag{6.1}
\]
If \(f\) has degree \(D\) and \(N>D\), this is an exact projector onto the
single coefficient \(f_aU^a\). Applied to (1.1), it preserves \(cU^n\) and
kills all higher layers.

This looks stronger than prime reduction, but it fails the accessibility
test. From
\[
\mathcal L(F_m)=0\qquad(m\ge1)
\]
one cannot infer
\[
\mathcal L(F_r(\zeta^sU))=0.
\]
The factorial functional is not invariant under \(U\mapsto\zeta U\).
Filtering the moment index \(m\) also does not help: it selects congruence
classes of already-zero moments without changing the mixed terms inside a
fixed moment.

Therefore roots-of-unity filters are:

* useful for an enlarged hypothesis containing rotated or scaled moment
  identities;
* useful as finite exact coefficient extractors in computations; but
* not a standalone replacement in the original theorem.

Combined with a q-deformation, they become the Fourier shadow of the
cyclotomic mechanism in Section 4.

## 7. Function-field valuations and Rees degenerations

Introduce a parameter \(t\) and the radial Rees deformation
\[
F_r(tU)=c\,t^nU^n+\sum_{j>n}c_jt^jU^j.           \tag{7.1}
\]
The \(t\)-adic initial form is \(cU^n\). Hence any identity forcing
\(\mathcal L(F_r(tU))=0\) in \(K(t)\) is impossible:
\[
\mathcal L(F_r(tU))
=c\,n!t^n+O(t^{n+1}).                            \tag{7.2}
\]
This is the clean valuation formulation of lowest-face preservation.
Valuations and their associated graded algebras are a standard source of
toric degenerations; see
[Kaveh--Manon](https://arxiv.org/abs/1610.00298).

Again, the deformation is not free evidence. Fixed-scale vanishing at
\(t=1\) does not imply the identity (7.2) vanishes in \(K(t)\). A useful
function-field attack must construct \(t\) from data already controlled by
the hypothesis, or explicitly pose a family version of the conjecture.

This mechanism does not need a separate mixed-term killer if it is applied
to \(F_r\) after the face power \(r\) has already produced \(c\). If it is
applied before that step, equal-valuation face mixtures must remain: they
are exactly the initial form whose constant term is nonzero.

## 8. Finite-difference and spectral projectors

Let
\[
G(t)=\mathcal L(F_r(tU))
=\sum_{j=n}^D j!c_jt^j.                          \tag{8.1}
\]
Because \(G\) is a polynomial, its \(t^n\)-coefficient can be recovered
from finitely many values by interpolation. More intrinsically, with
\(E=t\,d/dt\),
\[
\Pi_n
=\prod_{\substack{0\le j\le D\\j\ne n}}
\frac{E-j}{n-j}                                  \tag{8.2}
\]
satisfies
\[
\Pi_nG=n!c\,t^n.                                 \tag{8.3}
\]
A polynomial in a multiplicative shift \(G(t)\mapsto G(qt)\) gives a
finite-difference version of the same projector.

Standard forward differences are oriented the other way:
\(\Delta^{d+1}\) annihilates polynomials of degree at most \(d\). To kill
higher bounded layers while retaining the lowest one, one needs the
degree-aware projector (8.2), interpolation, or a limit \(t\to0\).

The obstruction is once more accessibility. Equation (8.2) needs scaled
values of the moment, not merely \(G(1)=0\). Finite differences therefore
provide an exact implementation of a scale-family argument, not a new
source of identities from the fixed-scale hypothesis.

## 9. A rigorous analytic scale theorem

There is an immediate characteristic-zero analogue. For \(R>0\), define
the circular Gaussian functional of precision \(R\) by
\[
\mathcal L_R(U^j)=j!R^{-j}.                      \tag{9.1}
\]
Equivalently, this is the radial moment functional for a circular complex
Gaussian whose covariance is \(R^{-1}\).

> **Scale-family lower-face theorem.** Let
> \(P=\sum_kT^kB_k(U)\). If
> \[
> \mathcal L_R(\operatorname{CT}_T(P^m))=0
> \]
> for every \(m\ge1\) and for an unbounded set of \(R\), then
> \(0\notin\operatorname{conv}(S)\).

**Proof.** If the support straddles zero, the exposed-face argument gives
an \(r\) and (1.1). Then
\[
\mathcal L_R(F_r)
=c\,n!R^{-n}+O(R^{-n-1}).                        \tag{9.2}
\]
Multiplication by \(R^n\) and passage to \(R\to\infty\) gives \(c\,n!=0\),
a contradiction. \(\square\)

This is also a special case of the general principle behind Watson's lemma:
the local lowest-order term controls the leading large-parameter
asymptotic of a Laplace transform; see
[DLMF §2.4(i)](https://dlmf.nist.gov/2.4).

The theorem shows exactly what an analytic analogue should say. It is
stronger in its hypothesis than the original GMC statement: vanishing at
one covariance does not imply vanishing along a covariance family.
However, it removes prime reduction, Frobenius, and factorial divisibility
entirely.

## 10. Comparison

| Mechanism | Preserves face term | Suppresses later mixed terms | Kills higher layers | Available from fixed-scale moments? | Verdict |
|---|---:|---:|---:|---:|---|
| Ordinary roots-of-unity filter | yes | only by residue class | yes if order exceeds degree | no | auxiliary projector |
| Cyclotomic \(q\)-divisibility | yes | yes, with q-binomial/shuffle product | yes via \(v_{\Phi_N}([m]_q!)\) | no, requires q-family | strongest algebraic analogue |
| Quantum Frobenius | yes | yes in the correct divided-power category | yes with cyclotomic normalization | no automatic commutative model | structural q-analogue |
| Function-field valuation | yes | unnecessary after the face step | yes by strict valuation | only with a deformed identity | best abstract language |
| Finite differences | yes | not independently | exactly, with degree bound | no, needs scaled samples | computational/exact projector |
| \(R\to\infty\) asymptotics | yes | unnecessary after the face step | yes by decay | only for a scale-family hypothesis | cleanest analytic analogue |

## 11. Recommended research program

1. **Record the scale-family theorem as the baseline analytic analogue.**
   It is already complete and clarifies the precise extra hypothesis.
2. **Build a q-lower-face toy model.** Define a q-shuffle power and
   \(\mathcal L_q(U^j)=[j]_q!\), then prove the cyclotomic isolation lemma
   using (4.2)--(4.4). Start in the three-level Bessel family before
   treating arbitrary support.
3. **Demand an intrinsic origin for q.** The decisive question is whether
   the q-moment identity comes from a natural q-Gaussian measure,
   braided Wick rule, or representation-theoretic trace.
4. **Use roots-of-unity and finite differences as audit tools.** They can
   verify coefficient isolation in bounded examples, but should not be
   advertised as proofs from the original hypothesis.
5. **Do not expect a function-field valuation alone to strengthen
   GMC(2).** The missing bridge is deformation invariance of moment
   vanishing, not valuation technology.

The central conclusion is therefore qualified but positive: primes can be
replaced. Cyclotomic specialization replaces them algebraically in a
q-deformed problem, and an asymptotic scale replaces them analytically in a
variance-family problem. Neither replacement currently recovers the same
fixed-variance theorem without adding structure.
