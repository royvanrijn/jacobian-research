# Pointwise self-sums for fixed prime cosets

Let `m` be fixed and even, `G=(Z/mZ)^*`, `K<=G`, and `S=aK` a
nonempty prime coset. Let `C_S` contain one and the positive integers
whose prime factors are outside `m` and have residues in `S`. Put

\[
 H=\langle K,a\rangle,\quad \rho=|K|/|G|,\quad
 d=\operatorname{ord}_{G/K}(aK),\qquad
 \Gamma_S(N)=\frac{m}{|H|^2}\#\{(u,v)\in H^2:u+v=N\bmod m\}.
\]

The prime density `rho` is distinct from the density of the generated
group `H`. This distinction is necessary for nonidentity cosets.

The [exact character check](verify_prime_coset_characters.py) and its
[receipt](prime_coset_character_checks.json) verify the quartic Fourier
conjugations, all ten target-residue factors, and two rejection controls
for the `3 mod 10` specialization. This finite algebra check does not
prove the analytic asymptotic or give a numerical onset.

## The pointwise result

Let `E_S` be the good primes outside `S`, and define

\[
 c_S=\frac1{\Gamma(\rho)}\lim_{s\downarrow1}(s-1)^\rho
                         \prod_{p\in\mathcal P_S}(1-p^{-s})^{-1},
 \quad C_{E_S}=\prod_{p\in E_S}(1-(p-1)^{-2}),\quad
 J_S(N)=\prod_{\substack{p\mid N\\p\in E_S}}\frac{p-1}{p-2}.
\]

**Theorem PrimeCosetSelfSums.** Under exactly the ZF and G inputs
stated in [the existing input record](PRIMITIVE_NORM_BINARY_IMPORTED_INPUTS.md),
uniformly on sufficiently large targets with `Gamma_S(N)>0`,

\[
 \#\{x+y=N:x,y\in C_S\}
 =c_S^2\Gamma_S(N)C_{E_S}J_S(N)
       \frac N{(\log N)^{2-2\rho}}
       \bigl(1+O((\log N)^{-\delta_S})\bigr).
\]

For `d>1`, one may take

\[
 g_S=\rho(1-\cos(2\pi/d)),\qquad
 \delta_S=\min\{1/(8192\varphi(m)),3g_S/10\}.
\]

For `d=1` use `delta_S=1/(8192 phi(m))`. The count is identically zero
when `Gamma_S(N)=0`, and positive `Gamma_S` is equivalent to solvability
in every finite quotient. The modulus/coset stay fixed. This is an
unfiltered numerical-value self-sum theorem; no arbitrary mixed-coset
or summand-state-filter assertion is included.

The full-group case `K=G` is elementary: use the exact unit-indicator
model and count the fixed residue classes. The decaying rough-model
bounds in the nontrivial application checks concern `rho<1`.

The standalone proof is [the flagship coset section](../../results/pointwise-binary-splitting-semigroups/cosets.tex),
Theorem `thm:coset`. Its analytic comparison is Theorem
`thm:model-comparison` in [the shared comparison section](../../results/pointwise-binary-splitting-semigroups/comparison.tex).
The proof uses no new energy estimate, graph construction, or
prime-distribution theorem. It is an application of the same exact
imported statements, not an independent validation of their proofs.

## Independently checkable application conditions

The checks below supply the shared comparison theorem's assumptions.
None assumes the desired binary cancellation or its residue-energy
conclusion.

1. **All character branches.** Fourier expansion of `aK` gives
   `F_psi(s)=G_psi(s) product_{chi in Kperp} L(s,psi chi)^(rho barchi(a))`,
   with a uniformly absolutely convergent correction to the right of
   one half. The zeta exponents are exactly
   `alpha_psi=rho psi(a)` for `psi in Kperp`; other fixed twists have
   no branch at one. The full moving-cutoff contour argument is in
   Lemma `lem:coset-centering`, including analytic-log bounds for complex
   powers. A fixed-function Selberg--Delange substitution is not used
   to assert that moving-cutoff estimate.
2. **A nonnegative invariant model.** Write `n=A v` with its small-prime
   record `A`. The rough model uses
   `D(n)=f(A)/V sum_{psi in Kperp} barpsi(v) Mpsi(v)`.
   The identities `M_{psi eta}=Mpsi` for `eta in Hperp` factor its exact
   support by `1_H(v)`. Conjugate branches make it real. On that support
   the positive principal term dominates all secondary terms by the
   relative factor `L^(-3g_S/5)`. Common extensions preserve these
   identities. Multiplication by an allowed small prime changes only
   `A`, so both `f` and `D` are exactly invariant.
3. **One fixed common bank class.** Choose `a0 in aK` and banks from
   primes `a0 mod m`, deleting target divisors. Equivalent fixed-modulus
   characters differ on a product of `s` slots by a unit phase
   `theta(a0)^s`. Both frame products have the same slot set, so the
   phase reproduces exactly and coefficient masses do not increase.
4. **Positive sparse-scale upper envelope.** The independent joint bank
   valuation law and its exact record probabilities are unchanged.
   Forbidden primes have density `1-rho`, so the upper sieve dimension
   is `2-2rho`. Domination by the positive leading model controls all
   four actual/model pair types. CRT gives the interval-length factor
   and power-smaller endpoints. Record tails follow from the same
   tilted Euler moment. These prove conditions (C1)--(C4), rather than
   assuming any additive cancellation.
5. **Main coefficient and local criterion.** Leading residues are uniform
   on `H`. Their fixed root factor is `Gamma_S`, while forbidden good
   primes are determined by `S`, not `H`. The fundamental sieve lemma
   evaluates the model product. A positive word in `S` realizing a
   desired `H` state, with a final CRT/Dirichlet prime, proves the finite
   image statement `lem:coset-finite-images`.

## Explicit checks

For `S={3} mod4`, `H=G`, `rho=1/2`, and the exponents are `1/2,-1/2`.
The general theorem gives exactly

\[
 A_-=\frac12\prod_{p\equiv1(4)}\left(1-\frac2{p(p-1)}\right),
 \qquad R_-(N)=A_-J_-(N)\frac N{\log N}
                       (1+O((\log N)^{-1/16384}))
\]

on all even targets, agreeing with the
[earlier special-case argument](PRIME_THREE_MOD_FOUR_SUMS.md).

For `S={3} mod10`, `H=G`, `rho=1/4`, and the four exponents are
`1/4,i/4,-1/4,-i/4`. The two complex branches have relative size
`L^(-3/20)` and cannot be dropped when arbitrary logarithmic centering
is needed. The forbidden good residues are `1,7,9 mod10`, and
`Gamma_S(N)=5/2` if `5|N`, `15/8` otherwise, on even targets.
Thus all even targets are congruence-admissible and have the predicted
pointwise scale `N J_S(N)/(log N)^(3/2)`, with explicit leading constant
and relative saving `1/32768`.

This generalization has a written application proof and a separate
source/argument audit. It does not claim external human review,
formal verification, an effective numerical onset, or historical priority.
