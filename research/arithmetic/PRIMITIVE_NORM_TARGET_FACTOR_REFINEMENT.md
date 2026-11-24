# Class-specific target factors for primitive norm sums

The mixed-field class-error estimate has a target factor because some
primes dividing the target cannot help change the ideal class of a
summand. That factor does not cancel against the local density of the
shadow count. It can, however, be restricted to the relevant obstruction
classes. The resulting density is usually half the obstruction-prime
density, with an explicitly characterized genus-character exception.

All forms, fields and orders below are fixed. The argument uses the
classical sieve and count comparison in the
[different-fields note](PRIMITIVE_NORM_DIFFERENT_FIELDS_CLASS_ERROR.md),
and the fixed-order prime-class theorem already used there. The additional
input is the classical description of genus characters in an enlarged
quadratic order. No joint independence of ring-class fields is assumed.

## 1. Statement

Let \(Q_i\) be primitive positive definite integral binary forms of
discriminants \(\Delta_i=f_i^2D_i\), where \(D_1,D_2<0\) are distinct
fundamental discriminants. Write
\[
 G_i=\operatorname{Pic}(\mathcal O_{f_i}),\quad B_i=2G_i,\quad
 S=\{p:p\mid2\Delta_1\Delta_2\},\quad
 L_N=\log(3N),\quad \ell_N=1+\log L_N.
\]
Let \(P_i\) and \(\Lambda_i\) be the positive primitive numerical
values and their local shadows. As before,
\(R_{A,B}(N)=\#\{a+b=N:a\in A,b\in B\}\) counts ordered pairs of
values, without representation multiplicity. Use
\[
 J_{12}(N)=\prod_{\substack{p\mid N,\ p\notin S\\
                     \chi_{D_1}(p)=\chi_{D_2}(p)=-1}}
                    \frac{p-1}{p-2}.
\]

For every proper subgroup \(H<B_i\), put
\[
 r_H=[B_i:H],\qquad d_H=\frac{r_H-1}{2r_H},\qquad K_H=r_H-2,
 \qquad A_{i,H}=\{g\in G_i:2g\in H\}.
\]
Let \(\Omega_{i,H}\) consist of good primes split in field \(i\)
whose inverse-pair ideal classes \(g,-g\) satisfy \(2g\notin H\).
Define the smaller prime set and target product
\[
 \mathcal E_{i,H}=\{p\in\Omega_{i,H}:\chi_{D_j}(p)=-1\},\qquad
 W_{i,H}(N)=\prod_{\substack{p\mid N\\p\in\mathcal E_{i,H}}}
                         (1-1/p)^{-1},\quad j\ne i.
 \tag{T.1}
\]
These products are one when their index sets are empty.

The density of \(\mathcal E_{i,H}\) has an exact finite-group
description. Set \(F_i=f_i|D_j|\),
\(G'_i=\operatorname{Pic}(\mathcal O_{F_i})\), and let
\(\pi_i:G'_i\twoheadrightarrow G_i\) extend ideals to the smaller
conductor. There is a nontrivial genus character
\(\lambda_i:G'_i\to\{1,-1\}\) with
\(\lambda_i([\mathfrak p])=\chi_{D_j}(p)\) at every good split prime
\(N\mathfrak p=p\). Its existence and the map \(\pi_i\) are justified
in Section 3. Put
\[
 \alpha_{i,H}=
 \begin{cases}
  1/4,&\lambda_i\text{ is trivial on }\pi_i^{-1}A_{i,H},\\
  (r_H-1)/(4r_H)=d_H/2,&\text{otherwise}.
 \end{cases}
 \tag{T.2}
\]
The first case can occur only when \(r_H\) is even.

**Theorem SubgroupTargetClassError.** For all integers \(N\ge3\),
\[
 \begin{split}
 0\le R_{\Lambda_1,\Lambda_2}(N)-R_{P_1,P_2}(N)
 &\ll_{Q_1,Q_2}\frac{NJ_{12}(N)}{L_N}
       \sum_{i=1}^2\sum_{H<B_i}
       W_{i,H}(N)\frac{\ell_N^{K_H}}{L_N^{d_H}},\\
 W_{i,H}(N)&\ll_{Q_1,Q_2}\ell_N^{\alpha_{i,H}}.
 \end{split}
 \tag{T.3}
\]
If \(B_i=0\), its subgroup sum is empty. In particular, if both
doubled groups are zero, the class error is exactly zero.

On sufficiently large locally admissible targets, division by the
classical shadow lower bound gives the same right-hand subgroup sum as
a relative error. The identical bound holds for replacing either fixed
form by another form in its genus, after adjusting the implied constant.
This uses MixedFieldValueCounts and MixedFieldGenusInvariance,
(D.19)--(D.20), with their existing classical premises. No new leading
asymptotic or numerical onset is asserted.

## 2. What the local normalization does and does not remove

At a good prime let \(I_i=\mathbf1_{\chi_{D_i}(p)=-1}\). The shadow
condition on the line \(a+b=N\) excludes \(a=0\pmod p\) if
\(I_1=1\), and \(a=N\pmod p\) if \(I_2=1\). Write
\(\beta_p(N)\) for its surviving density, and normalize by the two
one-variable densities:
\[
 \sigma_p(N)=\frac{\beta_p(N)}{(1-I_1/p)(1-I_2/p)}.
\]
Direct counting gives the complete table.

| Splitting types | \(\beta_p(N),\ p\nmid N\) | \(\beta_p(N),\ p\mid N\) | \(\sigma_p(N),\ p\nmid N\) | \(\sigma_p(N),\ p\mid N\) |
|---|---:|---:|---:|---:|
| split, split | \(1\) | \(1\) | \(1\) | \(1\) |
| split, inert, in either order | \(1-1/p\) | \(1-1/p\) | \(1\) | \(1\) |
| inert, inert | \(1-2/p\) | \(1-1/p\) | \(p(p-2)/(p-1)^2\) | \(p/(p-1)\) |

Thus the only target correction in the good-prime shadow density is
\(J_{12}\). In particular a split/inert target divisor does not increase
the shadow density. This table concerns the good primes; fixed bad-prime
profiles retain the treatment in (D.18)--(D.19).

For a marking prime in \(\Omega_{i,H}\), let \(b_{p,N}\) be the
conditional probability, after the shadow restrictions, of exact
valuation one in the marked summand. The three cases in (D.10) are
\((p-1)/p^2\), \(1/p\), and zero. Define a target-independent reference
\[
 c_p=\begin{cases}
 (p-1)/p^2,&\chi_{D_j}(p)=1,\\
 1/p,&\chi_{D_j}(p)=-1.
 \end{cases}
\]
For every \(0\le t\le1\) the exact finite product identity is
\[
 \prod_{\substack{p\le y\\p\in\Omega_{i,H}}}(1-tb_{p,N})
 =\prod_{\substack{p\le y\\p\in\Omega_{i,H}}}(1-tc_p)
   \prod_{\substack{p\le y,\ p\mid N\\p\in\mathcal E_{i,H}}}
                  (1-t/p)^{-1}.
 \tag{T.4}
\]
At full marking, \(t=1\), the second factor is exactly the truncated
\(W_{i,H}\). At a target-dividing opposite-inert prime, the baseline
already forbids divisibility of the marked summand, so its marking
increment really is zero.

Both possible values of \(c_p\) satisfy, uniformly for \(0\le t\le1\),
\[
 \log(1-tc_p)=-t/p+O(p^{-2}).
 \tag{T.5}
\]
For \(c_p=1/p\) this is the power series with \(t/p\le1/3\); for
\(c_p=(p-1)/p^2\), replacing \(t/p\) by \(t/p-t/p^2\) changes
the logarithm by \(O(p^{-2})\). The constants are independent of
\(t\), the target and its valuations. The retained class-prime Mertens
estimate therefore bounds the first product in (T.4) by
\(C(\log y)^{-td_H}\). The last product is at most \(W_{i,H}(N)\).

Apply the exact same finite-marking argument to MixedFieldMarkedSieve
(D.11). With \(z=N^{1/10}\), \(y=N^{1/160}\), the number of values
deficient for this subgroup is at most
\[
 C(1-t)^{-K_H}
 \left(\frac{NJ_{12}(N)W_{i,H}(N)}{(\log y)^{1+td_H}}
                       +N^{9/20}\right).
\]
For \(K_H>0\), use \(1-t=1/\log\log y\) once \(y\ge e^e\).
For \(K_H=0\), mark every prime deterministically. This gives the
subgroup term in (T.3), including the remainder
\(N^{9/20}(\log\log y)^{K_H}\). Its ratio to the stated main bound
is at most \(C N^{-11/20}(\log y)^{1+d_H}\), because both target
products are at least one. The remainder is absorbed, including after
survival normalization. The finite smaller targets change only the
constant. The all-proper-subgroup cover and the conductor-local
choices are exactly those in (D.7); none has been changed.

## 3. Joint density from one enlarged order

**Lemma QuadraticClassOverlapDensity.** The set in (T.1) satisfies
\[
 \pi_{\mathcal E_{i,H}}(x)
 =\alpha_{i,H}\operatorname{Li}(x)
       +O_{Q_1,Q_2}\bigl(xe^{-c\sqrt{\log x}}\bigr),
 \qquad
 \sum_{\substack{p\le x\\p\in\mathcal E_{i,H}}}\frac1p
 =\alpha_{i,H}\log\log x+C_{i,H}+o(1),
 \tag{T.6}
\]
for fixed data, where \(c>0\) may depend on that data. In particular,
the Mertens remainder is bounded.

**Proof.** Fix a slot and omit its index. Since \(f\mid F=f|D_j|\),
extension of invertible ideals induces
\[
 \pi:\operatorname{Pic}(\mathcal O_F)
              \longrightarrow\operatorname{Pic}(\mathcal O_f).
\]
This map is surjective. Every invertible ideal class of \(\mathcal O_f\)
has a representative prime to \(F\): multiply by a field element
which makes each of its finitely many local generators at \(p\mid F\)
a unit, using weak approximation, and clear remaining denominators by
an integer prime to \(F\). Away from \(F\), both orders are maximal;
extension to the maximal order and contraction to \(\mathcal O_F\)
preserve the ideal and its norm locally. At \(p\mid F\), take the
unit ideal. These local descriptions give an invertible
\(\mathcal O_F\)-ideal extending to the selected representative.
This is the prime-to-conductor ideal correspondence described on
page 13 of the primary source below, with the same approximation
argument at the additional finite primes.

The order of conductor \(F\) has discriminant \(\Delta'=F^2D_i\).
The integer \(D_j\) is a fundamental divisor of this discriminant:
it divides \(\Delta'\), and
\[
 \Delta'/D_j=f_i^2D_iD_j\equiv0\text{ or }1\pmod4.
\]
Proposition 3.4 and Theorem 3.5 of Ibukiyama,
[*Genus character \(L\)-functions of quadratic orders in an adelic way
and maximal orders of matrix algebras*, arXiv:2303.14983v1](https://arxiv.org/abs/2303.14983v1),
pages 12--14, give a genus character \(\lambda\) of \(G'\) such that
\[
 \lambda([\mathfrak p])=\chi_{D_j}(N\mathfrak p)=\chi_{D_j}(p)
 \quad(p\notin S,\ p\text{ split in field }i).
 \tag{T.7}
\]
Its reciprocal fundamental-discriminant pair is not the trivial pair
\(\{1,D_i\}\), since \(D_j\ne D_i\) and \(D_j<0\). Thus
\(\lambda\) is nontrivial. Narrow and ordinary ideal classes agree
in the imaginary quadratic field. The
[retained primary edition](sources/ibukiyama-genus-characters-2303.14983v1.pdf)
was inspected at the actual statements and their proof. No assertion
that a conductor-divisible ideal character is determined by its norm
is used; (T.7) is applied only at the displayed good primes.

Put \(A'=\pi^{-1}A_{i,H}\). Doubling onto \(B_i\) and surjectivity
of \(\pi\) give \([G':A']=r_H\). If \(\lambda|_{A'}\ne1\), exactly
half of \(A'\), and exactly half of \(G'\), have character value
minus one. Consequently
\[
 \#\{u\notin A':\lambda(u)=-1\}
     =\frac{|G'|}{2}-\frac{|A'|}{2}
     =\frac{|G'|(r_H-1)}{2r_H}.
\]
If \(\lambda|_{A'}=1\), the corresponding number is \(|G'|/2\).
In this second case \(\lambda\) descends to a nontrivial character
of \(G'/A'\), so its order \(r_H\) must be even.

The selected class set is closed under inverses. None of its classes
is ambiguous: \(2u=0\) would imply \(2\pi(u)=0\in H\). Its
inverse pairs therefore each have rational-prime density \(1/|G'|\).
Apply Bourgain--Fuchs, arXiv:1105.3698v2, Section 2, numbered part (5),
Lemma 2, equation (5.2), page 15, to these finitely many classes at the
single fixed discriminant \(\Delta'\). That exact statement, including
the ambiguous-class normalization, was reread in the
[retained edition](sources/bourgain-fuchs-1105.3698v2.pdf).
The selected number of classes divided by \(2|G'|\) is precisely
(T.2). This proves the prime-counting formula. Partial summation gives
the convergent Mertens remainder in (T.6), since
\(\int_2^\infty e^{-c\sqrt{\log t}}dt/t<\infty\).
Deleting the common fixed set \(S\) changes only constants. \(\square\)

This proof uses a single enlarged order, not a product-distribution
assumption for two ring-class fields. All enlargement parameters are
fixed. The higher conductor affects constants, not the target-uniform
exponents.

## 4. The compressed endpoint profile

For a nonzero \(B_i\), let \(q_i\) be the smallest prime divisor
of \(|B_i|\), and set \(\delta_i=(q_i-1)/(2q_i)\). Define
\[
 \beta_i=
 \begin{cases}
  \delta_i/2,&q_i\text{ odd},\\
  1/4,&q_i=2\text{ and }\lambda_i\text{ is trivial on }
                       \ker(2\pi_i:G'_i\to B_i),\\
  1/8,&q_i=2\text{ and that restriction is nontrivial}.
 \end{cases}
 \tag{T.8}
\]
The second case means exactly that the opposite-field quadratic
character descends to a nontrivial character of the doubled class
group. It is an alignment condition, not an assumption of generic
independence.

**Corollary CompressedMixedFieldEndpoint.** For all fixed data above,
\[
 R_{\Lambda_1,\Lambda_2}(N)-R_{P_1,P_2}(N)
 \ll_{Q_1,Q_2}\frac{NJ_{12}(N)}{L_N}
       \sum_{i:B_i\ne0}
          \frac{\ell_N^{q_i-2+\beta_i}}{L_N^{\delta_i}}.
 \tag{T.9}
\]
The same sum bounds the relative error on sufficiently large
admissible targets and the relative change under the fixed genus
replacements described after (T.3).

To prove the bound for \(W_{i,H}\) in (T.3), split its prime divisors
at \(x=\log N\), for large \(N\). Formula (T.6) bounds the small-prime
product by \(C(\log x)^{\alpha_{i,H}}\). At most
\(\log N/\log x\) distinct prime divisors exceed \(x\), so
\[
 \sum_{\substack{p\mid N\\p>x}}-\log(1-1/p)
 \le\frac{2\log N}{x\log x}=o(1).
 \tag{T.10}
\]
Here \(-\log(1-1/p)\le2/p\) suffices. This proves the asserted
uniform bound, with small \(N\) absorbed.

For the compression, \(d_H\) strictly increases with \(r_H\), and
the smallest possible proper-subgroup index is \(q_i\). Contributions
with larger index are absorbed by any fixed power saving in \(L_N\),
despite their fixed powers of \(\ell_N\). At index \(q_i\) the cap
is \(q_i-2\). If \(q_i\) is odd, (T.2) always gives
\(\alpha_{i,H}=\delta_i/2\).

If \(q_i=2\), put \(T_i=\ker(2\pi_i)\). A character trivial on
some \(\pi_i^{-1}A_{i,H}\) must be trivial on \(T_i\), since
\(T_i\subset\pi_i^{-1}A_{i,H}\). Conversely, when \(\lambda_i\)
is trivial on \(T_i\), it induces a nontrivial sign character of
\(B_i\); its kernel is an index-two subgroup attaining
\(\alpha_{i,H}=1/4\). If it is nontrivial on \(T_i\), every
index-two subgroup has \(\alpha_{i,H}=1/8\). These are exactly the
three cases in (T.8), proving (T.9).

## 5. Exact maximal order of the target product

The exponent of \(W_{i,H}\) itself cannot be decreased in a uniform
all-target bound. This statement concerns the product, not an optimal
lower bound for the actual class error.

**Theorem TargetFactorMaximalOrder.** For each fixed pair \((i,H)\),
there is a positive constant
\[
 C^*_{i,H}=\lim_{x\to\infty}
   \frac{\displaystyle\prod_{p\le x,p\in\mathcal E_{i,H}}(1-1/p)^{-1}}
        {(\log x)^{\alpha_{i,H}}}
\]
and
\[
 \limsup_{N\to\infty}
        \frac{W_{i,H}(N)}{\ell_N^{\alpha_{i,H}}}=C^*_{i,H}.
 \tag{T.11}
\]

Indeed, the convergent constant in (T.6), together with the convergent
sum of the higher terms in \(-\log(1-1/p)\), gives the positive
Euler-product limit. Splitting at \(x=\log N\) and using the
\(o(1)\) tail in (T.10) proves the upper limsup. For the reverse
inequality take
\[
 N_x=\prod_{\substack{p\le x\\p\in\mathcal E_{i,H}}}p.
 \tag{T.12}
\]
The prime-counting formula in (T.6) gives
\(\log N_x=\alpha_{i,H}x+o(x)\). Hence
\(\ell_{N_x}/\log x\to1\), and the Euler-product limit gives
the opposite inequality in (T.11).

## 6. Three concrete comparisons

The first slot in this table is the principal form of the stated
discriminant. The second slot has class number one, so its class-error
term vanishes.

| First discriminant | Second form | \(\delta_1\) | \(\beta_1\) | Relative class-error bound |
|---|---|---:|---:|---:|
| \(-23\) | \(x^2+y^2\) | \(1/3\) | \(1/6\) | \(O(\ell_N^{7/6}/L_N^{1/3})\) |
| \(-56\) | \(x^2+xy+2y^2\) | \(1/4\) | \(1/4\) | \(O(\ell_N^{1/4}/L_N^{1/4})\) |
| \(-56\) | \(x^2+y^2\) | \(1/4\) | \(1/8\) | \(O(\ell_N^{1/8}/L_N^{1/4})\) |

These statements concern sufficiently large admissible targets. For
the first row every positive target is locally admissible, as proved
in (D.21).

For discriminant \(-23\), the retained class-group calculation gives
\(G=C_3\). Thus \(q=3\), and odd index forces the half-density
case. The enlarged order can be taken to have discriminant \(-368\),
whose relevant reciprocal pair is \((-4,92)\). The previous coarse
relative bound \(O(\ell_N^{5/4}/L_N^{1/3})\) is therefore improved to
the first row.

For discriminant \(-56\), the
[retained class calculation](PRIMITIVE_NORM_CLASS_NUMBER_FOUR.md)
gives \(G=C_4\), with class parity detected by \(\chi_8\). Its unique
proper subgroup of \(B=C_2\) is zero; the marker primes are precisely
the split primes with \(\chi_8(p)=-1\). At good split primes,
\(\chi_{-56}(p)=\chi_8(p)\chi_{-7}(p)=1\), so
\(\chi_{-7}(p)=\chi_8(p)\). Every marker is therefore inert in
\(\mathbb Q(\sqrt{-7})\), giving the full quarter density in the
second row. This is a concrete counterexample to unconditional halving.

For the Gaussian second slot, the good prime
\(23=3^2+14\cdot1^2\) is principal in the first field and inert in
the second. Its class in the enlarged order has \(\lambda=-1\)
but projects to zero in \(G\). Thus \(\lambda\) is nontrivial on
\(\pi^{-1}A\), giving \(\alpha=\beta=1/8\). Equivalently, the
conditions \(\chi_8=-1\), \(\chi_{-7}=-1\), \(\chi_{-4}=-1\)
select exactly three of the 24 units modulo 56.

Finally, for the first row use the primorial targets (T.12). Every
prime divisor is split in the \(-23\) field and inert in the Gaussian
field, so \(J_{12}(N_x)=1\). All these targets are admissible, and
(D.19) gives
\(R_{\Lambda_1,\Lambda_2}(N_x)\asymp N_x/L_{N_x}\), whereas
\(W_{1,0}(N_x)\sim C^*\ell_{N_x}^{1/6}\to\infty\).
Consequently the target factor in the existing sieve bound cannot be
absorbed into a larger shadow-count normalization. This does not prove
that the actual class error attains that bound; a different class-error
argument could still improve it.

## 7. Exact controls and scope

The [checker](verify_primitive_norm_target_factor_refinement.py) directly
enumerates the stated prime-square residue spaces, verifies the
normalized local table and (T.4) for rational marking probabilities,
and counts character fibres in the listed finite abelian groups and
all their relevant subgroups. It checks the index-two alignment
criterion, the modulo-56 examples and the exact exponents in the table.
Negative controls reject universal density halving, an incorrect exact
quarter density in the odd-index calculation, a spurious shadow
factor at an exclusive target divisor, and a spurious marker gain
there. These are finite exact controls; the distribution, maximal-order
and all-target estimates are the written proofs with their explicit
primary premises.

Run from the repository root:

    python3 research/arithmetic/verify_primitive_norm_target_factor_refinement.py

The [receipt](primitive_norm_target_factor_refinement_checks.json) records
the actual replay. The [source ledger](primitive_norm_target_factor_refinement_source_ledger.json)
pins the proof, checker, primary editions and inherited results. The
classical lower comparison remains separate from the stronger analytic
imports used elsewhere for leading asymptotics.

No optimality of the actual class error, uniformity in varying forms or
conductors, numerical cutoff, external review, formal verification or
historical priority is asserted.
