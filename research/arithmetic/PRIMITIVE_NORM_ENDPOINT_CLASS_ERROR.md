# Endpoint logarithmic exponents for primitive norm class errors

The retained class-error sieve gives every exponent strictly below
\(\delta_*(2\operatorname{Pic}(\mathcal O))\). Two elementary refinements
reach the endpoint. Actual subset sums of class switches give a sharper
deficiency cap, and retaining the marking probability gives a uniform
endpoint estimate. If \(q\) is the smallest prime divisor of the nontrivial
doubled class group, the relative class error is
\[
 O_{\mathcal O}\left(
 \frac{(1+\log\log(3N))^{q-2}}
      {(\log(3N))^{(q-1)/(2q)}}\right).
\]
Every even-order doubled class group therefore has the pure endpoint
saving \((\log N)^{-1/4}\), with no loglog loss. A trivial doubled class
group has no class error at all.

This is a deduction from the retained ideal-genus interface, prime-class
Mertens estimate, and uniform prime-square sieve. No new zero-free or
weighted-dilation input is used. Orders are fixed; no optimality,
numerical onset, or varying-discriminant uniformity is asserted.

## 1. Statement and finite-group parameters

Fix an imaginary quadratic order \(\mathcal O\), of discriminant
\(\Delta=f^2D_0<0\), and put
\[
 G=\operatorname{Pic}(\mathcal O),\qquad B=2G,\qquad
 S=\{p:p\mid2\Delta\},\qquad \chi=\chi_{D_0}.
\]
Let \(P\) be the positive primitive values of the principal norm form,
and let \(\Lambda\) be the positive integers represented primitively
over every \(\mathbb Z_p\). For positive-integer sets \(A_1,A_2\), write
\[
 R_{A_1,A_2}(N)=\#\{(a,b)\in A_1\times A_2:a+b=N\},\qquad
 R_A=R_{A,A}.
\]
These are ordered numerical-value counts, without coordinate or ideal
multiplicities. Define
\[
 J_S(N)=\prod_{\substack{p\mid N\\p\notin S,\ \chi(p)=-1}}
                      \frac{p-1}{p-2},\qquad
 L_N=\log(3N),\qquad \ell_N=1+\log L_N.
 \tag{E.1}
\]
Both logarithms are positive for \(N\ge3\), and \(J_S(N)\ge1\).

For nontrivial \(B\), let \(\mathcal P_B\) be its finite set of **all
proper subgroups**, including the zero subgroup. These are not restricted
to maximal subgroups. For \(H\in\mathcal P_B\), set
\[
 r_H=[B:H],\qquad d_H=\frac{r_H-1}{2r_H},\qquad K_H=r_H-2.
 \tag{E.2}
\]
Here \(r_H\ge2\) may be composite, and \(K_H\ge0\). Put
\[
 \Phi_B(N)=\sum_{H\in\mathcal P_B}
              \frac{\ell_N^{K_H}}{L_N^{d_H}},\qquad \Phi_0(N)=0.
 \tag{E.3}
\]

**Theorem EndpointStabilizerClassError.** For every fixed imaginary
quadratic order, uniformly for every integer \(N\ge3\),
\[
 0\le R_\Lambda(N)-R_P(N)
       \ll_{\mathcal O}\frac{NJ_S(N)}{L_N}\Phi_B(N).
 \tag{E.4}
\]
If \(B=0\), then \(P=\Lambda\), so the left side is exactly zero.
On sufficiently large locally admissible targets, the retained classical
lower bound gives
\[
 R_P(N)=R_\Lambda(N)\bigl(1+O_{\mathcal O}(\Phi_B(N))\bigr).
 \tag{E.5}
\]
All target valuations are included. The implied constants and onset are
not evaluated numerically.

For nontrivial \(B\), let \(q\) be the smallest prime divisor of \(|B|\).
Then
\[
 \min_{H<B}r_H=q,\qquad
 d_* =\min_{H<B}d_H=\frac{q-1}{2q}=\delta_*(B),\qquad K_*=q-2,
 \tag{E.6}
\]
and
\[
 \Phi_B(N)\asymp_B
       \frac{\ell_N^{q-2}}{L_N^{(q-1)/(2q)}}.
 \tag{E.7}
\]
Every proper subgroup index is a divisor of \(|B|\) greater than one,
hence at least \(q\); a nonzero map \(B\to C_q\) supplies an index-\(q\)
kernel. For each larger index, the strictly larger \(d_H\) absorbs any
fixed extra power of \(\ell_N\). The index-\(q\) terms give both bounds
in (E.7). The remaining loglog power for odd \(q\) is not asserted optimal.

## 2. Exact inputs and the sharper stabilizer cover

The mathematical inputs from retained proofs are:

1. [The ideal-genus interface and FormIdealShadow](PRIMITIVE_NORM_ALL_ORDERS.md),
   (O.6), Section 2, and (F.2): a locally primitive value has a primitive
   integral invertible ideal in the target genus coset. Good split-prime
   factors can be conjugated independently while preserving norm,
   primitivity and conductor-local ideals. Conductor-divisible values
   are included. The earlier cover (O.7) remains a separate fallback.
2. [MixedPrimeSquareSieve](PRIMITIVE_NORM_MAXIMAL_ORDERS.md), (M.13):
   its bound is uniform in finite mark sets, target valuations and
   actual residue roots. Its proof pays for prime-square moduli.
3. Fixed-order prime-class Mertens (O.8), obtained by partial summation
   from Bourgain--Fuchs, Section 2, numbered part (5), Lemma 2,
   equation (5.2), printed/PDF page 15 of the retained
   [primary edition](sources/bourgain-fuchs-1105.3698v2.pdf).
   The nonambiguous inverse-pair coefficient is \(1/|G|\); their fixed-
   discriminant error \(O_{\mathcal O}(x e^{-c\sqrt{\log x}})\) supplies
   a bounded Mertens remainder. Qualitative density alone would not.
4. Relative assertions use the retained classical OrderValueCounts and
   [CommonFieldGenusInvariance](PRIMITIVE_NORM_LEADING_CONSTANT.md),
   (FM.8)--(FM.9), for a positive lower bound on sufficiently large
   admissible targets.

The proper-subgroup mechanism has classical antecedents in Bourgain--Fuchs.
Their squarefree exceptional-value theorem is not substituted for the
numerical-value convolution estimate here. No priority claim is made.

Choose one representative \(g\) of each inverse pair \(\{g,-g\}\) with
\(2g\ne0\), forming \(T_G\). Let \(\Omega_g\) be the good split primes
whose two prime-ideal classes are \(g,-g\), and set
\[
 \Omega_H=\bigcup_{\substack{g\in T_G\\2g\notin H}}\Omega_g,\qquad
 k_H(a)=\#\{p\in\Omega_H:v_p(a)=1\}.
\]

**Lemma StabilizerDeficiency.** If \(B\ne0\), every primitive form \(Q\)
of this discriminant satisfies
\[
 \Lambda_Q\setminus P_Q\subset
 \bigcup_{H\in\mathcal P_B}\{a:k_H(a)\le[B:H]-2\}.
 \tag{E.8}
\]

**Proof.** For any finite sequence \(v_1,\ldots,v_s\in B\), with
repetitions allowed, let \(T=\sum_i\{0,v_i\}\); for an empty sequence
put \(T=\{0\}\). Suppose \(T\ne B\), and put
\(H=\{h\in B:T+h=T\}\). This is a proper subgroup. Since \(T+H=T\),
the image \(\overline T\) in \(B/H\) has trivial stabilizer. Every
prefix subset sum in this quotient also has trivial stabilizer: a prefix
period would remain a period after adding the tail. Adding any nonzero
\(\overline v_i\) therefore strictly increases the prefix cardinality.
Otherwise the finite union of the prefix with its translate would equal
the prefix, making \(\overline v_i\) a nonzero period. Consequently
\[
 1+\#\{i:v_i\notin H\}\le |\overline T|\le[B:H]-1.
 \tag{E.8a}
\]
This proves the required finite cap without a general sumset theorem.

Now take \(a\in\Lambda_Q\setminus P_Q\). By (F.2), choose one
primitive integral invertible ideal \(I\) of norm \(a\), with
\([I]\in c_Q+B\). This is an existential choice: a conductor-divisible
norm need not determine a unique genus. Freeze its conductor-local
factors and all prime-power factors except good split primes occurring
to exponent one. Each of these has its factor in one of two conjugate
invertible prime ideals. Conjugation produces an increment
\(v_i=\pm2g_i\in B\). The finite local changes are independent and
glue to primitive integral invertible ideals of the same norm. A good
split prime still occupies only one local factor; all conductor-local
ideals are unchanged.

Thus attainable classes include \([I]+T\). If \(T=B\), that set
contains \(c_Q\), contradicting (F.2). For the proper sumset's
stabilizer \(H\), its number of increments outside \(H\) is exactly
\(k_H(a)\): the sign is irrelevant to membership, and ambiguous classes
give zero increments. Apply (E.8a). The same proof permits omitting any
fixed finite collection of extra good primes: freeze those factors too
and form \(T\) only from the remaining switches. \(\square\)

All proper subgroups are necessary for this cap. In \(B=C_4\), two
increments equal to \(1\) give \(T=\{0,1,2\}\) with zero stabilizer.
Both increments lie outside the unique index-two subgroup, contradicting
a putative index-minus-two cap restricted to maximal subgroups.

This is the two-element subset-sum case of the classical stabilizer
method. For primary context, Matt DeVos,
[*A short proof of Kneser's addition theorem for abelian groups*](https://arxiv.org/abs/1303.3539v1),
Theorem 1, page 1, states the general stabilizer inequality. That theorem
is context; (E.8a) is proved directly and imports no Kneser theorem.

The equal-size doubling fibres and disjoint inverse-pair count give
\[
 \sum_{\substack{p\le y\\p\in\Omega_H}}\frac1p
       =d_H\log\log y+O_{\mathcal O}(1).
 \tag{E.9}
\]
Indeed there are \(|G|(r_H-1)/(2r_H)\) inverse pairs outside the
doubling preimage of \(H\), each contributing \(1/|G|\). Deleting a
fixed finite prime set changes only the bounded error. This and the
lemma's final sentence permit two conductors in Section 5.

## 3. A product bound uniform through marking probability one

Put \(f_p=(p-1)/p^2\), the exact proportion of residues modulo \(p^2\)
with valuation one. Higher powers remain allowed.

**Lemma UniformMarkProduct.** Uniformly for \(0\le t\le1\) and large
\(y\),
\[
 \prod_{\substack{p\le y\\p\in\Omega_H}}(1-tf_p)
       \asymp_{\mathcal O}(\log y)^{-td_H}.
 \tag{E.10}
\]
The constants do not depend on \(t\).

**Proof.** For \(u=1/p\), define
\(F(t)=\log(1-t(u-u^2))+tu\). Its derivative is
\[
 F'(t)=\frac{u^2(1-t+tu)}{1-t(u-u^2)}\ge0.
 \tag{E.11}
\]
The denominator is positive. Since \(F(0)=0\), and
\(\log(1-x)\le-x\),
\[
 0\le\log(1-tf_p)+t/p\le t/p^2\le1/p^2.
 \tag{E.12}
\]
Sum this and use \(\sum_p p^{-2}<\infty\) and (E.9). The logarithm
of the product is \(-td_H\log\log y+O_{\mathcal O}(1)\), uniformly
in \(t\); the bounded Mertens error is multiplied by \(t\le1\).
Exponentiation proves the claim. \(\square\)

## 4. Retaining the survival loss on both sieve terms

Fix \(H\), and count \(1\le a<N\) for which both \(a\) and \(N-a\)
avoid good inert primes and \(k_H(a)\le K_H\). Dropping bad-prime and
genus conditions enlarges the hidden-failure count. Write this enlarged
count as \(D_H(N)\).

Take \(z=N^{1/10}\), \(y=z^{1/16}\), initially with \(y\ge e^e\).
For a finite mark set \(\mathcal M\subset\Omega_H\cap[1,y]\), let
\(S_{\mathcal M}\) count the inert-sieved integers with the additional
exclusions \(v_p(a)\ne1\) for \(p\in\mathcal M\). Uniform sieve
(M.13), with \(\epsilon=1/2\), gives
\[
 S_{\mathcal M}\le
 C_{\mathcal O}\frac{NJ_S(N)}{\log y}
        \prod_{p\in\mathcal M}(1-f_p)+C_{\mathcal O}N^{9/20}.
 \tag{E.13}
\]
Every constant is independent of the mark set.

Select each eligible prime independently with probability \(0\le t<1\).
An integer with \(r\) active primes up to \(y\) survives with probability
\((1-t)^r\). For a deficient integer \(r\le k_H(a)\le K_H\), so
its survival probability is at least \((1-t)^{K_H}\). Finite expectation
of (E.13) and (E.10) yield
\[
 D_H(N)\le(1-t)^{-K_H}\left(
 C_{\mathcal O}\frac{NJ_S(N)}{(\log y)^{1+td_H}}
       +C_{\mathcal O}N^{9/20}\right).
 \tag{E.14}
\]
The inverse survival factor multiplies both terms.

For \(K_H>0\), put \(v=\log\log y\ge1\), \(t=1-1/v\). This
choice is made after the uniform finite-mark-set bound: prime families,
root congruences and Mertens thresholds stay fixed. Then
\[
 (1-t)^{-K_H}=v^{K_H},\qquad
 (\log y)^{-td_H}=e^{d_H}(\log y)^{-d_H},
 \tag{E.15}
\]
and
\[
 D_H(N)\ll_{\mathcal O}
 \frac{NJ_S(N)(\log\log y)^{K_H}}{(\log y)^{1+d_H}}
       +N^{9/20}(\log\log y)^{K_H}.
 \tag{E.16}
\]
The ratio of the second term to the first is at most a constant times
\(N^{-11/20}(\log y)^{1+d_H}\), since \(J_S\ge1\). It tends to zero,
so the complete multiplied remainder is absorbed.

For \(K_H=0\), use the deterministic full mark set
\(\Omega_H\cap[1,y]\). Every deficient integer survives. Apply (E.13)
and (E.10) at \(t=1\), with no inverse division or loglog loss. The
same choice would be invalid for positive caps.

As \(\log y=(\log N)/160\), replace the logarithms by (E.1) up to
fixed constants. Sum over all proper \(H\) and the two possible failing
summands; (E.8) proves (E.4) for large \(N\). For nontrivial \(B\),
\(\Phi_B(N)>0\), so an enlarged constant covers the finite remaining
\(N\ge3\). For \(B=0\), the ideal-genus statement gives exact equality.
The retained positive lower bound gives (E.5). \(\square\)

## 5. Two forms with different conductors in the same field

Fix primitive positive definite forms \(Q_i\) of discriminants
\(\Delta_i=f_i^2D_0\), \(i=1,2\). The forms, conductors and common
imaginary quadratic field are fixed. Put
\[
 G_i=\operatorname{Pic}(\mathcal O_i),\quad B_i=2G_i,\quad
 S=\{p:p\mid2\Delta_1\Delta_2\},\quad
 J_S(N)=\prod_{\substack{p\mid N\\p\notin S,\ \chi_{D_0}(p)=-1}}
                  \frac{p-1}{p-2}.
\]
Write \(P_{Q_i}\), \(\Lambda_{Q_i}\) for positive primitive numerical
values and primitive local shadows.

**Theorem EndpointCommonFieldForms.** Uniformly for every \(N\ge3\),
\[
 0\le R_{\Lambda_{Q_1},\Lambda_{Q_2}}(N)-R_{P_{Q_1},P_{Q_2}}(N)
 \ll_{Q_1,Q_2}\frac{NJ_S(N)}{L_N}
                 \bigl(\Phi_{B_1}(N)+\Phi_{B_2}(N)\bigr).
 \tag{E.17}
\]
If \(Q_i'\) is in the same genus as \(Q_i\) in its own discriminant,
then on their common sufficiently large admissible targets,
\[
 R_{P_{Q_1},P_{Q_2}}(N)=R_{P_{Q_1'},P_{Q_2'}}(N)
 \left(1+O_{\mathbf Q}(\Phi_{B_1}(N)+\Phi_{B_2}(N))\right).
 \tag{E.18}
\]
If both \(B_i=0\), (E.17)'s left side is zero and the genus-pair counts
agree for every target, without dividing by a possibly zero count.

**Proof.** FormIdealShadow supplies the target \(c_{Q_i}\) and coset
\(c_{Q_i}+2G_i\), with an existential ideal at conductor-divisible norms.
Apply StabilizerDeficiency to each slot. Freeze additional good-prime
factors in the common finite set \(S\) before forming subset sums. This
preserves the cap, and finite deletion preserves (E.9)'s coefficient.
Both shadows forbid the same good inert primes, while marking primes
are split in that field. The one-sided proof (E.10)--(E.16) applies to
each slot; a union bound gives (E.17), without independence assumptions
on class partitions.

Genus equivalence gives identical local shadows in each slot. Apply
(E.17) to both pairs and use the retained classical lower bound
\(R_{P_{Q_1'},P_{Q_2'}}(N)\gg NJ_S(N)/\log N\) on sufficiently large
admissible targets. Their difference is at most the sum of the errors
from that common shadow, proving (E.18). A trivial \(B_i\) gives exact
equality of its individual shadow and value set. \(\square\)

No different-quadratic-field theorem follows from this proof. It also
gives no new error relative to the separate explicit leading singular
series.

## 6. The earlier familywise bound and concrete improvements

The retained maximal-subgroup cover provides a separate fallback. For
each maximal proper subgroup \(M<B\), define
\[
 d_M^{\rm fam}=\frac{[B:M]-1}{2[B:M]},\qquad
 K_M^{\rm fam}=\sum_{\substack{g\in T_G\\2g\notin M}}
                   (\operatorname{ord}(2g)-2),\qquad
 \Psi_G^{\rm fam}(N)=\sum_{M\text{ maximal}}
            \frac{\ell_N^{K_M^{\rm fam}}}{L_N^{d_M^{\rm fam}}}.
 \tag{E.19}
\]
Set this profile to zero if \(B=0\). The cover (O.7) has precisely
these caps: every inverse family outside \(M\) must have fewer than
\(\operatorname{ord}(2g)-1\) exponent-one primes. The same uniform
marking proof gives (E.4), (E.5), (E.17) and (E.18) with these fallback
profiles. Thus the endpoint optimization and the subset-sum cap are
separate improvements. Suppressing fixed subgroup multiplicities:

| Class group \(G\) | Familywise endpoint bound | Stabilizer endpoint bound |
|---|---|---|
| \(C_2\) | \(0\) | \(0\) |
| \(C_3\) | \(\ell_N/L_N^{1/3}\) | \(\ell_N/L_N^{1/3}\) |
| \(C_4\) | \(L_N^{-1/4}\) | \(L_N^{-1/4}\) |
| \(C_5\) | \(\ell_N^6/L_N^{2/5}\) | \(\ell_N^3/L_N^{2/5}\) |
| \(C_6\) | \(\ell_N^2/L_N^{1/3}\) | \(\ell_N/L_N^{1/3}\) |
| \(C_8\) | \(\ell_N^4/L_N^{1/4}\) | \(L_N^{-1/4}\) |
| \(C_3^2\) | \(\ell_N^3/L_N^{1/3}\) | \(\ell_N/L_N^{1/3}\) |
| \(C_4^2\) | \(L_N^{-1/4}\) | \(L_N^{-1/4}\) |

Every nontrivial even-order \(B\) gives the pure saving \(L_N^{-1/4}\):
index-two stabilizers have zero cap; higher indices have a strictly
stronger logarithmic exponent absorbing their loglog factors. Exponent
two is not required. The discriminant \(-56\) case already has this
endpoint in its [separate C4 proof](PRIMITIVE_NORM_CLASS_NUMBER_FOUR.md).
No optimality of the displayed loglog powers is asserted.

## 7. Exact controls and assurance

The [checker](verify_primitive_norm_endpoint_class_error.py) uses
[finite inputs](primitive_norm_endpoint_class_error_inputs.json) and
writes the [receipt](primitive_norm_endpoint_class_error_checks.json).
It enumerates complete valuation-one residues, all mark and active-set
choices in a four-prime toy space, and fixed finite abelian groups.
Independent subgroup enumeration is compared with maximal homomorphism
kernels. Both familywise and all-proper-subgroup parameters are recorded.
Bounded sequence enumeration checks actual subset sums, stabilizers,
quotient prefixes and translated class targets.

Negative controls reject a wrong local density, omitted inverse survival,
probability one with a positive cap, a corrupted group cap, and restricting
the new cover to maximal subgroups. These are finite algebra/probability
checks, not proofs of Mertens, the sieve, asymptotic bounds or all groups.
The unrestricted statements are proved above.

Run from the repository root:

    python3 research/arithmetic/verify_primitive_norm_endpoint_class_error.py --output research/arithmetic/primitive_norm_endpoint_class_error_checks.json

The [source ledger](primitive_norm_endpoint_class_error_source_ledger.json)
pins the reused proof chain, primary edition and finite controls. This
note has internal proof and exact-computation support. External human
review, formal verification and historical priority are not asserted.
