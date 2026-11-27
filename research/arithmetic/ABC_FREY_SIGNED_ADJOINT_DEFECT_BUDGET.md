# The signed Frey Hecke order and its adjoint Selmer deficit

Use the [fixed-level Frey family](ABC_FREY_FIXED_LEVEL_UNWEIGHTED_MONODROMY.md),
so the conductor is squarefree and equals \(R=\operatorname{rad}(abc)\).
Let \(\delta_1\) be the optimal classical modular degree and let
\(r_{\rm sign}\) be the congruence index of the exact \(R\)-new
space with all signs of the Frey form. Fix a coefficient prime
\(\ell\ge5\), possibly dividing \(R\), and put
\[
 A_\ell=\text{signed Hecke order localized at the residual Frey system},
 \qquad e_\ell=v_\ell(r_{\rm sign}),\qquad D_\ell=v_\ell(\delta_1).
\]
For the Frey augmentation \(\lambda\), put
\(I_\ell=\ker\lambda\) and
\(C_\ell=\operatorname{Hom}(I_\ell/I_\ell^2,\mathbb Q_\ell/\mathbb Z_\ell)\).
Use its [injection into adjoint cohomology](ABC_FREY_STRICT_ADJOINT_CONGRUENCES.md#the-conormal-module-and-its-galois-classes).
Write \(S_\ell\) for the classical torsion adjoint Bloch–Kato
Selmer subgroup of the same global cohomology, and
\(b_\ell=\operatorname{length}_{\mathbb Z_\ell}S_\ell\).

For \(q\mid R\) let \(n_q=v_q(\Delta_E^{\min})\), and let
\(s_{q,\ell}\) at \(q\ne\ell\) be the depth of scalar action
on \(E[\ell^j]\) over the local decomposition group. Put
\(s_\ell=\sum_{q\mid R,\ q\ne\ell}s_{q,\ell}\). In particular
\(s_{2,\ell}=0\) since \(2-1\) is a unit. The signed Wiles
defect has the uniform bound
\[
 0\le \operatorname{length}_{\mathbb Z_\ell}(I_\ell/I_\ell^2)-e_\ell
        \le 2s_\ell .
 \tag{SD1}
\]
Moreover the common subgroup \(K_\ell=C_\ell\cap S_\ell\)
satisfies
\[
 \operatorname{length}_{\mathbb Z_\ell}(S_\ell/K_\ell)
             \le 2s_\ell,\qquad
 \operatorname{length}_{\mathbb Z_\ell}(C_\ell/K_\ell)
             \le 3s_\ell.
 \tag{SD2}
\]
In particular if there is no scalar tame place, the signed
local Hecke ring is a complete intersection and the signed
cotangent dual equals the classical Selmer group as a subgroup
of global cohomology.

## The lower bound after removing the oldforms

Put \(t_{q,\ell}=v_\ell(n_q)\) for every \(q\mid R\).
Start with the full classical Hecke order at level \(R\),
whose congruence exponent at \(\ell\) is \(D_\ell\) by
squarefreeness. Apply the operator of (OA4) at every
\(q\mid R\). The operators remove all the \(q\)-old spaces,
and their Frey scalar has valuation at most \(t_{q,\ell}\)
for \(\ell\ge5\). This includes \(q=\ell\), using the ordinary
projector and the coefficient-prime Tate trace separator.
On the resulting newspace the integral sign projection has
unit Frey scalar at \(\ell\). Lift the signed Frey
idempotent multiple from the quotient and multiply by these
operators as in (OA5). This gives
\[
          D_\ell-\sum_{q\mid R}t_{q,\ell}\le e_\ell .
 \tag{SD3}
\]
[The all-prime oldform comparison, (OA4)–(OA5)](ABC_FREY_COEFFICIENT_PRIME_OLD_ANNIHILATOR.md).
[Agashe–Ribet–Stein, §2 and Theorem 3.6, squarefree congruence exponent](sources/agashe-ribet-stein-modular-degree-congruences.pdf).

Now the exactly calibrated modular degree formula (MD1) gives
\[
       D_\ell=b_\ell+\sum_{q\mid R}t_{q,\ell}+s_\ell .
 \tag{SD4}
\]
Here the \(t_{q,\ell}\) with \(q\ne\ell\) and the
\(s_{q,\ell}\) arise from the tame Euler comparison, including
the \(t_{2,\ell}\) term because this Frey curve is multiplicative
at \(2\). When \(\ell\mid R\), its \(t_{\ell,\ell}\) is the
coefficient-prime term of the ordinary calculation.
Consequently \(e_\ell\ge b_\ell+s_\ell\).
[The good and multiplicative coefficient calibration, (MD1) and (MD8)](ABC_FREY_MULTIPLICATIVE_ADJOINT_DEGREE.md).

## The tame signed local groups

For \(q\mid R,\ q\ne\ell\), let \(\mathcal L_{q,\ell}\)
be the subgroup of local divisible adjoint cohomology with
the fixed-determinant and signed Steinberg matrix conditions
of (SA2). Then
\[
                \mathcal L_{q,\ell}\simeq
                        (\mathbb Z/\ell^{s_{q,\ell}})^3 .
 \tag{SD5}
\]
Here the subgroup is described by the linearized matrix
conditions, irrespective of whether a local Hecke deformation
realizes each of its members. To verify the exact count,
remove the common unramified quadratic twist and write
\[
    F=\begin{pmatrix}q&u\\0&1\end{pmatrix},\qquad
    T=1+n_q e_{12},\qquad D=\mathbb Q_\ell/\mathbb Z_\ell.
\]
Normalize the first-order inertia change to \(x e_{21}\)
using a trace-zero change of basis. The Frobenius change has
trace zero, with diagonal entry \(a\), lower entry \(c\),
and an unrestricted upper entry. Equations (SA2) and the
determinants give
\[
 n_qx=(q-1)x=ux=0,\qquad
 n_q a=n_q c=0,\qquad (q-1)a+uc=0 .
\]
The upper Frobenius entry can be killed since multiplication
by \(q-1\) on \(D\) is surjective. The remaining action on
\((a,c)\) is translation by \((-u,q-1)\) times \(D[n_q]\).
If \(t=v_\ell(n_q)\) and \(s=s_{q,\ell}<t\), divide the row
\((q-1,u)\) by \(\ell^s\) and complete that primitive row to
an integral basis. In the new coordinates the first entry
of \((a,c)\) is in \(D[\ell^s]\). The second entry is in
\(D[\ell^t]\), modulo the subgroup \(D[\ell^{t-s}]\).
This gives two copies of \(\mathbb Z/\ell^s\), and the
\(x\)-coordinate gives the third. For \(s=t\) the row
condition and the translation action vanish on \(D[n_q]^2\),
giving the same formula.

These are sufficient local cocycle conditions. In the
square-zero ring the lifted \(T\) is \(1+\widetilde N\)
with \(\widetilde N^2=0\) by the trace and determinant
identities. Equations (SA2) imply
\(\widetilde F\widetilde N\widetilde F^{-1}=q\widetilde N\)
and hence the tame relation
\(\widetilde F\widetilde T\widetilde F^{-1}=\widetilde T^q\).
At each finite coefficient depth these matrices define an
action of the pro-\(\ell\) tame inertia and the unramified
quotient. The fixed sign and determinant are built in.
Passing to divisible coefficients gives the subgroup used
in (SD5).

The same proof works at the multiplicative prime \(2\)
because the local Tate inertia here is tame and pro-\(\ell\)
in the Steinberg representations. The equations with \(q=2\)
force \(x=0\) and make the two remaining directions disappear
after conjugacy. Thus the classes vanish on the full
decomposition group at \(2\) in this family. This is stronger
than just the unramified condition used for all dyadic
reduction types in the general Frey calculation.
[Local Tate matrices and the signed identities, (SA2)–(SA4)](ABC_FREY_STRICT_ADJOINT_CONGRUENCES.md#the-signed-steinberg-restriction-at-an-odd-prime).

## An exact error budget

Let \(S_\ell^{\rm st}\) enlarge the classical Selmer conditions
only at the primes \(q\mid R,\ q\ne\ell\), where it imposes
\(\mathcal L_{q,\ell}\). The classical condition is zero at
those multiplicative primes. Indeed the rational inertia
invariants are the upper Tate root and \(F-1\) is invertible
on that line over \(\mathbb Q_\ell\). The compact dual Selmer
group orthogonal to the classical conditions is zero by the
adjoint rank-zero result and the vanishing of the divisible
dual invariants. The same Poitou–Tate comparison as in
(BDG4) gives
\[
 S_\ell^{\rm st}/S_\ell \simeq
             \bigoplus_{q\mid R,\ q\ne\ell}\mathcal L_{q,\ell},
 \qquad
 \operatorname{length} S_\ell^{\rm st}=b_\ell+3s_\ell .
 \tag{SD6}
\]
The signed Hecke conormal \(C_\ell\) injects into this enlarged
group. Away from \(R\ell\) its classes are unramified. At
\(q\mid R,\ q\ne\ell\), the signed matrices give the conditions
defining \(\mathcal L_{q,\ell}\). At \(\ell\), (DA2) and
(BK2) place their divisible image in the classical condition,
also when \(E\) is multiplicative there.
[Classical coefficient-prime pullback and tame comparison, (BK2)–(BK3)](ABC_FREY_CLASSICAL_ADJOINT_COMPARISON.md).
[Vanishing of the compact dual and the local comparison, (BDG3)–(BDG4)](ABC_FREY_ADJOINT_DEGREE_CALIBRATION.md).

Write \(w_\ell=\operatorname{length}C_\ell-e_\ell\), the
Wiles defect. It is nonnegative by the augmentation
annihilator inequality. Equations (SD3), (SD4), and (SD6)
give the budget
\[
 w_\ell
 +\operatorname{length}(S_\ell^{\rm st}/C_\ell)
 +(e_\ell-b_\ell-s_\ell)=2s_\ell .
 \tag{SD7}
\]
All three summands are nonnegative, so (SD1) follows.
Equation (SD2) follows by embedding \(S_\ell/K_\ell\)
in \(S_\ell^{\rm st}/C_\ell\) and \(C_\ell/K_\ell\)
in \(S_\ell^{\rm st}/S_\ell\). This controls the actual
subgroups, as well as their orders. When \(s_\ell=0\)
the defect vanishes and \(C_\ell=S_\ell\); the local
order is finite flat, with smooth generic augmentation
because the Frey factor has dimension one over
\(\mathbb Q_\ell\). The Wiles–Lenstra criterion then gives
the complete-intersection assertion.
[Böckle–Khare–Manning, Propositions 3.26 and 3.27](sources/bockle-khare-manning-local-global-wiles-defect-2108.09729v2.pdf).

For all coefficient primes the total possible defect is small.
For each odd \(q\mid R\), take \(u_q^{\rm Tate}\in\mathbb F_q^\times\)
to be the residue of the unit part of a Tate parameter at \(q\),
and set
\[
   B_{\rm scalar}
     =\left(\prod_{\substack{q\mid R\\ q\ {\rm odd}}}
          \gcd\!\left(n_q,
              \frac{q-1}{\operatorname{ord}_{\mathbb F_q^\times}(u_q^{\rm Tate})}
              \right)\right)_{\ell\ge5}.
\]
At each \(\ell\ge5\) its valuation is \(s_\ell\). This follows
by the local Kummer test, and can equivalently be computed
by the square-root residue convention in (MD8). It follows
that
\[
 \prod_{\ell\ge5}\ell^{w_\ell}\ \mid\ B_{\rm scalar}^2,\qquad
     \log B_{\rm scalar}=O\!\left(\frac{\log c}{\log(\log c+2)}\right).
 \tag{SD8}
\]
No count of coefficient primes is used in this estimate.

The inequality \(b_\ell\le e_\ell\) also implies
\(S_\ell[\ell^{e_\ell}]=S_\ell\). The classical Selmer groups
in the degree identity therefore have no additional orders at
primes outside the support of the signed congruence index.

## Exact comparison at good coefficient primes

Suppose \(\ell
mid R\). Choose an even-cardinality set \(Q\)
of conductor primes that contains all primes with positive
\(s_{q,\ell}\). If \(\omega(R)\) is even, take every prime of
\(R\). Otherwise leave out a single prime \(q_0\) at which the
Frey residual representation is ramified. Such an odd prime
exists by the Frey discriminant argument of Pasten's Lemma 6.12.
Write \(N_0=R/\prod_{q\in Q}q\), so \(N_0=1\) or \(q_0\).
[Frey residual image and the ramified prime](ABC_FREY_JOINT_DIVISION_FIELDS.md#the-image-over-the-hecke-order).

Use the anemic image \(T^{\mathrm{st},Q}\) of Böckle–Khare–Manning,
§7, after extension to their coefficient ring. In that section
the Shimura cohomology is localized at the maximal ideal of
the full Hecke algebra determined by the Frey form. At the
primes in \(Q\) the maximal-compact invariants on the division
quaternion side give the Steinberg characters by
Jacquet–Langlands. The local uniformizer Hecke operator acts by
their sign, up to the common normalization, and the two signs
are distinct in the residue field. At the possible prime
\(q_0\) outside \(Q\), residual ramification excludes every
\(q_0\)-old constituent. Its exact-level signs are again
separated by \(U_{q_0}\). Thus after full localization the
good Hecke operators have the same characteristic-zero
eigensystems as on the classical signed space localized at the
Frey residual system. Their integral images agree since both
are torsion-free and embed in the product of these eigensystems.
These operators generate the full signed image by the
prime-to-level and prime-to-\(\ell\) coefficient argument.
Consequently \(T^{\mathrm{st},Q}\) is the base change of
\(A_\ell\). The residual image is full, and its cyclotomic
restriction contains \(\operatorname{SL}_2(\mathbb F_\ell)\),
so the residual hypotheses for the defect formula hold.
[Böckle–Khare–Manning, §7 and Theorem 7.4](sources/bockle-khare-manning-local-global-wiles-defect-2108.09729v2.pdf).
[Integral generation of the signed image](ABC_FREY_JOINT_DIVISION_FIELDS.md).

Theorem 7.4 gives twice the local scalar depth at the primes
in \(Q\), with zero contribution from \(q_0\) when present.
The depths and the \(\mathbb Z_\ell\)-lengths both scale by
the ramification index of the coefficient extension. Hence
\(w_\ell=2s_\ell\). The two other terms of (SD7) vanish:
\[
 e_\ell=b_\ell+s_\ell,\qquad C_\ell=S_\ell^{m st},
 \qquad C_\ell\cap S_\ell=S_\ell
       \qquad(\ell\ge5,\ \ell
mid R).
 	ag{SD9}
\]

## The abc endpoint

The subgroup \(K_\ell\) is exactly the strict signed conormal
group in this family: at \(q\mid R,\ q\ne\ell\) the classical
condition is zero, and the dyadic restriction already vanishes.
Consequently
\[
 (\delta_1)_{\{2,3\}}\prod_{\ell\ge5}|K_\ell|
 \quad\text{and}\quad
 (\delta_1)_{\{2,3\}}\prod_{\ell\ge5}|C_\ell|
\]
each have logarithm \(\log c+\log R+O(\log c/\log(\log c+2))\)
uniformly. This uses (SD2), (SD7), and the previous exact
Selmer length comparison (MD7). The sharp bound with exponent
\(2+\varepsilon\) for either product is equivalent to abc,
even restricted to the fixed \(n=60\) image.
In particular, replacing simultaneous congruence depth by
signed conormal length no longer requires a separate global
defect estimate. The remaining estimate is the size of the
matched global adjoint classes, including the joint budget at
\(2,3\).
