# Frey adjoint Selmer lengths and the classical modular degree

Let \(a+b=c\) be a primitive positive triple with at least two odd
radical primes. Write \(E/\mathbb Q\) for its Frey curve,
\(N=N_E\), \(R=\operatorname{rad}(abc)\), and \(\delta_1\) for
the classical optimal modular degree. Use the signed denominator
exponents \(e_\ell\), the classical torsion adjoint Selmer groups
\(S_\ell\), and the tame lengths \(T_\ell\) from the preceding
comparison. Here \(\ell\ge5\) throughout. Put
\[
 b_\ell=\operatorname{length}_{\mathbb Z_\ell}S_\ell,\qquad
 D_\ell=v_\ell(\delta_1).
\]
At every good coefficient prime there is an exact identity
\[
               b_\ell+T_\ell=D_\ell\qquad(\ell\nmid N).
 \tag{BDG1}
\]
This includes nonminimal residual conductor at the other primes
and additive reduction at \(2\). At a multiplicative coefficient
prime the local Greenberg condition enlarges the classical condition
by a cyclic quotient of order
\(\ell^{v_\ell(v_\ell(\Delta_E))}\). These two comparisons give
a mixed Selmer and modular-degree criterion with a converse to
the abc implication.
[Classical adjoint Selmer notation and the exact tame lengths, (BK1)–(BK8)](ABC_FREY_CLASSICAL_ADJOINT_COMPARISON.md).

## The weight-two integral pairing

Suppose \(\ell\nmid N\). Let \(f\) be the normalized newform for
\(E\), with trivial character and weight two. Use \(\eta_f^\emptyset\)
for the naive congruence ideal of Diamond–Flach–Guo, defined by the
image of the exterior square of their Betti lattice under the
pairing of §5.2. Then
\[
                   v_\ell(\eta_f^\emptyset)=D_\ell .
 \tag{BDG2}
\]
The integral normalization of this pairing matters when \(\ell\)
divides a tame diamond-group order.

For weight two and trivial character the level group in their §4.5
is \(U=U_0(N)\), with constant coefficients. To compare its
Betti lattice with the classical one, take \(U'=U\cap U_1(4)\).
It is sufficiently small, has full determinant, and its index
in \(U\) divides \(|\mathrm{GL}_2(\mathbb Z/4)|=96\).
At \(U'\) their Lemma 4.10 identifies parabolic Betti cohomology
with \(H^1(X_{U'}(\mathbb C),\mathbb Z_\ell)\). The inclusion of
the \(U\)-level Betti module into that at \(U'\) has a left inverse
after division by \([U:U']\), by the restriction and trace maps.
It therefore gives a saturated sublattice. Rationally this is
the pullback of \(H^1(X_0(N)(\mathbb C),\mathbb Q_\ell)\).
The geometric pullback of \(H^1(X_0(N)(\mathbb C),\mathbb Z_\ell)\)
is also saturated, since pushforward composed with pullback equals
the degree of \(X_{U'}\to X_0(N)\), a divisor of \([U:U']\).
The two lattices agree as saturated lattices in the same rational
subspace. The normalized pairing of Proposition 4.8 divides the
pairing at a common auxiliary level by the corresponding group
index. Geometric pullback multiplies cup product by the covering
degree. The group index and covering degree differ at most by the
ineffective central sign, so the normalized pairing at \(U\)
agrees with the usual cup product up to an \(\ell\)-adic unit.

Let \(J_0(N)\to E_0\) be the optimal quotient, and let
\(\phi:X_0(N)\to E_0\) be its parametrization. The connected
kernel of the quotient implies that the pullback of
\(H^1(E_0(\mathbb C),\mathbb Z_\ell)\) in
\(H^1(X_0(N)(\mathbb C),\mathbb Z_\ell)\) is saturated, since the
corresponding map on integral first homology is surjective.
This image has the \(f\)-eigenspace as its rational span, by
multiplicity one and the newform construction in §5.4.
It is therefore the localized Betti lattice used in
\(\eta_f^\emptyset\). The Fricke twist in the pairing of §5.2
acts by a sign on this eigenspace. Its cup product is \(\deg\phi\)
times the perfect elliptic cup product. This proves (BDG2).
[Diamond–Flach–Guo, §§4.2, 4.5 and 5.2–5.4, Proposition 4.8 and Lemma 4.10](https://arxiv.org/abs/2512.02348v1).

## The good coefficient identity

Let \(\Sigma\) contain every bad prime. It is disjoint from \(\ell\).
We may also include \(2\) if \(E\) has good reduction there.
The weight bound and good reduction give
\(\ell\nmid 2N\). The full residual image implies residual
absolute irreducibility on the quadratic subfield of
\(\mathbb Q(\mu_\ell)\), since the restricted image contains
\(\mathrm{SL}_2(\mathbb F_\ell)\). Stable lattices in the
two-dimensional rational representation are homothetic by residual
irreducibility, so the associated adjoint lattice in the theorem
identifies with that of \(E\). The Selmer condition at \(\ell\)
in Diamond–Flach–Guo is the classical one by their Corollary 7.2.
Outside \(\Sigma\) and \(\ell\), the unramified condition agrees
with the classical condition. Their Theorem 7.15 and Corollary 6.7
therefore give, for the group relaxed at \(\Sigma\),
\[
 \operatorname{length}S_\ell^\Sigma
   =v_\ell(\eta_f^\Sigma)
   =D_\ell+\sum_{q\in\Sigma} v_\ell
                \big(L_q^{\rm nv}(\operatorname{Ad}f,1)^{-1}\big).
 \tag{BDG3}
\]
This application permits residual level lowering at every prime
in \(\Sigma\).

The quotient of \(S_\ell^\Sigma\) by \(S_\ell\) is the direct sum
of the full local cohomology groups modulo the classical conditions
at \(\Sigma\). To see surjectivity, write \(W=T\otimes
(\mathbb Q_\ell/\mathbb Z_\ell)\) with \(T=\operatorname{ad}^0T_\ell E\),
and \(T'=\operatorname{Hom}(T,\mathbb Z_\ell)(1)\).
The compact Selmer group for \(T'\) orthogonal to the classical
conditions on \(W\) is zero. Local Tate duality for the rational
finite conditions identifies the compact orthogonal with the
preimage of the rational finite condition on \(T'[1/\ell]\).
In fact the pairing of a compact class against the image of a
rational class in \(W\) is their rational pairing modulo
\(\mathbb Z_\ell\). Vanishing against the entire rational
subspace forces that rational pairing to be zero.
Globally \(H^1(\mathbb Q,T')\to H^1(\mathbb Q,T'[1/\ell])\)
is injective because the divisible dual invariants are zero by
the full image and the scalar \(2I\). Both rational finite Selmer
spaces are zero by the adjoint rank-zero result and duality.
Thus the Selmer comparison sequence in Poitou–Tate gives the
asserted surjectivity. The same argument applies to any finite
set of enlargements of the classical local conditions, including
at the coefficient prime.
[Adjoint rational rank zero and its Tate dual](ABC_CUBIC_ADJOINT_TAMAGAWA_FACTORS.md#selmer-rank-zero-and-descent-from-the-cubic-field).
[Diamond–Flach–Guo, Lemma 8.11 and equation (59)](https://arxiv.org/abs/2512.02348v1).

For \(q\ne\ell\), the local length of the full group modulo the
classical condition equals
\[
 v_\ell\!\big(L_q(\operatorname{Ad}f,1)^{-1}\big)
   +\operatorname{length}\big(H^1(I_q,T)_{\rm tors}^{G_q}\big).
 \tag{BDG4}
\]
This is the local factor in the proof of Diamond–Flach–Guo's
Lemma 8.11 applied to their Tate twist of the adjoint, using
Tate duality and (57). The torsion length in (BDG4) is
\(t_q+s_q\) at odd multiplicative \(q\), \(t_2\) at potentially
multiplicative \(2\), and zero at good primes and potentially
good \(2\). This follows from the Tate inertia calculation (CT2)
and from prime-to-\(\ell\) inertia in the potentially good case.
For any \(\ell\ge5\) and \(q\ne\ell\), the local formula itself
applies regardless of the reduction type at \(\ell\).

The naive and primitive Euler factors coincide at every odd bad
prime. At \(2\) the primitive factor at \(1\) is always an
\(\ell\)-adic unit for the Frey curve. For potentially multiplicative
reduction its inverse is \(1-2^{-2}\). For potentially good
reduction with nonscalar inertia, the trace-zero inertia centralizer
has dimension zero or one. In dimension one Frobenius on it is
a sign, since Frobenius preserves or swaps the two inertial
character lines. For scalar inertia a quadratic twist makes the
Frey curve good with full rational \(2\)-torsion. The trace is
\(\pm1\) by the potentially good dyadic calculation in (SA3)'s
dyadic continuation. The inverse adjoint Euler factor is then
\((1-1/2)((1+1/2)^2-a_2^2/4)=1\).
The naive factor either coincides with the primitive one or is
one. Subtracting (BDG4) from (BDG3) yields exactly (BDG1).
In particular \(|V_{\ell,\infty}|=|U_{\ell,\infty}|=\ell^{D_\ell}\)
at every good coefficient prime, by (BK7).
[Diamond–Flach–Guo, §6.2, Corollaries 6.7 and 7.2, Theorem 7.15 and equations (57)–(59)](https://arxiv.org/abs/2512.02348v1).
[The dyadic twist and full \(2\)-torsion calculation](ABC_FREY_STRICT_ADJOINT_CONGRUENCES.md#unramified-adjoint-classes-at-the-prime-two).

## The local and global ordinary comparison at a multiplicative coefficient

Now let \(\ell\mid N\), and put \(n_\ell=v_\ell(\Delta_E)>0\) and
\(t_{\rm coeff}=v_\ell(n_\ell)\). In the adjoint \(W\) at \(\ell\),
let \(W^+\) denote the divisible upper root
\((\mathbb Q_\ell/\mathbb Z_\ell)(1)\) and \(W^-=W/W^+\).
The Greenberg ordinary condition consists of the classes whose
images in \(H^1(I_\ell,W^-)\) vanish. There is a canonical
identification of the excess over the classical local condition
with the following subgroup of the unramified diagonal quotient:
\[
 \frac{H^1_{\rm Gr}(\mathbb Q_\ell,W)}
      {H^1_f(\mathbb Q_\ell,W)}
    \ \simeq\
   (\mathbb Q_\ell/\mathbb Z_\ell)[n_\ell].
 \tag{BDG5}
\]
For finite coefficients \(M_m\) with the same definition on the
ordinary quotient, the excess over \(L_{\ell,m}\) has order
\(\ell^{\min(m,t_{\rm coeff})}\).

Indeed the quotient by the upper root has a trivial diagonal
line and a further quotient with inverse cyclotomic character.
The latter has zero inertia invariants, since some inertial
cyclotomic value differs from \(1\) by a unit. The unramified
cohomology in \(W^-\) is therefore the unramified cohomology of
the diagonal line. The obstruction to lifting such a class to
the full adjoint is the boundary to \(H^2(\mathbb Q_\ell,W^+)\),
given by cup product with twice the Tate Kummer class, up to
sign. This follows by conjugating the trace-zero diagonal in
the Tate extension. Under local reciprocity, the unramified
diagonal character evaluated on the Tate parameter is its value
on Frobenius times \(n_\ell\). Thus the obstruction on the
diagonal group is multiplication by \(2n_\ell\), with kernel
as in (BDG5). The kernel of the projection to the diagonal
cohomology is the image of the upper Kummer group, which is
the classical condition by the coefficient-prime calculation
in (BK2). The same argument with finite coefficients computes the
obstruction by the finite local Tate pairing and gives the
asserted finite-layer quotient.
[Maksoud, §3.2, equations (21), (22) and (25), the ordinary local condition](https://doi.org/10.1007/s40687-026-00598-y).
[The upper-root Bloch–Kato condition, (BK2)](ABC_FREY_CLASSICAL_ADJOINT_COMPARISON.md#the-coefficient-prime-pullback).

Let \(S_\ell^{\rm Gr}\) use the ordinary condition at \(\ell\)
and the unramified condition at every other finite prime.
The Poitou–Tate surjectivity used above and (BK3) give the exact
order
\[
                 |S_\ell^{\rm Gr}|=|S_\ell|\ell^{T_\ell+t_{\rm coeff}}.
 \tag{BDG6}
\]
If the conditions at every bad prime other than \(\ell\) are
fully relaxed, the additional length relative to (BDG6) is
\(\sum_{q\ne\ell,\ E\text{ potentially multiplicative at }q}
v_\ell(q^2-1)\). This follows from (BDG4), the local inertia
length already present in \(T_\ell\), and the dyadic unit
calculation. Each possible coefficient correction \(t_{\rm coeff}\)
contributes to the logarithm of the product by at most
\(\log(v_\ell(\Delta_E))\). Summed over multiplicative coefficient
primes these corrections are \(O(\log c/\log(\log c+2))\)
by the divisor estimate in (SA5).

The fixed-weight comparison computes the congruence ideal on
the depleted \(\Gamma_1\) cohomology. The tame determinant
cotangent contribution cancels the corresponding diamond-cover
degree. The resulting identity is
\(\operatorname{length}S_\ell+T_\ell+t_{\rm coeff}
=v_\ell(\delta_1)\).
[The integral ordinary modular-degree identity](ABC_FREY_MULTIPLICATIVE_ADJOINT_DEGREE.md).

## A mixed criterion with a converse

Define
\[
 Y(E)=(\delta_1)_{\{2,3\}}
        \prod_{\ell\ge5,\ \ell\nmid N}|S_\ell|
        \prod_{\ell\ge5,\ \ell\mid N}
                 \ell^{\min(D_\ell,b_\ell)}.
 \tag{BDG7}
\]
The bad coefficient term is the common part of the degree and
the classical group order in that characteristic. In particular
the definition does not use \(e_\ell\) or any division tower.
Then \(Y(E)\) divides \(\delta_1\) and
\[
             \log Y(E)=\log c+\log R+o(\log c)
 \tag{BDG8}
\]
uniformly on the whole family. For the upper bound use (BDG1)
and the definition of the common part. For the lower bound use
(DF7) and (BK1): \(e_\ell\le D_\ell\) and
\[
 b_\ell\ \ge\ e_\ell-\min(e_\ell,t_2)
                            -3\sum_{q\ne\ell}\min(e_\ell,s_q).
\]
The same lower bound applies to \(\min(D_\ell,b_\ell)\).
The error summed over all coefficient primes is lower order by
(SA5). Finally the signed denominator comparison (DF8) and its
following height argument give
\(\log((\delta_1)_{\{2,3\}}\prod_{\ell\ge5}\ell^{e_\ell})
=\log c+\log R+o(\log c)\), completing the squeeze.

The good coefficient factors in (BDG7) can be replaced by the
truncated orders \(|S_\ell[\ell^{e_\ell}]|\). For a finite
\(\ell\)-group of length \(b\), the length of its \(\ell^e\)-torsion
is between \(\min(e,b)\) and \(b\). In particular (BDG1) and (DF7)
bound the truncated good length below by \(e_\ell-T_\ell\), so
the same argument proves (BDG8).

Consequently a uniform \(Y(E)\ll_\varepsilon R^{2+\varepsilon}\)
for every \(\varepsilon>0\) is equivalent to full abc, including
the separate elementary small-support cases. In the forward
direction (BDG8) first gives \(\log c=O(\log R)\), after which
the lower-order error in (BDG8) is absorbed using the assumed
bound with a smaller exponent increment. In the reverse direction
the classical degree estimate gives the bound for \(Y(E)\) by
divisibility. A sequence on which \(Y(E)/R^{2+\varepsilon}\)
is unbounded for some \(\varepsilon>0\) would therefore disprove abc.
Since the good coefficient factors equal the corresponding
modular-degree factors after correction, a sequence with
unbounded good coefficient product times \((\delta_1)_{\{2,3\}}\)
over \(R^{2+\varepsilon}\) would already suffice for this
direction. A uniform estimate for those good factors still
requires a height bound on the degree, or another invariant that
can impose that estimate. The integral ordinary comparison
shows \(b_\ell\le D_\ell\) at every multiplicative coefficient
prime, so (BDG7) gives the entire classical adjoint order
there as well.
[Removal of the cap and the exact product, (MD1) and (MD7)](ABC_FREY_MULTIPLICATIVE_ADJOINT_DEGREE.md).
[Signed denominator and classical degree comparison, (DF7)–(DF8)](ABC_FREY_JOINT_DIVISION_FIELDS.md).
[Classical truncated length bound and local loss, (BK1) and (BK8)](ABC_FREY_CLASSICAL_ADJOINT_COMPARISON.md).
