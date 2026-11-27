# Classical dual adjoint classes in ordinary Cartan Kummer groups

Let \(E\) belong to the [oriented squarefree Frey family](ABC_FREY_FIXED_LEVEL_UNWEIGHTED_MONODROMY.md), with conductor \(R=\operatorname{rad}(abc)\). Fix a coefficient prime \(\ell\ge5\). Use
\[
 T=\operatorname{ad}^0(T_\ell E),\qquad T'=T^*(1),\qquad
 W=T\otimes\mathbb Q_\ell/\mathbb Z_\ell,\qquad
 W'=T'\otimes\mathbb Q_\ell/\mathbb Z_\ell,
\]
where the trace pairing identifies \(T^*\) with \(T\). Put \(M_m=W[\ell^m]\), \(M'_m=W'[\ell^m]\), \(S_\ell=H^1_f(\mathbb Q,W)\), and \(S'_\ell=H^1_f(\mathbb Q,W')\). The conditions on the divisible modules are images of the rational Bloch–Kato local conditions. For \(m\ge1\) there is an exact equality
\[
                   |S'_\ell[\ell^m]|=|S_\ell[\ell^m]|.
 \tag{CK1}
\]
Thus these two finite groups have the same elementary divisors, at every depth.

Let \(A_m\) be the nonsplit Cartan subfield of the projective \(\ell^m\)-division field \(P_m\). It has degree \(\ell^{2m-1}(\ell-1)\). The [Cartan Kummer projection](ABC_FREY_CARTAN_ORDINARY_KUMMER.md#projective-kummer-descent-and-the-compatibility-test) injects \(S'_\ell[\ell^m]\) into \((A_m^\times/A_m^{\times\ell^m})^-\). At a tame bad prime \(q\mid R,\ q\ne\ell\), put \(t=v_\ell(v_q(\Delta_E^{\min}))\), let \(s\) be the tame scalar depth from (SD1), and put \(\mu=\min(t,m)\). The simultaneous valuations at the places of \(A_m\) over \(q\) have image on the global group of order at most \(\ell^{\kappa_{q,\ell,m}}\), where
\[
        \kappa_{q,\ell,m}
          =(2\mu-m)_+ + (2\mu-m-t+s)_+ .
 \tag{CK2}
\]
Here \(x_+=\max(x,0)\). The exponent in (CK2) is exact for the possible valuations from the full *local* dual Bloch–Kato pullback. In particular every projected class has valuations zero modulo \(\ell^m\) over \(q\) for \(m\ge 2t\). This includes \(q=2\) with \(s=0\).

## Finite classical duality

At a finite place take the inverse images in finite cohomology of the two classical torsion conditions. They are the reductions of the integral pullbacks of the rational conditions. This follows from the multiplication sequence on \(W\) and \(W'\), as in the general argument for (BK2). At every place these finite local conditions are exact annihilators. Indeed rational Bloch–Kato conditions are orthogonal complements. If \(P\subset H^1(K,T)\) denotes the pullback of the rational condition, the annihilator in \(H^1(K,W')\) is the image of the rational dual condition. Indeed \(P\) contains the torsion and is saturated in \(H^1(K,T)\). Its annihilator by compact local Tate duality is divisible with corank equal to the codimension of the rational condition. The image of the rational dual condition annihilates \(P\) by compatibility of the pairings and has this same corank. A divisible group has no nonzero finite quotient, so the two groups agree. Reduction of \(P\) pairs with a finite dual class through its image in \(H^1(K,W')\). This proves the finite annihilator assertion by taking inverse images.
[Bloch–Kato, Proposition 3.8](sources/bloch-kato-tamagawa-1990.pdf).
[Wiles, proof of Proposition 1.4 and Proposition 1.6](sources/wiles-modular-elliptic-curves-1995.pdf).

The global invariants on \(M_m\) and \(M'_m\) vanish by the full Frey image. For \(M_m\) the centralizer of the full matrix group is scalar and trace zero forces a scalar to vanish. For \(M'_m\) the scalar \(2I\) acts as \(4\), whose difference from \(1\) is a unit. The finite global groups with the two preimage conditions therefore identify with \(S_\ell[\ell^m]\) and \(S'_\ell[\ell^m]\). Apply Greenberg–Wiles to these two finite structures. At each finite \(q\ne\ell\), the size of the primal local condition is \(|H^0(\mathbb Q_q,M_m)|\). In fact the rational unramified primal group has dimension \(\dim H^0(\mathbb Q_q,T\otimes\mathbb Q_\ell)\), by taking invariants and coinvariants of Frobenius on rational inertia invariants. If this dimension is \(h\), the divisible image has corank \(h\), while \(H^0(\mathbb Q_q,W)/\ell^m H^0(\mathbb Q_q,W)\) has order \(\ell^{-mh}|H^0(\mathbb Q_q,M_m)|\). The multiplication sequence multiplies these two orders. At \(\ell\) the local factor is \(\ell^m\) by (DA2) and (BK2), for both good and multiplicative reduction. It cancels the real factor \(\ell^{-m}\), since complex conjugation has a one-dimensional fixed space on trace-zero adjoints. This proves (CK1). Finiteness on the primal side follows from the adjoint degree formula, so (CK1) also gives finiteness on the dual side.
[The coefficient condition, (DA2)](ABC_FREY_DUAL_ADJOINT_SELMER.md).
[Its classical pullback, (BK2)](ABC_FREY_CLASSICAL_ADJOINT_COMPARISON.md).
[The all-coefficient classical degree identity, (MD1)](ABC_FREY_MULTIPLICATIVE_ADJOINT_DEGREE.md).

## The tame valuation calculation

Fix \(q\mid R,\ q\ne\ell\). For the primal adjoint write \(B_{q,m}\) for the inverse image of its classical local condition. This is the kernel of the map to \(H^1(\mathbb Q_q,W)\), because the rational unramified group is zero at a Tate prime. Thus the dual Bloch–Kato pullback is \(B_{q,m}^\perp\). The unramified group in \(W\) is the quotient
\[
 C_q=H^1_{\mathrm{unr}}(\mathbb Q_q,W)
       =D[\ell^t]^2\big/\langle\text{a primitive line of depth }t-s\rangle,
 \qquad D=\mathbb Q_\ell/\mathbb Z_\ell .
\]
The two coordinates are diagonal and lower root. Frobenius acts on their quotient by the rank-one map with column proportional to \((u_q,1-q)\), and its image has depth \(t-s\). The upper root is removed by the surjectivity of \(q-1\) on \(D\). Hence one can identify the displayed quotient with
\(D[\ell^t]\oplus(D[\ell^t]/D[\ell^{t-s}])\), in such a way that the image of the finite inertia invariants is the image of \(D[\ell^\mu]^2\). The integral basis change used to make the column primitive preserves \(D[\ell^\mu]^2\).
[The primal unramified quotient and its Tate coordinates, (BK3)](ABC_FREY_CLASSICAL_ADJOINT_COMPARISON.md).

Write \(L\) for the completion of \(P_m\) at a prime above \(q\), and \(L_0\) for its maximal unramified subextension of \(\mathbb Q_q\). The local projective image is triangular. Its intersection with each conjugate nonsplit Cartan is trivial because an element of a nonsplit Cartan fixing a primitive line must be scalar. Consequently each completion of \(A_m\) over \(q\) is a copy of \(L\). On \(M_m\), tame projective inertia \(J=\operatorname{Gal}(L/L_0)\) acts by conjugation with
\(1+x e_{12}\), \(x\in\ell^\mu(\mathbb Z/\ell^m)\). Its order is \(\ell^{m-\mu}\). In the adjoint matrix the terms in \(x\) and \(x^2\) have sum zero, since the sums of the first and second powers of \(0,\ldots,\ell^{m-\mu}-1\) are both divisible by \(\ell^{m-\mu}\) for \(\ell\ge5\). Thus the inertia norm on \(M_m\) is the scalar \(\ell^{m-\mu}\).

Let \(Y_0\) be the conjugate Cartan vector for a chosen completion, and \(u_L\) the unramified class with Frobenius value \(1\). By the local valuation pairing (OK5), for a local dual class \(\xi\) the projected valuation divided by \(\ell^m\) is
\[
 \big\langle\operatorname{cor}_{L/\mathbb Q_q}(Y_0u_L),\xi\big\rangle .
 \tag{CK3}
\]
In the totally ramified step the corestriction is unramified and has Frobenius value the inertia norm of \(Y_0\). Indeed the coset representatives can be taken in inertia, so the terms in the Frobenius evaluation each have degree one over \(L\). For the subsequent unramified step, coset representatives given by Frobenius powers show that corestriction sends a Frobenius value to the same value in the coinvariants for the base Frobenius. These formulas give the class of \(\ell^{m-\mu}Y_0\) in the base unramified cohomology. For \(\xi\in B_{q,m}^\perp\) it suffices to retain its image in \(C_q\).

As \(Y_0\) ranges through the conjugates, their \(\mathbb Z/\ell^m\)-span is \(M_m\) by residual irreducibility of the adjoint and Nakayama's lemma. A conjugate from the same completion gives the same image in \(C_q\). The Cartan acts trivially on its vector, inertia is absorbed by its norm, and Frobenius acts trivially on the coinvariants. It follows that the simultaneous valuations of local classes in \(B_{q,m}^\perp\) detect exactly the dual of the image of \(\ell^{m-\mu}D[\ell^\mu]^2\) in \(C_q\). The two quotient factors have orders \(\ell^{(2\mu-m)_+}\) and \(\ell^{(2\mu-m-t+s)_+}\). This proves (CK2) by perfect local duality. It is a local computation at all depths and makes no claim that every local pattern occurs globally.
[The valuation pairing and the nonsplit Cartan intersection argument, (OK5)–(OK6)](ABC_FREY_CARTAN_ORDINARY_KUMMER.md).

## The ordinary Cartan bound equivalent to abc

At good primes other than \(\ell\) the finite pullbacks are the finite unramified conditions, by trivial inertia and divisibility of the rational unramified image. After restriction to the unramified projective field the projected valuations are zero. Conversely for an ordinary Kummer class satisfying (OK4), all conjugate evaluations are represented by local units away from \(\ell\). An \(\ell^m\)-th root of a unit at residue characteristic different from \(\ell\) generates an unramified extension, by the étale power map on units. At a good rational prime the projective division field is unramified, so the recovered class is unramified already over the base. The local orthogonality test is therefore necessary only at \(q\mid R\) and at \(\ell\). At \(\ell\nmid R\), the valuations vanish by the good coefficient case of (OK5). At \(\ell\mid R\), (OK6) bounds the coefficient valuation order by \(\ell\), with a possible contribution only when the Tate parameter is an \(\ell^m\)-th power. In particular it vanishes if \(m>t_{\ell,\ell}=v_\ell(v_\ell(\Delta_E^{\min}))\).
[The good and multiplicative coefficient calculations, (OK5)–(OK6)](ABC_FREY_CARTAN_ORDINARY_KUMMER.md).

Let \(G^{f,\mathrm{ord}}_{\ell,m}\) be the subgroup of \(A_m(\varnothing,\ell^m)^-\) selected by (OK3)–(OK4) and by the following finite local tests on the recovered class:
- at each \(q\mid R,\ q\ne\ell\), it belongs to \(B_{q,m}^\perp\);
- at \(\ell\), it belongs to \(L_{\ell,m}^\perp\) of (DA2).

Here \(A_m(\varnothing,\ell^m)\) uses all finite valuations modulo \(\ell^m\), with no primes inverted. The preceding argument gives the exact identity between this group and the subgroup of \(S'_\ell[\ell^m]\) with ordinary Cartan valuations. Its index in \(S'_\ell[\ell^m]\) is bounded by
\[
     [S'_\ell[\ell^m]:G^{f,\mathrm{ord}}_{\ell,m}]
          \ \le\ \ell^{\,\epsilon_{\ell,m}+
                      \sum_{q\mid R,\,q\ne\ell}\kappa_{q,\ell,m}},
 \tag{CK4}
\]
where the index uses the injection into the Kummer group and \(\epsilon_{\ell,m}\) is the indicator in (OK1). The logarithm of the product of these bounds over all \(\ell\) is \(O(\log c/\log(\log c+2))\), independently of the depths. Indeed \(\kappa_{q,\ell,m}\le t_{q,\ell}+s_{q,\ell}\le 2t_{q,\ell}\), and a nonzero \(\epsilon_{\ell,m}\) requires \(\ell\mid v_\ell(\Delta_E^{\min})\). The product is bounded by \((\prod_{q\mid R}v_q(\Delta_E^{\min}))^2\), whose odd part in the place index has the divisor-height bound (SA5) and whose dyadic place contributes a valuation \(O(1+\log c)\).
[The divisor-height estimate (SA5)](ABC_FREY_STRICT_ADJOINT_CONGRUENCES.md).

Write \(e_\ell=v_\ell(r_{\mathrm{sign}})\) for the signed congruence exponent. In this family \(S_\ell[\ell^{e_\ell}]=S_\ell\) by (SD7), so (CK1) gives \(S'_\ell[\ell^{e_\ell}]=S'_\ell\) as well. The product
\[
       (\delta_1)_{\{2,3\}}
          \prod_{\substack{\ell\ge5\\ e_\ell>0}}
                          |G^{f,\mathrm{ord}}_{\ell,e_\ell}|
                 \ \ll_\varepsilon R^{2+\varepsilon}
                 \qquad(\varepsilon>0)
 \tag{CK5}
\]
is equivalent to abc, including on the fixed \(n=60\) images alone. This follows from (CK1), (CK4), and the exact height scale (MD7). In the direction from (CK5) to abc, the lower-order logarithmic error can be absorbed in \(\log c\) first to obtain \(\log c=O(\log R)\), and then in \(R^\varepsilon\).

For an exact product at every coefficient prime, one can use the depth
\[
 m_\ell=\max\!\left(e_\ell,\,2\max(\{0\}\cup
             \{t_{q,\ell}:q\mid R,\ q\ne\ell\}),\,b_\ell^\ast\right),
 \quad b_\ell^\ast=
       \begin{cases}1+t_{\ell,\ell}&\ell\mid R,\\1&\ell\nmid R.\end{cases}
\]
Then \(|G^{f,\mathrm{ord}}_{\ell,m_\ell}|=|S_\ell|\) for \(e_\ell>0\). The field degree increases from that at \(e_\ell\) by a factor \(\ell^{2(m_\ell-e_\ell)}\). The logarithm of the product of these factors is still \(O(\log c/\log(\log c+2))\), since \(m_\ell-e_\ell\le 2\sum_{q\mid R}t_{q,\ell}\).
[The exact signed depth comparison, (SD7)](ABC_FREY_SIGNED_ADJOINT_DEFECT_BUDGET.md).
[The classical Selmer product and height scale, (MD7)](ABC_FREY_MULTIPLICATIVE_ADJOINT_DEGREE.md).

By the ordinary ideal Kummer sequence, these groups have a unit intersection and an image in \(\operatorname{Cl}(A_m)[\ell^m]^-\), with exact product of the two orders. Both pieces must satisfy the global adjoint compatibility and the retained local tests. The bound in (CK5) concerns their actual compatible sizes and the \(2,3\)-primary degree. It has no \(q-1\) normalization. A family with unbounded excess in (CK5) for a fixed \(\varepsilon>0\) would disprove abc.
[Cremona–Lingham, §2, the ideal Kummer sequence](sources/cremona-lingham-good-reduction-kummer.pdf).
