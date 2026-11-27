# The adjoint modular-degree identity at multiplicative coefficient primes

Let \(E\) be the Frey curve of a primitive positive triple \(a+b=c\)
with at least two odd radical primes. Use \(N=N_E\), \(R=\operatorname{rad}(abc)\)
and the classical optimal modular degree \(\delta_1\). For
\(\ell\ge5\), write \(S_\ell=H^1_f(\mathbb Q,\operatorname{ad}^0
T_\ell E\otimes\mathbb Q_\ell/\mathbb Z_\ell)\) and retain the
prime-to-\(\ell\) tame length \(T_\ell\) in (BK3)–(BK7).
Set \(s_{\ell}^{\rm coeff}=v_\ell(v_\ell(\Delta_E))\) if
\(\ell\mid N\), and \(s_{\ell}^{\rm coeff}=0\) otherwise. Then
\[
 v_\ell(\delta_1)=
   \operatorname{length}_{\mathbb Z_\ell}S_\ell+
                         T_\ell+s_{\ell}^{\rm coeff}
                     \qquad(\ell\ge5).
 \tag{MD1}
\]
The good coefficient case is (BDG1). The multiplicative case
uses the weight-two ordinary Hecke ring. Its tame determinant
characters contribute the same length as the tame diamond-cover
degree, so the fixed-determinant comparison retains only the
classical modular degree and the explicit local factors.
[Good coefficient formula and the ordinary localization indices](ABC_FREY_ADJOINT_DEGREE_CALIBRATION.md).

## The fixed-weight ordinary rings

Fix \(\ell\mid N\) and work over \(O=\mathbb Z_\ell\).
Put \(\Sigma=\{q:q\mid N,\ q\ne\ell\}\). Let \(B\) be the
weight-two fixed-cyclotomic-determinant deformation ring of
\(\bar\rho_{E,\ell}\), ordinary at \(\ell\), unrestricted at
\(\Sigma\), and unramified outside \(\Sigma\cup\{\ell\}\).
The ordinary quotient is the one on the Tate module. Residual
inertia characters \(1,\bar\chi_\ell\) are distinct, so a
semistable lift of this residual representation is ordinary.
For a finite square-zero lift, its adjoint condition is the
kernel of restriction to inertia on the quotient by the upper
root. These are the definitions and local calculation of
Darmon–Diamond–Taylor §2.4, including Lemma 2.25(b).
The global residual image is full, so it is absolutely
irreducible on the quadratic subfield of \(\mathbb Q(\mu_\ell)\).
Diamond's Theorem 6.1 applied with \(\Sigma\cup\{\ell\}\)
therefore makes \(B\) finite flat and a complete intersection.
The prime-to-\(\ell\) finite-order determinant twist in his
deformation problem is trivial here because its reduction is
trivial. His theorem allows the dyadic residual types in the
Frey family.
[Diamond, §§3–6 and Theorem 6.1](https://www.wstein.org/people/kiming/boston/diamond.pdf).
[Darmon–Diamond–Taylor, §2.4 and Lemma 2.25](https://wstein.org/edu/Fall2003/252/references/darmon-diamond-taylor/Darmon-Diamond-Taylor-Fermats_Last_Theorem.pdf).

Let \(A\) be the specialization at weight two and trivial
wild coefficient character of the full ordinary algebra
of Emerton–Pollack–Weston, with the same residual ordinary line
and tame set \(\Sigma\). Their §2.4 and control theorem
realize \(A\) on the ordinary cohomology
of \(X_1(M)\), where
\[
 M=\ell\,N(\bar\rho_{E,\ell})
       \prod_{q\in\Sigma}q^{\dim_{\mathbb F_\ell}(\bar\rho_{E,\ell})_{I_q}} .
\]
Here \(N(\bar\rho_{E,\ell})\) is the prime-to-\(\ell\) Artin
conductor. Comparing the inertia invariant dimensions and the
equal Swan conductors of the lift and its reduction gives
\[
                   M=N\prod_{\substack{q\in\Sigma\\
                                        E\ {\rm multiplicative\ at}\ q}}q .
 \tag{MD2}
\]
The dimension of the invariant space of the lift is one at
multiplicative reduction and zero at additive reduction. In the
conductor comparison the coinvariant dimension equals the
invariant dimension for the residual two-dimensional representation:
the wild and prime-to-\(\ell\) tame inertia invariants are exact
direct summands, and the \(\ell\)-part of tame inertia is cyclic.
The augmentation of \(A\) is the eigensystem of the
\(\Sigma\)-depletion \(f_\Sigma\) of the Frey newform \(f\).
The algebra has \(U_q=0\) at \(q\in\Sigma\).
At the coefficient prime \(f\) is already new and ordinary.

Write \(\Delta\) for the \(\ell\)-primary quotient of
\((\mathbb Z/(M/\ell))^\times\), and put \(C=O[\Delta]\).
After the weight specialization the determinant of a lift to
the ordinary ring is \(\chi_\ell\) times a character through
\(\Delta\). This follows either from the diamond-operator
determinant in the ordinary family or from class field theory:
\(\chi_\ell^{-1}\det\rho\) is an \(\ell\)-primary character
unramified at \(\ell\) and outside \(\Sigma\), and at \(q\in\Sigma\)
such characters factor through the \(\ell\)-primary part of
\(\mathbb F_q^\times\). Conversely each such character has a
unique square root, since \(\ell\) is odd. Twisting an ordinary
lift by the inverse square root of its determinant character
gives a fixed-determinant lift. The square root is unramified
at \(\ell\), so the ordinary quotient and its unramifiedness
are preserved.

The resulting product decomposition can also be read directly
on the integral Hecke rings. Over \(\overline{\mathbb Q}_\ell\)
the ordinary depleted eigenpackets of \(A\) are pairs
\((g,\vartheta)\), where \(g\) belongs to Diamond's
fixed-determinant packet and \(\vartheta\) is a character of
\(\Delta\) giving the twist. At \(q\in\Sigma\) the twist
leaves the Swan conductor unchanged. The exponent of \(q\) in
the available level is that Swan conductor plus two, so the
twisted newform and its \(q\)-depletion both fit in the level.
Indeed the degree of the local depletion polynomial is the
inertia invariant dimension, which is precisely the gap between
the newform conductor exponent and this available exponent.
The ordinary line at \(\ell\) is unchanged by the unramified
twist, and \(\vartheta\) has trivial residual reduction.
These observations apply also to untwisting, giving the claimed
packet bijection.

In the product of these eigenvalue fields, the trace at a
prime \(r\nmid M\) is \(a_r(g)\vartheta(r)\). It lies in the
image of \(B\otimes_O O[\Delta]\), and these traces generate
\(A\) by Proposition 2.4.2. In the other direction the
diamond operators in \(A\) have eigenvalues
\(\vartheta(d)^2\). Squaring permutes \(\Delta\), so their
images give all of \(1\otimes O[\Delta]\). Dividing each
trace by its diamond-character unit recovers the corresponding
trace generator of \(B\otimes 1\); traces away from \(M\)
generate \(B\) by Diamond's definition. The two orders are
therefore equal inside the same product of fields, giving
\[
                      A\simeq B\otimes_O C .
 \tag{MD3}
\]
In particular it keeps the nonreduced tame-character directions
in the special fiber.
[Emerton–Pollack–Weston, Theorem 2.1.1, determinant in §2.2, Theorem 2.4.1, Proposition 2.4.2 and the ordinary universality at the end of §2.4](https://www.math.uchicago.edu/~emerton/pdffiles/iwasawa.pdf).

For an augmented ring \(A\) or \(B\) write
\(\Phi=\ker(\lambda)/\ker(\lambda)^2\) for its cotangent module
and \(\eta=\lambda(\operatorname{Ann}\ker(\lambda))\) for its
congruence ideal. From (MD3), \(\Phi_A\simeq\Phi_B\oplus\Phi_C\),
and \(\operatorname{length}\Phi_C=v_\ell(|\Delta|)\).
The group ring \(C\) is a complete intersection, using the
relations \(x_i^{|\Delta_i|}=1\) for cyclic generators.
The tensor in (MD3) is a complete intersection by combining
these relations with a presentation of \(B\). The equality
between conormal length and congruence valuation for a finite
flat complete intersection then gives
\[
                      v_\ell(\eta_A)
                 = v_\ell(\eta_B)+ v_\ell(|\Delta|).
 \tag{MD4}
\]
[Wiles, Appendix, Proposition 2](https://jontalle.web.engr.illinois.edu/Public/Wiles-Fermat.95.pdf).

## The integral pairing at the depleted level

The ordinary, localized \(H^1(X_1(M)(\mathbb C),O)\)
is free of rank two over \(A\) by ordinary cohomological
multiplicity one. Its Fricke-twisted cup product is perfect and
Hecke balanced. Indeed untwisted cup product pairs the ordinary
and anti-ordinary summands perfectly, and Fricke transports one
to the other. This remains true on the indicated localized
summand. Thus \(\eta_A\) is the ideal of this pairing on its
saturated \(f_\Sigma\)-eigenlattice. For the algebra check,
in a free \(A^2\)-basis the pairing is given by a linear
functional \(l:A\to O\) applied to the determinant. Perfection
identifies \(A\) with \(\operatorname{Hom}_O(A,O)\) via \(l\).
The saturated eigenlattice is \(uO^2\), where \(u\) generates
\(A\cap\mathbb Q_\ell e_{f_\Sigma}\). Under the identification
using \(l\), \(u\) maps to a unit times the augmentation, so
\(l(u)\) is a unit and \(l(u^2)=\lambda(u)l(u)\) generates
\(\eta_A\).
[Emerton–Pollack–Weston, Propositions 3.1.1 and 3.3.1](https://www.math.uchicago.edu/~emerton/pdffiles/iwasawa.pdf).
[Maksoud, Proposition 6, Proposition 7 and Theorem 8](https://doi.org/10.1007/s40687-026-00598-y).

First compare on \(X_0(N)\) and \(X_0(M)\). At each prime
\(q\) appearing in \(M/N\), let \(\alpha,\beta\) be the two
degeneracy maps at that step of raising the level from \(L\)
to \(Lq\). Here \(v_q(L)=1\). On differentials the pullbacks
of a cusp form \(g\) are \(g(z)\) and \(qg(qz)\).
The depletion map on the \(f\)-space is
\(\gamma_q=\alpha^*-(a_q(f)/q)\beta^*\).
The same construction after the other depletions uses the
unchanged \(U_q\)-eigenvalue \(a_q(f)=\pm1\).
Writing \(W_L\) for Fricke at level \(L\),
\[
                 \gamma_q^t W_{Lq}\gamma_q
                  =-a_q(f)(1-q^{-2})W_L
 \tag{MD5}
\]
on the \(U_q\)-eigenspace of the original form at that step.
Here \(t\) is the adjoint for untwisted cup product. Indeed
Fricke interchanges the two pullbacks composed with \(W_L\).
Each pullback composed with its own trace is multiplication
by \(q\). The cross compositions are \(U_q\) and its adjoint.
On the \(q\)-new subspace they agree, since \(U_q\) there
acts by minus the involution at \(q\). Thus the scalar on
the right is \(a_q(f)(1+q^{-2})-2a_q(f)\).
The product of (MD5) over all raised primes computes the
change in the twisted pairing.
[Degeneracy and Fricke computation, Diamond–Flach–Guo, §5.2 and Proposition 6.5 in the special local case](https://arxiv.org/abs/2512.02348v1).

The composite depletion map from the integral \(f\)-lattice
at level \(N\) has saturated image in the rational
\(f_\Sigma\)-space at level \(M\). For this assertion use the
pullback from the optimal elliptic quotient \(E_0\), which is
the saturated lattice at level \(N\). After clearing denominators
prime to \(\ell\), the composite is the cohomological pullback
of a homomorphism \(J_0(M)\to E_0\). The differential at
infinity of its induced modular parametrization has first
coefficient an \(\ell\)-adic unit: it is the original Manin
constant times a unit from these denominators, because normalized
depletion preserves \(a_1=1\). The Manin constant is an
\(\ell\)-adic unit by \(\ell^2\nmid N\).

If the pullback were zero modulo \(\ell\), the homomorphism
would kill \(J_0(M)[\ell]\) and would be divisible by \(\ell\)
in \(\operatorname{Hom}_{\mathbb Q}(J_0(M),E_0)\).
Dividing would give a parametrization with a nonintegral
first coefficient of its pulled-back Néron differential.
The formal neighborhood at the infinity cusp is
\(\operatorname{Spf}(\mathbb Z_\ell[[q]])\), in the smooth
locus of the integral modular curve. The Néron mapping
property gives integrality of that coefficient for every
such map. Thus the pullback is nonzero modulo \(\ell\).
By residual irreducibility it is injective modulo \(\ell\).
This proves saturation over \(O\), using comparison of Betti
and étale cohomology. Multiplicity one and the specified
\(U_q=0\) at the raised primes identify its rational image
with the full \(f_\Sigma\)-eigenspace.
[Česnavičius, Theorem 2.11, the Manin constant at exponent at most one](https://arxiv.org/abs/1703.02951v3).
[Česnavičius–Neururer–Saha, introduction and §5, the integral cusp parameter](https://ems.press/content/serial-article-files/42577).

The further pullback to \(X_1(M)\) is again saturated on
this eigenlattice. The kernel modulo \(\ell\) of pullback
on \(H^1(X_0(M),\mathbb F_\ell)\) injects by
inflation–restriction into the characters of the diamond
quotient by \(\Gamma_1(M)\). One can use the congruence groups
before compactification, since the cohomology of the compact
curve injects into that of the open curve and the elliptic
stabilizers have order prime to \(\ell\).
On the diamond characters \(T_r\) acts by \(1+r\) for
\(r\nmid M\): the two maps from the degree-\(r\)
correspondence induce the same map on the diamond quotient,
and transfer back has degree \(1+r\).
The kernel cannot contain the residual eigenlattice. The
full residual image and Chebotarev give a good \(r\) with
\(a_r(f)\not\equiv 1+r\pmod\ell\), by taking Frobenius in
the scalar conjugacy class of \(2I\).
This proves the required injectivity modulo \(\ell\) and
saturation on the two-dimensional lattice.

Fricke is compatible with the cover \(X_1(M)\to X_0(M)\),
whose degree is \(\varphi(M)/2\). Cup product on the
saturated pulled-back lattice is multiplied by this degree.
The \(\ell\)-part of the degree is \(|\Delta|\), since
\(\ell\) occurs to exponent one in \(M\). The starting
Fricke-pairing factor at \(X_0(N)\) is \(\delta_1\) up to sign,
by pullback from \(E_0\). Hence (MD4) and (MD5) give the
fixed-determinant formula
\[
 v_\ell(\eta_B)
       =v_\ell(\delta_1)+
         \sum_{q\mid M/N}v_\ell(q^2-1).
 \tag{MD6}
\]

## The Selmer order and the abc criterion

The cotangent dual of \(B\) at the Frey augmentation is
the ordinary adjoint Selmer group fully relaxed at \(\Sigma\).
This follows by applying the universal property to square-zero
lifts over \(O/\ell^m\) and passing to the direct limit.
At \(\ell\) the condition is precisely the kernel in (BDG5);
the residual ordinary characters are distinct, so such a lift
preserves the ordinary line after a strict change of basis.
Away from \(\Sigma\cup\{\ell\}\) both the deformation and
the Selmer conditions require unramifiedness. Consequently
(BDG6) and its subsequent relaxation formula give
\[
 \operatorname{length}\Phi_B
   =\operatorname{length}S_\ell+T_\ell+s_\ell^{\rm coeff}
         +\sum_{q\mid M/N}v_\ell(q^2-1).
\]
The only potentially multiplicative tame prime that may be
absent from \(M/N\) is an additive dyadic prime, whose
\(\ell\)-valuation of \(2^2-1\) is zero. Since \(B\) is a
complete intersection its cotangent length equals
\(v_\ell(\eta_B)\). Subtracting (MD6) proves (MD1).
[Diamond, equation (4)](https://www.wstein.org/people/kiming/boston/diamond.pdf).
[Darmon–Diamond–Taylor, §2.7](https://wstein.org/edu/Fall2003/252/references/darmon-diamond-taylor/Darmon-Diamond-Taylor-Fermats_Last_Theorem.pdf).
[Ordinary local and global indices, (BDG5)–(BDG6)](ABC_FREY_ADJOINT_DEGREE_CALIBRATION.md#the-local-and-global-ordinary-comparison-at-a-multiplicative-coefficient).

In particular the entire bad-coefficient excess in the
classical Selmer product is zero. Put
\[
 Z(E)=(\delta_1)_{\{2,3\}}\prod_{\ell\ge5}|S_\ell| .
\]
This is a finite product, since (MD1) implies each local
factor is supported on the modular degree. The formula
\[
 \delta_1=Z(E)\prod_{\ell\ge5}\ell^{T_\ell+s_\ell^{\rm coeff}}
 \tag{MD7}
\]
is exact. The correction can be calculated from the triple.
Write \(x_q\) for the term of \(abc\) divisible by the odd
prime \(q\), \(x_q'=q^{-v_q(x_q)}x_q\), and use \(y_q=b\)
if \(x_q=a\), \(y_q=a\) otherwise. Put
\(r_q=x_q'/(16y_q)\) in \(\mathbb F_q^\times\) and
\(n_q=v_q(\Delta_E)\). With \(n_2=1\) if \(E\) is
potentially good at \(2\), or \(n_2=-v_2(j(E))\) otherwise,
the correction in (MD7) is the prime-to-\(6\) part of
\[
 n_2\prod_{q\mid N,\ q\ {\rm odd}}
           n_q\gcd\!\left(n_q,\frac{q-1}
                         {\operatorname{ord}_{\mathbb F_q^\times}(r_q)}\right).
 \tag{MD8}
\]
In fact (SA4) characterizes \(s_q\) by a power test on \(r_q\),
so its contribution at \(\ell\ge5\) is the \(\ell\)-part of
the displayed gcd. The coefficient correction supplies the
\(\ell=q\) part of \(n_q\) whenever \(q\ge5\).
The logarithm of (MD8) is
\(O(\log c/\log(\log c+2))\) by the divisor estimate for
the tame Tate valuations and at the coefficient primes.
Thus \(\log Z(E)=\log c+\log R+o(\log c)\).
It follows that \(Z(E)\ll_\varepsilon R^{2+\varepsilon}\)
for every \(\varepsilon>0\) is equivalent to full abc.
The same is true of the product with the old signed-support
truncations. Its upper bound follows from (MD7), while (BK1)
and the signed denominator height comparison give the matching
lower height scale. An unbounded ratio for any fixed positive
\(\varepsilon\) in either criterion is a disproof certificate.
The remaining estimate is global. All \(\ell\ge5\) factors
are now calibrated with no ordinary determinant correction
left, and the two small coefficient primes enter through the
classical degree as before.
[Signed degree height comparison, (DF7)–(DF8)](ABC_FREY_JOINT_DIVISION_FIELDS.md).
[The finite-depth signed comparison, (BK1)](ABC_FREY_CLASSICAL_ADJOINT_COMPARISON.md).
[The scalar-depth unit test, (SA4)](ABC_FREY_STRICT_ADJOINT_CONGRUENCES.md).
