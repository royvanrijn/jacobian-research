# Removing bad-prime conditions from the Frey adjoint congruence space

Let \(E/\mathbb Q\) be the Frey curve of a primitive positive triple
\(a+b=c\) with at least two odd radical primes, \(R=\operatorname{rad}(abc)\),
and \(N=N_E\). For \(\ell\ge5\), let \(A_\ell\) be the local full-level
signed Hecke order in the [division-field comparison](ABC_FREY_JOINT_DIVISION_FIELDS.md).
Write \(\lambda:A_\ell\to\mathbb Z_\ell\) for the evaluation on \(E\),
\(I_\ell=\ker\lambda\), and \(\ell^{e_\ell}\) for its Frey projector
denominator. Set
\[
 D_\ell=\mathbb Q_\ell/\mathbb Z_\ell,\qquad
 W_\ell=\operatorname{ad}^0(T_\ell E)\otimes D_\ell,\qquad
 C_\ell=\operatorname{Hom}_{\mathbb Z_\ell}(I_\ell/I_\ell^2,D_\ell).
\]
There is a canonical injection \(C_\ell\hookrightarrow
H^1(\mathbb Q,W_\ell)\), obtained from the augmented Hecke
representation. Let \(C_\ell^{\rm str}\) be the subgroup whose classes
restrict to zero at every odd \(q\mid N\), \(q\ne\ell\), on the full
decomposition group. The classes in \(C_\ell^{\rm str}\) are
unramified at **every finite prime other than \(\ell\)**, including \(2\).

For \(q\) odd, \(q\mid N\), and \(q\ne\ell\), define
\(s_{q,\ell}\) as the largest integer \(s\) for which \(E[\ell^s]\)
has scalar action of the local decomposition group at \(q\).
The main comparison is
\[
     \ell^{e_\ell}\le |C_\ell|
        \le |C_\ell^{\rm str}|\,
                 \ell^{3\sum_{q\mid N,\ q\ {\rm odd},\ q\ne\ell}s_{q,\ell}} .
 \tag{SA1}
\]
In particular, controlling the strict groups is a sufficient
route to the full abc bound. The local factor in (SA1) has
total logarithm \(o(\log c)\) across all coefficient primes,
with no restriction on the number of odd conductor primes.

## The conormal module and its Galois classes

The algebra \(A_\ell[1/\ell]\) is a product of fields, with the
\(E\)-factor equal to \(\mathbb Q_\ell\). Consequently \(I_\ell/I_\ell^2\)
is finite over \(\mathbb Z_\ell\), and
\(\lambda(\operatorname{Ann}_{A_\ell}I_\ell)=(\ell^{e_\ell})\).
Indeed in the generic product the annihilator consists precisely
of the \(E\)-factor. If \(d_i\) are the elementary divisors
of \(I_\ell/I_\ell^2\) and \(x_i\in I_\ell\) are corresponding
lifts, the \(x_i\) generate \(I_\ell\) by Nakayama. Write
\(d_i x_i=\sum_j a_{ij}x_j\) with \(a_{ij}\in I_\ell\). The
adjugate shows that \(\det(d_i\delta_{ij}-a_{ij})\) annihilates
\(I_\ell\), and its augmentation is \(\prod_i d_i\).
Thus \(e_\ell\le\operatorname{length}(I_\ell/I_\ell^2)\).
Also \(\ell^{e_\ell}e_f\) annihilates \(I_\ell\), so on its conormal
module \(\ell^{e_\ell}\) acts by zero. This gives both
\(\ell^{e_\ell}\le |C_\ell|\) and
\(\ell^{e_\ell}C_\ell=0\).
[Congruence and relative cotangent modules, Böckle–Khare–Manning](https://arxiv.org/abs/2108.09729).

Choose the integral realization of the \(A_\ell\)-representation
whose augmentation is the Tate module of \(E\). For
\(d\in C_\ell\), apply the ring map
\[
 a\longmapsto \lambda(a)+\epsilon d(a-\lambda(a)),
 \qquad \epsilon^2=0,\quad \epsilon D_\ell\text{ the square-zero ideal},
\]
to its entries. This gives a continuous first-order deformation
with fixed determinant and hence a cocycle with values in \(W_\ell\).
If the cocycle is a coboundary then the deformation is conjugate
to the constant lift and all trace derivatives vanish. The traces
generate \(A_\ell\), so \(d=0\). This proves injectivity. The
deformation is unramified at every finite prime away from \(N\ell\)
because the same is true on every characteristic-zero factor of
the signed Hecke order. The injectivity of the generic product
ensures that each local matrix identity holding in all factors
holds in the coefficient ring, and hence on these deformations.
[Integral trace construction and full Frey residual image](ABC_FREY_JOINT_DIVISION_FIELDS.md#the-image-over-the-hecke-order).

## The signed Steinberg restriction at an odd prime

Fix \(q\) odd with \(q\mid N\) and \(q\ne\ell\). At \(q\) all
characteristic-zero constituents are Steinberg with the same
unramified quadratic character. Normalize Frobenius by this
sign, and choose a Tate basis and a generator of the
pro-\(\ell\) tame inertia such that
\[
 T=1+nU,\quad U=e_{12},\quad
 F=\begin{pmatrix}q&u\\0&1\end{pmatrix},\qquad n=v_q(\Delta_E).
\]
Then \(s_{q,\ell}=\min(v_\ell(n),v_\ell(q-1),v_\ell(u))\), allowing
\(v_\ell(0)=+\infty\). This is the depth of local projective
triviality, independent of the choices of tame generator and
Frobenius lift. Every deformation coming from \(C_\ell\)
has the same traces and determinants on \(F,T\) and satisfies
\[
     (\widetilde F-q)(\widetilde T-1)=0,\qquad
     (\widetilde T-1)(\widetilde F-1)=0.
 \tag{SA2}
\]
These identities follow in each signed Steinberg constituent by
the local Weil–Deligne description. All these constituents
factor through the quotient with inertia \(\mathbb Z_\ell\),
so the specified inertia and Frobenius determine the local
restriction.
[Carayol, §0.7, Theorem A](https://numdam.org/item/10.24033/asens.1512.pdf).
[Böckle–Khare–Manning, Definition 7.7 and Lemma 7.11](https://arxiv.org/abs/1910.08507).

Write \(\widetilde T=T+\epsilon X\),
\(\widetilde F=F+\epsilon Y\) with trace-zero \(X,Y\). The inertia
determinant gives \(nX_{21}=0\). A change of basis by \(1+\epsilon Z\)
adds the commutator \(n[Z,U]\) to \(X\). Since \(D_\ell\)
is divisible and \(2\) is a unit, this kills the diagonal and
upper entry of \(X\), leaving \(X=x e_{21}\). Write
\(Y=\left(\begin{smallmatrix}a&b\\c&-a\end{smallmatrix}\right)\).
Linearizing (SA2) and the two determinants now gives
\[
       nx=(q-1)x=ux=0,\qquad na=nc=0,\qquad
                          (q-1)a+uc=0.
 \tag{SA3}
\]
The remaining transformations preserve \(X\) precisely when
\(nZ_{21}=nZ_{11}=0\). The free entry \(Z_{12}\) kills \(b\)
because \(q-1\) acts surjectively on \(D_\ell\).
The action on \((a,c)\) is addition of a multiple of
\((-u,q-1)\) with multiplier in \(D_\ell[n]\).
The \(x\)-coordinate is exactly \(D_\ell[\ell^{s_{q,\ell}}]\).
With \(t=v_\ell(n)\), the kernel of \((q-1)a+uc\) on
\(D_\ell[n]^2\) has order \(\ell^{t+s_{q,\ell}}\), and the
remaining transformations on this kernel have image of order
\(\ell^{t-s_{q,\ell}}\). Thus the normalized solution space
modulo conjugacy has order \(\ell^{3s_{q,\ell}}\), bounding the
image of \(C_\ell\). In the case \(s_{q,\ell}=0\) every class
already restricts to zero at this prime. In particular, residual
noncentrality of Frobenius suffices even when inertia on
\(E[\ell]\) is trivial. This calculation agrees numerically
with the \(3s_{q,\ell}\) derived local invariant computed by
Böckle–Khare–Manning. Their resulting local Wiles defect is
\(2s_{q,\ell}\), since it subtracts a separate congruence term.
[Corollaries 5.15 and 5.17](https://arxiv.org/abs/2108.09729).

There is an explicit residual test for \(s_{q,\ell}\). Let \(x_0\)
be the letter among \(a,b,c\) divisible by \(q\), let
\(x_0'=q^{-v_q(x_0)}x_0\), and take \(y_0=b\) if \(x_0=a\),
or \(y_0=a\) otherwise. For \(\ell^s\mid \gcd(n,q-1)\),
one has \(s\le s_{q,\ell}\) exactly when
\[
                  \left(\frac{x_0'}{16 y_0}\right)^{(q-1)/\ell^s}
                              =1\quad\text{in }\mathbb F_q^\times .
 \tag{SA4}
\]
In fact the Tate parameter divided by \(q^n\) is congruent to
\(q^{-n}\Delta/c_4^3\bmod q\), from the first term of the Tate
\(j\)-series. For the standard Frey model the latter residue is
\((x_0'/(16y_0))^2\). An unramified quadratic twist does not
change \(j\) or the Tate parameter. Since \(\ell^s\mid q-1\),
the cyclotomic roots of order \(\ell^s\) are already present
locally. With \(\ell^s\mid n\), triviality of the local Kummer
class is exactly the indicated power condition on the unit.
The square factor has no effect for odd \(\ell\).
[Tate uniformization and \(j\)-series](https://people.math.harvard.edu/~ctm/home/text/others/tate/tate-unif/tate-unif.pdf).
[Earlier Tate inertia and Frobenius coordinates](ABC_CUBIC_ADJOINT_TAMAGAWA_FACTORS.md).

## Unramified adjoint classes at the prime two

For \(\ell\ge5\), a two-dimensional \(\ell\)-adic representation of
\(G_{\mathbb Q_2}\) with finite inertia image has inertia of order
prime to \(\ell\). To see this, let \(P\) be the image of wild
inertia, a \(2\)-group. The quotient of the inertia image by
\(P\) is cyclic. If there is a nontrivial \(\ell\)-Sylow subgroup,
it is cyclic and maps injectively to this quotient. Up to
inertial conjugacy, Frobenius carries a generator \(g\) of it
to \(g^2\), by the tame Frobenius action and conjugacy of Sylow
subgroups. The multiset of the two eigenvalues of \(g\) is then
stable under squaring. A root of \(\ell\)-power order in a
multiset of size two can have order only dividing
\(2^2-1=3\), a contradiction. In addition, reduction on a stable
lattice is injective on a finite inertia group of order prime
to \(\ell\), since the reduction kernel is pro-\(\ell\).

If \(E\) has potentially good reduction and its inertia is
nonscalar, some inertial element has distinct eigenvalues on
the Frey residual representation by this prime-to-\(\ell\)
property. No competitor can be a special representation at
\(2\), since its inertial eigenvalues would always coincide.
All competitors therefore have finite inertia of order prime
to \(\ell\), by local-global compatibility. The inertia image
over \(A_\ell\) embeds into the product of these finite groups.
Its first cohomology with divisible \(\ell\)-primary coefficients
vanishes, proving unramifiedness of the adjoint classes at \(2\).

If \(E\) has potentially good reduction and its inertia is
scalar, the determinant shows that its restriction to inertia
is a quadratic character. It extends to a character of
\(G_{\mathbb Q_2}\) since it is Frobenius invariant and the
unramified quotient is procyclic. Make this same twist on
all the local representations. The Frey curve now has good
reduction by the good-reduction criterion, and still has full
rational \(2\)-torsion. Its reduction is ordinary: for good
supersingular reduction the special-fiber group has odd order,
and the height-two formal group has no rational nontrivial
\(2\)-torsion. Indeed in \([2](T)\) the coefficients of
\(T^2,T^3\) are even, so \(2T\) is the unique term of lowest
valuation at \(0\ne T\in 2\mathbb Z_2\). Hence the twisted
Frey Frobenius trace is \(1\) or \(-1\). A special competitor,
after the same twist, would have unramified twisting character,
since its residual inertial eigenvalues would both be one.
Its Frobenius trace would then be \(3\) or \(-3\), which cannot
agree modulo \(\ell\ge5\). All remaining competitors have finite
prime-to-\(\ell\) inertia with trivial reduction after twisting,
so their inertia is trivial. This proves the desired adjoint
unramifiedness. It also covers good reduction without twisting.
[Serre–Tate, good and potential good reduction criteria](https://wstein.org/papers/bib/Serre-Tate-Good_Reduction_of_Abelian_Varieties.pdf).

Finally, suppose \(E\) has potentially multiplicative reduction.
Twist locally by a splitting quadratic character and use Tate
inertia, with the baseline matrices \(T,F\) of the previous
section at \(q=2\). All competitors after this twist are either
unramified or special with unramified character. For finite
inertia this follows from injectivity of reduction and the fact
that their residual inertial eigenvalues are one. For special
inertia the scalar semisimplification has order at most two
and is trivial on inertia after this twist, again by its
residual eigenvalues. The special competitors have the same
Frobenius sign as the twisted \(E\), since \(+3\) and \(-3\)
are distinct modulo \(\ell\ge5\). Both unramified and same-sign
special competitors satisfy (SA2) at \(q=2\), even if an
unramified competitor has variable Frobenius trace. The inertia
trace and determinant are still two and one. Normalize \(X\)
again to \(x e_{21}\) using these inertia constraints. The
lower-left entry of the first equation in (SA2) yields
\((1-2)x=0\), independently of the variation of the Frobenius
trace. Thus the restricted class on inertia is zero.

All the twists here are local and preserve the adjoint
representation. The outcome in all reduction types is that
\(C_\ell\) is unramified at \(2\), including when its local class
on the decomposition group need not vanish.

## Strict classes and an equivariant class-field group

The restriction kernel defining \(C_\ell^{\rm str}\) and the
odd-prime calculation (SA3) establish (SA1). Moreover
\[
 \prod_{\ell\ge5}\ell^{3\sum_{q\mid N,\ q\ {\rm odd},\,q\ne\ell}s_{q,\ell}}
       \ \le\ \prod_{q\mid N,\ q\ {\rm odd}}\gcd(v_q(\Delta_E),q-1)^3
       \ \le\ \exp\!\left(O\!\left(\frac{\log c}{\log(\log c+2)}\right)\right)
 \tag{SA5}
\]
Indeed \(v_q(\Delta_E)=2v_q(abc)\) at
odd \(q\), and the product of the valuations times the factor
\(2^{\omega(abc)}\) is bounded by a fixed power of the divisor
count of \(abc\). The usual elementary divisor estimate and
\(abc<c^3\) give (SA5).
[The same divisor bound in the Shimura degree comparison](ABC_CUBIC_ADJOINT_SHIMURA_PROJECTORS.md#uniform-shimura-transfer).

Let \(K_\ell\) be the field cut out by the projective
\(\ell\)-adic representation of \(E\), and let \(X_\ell\) be the
Galois group over \(K_\ell\) of its maximal abelian
pro-\(\ell\) extension unramified away from primes over \(\ell\),
with complete splitting at primes over every odd \(q\mid N\),
\(q\ne\ell\). These conditions use the respective inertia
and decomposition subgroups over \(K_\ell\). Conjugation
by \(\Gamma_\ell=\operatorname{Gal}(K_\ell/\mathbb Q)\)
acts on \(X_\ell\). Restriction gives
\[
                C_\ell^{\rm str}\hookrightarrow
             \operatorname{Hom}_{\rm cont,\Gamma_\ell}(X_\ell,W_\ell).
 \tag{SA6}
\]
To check injectivity, \(\Gamma_\ell=\operatorname{PGL}_2(\mathbb Z_\ell)\)
by (DF3) projected to \(E\). Inflation injects its \(H^1\) into that
of \(\operatorname{GL}_2(\mathbb Z_\ell)\) with these coefficients.
One has
\(H^1(\operatorname{GL}_2(\mathbb Z_\ell), W_\ell)=0\).
Modulo \(\ell\), inflation from
\(\operatorname{GL}_2(\mathbb F_\ell)\) has zero \(H^1\) by the
calculation in (DF3). Restriction to the first principal
congruence subgroup is a homomorphism factoring through its
first layer \(M_2(\mathbb F_\ell)\), since every element of the
second principal congruence subgroup is an \(\ell\)-th power
in the first, by the logarithm. The restriction is equivariant,
so it is a scalar map on the trace-zero summand and zero on the
scalar summand. On the \(\ell\)-th power of \(1+e_{12}\) it is
zero by the vanishing of the cyclic norm operator on the
residual adjoint. This power represents \(e_{12}\) in the first
layer, so restriction is zero. This proves \(H^1=0\) with
residual coefficients, and the coefficient exact sequences
give zero for \(W_\ell[\ell^j]\) for every \(j\). Continuous
cocycles into \(W_\ell\) have finite image, proving the
desired assertion. Inflation–restriction now gives (SA6).
The kernel of the action on \(W_\ell\) contains
\(G_{K_\ell}\), so the restricted cocycles are homomorphisms,
and the local conditions ensure the required factorization.
[Full Galois image and the finite residual adjoint calculation](ABC_FREY_JOINT_DIVISION_FIELDS.md#the-image-over-the-hecke-order).

For a finite class-field realization of the subgroup killed
by \(\ell^m\), let \(K_{\ell,j}\) be the projective
\(\ell^j\)-division field. Take \(j=m+\beta_\ell\), where
\(\beta_\ell\) is the maximum of zero, the numbers
\[
          v_\ell(v_q(\Delta_E))+v_\ell(q-1)
                \qquad(q\mid N,\ q\ {\rm odd},\ q\ne\ell),
\]
and \(v_\ell(-v_2(j(E)))\) if the reduction at \(2\) is
potentially multiplicative. A derivation of order dividing
\(\ell^m\) has a representative with values in \(W_\ell[\ell^m]\).
At an odd prime in this maximum, local strictness writes that
representative as a coboundary from a matrix \(Z\). The
diagonal and upper entries of \([Z,T]\) show that
\(Z_{11},Z_{21}\) are killed by \(\ell^{m+v_\ell(n)}\);
the upper entry of \([Z,F]\) then shows that \(Z_{12}\)
is killed by \(\ell^{m+v_\ell(n)+v_\ell(q-1)}\).
At \(2\) with finite prime-to-\(\ell\) inertia, averaging
gives a coboundary already in \(W_\ell[\ell^m]\).
For the potentially multiplicative case at \(2\), the inertia
gauge used above can take \(Z_{12}=0\) and the two other entries
killed by \(\ell^{m+v_\ell(-v_2(j(E)))}\). Thus restriction
to \(G_{K_{\ell,j}}\) already has the required unramifiedness and
splitting. The resulting homomorphisms are injective on the
indicated subgroup, since restriction further to the infinite
tower is injective.

At a prime \(v\mid\ell\) of \(K_{\ell,j}\) with ramification
index \(e_v\), the local logarithm is an isomorphism on
principal units of depth greater than \(e_v/(\ell-1)\).
Therefore all local characters of exponent dividing \(\ell^m\)
vanish on the principal units of depth
\(e_v m+\lfloor e_v/(\ell-1)\rfloor+1\). Global reciprocity
realizes the homomorphisms in the adjoint part of the ray class
group for the product of these powers of the primes above
\(\ell\), modulo the classes of primes required to split.
This finite modulus gives a concrete class-field target for
each depth \(m\), with no ramification modulus at other primes.

There is also a finite Selmer group containing
\(C_\ell^{\rm str}\). At \(\ell\), each weight-two constituent
of this Hecke order comes from a modular Jacobian with good
or semistable reduction, since \(v_\ell(N)\le1\). Thus
\(A_\ell^2\) is a stable \(\mathbb Z_\ell\)-lattice in a semistable
representation with the two Hodge–Tate weights of \(T_\ell E\).
A derivation of exact order \(\ell^m\) takes values in
\(D_\ell[\ell^m]\) and surjects onto that group. Identifying it
with \(\mathbb Z/\ell^m\), the square-zero ring map is surjective
onto \((\mathbb Z/\ell^m)[\epsilon]/(\epsilon^2)\), so the
associated first-order representation is a quotient of that lattice. The local classes
satisfy the torsion semistable condition of
Khare–Patrikis. Their Theorems 1.1 and 3.14 give finiteness
for the larger Selmer group using their semistable condition
at \(\ell\) and allowing all local classes at the other bad
primes. Here the residual cyclotomic restriction is irreducible
by full image. In the exceptional \(\ell=5\) hypothesis their
projective field does not contain \(\mathbb Q(\mu_5)\),
because a scalar residual matrix \(2I\) has determinant
\(4\ne1\bmod5\). The effective bound there still contains
the original Hecke cotangent and congruence terms. It does
not furnish an independent sharp radical-logarithmic bound.
[Conrad–Stein, §7.2, the modular Jacobian at a prime exactly dividing the level](https://wstein.org/papers/compgrp/compgrp.pdf).
[Khare–Patrikis, Definition 3.5, Theorem 1.1 and Theorem 3.14](https://arxiv.org/abs/2611.33361).

Let \(\delta_1\) be the classical optimal modular degree.
A uniform bound
\[
       (\delta_1)_{\{2,3\}}\prod_{\ell\ge5}|C_\ell^{\rm str}|
                    \ll_\varepsilon R^{2+\varepsilon}
                     \qquad(\varepsilon>0)
 \tag{SA7}
\]
would prove full abc. To see this without computing the joint
field, use the signed toric projector denominator \(B\) of
(HN1). For \(\ell\ge5\), the common signed Hecke order acts
integrally on each toric lattice, so \(v_\ell(B)\le e_\ell\)
by (DF6). The transfer estimate (SP9) and \(B\mid\delta_{D,M}\)
give \(B_{\{2,3\}}\le (\delta_1)_{\{2,3\}}\exp(O(1+\omega(D)))\).
Combining these with (SA1) and (SA5) bounds \(B\) by the left
side of (SA7) times \(\exp(o(\log c))\). The comparison
\(\log B=\log c+\log R+o(\log c)\) then gives the sharp abc
bound, adjoining the elementary low-support cases. This uses
neither the exact field intersection (DF4) nor the
linear disjointness of the joint fields (DF8).
[Signed projector comparison](ABC_FREY_FULL_LEVEL_CONGRUENCES.md).
[Exact transfer and low-support reduction](ABC_CUBIC_ADJOINT_SHIMURA_PROJECTORS.md).
If abc fails, the groups in (SA7) therefore have a sequence
exceeding the sharp bound. On the other hand,
\(\ell^{e_\ell}\) annihilates \(C_\ell^{\rm str}\) by the
conormal calculation. A converse estimate for their *order*
requires control of the number and size of their cyclic
summands, in addition to the largest exponent. This issue
already occurs for finite flat trace orders. For integers
\(b,k\ge1\), take
\[
 A=\{(a_0,\ldots,a_k)\in\mathbb Z_\ell^{k+1}:
            a_i\equiv a_0\pmod{\ell^b}\ \ (1\le i\le k)\},
 \qquad \lambda(a)=a_0.
\]
Its first projector denominator is \(\ell^b\), while its
conormal dual is \((\mathbb Z/\ell^b)^k\), by coordinatewise
multiplication in the augmentation ideal. The determinant-restricted
general linear group of this order has full residual image and
generates the order by traces: it contains determinant-one
matrices with each prescribed trace. Thus those algebra and
image hypotheses alone place no uniform bound on the excess
length of the conormal module.

In the actual signed Hecke order that excess length is its
Wiles defect at the Frey augmentation. A sublinear bound in
\(\log c\) for the sum of these defects weighted by
\(\log\ell\) would make (SA7) equivalent to abc using (DF7).
Böckle–Khare–Manning's local-global formula applies to a
squarefree Shimura level, finite flat at a coefficient prime
outside the level, with residual minimality at the other level
primes specified in their theorem. It does not give a uniform
formula for the signed classical orders here, which allow
\(\ell\mid N\) and higher conductor power at \(2\). The sufficient
implication (SA7) uses no such defect bound.
[Theorem 1.1](https://arxiv.org/abs/2108.09729).

The remaining arithmetic problem is to bound the equivariant
classes in (SA6) that come from the signed Hecke order, with
semistability at the coefficient prime, across all \(\ell\).
The finite ray class realization specifies the required splitting
and allowed conductor. The outstanding estimates are the size of
the semistable adjoint part and the contribution at \(\ell=2,3\) to
the modular degree.
