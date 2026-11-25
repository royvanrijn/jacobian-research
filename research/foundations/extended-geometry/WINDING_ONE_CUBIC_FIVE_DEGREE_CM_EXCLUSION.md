# Exclusion of the degree-five leading coefficient for winding-one septics

Use the sphere functional \(\Lambda\) and the winding-one presentation
\(f=x^{-1}(C+Bt+Dt^2+Et^3)\) over a characteristic-zero field.

**Theorem.** If \(f\) is a regular winding-one cubic of least sphere
degree seven, \(E(0)\ne0\), and \(\deg E=5\), then
\(\Lambda(f^m)\ne0\) for infinitely many positive orders \(m\).

By the [boundary theorem](WINDING_ONE_CUBIC_LOGARITHMIC_BOUNDARY.md),
it suffices to consider \(E=e(x-r)^4(x-s)\), with \(r,s\ne0\) and
\(r\ne s\), and its boundary divisibilities. Assume eventual
purity. The [marker and ratio argument](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_MARKERS.md)
makes \(\mu=r/s\) algebraic. Name the roots by extending the
coefficient field, then specialize to a number field over
\(\mathbb Q(\mu)\), preserving the required nonzero data and all the
putative vanishing moment identities. Rescaling the auxiliary
coordinate \(W\) allows \(e=1\) in the auxiliary curve without
changing the zero Cartier tests.

On the genus-two curve
\[
                 W^3=x^2(x-r)(x-s)
\]
put \(\eta=dx/((x-r)W)\). At all sufficiently large split primes
\(p=3k+1\), Cartier kills \(\eta\): its only possible Cartier
coefficient lies at \(x^{p-1}\), and the unmarked prime row is a unit
times this coefficient. The same computation gives
\[
 H_p'(\mu)=0,\qquad H_p(\mu)\ne0,\qquad
 H_p(Z)=(-1)^k\sum_{j=0}^k\binom{k}{j}^{\!2}Z^j .          \tag{1}
\]
The nonvanishing is the ordinary-point uniqueness calculation for
the hypergeometric equation in the marker and ratio argument. By the
[inert Cartier and CM theorem](WINDING_ONE_CUBIC_FIVE_DEGREE_INERT_CARTIER.md#theorem-complex-multiplication-restriction),
the genus-two curve has an elliptic quotient with CM field
\(\mathbb Q(\sqrt{-3})\). Its elliptic quotients below are ordinary
at the split primes in (1).

## The two elliptic quotients and their short models

Set \(c=r-s\), choose \(c^{1/3}\), and write
\[
 u=W/(x c^{1/3}),\qquad V=(x-rs/x)/c,\qquad
 a=2(r+s)/(r-s).
\]
Then
\[
 V^2=u^6+a u^3+1,qquad x=\frac{c}{2}(u^3+V)+\frac{r+s}{2}.        \tag{2}
\]
The involutions \((u,V)\mapsto(1/u,\pm V/u^3)\) give elliptic
quotients \(E_\pm\). The two involutions differ by the hyperelliptic
involution. The elliptic involution in the inert Cartier theorem is
one of these up to deck conjugacy, as its action on \(W/x\) is
proportional to \(x/W\). With \(z=u+1/u\) and
\(w_\pm=V(u\pm1)/u^2\), the quotient equations are
\[
              w_\pm^2=(z\pm2)(z^3-3z+a).
\]
For the plus sign take \(n=a-2\), \(z=-2+n/X\); for the minus
sign take \(n=a+2\), \(z=2+n/X\). In each case take
\(Y=w_\pm X^2/n\) and \(Z=X+3\). These give short Weierstrass
models \(Y^2=Z^3+A_\pm Z+B_\pm\), where
\[
 A_+=-6a-15,\quad B_+=a^2+14a+22,
 \qquad A_-=6a-15,\quad B_-=a^2-14a+22.                  \tag{3}
\]

There is a cyclic isogeny of degree three between these quotients.
Write \(J\) for the genus-two Jacobian, \(i_\pm:E_\pm\to J\) for
the pullbacks and \(p_\pm:J\to E_\pm\) for the norms. Write \(s_0\)
for the plus involution on \(J\), and \(g\) for the order-three deck
transformation. Since the hyperelliptic involution acts by minus one,
\(i_-p_-=1-s_0\) and \(p_\pm i_\pm=[2]\). The relation
\(1+g+g^{-1}=0\) holds on \(J\) because the deck quotient has genus
zero. The morphism \(\varphi=p_-g i_+\) has adjoint \(p_+g^{-1}i_-\).
Since \(s_0g s_0=g^{-1}\), we have \(g^{-1}s_0g=g s_0\) and
\(p_+ g i_+=p_+g^{-1}i_+=[-1]\). Thus
\[\varphi^\vee\varphi=p_+(1-gs_0)i_+=[3].\]
This gives the claimed degree and shows that both quotients have
CM by orders in \(\mathbb Q(\sqrt{-3})\).

## The Cartier condition determines the parameter over either j-invariant

Up to a nonzero scalar and an exact differential, (2) gives
\(\eta=(u^{-2}-u)\,du/V\). Indeed \(dx/(x(W/x))\) is proportional
to \(u\,du/V\), and \((1-u^3)^2-V^2=-(a+2)u^3\) rationalizes
the factor \(x-r\). Adding its transform by the plus involution gives
the pullback of \(-(z^2+z)\,dz/w_+\), up to the same scalar. For the
minus involution the corresponding form is \((z^2-z-2)\,dz/w_-\).
These identities follow from \(dz/w_+=(u-1)du/V\) and
\(dz/w_-=(u+1)du/V\).

In the plus cubic model before translating \(X\),
\(d(Y/X)=(X/2+3n/X-n^2/X^2)\,dX/Y\). This reduces
\(-(z^2+z)dX/Y\) to \(-(X/2+2)dX/Y\) modulo exact forms. In
the minus model \(d(Y/X)=(X/2-3n/X-n^2/X^2)\,dX/Y\), which
reduces \((z^2-z-2)dX/Y\) to \((X/2)dX/Y\). At \(Z=X+3\),
the two elliptic de Rham lines are consequently
\[
               [(Z+1)dZ/Y],\qquad [(Z-3)dZ/Y].             \tag{4}
\]

Cartier kills these elliptic classes at the split primes in (1).
This follows by pullback to the genus-two curve, functoriality of
Cartier, vanishing on exact differentials, and injectivity of the
pullback on holomorphic forms. The second-kind representatives
\(Z\,dZ/Y\) have just a double pole at infinity, so Cartier maps
these representatives to holomorphic forms at the good odd primes.

For a CM elliptic curve its de Rham cohomology splits into lines for
the two CM embeddings, one of which is the Hodge line. At all
sufficiently large split primes the non-Hodge line is killed by
Cartier. To see this, choose a nonrational CM endomorphism defined
over a finite number field extension. Its distinct eigenvalues
reduce to elements of \(\mathbb F_p\) at the split primes. Cartier
commutes with pullback by the endomorphism and is inverse Frobenius
semilinear, so only the Hodge eigenvalue maps to the holomorphic
eigenline. To reduce these eigenclasses, choose their representatives
as linear combinations of \(dZ/Y\) and \(Z\,dZ/Y\), and retain the
exact rational derivative relations for the endomorphism action over
the number field. Those identities reduce at all but finitely many
primes, and Cartier kills the exact terms. Ordinary reduction makes
Cartier nonzero on the Hodge line. Thus
each line in (4) is the non-Hodge CM eigenline in characteristic
zero. A nonzero coefficient on the Hodge line would remain nonzero
at all but finitely many primes, contradicting (1).

Let \(a'\) be any algebraic conjugate of \(a\) fixing \(j_+\).
The same eigenline criterion holds for \(a'\). An isomorphism of
the short models over the algebraic closure is given by
\((Z,Y)\mapsto(u_0^2Z,u_0^3Y)\) for some \(u_0\ne0\). It must
identify the non-Hodge CM eigenlines, because it preserves the
Hodge line and conjugates the endomorphism algebras. This means
\([u_0Z\,dZ/Y+u_0^{-1}\,dZ/Y]\) spans \([(Z+1)dZ/Y]\).
The classes \([dZ/Y]\) and \([Z\,dZ/Y]\) are independent. An
exact combination with nonzero second coefficient would give a
rational primitive with a unique simple pole on the elliptic curve,
while a nonzero holomorphic form is also nonexact. Thus \(u_0^2=1\)
and the coefficient \(A_+\) in (3) gives \(a'=a\). For \(j_-\)
the same argument uses the nonzero offset \(-3\) in (4) and the
coefficient \(A_-\). This proves
\[
                     \mathbb Q(a)=\mathbb Q(j_+)=\mathbb Q(j_-).
                                                                    \tag{5}
\]

## Comparing the CM conductors

Suppose first that the conductors are equal. The two \(j\)-invariants
are then Galois conjugate by the main theorem of complex
multiplication. Choose an automorphism \(\sigma\) of a finite Galois
extension with \(j_+(a^\sigma)=j_-(a)\). Compare the short plus
model at \(a^\sigma\) with the short minus model at \(a\), writing
\(Z_+=u_0^2Z_-\) and \(Y_+=u_0^3Y_-\). The two non-Hodge lines
(4) must coincide under this isomorphism. This makes the lines of
\((Z_-+u_0^{-2})dZ_-/Y_-\) and \((Z_--3)dZ_-/Y_-\) equal.
Therefore \(u_0^{-2}=-3\) and \(u_0^4=1/9\). By (3),
\[
          -6a^\sigma-15=(6a-15)/9,
          \qquad \sigma(a+2)=-(a+2)/9.
\]
An automorphism of a finite Galois extension has finite order, so
this forces \(a=-2\). That is impossible since \(r\ne0\).

If the conductors differ, they are \(f_0\) and \(3f_0\).
Indeed transporting endomorphisms through a cyclic \(3\)-isogeny
shows that three times either order is contained in the other.
For instance an endomorphism \(b\) of the source induces
\(\varphi b\varphi^\vee\) on the target, corresponding to \(3b\)
under the common CM field identification. In conductor coordinates
these two inclusions force the indicated ratio. The degree of a CM
\(j\)-invariant over \(\mathbb Q\) is the order class number. For
\(\mathbb Q(\sqrt{-3})\), the formula for the order class number
gives \(h(3f_0)=3h(f_0)\) when \(f_0>1\). The Euler factor at the
ramified prime three is one, and both orders then have the same
unit index. The unequal degrees contradict (5). The class
field facts used in both cases are the ideal-class action, the
irreducibility of the class polynomial and the order class-number
formula. [Deuring, §§7, 9, 12 and 15–16](https://arxiv.org/pdf/2006.04531)

In the remaining pair of conductors \(1,3\), the maximal order
quotient has \(j=0\). On that short model the automorphism
\(Z\mapsto\zeta_3 Z,\ Y\mapsto Y\) puts \([dZ/Y]\) and
\([Z\,dZ/Y]\) in different eigenspaces. Its non-Hodge CM line is
\([Z\,dZ/Y]\). Each line in (4) has a nonzero coefficient on
\([dZ/Y]\), contradicting the eigenline condition for the maximal
order quotient. This exhausts the possibilities and proves the
theorem.
