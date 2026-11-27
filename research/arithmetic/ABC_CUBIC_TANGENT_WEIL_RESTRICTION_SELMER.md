# The common cubic-twist Selmer group as one isogeny descent

The restriction of scalars from the pure cubic field
admits a degree-\(27\) isogeny from the three
Mordell twists. Its degree-\(3\) factor has
Selmer group equal to the full intersection of
the three common quadratic-algebra Selmer
groups in the tangent family. A second factor
has a zero-sum condition on the three dual
Selmer groups. The local equality at \(3\)
uses a point over the pure cubic completion.

## The isogeny and its three-primary kernel

Let \(b\in\mathbb Q^\times\), \(a\in\mathbb Q^\times\)
a noncube, \(r=\sqrt[3]a\) the real cube root,
\(L=\mathbb Q(r)\), and \(M=L(\zeta_3)\).
Write
\[
 E_j:y^2=x^3+a^{2j}b,\quad
 E_j':y^2=x^3-27a^{2j}b,\qquad j=0,1,2.
\]
Write \(\phi_j:E_j\to E_j'\)
for the usual three-isogeny and
\(\widehat\phi_j:E_j'\to E_j\) for its dual.
The common kernel modules are
\(W=\mathbb F_3(\chi_b)\) for \(\phi_j\) and
\(W^*=\mathbb F_3(\chi_{-3b})\) for its dual.
Put \(\mathscr A=\prod_j E_j\) and
\(\mathscr R=\operatorname{Res}_{L/\mathbb Q}(E_{0,L})\).

**Theorem WR1.** There are isogenies
\(\Phi:\mathscr A\to\mathscr R\) and
\(\Psi:\mathscr R\to\mathscr A\) over \(\mathbb Q\)
whose two compositions are multiplication by three.
Their ranks and their Selmer groups at primes
different from three decompose accordingly.
The kernel of \(\Phi\) has order \(27\) and is
an extension
\[
          0\longrightarrow W^2\longrightarrow
          \ker\Phi\longrightarrow W^*\longrightarrow0 .
\]
In the basis below this extension has the two
Kummer classes \([-4b]\) and \([a]\), up to
the choice of roots fixing the common generator
of \(\mu_3\).
The Galois field of this kernel is exactly
\(\mathbb Q(\sqrt b,\zeta_3,\sqrt[3]a,\sqrt[3]{-4b})\).

For the construction, use
\(\psi_j:E_{j,L}\to E_{0,L}\) given by
\((x,y)\mapsto(r^{-2j}x,a^{-j}y)\).
The three coordinates of \(\mathscr R_M\)
are indexed by \(r\mapsto\zeta_3^i r\).
In these coordinates take
\[
 \Phi_i(P)=\sum_{j=0}^2 [\zeta_3^{ij}]\psi_j(P_j),
 \qquad
 \Psi_j(Q)=\psi_j^{-1}\sum_{i=0}^2[\zeta_3^{-ij}]Q_i .
\]
The conjugates of \(\psi_j\) have precisely
the factors \([\zeta_3^{ij}]\), so the first map
descends by the restriction-of-scalars property.
The second descends as the sum of the conjugate
inverse maps. Matrix multiplication gives the
two compositions. The eigenpoint rank
decomposition over the cyclotomic field is
given by
[Mikić–Najman, Lemma 3.1](https://web.math.pmf.unizg.hr/glasnik/50.2/50%282%29-05.pdf).
Here the two maps descend to the restriction
from the nonnormal pure cubic field, and the
calculation also determines the kernel at three.
In particular
\(\operatorname{rank}E_0(L)=\sum_j\operatorname{rank}E_j(\mathbb Q)\).
At any \(\ell\ne3\) the mutual inverse isogeny
identifications on \(\ell\)-power torsion and
local Kummer images give
\(\operatorname{Sel}_{\ell^n}(E_0/L)
 \cong\bigoplus_j\operatorname{Sel}_{\ell^n}(E_j/\mathbb Q)\)
by restriction-of-scalars cohomology.

Over a splitting field put \(P=(0,\sqrt b)\)
and \(Q=(c,\lambda\sqrt b)\), with
\(c^3=-4b\) and \(\lambda=1+2\zeta_3\).
The three-division polynomial and addition formula
give \(P,Q\in E_0[3]\) and
\([\zeta_3]Q-Q=P\). The three-isogeny maps \(Q\)
to \((0,3\lambda\sqrt b)\) on \(E_0'\).
Write the aligned \(j\)-th coordinate of a kernel
point as \(U_jP+V_jQ\). Since the kernel is killed
by three, the three row equations are equivalent
to \(V_0=V_1=V_2=v\) and \(\sum_j U_j=0\).
The subgroup \(v=0\) is the zero-sum plane
in the product of the three \(W\)'s. The quotient
is the common diagonal in the three \(\phi_j(E_j[3])\)'s,
which identifies it with \(W^*\).

For a Galois element \(g\), write
\(\chi=\chi_b(g)\), \(g\zeta_3=\zeta_3^e\),
\(g c=\zeta_3^u c\), and \(g r=\zeta_3^z r\).
The \(E_0\) coordinates have
\(gP=\chi P\) and \(gQ=\chi e(Q+uP)\).
On the aligned \(E_j\) point there is also
the factor \([\zeta_3^{-jz}]\).
Thus
\[
 v'=\chi e\,v,\qquad U_j'=\chi U_j+\chi e(u-jz)v .
\]
Taking kernel coordinates \(U_0,U_0-U_1\)
gives extension cocycles \(u,z\) after
dividing by the action on the quotient.
The action is trivial exactly when
\(\chi=e=1\) and \(u=z=0\), proving
the assertion about the kernel field.

## Two factors and their local images

Let \(\mathscr B=\mathscr A/W^2\) using the
zero-sum subgroup of WR1, and factor
\(\Phi=\tau\rho\) with
\(\rho:\mathscr A\to\mathscr B\) and
\(\tau:\mathscr B\to\mathscr R\).
The product isogeny \(\prod_j\phi_j\)
factors through \(\rho\) by
\(\eta:\mathscr B\to\prod_j E_j'\).
Then
\[
                   (\textstyle\prod_j\widehat\phi_j)\eta=\Psi\tau .
\]
On \(\ker\tau=W^*\), \(\eta\) is the diagonal
map to \(\prod_j E_j'[\widehat\phi_j]\).

For any completion let \(I_j(v)\subset H^1(\mathbb Q_v,W^*)\)
and \(J_j(v)\subset H^1(\mathbb Q_v,W)\)
be the local images of \(\widehat\phi_j\) and \(\phi_j\).

**Proposition WR2.** For every \(a,b\) as in WR1,
the local image of \(\rho\) is
\[
 \{(c_0,c_1,c_2)\in\textstyle\prod_j J_j(v):
                           c_0+c_1+c_2=0\}.
\]
The local image of \(\tau\) is contained in
\(\bigcap_j I_j(v)\). In the tangent family
of [CS1](ABC_CUBIC_TANGENT_COMMON_SELMER_CLASS.md#three-local-images-in-a-common-quadratic-algebra)
with \(a=2,\ b=(3k)^2D\), this containment is
equality at every place. The same equality holds
for \(k=1\) with the special place \(k\) omitted.

The zero-sum subgroup \(W^2\subset W^3\)
is a direct summand of the Galois module.
The composition of its map to
\(H^1(\mathbb Q_v,W^3)\) and then
\(H^1(\mathbb Q_v,\mathscr A)\) tests the
three individual \(\phi_j\) local images.
Its kernel is the local Kummer image for
\(\rho\), giving the stated equality.
For \(\tau\), the isogeny diagram shows that
the Kummer images of the three coordinates
\(\Psi_j(R)\) of an \(\mathscr R\)-point \(R\)
are the same common class. This proves the
containment.

At a place where \(L\otimes\mathbb Q_v\)
has a rational factor, a point supported on
that factor has \(\Psi_0(R)\) equal to any
prescribed point of \(E_0(\mathbb Q_v)\).
This gives equality there. If the local cubic
field is unramified, the prime is different
from \(3\), and \(E_0\) has good reduction,
the norm from the cubic extension surjects
onto \(E_0(\mathbb Q_v)\). On the formal
kernel this follows by division by three
of rational points. On the finite fibre the
norm is \(1+F+F^2\), for Frobenius \(F\);
a preimage in the algebraic closure of a fixed
point of \(F\) is fixed by \(F^3\) because
\((F-1)(1+F+F^2)=F^3-1\).
This gives equality at the remaining good primes.

At \(2\), at the ramified quadratic places
\(q>3\) dividing \(R\), and at infinity,
the local intersection is zero by CS1.
At \(k\) there is a rational cubic factor if
\(2\) is a cube in \(\mathbb Q_k\).
In the other case the three lines have zero
intersection by CS1.
It remains to examine \(3\).

### A local trace at three

Use \(H=\mathbb Q_3(r)\) with \(r^3=2\),
and \(\pi=r+1\), which is Eisenstein at
\(3\). Let \(d_0=1+r+r^2=\pi^3/3\).
Its norm is one and
\(N(d_0)\operatorname{Tr}(d_0^{-1})=-3\)
since \(d_0^{-1}=r-1\).
Write \(b=27B\) and \(\sqrt b=3\lambda w\)
with \(w\in\mathbb Z_3^\times\) and
\(\lambda^2=-3\), as in CS1.

There is a point \(P_3\in E_0(H)\) with
\(x=\pi^2,\ y=\pi^3\sqrt{1+b/\pi^6}\), since
\(b/\pi^6=3B/d_0^2\).
Use the square-root branch congruent to one.
Then \(u=y/3=d_0+3B/(2d_0)\pmod{9\mathcal O_H}\).
Write \(e_i(u)\) for the elementary symmetric
functions in its three conjugates.
Both \(e_1(u)\) and \(e_2(u)\) vanish
modulo \(3\) since \(H\) is totally ramified.
Moreover \(e_2(u)\equiv e_2(d_0)=-3\bmod9\).
The first-order change under \(d_0\mapsto d_0+3z\)
is \(3(\operatorname{Tr}(d_0)\operatorname{Tr}(z)
-\operatorname{Tr}(d_0z))\).
Both terms are multiples of nine when
\(z=B/(2d_0)\).

The norm to \(\mathbb Q_3(\lambda)\) of
\(u+\lambda w\) has even part
\(N(u)+\lambda^2w^2 e_1(u)\), a unit,
and odd part
\[
 \lambda w e_2(u)+\lambda^3 w^3
   =3\lambda w\big(e_2(u)/3-w^2\big).
\]
The parenthesis is a unit, as its residue
is \(-1-1\bmod3\).
The quotient by the conjugate therefore
lies in \(U_3\) and is nontrivial modulo \(U_4\).
The discarded rational factor \(3\) in the
descent value \(y+\sqrt b\) has cubic norm.
By Kummer norm functoriality this class
is the image of the elliptic norm of \(P_3\),
which is the first coordinate \(\Psi_0(P_3)\).
As the common local intersection is exactly
the \(U_3\) line, this proves equality at \(3\).

## Global Selmer groups and the common obstruction

Write \(S_j=\operatorname{Sel}_{\widehat\phi_j}(E_j')\)
and \(T_j=\operatorname{Sel}_{\phi_j}(E_j)\).

**Corollary WR3.** In the tangent family,
\[
 \operatorname{Sel}_\tau(\mathscr B)=\bigcap_j S_j,\qquad
 \operatorname{Sel}_\rho(\mathscr A)
   =\{(c_0,c_1,c_2)\in\textstyle\prod_j T_j:
                         c_0+c_1+c_2=0\}.
\]
A common Selmer class maps from the associated
class in \(\Sha(\mathscr B)[\tau]\) to the
three Tate–Shafarevich classes of the original
coverings. In particular a class obstructed
on one of the twists gives a nonzero
class in \(\Sha(\mathscr B)[\tau]\).

The identities follow by enforcing the local
conditions of WR2. The morphism \(\eta\)
sends the kernel of \(\tau\) diagonally to
the three dual elliptic kernels. The image
of a \(\tau\)-Selmer class in \(H^1(\mathbb Q,\mathscr B)\)
therefore maps to the three corresponding
torsor classes. A nonzero image on a twist
proves the asserted nonzero class on \(\mathscr B\).
This applies to the three certified tangent
Hasse failures in
[CS5](ABC_CUBIC_TANGENT_COMMON_SELMER_CLASS.md#an-explicit-failure-of-the-hasse-principle)
and [TD9](ABC_CUBIC_TANGENT_PURE_CUBIC_TWO_SELMER.md#two-certified-prime-gap-hasse-failures).

**Corollary WR4.** Suppose the signed tangent
quartic is squarefree with \(|R|>3\).
Let \(\epsilon_k=1\) if \(k>3\) and either
\(k\equiv2\bmod3\) or \(k\equiv1\bmod3\)
with \(2\) a noncube modulo \(k\).
Set \(\epsilon_k=0\) otherwise, and put
\(\epsilon_D=\mathbf 1_{D<0}\).
The global Selmer ratio is
\[
 c(\tau)=3^{-1-\epsilon_k-\epsilon_D},\qquad
 \dim \operatorname{Sel}_{\tau^\vee}(\mathscr R^\vee)
      =\dim\bigcap_j S_j+1+\epsilon_k+\epsilon_D .
\]
In the imaginary prime subfamily with
\(k\equiv2\bmod3\), the dual group has
dimension at least four when the tangent
class is nonzero.

The local ratios of \(\tau\) are its local
Kummer image sizes divided by the local
kernel sizes. At \(2\) the size ratio is
\(1/3\). At \(3\) the common image and kernel
both have order three. At \(k\) the ratio
is \(1/3\) exactly when \(\epsilon_k=1\).
At infinity it is \(1/3\) exactly for \(D<0\).
The other factors are one, as follows from
good reduction or the ramified quadratic
calculation of CS1. The two global rational
kernels are trivial because \(D,-3D\)
are nonsquares under the stated conditions.
Now apply the
[Bhargava–Klagsbrun–Lemke Oliver–Shnidman, Proposition 2.1](https://lemkeoliver.github.io/papers/22-ShaElements.pdf).
The lower bound uses the common nonzero
class of CS2.

## The secondary lifting condition

**Proposition WR5.** Every common Selmer class
\(\alpha\) of WR3 has a global cohomology lift
\(z\) in \(H^1(\mathbb Q,\ker\Phi)\).
For each place \(v\), take a local \(\Phi\)-Kummer
lift \(z_v\) of \(\alpha_v\), and take
\(x_v\in H^1(\mathbb Q_v,W^2)\) mapping to
\(z_v-z|_v\). Write \(L_\rho(v)\) for the
local zero-sum group of WR2. Then \(\alpha\)
comes from \(\operatorname{Sel}_\Phi(\mathscr A)\)
if and only if
\[
  (x_v\bmod L_\rho(v))_v\ \in\
    \operatorname{im}\!\left[H^1(\mathbb Q,W^2)\longrightarrow
         \bigoplus_v H^1(\mathbb Q_v,W^2)/L_\rho(v)\right].
\]
Equivalently it lifts if and only if the sum of the local
Tate pairings \(\sum_v\langle x_v,\gamma_v\rangle_v\)
vanishes for every
\(\gamma\in\operatorname{Sel}_{\rho^\vee}(\mathscr B^\vee)\).
For the dual group, take tuples of local classes
in \(H^1(\mathbb Q_v,W^*)^3\) modulo the diagonal.
The dual local condition consists exactly of
the images of \(\prod_j I_j(v)\) in that quotient.

The extension of WR1 gives a connecting class
in \(H^2(\mathbb Q,W^2)\). It is everywhere
locally zero since each local \(\tau\)-Kummer
lift along \(\tau\rho\) is the image of a
local \(\Phi\)-Kummer class. Restrict to the
biquadratic splitting field of \(W\) and
\(\mu_3\). Restriction is injective on
three-torsion by corestriction, and the local
restriction map for \(H^2\) upstairs is
injective by the Brauer group local-global
sequence. Thus \(z\) exists.

The ambiguity in \(x_v\) from
\(H^0(\mathbb Q_v,W^*)\) belongs to the
Kummer subgroup for \(\rho\), since it
is the boundary of a point of \(\ker\tau\)
on \(\mathscr B\). A change of the local lift
\(z_v\) changes its coset by zero. Indeed a
local \(\Phi\)-Kummer class mapping to zero
under the \(\tau\) Kummer map comes from
a point in \(\tau(\mathscr B(\mathbb Q_v))\),
and its lift is then a \(\rho\)-Kummer class.
A change of the global lift \(z\) changes the
coset vector by a global localization, which
proves the first criterion. At all but finitely
many good primes the isogenies have their
unramified Kummer conditions, so the vector
has finite support.

For the pairing criterion use the
[Poitou–Tate comparison of Selmer structures](https://msp.org/ant/2018/12-10/ant-v12-n10-p02-p.pdf),
Theorem 3.7, by allowing full local classes
at the finite set supporting the discrepancy.
Its rightmost localization images are
orthogonal complements. The zero-sum
condition in WR2 lives in the direct
summand \(W^2\subset W^3\).
The Tate pairing therefore identifies its
dual with the quotient by the diagonal
\(W^*\). The annihilator of the intersection
with \(\prod J_j(v)\) is the image of the
annihilator \(\prod I_j(v)\) in this quotient.
This gives the asserted local test on
\(\gamma\).

A nonzero secondary pairing obstructs
rational realization by the simultaneous
\(\tau\)-covering. Its vanishing permits a
further \(\Phi\)-Selmer descent, whose
kernel is the image of the zero-sum
\(\rho\)-Selmer group.
Rational realizability by \(\tau\) implies
rational realizability on every original
twist by \(\Psi\). The prime-gap Hasse
question asks for failure on the third
twist, so a failure of simultaneous
realizability alone leaves that question open.

## The three elliptic pairings on the common group

**Theorem WR6.** For \(\alpha,\beta\in\bigcap_j S_j\),
the elliptic Cassels–Tate pairings, with the
common kernel orientations of WR1, satisfy
\[
             \sum_{j=0}^2 \langle\alpha,\beta\rangle_j=0 .
\]
In particular the displayed tangent class, which is
rational on the first two twists, pairs trivially
on the third twist with every common Selmer class.

If \(k=1\), \(k\equiv2\bmod3\), or \(k\) splits
and \(2\) is a cube modulo \(k\), choose
\(j_0\in\{0,1\}\) with \(j_0\ne e_*\) from CS1.
A common class lifts to the degree-\(27\)
Selmer group if and only if its two pairings
against the full \(S_j\) for \(j\ne j_0\) vanish.
In this range the displayed tangent class
therefore lifts precisely when it has a
full three-Selmer lift on \(E_2\).
If \(k\) splits and \(2\) is a noncube, the
full dual testing group of WR5 modulo the
image of \(\prod_j S_j\) has dimension
at most one.

For the pairing formula, project the
global and local lifts \(z,z_v\) of WR5
to \(H^1(E_j[3])\). The local projection
is the full three-Kummer class of
\(\Psi_j(R_v)\), since
\(\Psi\Phi=[3]\).
Under \(W^2\subset W^3\), the coordinates
of \(x_v\) sum to zero, and they are the
three individual discrepancies between
the global and local lifts.
The Weil pairing definition of
[Van Beek–Fisher, Definition 1.2](https://arxiv.org/pdf/1711.02432)
takes the sum of their local pairings
against the individual tester \(\beta\).
The kernel pairing is the same in each
coordinate by transport through \(\psi_j\)
and the orientation \(\phi_j(Q)\) of WR1.
Thus the sum across twists is zero place
by place. Rational classes pair to zero
because they have global Kummer lifts.

For the second assertion, \(I_{j_0}(v)\)
is contained in both other local images
at every place by CS1 and the condition
on \(k\). Fix the representative with
coordinate \(j_0\) zero in the quotient
by the diagonal. The image of
\(\prod I_j(v)\) then becomes simply
the product of the other two local
groups. Thus the full dual testing
group of WR5 is \(S_j\oplus S_{j'}\)
in those coordinates. The two elliptic
pairings give exactly the two required
vanishings, by the same local discrepancy
calculation. When \(2\) is a noncube at the
split \(k\), the same argument works outside
\(k\). At \(k\) the image of the three
distinct one-dimensional lines in the
diagonal quotient is three-dimensional.
It contains the direct sum of the two
non-anchor lines, of dimension two.
Localization modulo that two-dimensional
subspace has kernel \(S_j\oplus S_{j'}\)
and image of dimension at most one.
The image of \(\prod_j S_j\) contains the
kernel, giving the bound.

## Two explicit failures that pass the secondary descent

**Corollary WR7.** For the two imaginary
Hasse examples with \((p,k)=(17,1)\)
and \((101,5)\), the common class
lifts to \(\operatorname{Sel}_\Phi(\mathscr A)\).
Its nonzero threefold Tate–Shafarevich
class is in the image of
\(\Sha(\mathscr A)[\Phi]\) under \(\rho\).
On the third elliptic twist, the
obstruction class on \(E_2'\) is the
\(\phi_2\)-image of a class in \(\Sha(E_2)[3]\).
The three-torsion sizes of \(\Sha(E_2)\)
in these two examples are \(3^2\) and \(3^4\),
respectively.

For \((17,1)\) and \((101,5)\) the discriminants
of \(K_D\) are \(-186843\) and \(-256285803\).
Both are fundamental.
Their ordinary class numbers are \(96\) and \(2208\),
respectively, by enumerating positive reduced
binary quadratic forms as in the code below.
Each class number has \(3\)-valuation one.
The \(S_2\) local conditions in CS1 force virtual
units at all primes in both examples. At \(k=5\)
the image is zero, and at each split good prime
it is the unit line. The remaining quadratic
field completions have no nonzero valuation
class in the norm kernel. Thus \(S_2\) is
contained in the one-dimensional space of
virtual units. It contains the common
class, so \(S_2=\bigcap_j S_j\).

The common class pairs with all of \(S_2\)
trivially by WR6. It therefore lifts to full
three-descent on \(E_2\), and then to the
simultaneous group by the criterion of WR6.
Its projection under \(\phi_2\) gives the
asserted Tate–Shafarevich class on \(E_2'\).
The quotient diagram gives the simultaneous
assertion on \(\mathscr B\).
By CS3, including the omission of the \(k\)-factor
at \(k=1\), the dimensions of \(T_2\) are
\(1\) and \(3\), respectively. The rational
three-torsion kernels are zero.
The compositional descent sequence for
\([3]=\widehat\phi_2\phi_2\) thus gives
full three-Selmer dimensions \(2\) and \(4\).
The corresponding \(E_2(\mathbb Q)\) are
trivial by CS5 and TD9. Hence these full
three-Selmer sizes are exactly the claimed
three-torsion sizes.

The reduced-form enumeration uses
\(|b|\le a\le c\) and \(b\ge0\) on the tie
boundaries. Both discriminants are odd,
so \(b=0\) does not occur.

```python
from math import gcd, isqrt

def class_number(discriminant):
    maximum_a = isqrt((-discriminant) // 3)
    count = 0
    for b in range(1, maximum_a + 1, 2):
        product = (b*b - discriminant) // 4
        for a in range(b, min(maximum_a, isqrt(product)) + 1):
            c, remainder = divmod(product, a)
            if remainder == 0 and gcd(a, gcd(b, c)) == 1:
                count += 1 if a == b or a == c else 2
    return count

assert class_number(-186843) == 96
assert class_number(-256285803) == 2208
```

## A real example with one additional simultaneous tester

**Corollary WR8.** For the certified real
Hasse example \((p,k)=(167,19)\), all three
\(S_j\) equal the line of the common class.
That class lifts to full three-descent
on every individual twist, and
\(\left|\Sha(E_2)[3]\right|=3^2\).
The simultaneous dual testing group in WR5
has a one-dimensional quotient by the image
of \(\prod_j S_j\). An explicit tester in
this quotient is given below.

Here \(D=10179717\). Write
\(\omega=(1+\sqrt D)/2\).
Take \(\mathfrak q\) above \(19\) with
\(\omega\equiv4\bmod\mathfrak q\),
a generator \(g\) of norm \(19\), and
a fundamental unit \(\epsilon\).
Set \(\gamma=g/\bar g\).
The ordinary class group is cyclic of
order \(3\), verified by the class-group
certificate below. Thus the virtual unit
space has dimension two. Its basis consists
of the common class and \([\epsilon]\),
since the common class vanishes at \(2\)
and \(\epsilon\) has nontrivial residue there.
The space allowing valuations at \(19\)
has one extra direction \([\gamma]\), of
opposite valuations \(1,-1\). Indeed
subtracting its suitable multiple from
a norm class removes the valuations modulo
three, leaving a virtual unit.

At \(3\), \(\gamma\equiv\omega+2\bmod3\)
is non-scalar, and \(\epsilon\) is scalar.
The line \(U_3\) consists of classes whose
projective residue modulo \(3\) is scalar.
In fact this residue detects the \(U_1\)
coordinate in the anti-invariant space.
At \(2\) the class of \(\gamma\) is trivial.
At \(19\), encode a cubeclass by its valuation
modulo three and by the sixth power of its
normalized unit in \(\mathbb F_{19}^\times\).
The three local generators \(2^{j+1}s_0\)
of CS1 have coordinate pairs
\((1,11),(1,1),(1,7)\), respectively.
The coordinates of \(\gamma,\epsilon\)
are \((1,7),(0,7)\).
The common class is zero there.

Here \(e_*=1\). For the \(S_0\) condition,
the \(U_3\) line excludes the \(\gamma\)
direction and the place \(2\) excludes the
unit direction. For \(S_2\), the \(U_3\)
line again excludes \(\gamma\) and the
unit part at \(19\) excludes \(\epsilon\).
For \(S_1\), place \(2\) excludes \(\epsilon\)
and the \(19\)-line excludes \(\gamma\).
Hence all three groups equal the common
line. Alternation of the individual
Cassels–Tate pairing implies the full
three-descent lift. The ratio from CS3
gives \(\dim T_2=1\), since the quadratic
field is real, \(k\) splits, and \(e_*\ne2\).
The torsion and rank vanish by TD9,
so the descent dimensions give
\(\left|\Sha(E_2)[3]\right|=3^2\).

The triple \((0,[\gamma],[\epsilon^2])\)
modulo the diagonal is a tester for WR5.
At \(2\) the only nonzero component lies
in the allowed last line. At \(3\) the
second component is unrestricted and the
last lies in \(U_3\). At \(19\) take
one copy of the \(I_0\) generator,
two of the \(I_1\) generator and one
of the \(I_2\) generator and subtract
the first from all three. Their two
differences have the specified coordinates
\((1,7)\) and \((0,7^2)\).
At the other primes both entries are
virtual units and lie in the local images.
This tester is outside the image of
\(\prod S_j\), whose evaluations at
\(19\) all vanish. The upper bound
in WR6 proves the asserted dimension.

The computation of the local data uses
the following exact PARI/GP replay.
The relation \(N g=19\) gives
\(\gamma/19=1/\bar g^2\) at \(\mathfrak q\).

```gp
default(realprecision,600);
d=10179717;
B=bnfinit(x^2-x+(1-d)/4,1);
if(B.cyc!=[3]||bnfcertify(B)!=1,error("quadratic class certificate"));
P=idealprimedec(B,19)[1];
g=nfbasistoalg(B,bnfisprincipal(B,P,3)[2]);
if(norm(g)!=19||idealhnf(B,g)!=idealhnf(B,P),error("prime generator"));
z=g/conj(g);
e=quadunit(d);
u=Mod(component(e,2)+component(e,3)*x,B.pol);
print([lift(Mod(lift(z),2)),lift(Mod(lift(u),2))]);
print([lift(Mod(lift(z),3)),lift(Mod(lift(u),3))]);
v=nfmodprinit(B,P);
if(lift(nfmodpr(B,Mod(x,B.pol),v))!=4,error("prime orientation"));
y=nfmodpr(B,2*Mod(x,B.pol)-1,v);
w=nfmodpr(B,conj(g),v);
print(vector(3,j,lift((2^j*3*y)^6)));
print([lift(w^(-12)),lift(nfmodpr(B,u,v)^6)]);
```

The four output rows are \([1,x]\),
\([x+2,2]\), \([11,1,7]\), and \([7,7]\).
The individual descents leave the pairing
against this tester undetermined.
The following local calculation resolves it.

## The mixed pairing in the real example

**Theorem WR9.** For \((p,k)=(167,19)\)
the common class also lifts to
\(\operatorname{Sel}_\Phi(\mathscr A)\).
Its nonzero simultaneous obstruction in
\(\Sha(\mathscr B)[\tau]\) is in the image
of \(\Sha(\mathscr A)[\Phi]\).
With the local Kummer orientations below,
the pairing against the extra tester of WR8
has contributions \(2/3,0,1/3\) at \(2,3,19\),
respectively, and zero at other places.

For a point \((x,y)\) on \(y^2=x^3+b_j\)
write \(s_j^2=b_j\) and \(w^3=y+s_j\).
A dual-isogeny preimage has coordinates
\[
 U=3(x+w^2+x^2/w^2),\qquad
 V=\frac{27yU^3}{U^3+216b_j}.
\]
This is Cardano's formula for
\(U^3-9xU^2-108b_j=0\) and the
ordinate formula of the dual isogeny.
Changing the cube root translates the
preimage by dual-kernel torsion.

Use \(P_0=(p(p-h),u_-/2)\), and for \(P_1\)
use the negative of the third-colour point
\((p^2-h^2,p^3-3ph^2)\).
Put \(\alpha=y(P_0)+s_0\) and
\(\beta=(\sqrt D-p^2-h^2)/2\).
The identity \(\beta^3=\alpha_-\overline{\alpha_+}\)
of CT1 shows that the Kummer value at \(P_1\)
is \(\alpha\Delta^3\), where
\(\Delta=-(p-h)\beta/\alpha\).
Thus the respective Cardano roots \(w\)
and \(w\Delta\) give compatible preimages.
In WR1 the projection of \(\ker\Phi\) to
the first two coordinates is precisely
the fibre product of \(E_0[3]\) and \(E_1[3]\)
over their common \(W^*\). The two
rational points with these preimages
therefore specify a global lift as in WR5.

At \(2\) and \(19\) the common local class
is zero, so take zero as local lift.
Write \(u_0,u_1\) for \(V+3\lambda s_j\)
on the two Cardano preimages, where
\(\lambda^2=-3\). The local pairing against
the tester \((0,[\gamma],[\epsilon^2])\)
is then the sum of the pairings of
\(-[u_1]\) with \([\gamma]\) and
\([u_0u_1]\) with \([\epsilon^2]\),
using the zero-sum formula of WR5.
The tame Hilbert symbol and restriction
factor are as in
[Van Beek–Fisher, Lemma 1.3 and Proposition 1.5](https://arxiv.org/pdf/1711.02432).
At \(2\), the valuations of \(u_0,u_1\)
are \(0,4\) and the residue of \(\epsilon^2\)
is \(\bar\omega^2\) over the unramified
quadratic field. The restriction factor
is \(1/2\) in \(\tfrac13\mathbb Z/\mathbb Z\).
At \(19\) the valuations are \(2,0\);
the sixth-power normalized residues are
\(7,1\). WR8 supplies the corresponding
tester residues. These give the first
and third entries of the asserted sum.

At \(3\) the unit \(\epsilon\) is
\(8\bmod9\) in the integral quadratic
basis, as obtained by the same exact
unit calculation in WR8 with modulus
nine. Its local cubeclass is trivial
since \(U_1^3=U_4\) in \(\mathbb Q_3(\lambda)\).
Only the second coordinate of the
discrepancy can contribute. For the
local point of \(\mathscr R\), use the
negative of the point \(P_3\) in WR2.
In the notation of that proof,
\(y(P_0)\equiv2\bmod3\).
Its Kummer quotient by its conjugate
has leading term \(3\lambda w\)
modulo \(\lambda^4\), whereas the
elliptic norm of \(P_3\) has leading
term \(6\lambda w\).
Negating \(P_3\) gives the same
nonzero local class as \(P_0\).

Use the unramified cubic extension
\(U/\mathbb Q_3\) to take dual-isogeny
preimages \(Q_j\) of the three
\(\Psi_j(-P_3)\).
The true tuple of preimages from
\(\eta\) exists over \(U\) because
the \(U_3\)-line Kummer class becomes
trivial under the unramified cubic
extension. The checker selects the tuple by
two tests. First the sum of the
three aligned preimages equals
\(\phi_0(-P_3)\) in the first
restriction coordinate. This tests all
three rows since differences between
candidate preimages lie in the three
dual kernels and every row acts on
such points by the same sum.
Second the sum of the three individual
\(\phi_j\)-Kummer classes of \(Q_j\)
vanishes over \(U\). It is the
Kummer obstruction for \(\eta\), whose
kernel is the sum quotient of \(W^3\)
by \(W^2\).
Together the tests suffice. After lifting
to \(\mathscr B(U)\), the image under
\(\tau\) differs from \(-P_3\) by a
point of \(\operatorname{Res}_{H/\mathbb Q_3}W\).
There is no such nonzero point over
\(U_3\), since \(HU\) has odd
degree over \(\mathbb Q_3\) and \(W\)
has nontrivial ramified quadratic character.

The difference of the global and
local preimages in the second coordinate
has its \(\phi_1\)-Kummer value in
\(\mathbb Q_3^\times\). The replay gives
valuation \(3\) and normalized unit
\(-1\bmod9\) for one of the selected
tuples, and verifies the cube condition
for every diagonal shift.
Thus this discrepancy pairs trivially.
The checker isolates the Cardano roots
by Hensel estimates, and tests equality
of the first row using residues
that distinguish the two nonzero
torsion translations of the target.
All unit and valuation tests are
on finite residues. See
[the Sage replay](verify_cubic_tangent_mixed_selmer.py).

At every good prime outside \(2,3,19\)
the lifts and testers can be taken
unramified, so the local pairing is
zero. At the other primes dividing
\(R\) the local tester group is zero,
and the odd-order pairing at infinity
is zero. Thus the mixed pairing
vanishes globally. All individual
tests also vanish by WR8.
WR5 gives a full simultaneous lift,
while TD9 still gives a nonzero
obstruction on the third twist.
