# Pure cubic ray classes and the third tangent twist

The third cubic twist of the tangent progression has
a two-descent in a pure cubic field. For one of the
two possible local types at \(3\), its entire two-Selmer
group is a ray class character group with modulus
supported at \(2\). In the other type, a single
explicit paired condition at \(3\) is added.

## The exact local conditions

Let \(m\) be a signed cubefree integer with
\(m\equiv2\bmod8\) and \(3\nmid m\). Put
\[
 T=\mathbb Q(\rho),\qquad \rho^3=m,\qquad
 \theta=6\rho,\qquad E(m):y^2=x^3-216m .
\]
The cubic is irreducible by its valuation at \(2\),
so \(E(m)\) has no rational two-torsion.
Let \(\mathfrak p_2\) be the unique
prime above \(2\) in \(T\).

**Theorem TD1.** Identify the two-Selmer
group of \(E(m)\) in
\[
           \ker\bigl(T^\times/T^{\times2}
                          \xrightarrow{N}\mathbb Q^\times/\mathbb Q^{\times2}\bigr)
\]
using the class of \(x-\theta\).
Its local conditions have the following
complete description.

- At all finite primes other than \(3\),
  the valuations in \(T\) are even.
  At \(\mathfrak p_2\), after removing a
  uniformizer square to get a unit, the
  unit must reduce to \(1\bmod\mathfrak p_2^2\).
- If \(m\not\equiv\pm1\bmod9\), \(3\) is totally
  ramified and imposes no condition beyond the
  global norm. If \(m\equiv\pm1\bmod9\), put
  \(u^3=m\) in \(\mathbb Q_3\) and
  \(\theta_0=6u\).
  With \(T\otimes\mathbb Q_3=
  \mathbb Q_3\times\mathbb Q_3(\zeta_3)\),
  the required local class is either trivial
  or the class of
  \((3\theta_0^2,\theta_0(1-\zeta_3))\).

Use the odd-degree
[Schaefer descent](https://arxiv.org/pdf/1507.08325)
with \(p=2\) and a separable cubic.
The norm kernel of the three
geometric squareclass coordinates
is the Galois module of two-torsion,
and the evaluation is \(x-\theta\).
At a good odd prime the local
image is the norm-square unit
space. Indeed the evaluation has
even valuations, and the image
has dimension \(\dim E(m)(\mathbb F_q)[2]\).
If the local cubic algebra has
\(r\) factors, that dimension is
\(r-1\), which is also the dimension
of the norm-square residue-unit
space. At any odd prime dividing
\(m\), the cubic completion is
totally ramified by the root
valuation, and norm on its local
squareclasses is an isomorphism.
At infinity the real component
of a norm-square value is positive
and gives the trivial class.

At \(2\) the change \(x=4X,\ y=8Y+4\)
gives the good model
\[
 Y^2+Y=X^3-(27m+2)/8 .
\]
Its group over \(\mathbb F_2\) has order three.
As the local two-torsion is zero, the local
Kummer image has dimension one by
[the local Kummer index formula](https://arxiv.org/pdf/2001.02263)
(Lemma 1.3). A point at \(x=1\) exists by
the unit square criterion in \(\mathbb Q_2\)
and gives \(\gamma=1-6\rho\).
The difference
\(\gamma-(1+\rho^2)^2=-(6+m)\rho-2\rho^2\)
has valuation five, where \(v_\rho(2)=3\).
If the quotient \(\gamma/(1+\rho^2)^2\) were
\((1+c)^2\), the valuation of \(2c+c^2\)
would be even for \(v_\rho(c)<3\) and at
least six for \(v_\rho(c)\ge3\).
Thus \(\gamma\) is nonsquare and spans the image.
The norm-square unit squareclasses
of \(T_2\) have dimension two.
Reduction modulo \(\rho^2\) on
unit squareclasses is a character
of order two, nonzero even on
the norm-square subspace:
\((1+\rho)/N(1+\rho)\) has square
norm and nontrivial residue.
The class of \(1-6\rho\) has
trivial residue there. Its line
is therefore precisely the kernel.

At \(3\) the unit cube criterion in
\(\mathbb Q_3\) gives the product cubic
algebra exactly at residues \(1,-1\bmod9\). Otherwise the cubic
completion is a field and totally
ramified: even in an unramified
cubic extension a cube root of
a rational unit would have
residue in \(\mathbb F_3^\times\),
and its cube would be \(1\) or
\(-1\bmod9\).
In the split case there is one
nonzero rational two-torsion
point locally, \((\theta_0,0)\).
Its Kummer class under \(x-\theta\)
is the stated pair, using the
derivative \(3\theta_0^2\) at
the vanishing component.
It is nonzero by its odd
valuation in the quadratic
factor. Since the place is
odd, the image dimension
is one. This proves both the
necessity and sufficiency of
the listed conditions.

## A ray class description and an arithmetic rank-zero test

**Corollary TD2.** If \(m\not\equiv\pm1\bmod9\),
the two-Selmer group is canonically the
group of quadratic characters of the
ray class group of \(T\) of modulus
\(\mathfrak p_2^2\) with no real place
in the modulus. If \(m\equiv\pm1\bmod9\),
use the modulus obtained by also
including both primes above \(3\)
with exponent one. Impose the
paired local condition at those
two primes from TD1 on its
quadratic characters.

Use the ray class correspondence and
the conductor-discriminant formula, as stated
in [Milne's *Class Field Theory*](https://www.jmilne.org/math/CourseNotes/cft.html),
Chapter V, Theorems 3.6 and 3.25.
A quadratic character is identified with its
Kummer squareclass in \(T\). Away from \(2,3\), the
even valuations of TD1 give the
unramified condition. At the real
place the global norm square
gives positive sign and a split
local extension. In the nonsplit
case at \(3\) the character is
trivial there by total ramification
and norm. In the split case the
nonzero pair of TD1 corresponds
to tame ramification at both
primes, of exponent one.

At \(2\) the nonzero Kummer class
\(\gamma=1-6\rho\) has conductor
exponent two. To see this,
divide by the square
\((1+\rho^2)^2\) and call the
result \(d_0\).
Since \(m\equiv2\bmod8\),
\[
 \gamma-(1+\rho^2)^2=-(6+m)\rho-2\rho^2
\]
has valuation five.
The element \((\sqrt{d_0}-1)/\rho^2\)
has an Eisenstein quadratic
equation with discriminant
\(4d_0/\rho^4\) of valuation two.
The other generator of the
local characters of conductor
at most two is the unramified
class, represented by \(5\).
The norm of this class is
not a square in \(\mathbb Q_2\).
The group has order four by
local reciprocity, since it
has one unramified direction and
\((\mathcal O_{T_2}/\mathfrak p_2^2)^\times\)
has order two. Hence its norm-square
part consists exactly of the
two classes of TD1.
Conversely the stated ray
characters have norm a square:
outside the split case at \(3\)
their valuations are even,
in that case the two odd
valuations occur together with
residue degree one, and the
real sign is positive.
This proves the ray interpretation.

## The marked tangent family

For the squarefree signed tangent quartic
\(R=3p^4-6p^2h^2-h^4\), take \(h=6k\)
as in
[CS1](ABC_CUBIC_TANGENT_COMMON_SELMER_CLASS.md#three-local-images-in-a-common-quadratic-algebra).
The third twist \(E_2\) is \(E(m)\) with
\(m=2k^2(R/3)\). Allow \(k=1\) as well.
Here \(m\equiv2\bmod8\), \(3\nmid m\) and
\(R/3\equiv p^4\bmod9\), so the residues at
\(3\) are \(2,5,8\bmod9\).

**Corollary TD3.** Suppose this marking gives a
nonzero common tangent class, \(m\not\equiv-1\bmod9\),
and \(|R|>3\). If the ray class group of \(T\)
of modulus \(\mathfrak p_2^2\) has odd order, the
common-class covering fails the Hasse principle.

The rank of \(E_2(\mathbb Q)\) is zero by TD2.
Odd torsion injects in the good fibre at \(2\),
which has order three, and there is no
two-torsion. Choose a prime \(q>3\) dividing \(R\).
It has valuation one on \(m\), since \(q\ne k\).
The nonzero three-torsion abscissas in the
three-division equation satisfy \(x=0\) or
\(x^3=864m\). The valuations at \(q\) exclude
the required ordinate square at \(x=0\) and
the abscissa cube in the second case.
Thus the entire rational group is zero.
By CS2 and CS4 the nonzero shared class
gives an everywhere locally soluble curve
without a rational point.

The test holds if the ordinary class number
of \(T\) is odd and a fundamental unit has
nontrivial residue modulo \(\mathfrak p_2^2\).
This follows from the ray class exact sequence,
since the residue unit group has order two.
The example \(m=124562\) of CS5 has
class group \([3,3]\) and such a unit.

The condition \(p\equiv1,\ k\equiv5\bmod9\)
is compatible with the prime
and squarefree sieve and gives
\(m\equiv5\bmod9\).
It also lies in the range
of the Selmer lower count
of CS3, since \(k\equiv2\bmod3\)
and \(e_*=0\). Thus odd ray
class number for infinitely
many of these marked \(m\)
would produce the unbounded
Hasse failures in the same
prime family.

## The one-dimensional ramified choice at three

**Theorem TD4.** For the split type \(m\equiv\pm1\bmod9\)
in TD1 let \(V_0\) be the quadratic character group of
the finite ray modulus \(\mathfrak p_2^2\), and \(B\)
the one with both primes above \(3\) added. Then \(V_0\)
is the kernel of localization at \(3\) on \(B\), and
\([B:V_0]=2\). Its unique nonzero local image is ramified
at both primes above \(3\). The group
\(\operatorname{Sel}_2(E(m))\) is either \(V_0\) or \(B\).
The local image of \(B\) has the unique nonzero class
\(\big[(3\theta_0^2,\theta_0(1-\zeta_3))\big]\) or
\(\big[(3\theta_0^2,-\theta_0(1-\zeta_3))\big]\).

Every character in \(B\) has positive norm sign and
even norm valuations away from \(3\). Its rational norm
class is therefore \(1\) or \(3\). Its dyadic component
has conductor exponent at most two, whose local
norm classes are \(1\) and \(5\) by TD2. This excludes
\(3\) and proves that all of \(B\) is in the norm
kernel. The local description of TD1 now identifies
\(B\) with the relaxed two-Selmer group at \(3\):
every element of the local norm kernel at \(3\)
is allowed by the tame ray conductor bound.

For a norm-square class \(t\) define the ternary
form \(q_t(z)=\operatorname{Tr}_{T/\mathbb Q}(t z^2)/3\).
This is the cubic trace form of
[Shnidman–Siad, Proposition 4.2 and Corollary 5.1](https://arxiv.org/pdf/2506.13749).
At every place other than \(3\), a class in \(B\)
is a local Kummer value. Such a value can be represented
by \(t=x-6\rho\), since a local point can be changed
within its coset modulo doubling to avoid the exceptional
evaluations. The vector \(\rho^2\) is isotropic
because \(\operatorname{Tr}(\rho)=\operatorname{Tr}(\rho^2)=0\).
The same trace form is therefore locally isotropic
at \(3\) by reciprocity for the quaternion
even Clifford algebra of a ternary form.
The nontrivial unramified norm-square class at
\(3\) is \([(1,-1)]\) in the product algebra.
Its trace form, multiplied by \(3\), is
\(a^2-2b^2+6c^2\), using \(\sqrt{-3}\) in the
quadratic factor. A primitive isotropic triple
over \(\mathbb Z_3\) would first force \(3\mid a,b\)
modulo three and then \(3\mid c\) modulo nine.
This excludes the unramified line.

The strict group, defined with zero local image
at \(3\), is consequently \(V_0\).
The dimension difference of the relaxed and
strict groups is the dimension of the local Kummer
line, which is one by TD1 and
[Mazur–Rubin, Lemma 3.2](https://arxiv.org/pdf/0904.3709).
The norm-square space of the split cubic
algebra has dimension two, and its other
two lines are exactly the pairs in the
statement. The result follows by intersection
with the line of TD1.

## Transfer to a fixed cubic-twist family

**Theorem TD5.** For \(m\) in TD1 and \(s\in\{1,-1\}\),
define \(C_s(m):y^2=x^3-216s m^2\).
Its two-division field is the same pure cubic field
\(T\), and hence its two-Selmer group and that of
\(E(m)\) are in the same norm-square group.
If \(m\not\equiv\pm1\bmod9\), all three Selmer
groups coincide. If \(m\equiv\pm1\bmod9\),
the Selmer group of \(E(m)\) equals that of \(C_s(m)\)
for \(s\equiv m\bmod3\). In this latter case the
two groups of \(C_+(m)\) and \(C_-(m)\)
in some order are \(V_0\) and \(B\).

The Kummer evaluation for \(C_s(m)\) is
\(x-6s\rho^2\). At all odd places except the split
case at \(3\) the proof of TD1 gives the same
conditions. The primes dividing \(m\) have
zero local norm kernel, the good odd primes
impose even valuations, and the totally
ramified type at \(3\) imposes no further condition.
At infinity the norm-square kernel is again zero.

At \(2\) the Kummer image is one-dimensional
by the local index formula and absence of local
two-torsion. The point at \(x=1\) exists and
gives \(\gamma_s=1-6s\rho^2\). In the basis
\(1,\rho,\rho^2\), residues modulo \(\rho^7\)
are stored modulo \(8,4,4\), respectively.
Squaring \(a+b\rho+c\rho^2\) with \(a\) odd gives
\((1+4bc,2(b+c),b^2+2c)\) in these residues.
The four triples from the two coefficient parities
are \((1,0,0),(1,2,1),(1,2,2),(5,0,3)\).
The triple \((1,0,2)\) of \(\gamma_s\) is absent.
Since \(\gamma_s\) is a norm-square unit
congruent to one modulo \(\rho^2\), its line
is the dyadic line of TD1.

In the split case at \(3\) the nonzero
two-torsion evaluation is
\((3(6su^2)^2,6su^2(1-\zeta_3))\) with \(u^3=m\).
The quotient of this pair by the one of TD1
has squareclass \((1,su)\). The class of \(su\)
is trivial in \(\mathbb Q_3(\zeta_3)\)
exactly when \(s\equiv m\bmod3\), by its
residue in \(\mathbb F_3^\times\).
The two choices are therefore the two
ramified lines of TD4. In particular
one Selmer group equals \(B\) and the other \(V_0\).

## The quartic-field count and injective marking

**Corollary TD6.** Over \(T\) from TD1, the nonzero
classes of \(V_0=\operatorname{Hom}(
 \mathrm{Cl}_{\mathfrak p_2^2}(T),\mu_2)\)
correspond to \(S_4\) quartic fields with cubic
resolvent \(T\) and discriminants
\(\operatorname{Disc}(T)\) or \(4\operatorname{Disc}(T)\).
In the split type, the \(|V_0|\) remaining
nonzero classes of \(B\) correspond precisely to
the quartic fields over that resolvent with
discriminants \(9\operatorname{Disc}(T)\) or
\(36\operatorname{Disc}(T)\). The marked pure
cubic fields from distinct \((p,k)\) with \(p>6k\)
and squarefree signed \(R\) are distinct.

Use [Cohen–Thorne, Theorem 2.2](https://arxiv.org/pdf/1302.5728)
on the quadratic extensions \(T(\sqrt t)/T\)
of square norm. The quartic discriminant is
\(\operatorname{Disc}(T)N\mathfrak d_{T(\sqrt t)/T}\).
The dyadic conductor contributes \(1\) or \(4\).
In the split type both primes above \(3\) have
residue degree one, and the ramified line contributes
\(9\) or neither prime contributes. Conversely these
discriminant ratios bound the quadratic ramification
by the corresponding finite ray modulus, because
the norm of its relative discriminant is the ratio.
Square norm forces positivity at the unique real
place and paired ramification at \(3\).
The number at the two additional discriminants
is \(|B|-|V_0|=|V_0|\).

For injectivity, equality of two pure cubic fields
implies equality of their Kummer lines over
\(\mathbb Q(\zeta_3)\) by adjoining their normal closure.
The map from rational cubeclasses is injective
into that quadratic extension by the norm.
Since \(v_2(m)=1\) for both markings, the valuation at \(2\) selects the same
generator of the two Kummer lines. Cubefreeness
determines \(|m|\), and its unique prime marked by
exponent two is \(k\), with no such prime if \(k=1\).
Thus \(|R|\) is the same. The sign is determined by
\(R\equiv3\bmod8\), and \(R\) is strictly increasing
in \(p>h\) for fixed \(h=6k\).

## Root signs and the range of the average

**Corollary TD7.** In the marked tangent family of
TD3 with \(k\) prime, the root number of \(E(m)\) is
\[
 W(E(m))=-\operatorname{sgn}(R)
              \left(\frac{k}{3}\right)
              (-1)^{\mathbf 1_{m\equiv2\bmod9}} .
\]
In the split type the same odd-order criterion
in TD3 implies a Hasse failure if this root
number is positive. Each signed dyadic window
of CT4 has a compatible fixed root-positive
congruence family in the split type.

Indeed reduction at \(2\) is good, and \(W_\infty=-1\).
The coefficient \(-216m\) has valuation three at
\(3\), with \(3\)-unit part \(-8m\equiv m\bmod9\).
Its local sign is negative exactly at \(m\equiv2\bmod9\)
by [Várilly-Alvarado, Lemma 4.1 and Proposition 4.3](https://msp.org/ant/2011/5-5/ant-v5-n5-p04-p.pdf).
At \(k\), valuation two gives the factor
\((-3/k)=(k/3)\). The product over the simple primes
of \(R/3\) is \((-1/|R/3|)=\operatorname{sgn}(R)\),
since \(R/3\equiv1\bmod4\). For \(k=1\) omit
the \(k\)-factor. This proves the sign formula.
For odd ray order, TD4 bounds the split-type
two-Selmer dimension by one. The parity theorem
gives even dimension at positive root sign, since
there is no rational two-torsion. Use
[Alpöge–Bhargava–Shnidman, Theorem 2.8](https://arxiv.org/pdf/2210.10730)
for this form of the parity theorem. Thus the
two-Selmer group is zero and the torsion argument
of TD3 proves the asserted obstruction.
For the imaginary sign choose \(p\equiv1,\ k\equiv2\bmod9\).
For the real sign choose \(p\equiv1,\ k\equiv7\bmod9\).
Both give \(m\equiv8\bmod9\), the required
root sign, and admissible fixed congruences in CT4.

**Corollary TD8.** In any acceptable congruence
family of cubefree \(m\equiv2\bmod8\) with \(3\nmid m\)
and with fixed \(m\bmod9\), the average size of
\(\operatorname{Sel}_2(E(m))\) is three. In
particular on the unmarked positive squarefree
progression \(m\equiv2\bmod72\) the finite ray
group of modulus \(\mathfrak p_2^2\) has odd
order for a lower proportion at least \(1/3\).

Here acceptable means the local conditions
include all integers of valuation at most one
outside a finite set of primes.
Choose \(s\) at the fixed \(3\)-adic residue and use
TD5 followed by
[Alpöge–Bhargava–Shnidman, Theorem 1.4](https://arxiv.org/pdf/2210.10730)
on the fixed coefficient \(d=-216s\) and parameter \(n=m\).
For a one-sided count, retain the indicated sign of
the real invariant in the proof of their Theorems
4.1 and 5.3. This restricts the real fundamental
domain to one sign component and scales the real
mass and the radicand volume by the same factor.
For the squarefree progression all primes
\(q>3\) in \(m\) have valuation one and their
root sign product is \(+1\) since \(m/2\equiv1\bmod4\).
The \(3\)-adic sign is negative, and
\(W_2W_\infty=-1\), so the root number is positive.
The parity theorem makes the Selmer size \(1\)
or at least \(4\). An average of \(3\) forces
a lower proportion at least \(1/3\) of size \(1\).
This is the type with \(3\) totally ramified
and TD2 gives the asserted odd ray order.

[Shnidman–Siad, Theorem 1.2](https://arxiv.org/pdf/2506.13749)
give averages of \(2\) and \(3/2\) for the **ordinary**
two-torsion in the wildly and tamely ramified pure
cubic families, respectively, ordered by cubefree
radicand. They also identify the obstruction
at the unique ramified split cubic place used in TD4.
Both those averages and TD8 involve wider families.
In a marked dyadic window of prime scale \(Z\), there
are at most \(O(Z^2)\) marked values of \(m\) at scale
\(Z^6\). Their definition prescribes the squarefree
quartic part and its coupling to the squared prime.
They are outside the acceptable congruence setting
of TD8. A uniform marked mean size strictly below
\(4\) on either root-positive signed subfamily of
TD7 would give a positive proportion of marked
Hasse failures by parity and CS2. Establishing that
mean, or finding infinitely many individual zero
groups in the marked family, remains open.

## Two certified prime-gap Hasse failures

**Theorem TD9.** The marked triples
\((71,101,131)\) and \((53,167,281)\), with
\(k=5\) and \(k=19\) respectively, give smooth
plane cubics
\[
 \begin{aligned}
  5 W^3&=10(A^3-B^3)+101AB(A+B),\\
 19 W^3&=38(A^3-B^3)+167AB(A+B).
 \end{aligned}
\]
Each has points over all completions of \(\mathbb Q\)
and has no rational point. They occur on the
imaginary and real quadratic norm sides, respectively.

For \((p,k)=(101,5)\), the signed quartic is
\(R=256285803=3\cdot23\cdot673\cdot5519\),
and \(m=4271430050\). The common tangent class
is nonzero: an integral cube root of the norm
element would have norm \(71\cdot101<R/4\),
which forces it to be rational in the imaginary
quadratic field, an impossibility for a nonrational
norm element. For \((p,k)=(167,19)\),
\(R=-10179717=-3\cdot3393239\) and
\(m=-2449918558\). The common class is nonzero
by [RV1](REAL_COHERENT_NORM_VIRTUAL_CLASS_RANK.md),
since the doubled irrational coefficient
\(114\) is less than \(|R|/4\).
The local persistence proof of CS2 applies in
each case, and the covering is the plane model
of [TC1](ABC_CUBIC_TANGENT_COVER_SURFACE.md).

Both \(m\) have split type at \(3\) and positive
root sign by TD7. The pure cubic field ordinary
class groups are certified quotients of
groups of orders \(27\) and \(9\), respectively.
In each case an explicit computed integral unit
has nontrivial residue modulo \(\mathfrak p_2^2\).
This proves that both \(\mathfrak p_2^2\) ray
orders are odd. Now TD7 gives the failure of the
Hasse principle. The following
[PARI/GP](https://pari.math.u-bordeaux.fr/dochtml/ref/General_number_fields.html)
replay checks the class group quotient with
`bnfcertify(b,1)`. A computed unit with the
required residue suffices, so no certification
of the full unit group is needed. Its coefficient
of \(\rho\) modulo \(2\) detects that residue.

```gp
checkmarked(p,k,rr,ff,cc)={
  my(h=6*k,R=3*p^4-6*p^2*h^2-h^4,m,b,u);
  if(!isprime(p-h)||!isprime(p)||!isprime(p+h)||!isprime(k),error("prime marks"));
  if(R!=rr||factor(abs(R))[,1]~!=ff||!issquarefree(R),error("squarefree quartic"));
  m=2*k^2*R/3;
  if(m%9!=8||-sign(R)*kronecker(k,3)!=1,error("root sign"));
  b=bnfinit(x^3-m,1);
  if(b.cyc!=cc||b.no%2!=1,error("computed class order"));
  u=b.fu[1];
  if(abs(norm(u))!=1||denominator(nfalgtobasis(b,u))!=1,error("integral unit"));
  if(lift(Mod(polcoeff(lift(u),1),2))!=1,error("dyadic residue"));
  if(bnfcertify(b,1)!=1,error("class group quotient certificate"));
  print([p,k,m]);
};
checkmarked(167,19,-10179717,[3,3393239],[3,3]);
checkmarked(101,5,256285803,[3,23,673,5519],[3,3,3]);
```
