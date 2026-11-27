# Pure cubic ray classes and the third tangent twist

The third cubic twist of the tangent progression has
a two-descent in a pure cubic field. For one of the
two possible local types at \(3\), its entire two-Selmer
group is a ray class character group with modulus
supported at \(2\). In the other type, a single
explicit paired condition at \(3\) is added.

## The exact local conditions

Take a squarefree signed tangent quartic
\(R=3p^4-6p^2h^2-h^4\) with \(h=6k\)
as in
[CS1](ABC_CUBIC_TANGENT_COMMON_SELMER_CLASS.md#three-local-images-in-a-common-quadratic-algebra).
Allow \(k=1\) as well. Set
\[
       m=2k^2(R/3),\qquad T=\mathbb Q(\rho),
              \quad \rho^3=m,\qquad \theta=6\rho .
\]
Then \(m\) is cubefree,
\(v_2(m)=1,\ m\equiv2\bmod8\), and
\(E_2:y^2=x^3-216m\).
The cubic is irreducible, so there is
no rational two-torsion.
Let \(\mathfrak p_2\) be the unique
prime above \(2\) in \(T\).

**Theorem TD1.** Identify the two-Selmer
group of \(E_2\) in
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
- If \(m\not\equiv-1\bmod9\), \(3\) is totally
  ramified and imposes no condition beyond the
  global norm. If \(m\equiv-1\bmod9\), put
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
has dimension \(\dim E_2(\mathbb F_q)[2]\).
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

At \(2\) use the good model with
\(x=4X,\ y=8Y+4\) from CS5.
Its group over \(\mathbb F_2\)
has order three. In the formal
subgroup, the original abscissa
is odd integral or has negative
even valuation. The same argument
as in CS5 gives the local image
spanned by \([1-6\rho]\).
The class is nontrivial: a
local point with \(x=1\) has a
formal parameter of valuation
one in the good model, and cannot
be twice a local point.
Alternatively the square criterion
modulo \(\rho^7\) in CS5 applies
unchanged since \(m\equiv2\bmod8\).
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

At \(3\), \(m\) has residue in
\(\{2,5,8\}\bmod9\) because
\(R/3\equiv p^4\bmod9\).
The unit cube criterion in
\(\mathbb Q_3\) gives the split
cubic algebra exactly at residue
\(8\). Otherwise the cubic
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

**Corollary TD2.** If \(m\not\equiv-1\bmod9\),
the two-Selmer group is canonically the
group of quadratic characters of the
ray class group of \(T\) of modulus
\(\mathfrak p_2^2\) with no real place
in the modulus. If \(m\equiv-1\bmod9\),
use the modulus obtained by also
including both primes above \(3\)
with exponent one. Impose the
paired local condition at those
two primes from TD1 on its
quadratic characters.

Here a quadratic character is identified
with the corresponding squareclass
of \(T\). Away from \(2,3\), the
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

In the totally ramified case at
\(3\), odd ray class number of
modulus \(\mathfrak p_2^2\) is
therefore a sufficient global
insolubility test for the
common-class covering.
The rank of \(E_2(\mathbb Q)\)
is then zero by two-descent,
and its rational torsion is
zero by the good fibre at
\(2\), absence of two-torsion,
and the three-division test
at a prime \(q>3\) dividing
\(R\), which exists for all
sufficiently large \(|R|\).
By CS2 the nonzero shared
class survives locally, so
the covering has no rational
point. In particular the
test holds if the ordinary
class number of \(T\) is odd
and a fundamental unit of
positive norm has nontrivial
residue modulo \(\mathfrak p_2^2\).
This follows directly from
the ray class exact sequence,
whose unit residue group
has order two. The example
\(m=124562\) of CS5 has
class group \([3,3]\) and
such a unit.

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
prime family. Results averaging
two-torsion over all cubic
fields do not give this
thin-family ray class parity
input.
