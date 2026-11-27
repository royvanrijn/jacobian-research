# Selected coherent Frobenius at mu-ordinary nonsplit primes

The relative \(\mu\)-ordinary polygon can force a directed coherent
Frobenius entry to be nonzero even if the return map vanishes.
The following criterion works in arbitrary relative rank.

## The rank of a consecutive coherent segment

**Lemma.** On a coefficient Frobenius orbit of a
\(\mu\)-ordinary \(p\)-divisible group with unramified
coefficient action, the rank of a consecutive coherent
Frobenius composite equals the smallest dimension of the
coherent Hodge quotients along the segment, including its two
endpoints. In particular it is a nonzero scalar map between
lines when no intermediate quotient is zero.

**Proof.** Use contravariant \(H^1\) and write the Frobenius
steps as \(\varphi_i:M_i\longrightarrow M_{p i}\) on the
oriented crystalline lattices over an unramified extension.
Let \(d_i\) be the dimension of the coherent Hodge quotient,
and \(r\) its minimum along the selected segment. The case
\(r=0\) is immediate. There are \(d_i\) unit elementary divisors
at step \(i\) and the rest are \(p\). The reduction of
\(\varphi_i\) factors through the Hodge quotient, and the induced
map to the next Hodge quotient is the coherent Frobenius.
For a complete coefficient orbit normalize
\(\wedge^r\varphi_i\) by \(p^{\max(0,r-d_i)}\) at each step.
These normalized maps are integral. This divided exterior
Frobenius also underlies the \(\mu\)-ordinary Hasse construction.
[Goldring–Nicole, Section 3](https://arxiv.org/abs/1305.6956)
· [Lin–Mantovan–Singhal, Section 2.4](https://arxiv.org/html/2406.09632v1#S2.SS4).

The \(\mu\)-ordinary polygon averages the sorted Hodge slopes
around the orbit of length \(f\). The sum of its first \(r\)
slopes is therefore \(f^{-1}\sum_i\max(0,r-d_i)\), exactly
the shift just removed from the exterior return. Hence this
return has slope zero. No power of its reduction can vanish,
as is seen by iterating over a finite field of definition until
the return is linear, when its integral matrix has a unit root.

At all steps of the selected segment and at its endpoint,
\(d_i\ge r\), so the normalized exterior step is undivided.
If the coherent composite had rank less than \(r\), its
exterior composite would be zero. The crystalline exterior
product along the segment would have image in the kernel of
projection onto the exterior power of the endpoint quotient.
The next undivided exterior step kills that kernel, forcing
a sufficiently long normalized return to vanish in the residue
field. The unit root excludes this. The coherent rank is at
most \(r\) by the quotient dimensions, giving equality.
\(\square\)

## A fixed cyclic cover and a positive-density moment test

**Proposition.** For the curve
\[
        C:\quad w^5=z^2(z-1)^2(z-3)^2
\]
over \(\mathbb Q\), there is a set of rational primes of
natural density \(1/4\) where the directed coherent map
from character four to character one by two Frobenius steps
is nonzero, while the character-one return is zero.
These primes give nonzero prescribed moments on the sphere
by an explicit radial term.

**Proof.** The smooth normalization has genus four by
Riemann–Hurwitz, with four tame branch points including
infinity. Its Jacobian has the rank-two action by
\(E=\mathbb Q(\zeta_5)\) from the deck transformation.
Rosati is complex conjugation since deck transformations
preserve the canonical polarization. Galois preserves
the field action and acts by the cyclotomic permutation.
The coherent dimensions at characters \(1,2,3,4\) are
\(1,2,0,1\), from the usual Kummer eigensheaves. Thus there
is one balanced conjugate pair.

The Jacobian is not potentially CM. Here are the exact
finite-field data used to verify this. At a split prime
\(p\ne 2,3,5\), take a primitive order-five character on
\(\mathbb F_p^\times\), and lift it by norm to
\(\mathbb F_{p^2}^\times\). The negatives of the two
sums of the character of \(z^2(z-1)^2(z-3)^2\), with
zero contribution at its roots, give the relative trace
\(a\) and the second power trace \(a_2\). The relative
polynomial is \(X^2-aX+(a^2-a_2)/2\). The sums compute
the individual primitive character factors of the point
count by the usual character projector. Indeed, a deck translate
of Frobenius fixes five points over a nonbranch value \(z\) or
none, according to its fifth-power residue character. Projection
of the fixed-point formula to a nontrivial deck character
removes the branch contributions. At the two primes
we obtain, writing \(\zeta=\zeta_5\) and \(\delta=(a^2-a_2)/2\),
\[
\begin{array}{c|l|l|l}
p&a&\delta&N_E(a^2-4\delta)\\\hline
11&-\zeta^2-2\zeta&7\zeta^3+9\zeta^2-3\zeta+3&5^2\,11^2\,29^2\\
61&3\zeta^3+2\zeta^2+2&24\zeta^3+16\zeta^2+39\zeta+72&5^2\,61^2\,181^2
\end{array}
\]
Both norms of the relative polynomials are irreducible over
\(\mathbb Q\). The two Frobenius fields have degree eight
and contain \(E\): the relative discriminants are nonsquares
and the four conjugates of the relative traces are distinct.
At \(p=11\) there is a valuation-one discriminant place
above 29, since 29 has residue degree two in \(E\). At
\(p=61\) there are two valuation-one places above 181,
as the discriminant vanishes at \(42,125\) among the four
primitive fifth roots modulo 181. In particular the degree-eight
field at 11 ramifies above 29 and the field at 61 is
unramified there.

The reductions are absolutely simple. If a ratio of their
eigenvalues is a root of unity of order \(j\), then
\(\varphi(j)\le64\), since a pair of eigenvalues generates a
field of degree at most 64. Put
\(N=\operatorname{lcm}\{j:\varphi(j)\le64\}\).
The characteristic polynomial of the \(N\)-th eigenvalue
powers is squarefree modulo 103 at 11, and modulo 137
at 61. These are checks in the polynomial quotient ring
over the corresponding finite field. The search for \(j\)
is finite with \(j\le8192\) by \(\varphi(j)\ge\sqrt{j/2}\).
Thus the eight eigenvalues stay distinct under every power.
Irreducibility persists under every finite residue extension,
and Tate's theorem identifies each geometric endomorphism
algebra of the reduction with its degree-eight Frobenius
field. A CM algebra of dimension eight on the characteristic-zero
Jacobian would inject into both fields by specialization,
impossible because of their different ramification.
[Exact replay with SageMath 10.9](verify_cm_newton_mu_cyclic_path.py)
· [Tate, Theorem 1](https://numdam.org/item/SB_1968-1969__11__95_0/).

Now apply the [one-balanced-pair theorem](CM_NEWTON_SEMILINEAR_PROJECTIVE_BOUNDS.md#averaged-hodge-types-and-a-single-balanced-pair)
over \(\mathbb Q\). The reduction is \(\mu\)-ordinary
relative to \(E\) outside a density-zero set. At
\(p\equiv3\bmod5\), the segment is \(4\to2\to1\) with
dimensions \(1,2,1\), so the lemma gives the nonzero map.
The character one immediately maps into character three,
whose coherent dimension is zero, so its return is zero.
These primes have density \(1/4\). For
\(k=(4p^2-1)/5\) and \(P(z)=(z-1)^2(z-3)^2\), the entry is
\[
             C_p=[z^{3k}]P(z)^k\ne0\pmod p .
\]
Indeed \([w^4/z^4]\) represents a generator for character
four: in the normalized generator
\(e_4=w^4/[z(z-1)(z-3)]\) its coefficient at \(e_4/z\)
is one. Its second Frobenius image is
\(w(z^2P(z))^k/z^{4p^2}\), giving the displayed coefficient.
For example the entry at \(p=13\), \(k=135\), is \(11\bmod13\).
[Cyclic coefficient rule](CM_NEWTON_SELECTED_EIGENLINES.md#a-moment-detecting-a-transition-rather-than-a-return).

To transfer the entry, take any integral phase \(g\) with
highest Laurent weight \(j+2i=5\) and angular polynomial
\(z^{-3}P(z)\) there. A regular example is
\[
 g=x t^3-(8+44y+96y^2+72y^3)t^5
\]
on \(t^2+2xy=1\). Smaller weights can be arbitrary.
Put \(m=k+1\), \(S=5k+1=4p^2\), choose \(p^R>mS\),
and add the radial term \(t^J\) with \(J=p^R-S\).
If exactly \(l\) factors of \((g+t^J)^m\) use this radial
term, every angular-constant denominator has the form
\[
        l p^R+\Delta,\qquad
        \Delta=j+1-lS,\quad 0\le j\le5(m-l).
\]
Here \(|\Delta|<p^R\), and \(\Delta=0\) holds only for
\(l=1,j=5k\). That term has even height \(p^R-1\)
and coefficient \(m C_p\). As \(m\equiv4/5\bmod p\),
this gives
\[
           p^R\mathcal A((g+t^J)^m)\equiv m C_p\ne0\pmod p .
\]
This family uses the fixed subleading cyclic-cover face
and increasing dominant radial terms. The ordinary
return or full ordinarity of the primitive factor is
unnecessary for the selected moment. \(\square\)

## Blache's generic ordinarity and prescribed covers

Blache's Theorem 1.1 gives generic \(\mu\)-ordinarity on the
fixed-monodromy Hurwitz space of tame abelian covers in characteristic
\(p\). Theorem 4.5 gives the cyclic character-orbit Newton polygon
and its Hasse-polynomial test for distinct branch points.
[arXiv:2607.21033v1](https://arxiv.org/html/2607.21033v1)

Check the monodromy tuple, character orbit and polygon conventions
against the actual prescribed cover or subfamily. For the fixed
quintic cover above, start with branch exponents \((2,2,2,4)\),
including infinity, and the residue class \(p\equiv3\bmod5\).
Choose an admissible branch chart before restricting the Hasse
polynomial. Generic nonvanishing on the ambient Hurwitz space leaves
nonvanishing on the prescribed subfamily or at its specified point
to be proved. A prime-density application for a fixed
characteristic-zero cover also needs a separate estimate on the
reductions at varying primes. The fixed-cover proof above uses its
own arithmetic input.
