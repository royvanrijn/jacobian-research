# A coefficient-prime trace operator for the signed Frey space

Let \(E/\mathbb Q\) be the Frey curve from the
[fixed-level family](ABC_FREY_FIXED_LEVEL_UNWEIGHTED_MONODROMY.md).
Its conductor is squarefree and equals \(R=\operatorname{rad}(abc)\).
Use \(R=DM\) with \(M=110\) or \(770\), and \(n_r=v_r(\Delta_E^{\min})\).
The argument below isolates all the eigensystems of conductor \(R\)
with the signs of \(E\) in two settings. One is the classical
\(D\)-new Hecke order. The other is the single toric lattice at \(5\).
The loss is bounded by a constant, \(2^{\omega(R)}\), and the product
of the \(n_r\) at the fixed level primes that can be removed.

The new local step concerns an oldform at \(p\) when the coefficient
prime is also \(p\). The ordinary projector discards the nonordinary
oldforms. A single Galois element acts trivially on every remaining
good ordinary representation and detects the valuation of the Tate
parameter on \(E\). The construction includes \(p=2\). The local argument and (OA5)
apply to any rational elliptic curve of squarefree conductor and
any choice of \(D\mid R\) in the classical space, with \(L_0\)
in place of \(24\). The fixed Frey
family is used for the uniform height comparison and the fixed
quaternionic level.

## Separating the valuation from the unit Kummer classes

Fix a prime \(p\), put \(K=\mathbb Q_p^{\mathrm{ur}}\), the algebraic
maximal unramified extension, and let \(L=K(\mu_{p^\infty})\).
Write \(k_x\) for the \(p\)-adic Kummer character of \(x\in K^\times\)
on \(G_L\), after choosing compatible generators of the roots of
unity. There exists \(\tau\in G_L\) such that
\[
 k_u(\tau)=0\quad(u\in\mathcal O_K^\times),\qquad
 v_p(k_p(\tau))=v_p(2).
 \tag{OA1}
\]

To prove this, let \(m\ge1\) and consider the restriction
\[
 K^\times/K^{\times p^m}\longrightarrow L^\times/L^{\times p^m}.
\]
Its kernel is \(H^1(\operatorname{Gal}(L/K),\mu_{p^m})\) by Kummer
theory and inflation–restriction. The cyclotomic group is
\(\mathbb Z_p^\times\), because the finite cyclotomic extensions
of \(\mathbb Q_p\) here are totally ramified. In additive notation
a cocycle \(b\) on this abelian group satisfies
\((\chi(z)-1)b(g)=(\chi(g)-1)b(z)\) for each \(g,z\). Take
\(z\) with \(\chi(z)=-1\). Hence the kernel is killed by \(2\).
The class of \(p\) in the quotient by the images of the units
has order at least \(p^{m-v_p(2)}\) for \(m>v_p(2)\).
Indeed if \(p^a\) has the same image as a unit, the valuation
of the square of their quotient is \(2a\), which must be
divisible by \(p^m\).

Kummer duality over \(L\) therefore gives an element killing all
unit radicals of degree \(p^m\) while its value on \(p^{1/p^m}\)
has exponent of valuation at most \(v_p(2)\). The assertion also
holds for the full subgroup of unit classes by inverse limits of
the finite Kummer pairings. These conditions define nested
nonempty closed subsets of \(G_L\) for \(m>v_p(2)\), using the
same nonvanishing test modulo \(p^{v_p(2)+1}\). Compactness gives
\(\tau\) at every depth. At \(2\), \(\sqrt2\in L\), so the
valuation is exactly \(1\). This proves (OA1).

## The action on weight-two local representations

On the Tate module of a multiplicative elliptic curve, the restriction
to \(G_L\) is a unipotent extension given by the Kummer character of
its Tate parameter \(q=p^{n_p}u\). Any unramified quadratic twist
disappears on \(G_L\). In a Tate basis, (OA1) gives
\[
                 \rho_{E,p}(\tau)=1+k_p(\tau)n_p e_{12}.
 \tag{OA2}
\]
[Tate uniformization and its torsion extension](sources/tate-nonarchimedean-uniformization.pdf).

For a weight-two newform \(g\) of trivial character with conductor
prime to \(p\), consider any coefficient place above \(p\) at which
\(a_p(g)\) is a unit. The ordinary quotient of the covariant
representation \(\rho_{g,p}\) is unramified and the subrepresentation
is cyclotomic times an unramified character. This also holds at \(2\).
In cohomological convention the good-reduction crystalline Frobenius
has polynomial \(X^2-a_p(g)X+p\), by the Eichler–Shimura relation
and comparison with the special fiber. Its two slopes are \(0,1\).
In the filtered module for cohomology the unit-root line has
Newton number zero and its Hodge number is nonnegative. Weak
admissibility forces Hodge number zero, so the line corresponds to
an unramified subrepresentation. Dualizing gives the claimed
quotient, and the determinant gives the other character.
[Darmon–Diamond–Taylor, Theorem 1.29 and Theorem 3.1(b), (f)](sources/darmon-diamond-taylor-fermat.pdf).
[Colmez–Fontaine, Théorème A](sources/colmez-fontaine-weakly-admissible-crystalline.pdf).

For completeness, every such \(\rho_{g,p}\) is trivial on the
\(\tau\) of (OA1). A stable \(\mathbb Z_p\)-lattice in the underlying
crystalline representation is the Tate module of a Barsotti–Tate
group over \(\mathbb Z_p\), by the equivalence for weights in
\(\{0,1\}\) in Tate-module convention. The maximal étale quotient
has the full unramified quotient as its rational Tate module.
Indeed the latter has a stable quotient lattice coming from an
étale Barsotti–Tate group, and the quotient map lifts by full
faithfulness and factors through the maximal étale quotient. The
connected part has a cyclotomic-times-unramified rational Tate
module and so is of multiplicative type by the same equivalence
and Cartier duality. Over \(\mathcal O_K\) the étale quotient and
the dual of the connected part are constant. At finite depth
\(p^m\), a lift of each generator of the étale quotient is a torsor
under a product of \(\mu_{p^m}\). The torsor classes are unit
Kummer classes since
\(H^1_{\mathrm{fppf}}(\mathcal O_K,\mu_{p^m})
=\mathcal O_K^\times/\mathcal O_K^{\times p^m}\).
Thus \(\rho_{g,p}(\tau)=1\) at every depth.
[Liu, Theorem 2.2.1 on the Barsotti–Tate correspondence, including \(p=2\)](sources/liu-barsotti-tate-kisin-characteristic-two.pdf).

## Removing the whole old subspace at one prime

Work on a squarefree-level cuspidal Hecke space or its integral
character lattice, with a rational Frey eigensystem \(f\) of full
conductor. At a split level prime \(r\), \(U_r\) denotes the
Iwahori Hecke operator. At a coefficient prime \(\ell\), complete
the integral operator order over \(\mathbb Z_\ell\).
Write \(t(\gamma)\) for the operator acting as the trace of
\(\gamma\) on each system, with the same trace on both old copies.
These operators lie in the completed order generated by the good
\(T_s\), by approximating \(\gamma\) with Frobenius traces using
Chebotarev as in [(HN2)](ABC_FREY_FULL_LEVEL_CONGRUENCES.md#trace-differences-at-a-prime-exactly-dividing-the-level).

When \(r\ne\ell\), the tame Tate inertia operator of (HN2) kills
the \(r\)-old subspace and evaluates on the Frey line to
\(n_r(\rho_{E,\ell}(\sigma))_{21}\) for a suitable inertia element.
When \(r=\ell=p\), let \(e_{\mathrm{ord}}\) be the projector
onto the \(U_p\)-unit part. It is integral in the completed full
order, since the Fitting decomposition on the finite reductions
gives \(e_{\mathrm{ord}}=\lim_{j\to\infty}U_p^{j!}\).
On a \(p\)-old pair from a form \(g\) of level prime to \(p\),
the two \(U_p\)-roots satisfy \(X^2-a_p(g)X+p=0\). No unit
root exists when \(a_p(g)\) is a nonunit. When \(a_p(g)\) is
a unit, the \(p\)-old unit-root part has the good ordinary
representation just considered. Use the operator
\[
     S_p=e_{\mathrm{ord}}\bigl(t(\tau\sigma)-t(\sigma)\bigr).
 \tag{OA3}
\]
Then \(S_p\) annihilates the entire \(p\)-old summand, through
all powers of the coefficient prime. On the Frey line
\(e_{\mathrm{ord}}\) is the identity since \(a_p(f)=\pm1\).
The evaluation of (OA3) is
\(k_p(\tau)n_p(\rho_{E,p}(\sigma))_{21}\), by (OA2).

Put \(L_0=\operatorname{lcm}(1,\ldots,163)\). If all lower-left
entries in the Tate basis were divisible by \(\ell^b\), then
the first vector modulo \(\ell^b\) would give a
Galois-stable cyclic subgroup of order \(\ell^b\). The
rational cyclic-isogeny bound gives \(\ell^b\le163\).
Hence at each squarefree level prime \(r\) one can choose an
operator which kills the \(r\)-old space and whose Frey scalar
\(\lambda_r\) is nonzero with
\[
 v_\ell(\lambda_r)\le v_\ell(L_0n_r)
                             +\boldsymbol 1_{\ell=r=2}.
 \tag{OA4}
\]
On a Frey curve with full rational \(2\)-torsion, \(L_0\) can be
replaced by \(24\). For odd \(\ell\), a rational cyclic
\(\ell^b\)-subgroup together with the two independent rational
\(2\)-subgroups gives a rational cyclic \(4\ell^b\)-isogeny on an
isogenous curve, by the same construction as in the
[residual image proof](ABC_FREY_JOINT_DIVISION_FIELDS.md#the-image-over-the-hecke-order).
The cyclic-isogeny list allows \(b>0\) only for \(\ell=3,b=1\).
At \(2\), given a stable cyclic subgroup \(C\) of order \(2^b\),
choose a rational order-two subgroup \(C'\) disjoint from \(C\).
The dual isogeny \(E/C'\to E\) maps \((E/C')[2]\) onto \(C'\).
The inverse image of \(C\) thus has order \(2^{b+1}\) and just
one nonzero point of order two. It is cyclic, so \(b\le3\).
Consequently for every \(\ell\ge5\) the Frey scalar in (OA4)
has valuation exactly \(v_\ell(n_r)\).
There is one operator for all competing old packets, independent
of their number or coefficient fields.

## One classical signed space and one signed toric lattice

In the setting of the
[classical Fourier minors](ABC_FREY_CLASSICAL_D_NEW_CONGRUENCE_MINORS.md#the-congruence-index-from-fourier-minors),
let \(r_{\mathrm{sign}}\) denote the same full Hecke congruence
index, restricted to the exact \(R\)-new space with all signs
of \(f\). The same Sturm bound and the same saturated
coefficient-minor formula apply to this quotient. For the
original \(D\)-new index \(r_D\),
\[
 r_{\mathrm{sign}}\mid r_D\mid
     2^{\omega(R)+1}\,r_{\mathrm{sign}}\prod_{r\mid M}24n_r .
 \tag{OA5}
\]
For the upper divisibility, take the product of the annihilators
from (OA4) at \(r\mid M\). Its rational image is contained in the
exact newspace. On that space every \(U_r\) for \(r\mid R\) acts
by a sign, so multiplication by
\(\prod_{r\mid R}(U_r+a_r(f))\) cuts out the sign pattern and
costs at most \(2^{\omega(R)}\). The product \(S\) is an element
of the completed full Hecke order. Lift a multiple of the rational
Frey idempotent from the signed quotient and multiply by \(S\).
The result is that same multiple times the Frey scalar of \(S\)
on the full space and vanishes on the complementary space.
Bounding the valuation of that scalar gives the upper divisibility.
The lower one follows by restriction.

For the fixed toric lattice \(X_5\) at \(5\), let \(Y_5\) be
its saturated intersection with the rational exact \(R\)-new
space of the Frey signs. Write \(P\) for the Frey projector and
\(B_5=d_{Y_5}(P)\) for its denominator on that single lattice.
The Hecke action at \(r\mid M/5\) is the integral \(U_r\) on the
definite Eichler lattice, and the local systems and their
\(r\)-old roots agree with those in (OA3) by Jacquet–Langlands.
At \(5D\) the systems already have the required conductor.
Applying the same annihilators for \(r\mid M/5\), followed by the
integral Atkin–Lehner sign sum of (HN4), gives
\[
 B_5\mid d_{X_5}(P)\mid
       2^{\omega(R)+1}\,B_5\prod_{r\mid M/5}24n_r .
 \tag{OA6}
\]
Indeed the composite operator \(S\) maps into \(Y_5\) over
\(\mathbb Z_\ell\) and obeys \(PS=\lambda P\), where \(\lambda\) is the scalar of
\(S\) on \(f\).
Thus \(\ell^{v_\ell(B_5)}PS\) is integral on \(X_5\).
Restriction gives the lower divisibility. The relation
\(d_{X_5}(P)=H/i_5\) and the absolute bound for \(i_5\)
give the corresponding comparison with the primitive norm \(H\).
[Fixed-level toric lattice and monodromy identity](ABC_FREY_MONODROMY_NORMS.md#exact-integral-projector-on-the-torus).
[Integral Hecke and Atkin–Lehner actions in (HN2) and (HN4)](ABC_FREY_FULL_LEVEL_CONGRUENCES.md).

In (OA5) and (OA6) the logarithmic loss is
\(O(\omega(R)+\log(\log c+2))\) since there are at most four
level primes in \(M\) and \(n_r=O(\log c+1)\).
Both \(r_{\mathrm{sign}}\) and \(B_5\) consequently have
logarithm
\[
 \log c+\log R+O\!\left(\frac{\log c}{\log(\log c+2)}\right).
 \tag{OA7}
\]
Thus either single-space sharp bound with exponent \(2+\varepsilon\)
for every \(\varepsilon>0\) is equivalent to full rational abc,
also on the fixed \(n=60\) image. At a coefficient prime \(\ell\ge5\) the
comparison with the corresponding ambient index is exact whenever
the removed level primes have \(\ell\nmid n_r\).

The remaining estimate concerns congruences entirely among systems
of exact conductor \(R\) and the same signs. The local annihilator
bounds how much congruence can be lost by removing all the other
systems, but imposes no upper bound on depth within the retained
space. The integral Fourier relation test now needs to separate
\(f\) only from the other members of this single signed newspace.


## Shorter coefficient vectors for a fixed sign character

For squarefree \(R\) write \(g_R=g(X_0(R))\) and
\(c_R=2^{\omega(R)}\), the number of cusps. For every fixed
Atkin–Lehner sign character the coefficient lattice is already
detected by
\[
       B_{\mathrm{sign}}(R)=1+
               \left\lfloor\frac{2g_R-2}{c_R}\right\rfloor .
 \tag{OA8}
\]
Here \(g_R\ge1\) in the present family. Detection means that for a
rational form in that sign space, integrality of the first
\(B_{\mathrm{sign}}\) coefficients implies integrality of the
whole expansion. It also means that congruences between integral
forms through this bound are congruences of the whole expansion,
at every coefficient prime and power. The result applies to
rational subspaces of the fixed sign space, including the signed
newspace and its Frey complement.

Fix a coefficient prime \(\ell\). On the Deligne–Rapoport model
over \(\mathbb Z_\ell\), use the invertible relative dualizing
sheaf. A weight-two cusp form gives the generic differential
\(f(q)\,dq/q\). A cusp has a smooth integral formal parameter,
with the standard \(q\) at infinity. The Atkin–Lehner
involutions extend to the model. At squarefree level they act
transitively on the cusps. For \(\ell\mid R\), each of the two
special-fiber components has \(c_R/2\) cusps. The involutions
for \(q\mid R/\ell\) act as the corresponding ones at level
\(R/\ell\) on each component, and \(W_\ell\) exchanges components.
These actions follow from taking the quotients by the indicated
cyclic levels. Away from \(\ell\) they commute with the
Frobenius and Verschiebung maps defining the two components.
[Deligne–Rapoport, VI §6, Théorème 6.9 and VII §2](sources/deligne-rapoport-elliptic-moduli-1973.pdf).

If the first \(B_{\mathrm{sign}}\) coefficients of an integral
differential vanish modulo \(\ell\), its reduction has order
at least \(B_{\mathrm{sign}}\) at infinity and, by sign
equivariance, at every cusp. At good reduction a nonzero
section has zero divisor of degree \(2g_R-2\), so this is
impossible. At bad reduction let \(g_0=g(X_0(R/\ell))\)
and let \(s\) be the number of geometric nodes in the
special fiber. The restriction of the dualizing sheaf to
each of the two smooth components has degree
\[
                    2g_0-2+s=g_R-1 .
\]
The \(c_R/2\) required zeros exceed this degree on each
component. Thus the reduction vanishes on both components
and is zero on the reduced special fiber. This argument uses
the sign eigenvalue after reduction, which is still a unit
at \(2\).

A rational differential whose first \(B_{\mathrm{sign}}\)
coefficients are integral can be multiplied by a minimal
nonnegative power of \(\ell\) to lie in the lattice of
global sections of the dualizing sheaf. If that power is
positive, the resulting section has all these coefficients
zero modulo \(\ell\). By the preceding argument its reduction
is zero. The exact sequence for multiplication by \(\ell\)
on the relative dualizing sheaf makes it divisible by \(\ell\)
inside the lattice, contradicting minimality. Hence the
original form is integral and all of its coefficients are
\(\ell\)-integral. Applying the same argument to
\(\ell^{-k}(f_1-f_2)\) gives the congruence assertion.

Consequently \(B_{\mathrm{sign}}\) replaces the full Sturm
bound in the saturated Fourier-minor presentation of
\(r_{\mathrm{sign}}\). The dual coefficient relation can be
formed from the same prefix. Indeed the duality with the
signed Hecke order identifies the integral cusp lattice with
the dual of that order. Saturation of the shorter coefficient
map means exactly that the \(T_n\) through that prefix
span the quotient order additively. The bound is at most
\(1+R\prod_{p\mid R}(1+1/p)/(6c_R)\) by Riemann–Hurwitz
for the map to \(X_0(1)\).

### Finite check of the signed coefficient index

At level \(106=2\cdot53\), the exact newspace has four rational
packets. The two with \((U_2,U_{53})=(1,-1)\) have initial
coefficient vectors \((1,1,1,1,0,1)\) and
\((1,1,-2,1,3,-2)\). At the Sturm bound \(27\), PARI gives
\(r_{\mathrm{sign}}=3\) for the first packet by Smith normal form.
The full level index is \(12\), using the saturation of the
complementary span cut out by \(T_q-a_q(f)\) for good primes
\(q\le23\). This checks the restriction on a packet with a
nonzero remaining congruence.

Here (OA8) needs only six coefficients, and the resulting two-column
minor gcd is again 3.
