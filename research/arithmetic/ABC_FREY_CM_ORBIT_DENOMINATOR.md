# A single CM orbit detects the primitive Frey monodromy denominator

Use the [fixed-level Frey family](ABC_FREY_FIXED_LEVEL_UNWEIGHTED_MONODROMY.md).
Write \(R=\operatorname{rad}(abc)\), \(B\) for the definite quaternion
algebra of discriminant \(5D\), and \(\mathcal R\) for an Eichler order
of level \(m=M/5\) in \(B\), where \(M\) is \(110\) or \(770\).
Let \(X=B^\times\backslash\widehat B^\times/\widehat{\mathcal R}^{\times}\)
be its finite ideal class set and \(s=|X|\). The primitive integral
Frey eigenfunction \(v:X\to\mathbb Z\) has
\(\sum_x v(x)=0\), \(\gcd_x v(x)=1\), and
\(H=\sum_x v(x)^2\), with all ideal class weights one.

There is one imaginary quadratic field \(K\) of discriminant
\(|d_K|\ll R\), an auxiliary prime \(q\ll\log(3R)\) with \(q\nmid R d_K\),
and a CM orbit of conductor \(q^n\) that maps onto the whole set \(X\).
For the ring class group \(G\) of that conductor, the pullback \(u:G\to
\mathbb Z\) of \(v\) is primitive. Let
\[
 P_\chi=\sum_{t\in G}u(t)\chi(t),\qquad
 B_\chi=\frac{P_\chi P_{\chi^{-1}}}{H}\qquad(\chi\in\widehat G).
\]
If \(d_{\rm tor}\) is the least positive integer making every
\(d_{\rm tor}B_\chi\) an algebraic integer, then
\[
     d_{\rm tor}\mid H,\qquad
     H\mid d_{\rm tor}\big(h_K q^{n-1}(q-\chi_K(q))\big)^2 .
\tag{CO1}
\]
Here \(h_K=|\operatorname{Pic}(\mathcal O_K)|\) and \(\chi_K\) is the
quadratic Dirichlet character of \(K\).
Thus the full norm denominator is detected by scalar toric products
at every coefficient prime outside that explicit class number.
The coefficient array \(u(t)u(t')/H\) has exact common denominator
\(H\) at every prime. For the exact-conductor scalar denominator
\(d_{\rm prim}\), the stronger comparison (CO8) loses only
\(q(i_5h_K(q-\chi_K(q)))^2\), where \(i_5\) is a bounded component
image order. Choosing \(K\) by Burgess gives
\(H\ll_\epsilon d_{\rm prim}R^{1/4+\epsilon}\).
Equations (CO11) and (CO12) identify the Petersson and Néron
normalizations of the same central values. At qualitatively large
depth, (CO14) recovers \(H/\gcd(H,i_5^2)\) from the lcm of
two exact-conductor denominators.

## An effective surjective CM orbit

Use the Chinese remainder theorem to choose a negative integer
\(d\) of size \(O(R)\) with \(d\equiv1\bmod8\), with a nonsquare unit
residue at each prime of \(5D\), and with a square unit residue at
each odd prime of \(m\). Take \(K=\mathbb Q(\sqrt d)\). Its fundamental
discriminant \(d_K\) divides \(d\), with the same splitting conditions
at all these places. The quaternion embedding criterion gives an
embedding \(K\to B\): \(K\) is imaginary and all finite ramified
primes of \(B\) are inert in \(K\). Choose \(g_0\in\widehat B^\times\)
such that \(\widehat K\cap g_0\widehat{\mathcal R}g_0^{-1}=
\widehat{\mathcal O}_K\). Locally at a ramified prime of \(B\), the
maximal order contains \(\mathcal O_{K,p}\). At a prime of \(m\)
use the diagonal maximal split order in the local Eichler order.
At the other split places of \(B\), use the action of \(K_p\) on
itself with lattice \(\mathcal O_{K,p}\).
[Voight, Proposition 14.6.7 and §30.3](sources/voight-quaternion-algebras-author.pdf).

Choose the least prime \(q\nmid Rd_K\). It is \(O(\log(3R))\) by
the lower Chebyshev bound for the product of primes up to a cutoff.
At \(q\), the split quaternion algebra gives the tree of homothety
classes of rank-two \(\mathbb Z_q\)-lattices. Based at \(g_{0,q}\),
identify the reference lattice with \(\mathcal O_{K,q}\).
At distance \(n\ge1\) the lattices have unique representatives of
the form
\[
                  \mathbb Z_q u+q^n\mathcal O_{K,q}
\]
with the line of \(u\) primitive modulo \(q\). If \(q\) is inert,
every such \(u\) is a unit of \(\mathcal O_{K,q}\). If \(q\) splits,
omit the two initial directions corresponding to the coordinate
axes of \(\mathbb Z_q\times\mathbb Z_q\). The remaining \(u\)'s
again have both components units. These lattices form a single
orbit of \(\mathcal O_{K,q}^\times\), with representative
\(\mathbb Z_q+q^n\mathcal O_{K,q}\).
Modify \(g_0\) at \(q\) only to give this representative \(g_n\).
Its optimal intersection is
\(\mathcal A_n=\mathbb Z+q^n\mathcal O_K\). The map
\[
        G=\widehat K^\times/(K^\times\widehat{\mathcal A}_n^\times)
        \ \longrightarrow\ X,\qquad
        t\longmapsto [t g_n]
\tag{CO2}
\]
therefore contains the projection of the entire allowed sphere,
by taking ideles that are units at \(q\) and one elsewhere.
[Kohel, §3.2, ideal modules and cyclic Hecke correspondences](sources/kohel-hecke-module-structure.pdf).
[Chida–Hsieh, §2.2, optimal embeddings and Gross points](sources/chida-hsieh-anticyclotomic-special-values.pdf).

To see that the allowed sphere hits every class, let \(A=T_q\)
be the symmetric adjacency operator, counting local directions
with multiplicities. Its constant row sum is \(q+1\). On the
degree-zero class module its eigenvalues have absolute value
at most \(2\sqrt q\) by the weight-two Jacquet–Langlands
correspondence and the Ramanujan bound. Let \(S_j(x,y)\) count
the endpoints \(y\) of paths of length \(j\) in the tree rooted
at \(x\) without immediate reversals. Here \(S_0=I\), \(S_1=A\),
\(S_2=A^2-(q+1)I\), and \(S_j=A S_{j-1}-qS_{j-2}\) for \(j\ge3\).
On the nonconstant space the polynomial at degree \(j\ge2\)
is the difference of the Chebyshev terms
\(q^{j/2}U_j(A/(2\sqrt q))\) and
\(q^{(j-2)/2}U_{j-2}(A/(2\sqrt q))\).
For \(j\ge1\) this gives the entrywise estimate
\[
     S_j(x,y)=\frac{(q+1)q^{j-1}}s+O((j+1)q^{j/2}).
\]
For a specified initial direction \(x\to z\), the same count
has main term \(q^{j-1}/s\) and error
\(O((j+1)q^{(j-1)/2})\). Indeed it is the whole count of length
\(j-1\) from \(z\) minus the count starting with the reverse
direction. Iterating this subtraction bounds the error by a
geometric sum of the full-sphere errors, with the \(S_0\) error
bounded by one. These counts take place in the tree before
projection, so the argument includes loops and repeated edges
in the class graph.
[Kohel, §§3.3 and 4.3, cusp module and Jacquet–Langlands description](sources/kohel-hecke-module-structure.pdf).
[Charles–Goren–Lauter, §4.2, Brandt eigenvalues and the Ramanujan bound](sources/charles-goren-lauter-ramanujan-quaternion.pdf).

Deleting at most two directions leaves main term at least
\((q-1)q^{n-1}/s\) and error \(O((n+1)q^{n/2})\), uniformly in the
target class. Every count is positive if
\(q^{n/2}\) exceeds a sufficiently large absolute constant
times \(s(n+1)\). In particular take
\[
 n=2+\left\lceil\log_q\!\big(C s^2(\log(2+s))^2\big)\right\rceil
\]
with a sufficiently large absolute \(C\). Since \(s\ll R\),
this gives \(q^n\ll R^2(\log(3R))^5\) and
\(|\operatorname{disc}(\mathcal A_n)|\ll R^5(\log(3R))^{10}\).
The orbit (CO2) is surjective and \(u\) is primitive.

## Inversion and the scalar denominators

By Skolem–Noether there is \(j\in B^\times\) implementing conjugation
on \(K\). There is \(a\in\widehat K^\times\) such that
\(w=g_n^{-1}a^{-1}j g_n\) normalizes \(\widehat{\mathcal R}\).
To verify this locally, write \(a^{-1}j\) in the \(K_p\)-linear
lattice model at each split place of \(B\), where \(j\) is
multiplication by a \(K_p\)-unit or nonzero \(K_p\)-scalar followed
by quadratic conjugation. At places outside \(5Dm\), the chosen
lattice is either \(\mathcal O_{K,p}\) or \(\mathcal A_{n,p}\), both
stable under conjugation, so remove the scalar. At primes of \(m\),
\(K_p\) is split and the two local lattices can be written as
\(\mathbb Z_p^2\) and \(\mathbb Z_p\oplus p\mathbb Z_p\) in the two
eigenlines. The interchange of the eigenlines followed by the
diagonal matrix with entries \(1,p\) exchanges the two lattice
classes and normalizes their common Eichler order. Absorb the
difference between this operation and \(j\) in \(a_p\).
At the division places the unique maximal order is normalized
by every nonzero element. The choices are units at almost all
places and define an idele \(a\).

The induced right action of \(w\) commutes with all good Hecke
operators outside a finite set. By multiplicity one it preserves
the Frey line, and being an isometry of the unweighted lattice
it acts on this line by a sign \(\sigma\).
For the image of \(a\) in \(G\), also denoted by \(a\), the identity
\(j t g_n=\bar t\,a g_n w\) gives
\[
                    u(a t^{-1})=\sigma u(t).
\]
Here \([\bar t]=[t^{-1}]\) because rational ideles act trivially
in the ring class group over \(\mathbb Q\).
Consequently \(P_\chi=\sigma\chi(a)P_{\chi^{-1}}\).
This is the finite inversion relation needed here for the
unregularized toric sums.
[Chida–Hsieh, Theorem 4.8, related inversion for their regularized theta elements](sources/chida-hsieh-anticyclotomic-special-values.pdf).

Fix a prime \(\ell\) and a place over it in a cyclotomic field containing
all character values, with \(\operatorname{ord}_\ell(\ell)=1\).
All the \(P_\chi\)'s are integral. Fourier inversion at a coordinate
\(u(t)\) that is an \(\ell\)-adic unit gives
\[
    0\le e_\ell:=\min_\chi\operatorname{ord}_\ell(P_\chi)
                            \le v_\ell(|G|).
\]
The minimum is the same at all places over \(\ell\), since Galois
permutes the characters. The inversion relation above says that
\(P_\chi P_{\chi^{-1}}\) has valuation \(2\operatorname{ord}_\ell(P_\chi)\)
at the same place. Thus
\[
 v_\ell(d_{\rm tor})=
    \max\{0,\lceil v_\ell(H)-2e_\ell\rceil\}.
\tag{CO3}
\]
The class number of \(\mathcal A_n\) is
\(|G|=h_K(q-\chi_K(q))q^{n-1}\): the additional factor is the index
of the local units of \(\mathcal A_{n,q}\) in \(\mathcal O_{K,q}^\times\),
with no global unit quotient because \(K\) has units only \(\{\pm1\}\).
In fact either \(K\) embeds in an order with that unit group, or
the split condition at \(11\) alone excludes the two exceptional
quadratic unit groups. This proves (CO1). Finally
\(\gcd_{t,t'}u(t)u(t')=1\) by primitivity, proving the exact
denominator assertion for the coefficient array.

## Characters of exact conductor

The scalar comparison can be restricted to characters of exact
conductor \(q^n\), with one additional bounded factor. For a prescribed directed edge ending at a class
\(x\), the allowed paths of length \(n-1\) ending with that edge
are counted by first reaching its tail and excluding immediate
reversal. The first count has main term at least
\((q-1)q^{n-3}/s\). The discarded paths, before imposing the
initial restrictions, have main term \(q^{n-3}/s\) by reversal
and the single-direction formula. The total error is
\(O(nq^{(n-2)/2})\). All such final directions occur with the
choice of large \(C\).
Reversal here respects individual edges: an oriented edge is
a choice of one of the \(q+1\) local neighbors at a representative,
and transport of that choice to a different representative is
independent of the transporter because every class stabilizer
acts trivially on local directions. If a transporter fixes a
local lattice class, a rational power of \(q\) first removes its
local homothety factor, after which its residual ambiguity is a
unit of the associated Eichler order, hence central.

Let \(G_n=G\) and let \(U\) be the kernel of \(G_n\to G_{n-1}\).
It has order \(q\). For the local part of (CO2), its orbits
are exactly the \(q\) children of each depth-\((n-1)\) vertex in
the allowed sphere. Indeed the torus description by the lines of
units modulo \(q^n\) reduces to the same description modulo
\(q^{n-1}\); the local unit index from these two orders is \(q\)
and the global units are the same. Therefore every pair of local
directions out of \(x\) occurs among a set of sibling children:
choose an incoming directed edge whose reverse differs from the
two specified outgoing directions. The positivity just proved
supplies that incoming edge. The class graph is connected and
nonbipartite because its only eigenvalue of absolute value \(q+1\)
is the constant eigenvalue. Thus any two classes are joined by a
walk of even length. Differences of endpoint values of consecutive
two-edge pieces generate the ideal of all coordinate differences.
We obtain the exact equality
\[
       \gcd_{t\in G,\,k\in U}\big(u(t)-u(tk)\big)
               =\gcd_{x,y\in X}\big(v(x)-v(y)\big)=i_5 .
\tag{CO4}
\]
For the second equality use the unweighted monodromy pairing
on the degree-zero lattice. Its pairing gcd \(i_5\) is the
component image order in (MN1). It has prime support at most
\(163\) and is absolutely bounded here, by (MN2).
[Integral projector and component image](ABC_FREY_MONODROMY_NORMS.md#exact-integral-projector-on-the-torus).

Fourier inversion of the differences on the left of (CO4) only
uses characters nontrivial on \(U\), which are exactly the
characters of conductor \(q^n\). Consequently for the minimum
\(e^{\rm prim}_\ell\) over those characters,
\[
  0\le e^{\rm prim}_\ell\le v_\ell(|G|)+v_\ell(i_5).
\]
The same inversion symmetry holds for this character set. If
\(d_{\rm prim}\) clears its scalar products, its valuation is
given by (CO3) with \(e^{\rm prim}_\ell\). In particular,
\[
 d_{\rm prim}\mid H,\qquad H\mid d_{\rm prim} (i_5|G|)^2 .
\tag{CO5}
\]
Thus primitive characters suffice to detect the full denominator
outside \(|G|i_5\).

## The local factor at the auxiliary prime

The characters of exact conductor give a direct local calculation
for the unregularized spherical vector, for both split and inert
\(q\). Let \(T=K_q^\times/\mathbb Q_q^\times\) and give the image
\(T_0\) of \(\mathcal O_{K,q}^\times\) volume one. Let \(f_n\) be the
unit spherical vector translated to the lattice \(\mathcal A_{n,q}\).
Use the unitary local representation of the Frey newform and write
\(\Phi_j\) for the spherical matrix coefficient at distance \(j\).
For \(t\in T_0\), the relative position of \(\mathcal A_{n,q}\) and
\(t\mathcal A_{n,q}\) has distance \(2(n-k)\), where \(k\) is the
length, capped at \(n\), of their common initial path. For \(1\le i\le n\)
the condition \(k\ge i\) holds precisely on the image \(T_i\) of
\(\mathcal A_{i,q}^\times\). Hence the matrix coefficient is constant
on each \(T_i\setminus T_{i+1}\). A character of conductor \(q^n\)
is nontrivial on \(T_{n-1}\) and trivial on \(T_n\), so integrating
over \(T_0\) leaves \((1-\Phi_2)\operatorname{vol}(T_n)\).

When \(q\) is inert, \(T=T_0\). When \(q\) is split, a coset of \(T_0\)
outside \(T_0\) translates the apartment by some \(k\ne0\) edges.
The two vertices are off the apartment by \(n\) edges, so their
distance is \(2n+|k|\) independently of the unit part. These cosets
have zero integral against the primitive character. The local integral
thus is
\[
 I_q(\chi)=\int_T \langle \pi_q(t)f_n,f_n\rangle\chi_q(t)\,dt
     =\frac{(q+1)^2-a_q^2}{q(q+1)}
           \frac{1}{(q-\chi_K(q))q^{n-1}} .
\tag{CO6}
\]
Here \(\Phi_2=(a_q^2-q-1)/(q(q+1))\) by the spherical sum on the
sphere of radius two. The split integral is absolutely convergent
by the spherical decay bound \(O((j+1)q^{-j/2})\), which follows
from the same Chebyshev polynomial divided by the sphere size.

For the usual normalized local factor with this Haar measure,
multiply (CO6) by
\[
 \frac{L(1,\chi_{K,q})L(1,\operatorname{Ad}\pi_q)}
      {\zeta_q(2)L(1/2,\pi_{K,q}\otimes\chi_q)} .
\]
The last local factor in the denominator is one because
\(\pi_q\) is unramified and \(\chi_q\) is ramified. The other factors
follow from the Satake parameters:
\[
 L(1,\operatorname{Ad}\pi_q)^{-1}
     =(1-q^{-1})\,\frac{(q+1)^2-a_q^2}{q^2}.
\]
The normalized integral is exactly
\[
                          q^{-n}(1-\chi_K(q)/q)^{-2}.
\tag{CO7}
\]
In particular the local normalization introduces no dependence on
\(a_q\) or on an ordinary stabilization.
[Chida–Hsieh, §3.4, local toric integrals for stabilized vectors](sources/chida-hsieh-anticyclotomic-special-values.pdf).

## Removing the conductor power from the denominator loss

Write \(\eta=\chi_K(q)\). The kernel of \(G_n\to
\operatorname{Pic}(\mathcal O_K)\) has a cyclic \(q\)-Sylow subgroup
\(C\) of order \(q^{n-1}\), whose subgroup of order \(q\) is \(U\).
Indeed the local unit quotient is
\(\mathcal O_{K,q}^\times/\mathcal A_{n,q}^\times\).
On its principal units the \(q\)-adic logarithm identifies the
quotient with the additive cyclic group
\[
       q\mathcal O_{K,q}/(q\mathbb Z_q+q^n\mathcal O_{K,q}).
\]
Here \(q\) is odd and unramified in \(K\), including when the local
algebra is split. The unit quotient at the first level has order
\(q-\eta\), prime to \(q\).

Put \(b=v_q(i_5)\), choose a generator \(c_*\) of \(C\), and
write \(a_*=q^{n-2}\). For \(t\in G_n\) consider
\[
                 F_t(X)=\sum_{j=0}^{q^{n-1}-1}
                                 u(t c_*^j)X^j .
\]
Reducing by \(\Phi_{q^{n-1}}(X)=1+X^{a_*}+\cdots+
X^{(q-1)a_*}\), the coefficients for \(0\le i<a_*\) and
\(0\le j<q-1\) are
\[
             u(t c_*^{i+j a_*})-u(t c_*^{i+(q-1)a_*}).
\]
They are \(q^b\)-divisible, and some \(t\) has a coefficient of
valuation exactly \(b\) by (CO4). Indeed in a sibling set
containing a pair with difference of valuation \(b\), at least
one difference with the last child has valuation \(b\).
For a primitive \(q^{n-1}\)st root \(\zeta\), the powers below
\(\deg\Phi_{q^{n-1}}\) form a \(\mathbb Z_q\)-basis of
\(\mathbb Z_q[\zeta]=\mathcal O_{\mathbb Q_q(\zeta)}\).
Thus for this \(t\),
\[
                         v_q(F_t(\zeta))<b+1 .
\]
The value \(F_t(\zeta)\) equals the \(\zeta\)-weighted character
sum over \(C\) at the coset \(tC\). Extending its character from
\(C\) to \(G_n\) and applying Fourier inversion over \(G_n/C\)
expresses \(|G_n/C|F_t(\zeta)\) as a sum of the \(P_\chi\)'s of
those extensions multiplied by roots of unity. All extensions
are of exact conductor \(q^n\). Hence
\[
                 e_q^{\rm prim}<v_q(i_5h_K)+1.
\]
At primes \(\ell\ne q\) use the previous bound
\(e_\ell^{\rm prim}\le v_\ell(i_5h_K(q-\eta))\).
In (CO3) an upper bound \(e_\ell^{\rm prim}<k+1\) for integral
\(k\) contributes at most \(2k+1\) to the valuation deficit.
We obtain
\[
        d_{\rm prim}\mid H,\qquad
        H\mid d_{\rm prim}\,q\big(i_5h_K(q-\eta)\big)^2.
\tag{CO8}
\]
The integer on the right that multiplies \(d_{\rm prim}\) is
independent of the tower depth.

## A smaller discriminant with the same splitting conditions

For every \(\epsilon>0\) the field \(K\) can be chosen with
\[
                       |d_K|\ll_\epsilon R^{1/4+\epsilon}.
\tag{CO9}
\]
To prove this, let \(R_o\) be the odd part of \(R\), and count
\(1\le a\le X\) with \(a\equiv7\bmod8\), \((a,R_o)=1\), and the
prescribed Legendre symbols for \(d=-a\) at every prime of
\(R_o\). Multiply the count by \(2^{\omega(R_o)}\) and expand
the product of the sign indicators. The empty term is
\(X\varphi(R_o)/(8R_o)+O(2^{\omega(R_o)})\) by inclusion–exclusion.
For every nonempty term let \(r_0>1\) be the product of the
odd primes whose Legendre characters occur. Their product is a
primitive character of squarefree conductor \(r_0\). The remaining
coprimality conditions are removed by inclusion–exclusion over
divisors of \(R_o/r_0\). After writing \(a=bt\) for each such
divisor \(b\), the progression modulo \(8\) becomes an interval
sum with a shifted argument, since \(8\) is a unit modulo \(r_0\).
The Burgess bound at a fixed integer \(r\ge2\) gives the combined
nonempty error, for any \(\delta>0\),
\[
       \ll_{r,\delta}
       X^{1-1/r}R^{(r+1)/(4r^2)+\delta}.
\]
Here the factors from summing over subsets and divisors are
\(R^{o(1)}\) and are absorbed in \(\delta\). Empty intervals
contribute zero. Take \(X=R^{1/4+\epsilon}\) and then take \(r\)
large enough that \(1/(4r)<\epsilon\), with \(\delta\) smaller
than the resulting power saving. Since \(R_o/\varphi(R_o)=R^{o(1)}\),
the positive main term dominates. The finite initial range can
use the original Chinese remainder construction. Passing from
\(d<0\) to its fundamental discriminant preserves the required
splitting conditions, because \(d\) is an odd unit at the specified
odd primes and \(d\equiv1\bmod8\).
[Hasanalizade–Lin–Luna Martínez–Martin–Treviño, introduction and Theorem 1.1 on cubefree Burgess bounds](sources/hasanalizade-lin-luna-martin-trevino-burgess-2511.17778.pdf).

The original proof of surjectivity works with this \(K\).
The same choice \(q\ll\log(3R)\) and depth bound give ring
class discriminant \(O_\epsilon(R^{17/4+\epsilon})\), absorbing
logarithms with a smaller discriminant slack. By the quadratic
class number formula and (CO8), with \(d_{\rm prim}\) formed for
this choice of \(K\), for every \(\epsilon>0\)
\[
                    H\ll_\epsilon d_{\rm prim}R^{1/4+\epsilon}.
\tag{CO10}
\]
In particular a uniform bound \(d_{\rm prim}\ll_\epsilon
R^{7/4+\epsilon}\) for every positive \(\epsilon\) would suffice
for abc. The loss in the reverse comparison remains explicitly
supported on \(i_5h_K q(q-\eta)\).

## Global central values and the Néron period

Write \(c_K=q^n\) for the ring conductor and \(L_\chi=L(1,f_E,\chi)\)
for the finite Rankin central value in classical normalization.
Here \(R\) is squarefree, \((d_K c_K,R)=1\), primes of \(m\) split,
and exactly the primes of \(5D\) are inert among primes of \(R\).
The definite algebra ramifies at those inert primes and infinity,
an even number of places. These are the local torus signs in
Cai–Shu–Tian, Lemma 3.1: at inert conductor primes the special
representation has an unramified quadratic twist and the local
ring character is trivial. At the other finite places the algebra
is split, and at infinity the weight is two and the torus character
is trivial. The orders \(g_n\mathcal R g_n^{-1}\) adelically
have the optimal intersections used in their formula. At level
primes the order is the diagonal split Eichler order or the
unique maximal division order. At \(q\) it is the stabilizer
order of \(\mathcal A_{n,q}\), as also checked by (CO6).
Changing from \(\widehat{\mathcal R}\) to its conjugate by \(g_n\)
takes \(t\mapsto[tg_n]\) to the natural CM map in their formula.

Theorem 1.10 of Cai–Shu–Tian now gives, with their omitted places
consisting of infinity alone and the unit factor equal to one,
\[
             L_\chi=\frac{8\pi^2(f_E,f_E)_{\Gamma_0(R)}}
                                {\sqrt{|d_K|}\,c_K}
                              \frac{P_\chi P_{\chi^{-1}}}{H}.
\tag{CO11}
\]
Indeed the test-vector norm in their ideal class pairing is
\(H\) and the class-period factor is \(P_\chi\). The relative
class group over \(\mathbb Q\) is the ordinary ring class group,
and the quadratic unit quotient is trivial.
[Cai–Shu–Tian, Theorem 1.10, its ideal class formula, Lemma 3.1 and Proposition 3.12](sources/cai-shu-tian-explicit-waldspurger.pdf).

Let \(P(E)=\frac{i}{2}\int_{E(\mathbb C)}\omega_E\wedge
\overline{\omega_E}\) be the Néron area. Put \(\mathcal V_\chi=
\sqrt{|d_K|}\,c_K L_\chi/(2P(E))\), and let \(B_E\) be the
classical modular degree divided by the square of the differential
multiplier, as in (SP1). Then
\[
 \mathcal V_\chi
       = B_E B_\chi
       = \kappa_D n_5\left(\prod_{p\mid D}n_p\right)
                                       P_\chi P_{\chi^{-1}} .
\tag{CO12}
\]
The classical optimal Manin constant is one in absolute value
because the Frey conductor is squarefree.
[Česnavičius, Theorem 1.2](sources/cesnavicius-manin-semistable.pdf).
Here \(n_p=v_p(\Delta_E)\) for odd \(p\mid R\) and
\(\kappa_D=\Gamma_D(n_5^0/n_5)/i_5^2\) in the notation of (SP2),
(SP7), and (MN1). Indeed (SP1) gives the first equality,
the degree transfer gives \(B_E=\Gamma_D\delta_{D,M}
\prod_{p\mid D}n_p\), and (MN1) gives
\(H=\delta_{D,M}i_5^2/n_5^0\). The rational correction
\(\kappa_D\) is supported on primes at most \(163\) and has
logarithmic height \(O(1+\omega(D))\).
[Classical degree and Shimura transfer](ABC_CUBIC_ADJOINT_SHIMURA_PROJECTORS.md).
[Monodromy degree identity](ABC_FREY_MONODROMY_NORMS.md#exact-integral-projector-on-the-torus).

Thus ordinary integrality of the Néron-normalized central values
away from those bounded coefficient primes follows already from
the Tate exponents and the integral period products. A bound
for \(H\) needs information about the Petersson-normalized
denominator in (CO11), or about the archimedean size of the
nonzero Néron-normalized values. For example, from (CO12) the
minimum valuation over the exact-conductor family at a prime
\(\ell>163\) outside \(h_Kq(q-\eta)\) is exactly
\[
                       v_\ell\!\left(n_5\prod_{p\mid D}n_p\right).
\]
The min applies at any fixed place over \(\ell\) because of
(CO3), (CO4) and the inversion identity.

## Removing the odd class-number loss by simultaneous CM surjectivity

The quantitative comparison (CO10) uses a bounded conductor.
For arbitrary depth there is an exact denominator comparison at
odd coefficient primes. Fix a quadratic field \(K\) satisfying the
splitting conditions above, an odd auxiliary prime \(q\nmid R d_K\)
inert in \(K\), and an odd prime \(\ell\ne q\).
For all sufficiently large \(n\) in this tower,
\[
   e^{\rm prim}_\ell=v_\ell(i_5),\qquad
   v_\ell(d_{\rm prim}(K,q,n))
      =\max\{0,v_\ell(H)-2v_\ell(i_5)\}.
\tag{CO13}
\]
There is no restriction on the \(\ell\)-part of \(h_K\).
The depth threshold here is qualitative.

Let \(A\) be the \(\ell\)-primary torsion of the inverse limit of
the ring class groups \(G_n\). The successive maps beyond
level one have kernel of order \(q\). Thus \(A\) maps
isomorphically to the \(\ell\)-Sylow subgroup at every level
\(n\ge1\), even when \(\ell\mid q+1\).
Fix one lift of every element of \(A\) to
\(\operatorname{Gal}(K^{\rm ab}/K)\).
These lifts are pairwise inequivalent modulo the \(q\)-rational
subgroup in Cornut–Vatsal. To see this, consider a difference
which has a \(q\)-rational lift. Its idele can be taken to have
component \(1\) at \(q\), by multiplying by \(K^\times\) and
\(\mathbb Q_q^\times\); neither changes its image in the ring
class tower. If its image has finite order \(r\), its \(r\)th
power at every level is a global element \(\alpha_n\) times an
order unit idele. The valuations of \(\alpha_n\) at every prime
of \(K\) are fixed, so the \(\alpha_n\)'s differ only by global
units \(\{\pm1\}\). At \(q\) it follows that their common global
generator is in
\(\bigcap_n\mathcal A_{n,q}^\times=\mathbb Z_q^\times\), and
hence is rational. The \(r\)th power of the original fractional
ideal is therefore rational, making the ideal itself invariant
under conjugation. Its square is a rational ideal, since
split exponents pair up, inert prime ideals are rational and
the square of a ramified prime ideal is rational. As the idele
at \(q\) is \(1\), its square is trivial in every ring class
group. An odd-primary element of such order is trivial.

Fix a \(\mathbb Z_q\)-basis \((1,\tau)\) of the local lattice
\(\mathcal O_{K,q}\).
In these coordinates the reference depth-\(j\) lattice has
matrix \(\operatorname{diag}(1,q^j)\).
Use the two outgoing directions given modulo \(q\) by the
lines of \((1,0)\) and \((1,1)\), and denote the two corresponding
local matrices by
\[
      b_j=\begin{pmatrix}1&0\\j&1\end{pmatrix}
                \operatorname{diag}(1,q),\qquad j=0,1.
\]
Extend \(b_j\) to an adele trivial away from \(q\).
Write \(U^* \subset \widehat{\mathcal R}^{\times}\) for the open
subgroup which at \(q\) stabilizes both specified lines and
elsewhere is unchanged. Right translation by \(b_j\) induces
a map from the finite class set \(X^*\) of \(U^*\) to \(X\).
The reduced norm of \(U^*\) is \(\widehat{\mathbb Z}^{\times}\),
since its subgroup at \(q\) is the inverse image of a split
diagonal group in a suitable basis. Thus the reduced norm
component set over \(\mathbb Q\) is a singleton.

For \(n\ge2\), put \(\kappa=1+q^{n-1}\tau\) at \(q\) and \(1\)
elsewhere. In the basis of the depth-\((n-1)\) lattice,
multiplication by \(\kappa\) reduces to
\(\begin{pmatrix}1&0\\1&1\end{pmatrix}\) modulo \(q\).
It follows that the right translates by \(b_0,b_1\) of
\([t g_{n-1}]\) in \(X^*\) map to the sibling classes
\([t g_n]\) and \([t\kappa g_n]\) in \(X\), respectively.
The same \(\kappa\) works after every idele translation, since
the torus is commutative. Its image in \(G_n\) generates
the order-\(q\) kernel of \(G_n\to G_{n-1}\).

Write \(b=v_\ell(i_5)\) and normalize the differences by putting
\(w(x)=(v(x)-v(x_0))/\ell^b\), for any fixed \(x_0\in X\).
This is integral by (CO4), with nonconstant reduction modulo
\(\ell\) on \(X\). On \(X^*\) put
\[
                   F(y)=w(y b_1)-w(y b_0).
\]
Its reduction is nonconstant: over a class \(x\), the local
group acts transitively on ordered pairs of distinct outgoing
directions. If every pair of neighbors had equal values modulo
\(\ell\), the even-walk argument of (CO4) would make \(w\)
constant. A pair with unequal values gives opposite nonzero
differences when exchanged, distinct for odd \(\ell\).

Apply Cornut–Vatsal, Corollary 2.10, with their collection of
changed ramification sets equal to \(\{\varnothing\}\), since
\(B\) is already definite. Take their compact open subgroup
\(U^*\), their translating Galois group to be the full group,
and the finite set of translations to be the lifts of \(A\).
The simultaneous map
\[
                       t\longmapsto
          ([t\sigma g_{n-1}])_{\sigma\in A}\quad\text{in }(X^*)^A
\]
is surjective for all sufficiently large \(n\). The component
constraint in the corollary is void over \(\mathbb Q\) here,
and the growing conductor makes the sequence leave every finite
exceptional set in their CM quotient, even at the larger
compact open level \(\widehat{\mathcal R}^{\times}\).
Choose a tuple so that \(\sum_{\sigma\in A}F([t\sigma g_{n-1}])\)
is nonzero modulo \(\ell\). If \(A\) is nontrivial, take a constant
tuple with one entry changed between two distinct \(F\)-values;
the constant tuple contributes zero modulo \(\ell\). If \(A\) is
trivial, take one nonzero value.
[Cornut–Vatsal, Definition 2.8, Theorem 2.9 and Corollary 2.10](sources/cornut-vatsal-cm-points-quaternion-2005.pdf).

The resulting sibling difference of traces of \(w\) over \(A\)
is nonzero modulo \(\ell\). Fourier inversion on \(G_n/A\)
now has unit denominator at \(\ell\). Only characters
nontrivial on the image of \(\kappa\) appear in that difference.
Thus some character of exact conductor \(q^n\), trivial on \(A\),
has a \(w\)-period which is a unit at the chosen place over
\(\ell\). Every exact-conductor character annihilates the constant
subtracted in defining \(w\). Its \(v\)-period is therefore
\(\ell^b\) times its \(w\)-period, and (CO3) proves (CO13).

## The coefficient prime 2 for a prime discriminant

Suppose now that \(d_K=-r\) for a prime \(r\equiv7\bmod8\),
and retain the required splitting at \(R\). Then \(h_K\) is odd
by genus theory. In ideal language, every class of order two
is fixed by conjugation, and Hilbert 90 gives a representative
ideal fixed by conjugation. Modulo rational ideals the only
possible generator of such ideals is the ramified prime over
\(r\), which is principal. Choose \(q\) inert with \(q\equiv1\bmod4\).
The \(2\)-Sylow of \(G_n\), for \(n\ge1\), has order two, by
the formula for \(|G_n|\). Its nontrivial element \(\sigma\)
is represented by a uniformizer idele of \(K\) supported
over \(r\). Indeed its square is rational up to local units.
If the idele were trivial in \(G_n\), its global generator
would have to be \(\sqrt{-r}\) up to sign, by the valuations
at all primes. The reduction of this generator at inert \(q\)
is outside \(\mathbb F_q\), so it is not an
\(\mathcal A_{n,q}\)-unit.

Keep \(U^*\) and \(F\) as above with \(\ell=2\). The function
\(F\) is again nonconstant modulo \(2\). If all ordered pairs
of distinct outgoing directions gave the same difference,
three directions at one vertex would force that difference
to be zero by addition. Then the even-walk argument would
make \(w\) constant.

At \(r\), the local component of \(g_{n-1}\) is fixed across
the \(q\)-tower, and conjugating the uniformizer idele by
\(g_{n-1}\) gives an adelic element \(b_r\) supported at \(r\).
The base lattice and its translate by \((b_r)_r\) are
neighbors, since the quadratic uniformizer has norm of
valuation one. Let \(U^{**}=U^*\cap b_rU^*b_r^{-1}\).
Its reduced norm is still \(\widehat{\mathbb Z}^\times\),
since the additional stabilizer at \(r\) is an edge
stabilizer of the split local tree. The two maps from
\(X^{**}\) to \(X^*\) given by the identity and right
translation by \(b_r\) realize every ordered \(r\)-edge.
The \(r\)-neighbor graph on \(X^*\) is connected by the
Brandt spectrum bound. Indeed the level at \(q\) in \(U^*\)
is Eichler of exponent two: stabilizing the two specified
neighbor lattices stabilizes the path between them. The
Ramanujan bound at good \(r\) then gives a one-dimensional
\(r+1\) eigenspace. Since \(F\) is nonconstant modulo \(2\),
some \(r\)-edge has distinct \(F\)-values at its endpoints.

Apply Cornut–Vatsal, Corollary 2.10, at \(U^{**}\) with just
the identity translation. For all sufficiently large \(n\)
some \(t g_{n-1}\) realizes this \(r\)-edge. In the torus
trace over \(\{1,\sigma\}\), the two parent classes in \(X^*\)
are \([t g_{n-1}]\) and \([t g_{n-1}b_r]\).
Their \(F\)-values therefore sum to a unit modulo \(2\).
This is the sibling difference for the \(2\)-Sylow trace.
Fourier inversion on the odd-order quotient proves
\(e_2^{\rm prim}=v_2(i_5)\), exactly as in (CO13).

## An exact comparison using two auxiliary primes

Choose \(K\) of prime discriminant \(-r\), with
\(r\equiv7\bmod8\), preserving the prescribed splitting at
all primes of \(R\). This is possible by the Chinese remainder
theorem and Dirichlet's theorem, choosing \(r\nmid R\).
Choose two distinct inert auxiliary primes \(q_1,q_2\nmid Rr\),
both \(1\bmod4\); Dirichlet's theorem supplies them from the
nonsquare residues modulo \(r\).
For sufficiently large \(n_1,n_2\) put
\[
                   d_*=\operatorname{lcm}\big(
                     d_{\rm prim}(K,q_1,n_1),
                     d_{\rm prim}(K,q_2,n_2)\big).
\]
Every odd coefficient prime is covered by (CO13) in at least
one tower, and \(2\) is covered in both. In either tower the
lower bound \(e^{\rm prim}_\ell\ge v_\ell(i_5)\) holds by the
constant subtraction argument. For
\(\ell\nmid h_K q_j(q_j+1)\) the elementary Fourier bound
already gives equality. Thus beyond the coefficient prime
\(2\), only odd primes dividing \(h_K(q_j+1)\) other than
\(q_j\) require a depth threshold from (CO13). Therefore
\[
                      d_*=\frac{H}{\gcd(H,i_5^2)} .
\tag{CO14}
\]
The factor \(\gcd(H,i_5^2)\) is absolutely bounded for the
fixed-level Frey family. Hence the sharp \(R^{2+\epsilon}\)
bound for this lcm at sufficiently large depths is equivalent
to full rational abc. There is also a one-orbit version if one
permits the choice of \(q\) to depend on \(H\), by choosing
one inert \(q\nmid RHr\) with \(q\equiv1\bmod4\).
At sufficiently large depth the individual \(d_{\rm prim}\)
then satisfies (CO14) in place of the lcm.

This removes the class-number loss for the unbounded-depth
comparison. The effective comparison (CO10) still bounds the
discriminant and depth in terms of \(R\). The proof of (CO14)
gives no uniform quantitative bound for the depths in terms of
\(R\), and the sharp Petersson-normalized denominator bound
remains open. An unbounded ratio of \(d_*\) to
\(R^{2+\epsilon}\) at a fixed positive exponent slack would
disprove abc because \(d_*\mid H\).
