# A dual Selmer bound for the full-level Frey congruence denominator

Let \(E/\mathbb Q\) be the Frey curve of a primitive positive
\(a+b=c\) with at least two odd radical primes, \(R=\operatorname{rad}(abc)\),
and \(N=N_E\). For each \(\ell\ge5\), use the local signed Hecke order,
Frey projector denominator \(\ell^{e_\ell}\), and conormal groups
\(C_\ell,C_\ell^{\rm str}\) of the
[strict adjoint comparison](ABC_FREY_STRICT_ADJOINT_CONGRUENCES.md).
Write \(Q_\ell\) for the set of odd primes \(q\mid N\) distinct from
\(\ell\). Put
\[
 t_q=v_\ell(v_q(\Delta_E)),\qquad d_q=v_\ell(q-1),\qquad
 s_q=s_{q,\ell}.
\]
Here \(s_q\) is the scalar depth of the local Tate representation
in (SA2)–(SA4), so \(0\le s_q\le\min(t_q,d_q)\).
Let \(\beta_\ell\) be the maximum of zero, all \(t_q+d_q-s_q\) for
\(q\in Q_\ell\), and \(v_\ell(-v_2(j(E)))\) if \(E\) has potentially
multiplicative reduction at \(2\).

The bound uses the finite groups \(D^0_{\ell,m}\) defined below.
For each \(\ell\ge5\) with \(e_\ell>0\), take \(m=e_\ell+\beta_\ell\).
Then
\[
 \ell^{e_\ell+\sum_{q\in Q_\ell}d_q}
        \ \le\ |D^0_{\ell,m}|\,\ell^{4\sum_{q\in Q_\ell}s_q}.
 \tag{DA1}
\]
The group \(D^0_{\ell,m}\) is a group of equivariant twisted-adjoint characters of a
ray class group over the field
\(\mathbb Q(E[\ell^m])^{\{\pm I\}}\), cut out by local conditions at
\(2\) and \(\ell\). Its modulus is supported only above \(\ell\).
The factor on the right of (DA1) has total logarithm
\(O(\log c/\log(\log c+2))\) across coefficient primes by (SA5).
The factors \(\ell^{d_q}\) on the left give the full \(\ell\)-part
of each \(q-1\). In particular there is no local error at a prime
of scalar depth zero, even when \(\ell\mid v_q(\Delta_E)\).

## Local conditions at the coefficient prime

For \(m\ge1\) put \(B_m=\mathbb Z/\ell^m\),
\(M_m=\operatorname{ad}^0(E[\ell^m])\), where \(\operatorname{ad}^0\)
means trace-zero \(B_m\)-linear endomorphisms, and
\(M_m'=\operatorname{Hom}(M_m,\mu_{\ell^m})\simeq M_m(1)\)
by the trace pairing. Let \(L_{\ell,m}\subset H^1(\mathbb Q_\ell,M_m)\)
be the following local condition.

- If \(\ell\nmid N\), use the fixed-determinant weight-\(\{0,1\}\)
  Fontaine–Laffaille condition, with signs on the weights chosen
  to match the Tate module convention.
- If \(\ell\mid N\), \(E\) is multiplicative at \(\ell\).
  In a Tate basis modulo \(\ell^m\), take the image of
  \(H^1(\mathbb Q_\ell,U_m)\), where \(U_m=B_m e_{12}\simeq B_m(1)\)
  is the upper root line inside \(M_m\).

In both cases the signed conormal restriction lies in this local
condition at every finite depth at which a derivation is defined,
and
\[
                  |L_{\ell,m}|=\ell^m|H^0(\mathbb Q_\ell,M_m)|.
 \tag{DA2}
\]

For good reduction, the characteristic-zero packet at \(\ell\)
comes from the good-reduction modular Jacobian and is crystalline
with the two indicated weights. The representation over
\(A_\ell^2\) is a stable \(\mathbb Z_\ell\)-lattice in the direct
sum of these crystalline factors. A square-zero derivation map
\(A_\ell\to B_m[\epsilon]/(\epsilon^2)\) need not be surjective.
Nevertheless a finite direct sum of \(A_\ell^2\) maps onto the
associated rank-two deformation over \(B_m[\epsilon]/(\epsilon^2)\):
take a finite set of \(A_\ell\)-generators of the coefficient
ring and multiply the image of \(A_\ell^2\) by them. The quotient
has the required Fontaine–Laffaille type, using closure under
subquotients in this weight range and the exact, fully faithful
functor. The functor transports the coefficient-ring action
and preserves freeness and the ranks of the filtration steps.
[Fontaine–Laffaille, §6 and §8](https://www.numdam.org/item/10.24033/asens.1437.pdf).
[Booher, Facts 4.10 and 4.12 and Lemmas 4.15 and 4.18](https://arxiv.org/abs/1612.04237).

For the local count in good reduction, use a basis adapted to the
rank-one filtration of the Fontaine–Laffaille module, arranging
the two weights as \(0,1\) by duality if necessary. Frobenius on
the unfiltered generator and divided Frobenius on the filtered
generator give the two columns of an invertible matrix \(P\)
over \(B_m\). In a square-zero lift a lifted adapted basis also
exists. Its matrix is \(P+\epsilon X\). The determinant condition
is \(\operatorname{tr}(P^{-1}X)=0\), using the exterior square,
which remains in the Fontaine–Laffaille weight range. Since
the base local field is \(\mathbb Q_\ell\), there is no nontrivial
Witt-vector Frobenius on the base coefficients. The matrices
\(X\) form a free \(B_m\)-module of rank three. The trace-zero
infinitesimal changes of adapted basis have rank two. Scalar
changes of basis act trivially, and the stabilizer in the
trace-zero changes is \(H^0(\mathbb Q_\ell,M_m)\) by full
faithfulness. The orbit calculation proves (DA2). This is an
additive quotient of linear maps, so the local condition is a
subgroup. In particular the count does not require a
nonsplit ordinary residual representation or exclude supersingular
reduction.

For multiplicative reduction the \(\ell\)-new special factors
of the packet have the same unramified quotient character
\(\mu\), with Frobenius value \(U_\ell=\pm1\). Their diagonal
characters are \(\mu\chi_\ell,\mu\). To descend the line to
\(A_\ell\), choose an inertial element with \(\chi_\ell-1\)
a unit. Its eigenprojector for \(\chi_\ell\) is integral over
\(A_\ell\), since the eigenvalues in every factor are
\(1,\chi_\ell\). The image is a free rank-one summand. It
is stable and has the specified character by checking the
matrix identities in the faithful characteristic-zero product.
The quotient likewise has character \(\mu\). Thus after
choosing a basis on the two summands and augmenting, the
first-order classes vary only in the upper root direction.
The quotient \(M_m/U_m\) has a trivial diagonal line and an
inverse-cyclotomic quotient with zero invariants, so
\[
 |H^0(\mathbb Q_\ell,M_m/U_m)|=\ell^m,\qquad
 H^0(\mathbb Q_\ell,U_m)=0.
\]
Kummer theory gives \(|H^1(\mathbb Q_\ell,U_m)|=\ell^{2m}\)
from the valuation and principal units. The kernel of its map
to \(H^1(\mathbb Q_\ell,M_m)\) has order
\(\ell^m/|H^0(\mathbb Q_\ell,M_m)|\) by the invariant sequence,
proving (DA2) also in this case.
[Wiles, §2, the local ordinary characters and \(U_\ell\)](https://jontalle.web.engr.illinois.edu/Public/Wiles-Fermat.95.pdf).

## Finite strict conditions and the duality identity

Define \(S_{\ell,m}\) using \(M_m\), the condition \(L_{\ell,m}\)
at \(\ell\), zero at every \(q\in Q_\ell\), and unramified
conditions at every other finite prime. Let \(D_{\ell,m}\)
be its Tate-orthogonal dual with \(M_m'\). The dual condition
is all of \(H^1(\mathbb Q_q,M_m')\) at \(q\in Q_\ell\)
and unramified at every finite prime outside \(Q_\ell\cup\{\ell\}\).
At residue characteristic different from \(\ell\), the
unramified subgroups for \(M_m,M_m'\) are orthogonal by
inflation from the procyclic quotient. Their sizes are the
respective \(H^0\) sizes, so they are exact annihilators by
local Euler characteristic and Tate duality, including for
ramified coefficients.

For \(q\in Q_\ell\) use the Tate matrices \(T=1+n e_{12}\)
and \(F=\left(\begin{smallmatrix}q&u\\0&1\end{smallmatrix}\right)\)
from (SA2). If \(Z=\left(\begin{smallmatrix}x&y\\z&-x\end{smallmatrix}\right)\),
then commuting with both is equivalent to
\[
        nz=uz=(q-1)z=0,\qquad
        nx=0,\qquad 2ux-(q-1)y=0.
\]
The \(z\)-factor is killed by \(\ell^{s_q}\). The two Smith
valuations for the \(x,y\) system are \(s_q,t_q+d_q-s_q\):
the gcd of its entries has valuation \(s_q\) and its
determinant has valuation \(t_q+d_q\). Thus for all \(m\),
\[
 \big|H^0(\mathbb Q_q,M_m)\big|
       =\ell^{2\min(m,s_q)+\min(m,t_q+d_q-s_q)}.
 \tag{DA3}
\]
In particular, its saturated order is \(\ell^{t_q+d_q+s_q}\)
at \(m\ge t_q+d_q-s_q\). At this depth the subgroup killed
by \(\ell^{t_q}\) has order \(\ell^{t_q+2s_q}\).

For each \(b\ge1\) and \(m\ge b+\beta_\ell\), there is an
injection \(C_\ell^{\rm str}[\ell^b]\hookrightarrow S_{\ell,m}\).
Indeed \(H^0(\mathbb Q,W_\ell)=0\) by full residual image, so
the finite global representatives map injectively into
\(H^1(\mathbb Q,W_\ell)\). Choose the representative from the
derivation already killed by \(\ell^b\). At \(q\in Q_\ell\),
strictness over \(W_\ell\) expresses this representative locally
as a coboundary from \(Z\in W_\ell\), with \(\ell^b Z\) invariant.
By the Smith calculation the exponent of
\(H^0(\mathbb Q_q,W_\ell)\) divides \(\ell^{t_q+d_q-s_q}\),
so \(Z\) is killed by \(\ell^m\). This makes the representative
strict at finite depth. At \(2\), for finite prime-to-\(\ell\)
inertia the coboundary in the proof of dyadic unramifiedness
can be chosen by averaging at depth \(b\). For potentially
multiplicative reduction, twist locally to the Tate
representation. A matrix giving the inertia coboundary
can have its upper entry discarded since that direction
is invariant under the Tate inertia. The other two entries
are killed by \(\ell^{b+v_\ell(-v_2(j(E)))}\), from their
commutator with the Tate inertia generator. Unramifiedness
at every other prime away from \(N\ell\) holds already in
the signed Hecke representation.
[Strict local calculation and dyadic unramifiedness](ABC_FREY_STRICT_ADJOINT_CONGRUENCES.md).

Global \(H^0\) of \(M_m\) and \(M_m'\) is trivial. For \(M_m\)
use full image and trace zero, and for \(M_m'\) use the
scalar \(2I\), which acts by \(4\) on the twist. For complex
conjugation \(M_m\) has a rank-one invariant summand, since
the two eigenvalues on \(E[\ell^m]\) are \(1,-1\). The finite
Greenberg–Wiles formula then cancels the \(\ell^m\) of (DA2)
with the real place. Every unramified condition contributes
a ratio of one. This gives the exact identity
\[
      |S_{\ell,m}|=
       \frac{|D_{\ell,m}|}{\prod_{q\in Q_\ell}|H^0(\mathbb Q_q,M_m)|}.
 \tag{DA4}
\]
[Wiles, Proposition 1.6, the finite Greenberg–Wiles ratio](https://jontalle.web.engr.illinois.edu/Public/Wiles-Fermat.95.pdf).

## Descent at the same division layer and removal of tame ramification

Let \(F_{\ell,m}=\mathbb Q(E[\ell^m])^{\{\pm I\}}\), equivalently
the compositum of the projective division field and
\(\mathbb Q(\mu_{\ell^m})\). The quotient of
\(\operatorname{GL}_2(B_m)\) by \(\{\pm I\}\) is its full Galois
group \(\Gamma_{\ell,m}\). It acts on \(M_m'\) by
\(\det\otimes\operatorname{Ad}\). In particular the central
class of \(2I\) acts as \(4\); cohomology of the finite group
with these coefficients vanishes in all degrees since the action
of a central group element on cohomology is trivial while
\(4-1\) is invertible. Inflation–restriction therefore gives
\[
 H^1(\mathbb Q,M_m')\ \simeq\
    \operatorname{Hom}_{\rm cont,\Gamma_{\ell,m}}
                     (G_{F_{\ell,m}},M_m').
 \tag{DA5}
\]
The isomorphism uses both \(H^1(\Gamma_{\ell,m},M_m')=0\)
and \(H^2(\Gamma_{\ell,m},M_m')=0\).
[Full Frey Galois image, (DF3)](ABC_FREY_JOINT_DIVISION_FIELDS.md).

For \(q\in Q_\ell\), choose a prime of \(F_{\ell,m}\) above \(q\).
Tate uniformization shows that its relative pro-\(\ell\)
inertia in the absolute local group is generated in tame
coordinates by \(\tau^{\ell^r}\), where \(\tau\) has image
\(T\) and \(r=\max(m-t_q,0)\). The wild and prime-to-\(\ell\)
tame parts make no contribution to an \(\ell\)-primary cocycle
on inertia, because they act trivially. For any local
cocycle \(\xi\) with coefficients in \(M_m'\), put
\(v=\xi(\tau^{\ell^r})\). The norm formula gives
\[
 v=\left(\sum_{i=0}^{\ell^r-1}\operatorname{Ad}(1+i n e_{12})\right)
      \xi(\tau)=\ell^r\,\xi(\tau).
 \tag{DA6}
\]
The coefficients of \(i\) and \(i^2\) contain \(n\) and \(n^2\),
respectively. The sums of \(i\) and \(i^2\) are divisible by
\(\ell^r\), since \(\ell\ge5\), and \(n\ell^r=0\) in \(B_m\).
The cyclotomic twist is trivial on inertia.

By (DA5), the homomorphism over \(F_{\ell,m}\) is equivariant
under the original decomposition group. Inertia acts
trivially on the tame inertia quotient, so \(v\) is fixed
by \(\operatorname{Ad}(T)\). Frobenius sends the tame generator
to its \(q\)-th power, so it acts on \(v\) by multiplication
by \(q\). The determinant twist takes precisely this value,
leaving \(\operatorname{Ad}(F)v=v\) on the underlying adjoint
module after identifying the cyclotomic factor with \(B_m\).
Thus (DA6) bounds the image on this inertia generator by
\[
           v\ \in\ H^0(\mathbb Q_q,M_m)[\ell^{\min(m,t_q)}].
 \tag{DA7}
\]
At primes conjugate to the chosen one the values are determined
by the same homomorphism and equivariance. In particular,
the possible inertia at all primes over a given rational \(q\)
costs a single factor of size at most \(\ell^{t_q+2s_q}\)
when \(m\ge\beta_\ell\).

Define \(D^0_{\ell,m}\subset D_{\ell,m}\) by requiring that its
homomorphism in (DA5) be unramified outside the places over
\(\ell\). The only extra restrictions are the tame restrictions
at \(q\in Q_\ell\), so by (DA7), for \(m\ge\beta_\ell\),
\[
       [D_{\ell,m}:D^0_{\ell,m}]
             \ \le\ \prod_{q\in Q_\ell}\ell^{t_q+2s_q}.
 \tag{DA8}
\]
Consequently (DA4) gives
\[
           |S_{\ell,m}|\ \le\
                    |D^0_{\ell,m}|\ell^{\sum_{q\in Q_\ell}(s_q-d_q)} .
\]
Take \(b=e_\ell\) and \(m=b+\beta_\ell\). The conormal group
is killed by \(\ell^{e_\ell}\), so the finite strict injection
and (SA1) prove (DA1).

For a completely finite realization, let \(e_w\) be the
absolute ramification index at a place \(w\) of \(F_{\ell,m}\)
over \(\ell\), and use the ray modulus
\[
           \mathfrak m_{\ell,m}
           =\prod_{w\mid\ell}\mathfrak p_w^{e_w m+
                                           \lfloor e_w/(\ell-1)\rfloor+1}.
 \tag{DA9}
\]
The local logarithm on principal units of depth greater than
\(e_w/(\ell-1)\) shows that every character killed by \(\ell^m\)
vanishes on the displayed depth. Global class field theory
identifies \(D^0_{\ell,m}\) with the subgroup of
\(\operatorname{Hom}_{\Gamma_{\ell,m}}
(\operatorname{Cl}_{\mathfrak m_{\ell,m}}(F_{\ell,m}),M_m')\)
whose descended classes satisfy the two local tests. At
\(\ell\) use \(L_{\ell,m}^{\perp}\). In the multiplicative
case this is explicitly the kernel of the map
\[
 H^1(\mathbb Q_\ell,M_m')\longrightarrow
 H^1(\mathbb Q_\ell,U_m^D)=H^1(\mathbb Q_\ell,B_m),
\]
induced by restriction of the dual module to the upper root
line \(U_m\). The formula follows from the compatibility of
local Tate pairings with the inclusion of \(U_m\).
At \(2\) require
unramifiedness over \(\mathbb Q_2\) on the descended class,
since unramifiedness over the division field alone only
tests a subgroup of inertia. At the other good primes the
division field is itself unramified, and at \(Q_\ell\) the
dual condition over \(\mathbb Q_q\) imposes no restriction.

## The remaining uniform estimate

Let \(\delta_1\) be the classical modular degree and put
\(\mathcal E(E)=\{\ell\ge5:e_\ell>0\}\). Define the normalized
product
\[
    \mathcal D(E)
      =(\delta_1)_{\{2,3\}}
          \prod_{\ell\in\mathcal E(E)}
             \frac{|D^0_{\ell,e_\ell+\beta_\ell}|}
                  {\prod_{q\mid N,\ q\ {\rm odd}}(q-1)_\ell},
 \qquad (q-1)_\ell=\ell^{v_\ell(q-1)}.
 \tag{DA10}
\]
A bound \(\mathcal D(E)\ll_\varepsilon R^{2+\varepsilon}\)
would prove full abc by (DA1) and the signed projector degree
comparison, with the same low-support treatment as in (SA7).
In particular, an abc failure would force this normalized
product to violate the sharp bound. The numerator involves
the entire selected dual space, so large normalized values
alone would still need to be traced back to the Frey
congruence denominator to give a counterexample. The local
calculation supplies the factors \(q-1\) at the coefficient
primes in \(\mathcal E(E)\). Summing over a larger set of
coefficient primes or attributing the full product of \(q-1\)
to the credit would require a separate argument.

The global problem is an upper bound for this specific
determinant-twisted adjoint component of the ray class
group with the two descended local tests. Unrestricted ray
class sizes include many other representations and do not
isolate this component. The classical \(2,3\)-primary degree
also still appears in (DA10).

