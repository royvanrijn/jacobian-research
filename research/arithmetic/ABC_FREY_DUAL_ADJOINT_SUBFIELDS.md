# Smaller division subfields for the dual Frey adjoint bound

Use the Frey hypotheses of the [dual Selmer comparison](ABC_FREY_DUAL_ADJOINT_SELMER.md).
In particular \(N\) is the conductor, \(\ell\ge5\), \(m\ge1\),
\(B_m=\mathbb Z/\ell^m\), and
\[
 F=\mathbb Q(E[\ell^m])^{\{\pm I\}},\qquad
 \Gamma=\operatorname{GL}_2(B_m)/\{\pm I\},\qquad
 M_m'= \operatorname{ad}^0(E[\ell^m])(1).
\]
The group \(D^0_{\ell,m}\) from (DA8) consists of the dual Selmer
classes whose restriction to \(F\) is unramified away from \(\ell\).
For each odd \(q\mid N\), \(q\ne\ell\), keep
\(t_q=v_\ell(v_q(\Delta_E))\), \(d_q=v_\ell(q-1)\) and the local scalar
depth \(s_q\). Write \(Q_\ell\) for this set and \(r_q=(m-t_q)_+\).

For each primitive binary symmetric form \(S\) over \(B_m\) there
is a subfield \(L_S\subset F\) whose ray class characters detect
\(D^0_{\ell,m}\). Primitivity means that at least one matrix entry
of \(S\) is a unit. Put \(h=\min(m,v_\ell(\det S))\), with \(h=m\)
for determinant zero. There is a specified group \(\mathcal X_S\)
of one-dimensional ray class characters of \(L_S\), of modulus
supported only above \(\ell\), such that
\[
 |D^0_{\ell,m}|\ \le\ |\mathcal X_S|\,
              \ell^{\sum_{q\in Q_\ell}\kappa_q(h)},\qquad
 \kappa_q(h)=\big(\min(d_q,r_q)-(r_q-h)_+\big)_+.
 \tag{BD1}
\]
At \(h=0\) the descent preserves all tame factors in (DA1)
and lowers the field degree from order \(\ell^{4m}\) to order
\(\ell^{3m}\). At \(h=m\) it lowers the degree to order
\(\ell^{2m}\). The intermediate depths make the degree and
tame contribution explicit at every coefficient and division
level.

## One-dimensional characters from binary forms

Use the matrix
\(J=\left(\begin{smallmatrix}0&1\\-1&0\end{smallmatrix}\right)\).
Under a fixed Tate-module basis identify \(M_m'\) with trace-zero
matrices, with action \(g\cdot X=(\det g)gXg^{-1}\). The
symmetric matrix \(XJ\) then transforms as \(g(XJ)g^t\).
Define the \(B_m\)-linear functional
\[
             \pi_S(X)=\operatorname{tr}(S XJ).
\]
It satisfies \(\pi_S(g\cdot X)=\pi_{g^t Sg}(X)\).
Let \(P_S\subset\Gamma\) be the similitude group of \(S\) and
\(H_S\) the isometry group, both taken modulo \(\{\pm I\}\).
The multiplier \(\mu_S:P_S\to B_m^\times\) is well-defined and
unique by primitivity. Its kernel is \(H_S\). Put
\[
 L_S=F^{H_S},\qquad K_S=F^{P_S},\qquad
 \Delta_S=\operatorname{Gal}(L_S/K_S)=P_S/H_S.
\]
This is a cyclic extension. Its character \(\mu_S\) is the
multiplier on the quotient.

The restriction of a cocycle \(\xi\) to \(G_{L_S}\), followed by
\(\pi_S\), is a homomorphism \(y_S:G_{L_S}\to B_m\), because
\(H_S\) acts trivially on the functional. Coboundaries map to
zero on this subgroup. Conjugation by a lift of \(P_S\) acts
on \(y_S\) through \(\mu_S\): applying the cocycle identity to
\(gkg^{-1}\), \(k\in G_{L_S}\), cancels the two terms involving
\(\xi(g)\) after projection. This gives
\[
 H^1(\mathbb Q,M_m')\ \longrightarrow\
 \operatorname{Hom}_{\Delta_S}(G_{L_S},B_m(\mu_S)).
 \tag{BD2}
\]
The map is injective. If \(y_S=0\), its restriction to \(G_F\)
vanishes, and equivariance of the homomorphism on \(G_F\)
makes its projections by all \(\pi_{g^t Sg}\) vanish. These
functionals span the dual of \(M_m'\). Indeed, modulo \(\ell\)
the diagonal torus separates the three symmetric coordinates
with characters \(a^2,a,1\) for \(\operatorname{diag}(a,1)\).
The two elementary unipotent subgroups and the coordinate
swap generate all three coordinates from any nonzero
coordinate. Thus every nonzero symmetric form has spanning
orbit modulo \(\ell\). Primitivity, Nakayama's lemma and the
perfect trace pairing on symmetric matrices give the required
span over \(B_m\). Now (DA5) makes the original class zero.
[Full image and central cohomological descent, (DA5)](ABC_FREY_DUAL_ADJOINT_SELMER.md).

## Exact degrees

After basis change and scaling by a unit any such form is
\(S=\operatorname{diag}(1,\eta)\), with \(v_\ell(\eta)=h\) for
\(h<m\) or \(\eta=0\) for \(h=m\). To see this, a primitive
symmetric form in odd characteristic takes a unit value on
some primitive vector. Use this as the first basis vector
and complete the square. Scaling does not alter the groups.
Write \(\varphi_m=\ell^{m-1}(\ell-1)\). The table counts
isometries before quotienting by \(\{\pm I\}\). The last
column is the field degree, with the two factors of \(2\)
from \(\{\pm I\}\) in the Galois group and subgroup cancelling.

| Form | Number of isometries in \(\operatorname{GL}_2(B_m)\) | \([L_S:K_S]\) | \([L_S:\mathbb Q]\) |
|---|---:|---:|---:|
| \(h=0\), split | \(2\varphi_m\) | \(\varphi_m\) | \(\ell^{3m-2}(\ell^2-1)/2\) |
| \(h=0\), nonsplit | \(2\ell^{m-1}(\ell+1)\) | \(\varphi_m\) | \(\ell^{3m-2}(\ell-1)^2/2\) |
| \(0<h<m\) | \(4\ell^{m+h}\) | \(\varphi_m/2\) | \(\ell^{3m-h-3}(\ell-1)^2(\ell+1)/4\) |
| \(h=m\) | \(2\ell^m\varphi_m\) | \(\varphi_m/2\) | \(\ell^{2m-2}(\ell^2-1)/2\) |

For the split form use the hyperbolic coordinates \(xy\);
the isometries are diagonal or antidiagonal with product
of their two nonzero entries equal to one. For the nonsplit
form use the norm of the unramified quadratic algebra over
\(B_m\). Multiplication by norm-one elements and conjugation
are the isometries: after multiplying by the inverse image
of \(1\), an isometry fixes \(1\) and has sign \(+1\) or
\(-1\) on its orthogonal complement. The norm on units is
onto, by the finite-field norm and by square roots on
principal units, so the kernel has order \(\ell^{m-1}(\ell+1)\).
Likewise all multipliers occur for the split form.

For \(0<h<m\), put \(g=\left(\begin{smallmatrix}a&b\\c&d\end{smallmatrix}\right)\).
The condition \(g^t Sg=S\) first gives \(a^2+\eta c^2=1\),
so each \(c\) gives two unit values of \(a\).
The off-diagonal equation fixes \(b=-\eta cd/a\).
The last equation is \(\eta((d/a)^2-1)=0\), which has
\(2\ell^h\) possible \(d/a\), all units. Its determinant
is \(d/a\). If \(\eta=0\), the conditions are \(a=\pm1\),
\(b=0\), \(c\) arbitrary, and \(d\) a unit. In both cases
every multiplier is a square modulo \(\ell\), by the first
diagonal equation of the similitude. Conversely each
unit square is a scalar multiplier. These counts and
\(|\operatorname{GL}_2(B_m)|=\ell^{4m-3}(\ell-1)^2(\ell+1)\)
give the degrees. They use the full finite isometry group,
including lifts whose determinants satisfy the isometry
equations only modulo the relevant power of \(\ell\).
[Binary norm forms and orthogonal group constructions, Bitan–Schein, §4 and Appendix A](https://arxiv.org/abs/1702.05080).

## Exact tame contribution to the descended character

Fix \(q\in Q_\ell\). In a Tate basis the tame inertia
generator \(\tau\) acts by \(T=1+n e_{12}\) and arithmetic
Frobenius acts by
\(F_q=\left(\begin{smallmatrix}q&u\\0&1\end{smallmatrix}\right)\)
up to the common unramified sign. Here \(v_\ell(n)=t_q\).
Choose a prime above \(q\) and let \(S'=\left(\begin{smallmatrix}A&B\\B&C\end{smallmatrix}\right)\)
be the corresponding conjugate form in this basis. Write
\[
 h'=\min(v_\ell(A),v_\ell(B),m),\qquad
 k=\ell^{(m-t_q-h')_+}.
\]
A power \(T^j\) is an isometry of \(S'\) precisely when
\(jnA=jnB=0\) in \(B_m\). The off-diagonal equation gives
the first condition and the final diagonal equation then
gives the second because 2 is invertible. A lift of the
pro-\(\ell\) tame inertia generator over \(L_S\) is thus
\(\tau^k\) in the tame quotient. The Smith factors of \(S'\)
are \(1,\ell^h\), interpreting the second as zero for \(h=m\).
Since the first Tate basis vector is primitive, \(h'\le h\).

Let \(\xi\in D^0_{\ell,m}\) and put \(Y=\xi(\tau)=
\left(\begin{smallmatrix}x&y\\z&-x\end{smallmatrix}\right)\).
The cocycle is unramified after restriction to \(F\), so
(DA6) gives \(\ell^{r_q}Y=0\). Thus its entries lie in
\(\ell^{\min(t_q,m)}B_m\). Its lower-left entry satisfies
\((q-1)z=0\). In fact the relation
\(\operatorname{Frob}_q\tau=\tau^q\operatorname{Frob}_q\)
gives the cocycle identity
\[
 (q\operatorname{Ad}F_q)Y-N_q(\operatorname{Ad}T)Y
      =((\operatorname{Ad}T)^q-1)\xi(\operatorname{Frob}_q),
 \qquad N_q(U)=1+U+\cdots+U^{q-1}.
\]
The lower-left coordinate gives \(z-qz=0\), since
\(\operatorname{Ad}T\) preserves that coordinate and
\(\operatorname{Ad}F_q\) acts on it by \(q^{-1}\).

Project the inertia norm for \(\tau^k\) by \(\pi_{S'}\).
The forms \((T^i)^tS'T^i\) differ from \(S'\) by entries
linear or quadratic in \(i\), with coefficients involving
\(nA,nB,n^2A\). Their sums vanish because
\(\sum_{i<k}i,\sum_{i<k}i^2\) are divisible by \(k\) when
\(k\) is an \(\ell\)-power and \(knA=knB=0\). Consequently
\[
 y_S(\tau^k)=k(-Ay+2Bx+Cz)=kCz .
 \tag{BD3}
\]
The last equality uses \(Y\in\ell^{\min(t_q,m)} M_m'\)
and the defining depth \(h'\). Changes of the local
cocycle representative have no effect after projection on
this inertia subgroup. At all primes over a fixed \(q\)
the values in (BD3) depend on the same \(z\) in the original
Tate basis. They factor through \(z\mapsto
\ell^{(r_q-h)_+}z\), because \(h'\le h\). This single
parameter has order at most
\[
       \ell^{\kappa_q(h)},\qquad
       \kappa_q(h)=\big(\min(d_q,r_q)-(r_q-h)_+\big)_+.
 \tag{BD4}
\]
In particular a unit-determinant form introduces no tame
ramification. The same conclusion holds at a given \(q\)
for every form when \(d_q=0\) or \(r_q=0\).

If \(d_q>0\) and \(r_q>0\), the disappearance of the
projected inertia is equivalent to
\(h\le (r_q-d_q)_+\). Here when \(r_q<d_q\) the condition
forces \(h=0\). Thus the largest depth that gives zero tame loss
for all allowed local cocycles at the odd bad primes is
\[
 h_{\mathrm{free}}=\min\bigl(\{m\}\cup
             \{(r_q-d_q)_+:q\in Q_\ell,\ d_qr_q>0\}\bigr).
 \tag{BD4a}
\]
In particular, for
\(B_\ell=\max(\{0\}\cup\{t_q+d_q:q\in Q_\ell,\ d_q>0\})\)
the choice \(h=(m-B_\ell)_+\) is loss-free. If
\(m\ge B_\ell\), its field degree is \(O(\ell^{2m+B_\ell})\),
retaining all the \(q-1\) factors from (DA1). In particular,
at \(m=e_\ell+\beta_\ell\) the depth
\((e_\ell-\max_{q\in Q_\ell}s_q)_+\) is always loss-free,
with the maximum taken as zero for an empty set. This follows
from \(\beta_\ell\ge t_q+d_q-s_q\) for each \(q\).

The bound (BD4) is exact for the allowed **local** cocycles.
Take \(S'\) with the two diagonal entries exchanged, giving
\(h'=h\) and \(C=1\) after normalizing the first coefficient
of \(S\). Give \(Y\) only a lower-left entry \(z\), with
\(\ell^{r_q}z=(q-1)z=0\). The right-hand operator of the
displayed cocycle identity has image
\(\ell^{\min(t_q,m)}(B_m e_{12}+B_m\operatorname{diag}(1,-1))\):
this follows by applying it to the three trace-zero basis
matrices, since \(q\) and \(2\) are units. The left side
lies in this image by the condition on \(z\) and the
divisibility of \(Y\). Hence a Frobenius value exists.
These values define a continuous cocycle on the tame
semidirect product: the inertia generator maps to an
element of \(\ell\)-power order in the finite semidirect
product with the coefficient module, and the Frobenius
relation has just been checked. Its inertia restriction to
\(F\) vanishes by (DA6). Formula (BD3) then realizes
every value of order dividing \(\ell^{\kappa_q(h)}\) at the
chosen embedding. Thus a local argument from these
hypotheses alone cannot set the descended tame inertia to
zero when \(\kappa_q(h)>0\).
[Tate uniformization and division-point Galois action](https://people.math.harvard.edu/~ctm/home/text/others/tate/tate-unif/tate-unif.pdf).

## Ray class bound and discriminant

In (BD2) the image of \(D^0_{\ell,m}\) is already unramified
at \(2\) and all good primes away from \(\ell\), by the
original conditions over \(\mathbb Q\). The kernel of all
the tame evaluations at \(q\in Q_\ell\) has index at most
\(\prod_q\ell^{\kappa_q(h)}\) by (BD3)–(BD4). Its projected
characters are therefore unramified outside \(\ell\)
over \(L_S\). At \(w\mid\ell\) in \(L_S\), write \(e_w\)
for its absolute ramification index and put
\[
    \mathfrak n_S=\prod_{w\mid\ell}
     \mathfrak p_w^{\, e_w m+\lfloor e_w/(\ell-1)\rfloor+1},
 \qquad
 \mathcal X_S=
     \operatorname{Hom}_{\Delta_S}
       (\operatorname{Cl}_{\mathfrak n_S}(L_S),B_m(\mu_S)).
 \tag{BD5}
\]
The local logarithm bounds the conductor of a character
killed by \(\ell^m\) by the indicated depth, as in (DA9).
The modulus is \(\Delta_S\)-invariant. Global reciprocity
and injectivity of (BD2) prove (BD1).
The ambient group in (BD5) allows all its characters.
The actual image still satisfies the descended local
conditions and the compatibility of its projections under
the full Galois group.

One of these tests has a direct rank-one description. If
\(h=m\) and \(\ell\mid N\), take the embedding whose form in
an \(\ell\)-adic Tate basis is \(\operatorname{diag}(0,1)\).
The local quotient character of \(E[\ell^m]\) has values
\(\pm1\), so the corresponding form is fixed by the entire
local group and determines a prime of \(L_S\) with
completion \(\mathbb Q_\ell\). Projection by this form is
precisely the map to \(U_m^D\) in the multiplicative local
condition after (DA9). Therefore the projected character
is trivial on the whole decomposition group at that prime,
and at its \(\Delta_S\)-conjugates. In this case (BD5) can
omit these primes from the ray modulus and impose that
their ideal classes map to zero.

At \(q\in Q_\ell\), the tame discriminant has the exact
normalized exponent
\[
 \frac{v_q(\operatorname{disc} L_S)}{[L_S:\mathbb Q]}
       =1-\ell^{-r_q}
              \left(1+\frac{\ell-1}{\ell+1}\min(h,r_q)\right).
 \tag{BD6}
\]
Indeed the tame exponent is the field degree minus the
number of inertia orbits on its embeddings. The ratio
of the number of orbits to the degree is the mean of
\(k^{-1}=\ell^{-r_q+\min(h',r_q)}\).
For uniform \(g\in\operatorname{GL}_2(B_m)\) the first
column \(g e_1=(a,c)^t\) is uniform among primitive
vectors. The depth of \((g^t Sg)e_1\) is that of
\(Sg e_1\), namely \(\min(v_\ell(a),h)\) for diagonal
\(S\). Among primitive vectors the proportion with
\(v_\ell(a)\ge j\) is
\(\ell^{1-j}/(\ell+1)\) for \(1\le j\le m\):
the second coordinate must be a unit and the first is
any multiple of \(\ell^j\). The mean of
\(\ell^{\min(v_\ell(a),h,r_q)}\) is therefore
\(1+\min(h,r_q)(\ell-1)/(\ell+1)\) by summing successive
differences. Coset sizes are constant, so the average
over \(g\) is the average over embeddings. This proves
(BD6).
[Jones–Roberts, §4, the tame discriminant as degree minus inertia orbit count](https://msp.org/ant/2014/8-3/ant-v8-n3-p04-p.pdf).
[Other division-field discriminant comparisons](ABC_PROJECTIVE_DIVISION_FIELD_DISCRIMINANTS.md).

## Kummer groups in the Cartan fields for unit determinant

For \(h=0\) there is a further loss-free descent to degree
of order \(\ell^{2m}\). The determinant of the similitude
identity gives the quadratic orientation character
\[
      \psi_S=\det/\mu_S:P_S\longrightarrow\{\pm1\}.
\]
Its kernel is the split Cartan in the hyperbolic case
and the unit group of the unramified quadratic algebra
in the nonsplit case, with \(\{\pm I\}\) removed in both.
In particular this kernel has index two in \(P_S\).
Let \(A_S=F^{\ker\psi_S}\), the quadratic orientation
extension of \(K_S\). The degree formulas above give
\[
 [A_S:\mathbb Q]=
   \begin{cases}
   \ell^{2m-1}(\ell+1)&\text{split},\\
   \ell^{2m-1}(\ell-1)&\text{nonsplit}.
   \end{cases}
 \tag{BD6a}
\]
The multiplier on \(G_{K_S}\) is \(\chi_\ell\psi_S\),
with \(\chi_\ell\) taken modulo \(\ell^m\). Thus
\(L_S A_S=A_S(\mu_{\ell^m})\). The determinant on the
Cartan is surjective, by the product or norm formulas,
so \(A_S\) contains no nontrivial \(\ell\)-power roots
of unity.

Let \(T_S\) be the primes of \(A_S\) over \(\ell\).
Write \(A_S(T_S,\ell^m)\) for classes
\(a\in A_S^\times/A_S^{\times\ell^m}\) whose valuations
are divisible by \(\ell^m\) outside \(T_S\). The superscript
\({}^{-}\) denotes the minus eigenspace of the nontrivial
involution of \(A_S/K_S\). There is an isomorphism
\[
                    \mathcal X_S \simeq A_S(T_S,\ell^m)^- .
 \tag{BD6b}
\]
To prove it, inflation–restriction descends the equivariant
characters to \(H^1(K_S,B_m(\mu_S))\). Both end groups
\(H^1,H^2\) for the cyclic quotient \(\Delta_S\) vanish,
since the multiplier of \(2I\) is \(4\) and \(4-1\) is a
unit. On \(A_S\) the multiplier is cyclotomic, so
restriction through the quadratic orientation field
identifies this group with the minus part of
\(H^1(A_S,\mu_{\ell^m})\). The identification is
Kummer theory. To check ramification, the common
compositum \(L_S A_S\) is unramified over \(A_S\) away
from \(\ell\). Passing from \(L_S\) to this compositum
has degree at most two and detects unramifiedness for
an \(\ell\)-primary character. At primes of \(A_S\)
away from \(\ell\), a Kummer class is unramified exactly
when its valuation is divisible by \(\ell^m\), as
unit roots are obtained over the maximal unramified
field. At primes above \(\ell\), (BD5) already allows
every character of exponent dividing \(\ell^m\).
This proves (BD6b) in both directions.

Consequently the ambient group has the exact sequence
\[
 0\longrightarrow
 \bigl(\mathcal O_{A_S,T_S}^{\times}/
                \mathcal O_{A_S,T_S}^{\times\ell^m}\bigr)^-
 \longrightarrow \mathcal X_S
 \longrightarrow
 \operatorname{Cl}(\mathcal O_{A_S,T_S})[\ell^m]^-
 \longrightarrow 0 .
 \tag{BD6c}
\]
The ordinary Kummer sequence sends \(a\) to the ideal
class whose \(\ell^m\)-th power is \((a)\) in the
\(T_S\)-integer ring. Taking minus eigenspaces is exact
since \(2\) is invertible. In particular the exact
image of \(D^0_{\ell,m}\) can be studied through its
intersection with the minus units and its image in
the minus \(T_S\)-class torsion, both over \(A_S\).
The projected local and compatibility conditions
continue to select that image.
[Cremona–Lingham, §2, the \(m\)-Kummer sequence for \(S\)-integers](https://johncremona.github.io/papers/egros.pdf).

There is a substantial ambient unit contribution. Write
\(D=[K_S:\mathbb Q]\) and \(\varphi_m=\ell^{m-1}(\ell-1)\).
The minus rank of ordinary units is
\((D+r_1(A_S)-r_1(K_S))/2\), by Dirichlet's theorem and
descent of the invariant units. In the split case the
two real-embedding counts are \(2\) and \(1+\varphi_m/2\).
In the nonsplit case they are \(0\) and \(\varphi_m/2\).
Indeed complex conjugation on \(E[\ell^m]\) is conjugate
to \(\operatorname{diag}(1,-1)\), whose centralizer in
the full group has order \(\varphi_m^2\). The number of
fixed cosets of a subgroup is its intersection size
with this conjugacy class times the centralizer order,
divided by the subgroup order. The split Cartan
contains the two diagonal representatives and its
normalizer adds the \(\varphi_m\) antidiagonal
involutions. The nonsplit Cartan contains none, and
its normalizer adds \(\ell^{m-1}(\ell+1)\) norm-one
multiples of conjugation. This gives the stated
counts. The minus unit rank is therefore of order
\(D\) in both cases, and inverting primes above
\(\ell\) can only increase it. A sharp bound through
the exact projected image needs to control how much
of this unit subgroup meets the retained local and
compatibility conditions.

## The uniform estimate over smaller fields

Keep \(m=e_\ell+\beta_\ell\) on the coefficient support
\(\mathcal E(E)=\{\ell\ge5:e_\ell>0\}\). Combining (DA1)
with (BD1) gives, for every selected form \(S\),
\[
 \ell^{e_\ell+\sum_{q\in Q_\ell}(d_q-\kappa_q(h))}
                    \ \le\ |\mathcal X_S|\,
                                     \ell^{4\sum_{q\in Q_\ell}s_q}.
 \tag{BD7}
\]
One may choose the form separately at each \(\ell\).
There are at most two unit square classes of the second
diagonal entry for each \(h<m\) and one rank-one type with zero second
diagonal entry at \(h=m\); scaling and congruence leave the character
group up to isomorphism. In particular an estimate
\[
 (\delta_1)_{\{2,3\}}
   \prod_{\ell\in\mathcal E(E)}
    \min_S\frac{|\mathcal X_S|}
          {\prod_{q\in Q_\ell}\ell^{d_q-\kappa_q(h)}}
                \ \ll_\varepsilon\ \operatorname{rad}(abc)^{2+\varepsilon}
 \tag{BD8}
\]
would imply the normalized dual bound (DA10), and hence
full abc through the signed projector comparison.
An abc failure would force (BD8) to fail. Conversely an
oversize ambient ray class character group would need a
comparison with the subspace coming from the Frey
congruence denominator before it could yield a
counterexample.

The degree savings and the precise tame contribution
hold at every depth. A global bound for the selected
multiplier character at this range of degrees is still
needed. For large \(r_q\), (BD6) also shows why changing
the type of form mainly saves field degree, with a
normalized tame-discriminant difference of order at
most \(r_q/\ell^{r_q}\) at that prime.

