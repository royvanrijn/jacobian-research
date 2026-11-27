# Classical adjoint Selmer groups and the Frey denominator

Let \(E/\mathbb Q\) be the Frey curve of a primitive positive triple
\(a+b=c\) with at least two odd radical primes. Use \(N=N_E\),
\(R=\operatorname{rad}(abc)\), and the signed denominator exponents
\(e_\ell\) from the [strict congruence comparison](ABC_FREY_STRICT_ADJOINT_CONGRUENCES.md).
For \(\ell\ge5\), write \(T=\operatorname{ad}^0(T_\ell E)\),
\(\mathcal V=T\otimes\mathbb Q_\ell\), \(W=\mathcal V/T\), and
\(S_\ell=H^1_f(\mathbb Q,W)\) for the classical torsion
Bloch–Kato Selmer group. Its local condition is the image of
\(H^1_f(\mathbb Q_q,\mathcal V)\), defined using \(B_{\rm cris}\)
at \(\ell\) and unramified cohomology elsewhere. Use the notation
\(M_m=W[\ell^m]\), \(M'_m=M_m^*(1)\) with the finite Tate dual
convention of (DA2), and \(L_{\ell,m}\) from that equation.

At each odd bad prime \(q\ne\ell\) let \(n_q=v_q(\Delta_E)\),
\(t_q=v_\ell(n_q)\), \(d_q=v_\ell(q-1)\), and
\(s_q=\min(t_q,d_q,v_\ell(u_q))\), with \(u_q\) the upper Tate
Frobenius entry. At \(2\) with potentially multiplicative reduction
take \(n_2=-v_2(j(E))\) and \(t_2=v_\ell(n_2)\). Set \(t_2=0\) in
the remaining cases. Sums indexed by \(q\) below are over the odd
bad primes other than \(\ell\). Define
\[
 T_{\ell,m}=\min(m,t_2)+\sum_q\big(\min(m,t_q)+\min(m,s_q)\big),
 \qquad
 T_\ell=t_2+\sum_q(t_q+s_q),\qquad
 \gamma_\ell=\max(\{t_2,0\}\cup\{t_q\}_q).
\]
The sharpest direct comparison proved here is
\[
 e_\ell\ \le\ \operatorname{length}_{\mathbb Z_\ell}
             S_\ell[\ell^{e_\ell}]
             +\min(e_\ell,t_2)+3\sum_q\min(e_\ell,s_q).
 \tag{BK1}
\]
In particular a uniform product bound for the full classical groups at
the primes \(\ell\ge5\) with \(e_\ell>0\) and \(S_\ell[\ell]\ne0\),
together with the existing \(2,3\)-primary degree factor, would
prove abc. The reduction works in both good and multiplicative
coefficient characteristic.

## The coefficient-prime pullback

At all depths and for both reduction types,
\[
 L_{\ell,m}=\operatorname{preimage}_{H^1(\mathbb Q_\ell,M_m)}
                             H^1_f(\mathbb Q_\ell,W).
 \tag{BK2}
\]
First consider any \(\ell\)-adic lattice \(T\) with this definition of
the local condition on \(W\). The stated preimage consists of
reductions of integral classes whose rationalizations lie in
\(H^1_f(\mathbb Q_\ell,\mathcal V)\). Indeed if the image of a
finite class is represented by a rational finite class \(v\), then
\(\ell^m v\) is the image of an integral class by the exact
sequence for \(T\subset\mathcal V\). Its reduction has the same
image in \(H^1(\mathbb Q_\ell,W)\) as the original class. The
difference comes from the boundary of an invariant of \(W\) for
multiplication by \(\ell^m\). If \(Z\in\mathcal V\) lifts that
invariant, this boundary is represented by the reduction of the
integral cocycle \(gZ-Z\), which is rationally a coboundary.

In good reduction, a fixed-determinant Fontaine–Laffaille
deformation over \((\mathbb Z/\ell^m)[\epsilon]/(\epsilon^2)\)
lifts to \(\mathbb Z_\ell[\epsilon]/(\epsilon^2)\) by the adapted
matrix of (DA2) and the trace-zero lifting equation of (RC8).
The rationalized underlying extension of \(V_\ell E\) by itself
is crystalline, which is the \(B_{\rm cris}\) condition for its
adjoint cocycle. Conversely a rationally crystalline extension
from an integral fixed-determinant cocycle yields a lattice in a
crystalline representation with the weights of \(T_\ell E\).
Fontaine–Laffaille theory in these two weights and its
subquotient property give the finite lifts in \(L_{\ell,m}\).
[Fontaine–Laffaille, §§6 and 8](https://www.numdam.org/item/10.24033/asens.1437.pdf).
[The local deformation description, (DA2)](ABC_FREY_DUAL_ADJOINT_SELMER.md).
[The adapted-matrix lift, (RC8)](ABC_FREY_RESIDUAL_ADJOINT_CONTROL.md).

At a multiplicative \(\ell\), the rational adjoint contains the
upper root \(\mathbb Q_\ell(1)\). The quotient has a trivial
diagonal line and an inverse-cyclotomic quotient. Its invariant
line has boundary in \(H^1(\mathbb Q_\ell,\mathbb Q_\ell(1))\)
spanned by the Tate-parameter Kummer class, since conjugating the
trace-zero diagonal introduces twice the Tate extension. The
valuation of the parameter is positive. Thus the one-dimensional
unit subspace of that Kummer group has nonzero image in
\(H^1(\mathbb Q_\ell,\mathcal V)\). The dimension formula for
the Bloch–Kato condition gives dimension one, since the adjoint
has one filtered direction in \(D_{\rm dR}/\mathrm{Fil}^0\)
and no invariants. The latter follows also from the nonzero
boundary of the invariant diagonal. By functoriality the upper
unit subspace therefore maps onto \(H^1_f(\mathbb Q_\ell,\mathcal V)\).
The full upper rational Kummer space has the same image, because
it is spanned by the unit subspace and the Tate-parameter line.
Each finite upper Kummer class lifts integrally and maps into
that rational image after division by \(\ell^m\). This proves
the inclusion of \(L_{\ell,m}\) in the asserted preimage.

The image of the rational finite line in \(H^1(\mathbb Q_\ell,W)\)
is divisible of corank one. Its kernel is the intersection of
the line with the image of \(H^1(\mathbb Q_\ell,T)\), a full
lattice in rational cohomology. Also \(H^0(\mathbb Q_\ell,W)\)
is finite: it is the torsion of \(H^1(\mathbb Q_\ell,T)\), since
the rational invariants vanish. The multiplication sequence on
\(W\) now gives preimage order
\(\ell^m|H^0(\mathbb Q_\ell,M_m)|\), the order of \(L_{\ell,m}\)
by (DA2). This proves (BK2) at a multiplicative prime as well.
[Bloch–Kato, Corollary 3.8.4 and Example 3.9](https://virtualmath1.stanford.edu/~conrad/BSDseminar/refs/BKTamagawa.pdf).

## Unramified conditions away from the coefficient prime

Write \(C_q=H^1_{\rm unr}(\mathbb Q_q,W)/H^1_f(\mathbb Q_q,W)\)
for \(q\ne\ell\). At a good prime or a potentially good \(2\),
the inertia image on \(T\) has order prime to \(\ell\), so
\(\mathcal V^{I_q}\to W^{I_q}\) is surjective by averaging.
Taking Frobenius coinvariants gives \(C_q=0\), and inertia
invariants of \(T\) surject onto those of \(M_m\).

At a potentially multiplicative prime of residue characteristic
\(q\ne\ell\), work with the untwisted adjoint in the Tate basis.
The quadratic twist to a Tate curve disappears from its action.
Use \(n=n_q\), \(t=v_\ell(n)\), \(d=v_\ell(q-1)\), \(u=u_q\)
and \(s=\min(t,d,v_\ell(u))\), with \(d=s=0\) when \(q=2\).
For matrices \(\left(\begin{smallmatrix}x&y\\z&-x\end{smallmatrix}\right)\)
with entries in \(D=\mathbb Q_\ell/\mathbb Z_\ell\), the inertia
invariants are given by \(nx=nz=0\), with \(y\) unrestricted.
The same description holds with integral or finite coefficients.
The rational inertia invariants consist only of the upper root,
where Frobenius acts by \(q\), so the rational unramified group
is zero. In the divisible invariants the upper root has
surjective \(F-1\). On the quotient \(D[n]^2\) of the invariant
module, \(F-1\) sends \((x,z)\) to
\((uz/q,(1/q-1)z)\). Its image is cyclic of order \(\ell^{t-s}\)
and, after extracting \(\ell^s\) if \(s<t\), the generating
direction is primitive. Hence
\[
 C_q\simeq\mathbb Z/\ell^t\oplus\mathbb Z/\ell^s.
 \tag{BK3}
\]
The same Frobenius calculation shows that the quotient of \(W^{G_q}\)
by the invariants of the divisible upper root is killed by \(\ell^t\):
it embeds in \(D[n]^2\). The local index (BK3) is the untwisted
case of the Tate inertia computation in (CT2).
[Tate inertia and the Frobenius action on cohomology](ABC_CUBIC_ADJOINT_TAMAGAWA_FACTORS.md).

## Exact order of the dual tower and control at finite depth

At \(q\ne\ell\) the local Tate annihilator of the condition (RC7)
is
\[
 A_{q,m}=\operatorname{im}\big(H^1_{\rm unr}(\mathbb Q_q,T)
                                      \to H^1(\mathbb Q_q,M_m)\big).
 \tag{BK4}
\]
For the inclusion, pair an integral unramified cocycle with the
image in \(W'\) of a class satisfying (RC7). Both classes admit
unramified cocycle representatives, and their cup product uses
the procyclic unramified quotient of \(\ell\)-cohomological
dimension one. This pairing agrees with the finite pairing after
reducing the integral cocycle. At good primes and potentially
good \(2\), (RC7) is the finite unramified condition and
\(A_{q,m}\) is the full finite unramified condition. These are
exact annihilators.

In the Tate case set \(\mu=\min(m,t)\). The full finite primal
unramified condition has order
\(\ell^{2\min(m,s)+\min(m,t+d-s)}\), the order of the local
invariants in (RC3). The finite dual quotient of (RC7) by the
unramified condition has order \(\ell^{\mu+\min(m,s)}\), by the
local calculation in (RC9a). Thus the annihilator's length is
\[
                 \min(m,s)+\min(m,t+d-s)-\mu .
\]
To count \(A_{q,m}\), reduce the integral upper root in the
Frobenius coinvariants of the finite inertia invariants. For
\((F-1)\) of an invariant to have vanishing lower and diagonal
coordinates, one needs \(uz=(q-1)z=0\). The upper coordinates of
all such images generate \((q-1)(\mathbb Z/\ell^m)\)
plus \(2u(\mathbb Z/\ell^m)[\ell^\mu]\), since \(uz=0\).
This gives \(A_{q,m}\) length
\(\min(m,d,m-\mu+v_\ell(u))\). It equals the displayed length:
for \(m\le t\) both are \(\min(m,s)\). For \(m>t\), both are
\(\min(m,d)\) if \(s=t\), \(d\) if \(s=d<t\), and
\(\min(d,m-t+s)\) otherwise. This also covers \(q=2\) with
\(d=s=0\), including \(u=0\) with valuation \(+\infty\).
This proves equality in (BK4).
[Local invariants and the saturated dual inertia quotient, (RC3) and (RC9a)](ABC_FREY_RESIDUAL_ADJOINT_CONTROL.md).

Let \(A_{\ell,m}\) be the global primal Selmer group with
conditions (BK4) and \(L_{\ell,m}\). The globals have no
invariants by the full Frey image. The coefficient and real
factors in Greenberg–Wiles cancel by (DA2). At \(q\ne\ell\),
the ratio of the size of the full finite unramified condition
to that of \(A_{q,m}\) is the local factor just computed.
The dual of \(A_{\ell,m}\) is \(V_{\ell,m}\) from (RC7).
Consequently
\[
               |V_{\ell,m}|=\ell^{T_{\ell,m}}|A_{\ell,m}|.
 \tag{BK5}
\]
The analogous primal group with the full unramified conditions
away from \(\ell\) has order \(|U_{\ell,m}|\) by the same
formula.
[Finite duality, Wiles, Proposition 1.6](https://jontalle.web.engr.illinois.edu/Public/Wiles-Fermat.95.pdf).

The global maps into \(H^1(\mathbb Q,W)\) are injective and (BK2)
gives the exact coefficient condition there. The groups
\(A_{\ell,m}\) are nested under inclusion of the coefficients
and map into \(S_\ell[\ell^m]\), since the integral unramified
classes rationalize into the unramified condition. Conversely,
lift a class of \(S_\ell[\ell^j]\) globally to \(M_j\).
At any potentially multiplicative \(q\ne\ell\), its localization
is in the kernel of the map to cohomology of \(W\) by (BK3) and
the vanishing of rational unramified cohomology there. It is
therefore the boundary of \(b\in W^{G_q}\). After coefficient
inclusion at depth \(j+h\), it becomes the boundary of
\(\ell^h b\) at that depth. For \(h\ge t\), this element is
in the invariant divisible upper root. Its boundary is the
reduction of an integral unramified cocycle, obtained by lifting
it to the rational upper root and taking the coboundary.
At good primes away from \(\ell\) and potentially good \(2\),
the finite map on inertia cohomology is injective, by
trivial action or prime-to-\(\ell\) averaging on the inertia
image. Those localizations already belong to (BK4).
It follows from (BK5) that, for every \(m\ge1\),
\[
 |S_\ell[\ell^{(m-\gamma_\ell)_+}]|
       \ \le\ \ell^{-T_{\ell,m}}|V_{\ell,m}|
       \ \le\ |S_\ell[\ell^m]|.
 \tag{BK6}
\]
Here the zeroth torsion layer denotes the zero group.
In particular if \(\gamma_\ell=0\), all truncated orders are
exactly equal and the two finite limits have the same elementary
divisors.

There is also a comparison at the same depth in both directions.
At a potentially multiplicative \(q\ne\ell\), the pullback of
the classical condition on \(W\) is the kernel of the map to
\(H^1(\mathbb Q_q,W)\). Its order is \(|H^0(\mathbb Q_q,M_m)|\)
by the boundary calculation of (RC3). The same calculation
holds at \(2\), where the Tate matrix gives invariant order
\(\ell^{\min(m,t_2)}\). The index of \(A_{q,m}\) in that
pullback therefore equals its index in the finite unramified
condition. At the other primes these two conditions agree.
Thus \([S_\ell[\ell^m]:A_{\ell,m}]\le\ell^{T_{\ell,m}}\) by
localization. Equations (BK5) and (RC9a) give
\[
 \begin{split}
 |S_\ell[\ell^m]| &\le |V_{\ell,m}|
                 \le \ell^{T_{\ell,m}}|S_\ell[\ell^m]|,\\
 \ell^{-T_{\ell,m}}|S_\ell[\ell^m]| &\le |U_{\ell,m}|
                 \le \ell^{T_{\ell,m}}|S_\ell[\ell^m]|.
 \end{split}
 \tag{BK6a}
\]
In particular \(S_\ell[\ell]\ne0\) forces \(V_{\ell,1}\ne0\).

Let \(S_\ell^{\rm unr}\) use unramified cohomology in \(W\)
everywhere away from \(\ell\) and the classical condition there.
It is the union of the finite primal groups with the full
unramified conditions away from \(\ell\). Indeed at the finitely
many bad primes a finite localization that becomes unramified in
divisible coefficients becomes unramified already in some
finite layer. This uses the direct limit in continuous
cohomology on inertia. At good primes away from \(\ell\) the
finite maps on inertia cohomology are injective. The same argument
on the dual side using (RC7) and (RC8) gives
\(U_{\ell,\infty}=V_{\ell,\infty}\). The finiteness of this
tower is proved after (RC9a). The limit of \(A_{\ell,m}\) is
\(S_\ell\) by (BK6). Taking orders and using (BK5) proves
\[
 |U_{\ell,\infty}|=|V_{\ell,\infty}|
            =|S_\ell^{\rm unr}|=|S_\ell|\ell^{T_\ell},
 \qquad
 S_\ell^{\rm unr}/S_\ell\xrightarrow{\ \sim\ }
                  \bigoplus_{q\ne\ell} C_q .
 \tag{BK7}
\]
The last map is localization. It is injective by the definitions
and surjective by (BK3) and the order identity.
The equality between primal and dual orders uses duality and
asserts only a cardinality comparison.
[Pointwise finiteness and the coefficient inclusions, (RC2) and (RC9a)](ABC_FREY_RESIDUAL_ADJOINT_CONTROL.md).

## Signed classes and the uniform bound

At every \(q\ne\ell\), \(C_\ell^{\rm str}\) in the notation of
(SA1) is unramified in \(W\), including at \(2\). It is zero on
the full decomposition group at each odd bad \(q\ne\ell\).
At \(\ell\) it satisfies (BK2) by (DA2), so it injects into
\(S_\ell^{\rm unr}[\ell^{e_\ell}]\). The only possible
localization in (BK7) is at potentially multiplicative \(2\),
where the factor is cyclic of order \(\ell^{t_2}\). The kernel
of this localization lies in \(S_\ell[\ell^{e_\ell}]\).

At each odd bad \(q\ne\ell\), the quotient of the solution
space in (SA3) by the stated basis changes is killed by
\(\ell^{s_q}\) and has at most three cyclic factors. In addition
to the lower inertia coordinate of depth \(s_q\), the
Frobenius pair gives two factors of that depth. To see this,
divide the row \((q-1,u_q)\) by its common content and use a
basis with that primitive row as the first coordinate
functional. The vector \((-u_q,q-1)\), divided by the same
content, is a primitive vector on the second coordinate axis.
On \(D[n_q]^2\), the kernel condition and then quotient by
this vector give one factor of order \(\ell^{s_q}\) each.
If the row is zero modulo \(\ell^{t_q}\), both factors already
have depth \(t_q=s_q\). Since \(C_\ell\) itself is killed by
\(\ell^{e_\ell}\), the image at \(q\) has order at most
\(\ell^{3\min(e_\ell,s_q)}\). Together with the dyadic bound
and \(\ell^{e_\ell}\le |C_\ell|\), this proves (BK1).
[Signed conormal injection, annihilator and local matrix equations, (SA1)–(SA3)](ABC_FREY_STRICT_ADJOINT_CONGRUENCES.md).

Put \(n_2=1\) when \(2\) is not potentially multiplicative.
The product over all \(\ell\ge5\) of the local factors in (BK1)
divides \(n_2\prod_{q\mid N,\ q\ {\rm odd}}\gcd(n_q,q-1)^3\).
Its logarithm is \(O(\log c/\log(\log c+2))\) by (SA5).
In particular the **classical residual zero test** \(S_\ell[\ell]=0\)
eliminates these coefficient primes up to this one product.
It applies even if \(V_{\ell,1}\ne0\) due to tame contributions.
The remaining sufficient bound is
\[
 (\delta_1)_{\{2,3\}}
    \prod_{\substack{\ell\ge5:\ e_\ell>0\\ S_\ell[\ell]\ne0}}
        |S_\ell|
        \ \ll_\varepsilon\ R^{2+\varepsilon}
        \qquad(\varepsilon>0),
 \tag{BK8}
\]
where \(\delta_1\) is the classical optimal modular degree.
One can replace each \(|S_\ell|\) in (BK8) by
\(|S_\ell[\ell^{e_\ell}]|\). The deduction uses (BK1) and the
signed projector comparison in (SA7), with the same elementary
low-support cases. By (BK6a), the truncated bound and the
residual-support bound after (RC11) are equivalent. Indeed on
the common support the products differ by at most
\(\exp(O(\log c/\log(\log c+2)))\), and the extra dual
coefficients with \(S_\ell[\ell]=0\) cost at most the same
factor. In either direction, first use the abc consequence of
the assumed bound to get \(\log c=O(\log R)\). The comparison
factor can then be absorbed in \(R^\varepsilon\) by using the
assumed estimate with a smaller \(\varepsilon\).
The exact integral degree comparison at all coefficient primes
\(\ell\ge5\) gives the converse. Either the full or the
truncated product in (BK8) gives a sharp-bound criterion
equivalent to abc.
[Sharp signed degree implication, (SA7)](ABC_FREY_STRICT_ADJOINT_CONGRUENCES.md).
[Full degree calibration and the uncapped criterion](ABC_FREY_MULTIPLICATIVE_ADJOINT_DEGREE.md).

## The global length obstruction

Diamond–Flach–Guo, Theorem 8.10, supplies the adjoint special
value formula at good coefficient primes here. Its weight-two
hypotheses exclude primes dividing \(2N\) and require residual
absolute irreducibility over the quadratic subfield of
\(\mathbb Q(\mu_\ell)\). For \(\ell\ge5\) with \(\ell\nmid N\)
the latter holds by the full Frey image. Their equation (58)
uses the classical Selmer Fitting ideal after identifying the
Tate–Shafarevich module with the torsion Selmer group. Its
prime-to-\(\ell\) Tamagawa factors are the inertia torsion
invariants computed locally here: the long exact inertia
sequence identifies that torsion with
\(W^{I_q}/\operatorname{im}\mathcal V^{I_q}\), and the Tate
coordinates above give order \(\ell^{t_q+s_q}\) on its
Frobenius invariants. At good coefficient primes, the integral
pairing gives the exact length \(v_\ell(\delta_1)-T_\ell\).
The ordinary comparison at coefficient primes dividing \(N\)
gives a separate local correction there.
[Diamond–Flach–Guo, Theorem 0.1, equations (57)–(58) and Theorem 8.10](https://arxiv.org/abs/2512.02348v1).
[Integral degree calibration and the multiplicative coefficient comparison](ABC_FREY_ADJOINT_DEGREE_CALIBRATION.md).

For a globally semistable curve with good coefficient reduction
and the residual-minimal hypotheses of Kim, Corollary 1.4
identifies the classical adjoint length with the valuation of
the modular degree. The singular localization of the Flach
class at an auxiliary Flach prime is that same degree times
a coefficient-prime unit by §4.2. The explicit localization
therefore transmits degree bounds in that case and gives no
independent height estimate. Corollary 1.5 adds a cyclicity
criterion under an indivisibility hypothesis. A cyclic group
can still have arbitrarily large length, so its order remains
the quantity to control. These results have separate
coefficient and residual hypotheses from (BK1)–(BK8).
[Kim, §3, §4.2 and Corollaries 1.4–1.5](https://arxiv.org/abs/2403.01500v7).

The next quantitative step is to bound the order of the
classical adjoint groups on the surviving support by a
Frey-specific global argument, with a joint budget for the
\(2,3\)-primary factor. The computation of every tame
quotient in (BK7) isolates the part already explained by
Tate valuations. A method based on Euler or Stark classes
needs a singular invariant with an independently bounded
coefficient-prime valuation to improve the abc bound.
