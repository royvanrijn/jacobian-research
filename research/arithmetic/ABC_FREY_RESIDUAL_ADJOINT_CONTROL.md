# Residual control of the Frey adjoint denominator

Keep the Frey notation of the [finite duality comparison](ABC_FREY_DUAL_ADJOINT_SELMER.md).
For \(\ell\ge5\) the signed Frey projector denominator is
\(\ell^{e_\ell}\). Write \(R=\operatorname{rad}(abc)\),
\(M_j=\operatorname{ad}^0(E[\ell^j])\), and
\(M'_j=M_j(1)\). Put
\[
 T=\operatorname{ad}^0(T_\ell E),\qquad T'=T(1),\qquad
 W=T\otimes\mathbb Q_\ell/\mathbb Z_\ell,\qquad
 W'=T'\otimes\mathbb Q_\ell/\mathbb Z_\ell .
\]
For odd \(q\mid N\), \(q\ne\ell\), let \(t_q,d_q,s_q\)
be as in (DA1), and write \(Q_\ell\) for this set of primes.
Put \(t_2=v_\ell(-v_2(j(E)))\) when \(E\) is potentially
multiplicative at \(2\), and \(t_2=0\) otherwise.

There is a finite dual Selmer group \(U_{\ell,j}\) with the
condition \(L_{\ell,j}^{\perp}\) of (DA2) at \(\ell\), and
unramified conditions at every finite prime \(q\ne\ell\).
The depth needed to bound the denominator is now
\(j=e_\ell\), with no additional division depth:
\[
 \ell^{e_\ell}\ \le\
 |U_{\ell,e_\ell}|\,\ell^{\,\min(e_\ell,t_2)+
     \sum_{q\in Q_\ell}(3s_q+\lambda_{q,\ell,e_\ell})},
 \qquad
 \lambda_{q,\ell,j}
   =\min(j,t_q)+\big(\min(j,t_q)-t_q+s_q\big)_+ .
 \tag{RC1}
\]
The total logarithm of the local factor is
\(O(\log c/\log(\log c+2))\) across all coefficient primes.
The earlier dual bound (DA1) includes the \(q-1\) divisors.
Here the local comparison absorbs those factors to use the
smaller depth and unramified conditions over \(\mathbb Q\)
at every prime other than \(\ell\).

A second group \(V_{\ell,j}\) contains \(U_{\ell,j}\).
It uses the same coefficient-prime condition, and at
\(q\ne\ell\) allows the classes which become unramified
after mapping \(M'_j\) to \(W'\). The groups at all
depths satisfy
\[
              V_{\ell,j}\simeq V_{\ell,m}[\ell^j]
                            \qquad(1\le j\le m).
 \tag{RC2}
\]
Thus the residual group \(V_{\ell,1}\) controls the number
of generators at every depth. It is itself unramified
away from \(\ell\) and the bad primes for which the
positive Tate valuation is divisible by \(\ell\).
This gives a uniform elimination of coefficient primes
with zero residual group.

## Finite classes trivial in divisible coefficients

Work first at any finite depth \(j\ge1\). At each
\(q\in Q_\ell\), in place of the zero condition of
(DA3), use the kernel
\[
 B_{q,j}=\ker\big(H^1(\mathbb Q_q,M_j)
                              \longrightarrow H^1(\mathbb Q_q,W)\big).
\]
The exact sequence for multiplication by \(\ell^j\)
on \(W\) identifies this group with
\(H^0(\mathbb Q_q,W)/\ell^j H^0(\mathbb Q_q,W)\).
The invariant group here is finite. In the Tate
basis its three invariant divisors have depths
\(s_q,s_q,t_q+d_q-s_q\), by (DA3). In particular
\[
                    |B_{q,j}|=|H^0(\mathbb Q_q,M_j)|.
 \tag{RC3}
\]
The local Tate annihilator of \(B_{q,j}\) inside
\(H^1(\mathbb Q_q,M'_j)\) is
\[
 I_{q,j}=\operatorname{im}\big(H^1(\mathbb Q_q,T')
                                 \longrightarrow H^1(\mathbb Q_q,M'_j)\big).
 \tag{RC4}
\]
Indeed the sequence
\(0\to T'\xrightarrow{\ell^j}T'\to M'_j\to0\)
is dual to that for \(W\). The connecting map
\(H^1(\mathbb Q_q,M'_j)\to H^2(\mathbb Q_q,T')\)
is adjoint to the boundary from \(H^0(\mathbb Q_q,W)\)
which defines \(B_{q,j}\). This follows from the cup
product boundary identity and local Tate duality
between the compact lattice and the divisible module.
The degree-two pairing for the lattice is the inverse
limit of finite local Tate pairings with the degree-zero
invariants, whose direct limit is the invariant group of
the divisible module. Its kernel is the image in (RC4).
[Mazur–Rubin, §§1–2, local Tate duality for a lattice and its dual](https://jtnb.centre-mersenne.org/item/10.5802/jtnb.933.pdf).

Define \(S^{\rm sat}_{\ell,j}\) using \(B_{q,j}\) at
\(q\in Q_\ell\), \(L_{\ell,j}\) at \(\ell\), and the
unramified condition elsewhere, except that at \(2\)
allow all local classes if \(E\) has potentially
multiplicative reduction. Denote the dual Selmer
group by \(D^{\rm int}_{\ell,j}\). It uses \(I_{q,j}\)
at the odd bad primes and the zero condition at
\(2\) in the potentially multiplicative case.

The global invariant groups vanish as in (DA4).
The coefficient prime and real-place factors again
cancel, and (RC3) cancels the remaining odd factors.
At \(2\) the unrestricted primal condition contributes
\(|H^0(\mathbb Q_2,M'_j)|\) only in the potentially
multiplicative case. This order is
\(\ell^{\min(j,t_2)}\): after the splitting twist
Frobenius can be diagonalized with ratio \(2\),
since \(2-1\) is a unit. The twisted adjoint Frobenius
weights are \(4,2,1\). Only the lower root is invariant,
and commuting with Tate inertia kills this entry by
\(v_2(t_E)=-v_2(j(E))\). Thus
\[
             |S^{\rm sat}_{\ell,j}|
                  =|D^{\rm int}_{\ell,j}|\ell^{\min(j,t_2)} .
 \tag{RC5}
\]
[Wiles, Proposition 1.6, the finite Greenberg–Wiles ratio](https://jontalle.web.engr.illinois.edu/Public/Wiles-Fermat.95.pdf).

The strict signed conormal classes killed by \(\ell^j\)
inject into \(S^{\rm sat}_{\ell,j}\). Their natural
derivation cocycles have coefficients in \(M_j\).
By strictness their local restrictions at odd bad
primes are in \(B_{q,j}\). At \(2\) with potentially
good reduction their inertia cocycle is already a
coboundary at depth \(j\), by the prime-to-\(\ell\)
averaging argument used before (DA4). Potentially
multiplicative reduction uses the full condition at
\(2\), so imposes no finite-depth lift constraint.
The global class injection follows from
\(H^0(\mathbb Q,W)=0\). The local coefficient-prime
compatibility is (DA2).
[Signed conormal classes, strictness and the dyadic argument](ABC_FREY_STRICT_ADJOINT_CONGRUENCES.md).

## Cutting the integral dual condition to unramified classes

At an odd bad \(q\ne\ell\), put \(n=v_q(\Delta_E)\) and
write the Tate matrices
\[
       T_q=1+n e_{12},\qquad
       F_q=\begin{pmatrix}q&u\\0&1\end{pmatrix}
\]
after the unramified sign twist. On the pro-\(\ell\)
tame inertia generator the cyclotomic character is
trivial. Thus \(H^1(I_q,T')\) is the coinvariant
module for \(\operatorname{Ad}(T_q)-1\). Its diagonal
and upper entries have coefficient ring
\(\mathbb Z_\ell/\ell^{t_q}\), with lower entry free
over \(\mathbb Z_\ell\). Indeed the image of
\(\operatorname{Ad}(T_q)-1\) is \(n\) times the
span of the diagonal and upper vectors. The factors
of \(2\) in these matrices are units. The wild
and prime-to-\(\ell\) inertia kernels act trivially
and have no \(\ell\)-primary \(H^1\).

Frobenius acts on the inertia coinvariants by
\(\operatorname{Ad}(F_q)\). Its determinant twist is
cancelled when changing the tame generator by
the inverse \(q\)-power. On the lower coordinate
the resulting multiplier is \(q^{-1}\); it has no
nonzero invariant in the free \(\mathbb Z_\ell\)
coordinate. On the remaining two coordinates \((x,y)\),
the invariance condition is
\[
                 (q-1)y-2u x=0\bmod\ell^{t_q}.
\]
The invariants of \(H^1(I_q,T')\) are reached by
restriction from \(H^1(\mathbb Q_q,T')\), since the
unramified quotient has \(\ell\)-cohomological
dimension one. For compact coefficients this follows
from the finite coefficient assertion by inverse
limit, using finiteness of the degree-one terms.

After an invertible linear change over \(\mathbb Z_\ell\),
the last equation is multiplication by \(\ell^{s_q}\)
on one coordinate over \(\mathbb Z_\ell/\ell^{t_q}\).
The image on reduction to
\(\mathbb Z_\ell/\ell^{\min(j,t_q)}\) therefore has
order exactly \(\ell^{\lambda_{q,\ell,j}}\) as in (RC1).
This is the local image on inertia of the group in
(RC4). In particular
\(\lambda_{q,\ell,j}\le t_q+s_q\), and the locals
are already unramified if \(t_q=0\).
[Tate uniformization and the division-point action](https://people.math.harvard.edu/~ctm/home/text/others/tate/tate-unif/tate-unif.pdf).

Requiring \(D^{\rm int}_{\ell,j}\) to be unramified
at these primes costs an index at most
\(\prod_{q\in Q_\ell}\ell^{\lambda_{q,\ell,j}}\).
The resulting subgroup lies in \(U_{\ell,j}\), since
the condition at \(2\) was unramified or zero.
Use (RC5) at \(j=e_\ell\), where the entire conormal
group is killed by \(\ell^{e_\ell}\), and apply (SA1).
This proves (RC1). The proof only cuts rational
decomposition-group conditions and uses no bound
on the degree of a number field.

The product of these local factors is bounded by
\[
       \max(1,-v_2(j(E)))\!
                   \prod_{q\mid N,\ q\ {\rm odd}}
                                v_q(\Delta_E)^5 ,
 \tag{RC6}
\]
where the first factor is one in the potentially
good case at \(2\). This follows from \(s_q\le t_q\)
and \(\lambda_{q,\ell,j}\le t_q+s_q\).
The first factor is \(O(1+\log c)\) and the odd product
has the same divisor bound as in (SA5).
[Uniform divisor-height estimate, (SA5)](ABC_FREY_STRICT_ADJOINT_CONGRUENCES.md).

## Exact coefficient control and the residual conditions

For \(q\ne\ell\) define the local condition of \(V_{\ell,j}\)
to be the kernel of
\[
        H^1(\mathbb Q_q,M'_j)\longrightarrow H^1(I_q,W').
 \tag{RC7}
\]
The coefficient-prime condition is \(L_{\ell,j}^{\perp}\).
The compatibility of (RC7) with coefficient inclusions
is immediate from using the same target for every \(j\).
At \(\ell\) the reduction \(L_{\ell,m}\to L_{\ell,j}\)
is surjective for \(m\ge j\). In good reduction one
uses the adapted Fontaine–Laffaille matrix of (DA2);
the determinant-linear equation is the space of
trace-zero matrices after multiplication by the
inverse baseline matrix, and any such matrix lifts.
In the multiplicative case the upper-root Kummer
groups reduce surjectively, so their images do too.
Local Tate duality now gives
\[
 \iota_{j,m}^{-1}\big(L_{\ell,m}^{\perp}\big)
                               =L_{\ell,j}^{\perp},
 \qquad M'_j\xrightarrow{\iota_{j,m}} M'_m,
 \tag{RC8}
\]
where the preimage is on cohomology and the coefficient
map is multiplication by \(\ell^{m-j}\) in integral
coordinates. The global invariant groups \(H^0(\mathbb Q,M'_k)\)
vanish for all \(k\) by the scalar \(2I\). The
coefficient exact sequences consequently identify
\(H^1(\mathbb Q,M'_j)\) with the \(\ell^j\)-torsion
in \(H^1(\mathbb Q,M'_m)\). Pulling back all local
conditions proves (RC2).
[Mazur–Rubin, Proposition 3.3, control by compatible Selmer conditions](https://jtnb.centre-mersenne.org/item/10.5802/jtnb.933.pdf).

At an odd bad \(q\ne\ell\), the residual condition
(RC7) has the following explicit description.

- If \(t_q=0\), it is the unramified condition.
- If \(t_q>0\), a residual class may have on the
  pro-\(\ell\) tame generator a diagonal and upper
  entry \((x,y)\) satisfying
  \((q-1)y-2u x=0\bmod\ell\), with zero lower entry.
  Thus the extra ramification directions have
  dimension one if \(s_q=0\) and two if \(s_q>0\).

Indeed in \(W'\) the diagonal and upper directions
are already in the image of
\(\operatorname{Ad}(T_q)-1\), by divisibility. The
remaining lower direction detects restriction on
inertia. When \(t_q=0\), the first-layer inertia
coinvariants are already precisely this lower
direction. When \(t_q>0\), inertia is trivial on
the residual module, and (RC7) requires the lower
coordinate to vanish. Frobenius invariance of the
other two coordinates gives the stated relation.
All such pairs lift to a local cohomology class by
the unramified quotient's cohomological dimension.

The same calculation at \(2\) in the potentially
multiplicative case gives an extra one-dimensional
ramification line only if \(t_2>0\), since \(2-1\)
is a unit. At potentially good primes of residue
\(2\) the Frey inertia image has order prime to
\(\ell\). Restrict to its kernel and use
inflation–restriction and averaging: the map on
inertia cohomology induced by \(M'_1\to W'\) is
injective. The same is clear at good primes of
residue different from \(\ell\). Thus these cases
give only the unramified condition.

It follows that the only additional finite primes
at which the residual classes can ramify are
\[
 \{\,q\in Q_\ell:\ell\mid v_q(\Delta_E)\,\}
      \ \cup\ \{\,2:t_2>0\text{ and }E
                          \text{ is potentially multiplicative at }2\,\}.
 \tag{RC9}
\]
Here all tests at the stated support and the original
condition at \(\ell\) use the mod-\(\ell\) coefficient
group. Equations (OK2)–(OK4) describe the resulting
residual group in the minus Kummer group of a nonsplit
Cartan field of degree \(\ell(\ell-1)\). One needs to
allow valuations only at the places above (RC9) and,
in the exceptional multiplicative case, at \(\ell\)
when its Tate parameter is an \(\ell\)-th power.
The latter assertion uses (OK5)–(OK6), which depend
only on the dual condition at the coefficient prime.
The exact local tests on the recovered global class
are still imposed at the coefficient prime and at
its bad primes. In particular base unramifiedness at
a Tate prime with \(t_q=0\) is tested on the recovered
class over \(\mathbb Q_q\).
[The projective Kummer test](ABC_FREY_CARTAN_ORDINARY_KUMMER.md).

## Uniform comparison with the unramified group

The difference between \(V_{\ell,j}\) and \(U_{\ell,j}\)
is itself bounded independently of \(j\). At an odd bad
prime the inertia coinvariants at depth \(j\) have lower
coordinate free over \(B_j\) and two other coordinates
over \(\mathbb Z/\ell^\mu\), where \(\mu=\min(j,t_q)\).
The kernel of restriction to the divisible inertia
coinvariants is exactly the latter two coordinates.
Their Frobenius equation is
\((q-1)y-2ux=0\bmod\ell^\mu\), and its kernel has order
\(\ell^{\mu+\min(j,s_q)}\). At a potentially multiplicative
\(2\) the same computation gives \(\ell^{\min(j,t_2)}\),
and the remaining primes contribute one. Consequently
\[
 [V_{\ell,j}:U_{\ell,j}]
 \ \le\ \ell^{\,\min(j,t_2)+
        \sum_{q\in Q_\ell}(\min(j,t_q)+\min(j,s_q))}.
 \tag{RC9a}
\]
Over all coefficient primes this factor is again
\(\exp(O(\log c/\log(\log c+2)))\).

For a fixed Frey curve and \(\ell\), the limit
\(V_{\ell,\infty}\) is finite. To verify this, the dual
of \(U_{\ell,j}\) is the primal group defined by
\(L_{\ell,j}\) at \(\ell\) and the unramified conditions
at all other primes. These two finite groups have
equal order by Greenberg–Wiles, because every finite
factor away from \(\ell\) is one and the real and
coefficient factors cancel. The primal groups inject
at every depth into the finite torsion Selmer group
of Khare–Patrikis. At the coefficient prime in good
reduction, the adapted Fontaine–Laffaille matrix lifts
over \(\mathbb Z_\ell[\epsilon]/(\epsilon^2)\) by the
trace-zero description in (RC8), so its underlying
torsion representation is a quotient of a crystalline
lattice. In the multiplicative case, lift the upper
root class by Kummer theory to
\(H^1(\mathbb Q_\ell,\mathbb Z_\ell(1))\).
Over \(\mathbb Q_\ell\) the resulting square-zero
representation, without its common unramified sign,
is an extension of two trivial lines by two
cyclotomic lines. Every such local extension here
is semistable. The Kummer space for one trivial
and one cyclotomic line is spanned by parameters
of positive valuation: for a unit, take the
difference of the classes of that unit times
\(\ell\) and of \(\ell\). Positive-valuation classes
come from Tate curves. Direct sums, pullbacks and
pushouts give the extension with two lines on each
side. This proves the torsion semistable coefficient
condition used by Khare–Patrikis, with the weights
of \(T_\ell E\). At other primes their conditions
allow the classes in question. The residual image
is full, and for their exceptional \(\ell=5\)
condition a residual scalar of determinant \(4\)
shows that the projective field omits
\(\mathbb Q(\mu_5)\), as in the discussion after
(DA9). The injections into the divisible global
cohomology use \(H^0(\mathbb Q,W)=0\).
[Khare–Patrikis, Definition 3.5 and Theorems 1.1 and 3.14](https://arxiv.org/abs/2611.33361).

The finite groups \(U_{\ell,j}\) thus have bounded
orders for fixed \(E,\ell\), and (RC9a) gives the
same for \(V_{\ell,j}\). By (RC2) the latter tower
stabilizes. If its positive elementary divisor
exponents are \(a_{\ell,1},\ldots,a_{\ell,r_\ell}\),
then its order at depth \(j\) is
\(\ell^{\sum_i\min(j,a_{\ell,i})}\).
The finiteness result is pointwise in the curve,
and the known effective comparison has terms
depending on the original Hecke algebra.

## The denominator supported on nonzero residual groups

Write \(r_\ell=\dim_{\mathbb F_\ell} V_{\ell,1}\) and
\(V_{\ell,\infty}=\varinjlim_j V_{\ell,j}\) under the
coefficient inclusions. By (RC2),
\(V_{\ell,j}=V_{\ell,\infty}[\ell^j]\) under these maps.
The finite groups at any depth therefore have at
most \(r_\ell\) generators, and
\[
              |U_{\ell,j}|\le |V_{\ell,j}|
                               \le\ell^{j r_\ell}.
 \tag{RC10}
\]
Combining (RC1) with the elementary divisor formula
also gives
\[
 e_\ell\le \min(e_\ell,t_2)+
     \sum_{q\in Q_\ell}(3s_q+\lambda_{q,\ell,e_\ell})
          +\sum_i \min(e_\ell,a_{\ell,i}).
 \tag{RC10a}
\]
For residual dimension one, a denominator exponent
exceeding the local correction forces a single dual
class of order at least that excess power of \(\ell\).
In particular every coefficient prime with \(r_\ell=0\)
is removed from the denominator up to the explicit
lower-order factor in (RC6):
\[
       \prod_{\substack{\ell\ge5\\ r_\ell=0}}
                    \ell^{e_\ell}
        \le \max(1,-v_2(j(E)))\!
                 \prod_{q\mid N,\ q\ {\rm odd}} v_q(\Delta_E)^5 .
 \tag{RC11}
\]
This also gives a residue-level certificate for a
wide set of coefficient-prime eliminations when the
remaining valuation factors are coprime to \(\ell\).

The product of \(|U_{\ell,e_\ell}|\) over primes
with \(r_\ell>0\), together with the \(2,3\)-primary
classical degree, would suffice if bounded by
\(O_\varepsilon(R^{2+\varepsilon})\). Equation (RC1)
and the uniform error (RC6) then pass through the
signed degree comparison in (SA7) to full abc.
One may use \(V_{\ell,e_\ell}\) in that sufficient
estimate. The residual rank alone controls the
number of cyclic summands. When \(r_\ell>0\),
the bound \(\ell^{e_\ell r_\ell}\) still involves
the unknown \(e_\ell\), so a uniform estimate for
the length of the actual groups on the surviving
coefficient support remains necessary.

