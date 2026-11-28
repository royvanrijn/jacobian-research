# Local comparison on compact facets for label-indexed agreement kernels

Let \(T\) be the centered regular Gaussian tripod in the plane
and fix \(0<\rho<1\). On any of its facet rays \(e\), take a
fixed compact subinterval \(I_e\) away from the origin. The
intervals may be empty. Let \(E_h\) be the union of the
two-sided strips of normal width \(h\) about these intervals.
For sufficiently small \(h\) the strips are disjoint and the
reference labels on each strip are the two incident labels.

Consider a measurable matrix kernel \(P^{ij}(x,y)\) with
\(1\le i,j\le3\). Assume it has nonnegative entries, is
symmetric in the combined indices, and defines a positive
semidefinite integral operator on \(L^2(\gamma_2;\mathbb R^3)\).
Assume its pair marginals are a simplex-valued vector \(p(x)\)
independent of the second point, and that its same-label row
has mean \(1/3\). More explicitly,
\[
 \begin{gathered}
 P^{ij}(x,y)=P^{ji}(y,x),\qquad
 \sum_j P^{ij}(x,y)=p_i(x),\qquad
 \sum_i p_i(x)=1,\\
 \int \sum_i P^{ii}(x,y)\,d\gamma_2(y)=1/3 .
 \end{gathered}                                               \tag{LP1}
\]
The identities hold almost everywhere for the indicated
Gaussian measures. Suppose \(p\) is the hard tripod label
off \(E_h\) and uses only the two incident labels on each
strip.

With \(U_i=T_\rho1_{T_i}\) and
\[
 L=\int\sum_i p_i(x)(U_{T(x)}(x)-U_i(x))\,d\gamma_2(x),
\]
one has, for all sufficiently small \(h>0\),
\[
       \iint M_\rho(x,y)\sum_i P^{ii}(x,y)\,d\gamma_2^2
                    \le S_\rho(T)-cL,\qquad c>0.             \tag{LP2}
\]
The constants are uniform on compact positive-correlation
bands for fixed intervals \(I_e\). The comparison uses all
three compact strips. In this class equality of objectives
forces the genuine tripod pair kernel.

Retaining the label indices enforces the original signed
facet interactions even when the Gram features have additional
directions. Keeping only the sum of the diagonal label entries
does not enforce these signs. That scalar relaxation admits
a different local mode on several rays.
[Scalar agreement comparison on a single strip](GAUSSIAN_AGREEMENT_COMPACT_FACET_LOCAL_COMPARISON.md)
[Scalar agreement counterexample on three rays](GAUSSIAN_AGREEMENT_MULTI_FACET_COUNTEREXAMPLE.md)

## Exact incidence directions and balance

Write \(P^{ij}(x,y)=\langle V_i(x),V_j(y)\rangle\) in a real
Hilbert space. This feature representation follows from the
spectral theorem. The operator is Hilbert–Schmidt since the
entries are bounded. Its positive finite conditional
compressions have trace at most one. Indeed on an atom \(B\)
their trace contribution is
\(\gamma_2(B)^{-1}\iint_{B^2}\sum_i P^{ii}\,d\gamma_2^2\),
at most \(\gamma_2(B)\) by nonnegativity and (LP1).
Refining the compression shows that the operator is trace
class and gives square-integrable spectral features.

Set \(W(x)=\sum_i V_i(x)\). Its Gram kernel is constantly one.
Hence \(W(x)=v_*\) almost everywhere for some unit vector
\(v_*\). For example, the conditional means of \(W\) on atoms
of any finite partition all have norm one and mutual inner
products one, so they agree. Refining gives the constant
feature sum. The marginal identities then give
\(\langle V_i(x),v_*\rangle=p_i(x)\).
If \(p_i=0\) on a measurable set, \(V_i=0\) there almost
everywhere. In fact the joint entries \(P^{ij}\) vanish
from that set by nonnegativity and (LP1). The integral of
\(V_i\) on any subset then has squared norm zero.

Let \(m_i=\int V_i\,d\gamma_2\). By (LP1),
\(\sum_i\|m_i\|^2=1/3\), while \(\sum_i m_i=v_*\). Equality
in the norm inequality for the mean of three vectors yields
\(m_i=v_*/3\) for every \(i\). Put
\[
              w_i(x)=V_i(x)-1_{\{T(x)=i\}}v_*.
\]
Then each \(\int w_i\,d\gamma_2=0\), and the tuple \(w\)
vanishes outside the strips. If the strips have zero Gaussian
measure this already gives the conclusion. Otherwise, for an
oriented pair \(e=(i,j)\)
on its strip there is a Hilbert vector \(f_e(x)\) with
\[
       w(x)=(e_i-e_j)\otimes f_e(x),                          \tag{LP3}
\]
where \(e_i\) are the coordinate vectors of label space.
Indeed the other feature is zero and \(w_i=-w_j\).
Thus arbitrary Gram directions occur only as the Hilbert
coefficient of the fixed signed label incidence.

For \(q(x)=1-p_{T(x)}(x)\), the dot product of \(f_e\) at
two points on the same side of its strip is the joint
probability of transfer at both points. It lies between
zero and \(\min(q(x),q(y))\). Across its two sides the dot
product is the negative of a nonnegative joint entry and
so is nonpositive. Also on any two strips
\[
          |\langle w(x),w(y)\rangle|
                    \le 2\min(q(x),q(y)).                    \tag{LP4}
\]
In fact at any point with reference label \(i\) and other
incident label \(j\), its tuple is
\((e_j-e_i)\otimes V_j(x)\). The dot product is a label
incidence product of magnitude at most two times the joint
entry of the two transferred labels. That entry is at most
each of its marginals.

## Projecting the signed feature flux to the facets

In the coordinates of strip \(e\), write
\(d\gamma_2=d\nu_e(s)\omega(t)\,dt\), where
\(\omega(t)=e^{-t^2/2}\) and \(d\nu_e\) is Gaussian arclength
on the reference ray. Define
\[
             H_e(s)=\int f_e(s,t)\omega(t)\,dt
\]
on \(I_e\) and zero off \(I_e\). By (LP3), the signed tuple
of these fluxes has exact zero integral for each label. In
particular every Hilbert coordinate of the facet flux is
an ordinary balanced scalar speed.

The opposite-side dot products are nonpositive. On each
same side use its upper bound by \(\min(q(x),q(y))\) and
layer cake. A level set of weighted mass \(u\) on one
normal half-line has weighted first moment at least \(u^2/2\),
since \(\omega\le1\) and filling \([0,u]\) minimizes that
moment. It follows that
\[
   \|H_e(s)\|^2\le 2\int |t|q(s,t)\omega(t)\,dt.             \tag{LP5}
\]
With \(a_e(s)>0\) the magnitude of the raw pair-score
normal derivative, the regret on the strip is the integral
of \(q(s,t)(a_e(s)|t|+O(t^2))\). The Taylor estimate is
uniform on the fixed compact intervals. Therefore
\[
         \sum_e\int a_e\|H_e\|^2\,d\nu_e
                               \le (1+Ch)\,2L .             \tag{LP6}
\]
[The raw interface coefficient and its bounds, (B6)–(B7)](POSITIVE_CORRELATION_CONTINUATION.md#the-full-interface-form-at-a-fixed-positive-correlation)

The exact objective deficit is \(2L-Q(w)\), where
\[
        Q(w)=\iint M_\rho(x,y)\langle w(x),w(y)\rangle\,d\gamma_2^2.
\]
This follows by expanding about the reference features
\(1_{T=i}v_*\), since their linear term pairs with
\(\langle V_i,v_*\rangle=p_i\). Replace the Mehler kernel
in \(Q(w)\) by its restriction to the projected facets.
The two kernels differ by \(O(h)\) on the fixed compact
strips. The estimate (LP4) makes the resulting error
\(O(h)L\). For this, the pushforward of the Gaussian
measure of the strips to their signed transverse coordinate
has density at most \(W=\sum_e\nu_e(I_e)\). A level set
of total mass \(u\) on the strips has first absolute
normal moment at least \(u^2/(4W)\). Integrating the
layer-cake bound for \(\min(q(x),q(y))\) shows that the
integral of the absolute feature interaction without
the Mehler factor is at most \(8W\int |t|q\,d\gamma_2\).

The projected interaction contains the signed factors
\((e_i-e_j)\cdot(e_k-e_\ell)\) exactly, by (LP3).
Apply the planar simplex gap (S37) to every coordinate
of \(H\). In raw normalization the form is
\(\rho\mathfrak q_\rho\) of (B6). The compactly supported
facet speeds have distance squared from the planar
rotation space at least \(\delta\|H\|_{L^2(\nu)}^2\),
where \(\delta>0\) is the fraction of the rotation's
squared norm outside the union of the chosen intervals.
This follows from Cauchy–Schwarz against the rotation
restricted to the intervals. Since \(a_e\) is bounded
above on the intervals, the gap gives an upper bound
\((1-\xi)\sum_e\int a_e\|H_e\|^2\,d\nu_e\) for the
projected interaction, for some \(\xi>0\). The
coordinate sum converges by boundedness of the facet
operators. Combining with (LP6) and the projection error
gives (LP2) for small \(h\).
[The planar simplex spectral gap, (S37)](SIMPLEX_INTERFACE_SYMMETRY.md#reflected-kernels-close-the-three-cell-spectral-problem)

If \(L=0\), the strict reference score gaps almost
everywhere imply \(p\) is the hard tripod label.
The feature-sum identity and vanishing of zero-marginal
features then give the genuine tripod pair kernel.
The preceding proof concerns fixed compact facet strips.
The next comparison allows the tube to include the junction.

## A bounded tubular comparison including the junction

Fix \(0<R<\infty\), and let \(E_{R,h}\) be the part of the
\(h\)-neighborhood of the three reference rays inside \(B_R\).
Keep (LP1), nonnegativity and the label-indexed positive
semidefiniteness hypothesis. Assume only that \(p\) is the hard
reference label off \(E_{R,h}\). The marginals on \(E_{R,h}\)
may use all three labels, including arbitrarily near the
junction. There are \(h_0,c>0\) such that (LP2) holds for
\(0<h<h_0\), uniformly on compact positive-correlation bands.
The proof treats the small junction region as an interaction
error and uses the exact signed spectral comparison on the
remaining strips.

Keep \(w_i=V_i-1_{\{T=i\}}v_*\), which is balanced in every
feature coordinate and vanishes off \(E_{R,h}\).
For ordered transfers \(i\to j\), let
\(q_{ij}=1_{T_i}p_j\) and \(v_{ij}=1_{T_i}V_j\). One has
\[
       w=\sum_{i\ne j}(e_j-e_i)\otimes v_{ij},\qquad
       0\le\langle v_{ij}(x),v_{\ell m}(y)\rangle
                  \le\min(q_{ij}(x),q_{\ell m}(y)).
                                                               \tag{LP7}
\]
This follows from the nonnegative joint pair entries and
the constant feature sum. The estimate does not require the
support labels at either endpoint to be binary.

Let \(J_\delta=E_{R,h}\cap B_\delta\), where
\(0<\delta<\min(1,R/2)\)
will be chosen first. Write \(w_{\rm in}=1_{J_\delta}w\).
On \(T_i\), let \(t_{ij}\ge0\) be the normal distance from
the \(i,j\) supporting line. On \(B_\delta\) the reference
score margin bounds \(U_i-U_j\) below by \(c_0t_{ij}\),
with \(c_0>0\) uniform on the correlation band. In the
coordinate \(t_{ij}\) the pushforward of the Gaussian
measure on \(B_\delta\) has density at most \(C\delta\).
If \(m_{ij}(u)\) is the measure in \(J_\delta\) of the
level set \(q_{ij}>u\), the one-sided normal moment
bound therefore gives
\[
       \int_0^1 m_{ij}(u)^2\,du
                \le C\delta\int_{J_\delta} t_{ij}q_{ij}\,d\gamma_2.
                                                               \tag{LP8}
\]
Apply layer cake to (LP7), and Cauchy–Schwarz to pairs
of the layer masses. Since the Mehler kernel is bounded
on \(B_R^2\), this proves
\[
          Q(w_{\rm in})+\left\|\int w_{\rm in}\,d\gamma_2\right\|^2
                         \le C_R\delta L .                  \tag{LP9}
\]
[The reference score margin](GAUSSIAN_FINITE_SIMPLEX_GERM_PROPAGATION.md#a-comparison-in-a-growing-cylinder)

For \(h\) small compared with \(\delta\), the rest of the
tube lies in disjoint strips with tangential projections
in \([\delta/2,R+1]\) of their respective rays. On each
such strip the unused label \(k\) has a strictly positive
reference score gap, uniformly for the fixed \(\delta\)
and the correlation band. Consequently the sum of its
marginal masses over the exterior strips is at most
\(C_\delta L\).

Merge \(V_k\) into one of the two incident features on each
exterior strip, and set the features equal to the hard
reference features on \(J_\delta\). Let \(\widetilde w\)
be the resulting feature difference. Its pair entries are
nonnegative: merging takes sums of original joint entries,
and pairing with a hard reference feature uses a marginal.
Its sum of features is still \(v_*\). The mass error of the
merged tuple will be corrected on the facets, since the
merged features need not satisfy the exact same-label row
condition.

On the exterior strips the merge error
\(\eta=\widetilde w-w+w_{\rm in}\) is of the form
\((e_i-e_k)\otimes V_k\). Its absolute pair interactions
are bounded by twice the minimum of the corresponding
spectator marginals. The total Gaussian measure of the
strips is \(O_R(h)\). Using this bound and
\(\int p_k\le C_\delta L\) on these strips gives
\[
              Q(\eta)+\left\|\int\eta\,d\gamma_2\right\|^2
                            \le C_{\delta,R}hL.              \tag{LP10}
\]
For instance, the double integral of the minimum is
at most the strip measure times the integral of the
spectator marginal. Boundedness of the Mehler kernel on
the compact strips gives the interaction estimate.

The merged perturbation is binary on every exterior strip.
Its transfer probability \(\widetilde q\) is at most the
original total probability of changing the reference label.
With \(a_e\) the raw pair-score derivative, the reference
margin yields
\[
    2\sum_e\int a_e(s)|t|\widetilde q(s,t)\,d\gamma_2
                                 \le (1+Ch)\,2L .            \tag{LP11}
\]
The incident score gap is \(a_e(s)|t|+O(t^2)\) with
uniform constants on the compact rays. Once \(h\) is
sufficiently small for the fixed \(\delta\), the lower
bound for the spectator gap is larger than this pair cost.
Thus counting a merged spectator transfer as an incident
transfer cannot increase the bound in (LP11).

Project \(\widetilde w\) to the facets as in (LP5).
The projected flux \(H\) retains the exact incidence
vectors and \(\sum_e\int a_e\|H_e\|^2\,d\nu_e\le
(1+Ch)2L\). The compact kernel projection error is
\(O_R(h)L\), by the same layer-cake bound using
\(\widetilde q\). In particular \(Q(\widetilde w)=O_R(L)\)
with a constant independent of \(\delta\) once \(h\)
obeys the preceding restriction. The positive Mehler
form, (LP9), and (LP10) then give
\[
       Q(w)\le Q(\widetilde w)
                   +C_R\big(\sqrt\delta+C_{\delta,R}\sqrt h\big)L .
                                                               \tag{LP12}
\]
Here the mixed terms use Cauchy–Schwarz in the Mehler form
between \(\widetilde w\) and \(w_{\rm in}-\eta\).

The signed label means of the projected flux have total
Hilbert norm at most \(C_R(\sqrt\delta+
C_{\delta,R}\sqrt h)\sqrt L\), by the original exact
balance and (LP9)–(LP10). Remove this mass error using two
fixed compact regular facet bumps on a label spanning
tree. Their mass fluxes form a basis of the centered label
space, so the correction has the same bound in facet
\(L^2\), also for Hilbert-valued coefficients. The
corrected supports lie in fixed bounded ray intervals
independent of \(\delta\). Thus their gap to the planar
rotation is bounded below independently of \(\delta\),
by the rotation's squared norm outside this fixed
support. Apply (S37) to the corrected scalar coordinates
as before. Both the interaction and local form change
by at most \(C_R(\sqrt\delta+C_{\delta,R}\sqrt h)L\)
under the correction, by boundedness of the facet forms
and the flux norm bound from (LP11). Consequently
\[
         Q(w)\le (1-\xi)2L
                   +C_R\big(\sqrt\delta+C_{\delta,R}\sqrt h\big)L,
                    \qquad \xi>0 ,
\]
where \(\xi\) can be fixed before choosing \(\delta\).
Choose \(\delta\) small and then \(h\) small. This gives
the claimed strict comparison. Its equality statement is
the same as for (LP2).
[Fixed tree-facet repair](POSITIVE_CORRELATION_CONTINUATION.md#an-effective-local-radius-from-truncated-normal-flux)

The hard marginals outside a bounded tube are a support
hypothesis. Closeness of the marginals or of the pair
kernel in norm does not give that hypothesis. In the
convex unrestricted relaxation, a better feasible kernel
can be mixed with the tripod kernel in arbitrarily small
proportion.
