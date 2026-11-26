# Thin sets of relative Frobenius traces

Let \(A/K\) be an abelian variety over a number field with
\(\operatorname{End}^0_{\overline K}(A)=F\), where \(F\) is totally real of
degree \(g=\dim A\). Assume all endomorphisms are defined over \(K\).
At a good absolute degree-one place \(v\) over \(p\), use the relative
Frobenius polynomial \(X^2-a_vX+p\), with \(a_v\in\mathcal O_F\).
Choose an integral basis to view \(a_v\) as a point of \(\mathbb A^g(\mathbb Q)\).
For \(g=3\) put
\[
 b_v=\frac{N_{F/\mathbb Q}(a_v)}p,\qquad
 d_v=-\frac{N_{F/\mathbb Q}(a_v^2-4p)}p.
\]
At a split one-failure place both are integers, and \(d_v\) is the
[raw discriminant cofactor](CM_NEWTON_SPLIT_COFACTOR_NORM_AND_TRACE_COUNT.md).

A thin set here means a subset of a finite union of proper closed loci
and images of rational points under generically finite covers of degree
at least two with geometrically integral source. A \(\mathbb Q\)-irreducible
source that is geometrically reducible contributes only a proper closed
locus, since its rational points lie in intersections of geometric
components.
The sieve extends Serre's trace-and-determinant argument for elliptic
curves, with the quantitative version due to Zywina for tuples of
nonisogenous elliptic curves. Fixed Frobenius splitting fields are
treated by Burungale, Hida and Lai under more general monodromy hypotheses.
[Serre, §3.6, Theorems 3.6.2 and 3.6.4](https://www.ms.uky.edu/~sohum/ma561/notes/workspace/books/serre_galois_theory.pdf)
· [Zywina, Theorem 1.5 and §5](https://pi.math.cornell.edu/~zywina/papers/LargeSieve.pdf)
· [Burungale–Hida–Lai, Theorems 1.1 and 1.3](https://msp.org/ant/2026/20-3/ant-v20-n3-p05-p.pdf).

## Trace tuple and its effective count

**Theorem.** For a fixed thin \(\Omega\subset\mathbb A^{g+1}(\mathbb Q)\)
the good absolute degree-one places \(v\) with \((p,a_v)\in\Omega\)
have natural density zero among the primes of \(K\). For \(g=3\) their
count up to \(p\le X\) is
\[
 O_{A,K,\Omega}\left(
       \frac{X(\log\log X)^{31/30}}{(\log X)^{61/60}}\right).
\]
Under GRH for the Dedekind zeta functions of the finite torsion composita
used in the sieve, the count is
\[
 O_{A,K,\Omega}\bigl(X^{49/50}(\log X)^{2/25}\bigr).
\]

For a sufficiently large rational prime \(\ell\) split completely in
\(F\) and unramified in \(K\), the residual image on \(A[\ell]\) is
\[
G_\ell=\{(g_i)_{i=1}^g\in\mathrm{GL}_2(\mathbb F_\ell)^g:
                       \det(g_i)=\det(g_j)\}.
\]
The trace vector is \(a_v\bmod\ell\) under the split identification of
\(\mathcal O_F/\ell\mathcal O_F\), and the common determinant is \(p\bmod\ell\).
The residual image follows by reduction from Lombardo's large-image
theorem; discard primes dividing the index of the endomorphism order
in \(\mathcal O_F\) if needed.
[Lombardo, Theorem 1.4](https://arxiv.org/pdf/1411.1703).

These residual images are independent across any finite set of the
indicated primes \(\ell>3\). In each \(G_\ell\), the commutator subgroup
is \(S_\ell=\mathrm{SL}_2(\mathbb F_\ell)^g\), because \(S_\ell\) is
perfect and the determinant quotient is abelian. The joint commutator
projects onto each \(S_\ell\). Any nontrivial common quotient of one
\(S_\ell\) and a product of \(S_{\ell'}\) for distinct \(\ell'\) would
be perfect and have a nonabelian composition factor. The only such
factor on the first side is \(\mathrm{PSL}_2(\mathbb F_\ell)\), which
is absent on the second side. Goursat's lemma and induction give the
full product of the commutator groups. The joint determinant is
surjective by the cyclotomic character. The cyclotomic field of the
squarefree product of the auxiliary primes meets \(K\) trivially,
since their intersection is unramified at every finite prime of
\(\mathbb Q\). This proves independence.

At fixed determinant \(e\ne0\) and trace \(t\), at least
\(\ell(\ell-1)\) matrices occur. Indeed choose a first diagonal entry
\(x\) and nonzero upper-right entry \(y\), then set the other diagonal
entry to \(t-x\) and the lower-left entry to \((x(t-x)-e)/y\). Thus each
trace in a determinant fiber has probability at least \(1/(\ell+1)\).
The joint trace vector conditional on the determinant has product law.
Each point of \(\mathbb F_\ell^\times\times\mathbb F_\ell^g\) therefore
has probability at least \((\ell-1)^{-1}(\ell+1)^{-g}\) under \(G_\ell\).

Treat each thin component separately. For a proper closed locus, its
reductions have \(O_\Omega(\ell^g)\) points. For a cover component,
Serre's Theorem 3.6.2 gives a fixed finite extension such that at
sufficiently large \(\ell\) split in its normal closure the reduction
misses a fixed positive fraction of \(\mathbb F_\ell^{g+1}\).
The same assertion applies to reductions of rational points integral
at \(\ell\). Over a suitable open subset in the target, spread out the
cover as a finite morphism. An \(\ell\)-integral target point whose
reduction lies in this open subset has a lift integral over \(\mathbb
Z_\ell\) whenever it has a rational lift in the cover, by finiteness.
The excluded closed subset contributes \(O(\ell^g)\) further residues.
We further select \(\ell\) split in the normal closure of \(F\) and
discard determinant zero, a proportion \(1/\ell\). The trace atom
bound proves a uniform upper bound \(c<1\) for the proportion of
\(G_\ell\) admitted by this thin component. Independence and
Chebotarev for any \(m\) such primes give upper density at most \(c^m\).

For the effective estimates apply Zywina's Theorem 3.3 to these
auxiliary primes. The system is independent as proved above, and
each image is unramified away from a fixed set and the places over
\(\ell\). The admitted conjugation-invariant subsets have density
at most \(c<1\); enlarge any empty one to one conjugacy class if
needed. The product group has size at most \(Q^{3g+1}\) when the
product of the auxiliary primes is at most \(Q\).
The class number estimate in Zywina's Lemmas B.1–B.2 gives
\(|G_\ell^\#|\ll_g \ell^{g+1}\), where \(G_\ell^\#\) is the set of
conjugacy classes. Use just the singleton auxiliary subsets in the
sieve, so its denominator is \(\gg Q/\log Q\) since the auxiliary
primes run through a fixed completely split set. The size exponents
\(r=3g+1\) and \(s=g+1\) give the unconditional rate in Theorem
3.3(i) and the GRH rate in Remark 3.4(iv), with \(g=3\). Summing
over the finite union of thin components proves the result.
[Zywina, Theorem 3.3, Remark 3.4(iv), §5.2 and Appendix B](https://pi.math.cornell.edu/~zywina/papers/LargeSieve.pdf).

## Joint cofactor consequences

**Corollary.** For \(g=3\) the same counts hold for places with
\((p,b_v)\) in any fixed thin subset of \(\mathbb A^2(\mathbb Q)\), or
with \((p,b_v,d_v)\) in any fixed thin subset of \(\mathbb A^3(\mathbb Q)\).
In particular these exclusions hold on the split one-failure locus.

In applying the trace tuple theorem it suffices to work on the open
subset \(p\,N_F(a)\ne0\). A dominant morphism with geometrically integral
generic fiber pulls a thin set back to a thin set. Indeed it pulls
a proper closed locus back to a proper closed locus, and its regular
extension of function fields preserves the degree and geometric
integrality of a generically finite cover under base change over the
common open subset.

The map to \((p,b)\) has geometrically integral generic fiber. Over an
algebraic closure it fixes \(p\) and the nonzero product \(a_1a_2a_3\),
leaving a two-dimensional torus. Directly, for fixed nonzero residues of \(p\)
and \(b\) modulo a split \(\ell\), there are \((\ell-1)^{g-1}\) trace
tuples with the required product, so each such pair has probability
at least \((\ell-1)^{g-2}/(\ell+1)^g\).

For the triple map, work over an algebraic closure of the parameter
field \(\mathbb Q(p,b,d)\) and scale \(x_i=a_i/(2\sqrt p)\). The equations are
\[
   x_1x_2x_3=n=\frac b{8\sqrt p},\qquad
   (1-x_1^2)(1-x_2^2)(1-x_3^2)=\delta=\frac d{64p^2}.
\]
The map is dominant because at fixed generic \(p,n\) the second
product varies on the product-\(n\) torus. Eliminate \(x_3\) and put
\(x=x_1\), \(y=x_2^2\). The equation becomes
\[
 -(1-x^2)x^2 y^2
 +\big((1-x^2)(x^2+n^2)-\delta x^2\big)y
 -(1-x^2)n^2=0.
\]
Its discriminant in \(y\) is
\[
 \big((1-x^2)(x-n)^2-\delta x^2\big)
 \big((1-x^2)(x+n)^2-\delta x^2\big).
\]
Each quartic has simple roots because \(\delta\) is transcendental
over \(\mathbb Q(n)\) and the rational function
\((1-x^2)(x\pm n)^2/x^2\) is nonconstant and separable. The quartics
have disjoint roots: a common zero would occur at \(0,1,-1\), and
neither quartic vanishes there. Thus the discriminant is nonsquare
even over the algebraic closure of the parameter field. In the
quadratic field generated by \(y\), there is a place over \(x=1\)
with \(v(y)=1\). The linear coefficient is a unit there and the
constant coefficient has valuation one, so Hensel's lemma gives
the simple branch with this valuation. Hence \(y\) is nonsquare in
that quadratic field. This proves the quartic equation for \(x_2\)
irreducible over the geometric rational function field in \(x\).
Its coefficients have no common zero on \(x\ne0\), because at
\(x=\pm1\) the middle coefficient in \(y\) is nonzero. Gauss's lemma
on the torus proves geometric integrality of the generic fiber.

## The quotient by a moving coefficient prime

**Proposition.** For each ideal class of \(F\), choose once and for all an
ideal representative \(J\). For a split \(p\) and \(\mathfrak q\mid p\),
take \(J\) in its class and write \(\mathfrak q=\beta_{\mathfrak q}J\).
Let \(\Theta\subset F\), viewed in a fixed \(\mathbb Q\)-basis of \(F\),
be a fixed thin subset invariant under multiplication by the units
of \(\mathcal O_F\). The good absolute degree-one places \(v\) above
split \(p\) for which
\[
        a_v/\beta_{\mathfrak q}\in\Theta
        \quad\text{for some }\mathfrak q\mid p
\]
have natural density zero. The condition is independent of the chosen
generator \(\beta_{\mathfrak q}\). It applies in particular to the
quotient by the failed coefficient prime.

Enlarge \(K\) to contain the normal closure of \(F\). At each selected
degree-one place of the original field the completion is \(\mathbb Q_p\),
and the normal closure of \(F\) splits completely over \(\mathbb Q_p\).
Thus the base extension gives lifts of the same norm. At a chosen lift
\(v\), every
\(\mathfrak q\mid p\) is obtained by restricting \(v\) along some
embedding of \(F\) into the enlarged \(K\).

Choose a finite set of sufficiently large auxiliary primes \(\ell\)
split in \(F\) and avoiding the norms of all \(J\), and let \(M\) be
their product. Take the compositum \(R/K\) of the ray class fields
of all the corresponding embedded coefficient fields with modulus \(M\).
It is abelian over \(K\). At a place \(v\) away from its ramification,
its Frobenius determines the ray classes of the corresponding primes
\(\mathfrak q\) because \(v\) has residue degree one. The ray class
of \(\mathfrak qJ^{-1}\) determines the residue \(\beta_{\mathfrak q}\)
in \((\mathcal O_F/M)^\times\) up to the image of global units,
by the kernel of the ray-to-ordinary class group map.

Use the actual simultaneous reductions \(\Theta_M\) modulo \(M\) of
elements of \(\Theta\) integral at all selected auxiliary primes.
The set \(\Theta_M\) is stable under the image of the global units.
Thus the required residue condition is
\[
             a_v\bmod M \ \in\ \beta_{\mathfrak q}\Theta_M,
\]
which depends on the ray information alone for its choice of scalar
coset. If \(\Theta\) is contained in a finite union of basic thin sets
\(\Theta_j\), then \(\Theta_M\) is contained in
\(\bigcup_j\prod_{\ell\mid M}\Theta_{j,\ell}\), where each
\(\Theta_{j,\ell}\) is the reduction of the \(\ell\)-integral points of
that basic set. For each \(j\), choose \(m\) distinct primes where the
reduction misses a fixed positive fraction as above. The individual
\(\Theta_{j,\ell}\) need not be unit invariant.

Let \(H\) denote the joint image of the torsion representations and the
ray compositum. It surjects to the product of all \(G_\ell\). Because
the ray projection is abelian, \(H'\) contains the full product of
the \(S_\ell\) with trivial ray projection. Consequently, conditional
on any compatible values of the common determinants and the ray
Frobenius, the individual matrix blocks in their determinant fibers
are independent and uniform. Multiplication by an auxiliary unit
\(\beta_{\mathfrak q}\bmod\ell\) is a permutation of the residue trace
space. The trace atom bound now gives a conditional probability at
most \(c_j^m\), with \(c_j<1\), for the residue conditions of each
basic set after the chosen scalar multiplication. The bound holds
uniformly in the ray and determinant data. Add over the finite
number of pieces and embeddings. Chebotarev in the fixed compositum,
then \(m\to\infty\), gives the asserted upper density zero.

### Canonically balanced positive generators

With a fixed choice of generator, the quotient exclusion holds for
every fixed thin \(\Theta\subset F\). Choose a representative \(J\) of
each narrow ideal class and a half-open fundamental parallelepiped
\(P\) for the logarithmic lattice of the totally positive units.
For \(\gamma>0\) at all real embeddings, write
\[
 z(\gamma)=\bigl(\log|\sigma(\gamma)|\bigr)_\sigma
             -\frac{\log N_{F/\mathbb Q}(\gamma)}g(1,\ldots,1).
\]
Let \(\Lambda=z(\mathcal O_F^{\times,+})\) in the real sum-zero plane
\(V\). For each \(\mathfrak q\) use the unique totally positive
\(\beta_{\mathfrak q}\) with \(\mathfrak q=\beta_{\mathfrak q}J\) and
\(z(\beta_{\mathfrak q})\in P\). Then the places with
\(a_v/\beta_{\mathfrak q}\in\Theta\) for at least one coefficient
prime over split \(p\) have natural density zero.

For the proof, keep \(M\) fixed and use the narrow ray class fields
of modulus \(M\) in the compositum above. In each narrow ray class
\(c\) choose an ideal \(J_c\) prime to \(M\), and choose
\(\kappa_c>0\) with \((\kappa_c)=J_cJ^{-1}\), where \(J\) represents
the narrow class of \(J_c\). If \(\mathfrak q\) lies in \(c\), choose
\(\gamma_{\mathfrak q}>0\) with
\((\gamma_{\mathfrak q})=\mathfrak qJ_c^{-1}\) and
\(\gamma_{\mathfrak q}\equiv1\bmod M\). Set
\(\Lambda_M=z(\{\varepsilon\in\mathcal O_F^{\times,+}:
\varepsilon\equiv1\bmod M\})\). The angle
\(z(\gamma_{\mathfrak q})\bmod\Lambda_M\) is well defined. Reducing
\(z(\kappa_c)+z(\gamma_{\mathfrak q})\) into \(P\) selects the positive
unit \(\varepsilon\) with
\(\beta_{\mathfrak q}=\kappa_c\gamma_{\mathfrak q}\varepsilon\).
Thus \(\beta_{\mathfrak q}\bmod M=\kappa_c\varepsilon\bmod M\) is a
piecewise constant function of the angle and \(c\). Its boundaries
lie on finitely many translated faces of \(P\) in \(V/\Lambda_M\)
and have Haar measure zero.

Fix one embedding of \(F\) into \(K\) after the base extension above,
and let \(L/K\) be the finite Galois compositum of the torsion and ray
fields at the selected moduli. The angle is Haar equidistributed
conditional on each Frobenius conjugacy class in \(\operatorname{Gal}(L/K)\).
Indeed each nonzero Fourier mode \(\lambda\in\Lambda_M^\vee\)
defines the character \((\gamma)\mapsto e^{2\pi i\lambda(z(\gamma))}\)
on the narrow ray-principal subgroup. It extends to a unitary Hecke
character \(\psi_\lambda\) of \(F\) on ideals prime to \(M\), since
the quotient by that subgroup is finite and one can choose roots
in the unit circle. On a fixed ray class \(c\), its value at
\(\mathfrak q\) equals the angle mode times \(\psi_\lambda(J_c)\).

Expand the Frobenius class indicator in irreducible characters and
use Brauer induction. For an induced term
\(\operatorname{Ind}_{H}^{\operatorname{Gal}(L/K)}\eta\), the fixed
cosets of Frobenius correspond to primes \(w\) of \(E=L^H\) of
relative residue degree one over the prime of \(K\). The twisted
prime sum is then the prime sum in \(E\) for the unitary Hecke
character \(\eta\,(\psi_\lambda\circ N_{E/F})\), up to primes of
higher residue degree and finitely many ramified primes.
The radial archimedean exponent vector of \(\psi_\lambda\) is nonzero
with sum zero. On norm pullback the exponent at an archimedean
place, using the normalized absolute value, is inherited from the
place of \(F\) below it. The exponents still differ across places.
The finite-order factor \(\eta\) leaves these exponents unchanged.
For this fixed character the conductor and frequencies are fixed,
and the nonzero archimedean part excludes the real exceptional
character. The Hecke prime ideal theorem, with prime powers removed
and partial summation, gives \(o(X/\log X)\). The removed
higher-degree primes contribute \(O_{L}(\sqrt X)\), including the
restriction to absolute degree one in \(K\). The zero mode is
Chebotarev, so the Fourier criterion gives the conditional Haar law.
[Brauer, *On Artin's L-series with general group characters*](https://doi.org/10.2307/1969183)
· [Maynard, *Primes represented by incomplete norm forms*, Lemmas 4.4–4.6](https://doi.org/10.1017/fmp.2019.8).

Conditional on the ray data and determinants, the same product law
of the matrix blocks therefore holds after fixing the residue
\(\beta_{\mathfrak q}\bmod M\). For each admissible value of that
residue, multiplication by it permutes the auxiliary trace space.
The preceding missing-residue argument applies to the translated
sets \(\beta_{\mathfrak q}\Theta_M\), with the same bound uniformly in
the residue, ray data and determinants. Take the prime limit before
letting the number \(m\) of auxiliary primes per thin component
grow. Applying the argument to each of the finitely many embeddings
of \(F\) gives the stated conclusion.

## Scope

For \(p\) split in a cubic coefficient field, \(b_v\) and \(d_v\) being
integers impose divisibility by the residue characteristic. The thin
exclusions concern fixed algebraic restrictions on their values. The
local candidate space of full trace tuples \((p,a)\) with one failure
is itself nonthin. For each basic component of any proposed thin
superset, select a separate auxiliary prime \(\ell\) and a missing
residue tuple \((p\bmod\ell,a\bmod\ell)\) with \(p\bmod\ell\ne0\).
Chebotarev gives arbitrarily large \(p\) split in \(F\) with all the
specified \(p\)-residues, since the auxiliary moduli can be chosen
coprime to the discriminant of the normal closure. The
[moving-level lattice count](CM_NEWTON_SPLIT_COFACTOR_NORM_AND_TRACE_COUNT.md#counting-candidate-traces-with-a-moving-auxiliary-level)
gives a one-failure trace in the required trace residue classes for
each such sufficiently large \(p\). The images of this candidate
space in \((p,b)\) and \((p,b,d)\) are nonthin as well, by the
generic-fiber arguments. The sparsity of split one-failure primes of one
fixed \(A\) remains open.
