# Coupled local-global principles for splitting and primitive norm values

A coupled system must be locally solvable with the same shared variables.
Separate local tests for its individual equations do not establish this.
For fixed linear systems in abelian splitting-semigroup values, the correct
criterion consists of joint residue states and finitely many coefficient-prime
tests, together with a nondegenerate positive real region.

The main counting result below is a classical finite-complexity consequence
of Matthiesen's moving-target theorem and the Frobenian fixed-character
refinement, already scoped in the [Goldbach theory](SPLITTING_SEMIGROUP_LOCAL_GLOBAL.md#8-higher-goldbach-counts-the-classical-finite-complexity-range).
It does not use the 2026 binary graph or zero-free inputs.
Its explicit application gives an eventual local-global law for two shared
triple sums of separately primitive Gaussian norms: exactly 45 residue pairs
modulo eight are admissible. All counts here are of numerical value vectors,
rather than quadratic-coordinate representation weights.

The statements and source applications have internal written review.
No external human review, formal verification, historical priority claim,
numerical onset or complete small-exception classification is asserted.

## 1. Integer fibers and the geometric hypothesis

Fix an even modulus \(L\), subgroups \(H_i\le(\mathbb Z/L)^\times\),
filters \(U_i\subset H_i\), and splitting cores
\[
 C_i=\langle p:p\nmid L,\ p\bmod L\in H_i\rangle,\qquad
 \rho_i=|H_i|/\varphi(L),\qquad
 c_i=\frac{\lim_{z\downarrow1}(z-1)^{\rho_i}
                      \sum_{n\in C_i}n^{-z}}{\Gamma(\rho_i)}.
\]
All data are fixed. Let \(A:\mathbb Z^s\to\mathbb Z^r\) be surjective,
and put \(d=s-r\ge2\). Choose an integer right inverse \(R_0\), and a
basis matrix \(B\) for \(\ker_{\mathbb Z}A\). Then
\[
 x=R_0\boldsymbol n+Bt=:\ell_{\boldsymbol n}(t),
 \qquad t\in\mathbb Z^d                                  \tag{C.1}
\]
parametrizes the integer fiber \(Ax=\boldsymbol n\) bijectively.
Indeed every integer vector decomposes uniquely into a kernel vector
and \(R_0Ax\); hence \([B\ R_0]\) is unimodular. The same parametrization
is bijective modulo every integer modulus.

Assume the rows \(b_i\) of \(B\) are nonzero and pairwise linearly
independent over \(\mathbb Q\). Equivalently, no nonzero vector of
\(\operatorname{row}_{\mathbb Q}(A)\) is supported on fewer than three
coordinates: a relation between one or two coordinate restrictions
to the kernel is exactly such a row-space vector.
Thus the condition is independent of the chosen integer kernel basis.

This hypothesis is substantive. Two free parameters alone do not
guarantee pairwise independence; shared private sums can leave two
coordinate forms proportional.

## 2. A complete finite congruence criterion

For each pair \(i<j\), let \(D_{ij}>0\) be the gcd of the absolute
two-by-two minors of the matrix with rows \(b_i,b_j\).
At least one minor is nonzero. Put
\[
 {\cal P}_*=\{p:p\nmid L,\ p\le s\ \text{or}\ p\mid\prod_{i<j}D_{ij}\},
 \qquad M_*=L\prod_{p\in{\cal P}_*}p.
\]
These extra primes are used for local testing; they are not removed
from the original prime prescriptions.

For \(p\nmid L\), let \(I_p\) consist of the coordinates whose cores
forbid \(p\), and define
\[
 \begin{split}
 \gamma_L(\boldsymbol n)
 &=\frac{L^r}{\prod_i|H_i|}
   \#\{h\in\prod_iU_i:Ah=\boldsymbol n\bmod L\},\\
 \alpha_p(\boldsymbol n)
 &=p^{-d}\#\{t\bmod p:\ell_{\boldsymbol n,i}(t)\ne0
                                     \text{ for all }i\in I_p\},\\
 \sigma_p(\boldsymbol n)
 &=\frac{\alpha_p(\boldsymbol n)}
                {\prod_{i\in I_p}(1-1/p)}.                \tag{C.2}
 \end{split}
\]

**Theorem CoupledFiniteAdmissibility.** The following are equivalent:

1. \(\gamma_L(\boldsymbol n)>0\) and
   \(\alpha_p(\boldsymbol n)>0\) for every \(p\in{\cal P}_*\);
2. for every modulus \(m\), actual values
   \(x_i\in C_i(U_i)\) solve \(Ax=\boldsymbol n\bmod m\);
3. such a solution exists modulo the single modulus \(M_*\).

**Proof.** Any solution modulo \(M_*\) supplies the necessary root states
and forbidden-unit tests. Conversely, outside \({\cal P}_*\) each
nonzero forbidden form removes a proper affine hyperplane modulo \(p\).
There are at most \(s<p\) of them, so their union cannot cover
\(\mathbb F_p^d\). Thus a unit root exists at every good prime.

A root modulo \(p\) can be lifted to any \(p\)-power by lifting \(t\)
in (C.1); all forbidden coordinates remain units. At primes dividing
\(L\), fix the joint root modulo \(L\) and lift its parameter \(t\);
this preserves every required unit state. Combining these choices by
CRT gives a fiber residue vector modulo \(\operatorname{lcm}(m,L)\).
Each coordinate has an allowed \(L\)-state and is a unit at every
forbidden good prime. The finite-image lemma of the Goldbach theory
realizes it by an actual value of \(C_i(U_i)\). These values need only
satisfy the equation modulo \(m\). This proves all finite congruences.
Taking \(m=M_*\) proves the remaining implication. \(\square\)

This is a joint all-moduli criterion. Projecting the subgroup states
to one prime at a time can lose their coupling across the factors of
\(L\), as the modulus-twelve example in the Goldbach theory shows.

## 3. Pointwise counts in positive real regions

Let \(K\subset\mathbb R^d\) be a fixed bounded convex set of positive
volume. Consider targets \(\boldsymbol n=O(T)\) such that, throughout
\(t\in TK\),
\[
                0<\ell_{\boldsymbol n,i}(t)\le C T
                   \quad(1\le i\le s)                    \tag{C.3}
\]
for a fixed \(C\). Write \({\cal R}_K(\boldsymbol n;T)\) for the count
of these integer parameters with all \(\ell_i(t)\in C_i(U_i)\).
Set \(\kappa=s-\sum_i\rho_i\).

**Theorem CoupledSplittingCount.** Under the stated hypotheses,
\[
 {\cal R}_K(\boldsymbol n;T)=
 \operatorname{vol}(K)\left(\prod_i c_i\right)
 \gamma_L(\boldsymbol n)
 \left(\prod_{p\nmid L}\sigma_p(\boldsymbol n)\right)
 \frac{T^d}{(\log T)^\kappa}
 +o\!\left(\frac{T^d}{(\log T)^\kappa}\right).              \tag{C.4}
\]
The error is uniform over the targets satisfying (C.3).
On the finite-admissible classes the Euler product and its inverse
are uniformly bounded, and the displayed asymptotic is relative.
Failure of any finite test gives zero count identically.
No numerical logarithmic rate is claimed.

**Proof.** Use Matthiesen, *Linear correlations of multiplicative
functions*, arXiv:1606.04482v4, Theorem 2.1, pp. 6–7, at \(B_0=0\),
and Remark 2.2. The intercepts \(R_0\boldsymbol n=O(T)\) are explicitly
permitted by that theorem. Rescale the source parameter by a fixed
constant if necessary to put (C.3) in its range \((0,X]\); the product
of source volume and \(X^d\) remains \(\operatorname{vol}(K)T^d\).
The logarithms differ by \(O(1)\).

The class-\({\cal F}\) verification, fixed-filter character expansion,
and fixed-character \(W\)-refinement are exactly the ones proved in
Section 8 of the Goldbach theory, using Loughran–Matthiesen,
Proposition 3.3, Lemmas 3.4–3.5 and pp. 4799–4800 of the published
2024 edition. Choose \(w(T)=\log\log T\).
The pairwise-independent nonconstant rows verify the finite-complexity
hypothesis. Unlike the binary case, no graph transfer is used.

At \(L\), each marginal state has density \(L/|H_i|\) relative to its
unrestricted marginal. The number of parameter roots modulo \(L\)
equals the number of fiber roots in (C.2), by (C.1).
Its parameter density is that number divided by \(L^d\).
The product is exactly \(L^{s-d}/\prod_i|H_i|\) times the root count,
which is \(\gamma_L\).

At a good prime an allowed coordinate has local weight one; a forbidden
coordinate has weight \(p/(p-1)\) on the units. Haar averaging in the
parameter space gives exactly \(\sigma_p\). In a coupled system these
factors need not reduce to primes forbidden on every side: that
simplification used the geometry of a single sum equation.

Outside a fixed finite set the normals of every pair of forbidden
hyperplanes are independent modulo \(p\). Each hyperplane has
probability \(1/p\), and every intersection of at least two is either
empty or has probability at most \(p^{-2}\). Inclusion-exclusion gives
\[
 \alpha_p=1-|I_p|/p+O_s(p^{-2}),\qquad
                     \sigma_p=1+O_s(p^{-2}),              \tag{C.5}
\]
uniformly in the target. This proves uniform convergence of the Euler
tail. The finitely many small-prime factors depend only on target
residues modulo those primes; their positive values have a positive
minimum. The same is true of \(\gamma_L\).

For completeness, the source's exact small-prime prefix truncation is
uniform here even at coefficient primes. Put
\(e_i(p)=\min_jv_p(B_{ij})\), finite since \(b_i\ne0\).
For Haar \(t\in\mathbb Z_p^d\),
\[
 \mathbb P(v_p(\ell_i(t))\ge a)\le p^{e_i(p)-a}
       \ (a>e_i(p)),\qquad
 \mathbb E\,v_p(\ell_i(t))\le e_i(p)+1/(p-1).
\]
Condition on all forbidden coordinates being units. At each fixed
exceptional prime its positive probability is at least \(p^{-d}\).
For all sufficiently large primes it is at least \(1-s/p\ge1/2\),
and \(e_i(p)=0\). Conditional valuation expectations are consequently
\(O(1)\) at the finitely many exceptional primes and \(O_s(1/p)\)
elsewhere, uniformly in all target valuations. At \(p\mid L\) all
coordinates are units.

For the source \(W\le(\log T)^{B_1}\), each allowed prefix \(w_i\)
therefore has \(\mathbb E\log w_i=O(\log\log\log T)\).
Split \(\sum_{p\mid W}\log p/(p-1)\) at \(\log\log T\), as in the
Goldbach proof. Markov bounds the mass of
\(w_i>(\log T)^{B_2}\) by
\(O(\log\log\log T/\log\log T)=o(1)\).
Absolute filter tails are bounded by the unfiltered nonnegative model.

Thus, for every fixed source accuracy \(\varepsilon>0\), the finite
model becomes (C.4), independently of \(W,B_2\), up to uniform \(o(1)\).
The fixed-character refinement and Remark 2.2 let \(\varepsilon\)
tend to zero. This proves (C.4), its uniformity and positivity claims.
\(\square\)

**Corollary CoupledBulkLocalGlobal.** Let \(\Omega\) be a compact subset
of \(A(\mathbb R_{>0}^s)\). For all sufficiently large \(T\), every
\(\boldsymbol n\in T\Omega\cap\mathbb Z^r\) satisfying the finite tests
has a representation \(Ax=\boldsymbol n\) with \(x_i\in C_i(U_i)\).
Indeed it has \(\gg_\Omega T^d/(\log T)^\kappa\) distinct value vectors
in bounded positive regions.

**Proof.** At each \(y\in\Omega\), choose \(z>0\) with \(Az=y\), and
write \(z=R_0y+B u\). A small parameter box around \(u\) stays positive
for targets in a neighborhood of \(y\), with a common upper bound.
A finite cover of \(\Omega\) supplies finitely many fixed boxes \(K\)
to which (C.4) applies uniformly. Their volumes and positive local
factors have a positive minimum. The parametrization is injective.
\(\square\)

Real solubility without this positive-volume control gives no uniform
count: a target may approach a boundary where the available real
region collapses.

## 4. Image lattices and finite norm profiles

Surjectivity is a normalization, not an assumption to forget after
scaling variables. For full-row-rank \(A\), let \(Q\) be an integer
basis matrix for \(A\mathbb Z^s\). First require
\(\boldsymbol n\in Q\mathbb Z^r\). Then
\[
 A'=Q^{-1}A:\mathbb Z^s\longrightarrow\mathbb Z^r,\qquad
 \boldsymbol n'=Q^{-1}\boldsymbol n
\]
are integral and \(A'\) is surjective, with the same kernel.
The image-lattice condition is a finite cokernel congruence.
It can be found by Smith normal form.
All-moduli solubility for the original and normalized systems is
equivalent: multiplication by \(Q\) gives one direction; using modulus
\(|\det Q|m\) and its adjugate gives the other.
Real feasibility transforms by the same fixed invertible map.

**Corollary CoupledFiniteProfiles.** Replace each core by a finite
disjoint union \(P_i=\coprod_{u\in{\cal U}_i}uC_i(U_{i,u})\), with
positive fixed \(L\)-supported prefixes. Under the row-space and
positive-region hypotheses, eventual joint all-moduli admissibility
again implies global representability. Admissibility is a finite
congruence condition. On compact positive target regions the number
of value vectors is \(\gg T^d/(\log T)^\kappa\) on admissible classes.

**Proof.** Fix a prefix tuple. Its matrix is
\(A\operatorname{diag}(u_i)\). Normalize its finite-index image first.
Nonzero coordinate scaling and an invertible row transformation preserve
the row-space support hypothesis. Apply the preceding results.
The positive real image cone is unchanged by positive diagonal scaling.
There are finitely many prefix tuples. Each has a finite obstruction
modulus, including its image-lattice condition; a common multiple tests
them all. If every tuple failed, these finite obstructions together
would exclude the whole union. Thus all-moduli admissibility picks
one compatible stratum, to which the bulk theorem applies.
Thresholds and lower bounds are uniform over the finite list.
\(\square\)

This argument keeps profiles of shared variables consistent across all
equations. It does not treat unrestricted infinite conductor profiles.

## 5. A complete coupled primitive-Gaussian criterion

Let
\[
 P_G=\{x^2+y^2>0:\gcd(x,y)=1\}=C_4\sqcup2C_4,
\]
and count distinct ordered value vectors in
\[
 a+b+c=N,\qquad c+d+e=M,\qquad a,b,c,d,e\in P_G.            \tag{C.6}
\]
The value characterization follows from Gaussian unique factorization.
The Gaussian norm is Euclidean (rounding a quotient to a nearest lattice
point leaves remainder norm at most half the divisor norm).
Primitivity excludes every inert prime and permits at most one factor
of \(1+i\); conversely, choose just one prime above each split prime
to realize every value of \(C_4\) or \(2C_4\) primitively.

One primitive coordinate pair for \(c\) is shared by both equations.
The coefficient matrix is surjective; its columns for \(a,d\) form an
identity minor. With free parameters \(b,c,d\), the coordinate forms
\[
                  N-b-c,\ b,\ c,\ d,\ M-c-d
\]
have pairwise independent nonconstant parts.

Put \(S=\{1,2,5\}\subset\mathbb Z/8\), and define
\[
 R_8(r,s)=\#\{(a,b,c,d,e)\in S^5:
                    a+b+c=r,\ c+d+e=s\bmod8\}.            \tag{C.7}
\]

**Theorem CoupledPrimitiveGaussian.** Fix \(C\ge1\).
For all sufficiently large \(T\), and all integers \(N,M\in[T,CT]\),
system (C.6) is solvable if and only if
\[
 \begin{gathered}
 N,M\not\equiv2\pmod8,\\
 (N,M)\bmod8\notin\{(6,3),(6,7),(3,6),(7,6)\}.             \tag{C.8}
 \end{gathered}
\]
These 45 residue pairs are exactly the joint finite-congruence
conditions, and the same criterion holds for locally primitive
Gaussian coordinate solutions at every prime.
On these classes there are \(\gg_C T^3/(\log T)^{5/2}\)
distinct value vectors.

For fixed positive integers \(u,v\), put \(N=uT,M=vT\). On admissible
targets the full value count has asymptotic
\[
 R_G(N,M)\sim
 2R_8(N,M)c_4^5\,{\cal V}(N,M)
    \frac{\prod_{p=3\bmod4}\sigma_p(N,M)}{(\log T)^{5/2}},
                                                               \tag{C.9}
\]
where
\[
 \begin{split}
 m&=\min(N,M),\\
 {\cal V}(N,M)&=NMm-\tfrac12(N+M)m^2+\tfrac13m^3,\\
 \sigma_p(N,M)&=
 \frac{p^2}{(p-1)^5}
 \left[(p-1)(p-2)^2+(p-2)
       \bigl(\mathbf1_{p\nmid N}+\mathbf1_{p\nmid M}\bigr)
       +\mathbf1_{N=M\ne0\bmod p}\right].
 \end{split}                                                \tag{C.10}
\]
The prime product and its inverse are uniformly bounded.
The count is zero on excluded classes. No numerical onset or
logarithmic error rate is asserted.

**Proof.** A primitive Gaussian pair has norm \(1,2\) or \(5\bmod8\).
The pair-sum multiplicities for \(S\), in residue order \(0,\ldots,7\),
are
\[
                         q=(0,0,2,2,1,0,2,2).
\]
Hence \(R_8(r,s)=\sum_{c\in S}q(r-c)q(s-c)\).
The possible shared residues for a triple target \(r\) are
\[
\begin{array}{c|cccccccc}
r&0&1&2&3&4&5&6&7\\ \hline
c&S&\{2,5\}&\varnothing&\{1,5\}&S&\{1,2\}&\{2\}&\{1,5\}.
\end{array}
\]
Their intersections are nonempty exactly in (C.8), giving 45 pairs.

Every root modulo eight lifts to all higher powers using the original
matrix's integer right inverse; preserving the values modulo eight
preserves their allowed 2-adic value profiles.
At an odd inert prime every value must be a unit.
For any chosen shared unit, each private pair can sum to any prescribed
residue because \(p\ge3\); thus no odd-prime obstruction remains.
At a split odd prime there is no value restriction.
The finite-image lemma applied to the two profiles realizes these
residues by actual primitive norm values.
This proves the all-moduli criterion. Necessity for locally primitive
coordinate pairs follows from the same elementary mod-eight and
inert-prime calculations; sufficiency follows from the actual norm-value
representatives. At each prime the primitive solution sets modulo
successive powers form nested nonempty closed cylinders in the compact
space of local coordinates, so they supply compatible p-adic solutions. The finite-profile bulk theorem proves eventual
existence and the lower count, since the real cone is the whole
positive \((N,M)\)-quadrant.

For the full fixed-ratio asymptotic use the classical moving-target
theorem on the fixed positive polytope
\(b,c,d>0,\ b+c<N,\ c+d<M\), after scaling by \(T\).
Its volume is
\(\int_0^{\min(N,M)}(N-c)(M-c)\,dc={\cal V}(N,M)\).
The bounded multiplicative indicator of \(P_G\) is Frobenian on primes;
its mean exponent is \(1/2\).
The two exact prefixes at two must be retained. On every reduced
progression modulo a multiple of eight its mean is the odd-core mean,
with coefficient \(c_4\); the prefix two restores the even branch.
Each physical state \(1,2,5\bmod8\) has progression density \(4c_4\).
The finite coefficient is therefore
\((4c_4)^5 R_8/8^3=2R_8c_4^5\).
One must not assume that twisting the full indicator by \(\chi_{-4}\)
fixes its even branch.

For an inert odd prime, count unit solutions by conditioning on the
shared \(c\). A private unit pair has \(p-2+\mathbf1_{c=N}\)
possibilities, and likewise for \(M\).
Summing over \(c\ne0\) gives the bracket in (C.10).
Divide by \(p^3(1-1/p)^5\) to obtain \(\sigma_p\).
At a split prime the local factor is one.
The bracket is at least \((p-1)(p-2)^2>0\), and expansion gives
\(\sigma_p=1+O(p^{-2})\) uniformly in \(N,M\), including their diagonal.
This proves the product claims and (C.9).
\(\square\)

The full asymptotic (C.9) is stated along fixed target rays.
The uniform comparable-target existence and lower bound use finitely
many fixed positive boxes. No uniform moving-boundary asymptotic is
silently inferred.

**Proposition CoupledMarginalsFail.** There are infinitely many comparable
targets for which both equations in (C.6) are individually represented
by primitive Gaussian norms but the coupled system has no local
solution at two.

**Proof.** Put \(p=17^t\), \(N=2p+4\), \(M=p+2\), \(t\ge0\).
Both \(p\in C_4\) and \(2p\in2C_4\) are primitive Gaussian norm values.
The marginal witnesses are
\[
                       N=2+2+2p,\qquad M=1+1+p.
\]
But \(N=6\bmod8\) forces the shared norm to be \(2\bmod8\), whereas
\(M=3\bmod8\) forces it to be \(1\) or \(5\bmod8\).
This is impossible. The targets grow with \(N=2M\).
\(\square\)

The same shared-triple geometry also applies to \(C_{12}\).
Its criterion is simply \(N=M=3\bmod12\), with finite coefficient
\(144c_{12}^5\), exponent \(15/4\), the volume above and the product
(C.10) over \(p\nmid6,\ p\not\equiv1\bmod12\).
Every such value has both a primitive Gaussian and a primitive
Eisenstein realization: each prime \(1\bmod12\) splits in both Euclidean
rings, and choosing one prime above each rational prime keeps both
coordinate pairs primitive. Choosing one of each per value yields a
coupled system of seven quadratic equations in twenty coordinates,
with all ten coordinate pairs primitive. This conclusion counts
value vectors; it is not an asymptotic for coordinate multiplicities.

## 6. Coefficient primes and the actual analytic boundary

**Proposition CoupledCoefficientObstruction.** For \(C_{12}\) variables,
the system
\[
 x_1+5x_2-5x_3=N,\qquad 2x_2-x_3+x_4=2N                \tag{C.11}
\]
has eventual representability in fixed positive bulk regions exactly
when \(N=1\bmod12\) and \(5\nmid N\).

**Proof.** Parametrize by
\((N+5t_1,t_2,t_1+t_2,2N+t_1-t_2)\).
The rows are pairwise independent. The matrix is surjective, and the
only exceptional good prime in Section 2 is five.
The modulus-twelve condition is \(N=1\bmod12\).
If \(5\mid N\), its first variable is forced to be a multiple of five.
If \(5\nmid N\), the other three forbidden lines modulo five leave
13 of 25 parameter points, giving \(\sigma_5=325/256>0\).
All other good primes exceed four and have proper hyperplanes.
There is a fixed positive parameter box with \(t_1\) near \(N/4\)
and \(t_2\) near \(N/2\). Apply (C.4).
\(\square\)

Thus the infinitely many \(N=25\bmod60\) pass the original splitting
modulus but fail at a coefficient prime. Testing only \(L\) is not a
general coupled local criterion.

At the boundary, consider \(a+b=N,\ b+c=M\) in \(C_{12}\).
There is one free parameter and three dependent nonconstant forms
\(b,N-b,M-b\). The exact finite-admissibility criterion is
\(N=M=2\bmod12\): at a forbidden good prime at least five, avoid the
at most three residues \(0,N,M\).
Nevertheless \(N=26,M=38\) has no actual common choice of \(b\):
the first target forces \(b=13\), while the second permits only
\(b=1,37\). This is a small exception, not an eventual counterexample.

The diagonal \(M=N\) reduces exactly to the binary count and its
logarithmic exponent \(3/2\); treating it as three independent
marginals would predict the wrong exponent \(9/4\).
Off the diagonal, the current finite-complexity theorem does not
apply, and the binary theorem alone does not estimate the three-form
product. In the reflected graph proof the extra factor depends on
the edge-label product as well as the endpoint, so the existing
separable-endpoint comparison cannot simply absorb it.
A genuine joint correlation estimate is an additional obligation.
No impossibility theorem or historical unsolved-status claim is made.

## 7. Classical antecedents and scope

Matthiesen's *Linear correlations amongst numbers represented by
positive definite binary quadratic forms*, arXiv:1106.4690v3,
Theorem 1.3, already explicitly treats coupled quadratic equations
through an integer-kernel parametrization with the same row-space
support condition. There “primitive quadratic form” concerns the
gcd of its coefficients; it does not impose separate coordinate
primitivity. Its weighted coordinate correlations differ from the
unweighted primitive-value counts here.

Browning–Matthiesen, *Norm forms for arbitrary number fields as products
of linear polynomials*, Ann. Sci. ENS 50 (2017), Theorem 1.3,
proves rational Hasse-principle and weak-approximation results for
simultaneous norm equations. Its integral counting Theorem 5.2 and
Remark 5.3, printed p. 1401, allow moving \(O(T)\) affine constants;
Proposition 5.5 treats good local factors.
These are close classical antecedents, not claims supplied by the new
binary graph method.

The general coupled count in this note is also a classical consequence.
Its contributions within this research campaign are the explicit
joint finite criterion, correct image-lattice/profile handling,
the 45-class primitive-Gaussian application, and precise distinctions
between value counts, coordinate weights and degenerate correlations.
General nonlinear coupled varieties, arbitrary class-group norm-value
models, unbounded conductor profiles, effective small exceptions and
dependent-coordinate systems are not settled by this theorem.

General Euler indicators and nonsubgroup prime prescriptions are treated
in the [Hardy–Littlewood extension](MULTIPLICATIVE_HARDY_LITTLEWOOD.md),
including the explicit real-character and Gamma normalization conventions.
