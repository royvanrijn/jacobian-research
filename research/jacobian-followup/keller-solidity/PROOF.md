# Solidity and global regularization obstructions

14 November 2026. The nonzero-functional route admits a characteristic-free algebraic formulation and a short proof. The underlying finiteness theorem is classical: Datta--Smith, *Excellence in prime characteristic*, Proposition 3.7 with Lemma 3.6, [arXiv 1704.03628v3](https://arxiv.org/pdf/1704.03628v3), printed pages 8--10. Despite the paper's title, these statements do not assume positive characteristic. [Datta--Murayama, with Smith's appendix](https://arxiv.org/html/2007.10383v1), Proposition 6.1.2, explicitly records that implication and attribution.

The pairing proof below is self-contained. The finite-target and differential-operator consequences are elementary deductions recorded here for reuse; no novelty priority is claimed.

## A short algebraic finiteness proof

**FINSOL1.** Let A be a Noetherian domain, and let A be a subring of a domain R such that L=Frac(R) is a finite extension of K=Frac(A). Then

\[
 \operatorname{Hom}_A(R,A)\ne0
 \quad\Longleftrightarrow\quad
 R\text{ is a finite }A\text{-module}.
 \tag{1}
\]

No finite generation or Noetherian hypothesis on R is needed in this formulation.

First, S^{-1}R=L for S=A minus zero. Indeed S^{-1}R is a K-subalgebra of the finite-dimensional field L; every nonzero element has an inverse by its algebraic equation. Choose a K-basis b_1,...,b_n of L with each b_i in R.

If phi:R->A is nonzero and A-linear, its localization lambda:L->K is nonzero. The bilinear form lambda(xy) is nondegenerate: for x!=0 choose z with lambda(z)!=0 and take y=z/x. Therefore

\[
 r\longmapsto(\phi(b_1r),\ldots,\phi(b_nr))
 \quad\text{embeds }R\text{ in }A^n.
 \tag{2}
\]

Noetherianity of A makes this submodule finite. Conversely, if R is finite, clear one common denominator in a nonzero K-linear functional on finitely many A-module generators. This gives a nonzero A-linear R->A. This proves(1).

The proof uses no trace, separability, derivatives, normalization, grading, analytic continuation or affine-space recognition.

### An explicit common denominator

Let M=(phi(b_i*b_j)) and Delta=det(M), a nonzero element of A. For r=sum c_j*b_j write v_i=phi(b_i*r). Then Mc=v, so

\[
 \Delta(c_1,\ldots,c_n)^T=\operatorname{adj}(M)v\in A^n,
 \qquad
 \boxed{\Delta R\subseteq\sum_i A b_i.}
 \tag{3}
\]

One functional valid on all of R therefore supplies a single denominator for the entire source. Finite moment agreement does not certify that universal premise.

## Finite targets and differential operators

**FINSOL2.** Keep A,R as in FINSOL1 and let M be any nonzero finite torsion-free A-module. Then

\[
 \operatorname{Hom}_A(R,M)\ne0
 \quad\Longleftrightarrow\quad R\text{ is finite over }A.
 \tag{4}
\]

To prove the forward direction, embed M in A^q by choosing a basis of M tensor_A K and clearing the coordinate denominators of finitely many generators. Torsion-freeness makes this embedding injective. A nonzero map R->M has a nonzero coordinate R->A, so(1) applies. Conversely, compose any nonzero R->A from(1) with multiplication by a nonzero element of M.

This rules out bypassing a vanishing scalar dual by replacing the target with a finite free module, a vector bundle's finite module of sections, or another finite torsion-free coefficient module.

Now suppose A is a k-algebra. Regard R and M as A-modules. A k-linear differential operator relative to A has order zero when it is A-linear; an operator of order at most d has commutators

\[
 [D,a](r)=D(ar)-aD(r),\qquad a\in A,
\]

of order at most d-1. This is the [standard module definition](https://stacks.math.columbia.edu/tag/09CI).

For any A-modules N,M, the implication

\[
 \operatorname{Hom}_A(N,M)=0
 \ \Longrightarrow\
 \operatorname{Diff}^{\,d}_{A/k}(N,M)=0
 \quad\text{for every finite }d
 \tag{5}
\]

follows by induction. The order-zero case is the assumption. At the next order every commutator vanishes by induction, so the operator is A-linear and hence zero. No characteristic assumption enters this argument.

Combining(4) and(5), if R is not finite over A, then

\[
 \boxed{\operatorname{Diff}^{\,d}_{A/k}(R,M)=0}
 \quad\text{for every finite }d
 \text{ and every finite torsion-free }M.
 \tag{6}
\]

For each nonzero such M, existence of a nonzero finite-order operator is equivalent to finiteness, since an A-linear map is itself of order zero.

The module convention matters. A need not make M an R-module. If D is an ordinary finite-order source differential operator R->L whose image happens to lie in A, restriction of its commutators to elements of A puts it in Diff_{A/k}(R,A), so(6) applies as well.

## Application to complex Keller maps

For a polynomial map F:C^n->C^n with nonzero constant Jacobian, take A=C[F_1,...,F_n] and R=C[x_1,...,x_n]. The Jacobian condition supplies algebraic independence and a finite fraction-field extension.

FINSOL1 makes the existence of a nonzero A-linear functional equivalent to finiteness. A finite Keller map is etale and gives a finite analytic covering of simply connected C^n. The connected source makes its degree one. Normality of A then gives R=A and a polynomial inverse.

Thus the public campaign's plane criterion is a special case in every dimension. Only the final finite-etale-cover step uses complex affine space. This does not assert the existence of a functional, and any nonfinite higher-dimensional Keller examples must have zero global dual.

In characteristic p, keep the finiteness conclusion but do not replace it by invertibility. The finite etale Artin--Schreier map t->t^p-t has degree p and is not an automorphism. A purely inseparable extension can have identically zero field trace and still admit a nonzero functional; the trace-free proof covers it.

## A uniform obstruction for normal source completions

**WRIGHTSOL1, general form.** Let B be a Noetherian normal domain and let

\[
 B\subsetneq R\subset\operatorname{Frac}(B)
\]

be any proper intermediate domain. Then R cannot be finite over B: finite would imply integral, and normality would force R=B. Consequently

\[
 \operatorname{Hom}_B(R,M)=0,\qquad
 \operatorname{Diff}^{\,d}_{B/k}(R,M)=0
 \tag{7}
\]

for every finite torsion-free B-module M and finite order d. This requires neither a proposed target map nor any ramification index or degree bound. The normal-birational step is also [Stacks0AB1](https://stacks.math.columbia.edu/tag/0AB1).

For the actual [marked Wright model](../affine-complements/WRIGHT_MARKED_MODEL.md),

\[
 B=\mathbb C[Y,XY,X^2Y,X^3Y+X]
 \subset R=\mathbb C[X,Y].
\]

This B is the coordinate ring of the smooth marked affine normalization. Its fraction field is C(X,Y), since X=(XY)/Y. In the omitted-prime chart,

\[
 X=u^{-1},\quad Y=u^3v-u^2,\quad
 XY=u^2v-u,\quad X^2Y=uv-1,\quad X^3Y+X=v.
 \tag{8}
\]

All four generators are regular at u=0 but X has a pole. Hence X is not in B and the inclusion is proper. Equation(7) applies in all degrees to this fixed source inclusion.

In particular no nonzero ordinary finite-order differential operator on C[X,Y], even one with rational coefficients interpreted as an operator into C(X,Y), can send every source polynomial into B. Its restriction along B would contradict(7). Every finite matrix of such global operators is also zero entrywise. This excludes fixed finite-order differential regularization as an automatic construction from the ring alone.

### The pole proof for scalar functionals

There is also a direct degree-free proof for the scalar dual. A B-linear phi:R->B extends to an L-linear map L->L, so is multiplication by c=phi(1). If phi!=0 then c is a nonzero element of B with finite omitted valuation h>=0. But

\[
 \operatorname{ord}_{u=0}(c X^m)=h-m<0
 \qquad(m>h),
\]

contradicting cX^m in B. This identifies the missing universal input precisely: no single multiplier clears arbitrarily high source poles.

The same argument applies whenever an element of a source ring has a pole along a valuation nonnegative on the normalization. It is not special to the Wright generators.

## Controls showing the exact boundary

1. **Finite generic degree is essential.** For R=A[z], evaluation at z=0 is nonzero and A-linear but R is not finite over A. The fraction-field extension is transcendental.
2. **A domain generic fibre is essential.** For A=k[t] and a=t, with R=A times A[1/a], projection to the first factor is nonzero and A-linear although the second factor is nonfinite. One functional need not detect every factor.
3. **Torsion-free targets are essential.** For the Wright inclusion, evaluation at the source origin is B-linear R->C, where C is the residue module of that point. It is nonzero; C is torsion over B.
4. **Finite order is essential.** Evaluation E at X=0 maps R into C[Y], a subring of B. It is not B-linear. For b=XY, repeated commutators give
   ad_b^d(E)(1)=(-b)^d for every d>=0. These are nonzero, so E has no finite order relative to B. Its formal Taylor expression is pointwise finite on each polynomial but has no uniform order on R.
5. **Normality is essential for the proper-birational conclusion.** C[t^2,t^3] is a proper subring of C[t], but the latter is finite and multiplication by t^2 maps it nontrivially into the cusp ring.
6. **The base must be Noetherian for(1).** Let K/k be an infinite field extension, A=k+tK[t] and R=K[t]. Multiplication by t maps R into A, but R is not finite over A: reduction modulo tR would make K finite over k. Here tR/(tR)^2 is infinite-dimensional over k, so A is not Noetherian.

The checker verifies exact representative identities for these controls. Their unlimited scope and the algebraic finiteness arguments are written proofs.

## What changes in the research programme

The derivative/trace proof of the nonzero-functional criterion can be replaced by the short solidity argument. Finite torsion-free targets and finite-order differential operators do not weaken the global requirement. At the Wright normalization, no nonzero map of those types exists without changing the source or the type of target.

This is an architectural obstruction, not a proof that the Wright ring contains no Keller pair. Deriving a nonzero functional from the additional Keller equations would still be a valid contradiction strategy; the ring and its boundary alone cannot supply it. Infinite-order or nonlocal constructions and torsion targets are outside(6), and need their own global polynomial-value and linearity proofs.

The full quadratic/higher-jet systems, field degree, finiteness and general compatible completion remain unresolved. No r7, JC2, HC4 or dimension-minimum problem is closed.

The published solidity theorem retains its original authorship. The deductions here have no formal verification, independent whole-proof implementation, external review or novelty-priority claim.

