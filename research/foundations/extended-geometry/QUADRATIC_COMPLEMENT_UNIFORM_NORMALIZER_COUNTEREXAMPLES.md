# Uniform quadratic complements: explicit non-power counterexamples

Let K be ANY field of characteristic zero, let q be ANY nondegenerate quadratic form of dimension N>=3, and put M=N-2. Choose an orthogonal splitting over K

\[
 q=A(x_1,x_2)+B(z_1,\ldots,z_M),\qquad
 A=a_1x_1^2+a_2x_2^2,\quad B=\sum_{i=1}^M b_i z_i^2,
\]

with every coefficient nonzero. Let X={q!=0}, C=K[x_1,x_2,z_1,...,z_M,q^{-1}], and let Lambda be its unique normalized SO(q) x G_m-invariant trace. Write

\[
 U=q,\quad S=B/q,\quad V=A^2/q^3,\quad UV=(1-S)^2,
\]
\[
 R(U,S)=S(U+2)\{2-S(U^2+2U+2)\},
\]
\[
 P=(1+U)(V-R)=\frac{(q+1)\{q-(q+1)^2 B\}^2}{q^3},\quad
 F=P^M(1+P),\quad H=q^M.                                      \tag{1}
\]

These are actual regular functions; only q is inverted. In particular no angular coordinate or endpoint is inverted. All statements below hold over the given K and after EVERY characteristic-zero extension, without a rational-point, positivity, isotropy or parity assumption.

## 1. The exact all-order statements

Call f pure if Lambda(f^m)=0 for every m>=1. Call it safe if, for every fixed regular marker g, Lambda(g f^m)=0 for all sufficiently large m, with cutoff depending on g. A pure phase with one marker surviving infinitely often is unsafe and disproves the Mathieu property of ker Lambda.

For k>=0 put the nonzero positive rational

\[
 J_{M,k}=M\int_0^1 w^{M-1}(1-w^2)^{2k}\,dw
       =\frac{(2k)!}{(1+M/2)_{2k}}.                            \tag{2}
\]

Every integral in this proof is finite polynomial integration over Q, extended K-linearly. No integration on K or ordering of K is involved.

**Theorem 1 (exact moments and sharp pure powers).** For every k>=1 and every integer ell>=-k,

\[
 \Lambda(q^\ell P^k)=J_{M,k}\binom{k-M}{k-\ell},               \tag{3}
\]

with generalized binomial coefficients given by formal expansion of (1+U)^n for integer n, and binom(n,r)=0 for r<0. For ell>k the trace is zero independently by scalar weight. Formula (3) is NOT asserted for ell<-k.

For every positive integer a,

\[
 P^a\text{ is pure}\quad\Longleftrightarrow\quad a\geq M.     \tag{4}
\]

If a>=M and ell is ANY fixed integer, then the marker q^ell survives at every sufficiently large order of P^a if and only if ell>=M. If ell<M it vanishes eventually. More precisely, for every m>=1,

\[
 \Lambda(q^\ell P^{am})\ne0\text{ at EVERY }m
 \quad\Longleftrightarrow\quad M\leq\ell\leq a.               \tag{5}
\]

Here the left side means that the same fixed ell works at all positive orders. In particular Lambda(q^M P^k)=J_(M,k)!=0 for every k>=M.

**Theorem 2 (uniform absolute non-power witness).** The SAME explicit F,H in (1) satisfy, for EVERY m>=1,

\[
 \Lambda(F^m)=0,\qquad
 \Lambda(HF^m)=\sum_{r=0}^m\binom mr J_{M,Mm+r}\ne0.           \tag{6}
\]

In C tensor_K Kbar, F is not u g^e for ANY unit u, regular g, integer e>=2. Thus it is absolutely non-power even up to units. Consequently ker Lambda is not Mathieu for ANY such K,q,N>=3.

**Theorem 3 (direct ternary phase).** For N=3, M=1, P itself is pure, absolutely non-power even up to units, and the SAME H=q has Lambda(qP^m)=J_(1,m) at EVERY m>=1.

The normalizer circuit in (1) and its unweighted endpoint identity are Long's known method. The rank-two/rank-M splitting, direct trace calculation on X and uniform dimension formula are proved below. No previously constructed counterexample family or finite-cover failure-descent assertion is a premise.

## 2. The normalized trace over the original field

The group G=SO(q) x G_m is smooth connected reductive. Over Kbar, SO(q) is transitive on each nonzero norm level: an orthogonal map sending two vectors of the same nonzero norm can have its determinant corrected by reflection in the nondegenerate orthogonal complement. Scalar multiplication then makes G geometrically transitive on X. This requires no K-point on {q=1}. It also gives C^G=K by scalar extension, since X is geometrically integral.

Every rational representation is locally finite; same-field complete reducibility for connected reductive groups in characteristic zero splits off the trivial isotypic part. Therefore

\[
 C=K1\oplus\operatorname{Lie}(G)C
\]

as K-vector spaces, and projection to K1 is the unique normalized invariant trace Lambda. In particular Lambda(Df)=0 for every Lie derivation D. The invariant projection and this decomposition commute with every scalar extension. The classical input is Milne, Algebraic Groups (2017), Theorem22.42/Corollary22.43 with Section12l's locally finite/Reynolds discussion. The already retained geometrically homogeneous-form trace lemma supplies an alternate route at exactly this scope; its all-automorphism theorem is not needed.

Scalar multiplication acts on q^r g(S) with weight 2r. Consequently, for every integer r and every polynomial g over K,

\[
 \Lambda(q^r g(S))=\begin{cases}\Lambda(g(S)),&r=0,\\0,&r\ne0.\end{cases} \tag{7}
\]

This is a direct torus-weight identity on X. No double cover or passage to a norm-one real locus occurs.

## 3. The rank-two beta recurrence on X

For alpha=1,2 and i=1,...,M let

\[
 D_{\alpha i}=b_i z_i\partial_{x_\alpha}-a_\alpha x_\alpha\partial_{z_i}.
\]

They annihilate q and belong to Lie SO(q). For each j>=1, direct differentiation in C gives

\[
 \sum_{\alpha,i}D_{\alpha i}\left(S^{j-1}x_\alpha z_i/q\right)
 =2S^j-(M+2j-2)(1-S)S^{j-1}.                                \tag{8}
\]

Indeed the derivative of x_alpha contributes 2S^j; the derivative of z_i contributes -M(1-S)S^(j-1); and the derivative of S^(j-1) contributes -2(j-1)(1-S)S^(j-1). At j=1 the last contribution is zero and has no negative power of S.

Set mu_j=Lambda(S^j), mu_0=1. Applying Lambda to (8) yields

\[
 (M+2j)\mu_j=(M+2j-2)\mu_{j-1},\qquad
 \mu_j=\frac{(M/2)_j}{(M/2+1)_j}=\frac M{M+2j}.               \tag{9}
\]

All denominators are nonzero in characteristic zero. Thus, for EVERY polynomial g in K[S],

\[
 \Lambda(g(S))=M\int_0^1 t^{M-1}g(t^2)\,dt.                 \tag{10}
\]

Both (8) and (9) hold directly over K. The polynomial density t^(M-1) works for EVERY M>=1, so there is no parity restriction.

## 4. Regular circuit and Laurent coefficient extraction

The exact relation UV=(1-S)^2 gives

\[
 U(V-R)=\{1-(1+U)^2S\}^2,\qquad
 P=U^{-1}(1+U)\{1-(1+U)^2S\}^2.                            \tag{11}
\]

Substituting S=B/q gives precisely (1). Combining the scalar-weight identity (7) and height trace (10), for every k>=1 and integer ell,

\[
 \Lambda(U^\ell P^k)
 =M[U^{k-\ell}](1+U)^k\int_0^1
             t^{M-1}\{1-(1+U)^2t^2\}^{2k}\,dt.              \tag{12}
\]

The expression before coefficient extraction is a polynomial in U. All exchanges of integrals and coefficients are finite. No independence of the two formal variables U,S in C is required: this is evaluation of the finite displayed expression by (7) and (10).

## 5. Endpoint divisibility, exact formula and both sharp converses

For every k>=1 define

\[
 J_{M,k}(T)=M\int_0^T w^{M-1}(1-w^2)^{2k}\,dw.
\]

It is a rational polynomial divisible by T^M. Its derivative is divisible by (T-1)^(2k). Hence

\[
 J_{M,k}(T)-J_{M,k}(1)\in(T-1)^{2k+1}\mathbb Q[T].             \tag{13}
\]

Set T=1+U in (12), and substitute w=Tt in each monomial. It becomes

\[
 \Lambda(U^\ell P^k)=[U^{k-\ell}]T^{k-M}J_{M,k}(T).           \tag{14}
\]

Even if k<M, the displayed product is polynomial because T^M divides J. For endpoint comparison T^(k-M) is a unit power series at U=0. If ell>=-k then k-ell<=2k, so (13) permits replacing J(T) by J(1) at the required coefficient. This proves (3), including the upper-negative generalized binomial cases. If ell>k the required coefficient index is negative and the original polynomial has no such term. For ell<-k, endpoint error may contribute; no formula (3) is claimed there.

For k>=M, the upper binomial index k-M is nonnegative and less than k. Thus Lambda(P^k)=0 and Lambda(q^M P^k)=J_(M,k). Conversely, if 1<=k<M,

\[
 \Lambda(P^k)=(-1)^k\binom{M-1}{k}J_{M,k}\ne0.                \tag{15}
\]

If a<M, the first moment of P^a is therefore nonzero, so it is not pure. If a>=M, every positive power has exponent am>=M and vanishing trace. This proves BOTH directions of (4).

For fixed ell and large k=am, we have ell>=-k and ell<=k. The nonnegative upper binomial in (3) is nonzero precisely when 0<=k-ell<=k-M, that is ell>=M. This proves the eventual marker iff for ALL fixed integer ell, including negative ones. If M<=ell<=a, the same inequalities hold for every m>=1, giving (5)'s sufficient direction. If ell>a, its first marked moment is zero since ell>k=a; if ell<M it vanishes eventually. These give the necessary direction of (5).

Lastly (2) follows by finite beta integration, or by repeated integration by parts, and is a strictly positive rational number. Expanding F^m=P^(Mm)(1+P)^m gives (6), since every exponent Mm+r>=M. Its marked sum is nonzero in Q and remains nonzero under Q->K. This imposes no ordering on K and uses the ONE marker H=q^M at EVERY positive order.

## 6. Absolute non-power, even up to units

Work in the UFD D=Kbar[x_1,x_2,z_1,...,z_M,q^{-1}]. Choose a point xi with all z_i=0 and A=-1/2. Such a point exists over Kbar even when rank-two A is geometrically reducible. Nondegeneracy implies dA is nonzero there, since zero gradient would force x_1=x_2=0 and A=0. Also q=-1/2, so xi belongs to X.

On B=0, P=(q+1)/q=1+1/q. At xi, P=-1 and dB=0. Thus

\[
 d(1+P)|_\xi=-q^{-2}dq|_\xi=-4dA|_\xi\ne0.                 \tag{16}
\]

Consequently 1+P lies in the maximal ideal of the regular local ring at xi but not its square. In the UFD factorization of 1+P, exactly one prime factor through xi has multiplicity one: two vanishing factors or multiplicity>=2 would put it in the square. Since P(xi)=-1, this prime does not divide P. Its valuation on F=P^M(1+P) is exactly one.

Every unit has valuation zero, and the valuation of a unit times an e-th power is divisible by e>=2. Hence F is not u g^e over Kbar. No global irreducibility of 1+P or of rank-two A is asserted or needed. This proves Theorem2's full absolute non-power conclusion.

## 7. Direct ternary phase and scope boundary

When N=3, M=1. Theorem1 makes P itself pure with Lambda(qP^m)=J_(1,m) for every m>=1. Over Kbar, q+1 is irreducible: its homogenization q+T^2 has rank four, whereas a reducible quadratic has rank at most two. It remains a nonunit prime after q is inverted. In (11), the squared bracket is one modulo q+1; therefore ord_(q+1)(P)=1. The direct P is thus absolutely non-power even up to units. This proves Theorem3, including the positive-definite real norm-three complement.

The same construction proves failure for ALL nondegenerate q,N>=3 in characteristic zero; the present theorem makes no claim in dimensions one or two, on degenerate quadrics, or in positive characteristic. It is independent of earlier even-dimensional and isotropic constructions. It makes no classification claim about finite-cover descent: an arbitrary upstairs pure unsafe phase need not descend, and no such implication is used. No radial-primitivity, Image-Conjecture, formal-verification or worldwide-priority claim is made.

## 8. Exact primary attribution and replay boundary

The circuit R and P and its unweighted endpoint method are credited to Christopher D. Long, A Classification of Mathieu's Second Conjecture for Reductive Spherical Pairs, research draft dated22September2026, exact GitHub commit e9edac262b4b64eb68feb0e6e80496cefe4a1dee. The retained AUTHOR PDF has SHA256 b2bf25fbc2d6211333d4ff3e1f6de0b6562990e869c6078a2e219ef746f1026b (Git blob6d82f192462c8f342ee8577a64dfca1c84db96ea). The inspected scope is Section5, Proposition5.1, equations5.1--5.3 and the COMPLETE endpoint proof, printed pages5--6. Its versioned LaTex and author PDF were compared at those statements.

https://github.com/long-mathematics/mathieu-second-conjecture-spherical-pairs/blob/e9edac262b4b64eb68feb0e6e80496cefe4a1dee/mathieu_second_conjecture_spherical_pairs.tex

Long's statement concerns the complex normalizer pair SL2/N(T). This proof instead establishes the radial substitution, rank-two/rank-M moments over the given K, uniform endpoint formula, both sharp converses and the simple-zero consequence. All required circuit and endpoint steps are reproduced, so Long is credited classical method, not a necessary logical premise. No broader classification theorem or previous application of the paper is imported or certified.

Milne is reused only for same-field rational complete reducibility and local finiteness/Reynolds projection. The existing normalized homogeneous-form trace note offers an alternate narrow route to the normalized trace. Neither its all-automorphism principal nor its finite-etale machinery is a core premise here.

The compact exact checker authenticates the written proof, typed scope and retained sources; compares direct Laurent expansion to exact endpoint formulas; checks both sides of the sharp finite models, factor identity and simple-point values; and rejects a wrong phase. Bounded models do not prove the universal fields, forms, dimensions, orders, converse statements or divisor multiplicities. The written proof does. Every external-review, formal-verification, independent-replay and worldwide-priority flag remains false.
