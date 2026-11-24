# Explicit witnesses on seven standard compact classical products over characteristic-zero fields

The construction uses Christopher D. Long's interval–torus witness and normalizer circuit. Abstract failure of tensor-product stability already follows from his interval example, real positivity, and Laurent constant-term theory; complex existence on the complexifications below is already covered by his classification of reductive spherical pairs. The result established here gives explicit **real regular coordinates**, prescribed **rational models**, their witnesses over **every characteristic-zero field**, absolute nonpower witnesses, and genuine two-factor Mathieu instability over **every formally real field**. Section 15 gives exact source roles and pins. No worldwide-priority assertion is made.

## 1. Statement, conventions, and formulas

Let K be any field of characteristic zero. If X_Q is one of the seven rational affine homogeneous models defined in Section 2, put X_K=X_Q×_Q K, B_K=O(X_K), and let Lambda_K be its normalized invariant linear functional. On B_K[t,t^-1] put

    lambda_K = Lambda_K tensor CT_t.

Here CT_t is the coefficient of t^0. A phase f is **pure** if lambda_K(f^m)=0 for every integer m>=1. It is **safe** if, for every fixed h in the same coordinate ring, lambda_K(h f^m)=0 for all sufficiently large m. The kernel of lambda_K is a Mathieu subspace precisely when every pure phase is safe. The stronger witness constructed here has ONE fixed regular marker H with lambda_K(H f^m) nonzero at EVERY positive m. An **absolute unit-nonpower** is a function which, over an algebraic closure of K, cannot be written u g^e with u a unit and e>=2 an integer.

The seven families and their nontrivial ranges are:

| Family | Standard compact real orbit | Range |
| --- | --- | --- |
| BDI | rank-k real orthogonal projector orbit | n>=3, 1<=k<=n-1 |
| AIII | rank-k complex Hermitian projector orbit | n>=2, 1<=k<=n-1 |
| CII | rank-k quaternionic Hermitian projector orbit | n>=2, 1<=k<=n-1 |
| AI | SU(n)/SO(n) | n>=2 |
| AII | SU(2n)/Sp(n) | n>=2 |
| CI | Sp(n)/U(n) | n>=1 |
| DIII | chosen SO(2n)/U(n) complex-structure orbit | n>=2 |

Sp(n) denotes the compact quaternionic unitary group. The projector varieties are their **affine algebraic homogeneous models**, not complete projective Grassmannian schemes.

For the normalizer construction write

    p(S,t) = (1+t)t^-1 [1-(1+t)^2 S]^2,
    F_a(S,t) = p(S,t)^a (1+p(S,t)).                (1.1)

For the polynomial-CDF construction write

    v_(alpha,beta)(D)
      = [integral_0^D u^(alpha-1)(1-u)^(beta-1)du]
         /[integral_0^1 u^(alpha-1)(1-u)^(beta-1)du],
    E(v,t) = (1-t^-1)(1-v+vt),                    (1.2)

where alpha,beta are positive integers; consequently v_(alpha,beta) is a polynomial with rational coefficients. In all formulas an integral of a rational polynomial means its exact rational antiderivative and evaluation, so its value makes sense in every characteristic-zero K.

The concrete phase and marker choices are:

| Family | Regular height or CDF input | Phase f | Marker H |
| --- | --- | --- | --- |
| BDI | Q=P,r=k if k>=2; Q=I-P,r=n-1 if k=1; S=Q11 Q22-Q12 Q21 | F_(n-2)(S,t) | t^(r-1) |
| AIII | D=P11, v=v_(k,n-k)(D) | E(v,t) | t^-1 |
| CII | D=P11, v=v_(2k,2(n-k))(D) | E(v,t) | t^-1 |
| AI | C=gg^T, S=1-C11 bar(C11) | F_(n-1)(S,t) | t^(n-1) |
| AII | C=gJ0g^T, D=C12 bar(C12), v=1-(1-D)^(2n-2) | E(v,t) | t^-1 |
| CI | C=g(iI)g*, x=coefficient of i in C11, S=x^2 | F_(2n-1)(S,t) | t |
| DIII | C=gJ0g^T in the chosen orbit, x=C21, S=x^2 | F_(2n-3)(S,t) | t |

**Theorem 1 (seven standard compact classical products).** Every f,H is regular over Q on the indicated affine quotient times G_m. For EVERY char0 K and EVERY m>=1,

    lambda_K(f^m)=0,        lambda_K(H f^m) != 0.   (1.3)

Every f is an absolute unit-nonpower. Thus all these product traces fail the Mathieu property, with a unit marker H and all-order detection. For every formally real K, BOTH individual factor traces on B_K and K[t,t^-1] have the Mathieu property; these same formulas therefore give genuine two-factor instability over K, including Q and R.

The word “standard” has exact content: all seven all-field assertions refer to X_Q base changed to K, with its specified involutions and coefficient coordinates. There is one separately proved extension: the BDI construction works for the self-adjoint projector variety of ANY nondegenerate quadratic form over ANY characteristic-zero K of dimension n>=3. Formal-real factor positivity is asserted for the standard rational models, not for all such twisted forms.

## 2. The seven rational affine models and their geometry

Set E=Q[i]=Q[T]/(T^2+1), with Q-linear conjugation i->-i. Let D_H=Q+Qi+Qj+Qk be the Hamilton quaternion algebra, i^2=j^2=-1, ij=-ji=k, with canonical conjugation. A star means conjugate transpose in the relevant algebra. “Real part”, “imaginary part”, and “coefficient of i” below mean Q-linear coordinate projections in these fixed bases. They remain algebraic operations after base change; they do not require K to contain a scalar square root of -1.

### 2.1. BDI: real self-adjoint projectors

Take G_Q=SO_n for q0=sum x_i^2. Let X_Q be the closed scheme

    P^T=P, P^2=P, tr P=k.

Over an algebraically closed characteristic-zero field, an idempotent is diagonalizable with eigenvalues 0,1, so its trace fixes its rank k. Self-adjointness makes image and kernel orthogonal and nondegenerate. Conversely a nondegenerate k-plane defines its unique orthogonal projector. An orthogonal change of basis is transitive. A determinant-minus-one transformation in a nonzero image or kernel changes the determinant of a conjugating matrix, so SO_n itself is transitive. The stabilizer is S(O_k×O_(n-k)), reductive and possibly disconnected.

There is no extra infinitesimal structure in the displayed equations. At diag(I_k,0), differentiating P^2=P kills the two diagonal tangent blocks, and the off-diagonal blocks are arbitrary subject to transpose symmetry. Thus the tangent dimension is k(n-k), the orbit dimension. Conjugation gives the same calculation at every geometric point. Hence the closed scheme is smooth, reduced, and geometrically integral, with the stated affine homogeneous structure.

### 2.2. AIII: E-Hermitian projectors

Take G_Q=SU_n(E,identity Hermitian form). Define X_Q by

    P*=P, P^2=P, sum_i Pii=k,

in Mat_n(E). A diagonal entry of a Hermitian matrix belongs to the fixed scalar line Q, so the last equation is a scalar equation. G acts by conjugation, with stabilizer S(U_k×U_(n-k)) at diag(I_k,0).

After algebraic closure E splits as the product of two copies of the field, with conjugation interchanging the factors. A Hermitian matrix is then a pair (A,A^T), and the equations reduce to A^2=A,tr A=k. GL_n is transitive on these rank-k idempotents, and a determinant adjustment inside the block-diagonal centralizer makes SL_n transitive. The geometric stabilizer is S(GL_k×GL_(n-k)), reductive. At the standard projector the two off-diagonal tangent blocks are arbitrary, giving dimension 2k(n-k), equal to the orbit dimension. This proves smoothness and geometric integrality of the closed affine model. In particular P11 is a Q-regular scalar function.

### 2.3. CII: Hamilton-Hermitian projectors

Take G_Q={g in Mat_n(D_H):gg*=I}, the quaternionic unitary Q-group. Let X_Q be

    P*=P, P^2=P, sum_i Pii=k.

Each diagonal entry is again a scalar, since the fixed subspace of quaternion conjugation is the scalar line. The diagonal scalar trace is half the reduced matrix trace, so after a quaternion splitting extension it fixes the projector's ordinary rank to 2k. The stabilizer is Sp(k)×Sp(n-k) in quaternionic notation.

Over an algebraic closure D_H becomes Mat_2, and its conjugation becomes symplectic adjoint. Thus X becomes the symplectic-self-adjoint idempotents of rank 2k in dimension 2n. Their images and kernels are nondegenerate symplectic subspaces. Symplectic basis completion proves transitivity of Sp_(2n), with stabilizer Sp_(2k)×Sp_(2(n-k)). At the standard projector one tangent off-diagonal block is arbitrary and the other is its adjoint; its dimension is 4k(n-k), equal to the orbit dimension. This again gives a smooth geometrically integral affine scheme. The scalar diagonal coordinate P11 is Q-regular.

### 2.4. AI and AII Cartan matrix models

For AI let G_Q=SU_n(E) and H_Q=SO_n(q0), embedded by fixed-coefficient matrices. The map g->C=gg^T is right-H-invariant. Over an algebraic closure the quotient is the symmetric determinant-one matrices, SL_n/SO_n. Indeed every invertible symmetric matrix is gg^T; if det C=1 then det g is ±1, and an orthogonal reflection adjusts it to 1. This is a closed smooth affine hypersurface in the symmetric matrices, since the determinant differential is nonzero at every invertible matrix. It is one SL_n orbit and therefore geometrically integral.

For AII take G_Q=SU_(2n)(E), and fix the rational orthogonal skew matrix J0=diag(J2,...,J2), J2=[[0,-1],[1,0]]. Let H_Q be the subgroup preserving J0; its real compact group is Sp(n). The map C=gJ0g^T is right-H-invariant. Geometrically the quotient is the skew matrices of Pfaffian Pf(J0), an affine smooth SL_(2n) orbit with stabilizer Sp_(2n). Nonzero Pfaffian forces invertibility, and the Pfaffian differential is nonzero there. These descriptions also show that the rational Cartan maps define the entire affine homogeneous quotients, rather than only coordinate functions on an unrelated larger variety. Faithful base change verifies the corresponding identification over Q.

### 2.5. CI and DIII closed semisimple orbits

For CI use G_Q={gg*=I} over D_H and H_Q the centralizer of iI, namely the usual U_n(E). The Q-defined orbit map is C=g(iI)g*. After algebraic closure it is the closed semisimple orbit in sp_(2n) whose eigenvalues are +i,-i, each with multiplicity n. The symplectic form pairs the opposite eigenspaces nondegenerately, and its centralizer is GL_n. Symplectic basis completion proves transitivity. The two off-diagonal tangent blocks are symmetric, giving dimension n(n+1). Equivalently the Lie-algebra equation C^2=-I specifies that orbit geometrically; semisimplicity and the eigenspace decomposition show it is closed. The matrix coefficient of the quaternion basis element i in C11 is a Q-regular scalar coordinate.

For DIII let G_Q=SO_(2n)(q0), let J0 be the rational block matrix just used, and take its SPECIFIED conjugacy orbit X_Q with base point J0. Its stabilizer is the centralizer U_n(E). Concretely one may specify the chosen closed component by

    C^T=-C, C^2=-I, Pf(C)=Pf(J0).

Over an algebraic closure the ±i eigenspaces are maximal isotropic of dimension n; the centralizer in SO_(2n) is GL_n. Orthogonal basis completion, with the specified Pfaffian, gives one SO orbit. The equation without the Pfaffian qualification can contain the other orientation orbit, and is not our model. Tangent dimension n(n-1) agrees with the orbit dimension. This is a smooth geometrically integral affine semisimple orbit, with Q-regular coordinate C21.

Thus all seven X_Q are smooth geometrically integral affine homogeneous varieties for connected reductive Q-groups and reductive stabilizers. The real groups are the named compact classical groups. These descriptions are also valid when E or D_H splits after a later field extension. No compactness of K-points or classification of arbitrary K-forms is inferred.

## 3. Normalized trace, base change, and compact real positivity

Let B_Q=O(X_Q) for any model above. It is a locally finite rational G_Q-module. Rational representations of a characteristic-zero reductive group are completely reducible: on each finite-dimensional rational submodule there is a canonical projection to its trivial isotypic part. These projections are compatible on overlaps and define the Reynolds projection on B_Q. Geometric transitivity gives B_Q^G=Q, so the resulting Q-linear functional Lambda_Q is normalized by Lambda_Q(1)=1 and unique among invariant normalized linear functionals. Equivalently B_Q=Q direct-sum (Lie G_Q)B_Q. The same argument works over any characteristic-zero field.

These projections commute with field extension. One way to see this without choosing decompositions is that invariants are kernels of linear action equations and commute with field extension; the sum of the nontrivial isotypic summands has zero invariants after geometric scalar extension and the canonical projection remains the invariant projection. Therefore

    B_K = K tensor_Q B_Q,    Lambda_K = id_K tensor Lambda_Q.  (3.1)

The product trace is Lambda_K tensor CT_t, since G_m gives its character grading and the only trivial character is degree zero.

Over R, Lambda_R is normalized Haar integration on the compact group applied to the pullback of regular functions along G_R->X_R. Indeed averaging every finite-dimensional rational submodule is invariant under the compact group. Compact real points are Zariski dense in the real algebraic group: a smooth real neighborhood of the identity has full real dimension and is a totally real local chart in the complex group, so a regular function zero there is zero geometrically. Thus Haar averaging is algebraically invariant and equals the unique normalized trace. This argument uses the compact orbit through the chosen real base point, not a claim that every real point lies in that orbit.

For every nonzero f in B_R, its pullback to G_R is nonzero because the quotient map is dominant. Zariski density gives a compact real point where it is nonzero, and continuity gives an open set of positive Haar measure on which f^2 is strictly positive. Hence

    Lambda_R(f^2)>0 for every nonzero f in B_R.       (3.2)

The only pure phase on the single real factor is therefore zero. Section 13 upgrades this argument to all formally real fields by rational Gram matrices, without choosing an ordering or embedding the whole field in R.

The precise imported algebraic input here is characteristic-zero reductive complete reducibility and its normalized projection. It does not require the stronger automorphism-invariance or finite-cover conclusions of the separately retained reductive homogeneous trace theorem.

## 4. Polynomial square-root densities and the normalizer circuit

We prove the general all-order statement used in BDI, AI, CI, and DIII. Let X be smooth geometrically integral affine over a characteristic-zero K, with a normalized K-linear functional Lambda. Suppose S in O(X) satisfies the following identity of polynomial functionals:

    Lambda(g(S)) = integral_0^1 rho(z)g(z^2) dz        (4.1)

for EVERY polynomial g in K[z], where rho is a nonnegative rational polynomial on [0,1] of mass one. Write

    rho(z)=sum_(j=s-1)^b c_j z^j,
    c_(s-1)>0, c_b!=0, s>=1, a=b+1.

Individual higher coefficients c_j may have either sign. On X×G_m use p=p(S,t) from (1.1), and lambda=Lambda tensor CT_t.

### 4.1. Exact endpoint calculation and its domain

For integers k>=1 and ell put T=1+t and

    J_(j,k)(T)=integral_0^T w^j(1-w^2)^(2k)dw.

Taking Lambda coefficientwise, making the substitution w=Tz, and extracting the Laurent constant term gives

    lambda(t^ell p^k)
      = [t^(k-ell)] sum_j c_j T^(k-j-1)J_(j,k)(T).   (4.2)

The right side is interpreted in K[[t]], where T is a unit. This identity follows either as a polynomial identity after clearing powers of T or directly by substitution into the finite polynomial integral. Its left side is always a Laurent polynomial coefficient.

Since the derivative J'_(j,k)(T)=T^j(1-T^2)^(2k) vanishes to order 2k at T=1,

    J_(j,k)(1+t)-J_(j,k)(1) belongs t^(2k+1)Q[t].

Consequently, when ell>=-k, the coefficient in (4.2) is unchanged if every J_(j,k)(T) is replaced by J_(j,k)(1): the desired degree k-ell is at most 2k. This restriction matters for negative marker exponents; no unrestricted endpoint-replacement formula is asserted.

If k>=a, all exponents k-j-1 are nonnegative. For ell=0 each remaining polynomial has degree less than k, so its required coefficient is zero. For ell=s only j=s-1 has degree k-s, and that top coefficient is one. Thus for EVERY k>=a,

    lambda(p^k)=0,
    lambda(t^s p^k)=c_(s-1)J_(s-1,k)(1)>0.          (4.3)

The last expression is a strictly positive RATIONAL number. For F=p^a(1+p), expand F^m=sum_(h=0)^m binom(m,h)p^(am+h). Every exponent is at least a. Hence for EVERY m>=1,

    lambda(F^m)=0,
    lambda(t^s F^m)
      =c_(s-1)sum_(h=0)^m binom(m,h)J_(s-1,am+h)(1)>0.  (4.4)

No root-purity inference is involved: the seed is eventually pure, and the displayed polynomial in the seed is the pure phase.

### 4.2. What is sharp for general polynomial densities

For each fixed integer ell<s, choose k large enough that k>=a and k>=-ell. Every degree k-j-1 is at most k-s, strictly less than k-ell, so (4.2) is zero. It follows that t^ell is an eventually vanishing marker both for every pure power p^c with c>=a and for F; in contrast t^s detects every positive order of each. Thus s is the LEAST detecting integer exponent among the monomial markers t^ell in these cases. For a general multi-term rho this proof does not classify every ell>s, and a=b+1 is a sufficient cutoff, not a claimed least pure-power exponent.

### 4.3. The sharp monomial-density specialization

If rho(z)=s z^(s-1), then for every k>=1 and ell>=-k equation (4.2) gives

    lambda(t^ell p^k)
      =s J_(s-1,k)(1) binom(k-s,k-ell),             (4.5)

where binom(d,n)=[t^n](1+t)^d for n>=0, and is zero for n<0. Negative top integers use the formal binomial series. For a positive integer c, p^c is pure IFF c>=s. The forward direction follows from (4.3). If c<s, already its first pure moment is nonzero because

    binom(c-s,c)=(-1)^c binom(s-1,c)!=0.

For a pure power p^c with c>=s, a fixed integer monomial marker t^ell is unsafe IFF ell>=s: once k=cm is large, ell<s gives the degree zero above, whereas ell>=s gives a nonzero ordinary binomial coefficient binom(k-s,ell-s). Such a marker survives EVERY positive m IFF s<=ell<=c; if ell>c its first coefficient has negative lower index and is zero. These are genuine two-direction statements, proved from the exact identity rather than from finite probes.

### 4.4. Absolute nonpower

The seed p has order one along the prime divisor t+1: its bracket evaluates to 1 there. This does not say it is pure. If S vanishes at a geometric point x0 of X, set t=-1/2 at that point. Then p=-1 and

    partial_t(1+p)=-4.

In the smooth local ring of X_Kbar×G_m the function 1+p has a smooth hypersurface germ and a prime divisor of valuation one through this point. The factor p is a unit there. Thus F=p^a(1+p) has valuation one along that divisor and cannot be a unit times any e-th power with e>=2. No global irreducibility of 1+p is assumed or needed. This argument works over every algebraic closure of a characteristic-zero K.

## 5. Polynomial beta CDFs and the interval–torus witness

Let a regular D on an affine geometrically integral X satisfy, for positive INTEGERS alpha,beta,

    Lambda(g(D))=B(alpha,beta)^-1
         integral_0^1 u^(alpha-1)(1-u)^(beta-1)g(u)du   (5.1)

for every polynomial g. Define the polynomial v=v_(alpha,beta)(D) from (1.2). It has v(0)=0,v(1)=1 and derivative equal to the normalized beta density. The polynomial chain rule gives

    Lambda(v^j)= integral_0^1 v(u)^j v'(u)du
                =1/(j+1)                            (5.2)

for every j>=0. This is an identity of polynomial functionals; no analytic change of variables over K is required.

For E=(1-t^-1)(1-v+vt), direct polynomial integration in v gives, for m>=1,

    Lambda(E^m)
      =(t-1)^(m-1)(t^(m+1)-1)/[(m+1)t^m].           (5.3)

For example one integrates (1+(t-1)v)^m first and then multiplies by (1-t^-1)^m; the expression at t=1 is interpreted by its Laurent polynomial extension. Its constant coefficient is zero: the two summands have degree ranges [1,m] and [-m,-1]. Its coefficient of t is (-1)^(m-1)/(m+1). Therefore

    lambda(E^m)=0,
    lambda(t^-1 E^m)=(-1)^(m-1)/(m+1)!=0             (5.4)

at every positive m. The ideal (t-1) is geometrically prime because O(X_Kbar) is a domain. The factor 1-v+vt is 1 modulo t-1, so E has valuation one at this divisor. Thus E is an absolute unit-nonpower. We do not claim a least signed monomial marker exponent for this CDF witness; t^-1 is the single explicit marker needed here.

This is Long's interval–torus witness, with the uniform interval coordinate supplied by a regular polynomial CDF. Formula (5.2), rather than a finite distributional test, is the reason it applies.

## 6. BDI: a two-minor law from two Lie recurrences

Choose Q=P of rank r=k if k>=2, and Q=I-P of rank r=n-1 if k=1. Then 2<=r<n. Complementation is equivariant and preserves the normalized trace by uniqueness. In the standard symmetric model write

    A=Q11, B=Q22, C=Q12, T=A+B, S=AB-C^2,
    p_i=Q1i, q_i=Q2i, u_i=Bp_i-Cq_i, v_i=Aq_i-Cp_i  (i>2).

Let L_(1i),L_(2i) be the infinitesimal conjugation derivations from skew matrices with +1 in positions (1,i),(2,i) and -1 in the transposed positions. Matrix commutators give

    L_(1i)S=2u_i, L_(2i)S=2v_i,
    L_(1i)p_i=Qii-A, L_(2i)q_i=Qii-B,
    L_(1i)q_i=-C, L_(2i)p_i=-C,
    L_(1i)C=q_i, L_(2i)C=p_i.

The projector equations imply, with sums over i>2,

    sum p_i^2=A-A^2-C^2,
    sum q_i^2=B-B^2-C^2,
    sum p_i q_i=C(1-T),
    sum Qii=r-T.

Expanding and substituting these identities gives

    sum (u_i^2+v_i^2)=S(T-2S),
    sum (u_i p_i+v_i q_i)=S(2-T),
    sum [L_(1i)u_i+L_(2i)v_i]=(r-1)T-2(n-1)S,
    sum [L_(1i)p_i+L_(2i)q_i]=2r-nT.               (6.1)

Every invariant trace annihilates a Lie derivative. Put mu_j=Lambda(S^j), nu_j=Lambda(S^jT). For j>=1, apply Lambda to the sum of L_(1i)(S^(j-1)u_i)+L_(2i)(S^(j-1)v_i). Using (6.1) gives

    2(n+2j-3)mu_j=(r+2j-3)nu_(j-1).               (6.2)

Applying it to the sum of L_(1i)(S^j p_i)+L_(2i)(S^j q_i), for j>=0, gives

    (n+2j)nu_j=2(r+2j)mu_j.                       (6.3)

Combine (6.2) with (6.3) at j-1. With mu_0=1,

    mu_j/mu_(j-1)
      =[(r+2j-3)(r+2j-2)]/[(n+2j-3)(n+2j-2)],
    mu_j=(r-1)_(2j)/(n-1)_(2j).                  (6.4)

Here (x)_j=x(x+1)...(x+j-1) and (x)_0=1. The nonzero integer denominators are invertible in every characteristic-zero field. These are the even moments of the normalized polynomial density

    rho(z)=kappa z^(r-2)(1-z)^(n-r-1),
    kappa=(n-2)!/[(r-2)!(n-r-1)!].                (6.5)

Indeed direct beta integration gives its 2j-th moment as (r-1)_(2j)/(n-1)_(2j). Since monomials span the polynomial algebra, (6.4) proves the complete functional (4.1) for S. No distributional moment-determinacy theorem is needed.

The degree is b=n-3, its least exponent is s-1=r-2, and its least coefficient is kappa>0. Section 4 applies with

    a=n-2, F=p(S,t)^(n-2)(1+p(S,t)), H=t^(r-1).

For every k>=n-2 the marked seed moment is the positive rational

    lambda(Hp^k)=kappa integral_0^1 w^(r-2)(1-w^2)^(2k)dw
      =kappa (2k)!/[(r-1)((r+1)/2)_(2k)].          (6.6)

Every marked moment of F is the positive binomial sum (4.4). A coordinate projector of rank r avoiding the first coordinate axis exists because r<=n-1; at it S=0. Section 4.4 proves absolute unit-nonpower. The least detecting integer monomial exponent for F is r-1. In the codimension-one case r=n-1, (6.5) is the monomial density (n-2)z^(n-3), so the sharp pure-power and marker equivalences (4.5) also apply with s=n-2. No least pure-power assertion for the remaining two-minor densities is imported from that special case.

### 6.1. Arbitrary nondegenerate quadratic forms over K

This extension is independent of any compactness or positivity over K. Let q be ANY nondegenerate quadratic form on an n-dimensional K-space, n>=3. Diagonalize it to diag(a_1,...,a_n) with a_i!=0, and use its self-adjoint rank-k projector scheme. Its geometric transitivity and smoothness follow from Section 2.1 after scalar extension. In coordinates a_1 Q12=a_2 Q21, and the required regular scalar is

    S=Q11 Q22-Q12 Q21.                            (6.7)

Over the finite extension adjoining square roots d_i with d_i^2=a_i, conjugation by diag(d_i) changes q to the standard identity form. It preserves Q11,Q22 and the product Q12 Q21, so it preserves (6.7). Uniqueness of the extended normalized trace identifies its moments with (6.4). Injectivity of the field extension then proves the same rational moment identities over K. The same F,H and nonzero rational marked values follow. A diagonal coordinate projector avoiding the first axis is already K-defined and has S=0, so the geometric nonpower argument also applies. No transitivity on K-points is assumed.

## 7. AIII and CII: a uniform Hermitian-projector recurrence

We first calculate over R using the division algebra F=C or H of real dimension d=2 or 4. Let P be its rank-k Hermitian projector, and put D=P11, a real scalar. For each i>1 and each standard orthonormal basis unit e of F, take the skew-Hermitian infinitesimal matrix with entries e at (1,i) and -e* at (i,1). Write its conjugation derivation as L_(i,e), and put

    x_(i,e)=Re(P1i e*).

The off-diagonal generator has zero trace and belongs to the appropriate Lie algebra, including su(n) in the complex case. A direct commutator gives

    L_(i,e)D=2x_(i,e),
    L_(i,e)x_(i,e)=Pii-D.

The projector equation and the scalar trace give

    sum_(i,e) x_(i,e)^2=D-D^2,
    sum_(i,e)(Pii-D)=d(k-nD).                     (7.1)

Apply Lambda to the sum of L_(i,e)(D^(j-1)x_(i,e)), j>=1. With mu_j=Lambda(D^j) this yields

    (dn+2j-2)mu_j=(dk+2j-2)mu_(j-1),
    mu_j=(dk/2)_j/(dn/2)_j.                      (7.2)

All quantities in this calculation are rational coefficient coordinates in the E or D_H model. Alternatively the real identities in (7.2) descend to Q by Section 3 and extend to K. They are precisely the beta moments with

    alpha=dk/2, beta=d(n-k)/2.

Both are positive integers in the stated ranges. Therefore the polynomial-functional hypothesis (5.1) holds in every characteristic-zero K, with parameters (k,n-k) for AIII and (2k,2(n-k)) for CII. The CDF phase E(v,t) and marker t^-1 give exactly (5.4) in both families. The divisor t-1 supplies absolute unit-nonpower. This includes the usual complex and quaternionic projective-space orbits as rank-one cases, in their affine real algebraic projector models.

## 8. The compact Haar and Gaussian facts needed for AI–DIII

A Haar SU(m) matrix, m>=2, has first row uniform on the complex unit sphere. Complete any given unit row to a unitary basis and correct the determinant on the remaining rows; this proves SU(m) transitivity. Pushforward of Haar is the unique invariant probability on this compact transitive sphere. For m>=3, conditional on the first row u the second row is uniform in the unit sphere of u's Hermitian orthogonal complement: the stabilizer SU(m-1) acts transitively there. We use this conditional assertion only in AII, where m=2n>=4.

Likewise Sp(n) is transitive on the unit sphere in H^n by quaternionic orthonormal completion. Its Haar first row is uniform on the real sphere in R^(4n). This proof is independent of any particular random quaternionic matrix algorithm.

For a uniform complex unit vector u in C^m and nonnegative integers a_i,

    E product_i |u_i|^(2a_i)=product_i a_i!/(m)_(sum a_i).  (8.1)

To prove this, take independent standard complex Gaussians Z_i. Their direction u is uniform, independent of their squared radius R^2, by the polar-coordinate density. The latter has Gamma(m,1) density, so its A-th moment is (m)_A, while the product Gaussian moment is product a_i!. Division gives (8.1). The sphere law is invariant under independent coordinate phases, so any unbalanced powers of a coordinate and its conjugate have zero expectation.

More generally, split a uniform real sphere in dimension p+q into coordinate blocks of dimensions p,q. Independent real Gaussian blocks have independent squared radii with Gamma(p/2) and Gamma(q/2) densities (a common scale cancels). The change of variables from the two radii to their sum and ratio has Jacobian equal to their sum; the resulting density factors into gamma and Beta(p/2,q/2). Thus the first squared block norm has this beta law. A coordinate on the unit sphere in R^d has density proportional to (1-x^2)^((d-3)/2). These direct proofs supply all beta laws used below.

## 9. AI: SU(n)/SO(n)

For n>=2 use the Q-defined Cartan matrix C=gg^T and let u be the first row of g on real compact points. Put

    w=C11=sum_i u_i^2, D=w bar(w), S=1-D.

Expand |w|^(2j). Independent coordinate-phase invariance kills every pair of distinct multi-indices. Formula (8.1) then gives

    Lambda(D^j)=j!^2/(n)_(2j)
                sum_(a_1+...+a_n=j) product_i (2a_i)!/(a_i!^2).

The generating identity

    product_(i=1)^n sum_(h>=0) binom(2h,h)x^h=(1-4x)^(-n/2)

makes the finite sum 4^j(n/2)_j/j!. The duplication formula for (n)_(2j), which follows by separating its even and odd factors, yields

    Lambda(D^j)=j!/((n+1)/2)_j.                   (9.1)

These are the Beta(1,(n-1)/2) moments. Their complements S=1-D therefore satisfy

    Lambda(g(S))=(n-1) integral_0^1 z^(n-2)g(z^2)dz.  (9.2)

One may check this directly by setting S=z^2 in the Beta((n-1)/2,1) integral. Thus rho is the monomial (n-1)z^(n-2), and Section 4 applies with a=s=n-1. The phase is F_(n-1)(S,t), the all-order marker is t^(n-1), and all marked values are given by (4.4). The identity coset has C11=1 and S=0, giving absolute unit-nonpower.

In this family the sharp specialization also says p^c is pure IFF c>=n-1, and for each pure p^c a fixed integer monomial t^ell is unsafe IFF ell>=n-1. It detects every order IFF n-1<=ell<=c. The n=2 case is included; no conditional second-row law is used here.

## 10. AII: SU(2n)/Sp(n)

Let n>=2 and C=gJ0g^T, D=C12 bar(C12). Given the first row u of g, write its second row as v. Then C12=uJ0v^T. The vector conjugate(uJ0) has norm one because J0 is orthogonal. It lies in u's Hermitian orthogonal complement because uJ0u^T=0, by skew-symmetry. Conditional on u, v is uniform on the complex sphere in this complement, of complex dimension 2n-1. Thus C12 is a uniform unit-vector coordinate in that dimension. Formula (8.1) gives

    Lambda(D^j)=j!/(2n-1)_j,                     (10.1)

the Beta(1,2n-2) polynomial functional. Its regular polynomial CDF is

    v(D)=1-(1-D)^(2n-2).

Section 5 gives phase E(v,t), marker t^-1, all-order marked value (-1)^(m-1)/(m+1), and absolute unit-nonpower.

There is also a normalizer witness from S=1-D: its square-root density is (4n-4)z^(4n-5), so a=s=4n-4. Hence F_(4n-4)(S,t) with marker t^(4n-4) gives the other normalizer family and its sharp monomial-density pure-power/marker equivalences. At the identity coset D=1, so S=0 and its nonpower proof is valid too. This alternative is not needed for the simpler CDF witness. The excluded n=1 quotient is a point; its degenerate beta parameter is not used.

## 11. CI and DIII: even polynomial spherical densities

For h>=0 put

    c_h=integral_0^1(1-z^2)^h dz
       =4^h(h!)^2/(2h+1)!,                       (11.1)

a positive rational, as follows by elementary beta integration or its factorial recursion.

### 11.1. CI: Sp(n)/U(n), n>=1

In C=g(iI)g*, write the first quaternionic row q_l=z_l+w_l j, with z_l,w_l in the fixed complex subalgebra. Let x be the coefficient of i in C11. Quaternion multiplication gives

    x=sum_l (|z_l|^2-|w_l|^2).

The first row is uniform on the sphere in R^(4n), and these two real coordinate blocks both have dimension 2n. Therefore V=sum|z_l|^2 has Beta(n,n) law and x=2V-1. Its density on [-1,1] is proportional to (1-x^2)^(n-1). For S=x^2 this says

    Lambda(g(S))=c_(n-1)^-1
                 integral_0^1(1-z^2)^(n-1)g(z^2)dz.  (11.2)

Here the polynomial density has degree 2n-2 and positive constant coefficient c_(n-1)^-1. Section 4 gives a=2n-1,s=1: the phase is F_(2n-1)(x^2,t), the marker is t, and every marked moment is

    c_(n-1)^-1 sum_(h=0)^m binom(m,h)
                    integral_0^1(1-w^2)^(2((2n-1)m+h))dw>0.

A point with x=0 is obtained by taking g diagonal with first entry (1+j)/sqrt(2) and other entries 1: conjugation of i then has zero i coefficient. This point is defined over Qbar and proves the required absolute nonpower assertion. The cutoff 2n-1 is sufficient; for n>1 no smaller-power exclusion is asserted. The least detecting integer exponent for this displayed F is 1. For n=1, rho=1 is monomial and the sharper statements of Section 4.3 apply; this is a genuine sphere, not a point.

### 11.2. DIII: chosen SO(2n)/U(n) orbit, n>=2

Put C=gJ0g^T and x=C21. The vector Ce1 is a unit vector perpendicular to e1. Haar conjugation makes its law invariant under the stabilizer SO(2n-1) of e1, which is transitive on the unit sphere of e1-perpendicular. Thus it is uniform on the sphere in R^(2n-1), within the SAME chosen SO orbit. The coordinate x has density proportional to (1-x^2)^(n-2). For S=x^2,

    Lambda(g(S))=c_(n-2)^-1
                 integral_0^1(1-z^2)^(n-2)g(z^2)dz.  (11.3)

The polynomial degree is 2n-4 and the positive least coefficient is c_(n-2)^-1, so a=2n-3,s=1. The phase F_(2n-3)(x^2,t) and marker t satisfy (4.4). To obtain x=0 in the chosen orbit, use a determinant-one rational orthogonal rotation fixing e1 and sending e2 to e3, e3 to -e2. Since J0e1=e2, conjugating by this rotation gives Ce1=e3 and C21=0. This exists for 2n>=4. Hence absolute unit-nonpower follows. Again no least pure-power cutoff is claimed for the multi-term density; the least detecting integer monomial exponent for the displayed F is 1. At n=2 the density is constant and Section 4.3 supplies the extra sharp statements.

## 12. Why the same displayed pairs work over every char0 field

All groups, affine quotient models, coordinate projections, heights, polynomials v, phases, and markers in Sections 2 and 6–11 are defined over Q. Real averaging in the moment calculations is Lambda_R, the scalar extension of Lambda_Q by (3.1). Each equality is either an exact rational Lie identity or, when first calculated over R, an equality whose two sides belong to Q: for a polynomial b in the prescribed rational coordinates, Lambda_Q(b) is rational. The injection Q->R proves that equality over Q. Tensoring Q->K then proves it for the same polynomial over EVERY characteristic-zero K. There is no assumption that K embeds in R or C.

This applies separately to every moment order, and hence proves all-order identities. In the normalizer cases the marked values are positive rational polynomial integrals and finite binomial sums of them. In the CDF cases they are the rational numbers (-1)^(m-1)/(m+1). All remain nonzero under Q->K. Positivity here is used to certify a rational number before scalar extension; no order is asserted on a general K.

The nonpower arguments are also geometric and stable under this transfer. For the normalizer cases the specified S=0 points are Qbar-points, and the same equations have such a point over an algebraic closure of any K: BDI uses a coordinate projector, AI the identity, CI a quaternion with coefficients in Q(sqrt(2)), and DIII a rational rotation. At those points the divisor of 1+p has multiplicity one. For the CDF cases the geometrically prime divisor t-1 already gives multiplicity one. Every displayed f therefore remains an absolute unit-nonpower, not merely a nonpower in the original rational coordinate ring.

This proves (1.3) and its nonpower assertion. It is transfer of explicit rational trace identities for the prescribed models. It does not infer a global Mathieu-property ascent or descent from the existence of an arbitrary phase over an extension field.

## 13. Formal reality: the two individual factors really are Mathieu

We include the complete rational Gram argument so it applies to all seven new B_Q without importing a theorem stated only for earlier model families.

### 13.1. A generic rational Gram lemma

Let B_Q be any Q-algebra with a Q-linear functional Lambda_Q whose real scalar extension is strictly positive on the square of EVERY nonzero element of B_R. Let K be formally real: by definition -1 is not a finite sum of squares in K. Then

    Lambda_K(f^2)=0 implies f=0 in B_K.             (13.1)

To prove this, express f=sum_(i=1)^d c_i a_i with c_i in K and the a_i Q-linearly independent elements of B_Q. Their scalar extensions remain linearly independent over R, since field extension is flat. The rational symmetric Gram matrix

    G_ij=Lambda_Q(a_i a_j)

is positive definite over R: every nonzero real vector x gives the square of the nonzero element sum x_i a_i, whose trace is strictly positive by hypothesis. Exact rational Gaussian elimination gives an LDL^T decomposition with all diagonal pivots d_i positive rationals. In particular there is an invertible rational linear change of coefficient vector c, with coordinates ell_i(c), such that

    Lambda_K(f^2)=sum_i d_i ell_i(c)^2.             (13.2)

No irrational square root is used. Write d_i=A_i/B_i for positive integers A_i,B_i. It is the sum of A_i B_i copies of the rational square (1/B_i)^2. Thus the right side of (13.2) is a finite sum of squares over K.

In a formally real field a zero sum of squares has every term zero: otherwise divide by one nonzero square and rearrange to express -1 as a sum of squares. Hence every ell_i(c)=0. Invertibility of the rational change gives all c_i=0, proving (13.1). This proof uses only the definition of formal reality, not a chosen ordering, Artin–Schreier theory, an embedding of K into R, or a sum-of-four-squares theorem.

For each of the seven B_Q, (3.2) verifies the generic lemma's hypothesis. Consequently the only pure phase on B_K is zero whenever K is formally real: its second pure moment already forces it to vanish. Zero is safe, so ker Lambda_K is Mathieu on EACH individual compact-model factor over EVERY formally real K.

### 13.2. Scalar Laurent tori over arbitrary characteristic-zero fields

The classical Duistermaat–van der Kallen theorem says that a complex scalar Laurent polynomial whose every positive power has zero constant term has Newton polytope excluding the origin. We use the theorem for any finite Laurent rank; the theorem concerns scalar coefficients.

Let f be a pure Laurent polynomial over an arbitrary characteristic-zero field K. Its finitely many coefficients generate a finitely generated field K0 over Q. This field embeds into C: choose a transcendence basis, map it to algebraically independent complex numbers, and extend across the resulting finite algebraic extension into the algebraically closed field C. The embedding preserves all nonzero coefficients and constant-term identities. Thus DvK shows that 0 lies outside the same finite exponent polytope. If f=0 the conclusion is immediate; otherwise strict separation gives a rational linear form positive on every support exponent, and after multiplying by a denominator an integer linear form L with positive minimum delta.

For any fixed Laurent marker h, let -C be the minimum L-value in its finite support. Every exponent of h f^m then has L-value at least m delta-C, positive for large m; hence it cannot be exponent zero. So CT(h f^m)=0 eventually over K. This proves the scalar torus Mathieu property over every characteristic-zero field and every finite rank. Only the finitely generated coefficient field of the given phase is embedded; no embedding of the whole K and no function-valued DvK theorem is used.

### 13.3. Given-phase scalar extension and torus forms

For completeness, if A is a K-algebra with a K-linear functional lambda and L/K is ANY field extension, a GIVEN f in A is pure over K IFF it is pure after extension. This follows from injectivity of K->L on its scalar moments. It is also safe over K IFF it is safe after extension: for ascent, write a fixed h in L tensor_K A as a finite sum sum l_i h_i and take the maximum of the finitely many eventual cutoffs for h_i; descent follows by testing h in A and injectivity. Thus unsafety for a phase already defined over K has the same equivalence. This is a statement about that given phase. It is not general global Mathieu ascent, since an extended algebra can have new phases.

If a torus is a form of a split Laurent torus, choose a finite splitting extension. Its normalized trace becomes Laurent constant term, so the split Mathieu property and given-phase descent prove the Mathieu property for the original torus over K. We use this in the small BDI boundary. This descent alone suffices; no unproved general ascent is needed.

### 13.4. Genuine product instability

For every formally real K, Section 13.1 proves ker Lambda_K is Mathieu on B_K and Section 13.2 proves ker CT_t is Mathieu on K[t,t^-1]. Yet (1.3) gives a pure phase in their product algebra with the SAME fixed marker surviving every positive order. Therefore the product trace kernel fails the Mathieu property over every such K. This includes Q and R. This assertion does not require, and does not claim, compactness of K-points when K is not R. Individual-factor Mathieu assertions for non-formally-real K are not deduced.

The abstract tensor-instability phenomenon is known from Long's rational interval–torus example together with the same positivity and Gram reasoning. Here it is realized by the particular seven homogeneous coordinate algebras and the displayed regular formulas.

## 14. Exact boundaries of the result

**Points and the real toral boundary.** The projector cases k=0,n are points, so their product is G_m and is Mathieu. AI n=1, AII n=1, and the chosen DIII n=1 quotient are also points and excluded. CI n=1 is not a point and is included. For BDI n=2,k=1, write P11=a,P12=b,P22=1-a; the equation is

    (2a-1)^2+(2b)^2=1.

This is the standard norm-one circle torus, so its product with G_m is a torus and is Mathieu by Section 13.3. The same conclusion holds for any nondegenerate binary quadratic form's rank-one projector model: after a splitting extension it is a one-dimensional torus (one may use xy as split form and parameterize its self-adjoint rank-one projectors by a nonzero scalar). Its finite-splitting descent proves Mathieu. The BDI construction requires 2<=r<n and makes no contrary claim. There are no further nontrivial point/toral exceptions asserted within the stated seven rank ranges.

**Affine, non-toral models.** The seven failure families are the specific smooth affine homogeneous varieties in Section 2, with the independent split radial torus. Their complete complex projective Grassmannians or flag varieties have different coordinate rings and are not substituted. The chosen DIII component is retained throughout. No claim is made here for general nonsymmetric compact homogeneous spaces, exceptional Cartan families, disconnected ambient models, or arbitrary central forms.

**Fields and forms.** The all-field theorem is for the stated standard Q-model base changes, plus the separately proved arbitrary-q BDI extension in Section 6.1. It does not classify arbitrary real or K-forms of the other six types. In particular the positive-square Gram argument uses the compact real base change of the prescribed Q-model. The complex spherical classification supplies neither a real phase nor a descent to another real form. The fixed-coordinate rational transfer is proved directly instead.

**Finite quotients and descent.** An invariant trace or a finite-cover trace identity cannot make an arbitrary phase invariant under a finite group. A phase intended for a quotient must be a genuine regular function of that quotient (or satisfy a separately proved descent condition). Once a phase descends through a finite étale quotient in the reductive homogeneous category, normalized algebra trace can push a fixed upstairs marker down with identical marked moments; that is a separate retained theorem, not a premise needed here. No finite-cover or marker averaging is used in the displayed seven constructions, and no closure assertion for arbitrary finite quotients is made. In particular, no averaging of a phase is asserted to preserve purity.

**Sharp statements.** The polynomial-density lemma proves a sufficient pure-power cutoff a=b+1 and the least detecting monomial exponent s. The two-sided pure-power and full monomial-safety equivalences are asserted only for the monomial densities of Section 4.3 and their listed specializations. There is no root-purity converse outside those proved identities, no unproved optimal cutoff for the multi-term densities, and no finite scan used to fill such a gap.

## 15. Sources and attribution

The all-order coordinate calculations, density reductions, endpoint extraction, polynomial CDF substitution, valuations, rational-field transfer, and rational Gram proof are written in this document. The following retained primary sources have deliberately narrow roles.

### 15.1. Characteristic-zero reductive projection

J. S. Milne, *Algebraic Groups* (2017), retained physical file
`research/foundations/extended-geometry/sources/milne-algebraic-groups-2017.pdf`, SHA256
`503f1a7d511632af05541c52fb2bf6fc36f5a75b0323f2f9110f03dcce11fd67`.
The imported algebraic facts are complete reducibility, Theorems 22.42–22.43, printed 477–478, and the Reynolds construction in Section 12l, Proposition 12.53, printed 248–249. These provide Section 3's normalized invariant projection; the stronger full-automorphism theorem is unnecessary.

The generic rational Gram argument also has a complete proof in the local written source
`research/foundations/extended-geometry/FIELD_MATHIEU_PRODUCT_AND_QUADRATIC_COUNTEREXAMPLES.md`, SHA256
`8537845bffd21a60bcadcaaff190501b7d9d6869f395680c241b373ee443d97a`, Section 9 F8. Its reusable generic conditional facts are the finite rational Gram/LDL lemma and formal-real square anisotropy, after their positivity hypotheses are established here. Earlier model-specific factor conclusions are not imported outside their old scopes. Section 13.1 reproduces the generic proof rather than relying on those model-specific claims.

### 15.2. Scalar Laurent theorem

J. J. Duistermaat and W. van der Kallen, *Constant terms in powers of a Laurent polynomial*, Indagationes Mathematicae 9 (1998), 221–231, retained
`research/foundations/extended-geometry/sources/dvk-powers-1998.pdf`, SHA256
`bb002d12b539a1ff3fdd099f343d3accd921ccf60f8a149f29b07cd0f7e87e1b`.
The mathematical import in Section 13.2 is the scalar Laurent Newton-polytope constant-term theorem, including the arbitrary finite-rank treatment in Theorems 4–5 and their complete proofs. Rank-one Theorem 2 and its complete printed 222–223 proof also suffice for the independent split factor alone. The finite-coefficient-field transfer, fixed-marker support argument, and torus-form descent are proved here. No theorem with coordinate-algebra-valued Laurent coefficients is attributed to DvK.

### 15.3. Long's known constructions and complex comparison

Christopher D. Long, *Counterexamples to the xz-Conjecture and the Mathieu Conjecture for SU(2)*, [arXiv:2607.19012v1](https://arxiv.org/html/2607.19012v1), retained primary HTML
`research/foundations/extended-geometry/sources/long-2607.19012v1.html`, SHA256
`abeec76e95cba0e6244f4af31240634e28d657116c8afcf743b58db144dd7ebd`.
Theorem 2.1 gives the interval–torus witness, and Corollary 2.2 and Proposition 3.1 give related mixed constructions. Formula (5.3) is the same known interval mechanism, reproduced here. This paper already supplies the abstract real tensor-product counterexample after the elementary separate-factor arguments; the all-formally-real abstract example likewise follows from its rational identities and Gram reasoning. Neither phenomenon is claimed as a newly discovered abstract obstruction.

Long, *A Classification of Mathieu's Second Conjecture for Reductive Spherical Pairs*, author draft dated 22 September 2026, exact repository revision
`e9edac262b4b64eb68feb0e6e80496cefe4a1dee`, [pinned primary TeX](https://github.com/long-mathematics/mathieu-second-conjecture-spherical-pairs/blob/e9edac262b4b64eb68feb0e6e80496cefe4a1dee/mathieu_second_conjecture_spherical_pairs.tex).
Retained author PDF
`research/foundations/extended-geometry/sources/long-spherical-e9edac2.pdf`, SHA256
`b2bf25fbc2d6211333d4ff3e1f6de0b6562990e869c6078a2e219ef746f1026b`;
matching TeX SHA256
`934ad65de5ed4358fb33b29af15728ae86a1f7c5836e72050b00fc948c78df64`.
Section 5 Proposition 5.1 and its complete proof give exactly the normalizer circuit p=t^-1(1+t)[1-(1+t)^2S]^2 and endpoint method in its complex normalizer example. That known method is credited throughout; Section 4 reproves the precise polynomial-density form and its domain. Theorem 1.1 classifies the COMPLEX reductive spherical setting, and already covers existence of complex counterexamples for the seven nontrivial complex symmetric homogeneous products above. This comparison is context, not a logical premise of the real, rational, all-field, or nonpower construction.

### 15.4. Classical determinant and compact-Haar comparisons

Alain Rouault, *Asymptotic behavior of random determinants in the Laguerre, Gram and Jacobi ensembles*, ALEA 3 (2007), 181–230, [published primary PDF](https://alea.impa.br/articles/v3/03-09.pdf), retained PDF SHA256
`a148af1e0d1abe0e4c81e04277350700963e734bbd8a02c635e29386ed99b3fd`.
The actual Section 2.1, printed 184–185, and Section 2.2.3, formulas (2.15)–(2.17), including the complete proof of Proposition 2.4, printed 190–191, identify the BDI determinant law as classical Jacobi/matrix-beta theory. Its size-2 substitution, n1=r,n2=n-r,beta'=1/2, gives

    (r/2)_j ((r-1)/2)_j / [(n/2)_j ((n-1)/2)_j],

which is (6.4) by duplication. The read scope includes the singular codimension-one boundary. The retained erratum SHA256
`9256d1179b09def017f33a0dc0d6d829d0381e9a23cd29a201b14a387da480ca`
does not change these used Jacobi formulas. This is attribution and an independent comparison; the actual determinant law used in the construction is proved by the two Lie recurrences in Section 6.

Elizabeth Meckes, *The Random Matrix Theory of the Classical Compact Groups*, [author PDF](https://case.edu/artsci/math/esmeckes/Haar_book.pdf), retained
`research/foundations/extended-geometry/sources/meckes-haar-book.pdf`, SHA256
`12c9e49c9768f9f87951a59e6f4e960bf4a5bf6aa6f7fc5c8b442076e57d3dbe`.
Only Section 1.2 geometric orthonormal-column construction and complete invariance proof, printed/PDF 19–20, including its real/unitary and SO/SU adjustments, are credited. The retained draft's later symplectic-algorithm paragraph is not imported. The quaternionic sphere fact is proved independently by transitivity in Section 8, and the Gaussian moment calculations are supplied there.
