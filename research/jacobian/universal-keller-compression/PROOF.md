# Universal Keller maps with one fewer coordinate

14 November 2026. For every N>=5, the universal quadratic-gauge map can be compressed from A^N to A^(N-1), over the original field. The compressed map retains every rank-N finite etale algebra as a complete fibre, geometric degree N, symmetric monodromy, and absolute stable atomicity in characteristic zero. The scaling quotient has a global rational-power-free section on its clean torus. This removes a redundant parameter; it does not remove the intrinsic N-4-dimensional obstruction to descent through arbitrary primitive presentations.

The construction below gives a self-contained Jacobian and fibre calculation. It also proves Zariski-local universality for finite etale families over Q-schemes and a positive-characteristic extension with explicit sufficient bounds.

## The polynomial map

Write a=(a4,...,aN), and in three source coordinates put

\[
 t=1+xy,\quad q=t^2z+y^2(1+3t),
\]
\[
 \Pi=tq,\quad
 B=y+3xq+\sum_{j=4}^N j a_jt^2x^{j-2}q^j,\quad
 C=x(5-3t)-x^3z-\sum_{j=4}^N(j-2)a_j(xq)^j.
 \tag{1}
\]

Let F_a=(Pi,-B/2,C). These are polynomials over Z[1/2]. They are the normalized quadratic-gauge formulas in [URKC1/UAM1](../../foundations/verified/UNIVERSAL_RELATIVE_KELLER_MAP.md).

For N>=5 set

\[
 a_4=a_5=w,\qquad a_j=v_j\quad(j>=6).
 \tag{2}
\]

Promote w,v6,...,vN as unchanged source and target coordinates:

\[
 \boxed{\mathcal W_N=(w,v_6,\ldots,v_N,F_a):
             \mathbb A^{N-1}\longrightarrow\mathbb A^{N-1}.}
 \tag{3}
\]

There are N-4 promoted parameters and three moving coordinates. No parameter is inverted in the definition of this polynomial map.

## A short Jacobian and complete-fibre proof

Fix the parameters. On the rational open where x*t is nonzero, use coordinates (pi,S,t), where

\[
 \pi=tq,\qquad S=x/t.
\]

Writing b=-B/2 and c=C, direct substitution gives

\[
 b=\frac{t^{-1}-1-3\pi S^2-\sum_{j=4}^N j a_j\pi^jS^{j-1}}{2S},
 \quad
 c=2(S+bS^2+\pi S^3+\sum_{j=4}^N a_j\pi^jS^j).
 \tag{4}
\]

For fixed pi let D=1+2bS+3pi*S^2+sum j*a_j*pi^j*S^(j-1). Equation (4) says D=1/t. Differentiating the second identity gives, modulo d pi,

\[
 dc=2t^{-1}dS+2S^2db,\qquad b_t=-\frac1{2St^2}.
\]

Consequently

\[
 \det\frac{\partial(\pi,b,c)}{\partial(\pi,S,t)}
    =\frac1{St^3}.
\]

Also pi_z=t^3 and det partial(S,t)/partial(x,y)=x/t=S. Thus the other Jacobian is S*t^3 and their product is one. Since (1) is polynomial over Z[1/2], this rational identity proves the global polynomial identity

\[
 \det DF_a=1,\qquad \det D\mathcal W_N=1.
 \tag{5}
\]

The promoted variables contribute an identity diagonal block. This proof works over every field of characteristic other than two.

For a target (a,pi,b,c) with pi invertible define

\[
 E(S)=S+bS^2+\pi S^3+\sum_{j=4}^N a_j\pi^jS^j-c/2,\quad D=E'(S).
 \tag{6}
\]

Every source point over this target has t and q invertible, because tq=pi. Then S=x/t is defined, (4) holds without requiring x invertible, and E(S)=0, D=1/t. The apparently divided expression for b in (4) is replaced here by the polynomial identity below.

Conversely in the algebra (k[S]/E)[1/D] put

\[
 y=-2b-3\pi S-\sum_{j=4}^N j a_j\pi^jS^{j-2},
 \quad x=S/D,\quad
 z=\pi D^3-y^2(D^2+3D).
 \tag{7}
\]

The identity S*y=1-D gives t=1/D. Substitution gives q=pi*D, -B/2=b, C=c+2E=c. These formulas are inverse to S=x/t. The calculation is valid over any coefficient ring with 2 and pi invertible. Hence it is a scheme-theoretic isomorphism

\[
 \boxed{F_a^{-1}(\pi,b,c)
       \simeq\operatorname{Spec}(k[S]/E)[1/E'].}
 \tag{8}
\]

If the leading coefficient and discriminant of E are units, E' is already invertible modulo E. The complete fibre is the finite etale algebra k[S]/E. Repeated roots must otherwise be removed by the localization: dropping E' without squarefreeness would be false.

## Root-free compression and its inverse

For a nonzero scalar alpha define

\[
 a'_j=\alpha^{j+1}a_j,\quad
 (x',y',z')=(\alpha x,\alpha^{-1}y,\alpha^{-2}z).
\]

The exact covariance identity is

\[
 F_{a'}(x',y',z')
       =(\alpha^{-2}\Pi,\alpha^{-1}b,\alpha c).
 \tag{9}
\]

Indeed t'=t, q'=alpha^(-2)q, and every term of B,C has weight -1,+1 respectively. The inverse equation transforms by

\[
 E_{a',\alpha^{-2}\pi,\alpha^{-1}b,\alpha c}(S)
              =\alpha E_{a,\pi,b,c}(S/\alpha).
 \tag{10}
\]

On a4*a5!=0 choose

\[
 \alpha=a_4/a_5,\quad
 w=\alpha^5a_4=\alpha^6a_5=a_4^6/a_5^5,\quad
 v_j=\alpha^{j+1}a_j\ (j>=6).
 \tag{11}
\]

This uses no field extension or root extraction. Its inverse is

\[
 a_4=\alpha^{-5}w,\quad a_5=\alpha^{-6}w,\quad
 a_j=\alpha^{-(j+1)}v_j.
 \tag{12}
\]

Together with (9), it identifies the original promoted map on this parameter open with id_Gm times the compressed map on w!=0, by separate Laurent-polynomial isomorphisms of source and target. This is an open-set equivalence; it is not an equivalence of their entire affine spaces.

### The scaling quotient

On the clean torus where all a_j are units, (11)-(12) give an explicit isomorphism

\[
 (\mathbb G_m)^{N-3}
       \simeq\mathbb G_m\times(\mathbb G_m)^{N-4}.
 \tag{13}
\]

The scaling action sends alpha to alpha/s and fixes w,v6,...,vN. Thus a4=a5 is a global section of this torus quotient, over any base ring. The action is free because weights 5 and 6 are coprime.

The previous stable-moduli coordinates become

\[
 I_5=1/w,\quad J_6=v_6/w,\quad
 J_j=v_{j-2}v_j/v_{j-1}^2\quad(j>=7),\quad v_4=v_5=w.
 \tag{14}
\]

Conversely w=I5^(-1) and v_j=J_j*v_(j-1)^2/v_(j-2). These are Laurent formulas in both directions. This directly supplies the invariant-ring generators without an integer normal-form calculation.

At N=4 only weight 5 is available. Normalizing a4=1 over the original field can require a fifth root: a4=2 over Q is a countercontrol, by its 2-adic valuation. No dimension-three universal quartic map follows from this argument.

## Complete fibre universality and the new dimension bound

**UACOMP1.** For every characteristic-zero field k and N>=5, the map W_N has determinant one, geometric degree N, geometric and arithmetic monodromy S_N, and realizes every rank-N finite etale k-algebra as a complete fibre. It is absolutely and stably atomic.

Take a monic squarefree degree-N presentation P(T) of the algebra. Such a presentation exists over every infinite field, including disconnected finite etale algebras. Choose a in k so that

\[
 d_1d_3d_4d_5\ne0,\qquad d_j=P^{[j]}(a).
 \tag{15}
\]

All four Hasse derivatives are nonzero polynomials in a. Their product has degree 4N-13, so one of the 4N-12 integers 0,...,4N-13 is admissible in characteristic zero. This is a uniform finite translation rule once P is given; it does not choose P canonically.

An explicit compressed target is

\[
 \begin{aligned}
 w&=\frac{d_3d_4^6}{d_1^2d_5^5},&
 v_j&=\frac{d_jd_3d_4^{j+1}}{d_1^2d_5^{j+1}}\quad(j>=6),\\
 \pi&=\frac{d_1d_5^2}{d_3d_4^2},&
 b&=\frac{d_2d_5}{d_3d_4},&
 c&=-\frac{2d_0d_3d_4}{d_1^2d_5}.
 \end{aligned} \tag{16}
\]

With alpha=d3*d4/(d1*d5), its inverse polynomial is

\[
 E(S)=\frac{\alpha}{d_1}P(a+S/\alpha).
 \tag{17}
\]

It has degree N and is squarefree. Formula (8) therefore gives the complete desired fibre, with T=a+S/alpha.

### Monodromy without a Morse-slice argument

The original generic inverse (6), with all a_j independent, is the generic monic polynomial in a rational coefficient chart. To see this, divide a generic monic polynomial by its nonzero linear coefficient. This birationally identifies its coefficient field with that of a polynomial whose linear coefficient is one. Then h3=pi and h_j=a_j*pi^j recover (6).

The generic monic polynomial has splitting field k(r1,...,rN) over the symmetric rational functions, and permutation group S_N, in every characteristic. Thus the old generic cover has degree N and group S_N. Equations (9)-(12) identify it with the compressed generic cover after adjoining the independent variable alpha. Purely transcendental extension does not change a finite separable Galois group or its root action. Therefore the compressed cover also has degree N and group S_N, geometrically and arithmetically. This proof does not require generic Morse slices, a discriminant branch-cycle classification, or the weighted-seed monodromy theorem UM1.

The existing [primitive-monodromy atomicity proof](../../foundations/verified/PRIMITIVE_MONODROMY_ATOMICITY.md) now applies: a nontrivial polynomial factorization gives an intermediate field; the natural S_N action is primitive; a degree-one Keller factor is an automorphism. Adjoining independent identity variables preserves this field action, as does extending algebraically closed constants. This proves absolute stable atomicity.

In the programme's notation the result improves the existing upper bound to

\[
 \boxed{\operatorname{ed}_k(S_N)\le
        \operatorname{ktdim}_k(S_N)\le N-1\quad(N>=5).}
 \tag{18}
\]

The new argument proves the upper bound; the essential-dimension lower bound is the existing independent input. No optimality of N-1 is asserted. For N=3,4 the previous bounds are unchanged. N-4 is the minimum parameter dimension for domination of the displayed clean stable-class torus, by dimension alone; it is not a minimum for arbitrary field-versal Keller constructions.

## Zariski local universality for families

**UAFAM1.** Let S be any Q-scheme and A a finite locally free etale O_S-algebra of rank N>=5. There is a Zariski open cover S_i and morphisms y_i:S_i->A_Q^(N-1) such that Spec(A)|S_i is the complete pullback fibre of W_N along y_i.

At a point s, the residue field kappa(s) is infinite. Choose a primitive element in the rank-N residue algebra and an admissible translation satisfying (15). Trivialize the finite locally free module near s and lift these choices to local sections. The determinant of the power basis is nonzero at s, hence is a unit after shrinking. The characteristic polynomial is then a monic presentation of A. The four required Hasse coefficients are also units after shrinking. Equations (16)-(17), valid over rings with those units, give the target and complete-fibre isomorphism. Doing this at every point proves the assertion. Nilpotents on S cause no problem; the power-basis determinant and discriminant are units, not merely generically nonzero.

This is local existence with choices, not canonical global gluing. Schemes of generators and finite-etale algebras with generators are established constructions; see the [primary source manifest](sources.json). The new conclusion here is their realization by the compressed Keller map.

## Positive characteristic and a finite field bound

**UAPOS1.** The same W_N is a determinant-one degree-N map with monodromy S_N and absolute stable atomicity over every field of characteristic p>N. It realizes every rank-N finite etale algebra over an infinite such field. Over F_q it does so whenever q>4N-13.

The Jacobian and fibre calculations above apply, and the leading coefficients of all four Hasse derivatives in (15) remain nonzero when p>N. Infinite-field primitive elements and the nonzero-polynomial argument prove the infinite-field assertion. The generic-polynomial argument proves monodromy in this characteristic. The birational-etale/open-immersion/Hartogs proof of the degree-one Keller lemma works in this characteristic too, so the same primitivity proof gives absolute stable atomicity. This last statement concerns these maps; it does not assert the characteristic-p Jacobian conjecture.

For completeness, every rank-N etale F_q-algebra is monogenic if q>=N. Write it as a product of m_d copies of F_(q^d), with sum d*m_d=N. Distinct isomorphic factors need distinct monic irreducible degree-d polynomials. There are q such polynomials for d=1. For d>=2, elements of F_(q^d) not generating the full field lie in proper subfields and number at most

\[
 \sum_{i=1}^{\lfloor d/2\rfloor}q^i\le q^{d-1}.
\]

Thus the number of monic irreducibles is at least (q^d-q^(d-1))/d >= q/d >= m_d. Choose distinct ones and use the Chinese remainder theorem. This proves monogenicity. Since q>4N-13 implies q>=N for N>=5, choose P accordingly; at most 4N-13 translations fail (15), so an admissible a exists in F_q.

These bounds are sufficient, not sharp. For small characteristic, P=T^p-T can have identically zero third Hasse derivative, so the compiler can fail. For small finite fields, F_q^N is not monogenic when N>q. Neither example disproves the existence of some different universal map.

## What the compression does and does not make canonical

**UAQUOT1.** Equations (11)-(14) give a root-free global section of the clean scaling quotient. The generic Tschirnhaus non-descent obstruction remains after compression.

The [existing theorem GTND1](../../foundations/verified/GENERIC_TSCHIRNHAUS_NON_DESCENT.md) says that primitive presentations of the same split algebra move dominantly through the N-4-dimensional stable-boundary torus. Equations (13)-(14) merely give explicit coordinates on that very torus. Distinct values remain distinct after any extension of fields, or after dominant base change of the presentation relation. Thus removing the redundant scaling coordinate cannot turn the compiler into a presentation-free stable-class assignment.

Here the varying stable classes are those of the three-dimensional slices F_a (or their marked incidence data). The promoted maps U_N and W_N are each a single fixed map: changing a target does not change that map. The obstruction must not be misread as a non-descent claim for their unmarked constant families.

A separate, elementary categorical obstruction is worth keeping distinct. No finite-etale family E->T over an algebraic space T can have a functorial choice y(A) in T, together with a natural identification E_y(A) with A, for all rank-N algebras and their isomorphisms when N>=2. At the split algebra k^N, every S_N automorphism must map to the identity of its image in T. Naturality would then make that nontrivial automorphism act trivially on k^N, a contradiction. Equivalently, T->BS_N has no section as a morphism of stacks. This is the standard inertia obstruction, not a new general impossibility theorem for stack-valued receivers.

It does not prohibit a noncanonical choice on isomorphism classes, local presentations, a suitable stack with retained inertia, or another compiler. GTND1 is stronger in a different direction: even after forgetting inertia and considering only the displayed coarse stable classes, the particular compiler varies with primitive presentation. These are independent obstructions and should not be conflated.

## Proof assurance and remaining questions

The universal arguments are written proofs. The checker verifies exact polynomial Jacobians, covariance, inverse-equation transformations, rational and finite-field complete-fibre controls, quotient identities, and failure boundaries. Those controls do not independently verify the all-rank geometry, monodromy, or family descent. No full formalization, external review or novelty-priority claim is made.

The new universal dimension bound is N-1. Whether it is optimal, whether rank four can be realized universally in dimension three, whether a different compiler can reduce dimension further, and how to construct a presentation-free stack or correspondence receiver all remain open.

