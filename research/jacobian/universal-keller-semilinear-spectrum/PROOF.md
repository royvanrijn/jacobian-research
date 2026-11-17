# Semilinear isogenies and the exact prime-power degree spectrum

15 November2026. Written proofs with narrow exact controls. No formal verification, independent whole-proof replay, external review or publication-priority certificate is claimed. Classical arithmetic inputs are pinned in `sources.json`. OP-KDESC remains open.

The main new conclusion is a degree-spectrum theorem: in characteristic zero, a nontrivial strongly absolutely stably atomic finite etale self-map of a smooth affine geometric torus quotient has prime-power degree, and **every prime-power degree occurs**, including volume-preserving field-versal covers for specified finite etale groups. This does not assert full symmetric universality or an affine-space source.

## 1. UAISOTWIST1: allow a holonomy automorphism in the isogeny criterion

The commuting-isogeny construction UAISOCRIT1 extends to intertwiners. Let F be finite, alpha in Aut(F), S a split torus with action rho, and let phi:S->S be an isogeny of degree D satisfying

    phi rho(g)=rho(alpha(g)) phi.                             (1.1)

Assume char(k) does not divide D|F|. A free field-versal torus base B with an automorphism b satisfying b g=alpha(g) b can always be constructed. Start with B0 from UAFREEBASE1; if alpha has order t, use

    B=product_(j=0,...,t-1) B0^(alpha^j),
    g*(x_j)_j=(alpha^j(g)*x_j)_j,
    b((x_j)_j)=(x_(j+1))_j.                                 (1.2)

Each factor is free, and every F torsor twist has a point factorwise. The shift b is a k-automorphism and satisfies the required intertwining. No condition on the characteristic relative to t is needed: b is an automorphism, not a quotient by that cyclic shift.

The pair (b,phi) descends to a finite etale degree-D self-map of

    X=(B x S)/F.

With A=ker(phi), the action of A semidirect F on B x S is free, and quotienting by A followed by the pair (b,phi) identifies the target with the same X, relabelling F by alpha. For an F torsor Q/K, the isogeny twists to

    phi_Q:S_Q -> S_(alpha_*Q).                               (1.3)

Exactly as in UAISOCRIT1, the fixed A semidirect F torsor is field-versal if and only if the map in(1.3) has zero kernel on H1, for every Q and K. Here field-versal means every specified group torsor over every extension field, including finite ones, is realized; it supplies no natural choice of a realizing point. The possible affine S_Q-twist classes are the image of H1(K,A_Q), equal to that kernel. The adjugate isogeny intertwines alpha^-1 and composes with phi to [D]. Thus that kernel is killed by D and by |F|, so

    gcd(D,|F|)=1                                             (1.4)

is sufficient. More precise uniform cohomology annihilators may still replace |F|. The kernel's source F module determines primitive affine monodromy just as before. If the kernel is not constant, use its actual finite etale group scheme; no arbitrary constant identification is inferred.

For cyclic F=C_m the full cyclotomic etale algebra

    E=k[z]/Phi_m(z), B=Res_(E/k)Gm

is often smaller than(1.2). Translation by z gives the free C_m action. Every u in (Z/m)^* induces the automorphism z->z^u of E and hence b of B, intertwining g->g^u. This requires the FULL etale algebra: the smallest field k(zeta_m) need not have that automorphism over k. In particular a one-dimensional split translation base would allow only exponent +1 or -1, and cannot silently replace this B.

## 2. UACYCSEMI1: the ideal-class equation and its solution at prime conductor

Put L=Q(zeta_m), O=Z[zeta_m], and let sigma_u(zeta_m)=zeta_m^u. Choose an O ideal lattice I, viewed as the cocharacter lattice of a split torus S_I with cyclic holonomy multiplication by zeta_m. Every integral intertwiner with holonomy automorphism u has form

    T(v)=a sigma_u(v),  a in L*,  a sigma_u(I) subset I.       (2.1)

If its image is P I for a nonzero prime ideal P, its degree is Norm(P), and existence is equivalent to the class-group equation

    [P]=(sigma_u-1)[I] in Cl(O).                             (2.2)

Indeed a sigma_u(I)=P I exactly says that P I sigma_u(I)^-1 is principal. Conversely a solution class and a generator of that principal fractional ideal give the integral map(2.1). Its cokernel I/P I is a one-dimensional vector space over O/P. Thus the kernel is an irreducible cyclic module, up to the holonomy relabelling and torus duality, whenever O/P has characteristic p with p not dividing m. The resulting self-cover is strongly absolutely stably atomic, of degree p^f, f=ord_m(p), and is field-versal by(1.4).

This does not contradict UACYCIDEAL1. For its fixed lattice I=O, equation(2.2) requires P principal. Allowing a different ideal lattice and a nontrivial holonomy automorphism changes the exact obstruction. For a general conductor the image of sigma_u-1 can still be a proper subgroup; we make no claim that all such class equations are solvable.

### Prime conductor removes this class obstruction

Let m=ell be an odd prime, and let sigma generate Gal(Q(zeta_ell)/Q). Then

    Cl(O)^<sigma>=0,
    sigma-1:Cl(O)->Cl(O) is bijective.                       (2.3)

These are classical consequences of the ambiguous class number formula. In the ordinary class-group version, including infinite ramification, apply the formula to L/Q. Only ell and the real place ramify, with indices ell-1 and2. The units of Q are {+1,-1}; norms from the totally complex L are positive, so the unit-norm index is2. Since h(Q)=1, the formula gives

    |Cl(O)^Gal(L/Q)|=2(ell-1)/((ell-1)*2)=1.                 (2.4)

The class group is finite, so injectivity of sigma-1 is equivalent to bijectivity. We use the general cyclic formula in Li--Yu, Example2.5(1), not the prime-degree-only formulation in Lemmermeyer's shorter note. No class number one assumption is made.

It follows that for EVERY prime P above EVERY rational p different from ell, some ideal I and some a solve(2.2), using a generator sigma of the full cyclotomic Galois group. Together with the full cyclotomic base in Section1 this gives a fixed semilinear self-cover of degree p^ord_ell(p), regardless of whether P itself is principal. Matrices on a Z basis of I are integral; the source and target are the same quotient variety.

### A nonprincipal degree-2048 control

Let ell=23, let sigma(zeta)=zeta^5, and put

    tau=sum_(j a nonzero quadratic residue mod23) zeta^j,
    tau^2+tau+6=0,
    P=(2,tau),  I=P,  a=(2+tau)/2.                          (2.5)

There are two primes above2, each of residue degree11. Choose P with tau=0 modulo P. The other is sigma(P), since5 is a primitive root and is a quadratic nonresidue. In the quadratic subfield Q(sqrt(-23)), u=2+tau has norm8, is divisible only by the chosen prime above2, and generates its cube. Extending gives

    (u)=P^3,  (2)=P sigma(P),
    a sigma(P)=P^2 subset P.                               (2.6)

Hence a*sigma is an integral endomorphism of the lattice P of index2^11. In the standard degree22 cyclotomic basis, an explicit prime factor is

    g(z)=z^11+z^9+z^7+z^6+z^5+z+1 over F2.

It is irreducible, and P=(2,g(zeta)). The checker constructs its Hermite basis, the integral semilinear matrix, verifies(2.6), and obtains eleven Smith factors1 and eleven factors2. Its determinant is -2048 and it intertwines holonomy by the fifth power. The residue module has dimension11 and cyclic order23 action, hence is irreducible.

P is genuinely nonprincipal; this is not merely an example where a hidden principal generator was overlooked. In the quadratic subfield the chosen prime q above2 is nonprincipal: the norm equation

    A^2-A B+6 B^2=2

has no integer solution (complete the square), while q^3=(u). Its class therefore has order3. If P were principal, its relative norm q^11 would be principal, contradicting that order. Thus this is an explicit escape from the fixed-O principal-prime gate through a different ideal lattice and holonomy relabelling.

The full cyclic base has dimension22, as does S_P. The quotient self-map has dimension44 and monodromy F2^11 semidirect C23. Its logarithmic Jacobian is +2048: the cyclic-base automorphism has determinant -1, as does the semilinear sign above. Scaling a new nonzero volume coordinate by1/2048 gives a volume-preserving geometric torus-quotient self-map in dimension45. The model is defined over Z[1/46] and is field-versal for that specified constant group over every permitted field. This single control checks the semilinear mechanism; it is not a search over conductors or ideal classes.

## 3. UAPUNCTVOL1: keep volume lifting inside torus-quotient geometry

For any etale morphism f:X->Y between smooth d-folds, UAKANLIFT1 canonically lifts to the total canonical lines, preserving the tautological volume d(tau). Restrict this functor to the complements of the zero sections:

    K_X^times -> K_Y^times.                                 (3.1)

It is still the base change of f along K_Y^times->Y, with the same generic degree and geometric monodromy. It remains finite when f is finite, and preserves the restricted nowhere-zero volume form.

If X=T/F is a smooth geometric torus quotient, then K_X^times is another smooth affine geometric torus quotient. Over the torus cover choose its invariant logarithmic top form. Its F multiplier is det(rho(g)) in {+1,-1}. The punctured canonical line is

    (T x Gm)/F,
    g*(t,w)=(g*t,det(rho(g))^-1*w).                          (3.2)

The action remains free. Thus the operation adds one dimension inside the same category, rather than only producing an affine-line bundle. The generic extension is purely transcendental; primitive monodromy and UAFINATOM1 preserve strong absolute stable atomicity after every constant extension and identity stabilization.

For every field K, each point y of Y(K) has a nonzero K-covolume. Therefore pulling a fixed torsor back along K_Y^times->Y preserves its realization image on fields, and field versality persists. For general test schemes this map need not have a section, so its realization image can shrink; it never grows, as UAVERSLIMIT1 explains. No canonical choice of a nonzero covolume is asserted.

There is a precise marked-action boundary. If a finite group fixes y, a nonzero covolume over y is fixed by the canonical lifted action if and only if the determinant character of its tangent action at y is trivial. The character on the covolume line is its inverse. Thus puncturing preserves a fixed marked target exactly under that condition; the whole-line zero-section lift always preserves it. These alternatives must not be conflated when claiming canonical symmetry.

## 4. UAPPDEGREE1: the exact degree spectrum, in every degree

Let k be a field of characteristic zero and N>1 an integer. The following are equivalent:

1. N is a prime power.
2. Some smooth geometrically integral affine geometric torus quotient over k admits a finite etale degree-N self-map that is strongly absolutely stably atomic.
3. Such a self-map can be chosen to preserve a nowhere-zero global volume form and to arise from a fixed field-versal torsor for a specified finite etale group scheme, with its rank-N associated cover.

Necessity is UAINFRASELF1: the geometric monodromy has an abelian regular normal subgroup, and primitive affine monodromy forces an elementary abelian kernel. Strong atomicity in the stated finite normal-affine category is equivalent to primitivity by UAFINATOM1. Volume preservation does not weaken that necessity.

For sufficiency write N=p^e. When e=1, the power map on Gm is field-versal for mu_p, of prime degree, and its punctured-volume lift

    (x,w)->(x^p,w/p),   Omega=dlog(x) wedge dw

on Gm^2 supplies(3). For e>=2 use the following complete arithmetic construction. The Bang--Zsigmondy primitive-divisor theorem supplies an odd prime ell with

    ord_ell(p)=e                                             (4.1)

except precisely the relevant cases

    (p,e)=(2,6), or e=2 and p+1 is a power of2.              (4.2)

For the ordinary cases choose a prime P of Z[zeta_ell] above p. Its residue degree is e. The generator sigma of the full cyclotomic Galois group has sigma-1 bijective on the ideal class group by(2.3). Solve(2.2) and use the semilinear lattice isogeny a*sigma. Section1 gives a fixed self-cover with degree p^e, irreducible elementary abelian kernel and field versality because gcd(p^e,ell)=1. A uniform source uses the full etale cyclotomic base and has dimension2(ell-1). Section3 keeps volume preservation inside torus-quotient geometry with one extra dimension.

Both exceptional cases have small direct constructions:

- If e=2 and p+1 is a power of2, then p is odd and p=3 mod4. Use O=Z[i], F=C4, and multiplication by p. The quotient O/(p) is F_(p^2), since z^2+1 is irreducible modulo p. Use the full degree-two cyclotomic base. The source dimension is4 and its punctured-volume lift has dimension5. Coprimality p with4 gives field versality.
- For (p,e)=(2,6), use O=Z[zeta_9], F=C9 and multiplication by2. The polynomial Phi_9=z^6+z^3+1 is irreducible modulo2 because ord_9(2)=6. The source dimension is12 and the volume lift has dimension13. Its monodromy is F_64 semidirect C9 and is primitive, with field versality by gcd(64,9)=1.

This proves all prime powers, including every formerly remaining odd binary exponent. No assertion about infinitely many Mersenne primes, Artin's primitive-root conjecture, or cyclotomic class number one is used.

These particular spectrum models have trivial canonical line. The base action is by translations and the determinant of multiplication by zeta_m on S is1. Thus a global logarithmic top form Omega0 descends. In the generic prime-conductor construction sigma permutes the nonidentity roots in one cycle on both the full cyclic base and the cyclotomic fibre representation; their determinant signs cancel. The total logarithmic Jacobian is Norm(a)=N. In the scalar exceptional cases it is also N. Hence the volume map can be written explicitly on X x Gm as

    (x,w)->(f(x),w/N),  Omega=Omega0 wedge dw.                (4.3)

This uniform normalization makes no root choice and stays inside torus-quotient geometry. Field versality persists as in Section3; fixed marked actions still have that section's determinant boundary.

### Bounds, arithmetic formulation and limits

For e>=2 in the ordinary case, ell divides p^e-1, so the volume source dimension is at most2N-3. The exceptional dimensions5 and13 satisfy the same bound. The prime-degree source has dimension2. Thus one may choose

    dim(source)<=max(2,2N-3)                                (4.4)

for every prime-power N>1. This is an existence upper bound; no optimality or canonical marked symmetric action is claimed. The older quartic construction may have a smaller dimension and stronger full-symmetric universality.

For each chosen conductor and lattice, the construction has an integral model over Z[1/(p*m)]. Its geometric monodromy is F_p^e semidirect C_m (with m=ell,4 or9 in the cases above), with the indicated irreducible action. For p odd, constant-group torsor universality is stated after adjoining mu_p; over a general field use the actual finite etale kernel group scheme. For p=2 the kernel is already constant over every permitted field. The same fixed binary model is therefore field-versal for its specified constant group over every characteristic-zero field, and in the good positive-characteristic fibres of that model. The necessary prime-power classification is asserted in characteristic zero only.

These are universal realizations of the specified group torsors and their associated rank-N algebras. They are not universality for all rank-N etale algebras. In particular the symmetric self-universality ceiling4 from UAINFRASELF1 is untouched. Nonconstant units keep every constructed source outside affine space after affine identity stabilization. UAPUNCTVOL1's determinant condition remains necessary if one also requests a fixed marked target under a separately prescribed action.

## 5. Classical inputs, assurance and the surviving canonical problem

The arithmetic inputs are classical: the general cyclic ambiguous class formula, prime decomposition in cyclotomic fields and Bang--Zsigmondy's primitive-divisor theorem. Li--Yu Example2.5(1) explicitly includes the infinite ramification factor, which is essential to(2.4). The author-written Michels proof of Zsigmondy's theorem is pinned; its first page, containing the exact two exceptional cases, was visually inspected. The original1892 reference is identified without claiming a new primitive-divisor result.

The new repository deductions are the holonomy-twisted self-cover mechanism, its ideal-class equation, the prime-conductor class argument applied to self-isogenies, and the resulting complete degree spectrum. The exact matrix control at degree2048 confirms one nonprincipal semilinear realization; bounded controls do not prove the arbitrary-degree theorem. No publication-priority certificate has been established.

The original presentation-free whole-affine canonical problem remains open. Field versality still gives no natural choice of points and no fixed scheme receiver for all smooth affine families, by UAVERSLIMIT1. Effective simultaneous fixed-locus landing/gluing, actual compatible marked actions, atomicity when marked rank is smaller than generic degree, modular affine cases and A1-rich pseudo-plane mixing remain distinct open requirements. The prime-power degree-existence gap within finite torus-quotient geometry is closed by the written theorem; constructing arbitrary prescribed primitive affine groups or minimizing the dimensions is not claimed.
