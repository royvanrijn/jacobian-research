# Reduced branch pairings and symmetric trace matrix factorizations

17 November 2026. Characteristic zero over C, finite normal surface
covers of the affine plane. The actual polynomial Keller application
retains all hypotheses of the existing normalization and even-moment
theorems. JC2, HC4 and automatic full trace integrality remain open.

The trace defect has an exact reduced scheme structure: its annihilator
is the reduced branch equation, although its Fitting discriminant can
have arbitrarily high multiplicities. Multiplying the linking form by
that one reduced equation gives a perfect pairing over the branch curve.
At singular curve points this module need not be locally free, and
perfection does not give nondegenerate forms on closed-point fibres.

## 1. The exact annihilator and a simple-pole bound

Let B=C[P,Q], F=Frac(B), L/F finite separable and S the full normal
integral closure of B in L. Finite normal surface duality gives finite
B-projective S and its trace dual

    Lambda={z in L:Tr(z*S) subset B}.

As in LY0 and KN0, write Ram=sum_E a_E E, a_E=e_E-1, so that
Lambda=O_Z(Ram) for Z=Spec(S). For each target branch prime C choose
a reduced equation h_C, and put

    H=product_(C with some e_E>1 above C) h_C.

The product is squarefree; H=1 when the finite cover is etale. Set
D=Lambda/S. Then the exact global assertions are

    H*Lambda subset S,
    Ann_B(D)=H*B,
    H*Tr(Lambda*Lambda) subset B.                         (1)

At a source prime E above C, v_E(H) is e_E if C occurs in H, and
zero otherwise. Since e_E>a_E, H*Lambda has nonnegative valuations
at all such primes. At other primes Lambda=S. Normality gives the
first assertion at every finite point. It implies the third by
the definition of Lambda: H*lambda is in S and hence
Tr((H*lambda)*mu) is regular for all lambda,mu in Lambda.

For the reverse inclusion in the annihilator, localize at a branch
prime C. The module D_C is nonzero because it has length

    r_C=sum_(E above C) f_E*(e_E-1)>0.

It is killed by h_C, since all other factors of H become units.
It is therefore a nonzero vector space over B_(C)/(h_C).
An element of B annihilating it must vanish along C. This holds
for every factor of H, and factoriality of B gives H dividing that
element. This proves the second assertion, without an assertion that
D is locally free over the possibly singular ring B/(H).

The annihilator and the Fitting ideal must be distinguished. The
latter is generated locally by a trace Gram determinant, with order
r_C at C. The former has order one at every branch prime.

## 2. A perfect pairing on the reduced branch curve, at every point

LY0 proves the symmetric linking form

    beta:D times D -> F/B, beta([lambda],[mu])=Tr(lambda*mu) mod B

is perfect. Its values lie in H^-1*B/B by (1). Multiplication by H
identifies this target with A=B/(H), giving

    beta_H:D times D -> A,
    beta_H([lambda],[mu])=H*Tr(lambda*mu) mod H.             (2)

This is well-defined: changing lambda by s in S changes the numerator
by H*Tr(s*mu), which belongs to H*B. It is a perfect A-module pairing:

    D -> Hom_A(D,A) is an isomorphism.                    (3)

Indeed H kills D, so every B-linear map D->F/B takes values in the
submodule killed by H, exactly H^-1*B/B. Apply LY0's global perfection
and the displayed target identification. No generic-point gluing of
pairings is used, so (3) includes singular branch points.

The projective presentation S->Lambda of D also shows that D is a
maximal Cohen--Macaulay module on its reduced curve support: as a
B-module it has projective dimension one and depth one at closed
support points. It has no separate point-supported submodule.
None of these facts implies that D vanishes or is projective over A.

In a supplied local B-free basis of S, let G be its trace Gram matrix.
Then Lambda has coordinate lattice G^-1*B^d. Formula (1) says that

    N=H*G^-1 has regular entries, N^t=N,
    G*N=N*G=H*I_d.                                       (4)

Thus the trace matrix supplies a symmetric matrix factorization of the
REDUCED branch equation, despite its determinant having higher branch
orders. This is a local assertion in any free chart of the projective
B-module; no unproved global basis theorem is used.

The matrix-factorization framework is classical: Eisenbud,
[*Homological algebra on a complete intersection, with an application to
group representations* (1980), Section5, Proposition5.1, printed p49;
Theorem6.1, p52; Corollary6.3, p55](https://eisenbud.github.io/papers/pdfs/1980-002.pdf).
His definition applies directly to the pair (G,N) in (4). At a closed
point of H=0, the base B localized there is regular local and H is a
nonzero nonunit, giving the classical two-periodic resolution over
B/(H); after removing trivial summands, Corollary6.3 identifies the
associated maximal Cohen--Macaulay module. These statements and the
displayed Proposition5.1 proof were inspected in the original article.
The specific input H*Lambda subset S, the trace description of N and
its symmetry still come from the preceding argument. No new general
matrix-factorization or hypersurface-module theorem is claimed here.

## 3. Every higher-even-moment residual module has the same refinement

For every finite normal surface cover in Section 1 define the divisorial lattices

    T_m=O_Z(floor(Ram/m)), m>=2,
    T_m^#=O_Z(Ram-floor(Ram/m)),
    Q_m=T_m^#/T_m.

Each T_m is finite reflexive over the normal surface S. At a closed
base point the two regular base parameters form a system of parameters
at every point upstairs, and a regular sequence on these depth-two
reflexive modules. Thus they have B-depth two and are B-projective,
by the dimension-two regular-base projective-dimension formula.
At a source DVR, trace duality sends exponent -b_E to
-(a_E-b_E). The trace dual is again B-projective and S-reflexive;
intersecting all these DVR lattices proves the displayed dual formula
globally. These are algebraic divisor lattices; identifying them with
analytic moment spaces would require its own hypotheses.

Put b_E=floor(a_E/m) and c_E=a_E-b_E. Let H_m be the squarefree product
of target equations C for which at least one prime above C has
a_E-2*b_E>0. The same proof gives

    H_m*T_m^# subset T_m,
    Ann_B(Q_m)=H_m*B,                                    (5)
    Q_m -> Hom_(B/(H_m))(Q_m,B/(H_m)) is an isomorphism.

The pairing in (5) is H_m*Tr(z*w) modulo H_m for z,w in T_m^#.
At a prime above a factor of H_m the needed valuation inequality is

    e_E-c_E >= -b_E,

equivalently e_E>=a_E-2*b_E, which follows from a_E=e_E-1.
Above curves outside H_m the two lattices coincide. Reflexivity gives
the inclusion globally. The linking-form proof of LY0 applies to any
finite projective lattice T_m and its trace dual: lift maps through
F->F/B using projectivity, then use the nondegenerate field trace.
It therefore supplies perfection before the same H_m reduction.

The generic rank on a component C is exactly

    r_(m,C)=sum_(E above C) f_E*(a_E-2*floor(a_E/m)).        (6)

For EVERY m>=3, a-2*floor(a/m)>=a/3>0 when a>0, so H_m=H. Thus all these residual
modules are self-dual maximal Cohen--Macaulay modules over the SAME
reduced branch curve. Their ranks retain the different numerical moment
defects. In local free coordinates for T_m, its regular trace Gram
matrix G_m satisfies

    H_m*G_m^-1 regular, G_m*(H_m*G_m^-1)=H_m*I_d.          (7)

At m=2, H_2 keeps exactly target curves with some even ramification
index. This recovers LY0's parity distinction at the level of exact
annihilators, including curves over which odd and even primes coexist.

In a split tame chart t=u^e, the basis of a nonzero Q_m is

    u^-i,  b+1<=i<=e-1-b, b=floor((e-1)/m).

With H_m=t times a unit, its pairing is that unit times the
anti-diagonal matrix e*I_reverse: only i+j=e contributes. Its rank
e-1-2b is correct for every e,m. Unit twists and unramified residue
extensions must be retained outside a split chart. Higher moments
add hyperbolic pairs in these generic split forms; a Witt-class
vanishing argument still does not show Q_m=0.

## 4. A cusp control: perfect module duality is not fibre perfection

Take the finite normal degree-three cover

    S=C[p,v], B=C[p,q], q=v^3+p*v.

Its Jacobian 3*v^2+p vanishes, so this is not a Keller map.
The source basis 1,v,v^2 has trace Gram matrix

    G = [ 3,    0, -2*p
          0, -2*p,  3*q
       -2*p,  3*q, 2*p^2 ].

Its determinant is -H, with H=4*p^3+27*q^2, an irreducible reduced cusp.
The different is 3*v^2+p, and Lambda is its reciprocal times S.
The matrix N=H*G^-1 is regular and symmetric, and G*N=H*I_3 as in (4).

The module D=coker(G) has generic rank one on the cusp. At its origin
G has rank one, so D/(p,q)D has dimension two. Hence D is not locally
free there. All entries of N vanish at the origin. Thus the induced
bilinear form on this two-dimensional closed-point fibre is zero,
even though D->Hom_(B/(H))(D,B/(H)) is an isomorphism globally.
Base change of Hom fails for this nonprojective curve module.

This directly blocks upgrading (3), by an unsupported tensor-fibre
argument, into a nonsingular vector bundle or into zero defect.
The control is finite and normal on a smooth whole source, but its
Jacobian has zeros and its ramification meets that source. It is not
a counterexample to automatic perfection in the actual Keller setting.

## 5. Consequence for the open trace question

For an actual Keller source, full trace integrality is equivalently
the vanishing of the reduced-curve module D or of its perfect pairing
(2). Sixth-moment trace perfection is Q_3=0. Their exact support scheme,
generic ranks, pairings and local matrix factorizations are now
specified over the reduced branch curve, including its singular points.
The remaining implication must use additional whole-plane Keller
structure to force one of these modules to be zero.

The divisorial lattice identities and all conclusions above are algebraic;
they do not require an analytic L^(2m) characterization, cyclic quotient
singularities, or the differential-module source-layer theorem. The moment
terminology does not assert those separate geometric/analytic hypotheses.

The proofs use the existing tame inverse-different, projectivity,
global linking and moment-lattice theorems. No new classification or
external mathematical theorem is imported. The accompanying exact
controls verify the cusp matrix factorization, singular fibre dimensions,
and split tame forms. Finite controls do not certify the universal
proof or the retained boundary hypotheses. Independent mathematical
review, formal verification and external review remain absent.
