# Ramification-independent lifting obstructions for HC4

Claims HC4MYKJ0/1/2. [MATH_STATUS.json](../../MATH_STATUS.json) is the
authority. This note addresses a possible route to a counterexample to
full HC4, independently of the native graph-alignment questions. It closes
that route for the characteristic-two plane example and its cotangent
potential, with arbitrary correction degrees and arbitrary ramification.
It does not prove HC4 or exclude other characteristic-zero constructions.

## 1. The candidate and what was already proved

Over any field k of characteristic two, let

    P=x+x^2 y+x^4+x^6 y^2,
    Q=y+x^5+x^6 y+x^7 y^2+x^8 y^3,
    H0(x,y,u,v)=uP+vQ.

The plane map is credited to
[Mondello](https://arxiv.org/html/2608.02634v1), derived from Huq-Kuruvilla's
construction. Its Jacobian is 1 and

    (P,Q)(0,1)=(P,Q)(1,0)=(P,Q)(1,1)=(0,1).

The retained claim HKM2W1 already proves that the plane map has no Keller
lift through Z/4 in any degree, including its polynomial left--right orbit.
We do not replay its checker or claim that unramified obstruction anew.
Its proof is [retained here](../../foundations/verified/HUQ_KURUVILLA_PLANE_W2_OBSTRUCTION.md).

New issues are whether ramified coefficient extensions help, and whether
an arbitrary four-variable potential correction can lift H0 even though
lifting its two coefficient polynomials cannot. The second question is
strictly broader than insisting on cotangent form.

Let O be any complete discrete valuation ring of mixed characteristic
(0,2), uniformizer pi, residue field k, and

    e=val_pi(2)>=1.

No perfection assumption on k is required. Precision r means the
coefficient ring O/(pi^r). For this candidate the exact maximal precisions are

    plane constant-Jacobian lift: r=e;
    four-variable constant-unit-Hessian potential lift: r=2e. (1)

Thus neither admits a characteristic-zero polynomial lift over O, however
ramified it is. These are statements about every polynomial degree, not
successive searches at bounded degree.

## 2. Divided Jacobian coefficients: HC4MYKJ0

For arbitrary polynomials F=sum f_ab x^a y^b and G=sum g_cd x^c y^d
over Z, the coefficient of x^(r-1)y^(s-1) in their Jacobian is

    sum_(a+c=r,b+d=s) (ad-bc) f_ab g_cd
      =sum_(a+c=r,b+d=s) (as-br) f_ab g_cd.            (2)

Every summand is divisible by gcd(r,s). In particular for r=p i,s=p j,

    [x^(pi-1)y^(pj-1)] Jac(F,G)=p B_ij(F,G),
    B_ij=sum_(a+c=pi,b+d=pj) (a j-b i) f_ab g_cd,     (3)

where B_ij is an INTEGER polynomial in the coefficients. Constant-Jacobian
maps over a p-torsion-free ring satisfy B_ij=0. Their reductions satisfy
the reductions of these divided equations, even over a ramified ring.
This is the elementary coefficient form of an integer-saturated Keller
equation; no conjectural sufficiency of such equations is used.

More sharply, over a mixed (0,p) DVR with val_pi(p)=e, a constant-Jacobian
lift modulo pi^(e+1) forces B_ij=0 modulo pi, since p B_ij vanishes there.
For i,j>=1 the coefficient is nonconstant, so the value of the constant
Jacobian does not matter. Extra polynomial degrees cannot alter (2).

For the displayed P,Q, take p=2,i=j=1. Their coefficient vectors give

    B_11(P,Q)=1 in k.                                (4)

Only the pair x^2 y from P and y from Q contributes, with divided
weight 1. Thus there is no plane lift modulo pi^(e+1), over any O above.
For the integer representatives P0,Q0 their Jacobian is 1+2K with K
integral. They do give lifts modulo pi^e. This proves the sharp plane
part of (1), extending HKM2W1 beyond unramified Witt rings.

The same nonlifting holds for every polynomial plane left--right equivalent
representative over k. Indeed the plane tameness theorem retained in
HKM2W1 expresses each automorphism as affine and triangular factors.
Lift their coefficients and nonzero scalars to O; these factors and their
inverses are polynomial over O with constant unit Jacobian. A lift of an
equivalent map would compose to a forbidden lift of P,Q. The obstruction
therefore excludes arbitrary ramified plane rechartings of this orbit.

## 3. A universal cubic--quartic Hessian equation: HC4MYKJ1

Use ordered coordinates z=(x,y,u,v), q=xu+yv. After removing constant
and linear terms, write a potential with this quadratic part as

    H=q+C+D+terms of degree>=5,

where C is a completely general homogeneous cubic and D a completely
general homogeneous quartic. Put

    C=sum_(a+b+c+d=3) c_abcd x^a y^b u^c v^d,
    A=[x^2 y u]D, B=[x y^2 v]D,
    Delta=partial_x partial_u+partial_y partial_v.

The degree-one part of det Hess(H) is 2 Delta C. Write

    Delta C=ell_x x+ell_y y+ell_u u+ell_v v,
    ell_x=2c2010+c1101, ell_y=2c0201+c1110,
    ell_u=2c1020+c0111, ell_v=2c0102+c1011.

Define the following integral quadratic polynomial in the cubic coefficients:

    T=-c0102 c1200+3c0120 c3000-2c0201 c2010
                         +2c0210 c2001+3c0300 c1002-c1020 c2100.

The EXACT universal coefficient identity is

    [xy] det Hess(H)
      =4(A+B-T-ell_u c2100-ell_v c1200)+2ell_x ell_y. (5)

All twenty cubic coefficients are arbitrary. All quartic coefficients
other than A,B drop out of this one coefficient. Higher potential terms
cannot contribute to degree two of the determinant.

Here is a direct derivation. Let M=Hess(q), so M^2=1 and det M=1.
The degree-two determinant term is

    tr(M Hess(D))
      + ((tr(M Hess(C)))^2-tr((M Hess(C))^2))/2.

The first summand's xy coefficient is 4(A+B). In the second, expand the
four-by-four matrix, or its two-by-two horizontal/vertical blocks. Its
xy coefficient is

    -4T-4ell_u c2100-4ell_v c1200+2ell_x ell_y.

This is an identity over Z, including its displayed factors; there is no
division in a residue ring. A short complete polynomial check with all
twenty symbolic cubic coefficients is provided by the checker. It is a
universal jet identity, not a bounded search for potentials.

Over characteristic zero, constant Hessian determinant forces all ell's
zero, giving A+B=T. The divided equation survives reduction whenever the
quadratic part is q. In particular if C reduces to zero, every such
characteristic-zero reduction must have A+B=0.

## 4. Arbitrary ramification and arbitrary potential corrections

Suppose a polynomial H over O/(pi^r) reduces to H0 and has constant unit
Hessian determinant, with r=2e+1. Its quadratic part reduces to q and
its cubic part reduces to zero. First normalize its quadratic part by an
O-linear variable change congruent to the identity modulo pi.

This normalization is valid even in residue characteristic two. Lift the
quadratic part arbitrarily to O. An isotropic residue basis vector lifts
to an isotropic vector by solving Q(e_x+t e_u)=0 with t in (pi); its
derivative in t is a unit, so ordinary Newton iteration converges. Choose
its partner with polar pairing one and replace that partner w by
w-Q(w) v0 to make it isotropic, without dividing by two. Split off the
hyperbolic plane and repeat on its rank-two orthogonal complement.
Every basis change is the identity modulo pi. Thus the whole quadratic
part becomes q. This only multiplies the Hessian determinant by a constant
unit, and preserves the reduction of every homogeneous part.

For the normalized H, its constant Hessian determinant is 1. Its cubic
coefficients belong to (pi), and its linear determinant equations give

    2ell_i=0 modulo pi^r,
    ell_i in (pi^(r-e))=(pi^(e+1)).                   (6)

Since T is quadratic in cubic coefficients, T belongs to (pi^2). Since
H0 has quartic part x^2 y u, A+B is a unit, congruent to 1 modulo pi.
In (5), the error terms involving ell have valuation strictly greater
than r. The main term 4(A+B-T) has valuation EXACTLY 2e, so cannot
vanish modulo pi^(2e+1). This is a contradiction.

This proves nonexistence with arbitrary four-variable corrections at
every polynomial degree, without preserving cotangent form. Conversely
the integer cotangent representative H*=uP0+vQ0 satisfies

    det Hess(H*)=Jac(P0,Q0)^2=(1+2K)^2
                 =1+4K(1+K).

It therefore gives a constant-unit-Hessian lift modulo pi^(2e). This
proves the second sharp threshold in (1). In particular, over the
unramified two-adic integers the potential lifts modulo four but has
NO lift modulo eight, in any potential degree. The plane obstruction
alone would not prove that stronger statement.

More generally the nonlifting proof works for every characteristic-two
potential with quadratic part q, cubic part zero, and quartic coefficients
A+B!=0. It is not a theorem excluding all characteristic-two potentials
or all possible reductions of an HC4 counterexample.

## 5. Why this is a full-HC4 counterexample route, and its exact boundary: HC4MYKJ2

The block Hessian identity gives det Hess(H0)=1 in k. Its gradient at
(x,y,0,0) is (0,0,P,Q). Hence the three distinct source points

    (0,1,0,0), (1,0,0,0), (1,1,0,0)

all have gradient image (0,0,0,1). It is a characteristic-two unit-Hessian
noninjective gradient example. This does not contradict complex HC4.

If a characteristic-zero unit-Hessian lift over O had existed, it would
indeed give an HC4 counterexample over Frac(O). To see why noninvertibility
would follow, translate its gradient's constant output to zero. Its linear
Jacobian has unit determinant in O. The formal inverse is constructed
recursively with coefficients in O, using only that inverse linear matrix.
If the gradient had a polynomial inverse over Frac(O), that polynomial
would coincide with the formal inverse and would consequently be integral.
Reduction would give a polynomial inverse to grad(H0), contradicting the
three-point collision. This argument applies to every polynomial degree.
Passing to the finitely generated characteristic-zero field of its
coefficients and embedding that field in C preserves this conclusion:
the unique formal inverse would descend any polynomial inverse. Such a
lift would therefore also yield a complex HC4 counterexample.

KJ1 now excludes that entire ramified lifting route. It also excludes
affine variable rechartings of H0, multiplication by a nonzero residue
scalar and addition of an affine potential: lift the affine changes and
units over O, use Hessian congruence, and apply the theorem. Arbitrary
NONLINEAR source changes are not included: the global constant-Hessian
condition is not invariant under them. Nonintegral changes without this
reduction are also outside the statement. A different characteristic-two
potential passing the divided equations may still be a candidate.

The full HC4/JC2 decision from HCTANG2 remains open. Neither critical
Gorensteinness in the native quotient programme nor arithmetic support
uniqueness in the quantum programme supplies the missing characteristic-
zero counterexample or universal inverse theorem.

## 6. Literature and assurance

[Maubach--Rauf, Sections2--4](https://arxiv.org/html/1507.02946v1) construct
integer Keller equations by contraction of characteristic-zero equations
and define strong Keller maps. Those definitions and their motivating
coefficient divisions were inspected. We use only the elementary divided
identity (2)--(3), proved here; we import no strong-Keller invertibility or
lifting conjecture. The characteristic-two formula, Jacobian and collisions
were checked directly against Mondello's Sections1--2. Its degree-three
field theorem and formalization are not new theorem inputs here. The
unramified arithmetic work HKM2W1 is retained separately, with all its
original evidence unchanged. There are zero new external proved imports.

The [checker](verify_ramified_hc4_lifting.py),
[receipt](receipt-ramified-hc4-lifting.json), and
[context](sources/ramified-hc4-lifting-context.json) verify universal symbolic
coefficient identities, the cubic--quartic formula and exact ideal
divisibilities, the candidate's leading terms and collisions, and the
valuation thresholds. No old registered arithmetic checker is rerun.
Independent replay, formal verification and external review are not claimed.
