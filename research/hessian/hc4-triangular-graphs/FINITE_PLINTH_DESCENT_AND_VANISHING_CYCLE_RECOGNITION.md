# Finite plinth descent and vanishing-cycle recognition

Canonical claims: HC4MYKE0, HC4MYKE1, HC4MYKE2.
[MATH_STATUS.json](../../MATH_STATUS.json) is the authority. For every
remaining native frontier, radicality of the descended plinth is now
a consequence of the known finite vanishing cycles, not another action
hypothesis. The exact embedded kernel and descent of the entire plinth
ideal still have to be constructed for one action.

## 1. An arbitrary descended finite ideal is locally a complete intersection (KE0)

Let X be a smooth irreducible factorial affine complex fourfold, R its
ring, and delta a nonzero LND. Assume

    A=ker delta=B[h], B=C[b1,b2],
    pl(delta)=I A,  0<dim_C(B/I)<infinity.                (1)

The three displayed generators of A are algebraically independent and
the equality is of actual embedded subalgebras. The ideal I is proper.
It is not assumed reduced, locally a complete intersection, or generated
by two elements. The full ideal in (1), not only its radical, is required
to be extended from B. Let Z={b_i:1<=i<=n} be its reduced support.

Move Z by the actual JV plane interpolation to {(0,i)}. Write the
resulting coordinates as x,z, and C(z)=product(z-i). For some m,
x^m and C^m lie in I. Choose global local slices

    delta(p)=C^m, delta(q)=-x^m,
    xi=[p/C^m+q/x^m] in H1((A2 minus Z) times A1_h,O).

The kernel is factorially closed, as proved in KD0. Thus x,C are
coprime in the UFD R and R=R_x intersect R_C. The local-slice charts
and gluing argument from KD consequently give

    ann_A(xi)=pl(delta)=I A.                            (2)

This argument uses only the two powers in I, not generators of I:
a coefficient annihilates xi exactly when local primitives with
that delta-derivative glue to a global source function. The free
quotient over (A2 minus Z) times A1_h is the resulting actual Ga torsor.

Here is an elementary local inverse-system argument for the scheme I.
Put K=C(h). At b_i use parameters u=x,v=z-i and localize (2) there
after scalar extension to K. The class xi_i has a finite representative

    xi_i=sum_(a,b>=1) c_ab u^(-a) v^(-b),  c_ab in K.   (3)

Indeed expand the factors (z-j)^(-m), j!=i, at z=i and keep only
negative powers of u,v in the Cech quotient. Only finitely many terms
survive. Passing to formal parameters loses no annihilator information:
the class is killed by fixed powers of u,v, so both its relevant
coefficient module and the finite quotient by its annihilator are
unchanged by completion. Localization of the annihilator of this one
element is its localized annihilator. Hence ann(xi_i)=I_i over K.

On the finite local algebra T_i=K[u,v]_(u,v)/I_i define

    lambda_i(f)=coefficient_(u^-1 v^-1)(f xi_i).

The multiplication pairing (f,g) -> lambda_i(fg) is nondegenerate.
If f is nonzero in T_i, then f xi_i has a nonzero coefficient at some
u^-a v^-b. Multiplication by u^(a-1)v^(b-1) moves exactly that coefficient
to the residue position. Thus f pairs nontrivially with a monomial.
This is the classical finite Macaulay-duality mechanism; a direct
proof in the needed Cech model has just been given.

The socle of T_i is the orthogonal complement of its maximal ideal.
One inclusion is immediate. For the other, if lambda_i(s m)=0 for
every m in that ideal, then lambda_i(s m t)=0 for every t, so
nondegeneracy gives s m=0. Consequently the socle has dimension one.

One can deduce two generators without importing a classification of
Gorenstein ideals. Let L=K[u,v]_(u,v) and T=L/J any finite local quotient.
The Koszul complex of u,v with coefficients T has vector dimensions
(length T,2 length T,length T), H0=K and H2=Soc(T). Since u,v form
a regular sequence in L, this complex computes Tor^L(K,T). The exact
sequence 0 -> J -> L -> T identifies H1 with J/(u,v)J. Its Euler
characteristic is zero, proving

    minimal_number_of_generators(J)=1+dim_K Soc(T).    (4)

For T_i this is two. The dimension of I_(b_i)/m_i I_(b_i) is unchanged
by the flat extension C -> C(h), with the corresponding local rings.
It was already two over C. Nakayama then supplies generators F_i,G_i
in a principal base neighbourhood W_i of b_i containing no other
support point. They form a regular sequence, since their ideal has
height two in the regular plane. Thus I is locally a complete
intersection. No global two-generator assertion is needed here.

On W_i the full source has the presentation

    R_(W_i)=B_(W_i)[h,P,Q]/(F_i P+G_i Q-a_i(b,h)).     (5)

To verify it without assuming I reduced, use global local slices in
R_(W_i) with derivatives G_i,-F_i and their invariant numerator a_i.
The Cech colon identity gives I:a_i=I. Its residue

    alpha_i(h)=a_i(b_i,h)

is a nonzero polynomial. Otherwise a nonzero socle element of the
Artinian local algebra B_(b_i)/I would annihilate a_i, contradicting
that colon equality. The hypersurface in (5) is a domain and is smooth
off F_i=G_i=0. Its removed locus has dimension at most two, since
alpha_i is nonzero. It is therefore Cohen--Macaulay and regular in
codimension one, hence normal. Both source and presentation are the
intersection of their F_i and G_i localizations, proving (5).
This is the KD argument with codimension-two normality; no invertible
Jacobian of (F_i,G_i) at b_i was used.

## 2. Nonreduced holes have no isolated finite Morse cycle (KE1)

Continue (1), and suppose h has a finite critical set. At a support
point b_i let r_i be the rank of D_b(F_i,G_i)(b_i), and let c be a root
of alpha_i. In (5) put L=F_i P+G_i Q-a_i(b,h). A critical point of h
on the smooth hypersurface requires dh to be proportional to dL.
Thus it must have

    (a_i)_h(b_i,c)!=0,
    P dF_i+Q dG_i=(a_i)_b at (b_i,c).                 (6)

If the h derivative is zero, dh cannot be proportional to dL at a
smooth point. If it is nonzero, the second equation is a linear
system of rank r_i in P,Q. For r_i<2, either there is no solution
or an affine family of dimension at least one. The latter contradicts
the finiteness of the critical set. Hence a nonreduced support point,
where r_i<2, supplies no finite critical point at all.

At a reduced point r_i=2. Every root of alpha_i is simple: a repeated
root would make (5) singular at the unique solution of its base
derivative equations. Each root supplies one finite critical point.
In the local source chart (b1,b2,P,Q), its Hessian has mixed block
D_b(F_i,G_i)/(a_i)_h and zero vertical block, so

    det Hess(h)=J_b(F_i,G_i)(b_i)^2/(a_i)_h(b_i,c)^4 !=0. (7)

All finite critical points are therefore Morse. There are none off
the preimage of Z, since h is a coordinate in every local-slice chart.

Choose an open value set S excluding all roots of all alpha_i. The
entire family over S is the actual Ga torsor

    X_S -> (A2 minus Z) times S.                       (8)

The local presentations show that no source point over b_i occurs
at these values; the slice charts cover all other quotient points.
In particular every fibre in (8) has

    H3(X_s,Z)=Z^n,

identified with the primitive links e_i of the points b_i. These
identifications are over the fixed base A2 minus Z, so the local
system on S is constant.

At a reduced point and a simple root c, subtract analytic multiples
of F_i,G_i from P,Q on the nearby h=s fibres. Their local equation
becomes

    F_i P'+G_i Q'=alpha_i(s).

The right side is a local parameter at s=c. This is the ordinary
quadric smoothing, and its Morse vanishing sphere projects with
degree one to the link of b_i. This is exactly KA0's primitive local
sphere/link construction; a globally constant coefficient Jacobian
is not needed. Transport through (8) keeps that class equal to
plus or minus e_i, including collisions of critical values. Extra
excluded values of S can be avoided when choosing the paths.

Let V_fin be the subgroup of H3(X_s,Z) generated by all finite Morse
vanishing cycles transported to one generic fibre. Let

    g=count{i : I is reduced at b_i and deg alpha_i>0}.

Then the preceding description proves the exact integral statement

    V_fin=direct_sum_(such i) Z e_i,
    H3(X_s,Z)/V_fin=Z^(n-g).                           (9)

There is no hidden torsion or index in (9). In particular

    codim_Q span(V_fin) >= number of nonreduced support points. (10)

If the finite vanishing cycles span the full generic H3 even over Q,
then every support point is reduced, and each has a nonconstant
alpha_i. Thus the entire ideal I is radical. The theorem does not
assume reducedness, local complete intersection, or a generating pair
for I. The local algebra in Section 1 and finite-cycle exhaustion
force those conclusions.

If in addition X=A4, or more generally H3(X)=H4(X)=0, KD now applies
to the radical ideal. It forces exactly one reached value on every
line, constructs the global simple cotangent frame and its polynomial
inverse, and in the general case recognizes X as A4. No assertion
about multiplicity-weighted nontransverse normal maps is needed.

## 3. A genuine counterexample to automatic radicality from Morse points alone

In actual A4 coordinates x,z,P,Q put

    f=x(x+z), g=z(z-1),
    h=x+P f+Q g,
    delta=g partial_P-f partial_Q.                    (11)

The actual kernel is C[x,z,h], and its exact plinth is
I C[x,z,h], I=(f,g). For the latter equality, the Cech numerator is
h-x; it is monic in h, hence a nonzerodivisor on (B/I)[h]. The
regular-sequence colon calculation gives ann[(h-x)/(fg)]=I A,
and factorial source intersection identifies that annihilator with
the entire plinth. Coprimality gives the kernel as in KC.

The support consists of (0,0),(0,1),(-1,1). At the first point I
has length two and is not reduced; the other two points are simple.
The full critical ideal is exactly

    (z-1, x(x+1), P+2x+1, Q+x).

There are two Morse points, with values 0,-1 and critical Hessian
determinants both one. The nonreduced point supplies no critical
point because the x derivative of h is one there. Nevertheless,
the generic fibre is a Ga torsor over the three-punctured plane,
so H3=Z^3. Its two finite cycles span only the two simple-point
links, and their quotient is Z. Thus finite Morse points, a common
critical Hessian scalar, actual affine-space source and descended
finite plinth do not alone imply radicality. The finite-cycle
exhaustion used in Section 2 is essential. This example is not a
counterexample to HC4; its generic Betti count already fails the
necessary simple Keller-cotangent count for its two critical points.

## 4. The native action criterion no longer needs a radicality test (KE2)

The remaining normalized native functions have the finite-cycle
exhaustion required above: KA supplies an integral basis of finite
Morse cycles for h_D, K2 and every S_rho, including all critical-value
collisions. This property is invariant under actual polynomial
source equivalence. Therefore, for each such h, the following are
equivalent:

1. h has an actual simple cotangent frame and polynomial inverse.
2. One nonzero h-preserving LND has actual embedded kernel B[h],
   B a polynomial plane with the three generators independent,
   and exact plinth I B[h] for any proper finite-colength ideal
   I in B, with no reducedness or complete-intersection assumption.

For 2 -> 1, Sections 1--2 force I radical and KD supplies the global
frame. For 1 -> 2, the ordinary cotangent translation has the exact
radical plinth, so it meets the weaker requirement as well. The native
critical count forces the support size: two for K2/h1 and three for
S_rho/R_rho. The coordinate/unit-ideal branch stays separate.

Native common critical Hessian scalars and KB then make condition 2
equivalent to the explicit source alignment h~R_(N,c). Unit-Hessian
source existence remains exactly E_N AND condition 2. For every rho,
S_rho~R_rho is thus equivalent to this finite-ideal descent criterion;
its unit existence still separately requires E3. K2~h1 has the
analogous two-point statement. Neither action nor alignment is
constructed here.

The exact plinth equality in the same embedded B[h] remains a real
condition. A plinth whose radical has finite cylindrical support need
not itself be extended from B; this proof does not cover it. Nor does
it make B[h] an automatic kernel, make h a coordinate in an arbitrary
polynomial kernel, or straighten moving quotient holes. Collective
flexibility supplies none of these individual-action statements.
E2/E3, both endpoint inversion requirements and unrestricted HC4/JC2
remain open.

A local Cech control illustrates why full ideal descent cannot simply
be replaced by cylindrical support. For

    eta=h/(u^2 v)+1/(u v^2)

in C[h][u^-1,v^-1] with nonnegative powers discarded, its exact
annihilator is (u-hv,v^2), whereas its contraction to C[u,v] is
(u,v)^2. Its radical is the constant cylinder (u,v), but extension
of that contraction loses the linear relation u-hv. After using
u=hv, the two surviving Laurent classes are independent over C[h],
which verifies the annihilator equality. This is a local-cohomology
countercontrol, not a claimed native affine-space action or a proof
that moving action data can never be straightened.

## Sources and verification

[Khovanskii--Monin, Section 2.1](https://arxiv.org/html/2106.15562v1)
provides primary context for the classical annihilator/duality
construction. The finite Laurent-residue proof and Koszul generator
calculation used here are written in Section 1; no duality or
classification theorem is imported without proof.

[Masuda, Section 3](https://arxiv.org/html/2312.05455v1) was inspected
at Theorems 3.8/3.9 and Corollary 3.11 as a comparison. That source
assumes a principal plinth in the modification theorem, or smooth,
surjective equidimensional quotient hypotheses in its triviality
corollary. Those hypotheses do not supply (1). Its principal-divisor
conclusion is not transferred to our codimension-two plinth. No
restriction to h-fibres is made without a proved kernel/base change.

All structural inputs come from the unchanged KD, KA, KB, JV, JU
and KC evidence. [The checker](verify_finite_plinth_descent.py) checks
residue pairings and Koszul controls, local critical-rank and Morse
identities, the actual nonreduced negative example, primitive cycle
lattices and gate logic. It is not a search over degrees. [The context](sources/finite-plinth-descent-context.json)
pins prior inputs and identifies both primary comparisons as context
only; [the receipt](receipt-finite-plinth-descent.json) records the
implementation controls. Independent replay, formal verification
and external review are not claimed. The original goal is not complete.
