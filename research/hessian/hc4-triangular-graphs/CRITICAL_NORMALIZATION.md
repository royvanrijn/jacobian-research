# The entire finite polar locus is excluded; the remaining data live at the boundary

For `y=R(x,p)` in the original Meng--Yang family, let `L*N != 0` and
allow every actual polynomial automorphism of the retained three duals
over K[x]. Degrees and word lengths are unrestricted. Work first over C;
the same nonexistence statements descend to any characteristic-zero field
by embedding the finitely generated coefficient field of a candidate and
its dual inverse into C.

Put

\[
 F=3+2xR,\qquad
 D=-\frac{Np}{4x^2}-\frac{(M-13N/(2x^2))^2}{4L}.
\]

The [exact critical-curve formula](RELATIVE_CRITICAL_LOCUS.md) gives, on
every p-dependent component of F having multiplicity one,

\[
 K_f=-32LN^2x^6R_p^2\frac{d^2D}{dx^2}.                         \tag{1}
\]

These are precisely the components contributing to the critical locus
on `x*R_p != 0`; vertical and repeated factors contribute none there.
An allowed unit-Hessian pullback requires the **same nonzero constant**
K_f on every such component, by HC4MYGK0. The polynomial identity
holds in the component's function field, hence also on its normalization
at points outside the initial rational computation chart.

**HC4MYGM0/1. Under this hypothesis, R_p has no zero at any finite
point of any contributing component.** In particular such components
are smooth, mutually disjoint, disjoint from every other F component,
and etale over C* under x. The whole finite singular and polar locus is
excluded in every degree, extending [HC4MYGN3](NEWTON_CRITICAL_GATE.md)
from smooth ramification to singular branches and component intersections.

## 1. The local input is normalization, not a classification

We explicitly import two standard local statements from the
[retained primary sources](sources/critical-normalization-import.json):
Garcia Barroso--Ploski, [arXiv:1208.0913v1](https://arxiv.org/pdf/1208.0913v1),
Theorems 1.6/1.7; and Hefez--Rodrigues--Salomao,
[arXiv:1708.07412v1](https://arxiv.org/html/1708.07412v1), Theorem 2.1.
The former provides primitive normalization parameters and intersection
orders; the latter gives, for a local irreducible branch h,

\[
              \operatorname{ord}_t h_p=c(h)+\operatorname{ord}_t x_t. \tag{2}
\]

Here c is the conductor exponent of the branch ring in C[[t]]. It is
zero for a smooth branch and at least two for a singular branch:
if c were one, the conductor would contain t, making the branch ring
equal C[[t]] and smooth. No catalogue of singularities is assumed.

At a branch appearing with multiplicity one in F, write `F=h*g` in the
completed local ring, absorbing units into g. Let `I=ord_t g>=0`, the
weighted intersection contribution of all other local factors.
On h=0, `F_p=h_p*g`. If x0 is nonzero and `x=x0+t^e`, (2) gives

\[
                       j:=\operatorname{ord}_t R_p=e-1+c+I.      \tag{3}
\]

Each contributing component is nonvertical, so e is finite and positive.

## 2. All finite cases (HC4MYGM0)

Let `(x0,p0)` be any finite point of a contributing component. Necessarily
x0 is nonzero, since F(0,p)=3. Suppose R_p vanishes there, so j>0.

If the branch is smooth and e=1, x is a local parameter. D and D''
are regular at the point. Equation (1) has positive order at least 2j,
contradicting a nonzero constant.

If the branch is smooth and e>=2, then p-p0 is a local parameter.
With a primitive normalization parameter t,
`p=p0+beta*t+...`, beta nonzero. Since x-x0 has order e,

\[
 D=D_0-\frac{N\beta}{4x_0^2}t+O(t^2),\qquad
 \operatorname{ord}_t D''=1-2e.
\]

Here c=0, and (1)/(3) give

\[
                       \operatorname{ord}_t K_f=2I-1.            \tag{4}
\]

This odd integer is never zero. With I=0 it is the earlier simple pole;
with additional intersecting components it is a positive odd zero.

Finally suppose the branch is singular, so c>=2. The function D is
regular in t. Unless constant, its first nonconstant term has order at
least one, so `ord D'' >= 1-2e`; a zero derivative gives an immediate
contradiction to (1) as well. Thus

\[
                       \operatorname{ord}_t K_f\ge2c+2I-1\ge3.    \tag{5}
\]

Again a nonzero constant is impossible. These cases include arbitrary
singular multiplicity, arbitrary intersection order, and repeated other
components. A limiting point with R_p=0 is not mislabeled a finite
relative critical point: the obstruction is the valuation of (1) along
the genuine open critical component.

## 3. A polynomial Bezout test and the whole critical scheme (HC4MYGM1)

Let S be the product of the distinct p-dependent irreducible factors of F
that have multiplicity **one** in F. Choose their primitive polynomial
representatives in C[x,p]; p-independent factors and all higher
multiplicities are excluded from S. If there are none, set S=1.

Section 2 proves the single necessary condition

\[
                 (S,F_p)=(1)\quad\hbox{in }\mathbb C[x,p].        \tag{6}
\]

Equivalently there must be polynomials U,V with `U*S+V*F_p=1`.
The full F_p is essential: replacing it by S_p would lose the
intersection contributions of repeated and vertical other components.
This is a necessary polynomial test for all degrees, not a sufficient
Hessian condition or a proposed bounded search.

On the affine curve S=0, x is a unit, with explicit inverse `-2R/3`,
because S divides F. Condition (6) also makes F_p and hence R_p units.
It follows that every component is smooth and the x projection is etale;
F_p nonzero also shows that all these components are disjoint from
one another and every other F component. The projection need **not**
be finite: branches may escape to p=infinity at finite nonzero x.

The native critical formulas

\[
 T_c=\frac{13N/(2x^2)-M}{2L},\qquad
 q_c=-\frac{3T_c}{4x}+\frac1{16x^2R_p},\qquad
 r_c=\frac{(1+xR)^3p+3x(1+xR)^2q_c-T_c}{x^3}
\]

now define regular functions on the **whole** S=0. They give an
isomorphism from S=0 onto the entire open relative critical scheme,
with polynomial restrictions in both directions: replace the inverses
using (6) and `x^-1=-2R/3`. No omitted finite normalization points or
localized holes remain.

For the remaining nonautomorphic boundary class, R0 is nonzero with no
simple root and `R0(-M/(2L)) != 0`, so no critical points lie at x=0.
The entire relative critical scheme is therefore precisely this closed
smooth plane curve S=0. It must be nonempty by
[HC4MYGN0](NEWTON_CRITICAL_GATE.md), so S cannot equal one.

## 4. Exact escape law at finite nonzero base values (HC4MYGM2)

Let n=deg_p F=deg_p R, and set

\[
                         H(x,w)=w^nF(x,1/w).
\]

Take a normalization point of a contributing component with x tending
to x0 nonzero and p having a pole. Write `x=x0+t^e`, and
`w=1/p=eta*t^m+...`, eta nonzero, e,m positive. The active branch has
multiplicity one in H. Let c be its conductor in this reciprocal chart
and I its weighted intersection with **all other factors of H**, including
vertical or repeated factors if present. Both are nonnegative.

D has exact pole order m, so D'' has exact pole order m+2e, with
nonzero coefficient. From (1), a unit K_f requires

\[
                   \operatorname{ord}_t R_p=e+m/2.               \tag{7}
\]

In particular **m is even**. On H=0 direct differentiation gives
`F_p=-w^(2-n)*H_w`. Applying (2) in the reciprocal chart yields

\[
 \operatorname{ord}_t R_p=(2-n)m+c+e-1+I.
\]

Combining this with (7) gives the exact necessary equality

\[
                         2(c+I)=(2n-3)m+2.                     \tag{8}
\]

This includes all ramification indices e; e cancels, rather than being
assumed one. A smooth isolated reciprocal point has c=I=0: it is
impossible for n>=2, and for n=1 would require m=2. A reduced ordinary
node has c=0,I=1 on either branch, making (8) impossible for every
integer n>=1. Thus for n>=2 an admissible escape requires a boundary
singularity worse than a smooth point or ordinary node, with the
specified even pole order and conductor/intersection budget. In fact
`c+I >= 2n-2` follows immediately from m>=2.

These necessary conditions are not sufficient. A local cusp can meet
them. For example the reciprocal germ `(w-z)^2-z^3=0`, with
`z=t^2,w=t^2+t^3`, has e=m=2,c=2,I=0 and meets (8) when n=2.
Our checker also retains an actual original-family control of this
type whose critical Hessian has a nonzero limit at the escape but is
not constant. Passing this local test does not produce a unit Hessian.

## 5. The residual plane Keller-fibre obligation (HC4MYGM3)

Assume positive p degree and the nonautomorphic boundary conditions
from Section 3, and suppose an allowed pullback had unit Hessian.
Its gradient cannot be an automorphism: its critical curve is either
empty, or x is a nonconstant unit on each S component, so it cannot
be A1. Emptiness is already excluded by HC4MYGN0.

The generic dual value fibre has Euler1+r>=2. Hence the rank-zero
Tang form is impossible, as in [HC4MYGR1](CRITICAL_COMPONENTS.md).
The rank-one form supplies a plane Keller pair (A,C), a plane
polynomial E, and an algebraic isomorphism

\[
                         S=0\ \simeq\ \{C=0\}.
\]

Under this isomorphism the actual preserved critical target-line
function `tau=dD/dx` is A, while `E=D-x*tau` on the curve.
Indeed the rank-one normal form is `E(U,Z)+x A(U,Z)+V C(U,Z)`;
its critical equations solve for x,V, and its base derivative is A.
Consequently `dE=-x*dA` on C=0.

This is a necessary **abstract curve and function realization**, not
a conclusion that `(tau,S)` is a polynomial Keller pair in the
original (x,p) plane. A curve isomorphism is not silently extended
to a plane automorphism. Nor do (6), (8), or this realization alone
prove that a retained-dual automorphism producing unit Hessian exists.

## Verification and current frontier

[verify_critical_normalization.py](verify_critical_normalization.py) checks
the valuation identities on smooth, intersecting and singular branches,
the reciprocal differentiation factor, the even-pole law, exact Bezout
and whole-curve inverse controls, and the admissible-but-insufficient
cusp escape. [Replay receipt](receipt-critical-normalization.json).
The all-degree exhaustion and normalization arguments are written proofs;
the two local inputs are explicit byte-pinned imports, not independently
proved by the finite controls. Tang is used only in Section 5's residual
realization. No formal verification, independent external review or
novelty priority is asserted.

Finite affine singularities, component intersections and ramification are
no longer remaining loopholes. A candidate must instead supply the
nonempty S and polynomial Bezout certificate (6), satisfy the earlier
Newton constraints and (8) at every finite-base escape, and realize
the prescribed curve functions as an entire plane Keller fibre.
Later splitting at repeated initial faces and more complicated boundary
singularities have not been exhausted. The unit-ideal plane-surjectivity
and distinct singleton-zero properness endpoints, q/r-dependent omitted
bases, genuine free-base/dual mixing and general HC4/JC2 remain open.
