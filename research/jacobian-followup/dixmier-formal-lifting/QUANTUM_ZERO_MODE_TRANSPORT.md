# Quantum zero modes and canonical transport of a cokernel obstruction

15 November 2026. **QZM1.** A zero mode gives a coordinate-independent
cokernel obstruction to a missing quantum action. Normalized field trace
makes this obstruction survive every finite algebraic coefficient cover.
Transporting the representation as well as the operator preserves it
under actual finite algebra automorphisms and identity stabilization.

For the quadratic chart [a,v]=h/v, this excludes
S=v*(lambda*a+mu)+G(a) for every nonconstant polynomial G and
lambda nonzero, with arbitrary finite partner orders and algebraic
coefficient covers. Combined with QVL1, a separable operator
S=v*H(a)+G(a) with both H and G nonconstant has no finite mate.
The remaining pure-polynomial family has H constant; arbitrary negative
v corrections and deformed first images are not excluded.

The obstruction concerns a specified quantum operator and module.
It is not determined by the classical reduction alone, and it is not
being identified with QRF2's restricted order-seven cocycle.

## The vector space cokernel records a necessary commutator condition

Let A be an algebra over a characteristic-zero field K, V a left
A-module, and L,b in A. Suppose e in V satisfies
\[
       Le=0,\qquad
           \omega(L,V,e;b)=[be]\in V/LV .                    \tag{1}
\]
The quotient here is a K-vector-space cokernel; LV need not be an
A-submodule. If an operator T in A satisfies [L,T]=b, then
\[
                         be=L(Te),\qquad \omega=0.           \tag{2}
\]
Thus a nonzero class excludes every such T at once. No order,
coefficient-degree or matrix-size search is involved in this identity.

In particular take L=S-c and b=kappa, with kappa a nonzero
scalar constant. A Weyl mate [S,T]=kappa is excluded provided
[c,T]=0 and [e] is nonzero. The commutation with c is automatic
when c is scalar. When c is a function of another prescribed
image R, it is part of the full mixed-relation requirement.

## Trace makes the class survive every finite coefficient cover

Let F be a characteristic-zero differential field and E/F a finite
algebraic extension. Derivations extend uniquely. Each field embedding
into a normal closure therefore commutes with those derivations, so
normalized trace commutes with every finite differential operator over F:
\[
   \tau=\frac1{[E:F]}\operatorname{Tr}_{E/F},\qquad
          \tau P=P\tau,\qquad \tau|_F=\operatorname{id}.       \tag{3}
\]
Consequently F/PF -> E/PE is injective.

More generally, for any differential F-module V, the componentwise
map tau tensor id retracts E tensor_F V onto V and commutes with
the base-extended action. It therefore retracts the corresponding
vector-space cokernel. A nonzero class over F cannot acquire a
primitive on a finite algebraic coefficient cover.

For an inhomogeneous equation Pf=g over F, any algebraic solution
would yield the rational/base-field solution tau(f).
This is a linear statement. It does not assert that trace preserves
nonlinear quantum correction equations, specified spatial supports,
or a prescribed integral lattice.

## Transport the module together with the canonical change

For an actual algebra automorphism Phi:A->A, define V_Phi to have
the same vector space as V with action
a*v=Phi^(-1)(a)v. Then Phi(L) acts exactly as L on that vector
space and Phi(b)e acts as be. Thus
\[
 \omega(\Phi(L),V_\Phi,e;\Phi(b))
                          =\omega(L,V,e;b).                 \tag{4}
\]
For a finite Weyl automorphism this transports the class through
Fourier changes and polynomial Weyl shears. It is essential to
transport the module: after Fourier transformation the appropriate
module need not be the usual coefficient-function representation.
Changing the coefficient localization is a separate operation and
can kill a supported vector. In the Fourier control below, the new
position x annihilates e, so inverting that position would kill e.
The trace statement concerns finite extensions of an already valid
differential coefficient field; it does not justify this additional
localization.

Identity stabilization uses V tensor_K W for any nonzero module
W of the spectator Weyl algebra. Since L acts only on the first
factor, its cokernel is (V/LV) tensor_K W. A nonzero class remains
nonzero. Further actual finite canonical automorphisms, even ones
mixing the old and spectator pairs, are covered by (4).

This statement requires an actual algebra automorphism defined on
the algebra in question. A classical Poisson change without a proved
finite quantum action is not silently assumed to act. If an operator
belongs only to a rational coefficient localization, an arbitrary
polynomial Weyl automorphism need not preserve that localization.
For globally admitted operators, use the ambient Weyl module and
clear coefficient denominators as below.

## A cyclic module gives a deformation interpretation

For a domain differential-operator algebra D and L=S-c nonzero,
put M=D/DL and e=1 modulo DL. Right multiplication by L gives
the free left-D resolution
\[
 0\longrightarrow D\xrightarrow{\ \cdot L\ }D
             \longrightarrow M\longrightarrow0,\qquad
             \operatorname{Ext}^1_D(M,M)=M/LM .              \tag{5}
\]
The class [e] is the first-order eigenvalue-deformation class:
the relation (L-epsilon)e_epsilon=0 over epsilon^2=0 gives
L*e_epsilon=epsilon*e_epsilon. A finite commutator mate makes
that class zero by (2).

The standard relation between self-extensions and first-order
connection deformations is discussed in
[Arinkin, Section4.4](https://arxiv.org/html/0808.0699).
The particular identification (5) follows directly from the displayed
resolution. No rigidity or global lifting theorem is imported.

This is a necessary obstruction, not a complete mate criterion.
Vanishing of the class supplies neither a finite commutator identity
nor all the other right Weyl actions.

## A complete all-degree linear coefficient exclusion

Work over C(v), where C has characteristic zero, h is a nonzero
constant and a acts as (h/v)*partial_v. Let
\[
 S=v(\lambda a+\mu)+G(a),\qquad
 \lambda\ne0,\qquad G\in C[a],\quad\deg G=m\ge1 .             \tag{6}
\]
Then S has no finite scalar T with [S,T]=kappa a nonzero
constant, over any finite algebraic function-field extension.

Put alpha=-mu/lambda and apply the coefficient-preserving Ore
automorphism a->a+alpha. Subtract c=G(alpha). The transformed
operator is
\[
 P=\lambda va+\sum_{j=1}^m c_j a^j,\qquad c_m\ne0,\qquad
                              P1=0 .                       \tag{7}
\]
A mate would therefore give an algebraic coefficient function
\[
                        f=T1,\qquad Pf=\kappa.              \tag{8}
\]
By (3) it suffices to exclude rational f in C(v).

If m=1, equation (8) is
df=(kappa/h)*v/(lambda*v+c1) dv. At v=-c1/lambda its residue is
\[
                         -\frac{\kappa c_1}{h\lambda^2}
                                      \ne0 .                \tag{9}
\]
A rational derivative has zero residues, a contradiction.
The trace argument has already covered every finite algebraic cover.

Suppose m>=2. Away from v=0 the ordinary differential leading
coefficient of P is c_m*h^m*v^(-m), which is nonzero.
A rational pole of f there would give a uniquely highest pole
from the m-th derivative; it cannot match a constant forcing.
At v=0, if the first pole term is b*v^beta with beta<0, use
\[
 a^j v^\beta=
 h^j\beta(\beta-2)\cdots(\beta-2j+2)v^{\beta-2j}.              \tag{10}
\]
The j=m term has uniquely lowest exponent and nonzero coefficient.
No pole at zero is possible either. Thus f must be a polynomial.

For a nonconstant polynomial f of degree N, lambda*v*a lowers
its degree by one and every a^j term lowers it by at least two.
Equation (8) forces N=1 and f=(kappa/(lambda*h))*v+constant.
Substitution leaves
\[
 Pf-\kappa=
 \frac{\kappa}{\lambda h}
  \sum_{j=1}^m c_jh^j
       \prod_{t=0}^{j-1}(1-2t)\,v^{1-2j}\ne0 .               \tag{11}
\]
The j=m term is nonzero and has strictly lowest v power.
This proves [1] nonzero in C(v)/PC(v), and proves the theorem
for every polynomial degree m and every finite mate order.

The result also includes arbitrary constant-field dependence,
such as coefficients in K(r)((h)). The shift alpha can have an
h pole; no regular formal normalization is inferred.

## The remaining separable polynomial class becomes smaller

[QVL1](ALGEBRAIC_INFINITY_LIMITS.md) excludes
v*H(a)+G(a) whenever deg H>=2, regardless of G, because G's
coefficients are bounded at infinity. Equation (6) excludes
degree-one H whenever G is nonconstant. Hence
\[
  H,G\text{ both nonconstant}\quad\Longrightarrow\quad
       vH(a)+G(a)\text{ has no finite Weyl mate}.             \tag{12}
\]
For nonconstant H this gives a complete classification:
the sole positive case is degree-one H with G constant,
where T=v/lambda gives [S,T]=h.
The degree and nonconstancy conditions are therefore material.

In particular, a bounded full exact-R scalar Long lift of the
separable form v*H_h(a)+G_h(a), with H0=2 and
G0=4a-3ra^2, must keep H_h independent of a at generic h.
Adding arbitrarily high finite momentum powers to that coefficient
cannot supply the missing action.

The leading-linear rank-four example
S=v*(2+h*a)+4a-3ra^2+h^2*a^3 is excluded by (6).
Even tuning G_h'(alpha)=0 does not repair it. For example
G_h=4a-3ra^2-(r*h+h^2/3)*a^3 has derivative zero at alpha=-2/h,
but its shifted operator still has the nonzero top-pole term
in (11). The vanishing of that first derivative does not make the
cokernel class exact.

The family with H constant and nonlinear G_h remains outside
this new theorem. Its quadratic case is covered by QSO2;
higher generic orders and additional negative-v coefficient
corrections are not settled here.

## Relative parameters and globally admitted canonical pairs

The auxiliary value c may depend on a prescribed exact first image
R=r. The rational chart has [a,q]=[d,r]=h, q=v^2/2.
For alpha=alpha(r), the coefficient representation on K(r,v)
extends to the full chart by
\[
 a\mapsto(h/v)\partial_v+\alpha(r),\qquad
 d\mapsto h\partial_r+q\alpha_r .                             \tag{13}
\]
The mixed commutator vanishes because
alpha_r=(1/v)*partial_v(q*alpha_r).

The action of the original global generators uses the
[retained rational canonical chart](QUINTIC_PHASE_AND_BRANCH_OBSTRUCTION.md).
Suppose R,S are globally admitted finite Weyl operators and their
localized form is (6), with coefficients rational in r.
Write c(r)=g(r)/f(r), with nonzero polynomial f and polynomial g
over the generic scalar field. Use the finite polynomial operator
\[
             L=f(R)S-g(R),\qquad b=h f(R).                   \tag{14}
\]
Its coefficient zero mode is e=1. A full partner satisfying
[S,T]=h and [R,T]=0 would have [L,T]=b.
But its class [b e] is nonzero: multiplication by f(r) is
invertible in the coefficient field and the cokernel computation
is linear over that field of differential constants.

Consequently this actual quantum pair obstruction survives
every finite target Weyl automorphism and identity stabilization
by (4). This conditional global statement requires admission of
R,S themselves. It does not claim that every rational-chart
formula is globally polynomial, nor normalize every deformed R.

## Controls and scope

A positive case is G constant, with [v*(lambda*a+mu),v/lambda]=h.
For the canonical transport control take L=xD in the usual
coefficient module, [D,x]=h. The actual residue functional
ell(f)=Res_(x=0)(f*dx/x) has ell(Lf)=0 and ell(1)=1,
so it certifies the nonzero class. The same functional on the
transported underlying vector space continues to certify it.
Fourier transformation gives
L'=-xD-h. Its transported module has x acting as -D_old and D
as x_old, so L' acts exactly as x_old*D_old and still kills 1.
In the unchanged coefficient module L' would not kill 1.
A polynomial Weyl shear gives another exact transport control.

[Exact controls](check_zero_mode_transport.py) check trace and
differentiation on finite rational-function covers, the universal
monomial action, shifted candidate operators, the highest-pole
defects, the degree-one residue, the positive mate, Fourier/shear
module transport, a compatible parameter connection and clearing
of its auxiliary denominator. The all-degree pole/degree
exclusion, arbitrary finite covers, cokernel transport and
stabilization use the written proofs.

The nonzero class concerns an actual quantum operator on an
actual module. Changing its quantum correction, adding terms
outside (6), or changing a first image can change that class.
No identification with QRF2's cocycle or a theorem about every
Long representative is asserted. No finite-left Long carrier,
bounded full Long lift, scalar or matrix counterexample, independent
theorem replay, formal verification, external review, novelty
priority or DC2 solution is claimed.
