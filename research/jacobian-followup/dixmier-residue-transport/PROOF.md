# Strong quantum residue transport through polynomial Weyl shears

14 November 2026. Work over a characteristic-zero field. A strong restricted
residue transports through every positional polynomial Weyl shear, provided
its correction spaces, affine lower lifts and output functional are transported
together. The retained order-five residue has an exact symbolic example
under a nonlinear shear. Reusing the old Bernstein bounds in the new chart
fails: one admissible native correction pairs to -1296 lambda with the
transported functional.

The transport theorem is functoriality under an explicit algebra automorphism;
no novelty priority is claimed. Its application and support comparison are
recorded as QRT1. It does not give an unrestricted obstruction for DC2.

## The source residue and its domain

Use the retained seed (kappa,tau)=(0,1), or (a,tau)=(-1/2,1), and the
canonical pair S,T in A=K[X,Q,Z], with

\[
\delta=3X^2\partial_X+(2-6XQ)\partial_Q,\qquad
\Pi=\partial_Z\otimes\delta-\delta\otimes\partial_Z,\qquad
\{S,T\}=1.
\]

The Weyl product is m exp(hbar Pi/2). The two displayed derivations commute.
The [retained strong residue](../../archive/non-elliptic/extended-geometry/RESTRICTED_QUANTUM_DEFORMATION_COCYCLES.md)
is a 16-term coefficient functional Lambda. In the declared parity branch,
write the attempted lifts as S+hbar^2 s+hbar^4 u and T+hbar^2 t+hbar^4 v.
The order-five defect before the current correction is

\[
O_5(s,t)=\{s,t\}+\Pi^3(s,T)/24+\Pi^3(S,t)/24+\Pi^5(S,T)/1920.
\]

The current linear differential is d1(u,v)={u,T}+{S,v}. Its source
Bernstein weights are (1,1,3) on X,Q,Z. The correction bounds are
(degree,Z-order)=(25,3),(21,2) at order two, and (21,1),(17,0) at order
four. The retained proof also allows all finite negative X powers under
those bounds. It proves

\[
\Lambda d_1=0,\qquad
\Lambda O_5(s,t)=C=-47547660815739/190658\ne0
\tag{1}
\]

on every allowed lower lift, including Hamiltonian gauge. The new replay
checks the complete bounded polynomial lower-lift space of dimension 42;
the full finite-Laurent strength is an input from the retained proof.

## Exact quantum transport for arbitrary polynomial potentials

Let H belong to K[X,Q], independent of Z and hbar, and define

\[
\alpha_H=\exp(-\operatorname{ad}_{\star}H/\hbar).
\]

The star commutator is divisible by hbar. Its normalized negative adjoint is

\[
D_H=\sum_{r\ge0}\frac{\hbar^{2r}\delta^{2r+1}H}
 {2^{2r}(2r+1)!}\partial_Z^{2r+1}.
\tag{2}
\]

All coefficients are independent of Z, so the summands commute and every
application of D_H lowers Z degree. The exponential is therefore finite
on each polynomial, and on each finite-X-Laurent polynomial. It is an
exact star-algebra automorphism, with inverse alpha_(-H). No convergence
or truncation hypothesis is needed for this assertion.

Set g=delta H, tau_H(f)=f(X,Q,Z+g), A2=delta^3 H/24 and
A4=delta^5 H/1920. Then

\[
\alpha_H=\tau_H\left(1+\hbar^2 A2\partial_Z^3+
\hbar^4(A4\partial_Z^5+\tfrac12 A2^2\partial_Z^6)+\cdots\right).
\tag{3}
\]

Thus the classical images are S'=tau_H S,T'=tau_H T, but their lower
corrections are affine transforms

\[
s'=\tau_H(s+A2 S_{ZZZ}),\qquad
t'=\tau_H(t+A2 T_{ZZZ}).
\tag{4}
\]

The order-four correction space above this lower lift is the affine space

\[
(u',v')=(q_S,q_T)+(\tau_H u,\tau_H v),
\tag{5}
\]

where

\[
q_S=\tau_H(A2 s_{ZZZ}+A4 S_{Z^5}+\tfrac12 A2^2 S_{Z^6}),
\]

and q_T has the analogous expression. Omitting these affine shifts is
not the exact quantum coordinate change.

The Poisson map tau_H gives chain isomorphisms on the transported
Hamiltonian, correction and defect spaces:

\[
d'_0\tau_H=\tau_Hd_0,\qquad
d'_1(\tau_H u,\tau_H v)=\tau_Hd_1(u,v).
\tag{6}
\]

Define Lambda_H=Lambda composed with tau_H inverse. Applying alpha_H to
the star commutator and taking its order-five coefficient, after the
order-three equation is satisfied, yields

\[
O'_5(s',t')+d'_1(q_S,q_T)=\tau_H O_5(s,t).
\tag{7}
\]

Consequently Lambda_H annihilates every transported current variation,
and its value on the affine corrected defect in (7) is C, independent
of all transported lower lifts and gauges. This proves strong residue
transport for every H in the stated polynomial class, not just a finite
set of shear degrees or coefficients. These shears form an additive group;
arbitrary finite words in this same positional class are covered as well.

For the boundary valuation nu(X,Q,Z)=(1,-1,-2), the subclass H in K[X]
preserves the valuation filtration in both directions. Each nonconstant
term delta^(2r+1)H has X valuation at least 2r+2, while differentiation
in Z raises valuation by two per derivative. Every nonidentity term in
(2) strictly raises valuation. This does not imply preservation of the
separate Bernstein degree bounds.

## Ambient Weyl algebra and all six relations

In the second Weyl algebra take positions X,Q and momenta p_x,p_q, with
[p_x,X]=[p_q,Q]=hbar. The Ore generator

\[
Z=3X^2p_x+2(1-3XQ)p_q
\]

has commutators [Z,f]=hbar delta f with every position polynomial f.
Its normal and symmetric ambient expressions coincide because delta has
zero divergence. The exact Weyl shear induced by the same H is

\[
p_x\mapsto p_x+H_X,\qquad p_q\mapsto p_q+H_Q,
\qquad X,Q\text{ fixed}.
\]

It preserves all six ambient Weyl relations and induces Z -> Z+delta H.
Formula (3) uses the fiber Weyl symbol convention of the retained residue;
it is not plain commutative substitution in that convention.

More generally, apply this source algebra automorphism to any declared
four-image deformation and its six commutator defects. Every defect is
carried to its alpha_H image componentwise, including the coefficient
mixing in (4)--(7). Extend a functional on the necessary S,T relation by
zero on the other five components. It annihilates the image of the full
transported correction differential whenever the S,T projections obey
the source restrictions. Any solution of all six relations would solve
this necessary pair relation and contradict its nonzero period. This
does not construct a complete four-image endomorphism, nor allow
unrestricted S,T corrections outside the transported category.

## Symbolic nonlinear example and a native correction failure

Choose H=lambda X^3/3. In ambient variables this is the nonlinear shear
p_x -> p_x+lambda X^2. In the fiber convention,

\[
g=3\lambda X^4,\qquad
A2=\tfrac{45}{2}\lambda X^6,\qquad
A4=\tfrac{1701}{16}\lambda X^8.
\]

For a defect of Z degree at most four, the complete transported dual is
given by the following finite formula. If Lambda has coefficients c_(i,j,k),
then

\[
\Lambda_H(X^I Q^J Z^K)=
\sum_{k=0}^K {K\choose k}(-3\lambda)^{K-k}
c_{(I+4(K-k),J,k)}.
\tag{8}
\]

The saved table has 60 terms over K[lambda], including possible negative
X exponents. The bound K<=4 covers every order-five defect of this
declared lower family: S,T have Z degrees five,four and s,t have degrees
at most three,two. It is not a cutoff on the original X-Laurent pole order.

All 614 current polynomial columns and all 945 nonconstant coefficients
of the complete affine-42 lower-lift family pass the symbolic transport
test. Direct expansion verifies (7) at the selected lower lift and at a
second point differing by a nonzero kernel vector. Thirty-five gauge
controls and all six ambient shear relations also pass; the universal
identities above give the infinite-class coverage.

The transported current space differs from the old numerical box. For
example,

\[
\tau_H(X^{18}Z)=X^{18}Z+3\lambda X^{22},
\]

whose second term exceeds the old S4 bound 21. Conversely the native
monomial m=X^12 Q^6 Z has Bernstein degree 21, but

\[
\tau_H^{-1}m=m-3\lambda X^{16}Q^6.
\]

The source functional kills d1(m,0), whereas exact coefficient extraction
gives Lambda(d1(X^16 Q^6,0))=432. Therefore

\[
\boxed{\Lambda_H(d'_1(m,0))=-1296\lambda.}
\tag{9}
\]

For nonzero lambda, this functional is not a cocycle for the unchanged
native correction box. Equation (9) is the unique nonzero pairing among
its 614 polynomial current columns. Nineteen transported source basis
vectors leave that box. This explicitly identifies a correction direction
which a naive support-preserving argument misses.

An additional exact linear solve adjoins X^16 Q^6 to the original 614
current directions. At the selected lower lift, the full order-five
equations remain inconsistent. Killing the scalar period is therefore
not enough to repair that lift. No universal assertion about all enlarged
lower lifts or all possible dual supports follows from this last check.

## Assurance and remaining question

[verify_transport.py](verify_transport.py) retains the symbolic dual table,
input packet, exact counts, hashes and the native failure in
[receipt.json](receipt.json). [check_independent.py](check_independent.py)
uses separate SymPy differentiation, without importing the producer or
legacy arithmetic, to replay the scalar period, native failure, quantum
coefficients, dual coefficient values and six ambient shear relations.
It does not independently replay all 945 producer columns or the original
finite-Laurent strong-residue proof. No formal verification or external
review is claimed.

The next open step is a cocycle for the sum of the transported and native
correction spaces, with its complete lower-lift scheme. The current Lambda_H
cannot supply it by (9). The all-degree transport theorem establishes an
orbit-covariant restricted obstruction; it does not make the restrictions
coordinate independent or settle DC2.
