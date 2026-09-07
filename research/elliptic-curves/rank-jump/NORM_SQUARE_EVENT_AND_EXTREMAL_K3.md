# A square norm reduces the torsion field but adds no generic directions

In the fixed-cubic pencil, the point-independent condition
\[
D(u)=N(1-u\theta)\in\mathbf Q^{\times2}
\]
does not produce a new generic Mordell–Weil block. The quadratic twist by
\(D(u)\) is an extremal K3 surface with geometric generic rank zero.
Consequently adjoining \(\sqrt{D(u)}\) adds no independent generic
Mordell–Weil direction. This rank statement does not settle finite-index
changes in the section lattice.

There is also an exact coefficient-defined control:
\[
\boxed{u_*=-A/B,\qquad D(u_*)=1,\qquad
 1-u_*\theta=-\theta^3/B.}
\]
For the retained anchor, this element has two negative real conjugates.
It is not square even in the totally real cubic splitting field.
The Jacobian two-torsion field has degree 24 and its module is still
indecomposable. A square norm therefore does not imply that the three
branch square roots split together.

These conclusions concern generic rank and torsion structure. The
specialized rank at \(u_*\) is **UNKNOWN**. No point search or descent was
performed on that fibre, and no individual norm-square fibre is excluded
as a high-rank curve.

## The K3 twist supplies an exact generic-rank bound

Use
\[
E_u:y^2=x^3+a_2x^2+a_4x+a_6,
\]
where
\[
\begin{aligned}
a_2&=2Au,\\
a_4&=A+3Bu+A^2u^2,\\
a_6&=B+ABu^2-B^2u^3.
\end{aligned}
\]
Assume \(B\ne0\) and \(\delta=-4A^3-27B^2\ne0\). The
[reassessment](../notes/RANK_JUMP_REASSESSMENT_2026-09-05.md) proves that
\(E_u\) has three geometric \(I_2\) fibres at \(D=0\), an \(I_0^*\)
fibre at infinity, and geometric generic rank one.

Its \(D\)-twist is the explicit Weierstrass equation
\[
E_u^{(D)}:
y^2=x^3+D a_2x^2+D^2a_4x+D^3a_6.
\]
Direct symbolic calculation gives
\[
\begin{aligned}
\Delta(E_u^{(D)})&=16\delta D^8,\\
c_4(E_u^{(D)})&=
 16D^2(A^2u^2-9Bu-3A).
\end{aligned}
\]
The discriminant of \(D\) is \(\delta\), and its roots are disjoint
from the second factor of \(c_4\). Each finite bad fibre therefore has
valuations \((v(c_4),v(\Delta))=(2,8)\), hence type \(I_2^*\).
This also follows by twisting \(I_2\) once at each simple root of \(D\).

The coefficient degrees are \(4,8,12\). At infinity put
\(v=1/u,\ X=v^4x,\ Y=v^6y\). The transformed discriminant has constant
term \(16\delta B^8\ne0\); infinity is good. Thus the minimal elliptic
surface over \(\mathbf P^1\) has exactly
\[
3I_2^*,\qquad e=24,\qquad \chi(\mathcal O)=2.
\]
It is a K3 surface. The three root lattices are \(D_6\), so its trivial
lattice already has rank
\[
2+3\cdot6=20.
\]
The K3 Picard bound and Shioda–Tate formula force
\[
\boxed{\rho=20,\qquad
\operatorname{rank}E_u^{(D)}(\overline{\mathbf Q}(u))=0.}
\]
These are the standard twisting, fibre and Picard tools in
[Schütt–Shioda, Sections 5, 6 and 12](https://arxiv.org/pdf/0907.0298).
The application to the displayed pencil is the calculation above.

The finite invariant check does not mistake a high Picard number for
Mordell–Weil directions: all twenty divisor dimensions are already
accounted for by the zero section, fibre and reducible components.

## The norm double cover has no new character space

Let \(\mathcal C_D\) be the smooth genus-one curve
\[
v^2=D(u).
\]
It has the rational point \((u,v)=(0,1)\). The usual invariant and
anti-invariant decomposition over its quadratic function-field extension
gives, after tensoring with \(\mathbf Q\),
\[
E\bigl(\mathbf Q(\mathcal C_D)\bigr)\otimes\mathbf Q
 \simeq
 \left(E(\mathbf Q(u))\otimes\mathbf Q\right)
 \oplus
 \left(E^{(D)}(\mathbf Q(u))\otimes\mathbf Q\right).
\]
The second summand is zero even geometrically. Hence
\[
\operatorname{rank}E(\overline{\mathbf Q}(\mathcal C_D))=1.
\]
For the retained nonsquare \(B\), the original arithmetic generic rank
is zero, and therefore
\[
\boxed{\operatorname{rank}E(\mathbf Q(\mathcal C_D))=0.}
\]

This closes one proposed implication: the norm-square base change has
**no new generic rank contribution**. It does not say that its individual
rational fibres have rank zero. In particular the rank-at-least-20
anchor at \(u=0\) remains a specialization of this generic-rank-zero family.

## A coefficient-only norm-square control

Writing
\[
x_D=Bu,\qquad y_D=Bv
\]
identifies \(\mathcal C_D\) with
\[
y_D^2=x_D^3+A x_D^2+B^2.
\]
On this auxiliary curve, \(Q=(0,B)\) has tangent slope zero, giving
\[
2Q=(-A,-B).
\]
Its base parameter is precisely \(u_*=-A/B\). This is a deterministic
algebraic construction from the anchor coefficients, not a parameter
search. The auxiliary point is not a point on \(E_{u_*}\).

The cubic relation yields
\[
1+\frac AB\theta=-\frac{\theta^3}{B}.
\]
For the retained anchor, the certificate verifies \(B>0,\ \delta>0\).
The three distinct real roots have one negative and two positive signs:
their product is \(-B<0\), and their sum is zero. The displayed
\(\gamma_*=1-u_*\theta\) consequently has one positive and two negative
conjugates. It cannot be square in the totally real splitting field \(L\).

Nevertheless \(N(\gamma_*)=D(u_*)=1\). The three conjugate squareclasses
have the product relation but are not all zero. Their quotient of the
three-coordinate permutation module by the all-ones relation is the
irreducible two-dimensional \(\mathbf F_2[S_3]\)-module. They therefore
have exactly two independent classes over \(L\):
\[
[L(\sqrt{\gamma_{*,1}},\sqrt{\gamma_{*,2}},\sqrt{\gamma_{*,3}}):L]=4.
\]
Since the retained cubic has Galois group \(S_3\), the full degree is 24.

## Norm collapse is weaker than splitting the extension

On the norm-square locus, the signs in the
[six-branch-point model](TORSION_DIFFERENCE_AND_CT.md) have even parity.
The resulting signed-permutation group is
\[
(\mathbf Z/2)^2\rtimes S_3\simeq S_4.
\]
The new certificate checks all 24 actions on all 16 Jacobian
two-torsion vectors. Its commuting algebra is still
\(\mathbf F_2[\epsilon]/(\epsilon^2)\), whose only idempotents are zero
and one. Thus this degree-24 action has an indecomposable, nonsplit
two-torsion extension.

For comparison, \(u=0\) has \(D=1\) and all three \(\gamma_i=1\).
That is the separately proved
[split-correspondence event](LINEAR_TWIST_SOLUBLE_BLOCKS.md), which does
transport the known rational block. The two controls have the same
norm-square result but different branch splitting:

| Control | \(D\) | Independent branch squareclasses over \(L\) | Extension |
|---|---:|---:|---|
| \(u=0\) | 1 | 0 | split stable limit |
| \(u=-A/B\) | 1 | 2 | nonsplit degree-24 torsion action |

The second row does not determine a CT value or rational solubility.
It disproves the intermediate inference “square norm implies split
branch extension,” independently of that missing arithmetic.

## Mechanism ranking and next gate

1. **Solubility:** the full labelled extension and its cup obstruction
   remain stronger candidates than the norm. The
   [scalar calculation](INDEPENDENT_SCALAR_CUP_AND_TWIST_BLOCKS.md)
   already computes a nontrivial part independently.
2. **Weak incidence proxy:** a square \(D\) reduces a torsion-field
   character but creates no generic Mordell–Weil character space on
   its double cover. It is not an automatic source of independent points.
3. **Missing computation:** the original non-scalar \(1+\theta\) cup
   evaluation remains open. The present theorem does not compute it.
4. **Missing positive implication:** a useful block condition must
   produce rational classes beyond a norm identity or a smaller
   Galois group.

For Agent 1, the useful information is a theorem gate: do not justify a
norm-square selector by claiming a generic rank gain or complete branch
splitting. Its possible statistical effect on exceptional specializations
is unmeasured, so this result supplies no veto on individual candidates.
No selector or search policy was changed.

The next productive arithmetic target is still an independent non-scalar
cup value or a construction solving several full covers together.
Increasing a norm-square parameter budget would not resolve either gap.

## Replay and scope

The [protocol](NORM_SQUARE_RANK_PROTOCOL.json) freezes the structural
control and finite calculation. The
[certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_norm_square_rank_v1.json)
retains the coefficient identities, sign argument, rank accounting and
all 384 finite root-action checks.

```sh
python3 elliptic-curves/rank-jump/norm_square_rank.py check
sage -python elliptic-curves/rank-jump/verify_norm_square_rank.py
python3 -m unittest discover -s elliptic-curves/rank-jump -p test_norm_square_rank.py
```

The Sage verifier independently checks the universal invariants, absence
of a singular fibre at infinity, irreducibility and real signs of the
retained cubic, and the auxiliary-curve doubling identity. No class group,
point search, new-fibre rank calculation or Selmer computation is invoked.
Active-search files and mathematical-status entries remain unchanged.
