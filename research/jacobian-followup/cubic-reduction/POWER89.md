# An89-variable integral cubic-linear Keller collision

There is an explicit integer matrix \(A\in\operatorname{Mat}_{89}(\mathbb Z)\)
such that

\[
F(Y)=Y+(AY)^{\circ3}
\]

has Jacobian determinant identically one and two distinct rational points
with the same image. Its linear part is the identity, its nonlinear part
is homogeneous cubic, and each nonzero nonlinear coordinate is the cube
of a homogeneous linear form.

The [complete witness](power89/witness.json) contains all7921 integer
entries of \(A\), both89-coordinate points and their common image.
JSON padding and pivot indices are zero-based.
This improves the audited
[90-variable construction](../de-bondt-2610.24889/README.md) associated with
[de Bondt's preprint](https://arxiv.org/abs/2610.24889) to89 variables.
It is an upper bound, not a minimum dimension or priority claim.
The original counterexample and the economical12/18-variable base maps
are not new here.

## An exact89-form shared-cube identity

Let \(K(X)=X+Q(X)+C(X)\) be the verified five-cubic twelve-variable map,
using the source parameter scale11=2. Its cubic part is supported in the
first five outputs. In eighteen variables \(Z=(X,W,\tau)\), define

\[
H(Z)=(\tau Q(X)-\tau^2\iota W,\ C_1(X),\ldots,C_5(X),\ 0),
\]

where \(\iota:\mathbb A^5\to\mathbb A^{12}\) includes the first five
coordinates. Put \(G_\lambda=I+\lambda H\).

The retained decomposition has89 integer linear forms, the rows of
\(L\in\mathbb Z^{89\times18}\), and a rational
\(\overline B\in\mathbb Q^{18\times89}\), satisfying

\[
\overline B(LZ)^{\circ3}=H(Z)+e_{18}\tau^3.
\]

Every coefficient of all eighteen identities is independently checked.
One selected form is exactly \(\tau\); call its index \(p\).
The last row of \(\overline B\) is the coordinate vector selecting \(p\).

The denominator least common multiple in the first seventeen rows is24.
Scale those rows by24, leaving the last row unchanged, to obtain
\(B\in\mathbb Z^{18\times89}\). Let \(L_0\) be \(L\) with row \(p\)
replaced by zero, and set

\[
A=L_0B,\qquad b=Be_p-e_{18}.
\]

Then \(b_{18}=0\), and the cube identity gives

\[
B(L_0Z)^{\circ3}=24H(Z)-b\tau^3.
\]

The right-inverse certificate in the witness verifies \(\operatorname{rank}B=18\).

## Determinant proof

The companion block elimination gives, for every scalar \(\lambda\),

\[
\det DG_\lambda
=\det\!\left(I+\lambda\tau DQ+
                    \lambda^2\tau^2DC\right)
=\det DK(\lambda\tau X)=1.
\]

This uses the exact homogeneous degrees two and three of \(Q,C\).
The twelve-variable base determinant is recomputed independently.

Because \(b_{18}=0\), the map \(S(Z)=Z-b\tau^3\) is an elementary target
shear of determinant one. The previous identity implies

\[
T(Z):=Z+B(L_0Z)^{\circ3}=S(G_{24}(Z)),
\]

so \(\det DT=1\). The shear is retained explicitly; dropping it would
give an incorrect projection formula.

Using \(A=L_0B\), Sylvester's determinant identity now gives

\[
\begin{aligned}
\det DF(Y)
&=\det\!\left(I_{89}
       +3\operatorname{diag}((AY)^2)L_0B\right)\\
&=\det\!\left(I_{18}
       +3B\operatorname{diag}((L_0BY)^2)L_0\right)\\
&=\det DT(BY)=1.
\end{aligned}
\]

The identity \(\det(I+UV)=\det(I+VU)\) follows by eliminating either
identity block in
\(\begin{pmatrix}I&U\\-V&I\end{pmatrix}\), so it applies over the
polynomial coefficient ring. This is a symbolic factorization proof.
An89-variable polynomial determinant is not expanded.

## Explicit rational collision

Let \(p_0,q_0\) be the two verified base collisions of \(K\), and write
\(C_*=(C_1,\ldots,C_5)\). The points

\[
\widehat p=(p_0,-24C_*(p_0),1/24),\qquad
\widehat q=(q_0,-24C_*(q_0),1/24)
\]

satisfy \(G_{24}(\widehat p)=G_{24}(\widehat q)\). The same holds after
the target shear \(S\).

The witness supplies a rational right inverse \(R\) with \(BR=I_{18}\).
Set

\[
u=R\widehat p,\quad v_0=R\widehat q,\quad
\eta=u+(Au)^{\circ3},\quad
v=\eta-(Av_0)^{\circ3}.
\]

The projection identity \(BF(Y)=T(BY)\) yields \(Bv=\widehat q\).
Consequently \(Av=Av_0\), so \(F(v)=\eta=F(u)\).
The points are distinct because their \(B\)-images are distinct.
Every coordinate of this construction is checked in exact rational
arithmetic and retained in the full witness.

## Independent verification and provenance

- [Exact shared-cube candidate](power-basis-anneal/candidate.json)
- [Integer matrix and rational points](power89/witness.json)
- [Primary constructor](build_power89.py)
- [Independent checker](verify_power89.py)
- [Independent receipt](power89/validation.json)
- [Scoped integration receipt](power89/integration.json)

The primary constructor uses SymPy. The independent checker uses only
standard-library Fractions and sparse polynomial arithmetic. It
reconstructs the12-variable base from the pinned input by an explicit
diagonal conjugation, verifies its full polynomial determinant, constructs
\(Q,C,H\), and checks every cube and projection identity. It also verifies
the full integer factorization \(A=L_0B\), the finite right-inverse
certificate and both complete89-coordinate collision evaluations.

The determinant trust boundary is the proved companion block identity
and Sylvester identity, together with these exact coefficient and base
determinant checks. The verification uses exact polynomial identities and the stated determinant
lemmas. No formal proof assistant or external review is claimed.

From the repository root:

~~~sh
.venv/bin/python research/jacobian-followup/cubic-reduction/build_power89.py
.venv/bin/python research/jacobian-followup/cubic-reduction/verify_power89.py
.venv/bin/python research/jacobian-followup/cubic-reduction/check_integration_power89.py
~~~

These commands replay the saved exact identity; discovery need not be
rerun.

## Discovery scope and limits

The finite candidate space contains178 forms obtained from all five
polarization identities on the existing grouped eighteen-variable tensor.
Its cube span has rank121. The original90 forms are retained among them.

A60-second exact linear-arithmetic support search for89 forms timed out
with UNKNOWN; this is not an exclusion. A bounded plateau basis search
stayed at90. The successful
[basis-exchange search](power_basis_anneal.py) allowed temporary support
growth through94 forms before finding89, with seed261115. The
[checkpoint](power-basis-anneal/checkpoint.json) records the run and pins
the source inputs.

The earlier90-form independence result remains valid: the89-form identity
uses a different set and a different coefficient representation. It is
not deletion from that fixed independent list. The search does not prove
optimality even in the finite178-form space, much less over all linear
forms or all Keller maps.

The ordinary cubic12-variable,36-term example remains a separate
construction. Cubic-linear dimension minimization, other special-form
minima and the normalized three-variable support problem remain open.

