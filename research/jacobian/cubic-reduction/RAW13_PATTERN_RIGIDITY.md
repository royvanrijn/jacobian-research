# Rigidity of the established raw13-term collision pattern

The later [weighted tail classification](WEIGHTED_TAIL_CLASSIFICATION.md)
also allows all six extra degree-seven tail atoms with the first output fixed.

The established three-variable degree-seven counterexample has nonlinear
row counts(2,5,6), totaling13, in its original tangent-normalized
coordinates. This is an existing construction, not a new sparsity record.
It is distinct from the ordinary cubic12-variable36-term example and
from the collision-axis-normalized three-variable support interval[10,48].
Those normalized lower bounds must not be transferred to raw coordinates.

Consider all coefficient specializations of precisely these13 nonlinear
monomial atoms, retaining the identity linear part. Require determinant
one and a collision between (0,1,c) and(-1,1,d), where c,d are free.
Over every characteristic-zero field, **the entire family is a
one-parameter diagonal conjugate of the known representative**.
All13 nonlinear coefficients are nonzero. In particular, deleting any
atom cannot give a12-term collision within this support pattern and
collision family. This does not exclude other raw12-term patterns or
other collision configurations.

## Explicit classification

Put v=1+2ab and let

\[
\begin{aligned}
F_1&=a-3a^2b-a^3c,\\
F_2&=3v^2ac+(9v^2-6v-2)b,\\
F_3&=v^3c+(6v^2+2v)b^2.
\end{aligned}
\]

For every nonzero parameter w, define

\[
D_w=\operatorname{diag}(1,1,-1/w),\qquad
F_w(X)=D_w F(D_w^{-1}X).
\]

These and only these maps occur in the stated coefficient/collision
model. Their collision points and common image are

\[
p_w=(0,1,8/w),\qquad q_w=(-1,1,-4/w),\qquad
F_w(p_w)=F_w(q_w)=(0,1,0).
\]

The parameter w is the coefficient of a^3c in the first output.
Every coefficient is the corresponding nonzero coefficient of F times
(-w)^(e_c-delta_(output,3)). Hence all13 atoms remain present.
The determinant and collision are checked as Laurent polynomial
identities in w, not merely at sampled parameter values.

## Normalization and exact certificates

Write B for the first-output a^3c coefficient and A for its a^2b
coefficient. Two explicit unit identities show that neither B nor A
can vanish in the raw coefficient/collision equations.
Conjugate by

\[
D=\operatorname{diag}(1,-3/A,-1/B).
\]

The first two nonlinear coefficients become -1 and -3. The collision
points now share the second coordinate beta=-A/3, which is nonzero;
it must be retained as a variable. The second point's third coordinate
is1+3beta, forced by the first collision equation. Introduce sat with
sat*beta=1.

The normalized determinant and second collision equations have18
generators in13 variables. For each of these13 variables x_i, there is
an explicit membership identity for (x_i-x_i^known)^2 in that ideal.
Over a field, these squares force every normalized coefficient to its
known value, and beta=sat=1. Thus A=-3 in the original collision frame,
and only the scaling of B remains. The third collision equation then
recovers c and d as displayed above.

[Primary normalized calculation](raw13_normalized.py) retains the full
coefficient ideal and its zero-dimensional basis. The earlier
[unscaled boundary searches](raw13_specialization.py) had seven5-second
timeouts; those remain inconclusive artifacts. The normalized route
and its exact certificates close the stated family by a different
calculation, not by relabeling the timeouts.

[Certificate production](certify_raw13_rigidity.py) completed the
Singular calculations. Its generic symbolic parser was stopped during
slow postprocessing, with all completed logs retained. The
[direct sparse parser](parse_raw13_certificates.py) reads those logs
without redoing the algebra. [The certificate manifest](raw13-rigidity/certification-direct.json)
pins the two unit certificates and thirteen squared-coordinate
memberships, totaling94,324 multiplier monomials with maximum degree12.

[The independent checker](verify_raw13_rigidity.py) uses standard-library
Fraction polynomial arithmetic. It reconstructs the original map, its
full three-variable determinant and collision, the universal13-atom
determinant equations, and the normalized equations with the transformed
collision coordinate beta. It then multiplies out all15 certificates
and verifies the complete Laurent family. [Validation](raw13-rigidity/validation.json).
No formal verification, external review, global support optimum or new
counterexample construction is claimed.
