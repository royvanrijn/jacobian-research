# Canonical energy integrality and power traces

Claim HC4MYKO3. This is an all-degree inverse criterion, not an automatic
integrality theorem. The universal HC4 application USES the retained Tang
classification HCTANG2. The preceding KO0/1/2 criteria apply directly to any
Keller map without that classification.

## 1. Every cotangent Keller lift, in every even dimension

Let K:C^m -> C^m be any complex polynomial Keller map, let D be any polynomial,
and put h(z,w)=D(z)+w dot K(z). Its gradient outputs are

    eta=K(z), xi=grad(D)(z)+DK(z)^t*w.

Because det(DK) is a scalar unit, (z,w) -> (z,xi) is a polynomial source
automorphism. The gradient field extension is L0(xi)/F0(xi), where
L0=C(z), F0=C(K), so its degree is d=[L0:F0]. Its source energy is

    G=(z,w) dot grad(h)-h = xi dot z-D(z).

The following conditions are equivalent:

1. K, and hence grad(h), has a polynomial inverse.
2. G is integral over the gradient target C[eta,xi].
3. Tr_(L0(xi)/F0(xi))(G^k) is in C[eta,xi] for 1<=k<=d-1.

For d=1 the trace list is empty; the retained birational Keller theorem
already supplies the inverse.

Proof for d>=2. Condition 1 gives equality of source and target rings, hence
2. Integrality gives 3 by normality of the polynomial target. Since the
variables xi are algebraically independent over L0, field trace extends
coefficientwise, and

    [xi_i^k] Tr((xi dot z-D)^k) = Tr_(L0/F0)(z_i^k).

Thus 3 supplies the first d-1 power traces of EVERY plane/base source
coordinate. KO1 forces the Keller inverse, giving 3 -> 1. This argument
works for all m, all D, all map degrees and all extension degrees. No
finite-free normalization basis or monogenic order is assumed. It does
not replace power traces by ordinary topological periods.

## 2. The energy has a pole at every omitted normalization divisor

Let S0 normalize C[eta] in L0. The gradient-target normalization is
S0[xi]: it is finite and normal, and has the required full fraction field.
In the source coordinates (z,xi), its source open is X0 times A^m_xi.
Every omitted prime divisor is therefore D_i times A^m_xi for a base
omitted divisor D_i. Some z_j has negative valuation there by KO1. In

    G=xi dot z-D(z)

the independent xi variables prevent cancellation of its leading negative
valuation: the minimum of the valuations of D and all z_j is negative,
and its leading residue is a nonzero affine-linear polynomial in xi.
Consequently G has a pole on EVERY omitted prime divisor. It is a
canonical boundary detector, without choosing a generic linear direction.
This statement is invariant under the constant linear cotangent frame.

## 3. Every actual HC4 gradient

For a four-variable polynomial h with constant nonzero Hessian determinant,
the retained HCTANG2 alternative says either grad(h) is already invertible,
or an invertible constant linear source change gives the cotangent form of
Section 1 with m=2. Under h_A(v)=h(A*v), its energy is G_h(A*v), and its
gradient target is obtained by the invertible linear output change A^t.
Thus integrality and regular power traces are unchanged. Affine translations,
if used in normalization, only add a target polynomial to G; binomial
expansion makes the same trace criterion invariant under this addition.

It follows for EVERY HC4 polynomial that

    grad(h) invertible
      iff its canonical energy G is integral over C[grad(h)]
      iff Tr(G^k) is polynomial in grad(h) for 1<=k<=d-1.

Under hypothetical nonautomorphy G has a pole on EVERY omitted divisor
of the full finite normalization. In particular its primitive implicit
equation cannot have a scalar-unit leading coefficient in G. A monic
energy equation would force an inverse EVEN IF the finite order it defines
were nonnormal. This is stronger for this specific observable than asking
the primitive image to be normal; it does not impose monicity on arbitrary
primitive helpers or identify their orders with the normalization.

The retained IV1 energy primitivity is compatible with this conclusion but
is not needed for the trace-coefficient proof. In general a source
polynomial may be integral in a noninvertible Keller map, as the actual
three-dimensional KO2 control demonstrates. The special energy statement
for all HC4 polynomials uses Tang, not arbitrary-dimensional Keller
etaleness. For specified cotangent lifts it needs no classification.

## 4. What remains missing

No proof that G is integral, no monic energy equation, and no automatic
regularity of the listed traces is supplied. Tang's form is compatible with
every hypothetical nonautomorphic plane Keller pair, whose energy then
necessarily has the poles above. This equivalence gives a single intrinsic
observable for the remaining global gate, not a resolution of JC2 or HC4.
The earlier nontorsion-canonical and degree-six packets are not excluded.

The [exact checker](verify_canonical_energy_traces.py) checks coefficient
extraction, target-polynomial energy shifts, the cotangent energy identity,
and the generic trace base-extension representation. Universal claims use
the proofs above. [The context](sources/canonical-energy-trace-context.json)
pins retained artifacts without rerunning a registered checker. There are
no new external theorem imports, no independent replay, no formal
verification and no external review.
