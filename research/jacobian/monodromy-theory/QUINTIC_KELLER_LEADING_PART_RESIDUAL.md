# The leading-part alternatives for a quintic Keller map

**Theorem QuinticKellerLeadingPartResidual.** Let \(F\) be a
degree-five three-variable Keller map over a characteristic-zero
field. If the entries of its degree-five part are linearly
independent, then over the algebraic closure and up to invertible
linear source and target changes its degree-five part is
\[
              \bigl(z^5,\ z^3(x^2+yz),\ z(x^2+yz)^2\bigr).
              \tag{1}
\]
If the entries are dependent and \(F\) is noninvertible, every
nonzero target-linear combination which cancels the degree-five
part has degree four. The leading plane curve of each such
quartic combination is nonreduced or has a singularity outside
the ordinary node, cusp and tacnode types.

To prove the first assertion, work over \(\mathbb C\) by a
finitely generated coefficient-field embedding, or over the
algebraic closure of that field for the linear normal forms.
The leading part \(H_5\) has singular Jacobian by the highest
equation of \(\det JF\). [Wang and Yan, Lemma 2.4](https://arxiv.org/pdf/2610.36495v1)
give \(H_5=g\Phi(p,q)\) with \(\deg p=\deg q=m\),
coprime homogeneous arguments, coprime binary outer entries of
common degree \(e\ge2\), and \(\deg g+em=5\).
The linear-argument case \(m=1\) is
[excluded](QUINTIC_KELLER_LINEAR_ARGUMENT_LEADING_EXCLUSION.md).
The only other possibility is \(m=e=2,\ \deg g=1\).
The three outer quadratics are a basis of the binary quadratic
space, so a target change gives \(H_5=\ell N(p,q)\),
where \(N(u,v)=(u^2,uv,v^2)\) and \(\ell\) is linear.

The [conic pencil argument](QUINTIC_KELLER_CONIC_PENCIL_FIRST_INTEGRALS.md)
excludes a pencil with a nonsingular member unless \(\ell^2\)
belongs to the pencil. It also proves that a pencil with only
singular members is binary. A binary pencil with \(\ell\) in
its binary span is excluded by the
[independent binary theorem](QUINTIC_KELLER_INDEPENDENT_BINARY_LEADING_EXCLUSION.md).
An independent \(\ell\) in the binary-pencil case is excluded
by the [external linear factor theorem](QUINTIC_KELLER_EXTERNAL_BINARY_CONIC_PENCIL_EXCLUSION.md).

It follows that the pencil has a nonsingular member and contains
\(\ell^2\). Change basis to \(p=\ell^2=z^2\). The restriction
of \(q\) to \(z=0\) has rank one or two, since rank zero would
contradict coprimality. Rank two is excluded by the
[formal slice obstruction](QUINTIC_KELLER_CONIC_PENCIL_FORMAL_SLICE.md).
If the rank is one, change the binary coordinates to write
\(q=x^2+z(a x+b y+c z)\). Completing the square in \(x\) and
adding a multiple of \(p\) to \(q\) removes the \(xz,z^2\)
terms. If \(b=0\) the pencil is binary; otherwise rescale
\(y\) to obtain \(q=x^2+yz\) and (1). All these transformations
use coefficients in the algebraic closure of the original
field.

In the dependent case the [cubic component theorem](QUINTIC_KELLER_CUBIC_COMPONENT.md)
makes \(F\) invertible whenever there is a nonconstant
target-linear combination of degree at most three.
A nontrivial linear combination of the coordinates cannot
be constant because the full Jacobian is invertible. Thus a
combination canceling the degree-five part must have degree
four in a noninvertible case. The
[quartic boundary obstruction](KELLER_QUARTIC_COMPONENT_BOUNDARY_OBSTRUCTION.md)
excludes reduced leading curves with only the three stated
singularity types.

The [lower determinant equations](QUINTIC_KELLER_TERNARY_PENCIL_TRANSVERSE_PARAMETER.md#independent-degree-five-leading-parts)
exclude every Keller completion of (1). The [degree-five inversion theorem](THREE_VARIABLE_KELLER_INVERSION_THROUGH_DEGREE_FIVE.md)
completes the dependent case.
For dependent leading parts, the
[refined quartic boundary analysis](KELLER_QUARTIC_TRIPLE_BOUNDARY_AND_LINEAR_COMPONENTS.md)
adds ordinary triple points and transverse repeated support to the
exclusions. It also proves invertibility in the stated squarefree
affine-linear quartic component class.

The [refined boundary tests](KELLER_QUARTIC_QUADRATIC_SUSPENSIONS_AND_HIGHER_CONTACT.md#necessary-singular-point-in-a-reduced-leading-curve)
constrain the projective quadratic jet at a required singular point of
a reduced leading quartic. The [affine-linear classification](KELLER_LINEAR_SOURCE_COEFFICIENT_CUBIC_CLASSIFICATION.md)
and [pure cyclic completion](KELLER_PURE_CYCLIC_SUSPENSION_COMPLETION.md)
give additional component tests. In top rank two the
[quartic pencil argument](QUINTIC_KELLER_RANK_TWO_QUARTIC_PENCIL.md)
and [cubic-factor closure](QUINTIC_KELLER_RANK_TWO_CUBIC_PENCIL_CLOSURE.md)
leave only leading factors of degree at most two in a noninvertible case.
