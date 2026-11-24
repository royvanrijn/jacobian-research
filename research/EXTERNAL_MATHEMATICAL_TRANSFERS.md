# External mathematical transfers

The external result catalogue is pinned to
[solved-conjectures main, 758066055ef4](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/CONTENTS.md).
The research notes are pinned to
[research, dab1ac8d1ece](https://github.com/1337isnot1337/solved-conjectures/tree/dab1ac8d1ecea8dc34c38932fcd2a4b3f56ce0ef/research).
The repository's claim IDs refer to [the current mathematical DAG](MATH_STATUS.json).

| External result | Mathematical transfer and exact DAG target |
| --- | --- |
| Hajdu–Sárközy shifted powers and shifted products, Erdős perturbed squares, and Sárközy smooth fixed gaps | The [determinant and correlation audit](arithmetic/SHIFTED_PRODUCT_AND_CORRELATION_TRANSFER_AUDIT.md) reconstructs the four-cycle exclusion, the summable-coefficient and no-log argument, and the mechanism for eliminating the unknown third-order term in a fixed-gap test. The square-full determinant forces the singular case at large prime targets. The six-cycle and pointwise kernel-pair nodes still require independent bounds. For the Gaussian-norm frontier \(\texttt{GAP-ARITH-CK-AVERAGED-CENTERED-TRIPLES}\), pair limits plus the fourth-moment upper sieve yield triple occurrence for windows longer than \((3+\sqrt3)/6\) mean spacings, with the retained ZF and G inputs. The first moment alone, with the same upper sieve, suffices unconditionally beyond two mean spacings. |
| Karp–Vishnyakova–Zhang, Nuttall, Csordas–Vishnyakova and Cardon et al. | The [factorial audit](factorial/image-rigidity/GAMMA_QUOTIENT_TRANSFER_AUDIT.md) independently derives the all-order gamma-quotient signs and a product determinant on unit-spaced columns. It records the scope of the external Lean receipt. The actual unweighted linked beta array has opposite \(2\times2\) minor signs. An exact contiguous subtraction proves all-order consecutive-moment coprimality for the axis-supported monomial-linear family in \(\texttt{GAP-FACTORIAL-BINARY-LINKED-BETA-COPRIMALITY}\). Distinct positive exponents and the uniform Hurwitz determinant obligations remain open. The zero-geometry comparisons include a negative cubic test for the actual gamma multipliers. |
| Sokal Apéry threshold, Hamburger classification and Ferreira–Simon beta-prime convolution | The [Gaussian audit](moments/vanishing/GAUSSIAN_W_AFFINE_ROTATIONAL_DIVERGENCE.md) collapses a rotation orbit of dominant singularities to one actual-kernel singularity. This covers a symmetry slice beyond \(\texttt{GAUSSIAN-W-AFFINE-SELECTED-DOMINANT-DIVERGENCE}\). It also gives an exact positive-measure classification for normalized coefficients of \(W(1+Z^m)Y\), with separate real-line and nonnegative supports. General nonsymmetric dominant-point cancellation remains open. |
| Hecke bootstrap, CM descent, Artin uniform and average estimates, and Goldbach spectral bridge | The [analytic audit](arithmetic/EXTERNAL_HECKE_ARTIN_GOLDBACH_TRANSFER_AUDIT.md) extends the conditional window interface to the actual angular characters and norm twists, checks the Artin weighted family, and gives an unconditional CM square-norm inclusion and fixed-field splitting density. The gcd-stratified Goldbach projection transfers exactly. The operator-norm interface for shrinking windows loses two powers of the relative width and still requires an estimate for the actual weights. The Eisenstein central-window, sharp four-form energy and switching nodes retain their analytic obligations. |
| Ma–Tang finite-group random subset sums | [Local intensities](arithmetic/RANDOM_PRIME_SEMIGROUP_LOCAL_INTENSITIES.md) proves an unconditional strict improvement of the critical gap-density lower bound and an explicit parity obstruction to any residue-blind subcritical pointwise limit. The exact Walsh certificate in the pinned source uses a binary affine fibre. In the semigroup model the remaining conditions for a critical density are an empirical intensity law and a joint-support bound. |
| Féray conjugate-power words | The [Hurwitz census](arithmetic/HURWITZ_EXACT_CENSUS_AND_REMOVABLE_FACTORS.md) uses conditional class Fourier calculations to give exact absolute totals and to count directly removable last factors, with limiting fraction \(e^{-2}\). The few-cycle mixing hypothesis fails strongly for double transpositions, whose two-factor product stays at total variation \(1-o(1)\) from the parity-uniform law. The exact census targets \(\texttt{GAP-ARITH-HURWITZ-UNIFORM-REMOVABLE-DT-REPRESENTATIVE}\) and the uniform fixed-product fibre connectivity nodes. |
| Stone nonproper-locus obstruction | The degree criterion and the polynomial comparison below locate the missing hypothesis for \(\texttt{GOAL-JC2}\) and \(\texttt{GAP-JACOBIAN-1}\). A known three-dimensional Keller map gives an exact obstruction to inferring generic uniqueness from even a full target line of singleton fibres. |

## Further consequences and remaining estimates

- In the pure and marked mixed Hurwitz domains, the removable-factor
  count has an exact finite probability law and a normal limit. The
  exceptional fraction decays faster than every fixed inverse power.
  For each fixed \(C>0\), any orbit with relative size at least \(r^{-C}\)
  therefore meets the removable locus for all sufficiently large \(r\).
  The unresolved step for \(\texttt{GAP-ARITH-HURWITZ-UNIFORM-REMOVABLE-DT-REPRESENTATIVE}\)
  concerns orbits smaller than these thresholds. [Joint law and orbit criterion](arithmetic/HURWITZ_EXACT_CENSUS_AND_REMOVABLE_FACTORS.md#the-joint-removable-factor-law)
- In native increasing row and column order, a positive-exponent
  raw linked beta matrix is weakly sign regular exactly when \(u/r\)
  is a positive integer. The arbitrary-minor zero criterion gives
  an explicit test for which minors survive at each order.
  The physical moment polynomials also contain binomial row
  factors, and the remaining linked-beta coprimality problem
  requires control of their common complex roots. On the axis,
  the shifted-shape theorem covers \(\Gamma(1)\times\Gamma(\beta)\)
  for every \(\beta>0\).
  [Sign classification and shifted-shape proof](factorial/image-rigidity/GAMMA_QUOTIENT_TRANSFER_AUDIT.md)
- The finite-modulus critical gap bounds have a sharp truncation
  asymptotic along \(\operatorname{lcm}(1,\ldots,\lfloor y\rfloor)\),
  and the all-powers-of-two bound is explicit. The remaining
  critical-density route requires an empirical conditional intensity
  law and control of overlapping unexposed prime supports. The
  quantitative truncation concerns the proposed law itself.
  [Truncation and density criterion](arithmetic/RANDOM_PRIME_SEMIGROUP_LOCAL_INTENSITIES.md#strict-two-adic-improvement-and-sharp-finite-modulus-error)
- Gaussian divergence survives when the inverse equation rotates
  and the numerator breaks rotation symmetry. The full numerator
  separates into the \(m\) residue classes of exponents and each
  component has the same quotient singularity majorant.
  The next case is an inverse with tied dominant points and no such
  rotation of its defining equation.
  [Actual source proof](moments/vanishing/GAUSSIAN_W_AFFINE_ROTATIONAL_DIVERGENCE.md)
- Optimizing the centered Goldbach Hardy test shape has sharp
  relative norm cost \(e^2\pi/(32\Delta^2)\).
  This quantifies the operator-norm route. For the unresolved
  short-window estimate, a direct bound on the moving quadratic
  form could avoid the full norm cost.
  [Sharp cost and moving-window version](arithmetic/EXTERNAL_HECKE_ARTIN_GOLDBACH_TRANSFER_AUDIT.md#sharp-hardy-cost-over-all-supported-test-shapes)

## Stone's criterion and the location of a singleton fibre

In the [pinned Stone note](https://github.com/1337isnot1337/solved-conjectures/blob/dab1ac8d1ecea8dc34c38932fcd2a4b3f56ce0ef/research/stone-nonproper-locus/README.md),
the complementarity map on each orthant is linear in nonnegative
coordinates with matrix \(B_I\). For a sequence with images converging
and norms escaping, take a fixed orthant along a subsequence and
normalize its nonnegative coordinates by their norms. The limiting
nonzero vector belongs to the nonnegative kernel of \(B_I\).
The limit of the images lies in the closed image cone of that
orthant. In the reverse direction one can add an arbitrary
nonnegative multiple of the kernel vector to any preimage in the
cone. This verifies the identification of the nonproper locus with
the union of the stated polyhedral cones and, in this setting,
with the locus of unbounded fibres.

The other necessary fact in the source is uniqueness at an interior
point of a nonsingular complementary cone. Pivoting on that basis
takes the strictly complementary solution to the solution with
zero complementary input and positive right-hand side. Another
solution would have a nonzero nonnegative input with strictly
negative matrix action on its support, forbidden by semimonotonicity
of the pivoted matrix. Q0 identifies the image with the convex
feasible cone. Along a compact path in its interior avoiding the
nonproper locus all preimages are uniformly bounded, so Brouwer
degree on one large ball is constant. At a unique preimage the
degree is the sign of the corresponding principal minor. At a
strictly positive right-hand side it is \(+1\). If the exceptional
cones intersect the interior only in codimension at least two,
the two-segment general-position argument supplies the path. These
are the two separate hypotheses that yield the degree obstruction.

For complex polynomial Keller maps, the complement of the **nonfinite**
locus carries the finite étale cover of generic degree \(d\). A
single fibre of size one over that complement forces \(d=1\).
The source being open in its finite normalization is essential
here. If the nonfinite locus has **complex** codimension at least two,
the full finite normalization has no height-one ramification. Purity
applies because the normalization is normal, the affine target is
regular, the finite map is quasi-finite, and local dimensions agree by
going-down. The result is a connected finite étale cover of complex
affine space and therefore has degree one. A Keller morphism of
degree one is a birational open immersion by the quasi-finite
form of Zariski's main theorem. Its image meets every target prime
divisor generically: for an irreducible polynomial \(q\) cutting out
the divisor, \(q(F)\) is nonconstant by dominance, and a prime
component of its zero set dominates the divisor by dimension and
quasi-finiteness. The open immersion image thus omits at most a
codimension-two closed set. The coordinate functions of its inverse
extend across that set by normality of the affine target. Both
polynomial composition identities follow on dense opens and hence
everywhere. [Affine inertia and cover proof](jacobian/monodromy-theory/THEORY.md#exact-purity-and-affine-cover-applications)

Real codimension two of a complex exceptional hypersurface can
provide path connectivity and constancy of \(d\), with no conclusion
about its value. Fibre cardinality at an exceptional point need
not give \(d\). For the retained three-dimensional Alpöge–Gallagher
map, put \(U=1+xy\), and write its three outputs as
\[
 A=U^3z+y^2U(4+3xy),\quad
 B=y+3xU^2z+3xy^2(4+3xy),\quad
 C=2x-3x^2y-x^3z .
\]
Its Jacobian determinant is the constant \(-2\).
At a point with \(A=B=0\) and \(U\ne0\), the identity
\(B=y+3xA/U\) gives \(y=0\), then \(z=0\). On \(U=0\) one has
\(B=-2y\ne0\). Thus the fibre of \((0,0,c)\) is the singleton
\((c/2,0,0)\) for every \(c\in\mathbb C\). The arc
\((x,y,z)=(t^{-1},0,2t^2-ct^3)\) escapes with limiting image
\((0,0,c)\). In contrast \((0,1,0)\) has the three distinct preimages
\((0,1,-4),(-2,1,2),(2,-1/2,5/4)\). Their local analytic inverse
branches show \(d\ge3\). These identities and the determinant
were checked in exact symbolic arithmetic. The example is a
control on the inference in arbitrary dimension. The plane
Jacobian obligation and the classification of its actual omitted
divisors remain open. [Retained polynomial control](jacobian/affine-complements/NORMAL_FORM_COVERAGE_AND_COMMON_LEADING_DIVISORS.md#a-mandatory-three-dimensional-control-for-integrality)
