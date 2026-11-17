# Author-side obligation audit

This is not independent external review. The exact calculations support
the written arguments; they do not formalize those arguments.

| Obligation | Evidence and limit |
|---|---|
| Scaling covers the entire Laurent domain | z=x^3 with x nonzero; all powers after source/target scaling are multiples of three. The scaling is applied to the whole target when used in global matching. |
| Collision projection artifacts are kept separate | The v^10 factor is stripped only with the independent v0 chart retained. The degree23 squared factor is checked by a modular gcd, a quadratic common fibre, its discriminant and Jacobian norm. W has its own direct regular chart. |
| The degree27 factor gives one ordinary tacnode | Its discriminant exponent is one, all other factors are coprime, and fixed Sylvester degrees force corank one. Linear recovery is regular and no triple/nonimmersion occurs there. |
| Triple necessity loses no chart | S0 is eliminated separately with leading-coefficient and resultant checks. S nonzero uses only constant and S divisions; two literal integer identities force D50 or z=-1/3, and a separate unit excludes the latter. |
| Every D50 root has an actual triple | Four original cubic-remainder equations vanish at the exact point over Q[z]/D. Its denominator is nonzero, and immersion forces its cubic roots distinct. |
| The finite-field rank witness proves a characteristic-zero statement | Two verified modular factorization patterns prove D irreducible. A nonzero modular minor at a root of D proves its determinant is nonzero in the number field by Gauss's lemma. The exact point supplies the matching upper rank bound. This is not inference from a sampled curve. |
| Uniqueness holds for every conjugate root | The truncated coefficient matrix has a one-dimensional kernel over the degree50 field and hence over its algebraic closure. Unit and coordinate monomials fix the normalized evaluation vector, ruling out any second unordered triple or four-branch value. |
| Both cusp rows are inherited correctly | Exact scaled p,q comparisons identify class I directly and class II after t maps to -t. Their singularity types remain explicit PF2R7CUSP1 inputs. |
| Every affine singularity is counted | The family is finite birational and immersed outside the known two cusp values. Every other singularity is a collision. The pair charts, tangencies and all triples are exhausted. |
| Counts do not become an affine-moduli count | 237 counts Laurent x parameters, including the cube-root choices. The coefficient union count268 subtracts its three actual L/K intersection points. No uniqueness modulo all target transformations is asserted. |
| No global closure is imported | Triple fibres exhaust degree six and forbid further omitted-image intersections, but source surjectivity is not assumed. No proper infinity comparison, braid transport, second image or Keller realization is supplied. |
