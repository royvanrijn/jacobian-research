# Fixing the resolvent polynomial tests solubility, not class creation

There is an exact false-negative mechanism for a tempting class-creation
constructor. Searching for unramified S4 extensions with a given **cubic
field** is the incidence problem. Requiring their standard trace-free
resolvent to equal a specified **cubic polynomial** imposes a rational
elliptic 2-cover condition as well.

On the retained MW16-05 reference, the distinction excludes an entire
verified six-dimensional strict block:

| Constructor condition on the same cubic field | Classes admitted from the marked block |
|---|---:|
| Arbitrary cubic generator; strict unramified S4 class | all 63 nonzero classes |
| Trace-free resolvent exactly f(Z)=Z^3+AZ+B | none |
| Trace-free resolvent exactly -f(-Z)=Z^3+AZ-B | all 63 nonzero classes |

The exclusions hold for **every rational quartic coefficient choice**, not
just a bounded coefficient box. The second and third cubics define the
same field, through theta -> -theta. The six-dimensional class block and
its strict local conditions are unchanged.

We construct explicit quartics for six independent classes using the
first row, with quadratic expressions in theta as their resolvent
generators. Independent algebra and prime witnesses verify these
quartics. Their strictness and the obstruction use the already certified
**generic** control; no exceptional point, class or rank label enters.
These are existing classes, not newly created specialization directions.

The quartic number fields here encode unramified double covers of the
specialized cubic field K_t. They are different objects from the preceding
[degree-four scheme inside the fixed root Jacobian](A_QUARTIC_GOVERNS_THE_LAST_GLOBAL_CLASS.md).

## General quartics from norm-square classes

Let K=Q(theta), with theta^3+A theta+B=0 irreducible and nonsquare
discriminant. Take a nontrivial norm-square class represented by

    alpha=a+b*theta+c*theta^2,     N(alpha)=n^2 != 0.

Let T=Tr(alpha) and S be the second elementary symmetric function of its
three conjugates. The polynomial

\[
 Q_\alpha(Z)=Z^4-\frac T2Z^2-nZ+\frac{T^2-4S}{16}
\tag{1}
\]

defines its S4 quartic field. To see the construction directly, choose
square roots u_i of the conjugates alpha_i with u_1*u_2*u_3=n. Its four
roots are

    (u1+u2+u3)/2,  (u1-u2-u3)/2,
    (-u1+u2-u3)/2, (-u1-u2+u3)/2.

Their product polynomial is (1). A nontrivial norm-square class over a
non-Galois cubic gives the S4 extension; the norm-square correspondence
and relative-discriminant identity are also given by
[Cohen–Thorne, Theorem 2.2 and Proposition 2.3](https://link.springer.com/article/10.1007/s40993-015-0001-y).

For a depressed quartic Q(Z)=Z^4+pZ^2+qZ+w, use the standard resolvent

\[
 R_Q(U)=U^3-pU^2-4wU+4pw-q^2.
\tag{2}
\]

Its roots are the three sums of products associated with partitions of
four roots into pairs. For (1), these are alpha_i-T/2. Therefore

    R_Q(U)=minpoly_alpha(U+T/2),
    trace-free resolvent generator = alpha - T/3.

For alpha=a+b theta+c theta^2 this generator is

\[
 b\theta+c\left(\theta^2+\frac{2A}{3}\right).
\tag{3}
\]

A general class naturally produces a quadratic expression in theta.
Allowing that expression is essential to the field-incidence question.

The exact change-of-basis determinant from (1,theta,theta^2) to
(1,alpha,alpha^2) is

\[
 \delta=b^3+Abc^2+Bc^3=N(b-c\theta).
\]

Both symbolic and numeric checks verify

\[
 \operatorname{disc}(Q_\alpha)
 =\operatorname{disc}(\operatorname{minpoly}_\alpha)
 =\operatorname{disc}(f)\,\delta^2.
\tag{4}
\]

These are **polynomial** discriminants. For a strict class, the **field**
discriminant of the quartic still equals Disc(K), by the unramified
relative-discriminant identity. A large factor delta^2 in (4) need not
be new field ramification. With integral power bases the corresponding
changes are order indices; for arbitrary rational generators, delta is
simply the rational basis determinant. No fresh maximal-order calculation
is claimed here.

## Fixing a generator chooses an elliptic twist

Fix a nonzero rational scalar lambda, and require that the trace-free
resolvent generator under the identification with K be exactly
lambda*theta. Its polynomial is then

    U^3 + lambda^2*A*U + lambda^3*B.

For a quartic representing the class [alpha], equation (3) means that
some representative in its squareclass must have the form

\[
 \alpha z^2=s+\lambda\theta,\qquad z\in K^*,\ s\in\mathbb Q.
\tag{5}
\]

Taking norms gives

\[
 v^2=s^3+\lambda^2As-\lambda^3B.
\tag{6}
\]

This is exactly the elliptic curve E^(-lambda), for
E:y^2=x^3+Ax+B. Its two-division root is -lambda*theta, and its Kummer
class at the point (s,v) is s+lambda*theta. Conversely a rational point
of E^(-lambda) in class [alpha] supplies z in (5). Thus, for a nontrivial
class,

\[
 \boxed{\text{quartic presentation with prescribed oriented resolvent}
 \iff [\alpha]\in\delta E^{(-\lambda)}(\mathbb Q).}
\tag{7}
\]

The associated quartic is explicitly

\[
 Z^4-\frac{3s}{2}Z^2-vZ-\frac{3s^2}{16}-\frac{\lambda^2A}{4}.
\tag{8}
\]

Substitution into (2) leaves precisely (6) as its coefficient condition.
Translation of a general monic quartic removes its cubic coefficient
without changing the centered pair-product resolvent, so the statement
is not limited by the depressed-quartic normalization.

The field identification matters. A non-Galois cubic has no nontrivial
Q-automorphism, so prescribing its actual trace-free polynomial and the
identification cannot silently permute the original theta into another
rational generator. Allowing arbitrary isomorphic cubic polynomials
removes the restriction; that is the first row of the opening table.

Equivalently, using the
[general 2-cover equations](ONE_COMMON_CLASS_FOR_A_LARGE_NATIVE_BLOCK.md),
write alpha*z^2=Q0(z)+Q1(z)theta+Q2(z)theta^2. The fixed-generator problem
is the genus-one cover

    Q2(z)=0,       Q1(z)=lambda*w^2.

A solution of its common conic alone does not provide the required
square value. Rational solubility of this particular cover, not just
existence of the class field, is what (7) tests.

## Exact six-class experiment

The [frozen protocol](FIXED_RESOLVENT_SOLUBILITY_PROTOCOL.json) projects
only the cubic equation, its certified field discriminant, the six
strict generic representatives and the unconditional CT matrix from the
[fixed-incidence switch](FIXED_INCIDENCE_SIX_DIRECTION_SOLUBILITY_SWITCH.md).
The worker forms (1) for each representative and verifies its resolvent,
change-of-generator determinant and discriminant identity. All six
representatives have a nonzero theta^2 coefficient in (3).

| Generic mask | Prime with irreducible quartic reduction | Prime with irreducible original cubic f reduction |
|---:|---:|---:|
| 343 | 11 | 23 |
| 2101 | 7 | 23 |
| 4212 | 11 | 23 |
| 8942 | 7 | 23 |
| 16624 | 7 | 23 |
| 34044 | 11 | 23 |

Every indicated reduction is separable. Quartic irreducibility supplies
a 4-cycle. The invertible generator change identifies the resolvent field
with K, whose original cubic f is irreducible modulo 23; this supplies
the transitive resolvent action. Together these certify S4. The prime search was bounded
by 503 and stops on these exact witnesses. It does not search parameters
or rational points.

Each quartic has field discriminant

```text
128900477062442043600727490102612931938219670661531295245188752203875468
```

This equality is **deduced from the retained strict-class certificates**
and the standard correspondence. The independent replay checks the new
polynomial arithmetic, not a new maximal-order discriminant calculation.
The existing strictness certificates also imply local quartic algebra
Q_p × (K tensor Q_p) at every strict bad place, and complete real splitting.

The twist CT matrix on these classes has rank six. For every nonzero
six-bit word, the new artifact supplies an explicit basis class with
CT pairing one. Hence all 63 words are nonrational on E^(-1). Equation
(7), with lambda=1, excludes the fixed f resolvent presentation for every
one of the corresponding unramified S4 fields.

On the original E all these words are rational generic Kummer classes.
Equation (7), with lambda=-1, therefore supplies the reflected-resolvent
presentations. The proof does not search for these points or recompute
their coordinates. The norm-square/S4 correspondence is injective here
because K is non-Galois; the 63 nonzero classes represent distinct quartic
field isomorphism classes, not multiple generators of one field.

All statements about the six-dimensional block and its exclusions are
unconditional. No full class-group dimension, GRH assumption, exceptional
point or completed rank label is used in this experiment.

## Consequences for the class-creation question

This rules out using a fixed-oriented-resolvent count as an independent
incidence feature. On this control it can return zero for an entire
six-dimensional block which is already known to exist. Increasing its
coefficient bound cannot fix the problem.

The useful distinctions for future constructors are now precise:

- **Incidence:** independent norm-square classes defining strict
  unramified quadratic extensions of K_t, or their S4 quartic fields,
  allowing arbitrary cubic generators and testing maximal-order
  ramification. This is the still-missing object for the largest jumps.
- **Solubility:** realizing those fields by quartic polynomials with the
  prescribed oriented trace-free resolvent. The equivalence (7) is exact
  and depends on the elliptic twist, while the cubic field stays fixed.
- **Visibility:** small quartic coefficients, small generator determinant
  delta, or recovery within a bounded monogenic representation search.
  None is promoted to a class-rank predictor.

For the fresh +10/+11 and historic +12/+14 examples, the next required
incidence computation remains a construction of strict classes outside
the whole inherited pool. It must allow (3) and distinguish (4) from
field ramification. This turn establishes a concrete failure mode of a
restricted constructor and supplies six explicit control fields; it
constructs no extra specialization class and proves no condition on the
original family parameter t. The largest-jump mechanism remains open.

The present priority is therefore the full cubic S-class problem or an
ideal-square constructor with arbitrary cubic generator. The verified
CT switch can then test rational solubility. Fixed-polynomial searches
and coefficient-size scores cannot substitute for that incidence step.

## Replay

The [main certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_fixed_resolvent_solubility_v1.json)
contains the six quartics, their cubic generators and all 63 CT witnesses.
The [independent verifier](../../artifacts/generated-results/elliptic-curves/rank_jump_fixed_resolvent_solubility_verification_v1.json)
expands the four square-root sums in SymPy, reconstructs each quartic in
an independent rational cubic algebra, and rechecks twelve finite-field
irreducibility witnesses. It binds the retained unconditional strict/CT
verification without rerunning its number-field computations.

```sh
timeout 30 sage -python elliptic-curves/rank-jump/fixed_resolvent_solubility.py check
timeout 30 sage -python elliptic-curves/rank-jump/verify_fixed_resolvent_solubility.py check
```

All outputs are new rank-jump artifacts. Active search protocols, candidate
populations, worker limits, scoring policies and mathematical-status
entries are untouched.
