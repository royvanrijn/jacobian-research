# Elementary exponent decision with full jets

The fixed-coordinate odd-height fibre algorithm eventually asks whether
\([R]^j=b\) in a finite number-field algebra, with \(j\ge0\).
The retained proof imports a general point-orbit decision theorem for
this step. Here we supply an independent terminating proof and an exact
implementation, retaining every nilpotent direction. A repeated factor
of the cyclic minimal polynomial already determines the only possible
exponent from a derivative.

**This is a rederivation of a classical decidability result, not a new
Orbit theorem.** Kannan–Lipton proved the stronger polynomial-time
point-to-point result. Chonev–Ouaknine–Worrell, Theorem 2.2, contains the
stronger one-dimensional-target complexity statement used in the
retained route. The present algorithm makes no polynomial-time claim.
It does not decide the vanishing of a scalar observation of an orbit,
general recurrence zero sets, or higher-dimensional target spaces.

## 1. The complete exponent algorithm

**Theorem ElementaryRationalPointOrbit.** Given a rational square matrix
\(M\) and rational column vectors \(v,w\) of the same positive dimension,
there is a terminating procedure returning the exact set
\[
                       \{j\ge0:M^jv=w\}.                 \tag{O.1}
\]
Its output is empty, a singleton, or one progression
\(j_0+d\mathbb Z_{\ge0}\), with explicit integers \(j_0\ge0,d\ge1\).
The procedure uses rational linear algebra, rational polynomial
factorization, integer arithmetic and finite searches with the proved
bounds below. It uses no Orbit oracle, logarithmic approximation,
complex root-isolation routine, or unproved exponent cutoff.

This includes multiplication by an element of any finite-dimensional
number-field algebra after restriction of scalars. The input need not
be invertible, diagonalizable, or semisimple. In fact the point-orbit
formulation does not require an ambient algebra at all.

The [implementation](verify_elementary_reciprocal_exponents.py) has functions
`orbit_solutions` and `quotient_power_solutions` implement the stated
procedure. The other routines provide exact arithmetic and controls.
This is an implementation of the exponent subproblem, not an end-to-end
implementation of the entire fibre classifier.

## 2. The cyclic quotient and the zero-primary part

If \(v=0\), every exponent works when \(w=0\), and none works otherwise.
Assume \(v\ne0\). Find the first rational linear dependence in
\[
                  v,Mv,M^2v,\ldots.
\]
If it occurs at index \(r\), the earlier vectors are independent and
span an invariant space \(V\). Write
\[
 M^rv=\sum_{i<r}c_iM^iv,\qquad
                  f(T)=T^r-\sum_{i<r}c_iT^i.             \tag{O.2}
\]
The map \(q\mapsto q(M)v\) gives an isomorphism of rational
\(\mathbb Q[T]\)-modules from \(\mathbb Q[T]/f\) to \(V\).
If \(w\notin V\), the answer is empty. Otherwise solve for the unique
polynomial \(q\) of degree less than \(r\) representing \(w\).
The exact problem becomes
\[
                         T^j\equiv q(T)\pmod f.          \tag{O.3}
\]
This reduction concerns the entire target vector; it discards no jet
or field-coordinate condition.

Write \(f=T^e g\), with \(g(0)\ne0\); \(g=1\) is allowed. Test
\(0\le j<e\) directly modulo \(f\). If one matches, it is the only
match: the residues \(1,T,\ldots,T^{e-1}\) modulo \(T^e\) are distinct
and nonzero, whereas every later power is zero there. If none matches
and \(q\not\equiv0\pmod{T^e}\), there are no later matches either.

When \(q\equiv0\pmod{T^e}\), the later condition is exactly
\[
                 T^j\equiv q\pmod g,\qquad j\ge e.       \tag{O.4}
\]
If \(g=1\), every \(j\ge e\) works. Otherwise \(T\) is a unit modulo
\(g\). If \(\gcd(q,g)\ne1\), no exponent works. It remains to solve
the unit problem for \(j\ge0\) and intersect its answer with
\([e,\infty)\). These checks cover zero targets and nilpotent inputs.

## 3. A full jet determines the exponent

**Lemma RepeatedFactorExponent.** Suppose \(g(0)\ne0\),
\(\gcd(q,g)=1\), and an irreducible \(h\in\mathbb Q[T]\) satisfies
\(h^2\mid g\). Every solution of (O.4) satisfies
\[
            j\equiv Tq'(T)q(T)^{-1}\pmod h.              \tag{O.5}
\]
Thus the right side must be a constant rational polynomial representing
a nonnegative integer. If it does, that integer is the only candidate;
checking it modulo the **whole** polynomial \(g\) decides the problem.

Indeed, differentiation of \(T^j-q\), which is divisible by \(h^2\),
gives \(jT^{j-1}-q'\equiv0\pmod h\). Multiply by \(T/q\), using
\(T^j\equiv q\pmod h\), to get (O.5). This is valid also for \(j=0\).
All divisions are by units in the residue field.

Differentiate the full target representative before reducing modulo
\(h\). Replacing it by another representative modulo \(g\) does not
change (O.5), because \(h^2\mid g\) implies \(h\mid g'\). Replacing
it by its value modulo \(h\) **before differentiation** can destroy
the derivative information. Higher powers of \(h\), and all other
irreducible factors, are retained by the final test modulo \(g\).

A nonreduced ambient algebra can nevertheless have a squarefree cyclic
minimal polynomial for its particular multiplier. For example, a
multiplier congruent to 1 on \(\mathbb Q[x]/(x-1)^2\) has cyclic
polynomial \(T-1\). The preceding lemma applies to repeated factors
of the cyclic polynomial (O.2), not automatically to every repeated
factor of the original jet modulus.

## 4. A finite denominator bound for a nonintegral residue

After Section 3, assume \(g\) is squarefree. Let \(h\) be one of its
monic irreducible factors and suppose its primitive integer multiple
has leading coefficient greater than 1. Write \(q_h=q\bmod h\),
\(m=\deg q_h\), and choose a positive common denominator \(D\) for
the coefficients of \(q_h\). The unit hypothesis makes \(q_h\ne0\).

**Lemma NonintegralExponentBound.** Every solution of (O.3) satisfies
\[
                    j\le m+(\deg h)\lfloor\log_2D\rfloor.
                                                               \tag{O.6}
\]
No prime or ideal has to be identified by the algorithm. For \(D=1\),
the right side is simply \(m\).

**Proof.** Let \(H=A T^d+\cdots\in\mathbb Z[T]\) be the primitive
integer multiple of \(h\), with \(A>1\), and choose a prime \(p\mid A\).
Use a valuation on an algebraic closure of \(\mathbb Q_p\), normalized
by \(v(p)=1\). Existence of this extended valuation is the basic
input from Milne, *Algebraic Number Theory* v3.08, Theorem 7.38
(printed p.123) and Corollary 7.40 (p.124). Some root \(\alpha\) of \(H\) has \(v(\alpha)<0\).
Otherwise the elementary symmetric functions of all roots would have
nonnegative valuation, making every coefficient of the monic
polynomial \(H/A\) integral at \(p\). That would force \(p\) to
divide every coefficient of \(H\), contrary to primitivity.

In \(H(\alpha)=0\), at least two nonzero terms attain the least
valuation: a uniquely least term cannot sum to zero. For some
\(i\ne k\), therefore,
\[
 v(\alpha)=\frac{v(H_i)-v(H_k)}{k-i}.
\]
The numerator is an integer and the denominator has magnitude at most
\(d\). Since the value is negative,
\(v(\alpha)\le-1/d\). This proves the required denominator bound
directly; no Newton-polygon algorithm is needed.

The congruence implies \(\alpha^j=q_h(\alpha)\). Since
\(Dq_h\in\mathbb Z[T]\) and \(v(\alpha)<0\), the ultrametric
inequality gives
\[
 jv(\alpha)=v(q_h(\alpha))\ge m v(\alpha)-v_p(D).
\]
Dividing by the negative value yields
\(j\le m+d\,v_p(D)\le m+d\lfloor\log_2D\rfloor\).
This proves (O.6). \(\square\)

The code takes \(D\) to be the least common multiple of the rational
coefficient denominators and computes the logarithmic bound as
`D.bit_length()-1`, exactly. It tests every exponent from zero through
the displayed bound against the full modulus \(g\).

This branch is essential even when every complex conjugate lies on the
unit circle. The roots of \(5T^2-6T+5\), for example, have modulus one
and are not algebraic integers; an archimedean growth-only argument
would miss this case.

## 5. An exact coefficient-growth bound for an integral nontorsion residue

Suppose instead that \(h\in\mathbb Z[T]\) is monic irreducible, has
nonzero constant term, and its roots are not roots of unity. Write
\(d=\deg h\). Let \(C_h\) be multiplication by \(T\) in
\(\mathbb Q[T]/h\); it is an integer matrix. For \(s\ge1\), compute
\[
 \chi_s(Z)=\det(ZI-C_h^s)
          =Z^d+\sum_{k=1}^d c_{s,k}Z^{d-k}.              \tag{O.7}
\]

**Lemma IntegralGrowthCertificate.** A finite search finds \(s,k\)
with \(|c_{s,k}|>\binom dk\). From it one can find rational numbers
\(L>1\), \(B\ge1\) such that every solution has \(L^j\le B\).
Consequently every solution satisfies \(0\le j<J\), where \(J\) is
the first positive integer with \(L^J>B\), computed by rational
exponentiation and integer bisection.

**Proof.** If every coefficient in (O.7) obeyed its unit-disc bound for
every \(s\), only finitely many monic integer polynomials \(\chi_s\)
could occur. Their union of roots would be finite. For a fixed root
\(\alpha\ne0\) of \(h\), every \(\alpha^s\) belongs to that finite
union. Two powers would agree, making \(\alpha\) a root of unity,
a contradiction. This is the bounded-coefficient proof of Kronecker's
lemma, and proves that the search terminates.

Choose \(L=1+2^{-b}\), increasing \(b\) until
\[
                  \binom dk L^{sk}<|c_{s,k}|.            \tag{O.8}
\]
This also terminates. At least one complex root \(\lambda\) of \(h\)
has \(|\lambda|>L\). Otherwise every root of \(\chi_s\) has modulus
at most \(L^s\), so the elementary symmetric coefficient bound
contradicts (O.8).

Let
\[
 A=1+\max_{i<d}|h_i|,\qquad
 B=\max\left(1,\sum_i |(q_h)_i|A^i\right).              \tag{O.9}
\]
Cauchy's elementary root bound gives \(|\lambda|\le A\) for every
root, and hence \(|q_h(\lambda)|\le B\). A solution satisfies
\(\lambda^j=q_h(\lambda)\) for the root with \(|\lambda|>L\), so
\(L^j\le B\). All numbers in (O.8)–(O.9) are rational. This proves
the finite bound without numerical root isolation. \(\square\)

The implementation tests the resulting finite range against \(g\),
not merely against the selected factor. A non-torsion residue ensures
that at most one exponent can work.

## 6. The torsion branch and termination

If every factor of the squarefree \(g\) is cyclotomic, let \(d\) be
the least common multiple of its cyclotomic indices. Then
\(T\bmod g\) has exactly order \(d\). Check \(0\le j<d\); if a
representative \(j_0\) works, all and only the solutions are
\(j_0+d\mathbb Z_{\ge0}\).

The cyclotomic test is finite. If an irreducible polynomial of degree
\(r\) has a root of unity of order \(h\), then \(\varphi(h)=r\)
and
\[
                              h\le2r^2.                 \tag{O.10}
\]
Indeed \(\varphi(h)^2/h\ge1/2\): each odd prime-power contribution
is at least 1, and the contribution from 2 is at least \(1/2\).
Testing the finitely many \(\Phi_h\) in this bound determines every
cyclotomic factor and its order. If a factor is not cyclotomic, either
Section 4 or Section 5 applies. This exhausts the squarefree case.

Rational polynomial factorization and the Krylov linear algebra are
finite exact operations. Each of the remaining searches has its own
proved stopping condition. Finally intersect with \(j\ge e\) from
Section 2: a singleton is kept or discarded, and a progression is
shifted to its first member at least \(e\). This proves Theorem
ElementaryRationalPointOrbit, including all zero and repeated-root
cases. Nothing here asserts a polynomial bound on runtime or integer
bit sizes.

## 7. Application to reciprocal exponents

**Corollary ReciprocalExponentDecisionWithoutOrbitImport.** The exponent
step in the effective fixed-coordinate fibre classification can be
carried out by the preceding algorithm, without its imported general
Orbit decision theorem.

Let \(K\) be a presented number field and
\(\mathcal A=K[x]/J_{a,u}\). Handle \(J=1\) as the vacuous congruence.
For each admissible pole-twist base \(\eta_0\), all its denominators
are invertible modulo \(J\). If \([\eta_0]\) is not a unit, no exponent
can work. Otherwise the exact condition is
\[
                      [R]^j=-[\eta_0]^{-1}.              \tag{O.11}
\]
Choose a rational basis of \(K\) and the resulting rational basis of
\(\mathcal A\). Multiplication by \([R]\) is a rational matrix;
apply the theorem to the complete coordinate vectors of 1 and the
target. Restriction of scalars retains both the number-field and jet
directions. The general zero-primary branch also handles a nonunit
\([R]\), although the unit target in (O.11) already rules out positive
exponents in that case.

The bounded generator search, finite reciprocal box, and explicit
pole-twist box in Sections 3–5 of
[the effective fibre note](EFFECTIVE_ODD_START_FIBRE_CLASSIFICATION.md)
are independent of the
Orbit decision import. Combining those constructions with (O.11)
gives a second complete proof route for that theorem's termination and
output assertions. We import those bounded constructions and the
retained fixed-coordinate theorem, not the whole effective theorem
as a premise of its own alternative proof.

The algorithm distinguishes multiplicity conditions from point
values. If the multiplier actually has a nontrivial nilpotent component
in its cyclic algebra, Section 3 fixes at most one possible exponent.
If the cyclic polynomial is squarefree torsion, a periodic family is
possible even when the ambient quotient \(K[x]/J\) is nonreduced.
Whether the target belongs to the cyclic span is tested before either
conclusion. This is the precise reason a torsion-value test on the
roots of \(J\) alone is insufficient.

## 8. A coefficient field uniform in all odd heights

**Corollary UniformFibreCoefficientField.** Let \(K\) be any
characteristic-zero field and \(u\in K(x)\) satisfy \(u(0)=0\),
\(u'(0)=1\). Put \(N=\deg u\). There is a finite Galois extension
\(L/K\), of degree at most \((2N)!\), containing the coefficients of
every geometric pair in the fixed-coordinate fibre of \(u\),
**simultaneously at every positive odd height**.

The retained fixed-coordinate proof, (F15)–(F16), puts the finite
zeros and poles of every normalized reciprocal unit in the fixed set
\[
 S=u^{-1}(\infty)\cup u^{-1}(-u(\infty))\cup\{\infty\}.
\]
Write \(u=P/Q\) in coprime form, with \(P,Q\in K[x]\),
\(\max(\deg P,\deg Q)=N\), and \(Q(0)\ne0\). If
\(c=u(\infty)\) is finite, the finite support is contained in the roots
of \(Q(P+cQ)\), of degree at most \(2N\). The polynomial \(P+cQ\)
is nonzero since \(u\) is nonconstant. If \(c=\infty\), both finite
fibres in the union are the poles of \(u\), so \(Q\) alone suffices.
Remove the factor \(x\) if present, because normalization excludes a
zero or pole at the marked point, and remove repeated factors. Let
\(F\in K[x]\) be the resulting squarefree polynomial. Infinity is
already a \(K\)-rational point and requires no field extension.

Take the splitting field \(L\) of \(F\), using \(L=K\) if \(F\)
is constant. If there are \(s\) distinct roots,
its Galois group acts faithfully on them, giving
\([L:K]\le s!\le(2N)!\). Every such unit is exactly
\[
                    \eta=\prod_i(1-x/\alpha_i)^{n_i},
                    \qquad n_i\in\mathbb Z,
\]
so it belongs to \(L(x)\). There is no additional scalar parameter,
because \(\eta(0)=1\). The reconstruction
\(X=u(1+\eta)/(2xu')\), \(C=xX^a/u\), \(D=-xX^{a+2}/u\) stays in
\(L(x)\), and polynomiality makes \(C,D\in L[x]\). The same field
works for every odd \(a\). No rationality of the original return
branch is assumed. For a presented number field \(K\), the extension
is effectively constructible by finite polynomial splitting.

This strengthens the coefficient-field scope; it does not classify all
endpoint coordinates at once or give a degree bound on their pairs.

## 9. Antecedents, implementation and evidence

The point-to-point decidability problem is classical. The primary
comparison retained in this repository is Chonev–Ouaknine–Worrell,
*On the Complexity of the Orbit Problem*, Theorem 2.2 and Sections 3–4.
Its introductory statement attributes the polynomial-time point result
to Kannan–Lipton, *Polynomial-time algorithm for the orbit problem*,
JACM 33(4), 808–821 (1986), DOI 10.1145/6490.6496. The original ACM
paper was not available from its DOI during this check; no unseen
statement from it is used as a proof premise here.

Richard Lipton's author account,
[The Orbit Problem](https://rjlipton.com/2009/09/02/the-orbit-problem/)
(2 September 2009, sections “The Orbit Problem” and “The Proof”),
explicitly describes the algebraic-integer/Kronecker and nonintegral
ideal-growth split. This is historical context, accurately identified
as an author exposition, not a substitute for the written proofs above.
The current construction uses a direct denominator bound, an exact
coefficient witness, and a derivative test to make its stopping
conditions explicit. No historical priority is claimed.

The exact checker implements the algorithm on arbitrary rational
matrix/vector inputs through its Python functions. Its retained
controls include a target outside the cyclic span, both zero-start
outcomes, nilpotent early hits and tails, mixed zero/cyclotomic factors,
repeated jets that force an exponent, a corrupted higher jet, a
nonintegral unit-circle example, integral-unit growth, mixed finite
orders, number-field restriction of scalars, and a rational basis
change. It compares the complete predicted exponent descriptions with
independent direct powers on the control prefixes. Those bounded
comparisons test the implementation; the written argument proves the
full algorithm's termination and completeness.

There is no end-to-end implementation of the arbitrary-number-field
fibre classifier in this packet, no general recurrence zero-set
claim, no external mathematical review, and no formal verification.


Replay the exact controls from the repository root:

```sh
python3 research/foundations/extended-geometry/verify_elementary_reciprocal_exponents.py \
  --repo-root . \
  --output research/foundations/extended-geometry/elementary_reciprocal_exponents_verification.json
```

The [source ledger](elementary_reciprocal_exponents_source_ledger.json)
separates the required basic valuation infrastructure, the retained
finite fibre constructions, and the classical Orbit antecedents used
only for context. The [receipt](elementary_reciprocal_exponents_verification.json)
records the exact replay and its finite scope.
