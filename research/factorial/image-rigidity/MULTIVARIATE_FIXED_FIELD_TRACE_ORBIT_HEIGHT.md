# Fixed-field height bounds for reciprocal trace orbits

A common factor whose root field remains in a finite set has a
stronger projective denominator requirement near an algebraic
diffuse contact. The estimate works in every factor degree and
permits the primes in that denominator to change. It applies
separately to each unordered Galois orbit of root pairs.
The same field data sharpen the local upper bound by restricting
which root valuations can be assigned to the embeddings.

Fix consecutive ordinary moment polynomials and an irreducible
monic rational common factor \(H\) with nonzero roots
\(\alpha_1,\ldots,\alpha_d\). Use the common rational argument
normalization at a fixed positive diffuse reference. It
does not change the fields or root quotients. Let \(\eta\to0\)
along a sequence, and take a matched algebraic reference factor
selection as in the
[projective factor comparison](MULTIVARIATE_DIFFUSE_PROJECTIVE_FACTOR_SIEVE.md#consequences-at-finite-physical-diffuse-contacts).
The roots and reference roots lie in a fixed compact annulus.
Write
\[
         |\alpha_i-\alpha_{i,*}|=O(\eta^{\gamma_i}),\qquad
                              \gamma_i>0 .
\]
The result also applies if this matching and its rates are
supplied directly by other means.

## Denominators by Galois orbit

Let \(L=\mathbb Q(\alpha_1)\) and let \(E\) be its normal closure.
Identify the \(d\) roots with the images of \(\alpha_1\) under
the embeddings of \(L\) into \(E\). Thus the action on the
index set depends on the field \(L\) and its fixed embedding
labels. For an orbit \(O\) on unordered pairs of distinct
indices, put
\[
 t_{ij}=\frac{\alpha_i}{\alpha_j}+\frac{\alpha_j}{\alpha_i},
 \qquad T_O(Y)=\prod_{\{i,j\}\in O}(Y-t_{ij}).
\]
This is a rational monic polynomial. Let \(D_O\) be its least
positive denominator. The full trace denominator \(D_H\) of the
[moving-prime trace resultant](MULTIVARIATE_MOVING_PRIME_TRACE_RESULTANT.md)
satisfies \(D_H=\prod_O D_O\). Indeed, each \(D_O T_O\) is
primitive, so Gauss's lemma applies to their product, including
when some trace values repeat.

Fix an embedding into \(\overline{\mathbb Q}_p\) and write
\(s_i=v_p(\alpha_i)\). If \(s_i\ne s_j\),
\(v_p(t_{ij})=-|s_i-s_j|\). If \(s_i=s_j\) the trace is
integral at \(p\). The minimum valuation among the coefficients
of a rational monic polynomial is the sum of the negative
valuations of its roots. To see equality in the usual
elementary symmetric bound, take the coefficient with exactly
as many root factors as there are negative root valuations.
The product on that entire set is the unique minimizer.
Consequently
\[
                  v_p(D_O)=\sum_{\{i,j\}\in O}|s_i-s_j|.           \tag{1}
\]

Choose one ordered Galois orbit \(U\) above \(O\), and one
\((i,j)\in U\). Put \(K=\mathbb Q(\alpha_i,\alpha_j)\) and
\(\beta=\alpha_i/\alpha_j\). The embeddings of \(K\) correspond
bijectively to \(U\), even if some of the images of \(\beta\)
coincide. Set \(e=|O|\), and let \(c=1\) or \(2\) be the number
of orientations in \(U\) above one unordered pair. Thus
\([K:\mathbb Q]=ce\). The monic quotient polynomial on \(U\)
has least denominator
\[
 \log D_U=\sum_p\sum_{(i,j)\in U}\max(0,s_j-s_i)\log p.
\]
If \(c=2\) this is exactly \(\log D_O\).
If \(c=1\), the reversed orbit accounts for the other
orientations. Transitivity on the root labels gives equal
in-degree and out-degree at every label in \(U\).
In particular \(\prod_{(i,j)\in U}(\alpha_i/\alpha_j)=1\),
so the sum of the positive valuations matches the sum of
the negative magnitudes at each prime. Hence
\(D_U^2=D_O\) when \(c=1\).

Let \(h\) denote the absolute logarithmic Weil height.
The finite-place part of \([K:\mathbb Q]h(\beta)\) is exactly
\(\log D_U\). At a finite place above \(p\) the summand with
absolute value extending \(|p|_p=p^{-1}\) occurs once for each
of its \([K_v:\mathbb Q_p]\) embeddings into
\(\overline{\mathbb Q}_p\). This is precisely the multiplicity
in the ordered orbit sum. The archimedean part is bounded
uniformly by the compact annulus for the roots. We obtain
\[
                       \log D_O=2e\,h(\beta)+O(1).                \tag{2}
\]
The passage to the quotient parameter in (2) doubles the
height coefficient without increasing the rate at which
the root matching controls the approximation error.
Repeated quotient or trace values do not affect (1) or (2).

## Fixed-field approximation

Suppose the root fields \(L\) belong to a fixed finite set of
number fields up to isomorphism. After passage to a subsequence
we can fix an embedded \(L\), the normal closure \(E\) in a
fixed algebraic closure, the orbit \(O\), and the labels
in the matching. There are only finitely many such choices:
each field has finitely many embeddings in the algebraic
closure and the reference selection and root labels are finite.
Then \(K\) is fixed, since it is the compositum of the two
fixed images of \(L\) in \(E\). In particular the hypothesis
holds when the degrees and absolute field discriminants of
\(L\) are bounded, by the
[Hermite finiteness theorem, Theorem 8.43](https://www.jmilne.org/math/CourseNotes/ANT.pdf).
A fixed finite set of possible ramified primes and bounded
degree also suffice, by Theorem 8.42 of that source.

For the lower bounds alone, it suffices that the quotient
fields \(\mathbb Q(\beta)\) belong to a fixed finite menu.
After fixing that field, each of its embeddings occurs
\([K:\mathbb Q(\beta)]\) times in \(U\). The repeated
quotients have the same limit and the best rate among
their matched pairs. Thus one can sum the individual
rates with these multiplicities in the argument below,
even when the larger pair fields vary.

Put
\[
                   \Gamma_O=\sum_{\{i,j\}\in O}
                                         \min(\gamma_i,\gamma_j).
\]
For every \(\epsilon>0\), for sufficiently small \(\eta\), either
all ordered quotient limits on \(U\) are exactly attained or
\[
       \log D_O\ \ge\ \frac{2\Gamma_O}{2+\epsilon}
                                 \log(1/\eta)-O(1).               \tag{3}
\]
If the limiting ordered quotient tuple consists of the
embeddings of one element \(\beta_*\in K\), the alternative
when it is not attained strengthens to
\[
                     \log D_O\ge 2\Gamma_O\log(1/\eta)-O(1).      \tag{4}
\]
Constants may depend on the fixed fields, reference, matching
and \(\epsilon\). The same inequalities apply to all allowed
finite-field branches, by taking the largest constants and
smallest neighborhoods across their finite list.

To prove (3), use precisely the square-form
[number-field Subspace Theorem, equation (1.1)](https://arxiv.org/pdf/math/0408381)
on \(\mathbb P^1(K)\), with the absolute-value normalization
in its section 1.4. At each archimedean place \(v\) choose an
embedding \(\sigma\) and use the independent forms
\(X_0\) and \(X_1-\theta_v X_0\), where \(\theta_v\) is the
matched algebraic limit of \(\sigma(\beta)\).
Since proper linear subspaces of \(\mathbb P^1\) are points,
outside a finite set of \(\beta\) the conclusion is
\[
 \sum_{v\mid\infty}[K_v:\mathbb R]
           \bigl(-\log|\sigma(\beta)-\theta_v|\bigr)
                 \le (2+\epsilon)[K:\mathbb Q]h(\beta).
\]
In the theorem the norm of \((1,\beta)\) in the denominator
is at least one at each place, so discarding it gives this
inequality. The complex local degree is two and thus counts
both conjugate embeddings. Their matched errors have equal
magnitudes. At the embedding that sends the base pair to
\((i',j')\), the numerator error is
\(O(\eta^{\min(\gamma_{i'},\gamma_{j'})})\). Summing with
multiplicity gives \(c\Gamma_O\log(1/\eta)-O(1)\).
Both \(c\) and \([K:\mathbb Q]=ce\) cancel in (2).
A member of the finite exceptional set cannot occur
arbitrarily close to the limiting tuple unless it is
exactly that tuple. This proves the alternative in (3).

For (4), put \(w=\beta-\beta_*\ne0\). Its archimedean
sum of negative logarithms is at least
\(c\Gamma_O\log(1/\eta)-O(1)\). By the product formula,
the finite-place positive part of \([K:\mathbb Q]h(w)\)
is at least that large. The elementary height inequality
\(h(w)\le h(\beta)+h(\beta_*)+\log 2\) and (2) give (4).
When several unordered orbits are nonattained, their
lower bounds add in \(\log D_H\).

One can also include any fixed finite set \(V\) of finite
places of \(K\) at which \(\beta\) is a unit. Suppose for
each \(v\in V\) above \(p\) there is a fixed finite list of
unit algebraic targets and one of them satisfies
\(v_p(\beta-\theta_v)\ge k_v\). The numerator of the
parenthesis on the right of (3), before multiplication
by \(2/(2+\epsilon)\), then acquires the term
\[
                    \frac{1}{c}\sum_{v\in V}
                      [K_v:\mathbb Q_p]\,k_v\log p.              \tag{6}
\]
Use the same two forms at these places in the
Subspace Theorem. The unit hypothesis makes the vector norm
one, and the normalized absolute value supplies the local
degree. A finite union over the target lists preserves the
exceptional-point alternative. All other primes remain
free to vary in the denominator in (2). In the coherent
case where these targets are the local images of the same
\(\beta_*\), the product formula likewise adds (6) to the
\(\Gamma_O\log(1/\eta)\) term in (4) before multiplication
by two.

## Upper bounds with moving primes

At a prime \(p\) choose a prime of \(E\) above it and let \(G_p\)
be its decomposition group acting on the \(d\) labels.
The assignment of values \(s_i\) is constant on the orbits
of \(G_p\) by invariance of the chosen valuation. It must
also respect the shared Newton edge widths and any
[truncated Newton width bounds](MULTIVARIATE_PHYSICAL_TRUNCATED_NEWTON_FACES.md#common-root-widths-at-higher-precision)
from the parent polynomials. Define \(B_{p,O}\) as the maximum
of the right side of (1) over these constrained assignments.
An empty assignment set already excludes that field.
Otherwise the candidate bound is
\[
                     \log D_O\le\sum_p B_{p,O}\log p.             \tag{5}
\]
Conjugating the chosen decomposition group relabels the
assignments and does not change the maximum. The sum is
finite at each profile. A prime outside a primitive parent's
endpoint support gives unit roots in its original variable.
A rational rescaling translates all those valuations by the
same amount, preserving their zero dispersion.
At an unramified prime in \(E\), the orbits of \(G_p\) on
the labels are the cycles of its Frobenius permutation.
One can also use the local embedding action at ramified
primes. In particular if \(G_p\) carries one end of a pair
to its other end, that pair contributes zero to (1),
even when other pairs or orbits have nonzero dispersion.

For a nonattained orbit, (3) and (5) are incompatible if the
right side of (5) stays below
\((1-\delta)\Gamma_O\log(1/\eta)\) for any fixed \(\delta>0\).
In the coherent case of (4), the analogous coefficient
is \(2\Gamma_O\). These are orbit tests. A budget at a
different orbit cannot compensate for a shortage here.
The field menu and the nonattainment check are required
for this exclusion. Bounds on the ramification support
of the candidate fields can come from the
[common-factor discriminant ledger](MULTIVARIATE_COMMON_FACTOR_FIELD_DISCRIMINANT.md).

For a quadratic \(H\) there is one unordered pair.
On a nonattained branch with fixed quadratic field,
(3) gives \(D_H\gg\eta^{-2\gamma/(2+\epsilon)}\) when both
matched rates are \(\gamma\). This doubles the exponent
from Roth applied directly to the rational trace at rate
\(\gamma\). More generally this comparison applies to the
same pair-error rates in higher degree. If a trace has
extra vanishing, for example at a coincident limiting
quotient \(1\), it can have a faster rate than the quotient.
Each bound uses its own certified rates.



## Exact orbit denominator fixtures

For the cyclic quartic generated by \(\zeta_5\), take
\(\alpha=\zeta_5+2\), ordered using the exponent cycle
\(1,2,4,3\). Adjacent pairs form an orbit whose orientations
split. Its denominator is \(121\), and either ordered
polynomial has denominator \(11\). Opposite pairs form a
swap orbit with denominator \(11\). The full denominator
is \(1331\). All nonzero valuations here come from a split
prime, as in (5). The same arithmetic checks apply to
shifts \(\zeta_5+a\) with arbitrary integer \(a\).

The following exact replay tests both orientation patterns
over 25 shifts. It also tests the single unordered pair orbit
in the \(S_4\) quartic \(x^4-x+1\) over 16 shifts. At each
tested prime other than \(229\) dividing the constant term
of that shifted quartic, exactly one of the four roots has
positive valuation, and its valuation is the valuation of
the constant term. The primitive trace denominator therefore
has three times that valuation. The fixtures check the
denominator calculations for algebraic polynomials.

```python
from math import comb
from sympy import Poly, QQ, Symbol, factorint, ilcm, oo, resultant
from sympy.polys.numberfields import galois_group

x = Symbol("x")
Y = Symbol("Y")
Z = Symbol("Z")

def denominator(poly):
    return ilcm(*[int(c.q) for c in poly.all_coeffs()])

def polynomial_orbit_product(values, modulus):
    result = Poly(1, x, domain=QQ[Y])
    for value in values:
        result = (result * Poly(Y - value.as_expr(), x, domain=QQ[Y])).rem(
            Poly(modulus.as_expr(), x, domain=QQ[Y]))
    assert result.degree() == 0
    return Poly(result.nth(0), Y, domain=QQ)

def ordered_quotient(H):
    scaled = sum(c*Z**(H.degree() - i)*x**i for (i,), c in H.terms())
    out = Poly(resultant(H.as_expr(), scaled, x), Z, domain=QQ).monic()
    return out.exquo(Poly((Z-1)**H.degree(), Z))

def trace_from_ordered(Q):
    N = Q.degree() // 2
    remaining = Q
    answer = 0
    for k in range(N, -1, -1):
        leading = remaining.nth(N+k)
        answer += leading * Y**k
        basis = Poly(sum(comb(k,j)*Z**(N+k-2*j)
                   for j in range(k+1)), Z)
        remaining -= leading * basis
    assert remaining.is_zero
    return Poly(answer,Y,domain=QQ)

modulus = Poly(x**4+x**3+x*x+x+1, x, domain=QQ)
adj = 0
opp = 0
for a in range(25):
    alpha = [(Poly(x**k+a,x,domain=QQ)).rem(modulus) for k in [1,2,4,3]]
    reciprocal = [value.invert(modulus) for value in alpha]
    edge = [(alpha[k]*reciprocal[(k+1)%4]).rem(modulus) for k in range(4)]
    redge = [(alpha[(k+1)%4]*reciprocal[k]).rem(modulus) for k in range(4)]
    edge_trace = [u+v for u,v in zip(edge,redge)]
    diagonal = [(alpha[k]*reciprocal[(k+2)%4]).rem(modulus) for k in range(4)]
    diag_trace = [diagonal[k]+diagonal[k+2] for k in range(2)]

    P_edge = polynomial_orbit_product(edge,modulus)
    P_reversed = polynomial_orbit_product(redge,modulus)
    T_edge = polynomial_orbit_product(edge_trace,modulus)
    P_diag = polynomial_orbit_product(diagonal,modulus)
    T_diag = polynomial_orbit_product(diag_trace,modulus)
    assert denominator(T_edge) == denominator(P_edge) ** 2
    assert denominator(P_edge) == denominator(P_reversed)
    assert denominator(T_diag) == denominator(P_diag)
    H_a = Poly(modulus.as_expr().subs(x,x-a), x)
    T_H = trace_from_ordered(ordered_quotient(H_a))
    assert T_H == T_edge*T_diag
    assert denominator(T_H) == denominator(T_edge)*denominator(T_diag)
    assert all(p%5==1 for p in factorint(denominator(T_H)))
    if a == 2:
        assert (denominator(P_edge),denominator(T_edge),denominator(T_diag),denominator(T_H)) == (11,121,11,1331)
    adj += 1
    opp += 1

non_cm = Poly(x**4-x+1,x)
G, _ = galois_group(non_cm)
assert G.order() == 24
assert non_cm.count_roots(-oo, oo) == 0
assert non_cm.discriminant() == 229
s4 = 0
for a in range(16):
    H_a = Poly(non_cm.as_expr().subs(x,x-a),x)
    Q_H = ordered_quotient(H_a)
    T_H = trace_from_ordered(Q_H)
    assert T_H.degree() == 6 and denominator(Q_H) == denominator(T_H)
    for p,e in factorint(abs(H_a.nth(0))).items():
        if p != 229:
            exponent = factorint(denominator(T_H)).get(p,0)
            assert exponent == 3*e
    s4 += 1
print(f"{adj} split-orientation orbit identities, {opp} swap-orbit identities, {s4} S4 orbit and good-prime dispersion identities")
```
