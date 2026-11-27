# Affine reflections from genuine higher-prime progressions

The two norm classes in the tangent-square prime-progression
construction can both be made local powers at the exponent
prime. This preserves the exponent in the count of distinct
base fields. Their reflected degree-\(\ell\) fields then have
the smallest relative conductor in the prescribed mirror
character, and the two reflected directions give a common
degree-\(\ell^2\) field with an exact discriminant relation.

## Local powers while retaining prime progression witnesses

**Proposition PA1.** Fix an odd prime \(\ell\ge5\). In the
imaginary family of
[HP5](ABC_COHERENT_HIGHER_PRIME_PROGRESSION_FIELDS.md#removing-a-tangent-square-from-the-field-parameter)
one may further impose \(h\equiv0\bmod\ell^2\) and keep
any fixed unit residue for \(p\bmod\ell^2\). The two independent norm classes
then become \(\ell\)-th powers at the prime of the
imaginary quadratic field \(K\) above \(\ell\). This prime
is ramified and its local squareclass is \(2\ell\).
The lower count of distinct \(K\) of absolute fundamental
discriminant at most \(X\) remains
\(\gg_\ell X^{1/(\ell-1)}/(\log X)^{B_\ell}\),
for a fixed \(B_\ell>0\).

Put \(A=(p-h)^\ell,\ B=4p^\ell,\ C=(p+h)^\ell\),
and use \(E=h^2R_\ell(p,h)\) as in HP5. The polynomial
\(R_\ell\) is even in \(h\) and its constant term in \(h\)
is \(\ell(\ell-2)p^{2\ell-2}\). Thus
\[
  R_\ell(p,h)\equiv\ell(\ell-2)p^{2\ell-2}\pmod{\ell^2},
\]
of valuation exactly one. The unit factor of
\(-R_\ell/\ell\) is \(2p^{2\ell-2}\bmod\ell\), so its
squareclass is that of \(2\ell\) after restoring the
factor \(\ell\).
The norm generators are
\[
 \alpha_1=p^\ell+(A-C)/4+h\sqrt{-R_\ell}/2,\qquad
 \alpha_2=p^\ell-(A-C)/4+h\sqrt{-R_\ell}/2 .
\]
The first correction has \(\ell\)-adic valuation at least
three, since \(A-C\) has linear term
\(-2\ell p^{\ell-1}h\) and all other terms contain \(h^3\).
The second has valuation at least \(5/2\) in the
valuation normalized on \(\mathbb Q_\ell\).
Consequently each \(\alpha_i/p^\ell\) differs from one
by an element of valuation greater than two.
Hensel's lemma applied to \(Z^\ell-\alpha_i/p^\ell\) at
\(1\) gives both local powers.

At \(\ell\) the three prime forms remain units, and the
new congruence is independent of the two auxiliary
prime congruences and the modulus at \(2,3\). Apply the
same prime-form and fixed-field multiplicity arguments
as HP5 in the restricted congruence slice. Both norm
generators still come from exact oriented split-prime
ideal powers and the genuine third generator identity.

## The real mirror, fundamental discriminants and local types

**Theorem PA2.** For a field \(K\) in PA1 with
\(K\not\subseteq\mathbb Q(\zeta_\ell)\), let \(D\) be its
negative fundamental discriminant and let
\[
 B=K(\zeta_\ell),\qquad M=B^+,\qquad m=\ell-1,
 \qquad e_\ell=
 \begin{cases}
 m,&\ell\equiv1\bmod4,\\
 m/2,&\ell\equiv3\bmod4 .
 \end{cases}
\]
The two marked classes give an unramified elementary
abelian extension \(W/M\) of degree \(\ell^2\), split
at every prime above \(\ell\). The \(m\)-element
quotient \(\operatorname{Gal}(M/\mathbb Q)\) acts as
all nonzero scalars on \(\operatorname{Gal}(W/M)\).
In particular there are \(\ell+1\) distinct totally real
degree-\(\ell\) fields, with Galois group
\(\mathbb F_\ell\rtimes\mathbb F_\ell^\times\) on their
normal closures, each of discriminant
\[
 D_1=\ell^{\,m-m/e_\ell}(|D|/\ell)^{m/2}.
 \tag{1}
\]
The degree-\(\ell^2\) fixed field of a complement in \(W\)
has discriminant \(D_1^{\ell+1}\).
For every one of the degree-\(\ell\) fields the types
\((e,f)\) over \(\ell\) are
\[
 \begin{array}{c|l}
 \ell\bmod8 & (e,f)\\ \hline
 1,5 &(1,1),\ (m,1)\\
 3 &(1,1),\ (m/2,1),\ (m/2,1)\\
 7 &(1,1),\ (m/2,2).
 \end{array}
 \tag{2}
\]
The two auxiliary split-prime rows from HP5 give
independent prime classes modulo \(\ell\)-th powers
in \(\operatorname{Cl}(M)\).

The underlying reflection character and its minimum
discriminant are those of
[Cohen–Rubinstein-Salzedo–Thorne, Lemma 2.11 and Propositions 3.6–3.7](sources/cohen-rubinstein-salzedo-thorne-ohno-nakagawa.pdf).
Here the local power calculation in PA1 achieves complete
splitting at \(\ell\) on the prescribed prime-path classes.

Adjoin the two \(\ell\)-th radicals to \(B\), giving an
extension \(N/B\) of degree \(\ell^2\). Independence
persists because \([B:K]=m\) is prime to \(\ell\).
Away from \(\ell\) the virtual valuations imply
unramifiedness. At \(\ell\) both radicals already
exist locally. For \(g\in\operatorname{Gal}(B/\mathbb Q)\)
the action on the marked radicals modulo powers is by
the quadratic character \(\chi_K(g)\), since their
quadratic norms are powers. Hence the action on the
Kummer Galois group is by the scalar
\(\chi_\ell(g)\chi_K(g)\), where \(\chi_\ell\) is the
cyclotomic character. In particular simultaneous
conjugation on \(K\) and \(\zeta_\ell\) acts trivially.
There is a unique order-two lift of that element,
because its preimage is abelian of order \(2\ell^2\).
It is central, and its fixed field \(W\) is totally real.
The induced \(W/M\) remains unramified and split above
\(\ell\) by the coprime quadratic descent.
The quotient by the simultaneous involution acts by
the faithful scalar character. A lift of a generator
of this cyclic quotient has order \(m\), since its
\(m\)-th power lies in the scalar module and is fixed
by a nontrivial scalar. Thus the affine semidirect
product and its indicated complements exist.
The \(\ell+1\) hyperplanes in the two-dimensional
normal subgroup give distinct normal closures and
hence distinct degree-\(\ell\) fields. The complement
on the rank-two affine action is core free, so its
fixed field has normal closure \(W\).

At \(p\ne\ell\) inertia acts exactly as the quadratic
inertia from \(K\), since \(B/K\) and \(N/B\) are
unramified there. If it is nontrivial, the affine
involution has one fixed point and otherwise two-cycles.
The same order-two subgroup controls each nontrivial
higher ramification group at \(p=2\). More precisely the
upper filtration maps isomorphically to that of the
quadratic extension on inertia. Hence the discriminant
exponent is \((\ell-1)/2\) times the quadratic exponent
on the degree-\(\ell\) action, and \((\ell^2-1)/2\) times
it on the rank-two action.

At \(\ell\), the prime has local squareclass \(2\ell\)
in \(K\). The only tamely ramified quadratic subfield of
the cyclotomic extension has squareclass
\(\ell^*=(-1)^{(\ell-1)/2}\ell\).
On inertia the two quadratic characters therefore have
the same restriction, the order-two character of
\(\mathbb F_\ell^\times\). The scalar on the affine part
is \(r\mapsto r^{1+m/2}\), with image of order \(e_\ell\).
Its number of orbits on the rank-\(j\) affine space,
\(j=1,2\), is \(1+(\ell^j-1)/e_\ell\).
The tame conductor exponent is
\((\ell^j-1)(1-1/e_\ell)\).
This gives (1) and the exponent ratio on passing
from degree \(\ell\) to degree \(\ell^2\).

For (2), when \(\ell\equiv1\bmod4\) the mirror inertia
already has order \(m\). When \(\ell\equiv3\bmod4\),
\(K_\ell\) coincides with the local cyclotomic quadratic
field exactly if \((-2/\ell)=1\), i.e. if
\(\ell\equiv3\bmod8\).
In that case the mirror decomposition group has order
\(m/2\), since simultaneous conjugation lies in the
local cyclotomic group realized as the graph of its
quadratic character. Otherwise the composite also
contains the unramified quadratic extension, and the
mirror decomposition group has order \(m\).
There is no change in the decomposition group on
passing to the unramified extension \(W/M\), because
these primes split there. Each resulting scalar group
fixes the zero point and acts freely on the nonzero
points of the degree-\(\ell\) permutation set.
This yields (2).

At either auxiliary prime \(q_i\) from HP5, \(B\)
splits completely. The two residue characters are the
Frobenius coordinates at a prime above it in \(W/M\).
Both rows together span its Galois group, so the
corresponding split prime classes generate that
elementary abelian quotient of the mirror class group.

## Polynomial representatives and a quintic fixture

**Proposition PA3.** For any nonzero line in the
\(\mathbb F_\ell\)-span of \(\alpha_1,\alpha_2\), choose a
representative \(\gamma=\alpha_1^a\alpha_2^b\) with
\(0\le a,b<\ell\). Define \(n=(p_1p_2)^a(p_2p_3)^b\).
A defining polynomial of its degree-\(\ell\) field in
PA2 is
\[
 \mathcal D_\ell(T,n)-\operatorname{Tr}_{K/\mathbb Q}(\gamma),
 \qquad
 \mathcal D_0=2,\quad \mathcal D_1=T,\quad
 \mathcal D_i=T\mathcal D_{i-1}-n\mathcal D_{i-2}.
\]
Indeed for \(\rho^\ell=\gamma\) the roots are
\(\zeta_\ell^j\rho+\zeta_\ell^{-j}n/\rho\).
They are distinct: equality for two indices implies
\(\rho^2=n\zeta_\ell^{-i-j}\), hence
\(\gamma^2=n^\ell\) and \(\gamma\in\mathbb Q\).
A rational element whose square is an
\(\ell\)-th power is itself an \(\ell\)-th power,
contrary to independence.
The recurrence gives the trace identity, and the affine
action on the \(\ell\) distinct roots is transitive.
The polynomial is therefore irreducible and gives the
claimed fixed field.

For \(\ell=5,\ (p,h)=(419,150)\), the three primes
are \(269,419,569\) and \(D=-16290203579467969335\).
The two rows at \(11,71\), in the indicated bases
\(4,54\), are \((1,3),(2,4)\) modulo five.
For the first mark the polynomial is
\[
 T^5-563555T^3+63518847605T+3288814754552 .
\]
Its field discriminant is
\(125(|D|/5)^2\), and its two primes above \(5\)
have types \((1,1),(4,1)\). The
[exact checker](verify_higher_prime_path_affine_reflection.py)
checks the genuine marks, residue rows, primary local
formula and polynomial discriminant calculation.

## Counts by discriminant

**Corollary PA4.** For every fixed odd prime \(\ell\ge5\)
there is a fixed \(B_\ell\) such that the number of
degree-\(\ell\) fields of discriminant at most \(Y\)
arising in PA2 is at least
\[
 \gg_\ell Y^{2/(\ell-1)^2}/(\log Y)^{B_\ell}.
\]
The corresponding degree-\(\ell^2\) fields number at least
\[
 \gg_\ell Y^{2/((\ell-1)(\ell^2-1))}/(\log Y)^{B_\ell}.
\]
For the upper count of these degree-\(\ell^2\) fields,
a bound \(O_{\ell,\epsilon}(X^{\Theta+\epsilon})\) for
the negative quadratic \(\ell\)-torsion second moment
gives \(O_{\ell,\epsilon}(Y^{2\Theta/(\ell^2-1)+\epsilon})\).
At \(\ell=5\) the retained second moment exponent
\(3/2\) gives \(1/8+\epsilon\) in this upper count.

The lower counts follow by substituting the discriminant
laws into PA1. Distinct \(D\) in this construction give
distinct discriminants, since \(\ell\mid D\) throughout.
For the upper bound fix \(K\). The number of marked
two-dimensional subspaces is at most
\(|\operatorname{Cl}(K)[\ell]|^2\) by the virtual-unit
map, since the units of an imaginary quadratic field
are all \(\ell\)-th powers for \(\ell\ge5\).
Each subspace defines at most one normal closure and
one isomorphism class of fixed fields of complements.
The latter follows by conjugating the complements by
translations in the affine group. The discriminant law
then makes the hypothesis apply at
\(X\asymp_\ell Y^{2/(\ell^2-1)}\).
The [higher-prime moment theorem](ABC_COHERENT_PRIME_PATH_HIGHER_ODD_TORSION.md)
supplies the stated exponent for \(\ell=5\).

## Frobenius at the moving prime marks

**Theorem PA5.** Impose \(p\equiv1\bmod\ell\) in PA1.
At the middle prime \(p_2\) the reflected \(W/M\)
splits completely. At the two outer prime marks, the
Frobenius rows on the marked radicals are given by
the \(\ell\)-th power characters of the following fixed
rational numbers, up to simultaneous scaling of each
row by a nonzero scalar:
\[
 \begin{array}{c|cc}
 &\alpha_1&\alpha_2\\ \hline
 p_1&2/(4-2^\ell)&1/2\\
 p_2&1&1\\
 p_3&1/2&2/(4-2^\ell).
 \end{array}
 \tag{3}
\]
For either outer prime prescribe independently a zero
row or any projective row in
\(\mathbb P^1(\mathbb F_\ell)\).
The lower bounds of PA1 and PA4 still hold on the
restricted family. For a line \((a,b)\) the associated
degree-\(\ell\) field splits completely at a mark exactly
when its row pairs to zero with \((a,b)\), and is inert
otherwise.

The three primes split completely in \(B\), since
they split in \(K\) and are \(1\bmod\ell\).
Write \(u_1=(A+B-C)/2,\ u_2=(B+C-A)/2\).
At \(p_1\) choose the split branch on which
\(\alpha_1\equiv0\); its conjugate is \(u_1\).
The normalized unit \(\alpha_1/p_1^\ell\) reduces to
\(p_2^\ell/u_1\), while \(\alpha_2\) reduces to
\((C-A)/2=C/2\) there. The expression for \(u_1\) at
that prime is \((4-2^\ell)h^\ell/2\).
For \(p_3\) interchange the consecutive roles and
replace \(h\) by \(-h\), giving the other outer row.
At \(p_2\) one has \(C=-A\) and \(u_1=A,\ u_2=-A\).
On the branch where \(\alpha_1\) vanishes,
\(\alpha_1/p_2^\ell\) reduces to \(1\) and \(\alpha_2\)
reduces to \(-A=C\). These are both powers. Changing
the split branch inverts its row, and changing the
generator of roots of unity scales its exponents.
This proves (3), including the normalization at the
prime where the original generator vanishes.

Set \(a_\ell=2/(4-2^\ell)\) and \(b_\ell=1/2\).
Their rational classes modulo \(\ell\)-th powers are
independent: \(a_\ell/b_\ell=(1-2^{\ell-2})^{-1}\)
has an odd prime valuation not divisible by \(\ell\),
since \(1<2^{\ell-2}-1<2^\ell\),
while \(b_\ell\) has valuation \(-1\) at \(2\).
The Kummer extension adjoining both radicals to
\(\mathbb Q(\zeta_\ell)\) therefore has scalar Galois
action and maximal abelian subfield
\(\mathbb Q(\zeta_\ell)\). Each zero row or projective
row is exactly a conjugacy class in the scalar kernel.
Impose those two classes respectively on the first
and last prime forms, switching the order of the
characters for the last one as in (3).
Use [Frei–Tinková, Theorem 1.4 and Corollary 1.6](https://arxiv.org/pdf/2408.17441).
The restrictions of both conjugacy classes to the
maximal abelian subfield are trivial, and the
congruence \(p=1,\ h=0\bmod\ell\) already gives that
restriction. All fixed-prime Euler factors are positive
by the existing CRT choices, and at every other prime
the forms are units under \(p=1,h=0\) modulo that prime.
The same positive-density prime form count and the
fixed-field multiplicity and prime-factor estimates of
HP5 apply. In the affine actions a nonzero translation
is an \(\ell\)-cycle, and the zero translation acts
trivially. This gives the degree-\(\ell\) splitting claims.

For example, with \(\ell=5,\ (p,h)=(421,150)\) the three
marks are \(271,421,571\). The first reflected field
has polynomial
\(T^5-570455T^3+65083781405T+3167699018548\).
The mark \(271\) is inert and both \(421,571\) split
completely in this field, as in (3).
