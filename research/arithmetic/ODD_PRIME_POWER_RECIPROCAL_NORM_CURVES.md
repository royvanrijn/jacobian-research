# Two norm marks at every odd prime power

The two moving-prime directions and the rank on
successive reciprocal Jacobian factors extend from
powers of three to every odd prime power. The
cyclotomic mirror has degree \(\varphi(m)\).
At prime powers whose prime is \(1\bmod4\), its
unique quadratic subfield is the cyclotomic
quadratic field, so the discriminant count must
use the full mirror rather than its quadratic
subfield.

Fix an odd prime \(\ell\), \(n\ge1\) and
\(m=\ell^n\). Put \(L=2n+3\),
\(a=2\ell^L\), and \(A=a^m\).
Take positive coprime squarefree integer
coefficients \(p>q>1\) with
\[
 p-q=4A,\qquad p,q\equiv1-\ell^{2n+2}\pmod{\ell^L},
 \qquad p,q\text{ odd}.
\]
Write \(B=-(p+q)/2\) and use the even norm
curve \(R^2=F(u)=A^2u^{2m}+Bu^m+1\) with
\(f_\pm=Au^m\pm1+R\), norms \(p u^m,q u^m\)
and ratios \(\gamma_\pm=f_\pm/\bar f_\pm\).
The roles of \(p,q\) are fixed coefficients.
The moving numerator and denominator
\(P,Q\) are prime.

## Independent prime and mirror directions

**Theorem OPR1.** For each fixed \(\ell,n\)
and either real sign of \(F\), genuine ratios
\(u=P/Q\) can be taken with \(2\) either
split or inert and \(\ell\) split in
\(K=\mathbb Q(\sqrt{F(u)})\).
The two \(m\)-virtual ratios are local
\(m\)-th powers above \(\ell\).
Both moving primes split in \(K\), and
their two \(m\)-power residue vectors can
be prescribed independently up to a unit
scalar on each. A basis prescription gives
independent order-\(m\) virtual classes.
On the imaginary side, for \(|D|>4pq\),
the actual two prime classes have order
\(4m\) and jointly generate
\(C_{4m}\times C_{2m}\), with independent
\(\ell\)-primary parts and the exact
ramified ideal relations of TP3.
For each fixed choice of coefficients,
either sign and either dyadic type gives
\(\gg_{\ell,n,p,q} X^{1/m}/(\log X)^C\)
distinct fundamental discriminants
\(|D|\le X\).

The squarefree coefficient choice and
arbitrary additional good-prime residues
are available by the elementary sieve
of [ENJ4](EVEN_NORM_CURVE_PRIMITIVE_JACOBIAN_TWISTS.md).
At two use the odd residue \(u\) for which
\(F(u)\equiv4\) or \(20\bmod32\).
At \(\ell\) choose \(u\equiv1\bmod\ell^L\).
Then \(F(u)\equiv\ell^{2n+2}\bmod\ell^L\)
is a nonzero local square, with
\(v_\ell(R)=n+1\). Both norm ratios lie in
\(1+\ell^{n+1}\mathbb Z_\ell\), the
\(m\)-th powers in \(1+\ell\mathbb Z_\ell\).
At each prime dividing \(p\) or \(q\)
use \(u\equiv-1/a\) or \(1/a\),
respectively. The relevant norm identity
gives ramification of valuation one.
All other valuations follow as in TP1:
outside the fixed coefficient support
one conjugate in each norm has valuation
zero whenever \(u\) has a zero or a pole,
and the other valuation is divisible by \(m\).
At two the coefficient norms are odd,
and at the ramified coefficients norm-one
values are units.

For the genuine ratios the two frames
are still
\[
              (4/p,q/4),\qquad (4A^2/p,4A^2/q).
\]
A prime above a divisor of each coefficient
gives a unimodular valuation matrix over
\(\mathbb Q(\zeta_m)\). Each governing
Kummer field has translation group \(C_m^2\)
and maximal abelian field \(\mathbb Q(\zeta_m)\),
because complex conjugation inverts the
translations and \(m\) is odd. As in
[TP2](TWO_MOVING_PRIME_CLASSES_FROM_EVEN_NORM_CURVES.md),
the frame prescriptions are therefore
compatible with all fixed rational
progressions, including \(P,Q\equiv1
\bmod\ell^L\). Chebotarev supplies a
positive prime density in each frame.
The split-prime residues give independent
classes modulo \(\ell\)-th powers.
The same valuation support gives (6)–(7)
of TP3, with auxiliary ramified ideal
products of squarefree norms \(p,q\).
Their independence for large imaginary
\(|D|\) follows from the imaginary norm
bound. Units are \(m\)-divisible there.
The squarefree binary form has degree
\(2m\), and the fixed-genus counting
argument of TP2 applies to the same
fixed rectangles and the two real sign
intervals. The statement concerns fixed
coefficients and fixed \(m\).

**Theorem OPR2.** Write \(D=\operatorname{disc}K\)
on this slice, and put \(C=\mathbb Q(\zeta_m)\).
In \(KC\) let \(\tau\) act as the nontrivial
quadratic involution on \(K\) and inversion
on \(\zeta_m\). The fixed field \(E=(KC)^\tau\)
is cyclic of degree \(\varphi(m)\) over
\(\mathbb Q\). It supports an everywhere
unramified \(C_m^2\) quotient of its
ordinary class group, split above \(\ell\),
with the prescribed Frobenius vectors of
the two actual moving primes. Its
discriminant is
\[
 |\operatorname{disc}E|=
   \ell^{(n(\ell-1)-1)\ell^{n-1}}
        |D|^{\varphi(m)/2}.
 \tag{1}
\]
The map from these \(D\) to \(E\) is injective,
and each sign and dyadic type thus gives
\(\gg Y^{\,2/(m\varphi(m))}/(\log Y)^C\)
distinct \(E\) of discriminant at most \(Y\).

The quadratic field \(K\) is disjoint from
\(C\) since \(K\) splits at \(\ell\) and
the unique quadratic cyclotomic subfield
ramifies there. If \(\gamma\in K\) is a
non-\(\ell\)-th power, it stays so in \(KC\).
Indeed over \(K(\zeta_\ell)\), where the
degree over \(K\) is prime to \(\ell\),
it is still a nonpower by taking norms.
Then \(X^\ell-\gamma\) is irreducible
there and over \(K\). A root field of
degree \(\ell\) over \(K\) would have to be
normal if it were contained in the
abelian extension \(KC/K\). It would
then contain \(\zeta_\ell\), which has
degree \(\ell-1>1\) over \(K\).
It follows that the two global marks
give all \(m^2\) translations in the
Kummer extension of \(KC\).

The involution \(\tau\) lifts by inverting
the two \(m\)-th roots of the ratios.
It centralizes translations, since it
inverts both those roots and \(\zeta_m\).
Its fixed field gives the desired
relative \(C_m^2\)-extension of \(E\).
Virtuality at primes away from \(\ell\)
gives unramifiedness there, and the
local powers in OPR1 give complete
splitting at \(\ell\). Descent through
the degree-two base change preserves
these properties for the odd extension.
The moving primes split completely in
\(KC\) and retain their prescribed
relative Frobenius after descent.

The quotient group of \(KC/\mathbb Q\)
by this involution maps to
\((\mathbb Z/m)^\times\) by multiplying
the cyclotomic coordinate by the
quadratic sign, so is cyclic.
At \(\ell\), the two places of \(KC\)
are exchanged by \(\tau\), with
completion \(\mathbb Q_\ell(\zeta_m)\).
There is one place of \(E\) above \(\ell\)
with the full cyclotomic completion.
At an odd ramified prime of \(K\) away
from \(\ell\), inertia stays quadratic
in \(E\). At all other primes away from
\(\ell\) it is unramified. This gives (1).
Moreover \(EC=KC\). Of the two quadratic
subfields in this compositum outside
\(C\), exactly one is split at \(\ell\)
because the other is its product with
the quadratic cyclotomic subfield.
This recovers \(K\) from \(E\) on the
stated slice. The unique quadratic
subfield of \(E\) is
\(\mathbb Q(\sqrt{\ell})\) for
\(\ell\equiv1\bmod4\), and
\(\mathbb Q(\sqrt{-\ell D})\) for
\(\ell\equiv3\bmod4\).
At \(m=\ell\) this is the classical
mirror field of degree \(\ell-1\) in
Cohen–Rubinstein-Salzedo–Thorne §2 and
Lemma 2.11. Here the exact
prime-power extension carries two
independent moving-prime vectors.
[Cohen–Rubinstein-Salzedo–Thorne](https://thornef.github.io/ohno-nakagawa3.pdf).

## The geometric criterion for new rank

**Theorem OPR3.** Use the models \(C_k,J_k,N_k\)
and reciprocal factors \(B_{k,\pm}\) of ENJ1
with \(r_k=\ell^k\) in place of \(3^k\).
The two factors are simple over \(\mathbb Q\)
of dimension \(d_k=\varphi(r_k)/2\).
Each remains simple over any number field
having a dyadic place of ramification
degree prime to \(\ell\).
The \(\ell\)-primary torsion of \(N_k\)
over a field with such a place is
exactly the rational \(C_\ell^2\); in
particular every nontrivial quadratic
twist has zero rational \(\ell\)-torsion.

For any nontrivial quadratic specialization
of \(C_n\), each \(\gamma_\pm\) that is
a non-\(\ell\)-th power gives positive
rank on every corresponding
\(B_{k,\pm}^{(D)}\). This criterion holds
without a dyadic splitting requirement
on the quadratic specialization.
Thus the two genuine frames in OPR1
give rank at least \(2n\) on the
whole twisted Jacobian on both dyadic
types and both signs. The untwisted
\(J_n\) has full rational torsion
\(C_m^2\), and every nontrivial quadratic
twist of \(J_n\) has trivial rational
torsion.

The dyadic model of TP4 holds for odd
\(r_k\) and \(r_k\mid v_2(A)\).
Its two genus-\((r_k-1)/2\) components
have primitive inertia. On each, invariants
of the subgroup generated by the
\(\ell^j\)-th power have dimension
\(\ell^j-1\) by taking the quotient
curve of degree \(\ell^j\) over the
\(x\)-line. The successive rational
\(\ell\)-adic representation consequently
has characteristic polynomial
\(\Phi_{r_k}^2\). The two reciprocal
eigenspaces split its dimension evenly
by the same Dickson quotient formula
as in ENJ1. Over a number field with the
specified dyadic place, inertia has the
same image because its index divides
the dyadic ramification degree.
Irreducibility of \(\Phi_{r_k}\) over
\(\mathbb Q_\ell\) proves simplicity
over that field.

The pullback of the preceding Jacobian
is injective by the ramified prime-degree
argument of ENJ1, with \(\ell\) in place
of three. The inertia-fixed geometric
\(\ell\)-primary torsion on \(J_k\) is
the rational boundary group \(C_{r_k}^2\)
by TP4. Its image in \(N_k\) has order
\(\ell^2\), and exhausts the fixed group
by the determinant \(\Phi_{r_k}(1)^2\).
It is \(C_\ell^2\) by pullback of the
boundary divisors. A quadratic field
has dyadic ramification degree at most
two, so its untwisted \(\ell\)-torsion
on \(N_k\) is still all rational.
Consequently that torsion has no
quadratic anti-invariant part, even if
the quadratic field splits at two.
This proves the global torsion
exclusion required for every quadratic
specialization.

Evaluation of \(f_{\pm,k}\) gives a
homomorphism modulo \(\ell\)-th powers,
since its divisor is \(r_k\) times a
boundary divisor. The functions are
pullbacks from the preceding model
and evaluate to \(\ell\)-th powers on
pulled-back divisor classes. The
reciprocal sum and difference of
\([P_k-\iota P_k]\) give \(\gamma_\pm^2\)
under the respective evaluation as in
ENJ2. If a corresponding image in the
new quotient were torsion, one could
multiply by its order, prime to \(\ell\)
by the global torsion exclusion, to land
in the preceding Jacobian. This
contradicts the nonpower evaluation.

To check the full torsion assertion,
it remains to exclude rational
two-torsion on \(J_n\).
Over \(\mathbb Q(\sqrt{pq})\) put
\(\beta=(-B+\sqrt{pq})/(2A)\).
If \(\theta^\ell=\beta\) there,
\(N(\theta)=1\) and \(c\theta\) is
integral for \(c=a^{m/\ell}\).
Its trace \(T\) is integral and
\(D_\ell(T,c^2)=-B\), where
\(D_\ell(c(w+w^{-1}),c^2)=
c^\ell(w^\ell+w^{-\ell})\).
The identity
\[
          q=D_\ell(T,c^2)-2c^\ell
           =(T-2c)\,G_\ell(T,c)^2
\]
has \(G_\ell(T,c)\in\mathbb Z[T,c]\):
upon putting \(T=c(w+w^{-1})\),
it becomes
\(c^{(\ell-1)/2}\sum_{i=-(\ell-1)/2}^{(\ell-1)/2} w^i\).
Since \(q>0\), \(T>2c\). Choosing
\(w>1\) proves \(|G_\ell(T,c)|>1\),
contradicting squarefreeness of \(q\).
Thus \(\beta\) is a non-\(\ell\)-th
power. The binomial \(x^{\ell^n}-\beta\)
is irreducible over the quadratic field
by induction: the norm of a root at
each irreducible odd level is \(\beta\)
and hence that root is still a
nonpower. At each prime step,
irreducibility follows by adjoining
\(\zeta_\ell\) of prime-to-\(\ell\)
degree and applying Kummer theory.
The branch polynomial over \(\mathbb Q\)
is therefore irreducible of degree \(2m\).
No invariant even subset of its
branch points represents nonzero
rational two-torsion, by TP4. This
eliminates two-torsion also on every
quadratic twist. The dyadic
inertia-fixed odd torsion is all rational
on the untwisted Jacobian and remains
the full odd torsion over any quadratic
field, since the inertia image has
odd order. It cannot contribute
anti-invariant torsion to a nontrivial
twist.

## Nonisogeny and exact common conductor

**Theorem OPR4.** The fixed squarefree
coefficients can be chosen so that all
\(2n\) simple factors are pairwise
nonisogenous over \(\mathbb Q\).
On every genuine slice of OPR1 the
two factors at level \(k\) have the
common \(D\)-twist conductor
\[
 N(B_{k,\pm}^{(D)})
 =2^{2d_k}\ell^{(k(\ell-1)-1)\ell^{k-1}+1}
                     |D|^{2d_k}.
 \tag{2}
\]
The untwisted conductor uses
\((pq)^{d_k}\) in place of \(|D|^{2d_k}\).
For fixed \(\ell,n,k\) and coefficients,
both signs and both dyadic types give
\(\gg Z^{1/(2md_k)}/(\log Z)^C\)
distinct common conductors up to \(Z\),
with simultaneous positive rank and
the mirror quotient.

For the finite trace test in ENJ4 take
a distinct good prime \(s_k>4r_k+3\)
with \(s_k\equiv-1\bmod r_k\).
Its norm-one group has the nested
odd-index subgroups \(U^{r_k}\) and
\(U^{r_k/\ell}\). The primitive trace
difference is the difference of (3)
for \(j=r_k\) and \(j=r_k/\ell\).
The two ratios to the quadratic
character difference in the
orthogonality argument are now
\(-r_k+r_k/\ell\) and \(r_k/\ell\).
They are distinct, so some
\(\lambda_k\ne\pm2\) separates the
two traces. Use these finitely many
residues of \(B/A\) in the fixed
squarefree coefficient sieve. Different
levels have different dimensions.
Simultaneous twisting preserves
rational nonisogeny.

At two the conductor computation
in ENJ3 is tame. At all odd primes
different from \(\ell\), its
totally toric charts at the prime
divisors of \(pq\), and good Jacobian
charts elsewhere, apply with the
new exponent \(r_k\). At \(\ell\)
the unit root \(\alpha\) of the
quadratic in \(u^{r_k}\) is
\(1\bmod\ell^{2n+2}\). Its two branch
scales have \(r_k\)-th roots over
\(\mathbb Q_\ell\), and each block of
branch roots is a scalar times
\(\mu_{r_k}\). In the cyclotomic
splitting field, the within-block
depths are \(\ell^{k-j}\), with the
same even shift for the other block
as in ENJ3. The proper clusters of
size \(\ell^j\) have \(\ell\) odd
children, and their \(\nu\)'s are
odd. The last ramified quadratic
extension gives compact type and
an inertia element acting by
negation on every component.
Thus the conductor tame term on
each full reciprocal half is
\(r_k-1\). The wild term from one
block is
\(\sum_{j=1}^k
 [((j-1)(\ell-1)-1)\ell^{j-1}+1]\),
by the cyclotomic root-field sum of
Dokchitser–Dokchitser–Maistret–Morgan.
The reciprocal map swaps the two
blocks, so the two full factors
have the same conductor at \(\ell\).
Away from \(\ell\), the isogeny
given by the difference of the two
cyclic rotations is defined over a
cyclotomic extension unramified
at the prime, giving the same
conductor equality. Subtraction
at the preceding level gives the
exponent of \(\ell\) in (2), namely
\[
       (r_k-r_{k-1})
       +((k-1)(\ell-1)-1)\ell^{k-1}+1
           =(k(\ell-1)-1)\ell^{k-1}+1.
\]
The twist is unramified at two,
trivial at \(\ell\), and ramified
quadratic at every other prime
dividing \(D\). At those primes
the full \(2d_k\) is the tame
conductor exponent, including at
the primes dividing \(pq\).
This gives (2) and the count by
substitution in OPR1.

The prime-power extension uses fixed
coefficients for every \(\ell,n\).
The [exact checker](verify_odd_prime_power_reciprocal_norm_curves.py)
replays higher-prime local and moving
residues, cyclotomic and trace identities,
and the Dickson square factor used
in the branch irreducibility proof.

## Rational endomorphisms and minimum cross-isogeny degrees

**Theorem OPR5.** The nonisogeny over
\(\mathbb Q\) in OPR4 holds for every
coefficient pair in OPR1, without the
extra trace prescriptions. Fix a level
\(k\), set \(C_k^{\mathrm{cyc}}=
\mathbb Q(\zeta_{r_k})\), and write
\(C_k^+\) for its real subfield. For a
number field \(F_0\) with some dyadic
place of ramification degree prime to
\(\ell\), put \(F_1=F_0\cap C_k^{\mathrm{cyc}}\).
Then
\[
 \begin{split}
 \operatorname{End}_{F_0}^0(B_{k,\pm})
                &= F_1\cap C_k^+,\\
 \operatorname{Hom}_{F_0}^0(B_{k,+},B_{k,-})
                &\simeq \{\theta\in F_1:
                                    \bar\theta=-\theta\}.
 \end{split}
 \tag{3}
\]
In particular the opposite factors are
isogenous over \(F_0\) exactly when
\(F_1\) is nonreal. In that case the
minimum isogeny degree in either
direction is
\[
                        \ell^{[C_k^{\mathrm{cyc}}:F_1]} .
 \tag{4}
\]
Each factor has a rational cyclic
degree-\(\ell\) quotient from its
boundary torsion line. Over the
full cyclotomic field this quotient
is isomorphic to the other factor.
All conclusions can also be obtained
with both fixed coefficients prime,
allowing \(a\) to be a positive
multiple of \(2\ell^{2n+3}\).

Use the cyclic automorphism \(\rho\)
on the normalized curve at level \(k\).
Over \(C_k^{\mathrm{cyc}}\) it
descends to \(N_k\). Its rational
\(\ell\)-adic characteristic polynomial
on \(N_k\) is \(\Phi_{r_k}^2\):
the fixed-space dimensions for
subgroups follow by taking the cyclic
quotient curves and subtracting the
preceding level. Thus the rational
endomorphism \(\rho+\rho^{-1}\)
generates the real cyclotomic field
on each reciprocal factor. Its
minimal polynomial there is the
real cyclotomic polynomial, since
\(\Phi_{r_k}\) has no zero root and
the substitutions \(t+t^{-1}\)
pair its \(2d_k\) roots. In
particular both factors carry
\(C_k^+\) over the full cyclotomic
field. This cyclic correspondence
is the one in
[Tautz–Top–Verberkmoes, §3.1](https://research.rug.nl/files/1054394337/explicit-hyperelliptic-curves-with-real-multiplication-and-permutation-polynomials.pdf).

Over \(F_0 C_k^{\mathrm{cyc}}\) the dyadic
ramification degree is still prime
to \(\ell\), since an odd
cyclotomic field is unramified at
two. Each simple factor therefore
has a division endomorphism algebra
by OPR3. That algebra commutes
faithfully with dyadic inertia on
the rational \(\ell\)-adic Tate module.
The centralizer of its dyadic
generator is a field, by the
irreducible characteristic
polynomial in OPR3, so the
endomorphism algebra is commutative.
At a prime dividing \(p\) or \(q\)
both factors have totally toric
semistable reduction by ENJ3 and
OPR4. Their endomorphism field
acts faithfully on the rational
character lattice of the reduction
torus of dimension \(d_k\).
Indeed any nonzero endomorphism is
an isogeny and has an inverse up
to multiplication. Both maps extend
to the Néron model, and their
composition induces that same
nonzero multiplication on the
torus. Thus the endomorphism
field has degree at most \(d_k\).
It already contains \(C_k^+\) of
degree \(d_k\), so equals it.
Every such endomorphism is
already defined over
\(C_k^{\mathrm{cyc}}\).

The map \(\delta=\rho-\rho^{-1}\)
interchanges the two factors.
On \(N_k\), its kernel is that of
\(\rho-1\) since \(\rho+1\) has degree
one, from the cyclotomic polynomial
at \(-1\). The kernel has order
\(\ell^2\) from its value at one.
The two rational boundary directions
survive as \(C_\ell^2\) on \(N_k\)
and are fixed by \(\rho\), so give
the exact kernel. Their signs under
the reciprocal involution are plus
and minus as in ENJ2. Consequently
the restrictions of \(\delta\) in
both directions have cyclic kernel
of order \(\ell\), defined over
\(\mathbb Q\). This proves the
statement about the rational
degree-\(\ell\) quotients.

All rational cross-homomorphisms
over \(F_0 C_k^{\mathrm{cyc}}\)
therefore have the unique form
\(\delta\circ b\) for \(b\in C_k^+\).
Writing \(\zeta=\zeta_{r_k}\),
identify this line with the
anti-real part of
\(C_k^{\mathrm{cyc}}\) by
\(\delta\circ b\mapsto
(\zeta-\zeta^{-1})b\).
A cyclotomic automorphism
\(\sigma_a\) takes \(\delta\)
to \(\rho^a-\rho^{-a}\), which
multiplies it by the real element
\((\zeta^a-\zeta^{-a})/
(\zeta-\zeta^{-1})\). Thus the
identification is equivariant.
The same action on
\(\rho+\rho^{-1}\) gives the usual
action on the endomorphism field.
Taking invariants under
\(\operatorname{Gal}(F_0 C_k^{\mathrm{cyc}}/F_0)\)
proves (3).

For a nonzero \(F_0\)-isogeny written
as \(\delta\circ b\), its degree is
\(\ell\,|N_{C_k^+/\mathbb Q}(b)|^2\).
This follows by taking determinants
on rational homology, which has
dimension two over \(C_k^+\).
At the unique prime above \(\ell\)
in \(C_k^{\mathrm{cyc}}\), the
valuation of its image
\(\theta=(\zeta-\zeta^{-1})b\)
is therefore \(v_\ell(\deg(\delta\circ b))\).
Here \(v(\zeta-\zeta^{-1})=1\)
and \(v(b)=2v_\ell(N_{C_k^+/\mathbb Q}(b))\).
The degree is a positive integer,
so this exponent is positive and
odd. Since \(\theta\in F_1\)
and the cyclotomic prime is
totally ramified, the exponent
is divisible by \(e=[C_k^{\mathrm{cyc}}:F_1]\).
This proves the lower bound (4).

For attainment, a nonreal \(F_1\)
has odd index \(e\). Take
\(\theta=N_{C_k^{\mathrm{cyc}}/F_1}
(\zeta-\zeta^{-1})\), which is
anti-real and has valuation \(e\)
and no prime divisor away from
\(\ell\). Then
\(b=\theta/(\zeta-\zeta^{-1})\)
is real and integral, a product
of the nonidentity conjugates in
that norm. It lies in
\(\mathbb Z[\zeta+\zeta^{-1}]\):
write any element of
\(\mathbb Z[\zeta]\) in the form
\(f(s)+\zeta g(s)\) with
\(s=\zeta+\zeta^{-1}\) and compare
with its conjugate. Hence \(b\)
acts integrally by the established
cyclic endomorphism. The map
\(\delta\circ b\) is an actual
isogeny, descends to \(F_0\) by
(3), and has degree exactly
\(\ell^e\). This also explains
why the two rank directions at
each level have the same local
conductor but belong to distinct
rational isogeny classes.

For prime coefficients use the
positive-relative-density progression
\(1-\ell^{2n+2}mod\ell^{2n+3}\)
and the two polynomials \(0\) and
\(4(2\ell^{2n+3}t)^m\) in
[Tao–Ziegler, Theorem 5](sources/tao-ziegler-narrow-progressions-primes.pdf).
The coefficient arguments of OPR1–OPR4
only require the stated divisibility
of \(a\); if additional primes divide
\(a\), impose that \(u\) is a unit
at those primes in OPR1. Both fixed
coefficients are coprime to \(a\).
The good Jacobian calculation at
an odd prime dividing \(A\) already
uses only \(A=a^m\). The nonisogeny
for this prime-coefficient choice
follows from (3), without prescribing
auxiliary good-prime reductions.

## Absolute simplicity and unrestricted isogeny descent

**Theorem OPR6.** The primitive reciprocal
factors are absolutely simple and their full
geometric endomorphism rings are
\(\mathbb Z[\zeta_{r_k}+\zeta_{r_k}^{-1}]\).
In either direction the integral geometric
cross-homomorphisms form a free rank-one
module over this ring, generated by the
degree-\(\ell\) isogeny of OPR5.
Equations (3)–(4) hold over every number
field \(F_0\), with no restriction on its
dyadic ramification. The same results hold
after simultaneous quadratic twisting.

The argument gives a general criterion.
Let \(V\) be a \(d\)-dimensional abelian
variety with potentially totally toric
reduction at one finite place. At another
place suppose its rational Tate
representation has a finite odd cyclic
inertia image with irreducible
characteristic polynomial for a
generator. Suppose further that
\(\operatorname{End}^0(V_{\overline{\mathbb Q}})\)
contains a totally real field \(H\) of
degree \(d\) fixed pointwise by that
inertia group. Then
\(\operatorname{End}^0(V_{\overline{\mathbb Q}})=H\).

For the criterion write \(D_0\) for the
geometric endomorphism algebra and take a
finite field of definition for all its
elements, enlarged so the first place
has totally toric reduction. There is
a faithful action of \(D_0\) on the
\(d\)-dimensional rational cocharacter
space of that torus. To see faithfulness,
the kernel is a two-sided ideal of the
semisimple algebra \(D_0\). A nonzero
kernel would contain a nonzero rational
idempotent. Its abelian isogeny factor
has totally toric reduction since
toric and abelian ranks are preserved
under isogeny and add on products.
The idempotent induces the projection
on the torus of this nonzero factor,
contradicting a zero action.
[Grothendieck, SGA 7 I, Exposé IX, Proposition 2.2.6](https://grothendiecksga.com/read/sga7/fr/9-2.html).

The cocharacter space has dimension
one over \(H\), so the centralizer of
\(H\) in \(D_0\) is \(H\).
In particular every central idempotent
is in the field \(H\), so \(D_0\)
is simple. Its center \(Z_0\) is a
subfield of \(H\). The action of the
inertia generator on \(D_0\) fixes
\(H\) pointwise and hence fixes \(Z_0\).
By Skolem–Noether it is conjugation by
an element \(b\) of \(D_0^\times\).
That element commutes with \(H\), so
\(b\in H^\times\).
[Colliot-Thélène–Skorobogatov, Theorem 1.3.3](https://www.imo.universite-paris-saclay.fr/~jean-louis.colliot-thelene/BGgroup_book.pdf).

If the inertia order is \(h\), then
\(b^h\in Z_0^\times\), since conjugation
on endomorphisms is induced by conjugation
on the Tate representation. Each embedding
of \(H\) fixing \(Z_0\) has real image and
takes \(b\) to a real \(h\)-th root of
\(b^h\). This root is unique because
\(h\) is odd. Hence \(b\in Z_0\) and
inertia fixes \(D_0\) pointwise.
It follows that \(D_0\) is commutative
since the irreducible inertia generator
has a field centralizer on the rational
Tate module. Now \(D_0\) centralizes
\(H\) and must itself be \(H\).
This proves the general criterion.

For \(B_{k,\pm}\), use any prime
dividing \(pq\) for the totally toric place,
two for the irreducible odd inertia,
and \(H=C_k^+\) for the degree-\(d_k\)
field, defined over the odd cyclotomic
field unramified at two. The criterion
gives the full geometric endomorphism
algebra and therefore absolute
simplicity. Integral endomorphisms
lie in the ring of integers by their
action on singular integral homology.
The whole ring
\(\mathbb Z[\zeta_{r_k}+\zeta_{r_k}^{-1}]\)
already acts, by the cyclic
correspondence and the invariant-ring
argument in OPR5.

For the cross-homomorphisms, write
any geometric map as \(\delta\circ b\)
with \(b\in C_k^+\). Since \(\delta\)
has kernel killed by \(\ell\),
\(\ell\delta^{-1}\) is an integral
map, so \(\ell b\) is an integer
of \(C_k^+\). The only possible
denominator of \(b\) is at the unique
prime above \(\ell\). The degree
\(\ell |N_{C_k^+/\mathbb Q}(b)|^2\)
of a nonzero integral map has
nonnegative \(\ell\)-adic valuation,
which forces \(b\) integral at that
prime as well. Thus the integral
module is exactly the stated module
generated by \(\delta\).
All maps are already defined over
\(C_k^{\mathrm{cyc}}\). The equivariant
description and minimum-degree
argument of OPR5 now apply over
arbitrary \(F_0\). The same descent
calculation applies to a simultaneous
quadratic twist because the two
negation automorphisms commute with
every homomorphism.
