# Primitive translation observability in the binary projected return semigroup

## 1. Outcome

**Current status.**  The observability and elimination results below remain
valid, but their Hall/carry promotion programme is now parked.  The later
[Hall-envelope theorem](BINARY_GVC_ENVELOPE_CLOSURE.md) proves unrestricted
binary GVC before any primitive packet must be extracted.

The first proposed counterexample search at the primitive-semigroup level has
an all-span answer.  In the projected two-colour return configuration, every
genuinely mixed Graver move transfers a nonzero weighted level from one colour
to the other.  Translating either colour separately turns that transfer into a
strict degree gap between the two endpoint monomials.  Consequently no mixed
primitive packet can vanish along its complete one-colour Taylor orbit while
all coefficient channels remain nonzero.

This is independent of prime-power valuations and of the factorial weights.
Those invariants are still useful before a primitive move has been exposed,
but there is no adelic collision to find after exposure of one mixed Graver
binomial.

Two further all-scale reductions close the first linear shells.  Scales one
and two force every three-state orbit ratio to be constant, so no transferring
three-state shell survives.  More generally, a confluent-Wronskian argument
separates every fixed finite family of affine orbit rays

\[
 h_j(t)f_j(t)^N
\]

for all sufficiently large \(N\), even when its scalar coefficients depend
arbitrarily on \(N\).  Cancellation therefore splits by proportionality of
the exponential bases \(f_j\); nonexponential carry coefficients cannot
couple distinct affine rays in a fixed template.  They may still act inside
one common-base low-correction space.  For fixed signs, factorial trace
independence and power-sum partial fractions sharpen this to equal
opposite-sign orbit pairs.  The first six four-state pairs surviving both
Taylor directions are exactly the safe product Veronese locus.  Independently,
the first three rows already classify a four-state shell with arbitrary
nonzero row coefficients: every positive-dimensional ratio component is a
sum of two proportional pair blocks.  Once a common base has been exposed,
the entire infinite Cartesian ideal collapses exactly to the finite span of
its scale coefficient vectors; on a primitive orbit circuit, survival forces
projective constancy or a further rank drop.

The factorial part of that projective alternative is also exact.  Integer
affine factorial profiles are proportional only through
same-rational-boundary transfers, and finite nonzero carry or character
alphabets can vary only by one common projective scalar.  If the coefficient
is a finite additive factorial sum with an eventually periodic law, residue
restriction and factorial-trace independence reduce it to finitely many
rational-function identities.

The ordinary Cobham shortcut does not promote the remaining carry law.
Already for the central-binomial partition, the one-bit Kummer carry
indicator is two-state \(p\)-automatic and stationary on every ray \(qp^e\),
but is not ultimately periodic.  Therefore no common finite state automatic
in two unrelated bases can retain the complete prime-specific carry data.

The theorem does **not** prove affine-carry promotion.  A Graver decomposition decomposes binomial lattice
relations; it does not decompose a signed linear cancellation into vanishing
primitive summands.  The Wronskian theorem removes arbitrary scalar carry
weights and even changing support **between distinct bases** inside one fixed
finite affine-ray universe.  It leaves a finite common-base low-correction
circuit only when its carry law is genuinely nonperiodic automatic,
infinite-output, additive with a growing alphabet, or lies on a further
orbit-rank locus.  The positive-density Hall face may also have an unbounded
family of affine rays as the scale grows.  Classifying that residual bounded
law, or uniformly extracting a fixed finite quotient packet in the second
case, is the gap left inside the parked route.  This is narrower than the
former search for indistinguishable primitive packets, but it is not a gap in
the independent Hall-envelope proof.

## 2. The two-colour configuration

Fix \(s\geq1\).  The projected return matrix has columns

\[
 R_i\longmapsto(1,0,i),\qquad
 B_i\longmapsto(0,1,i),\qquad 0\leq i\leq s.
\tag{2.1}
\]

For a signed vector \(g=(g^R,g^B)\in\ker_{\mathbb Z}A_s\), put

\[
 m_C(g)=\sum_i g_i^C,\qquad
 \ell_C(g)=\sum_i i g_i^C\quad(C=R,B).
\tag{2.2}
\]

Thus

\[
 m_R(g)=m_B(g)=0,\qquad
 \ell_R(g)+\ell_B(g)=0.
\tag{2.3}
\]

Call \(g\) **genuinely mixed** when both colour restrictions are nonzero.
This is the relevant case after pure one-colour relations and support loss
have been removed.

> **Lemma 2.1 (a mixed primitive move transfers level).**  If \(g\) is a
> genuinely mixed Graver element of \(A_s\), then
> \[
>  \ell_R(g)=-\ell_B(g)\ne0.
> \tag{2.4}
> \]

### Proof

If \(\ell_B(g)=0\), then (2.3) also gives \(\ell_R(g)=0\).  Hence

\[
 (g^R,0),\qquad(0,g^B)
\tag{2.5}
\]

are two nonzero lattice relations for \(A_s\).  They are sign-compatible
with \(g\), and their sum is \(g\).  This is a proper conformal
decomposition, contrary to Graver primitivity.  The same argument applies
with the colours exchanged.  \(\square\)

## 3. Translation-degree separation

Use normalized Taylor variables \(B_0,\ldots,B_s\) and the locally nilpotent
derivation

\[
 D_B(B_i)=(i+1)B_{i+1}\quad(0\leq i<s),
 \qquad D_B(B_s)=0.
\tag{3.1}
\]

Its exponential is

\[
 B_i(t):=\exp(tD_B)B_i
 =\sum_{j=i}^s\binom ji B_jt^{j-i}.
\tag{3.2}
\]

On the coefficient torus \(B_s\ne0\), so

\[
 \deg_t B_i(t)=s-i.
\tag{3.3}
\]

Write \(g=g^+-g^-\), and let \(w_+,w_-\) be arbitrary nonzero scalars.
They may include the two multinomial weights, a common radial factorial,
torsion phases, and nonzero evaluations of the untranslated colour.  The
associated endpoint binomial is

\[
 F_g=w_+R^{(g^R)^+}B^{(g^B)^+}
     -w_-R^{(g^R)^-}B^{(g^B)^-}.
\tag{3.4}
\]

> **Theorem 3.1 (primitive one-colour translation observability).**  Let
> \(g\) be a genuinely mixed Graver element of \(A_s\).  There is no point
> of the coefficient torus at which
> \[
>  D_B^kF_g=0\qquad(k\geq0).
> \tag{3.5}
> \]
> The same statement holds with \(R\) and \(B\) exchanged.

### Proof

Because \(D_B\) is locally nilpotent, (3.5) is equivalent to the polynomial
identity

\[
 w_+R^{(g^R)^+}\prod_iB_i(t)^{(g_i^B)^+}
 =
 w_-R^{(g^R)^-}\prod_iB_i(t)^{(g_i^B)^-}.
\tag{3.6}
\]

The two \(B\)-selection counts agree; denote their common value by \(q_B\).
By (3.3), the degrees of the two sides are

\[
 sq_B-\sum_i i(g_i^B)^+,qquad
 sq_B-\sum_i i(g_i^B)^-.
\tag{3.7}
\]

Their difference is \(-\ell_B(g)\), which is nonzero by Lemma 2.1.  Both
leading coefficients are nonzero on the coefficient torus, so (3.6) is
impossible.  Exchanging the colours proves the second assertion.
\(\square\)

Let

\[
 I_B(g)=(F_g,D_BF_g,D_B^2F_g,\ldots).
\tag{3.8}
\]

The list is finite because \(D_B\) is locally nilpotent.

> **Corollary 3.2 (exact torus certificate).**  Over every algebraically
> closed field of characteristic zero,
> \[
> I_B(g):\left(\prod_{i=0}^sR_iB_i\right)^\infty=(1),
> \tag{3.9}
> \]
> and likewise for \(I_R(g)\).  Thus a primitive mixed packet hidden from
> either one-colour Taylor tower must lose a coefficient channel.

This follows from Theorem 3.1 and the Nullstellensatz.  Notice that no
equality or positivity assumption on \(w_+,w_-\) was used; only their
nonvanishing matters.  In particular the statement survives every
factorial normalization for which the two endpoints belong to one radial
fibre.

The proof is really graded rather than one-dimensional.  Suppose channels
\(X_\nu\) have additive weights \(\omega_\nu\), the return matrix records
the two colour counts and total \(\omega\)-weight, and one Cartesian
translation has

\[
 \deg_t\exp(tD)X_\nu=h_C-\omega_\nu
\tag{3.10}
\]

with nonzero leading coefficient on the active torus.  Lemma 2.1 and the
degree comparison prove the same conclusion verbatim.  For the Taylor
simplex of a translated monomial \(z^\gamma\), take
\(\omega_\beta=|\beta|\) and translate along a generic line with every
direction coordinate nonzero; then

\[
 \deg_t\frac{\partial^\beta(z+tv)^\gamma}{\beta!}
 =|\gamma|-|\beta|.
\tag{3.11}
\]

Thus a single genuinely mixed primitive move on a rank-one Cartesian
monomial face is also observable.  The qualification “single” remains
essential: a signed sum of several orbit rows can cancel its leading degrees
linearly, as Section 5 shows.

## 4. Why common translation is weaker

The separate one-colour operation is essential.  Under the diagonal
derivation \(D_R+D_B\), all monomials in one return fibre have the same total
translation degree, so Lemma 2.1 gives no degree separator.

There are already two diagonal survivors at span four:

\[
 R_3B_0B_3=R_1B_1B_4,
 \qquad
 R_4B_0B_3^2=R_1B_1B_4^2.
\tag{4.1}
\]

They are not new obstructions.  Put

\[
 q_i(t)=\binom4i b^i(a+bt)^{4-i}\qquad(ab\ne0).
\tag{4.2}
\]

Setting \(R_i(t)=B_i(t)=q_i(t)\) makes both identities in (4.1) hold for
every \(t\).  This is the quartic Veronese, hence the already-safe pure-power
or homogeneous locus.  Translating either colour alone destroys both
identities, exactly as Theorem 3.1 predicts.

## 5. A minimal linear-orbit obstruction

Pairwise observability does not make a signed shell inherit a primitive
binomial.  There is a full-torus three-state control already at span two.
All three states have colour counts \((2,2)\) and total level four:

\[
 S_1=R_1^2B_1^2,\qquad
 S_2=R_0R_2B_0B_2,\qquad
 S_3=R_0^2B_2^2.
\tag{5.1}
\]

Their two-colour multinomial weights are \(1,4,1\).  Mark \(R_0,R_2\)
by \(1\in C_4\) and every other channel by zero.  The resulting character
phases on (5.1) are \(1,-1,-1\).  At the coefficient-torus point

\[
 R_0=R_1=R_2=1,\qquad
 (B_0(t),B_1(t),B_2(t))=(t^2+3t+2,2t+3,1),
\tag{5.2}
\]

one has the complete translation identity

\[
 S_1(t)-4S_2(t)-S_3(t)
 =(2t+3)^2-4(t^2+3t+2)-1=0.
\tag{5.3}
\]

This is the quadratic discriminant circuit assembled linearly from safe
beta/centered moves.  It is not an all-scale packet.  There are two phase
conventions which must not be conflated.  If the displayed signs are fixed
external coefficients, the scale-two middle factorial weight is
\(\binom42^2=36\), and the next row at (5.2) is

\[
 (2t+3)^4-36(t^2+3t+2)^2-1
 =-64-216t-252t^2-120t^3-20t^4\ne0.
\tag{5.4}
\]

If instead the signs come from a fixed character, the phases themselves are
raised to the second power and the row is

\[
 (2t+3)^4+36(t^2+3t+2)^2+1
 =226+648t+684t^2+312t^3+52t^4\ne0.
\tag{5.5}
\]

Thus (5.3) is an exact obstruction to the inference “all primitive pairs are
observable, therefore the signed shell promotes.”  It is simultaneously a
safe calibration: imposing the next scaled-factorial row kills it under
either phase law.  The next theorem shows that this is a unit-ideal
obstruction on the whole coefficient torus, not merely failure at the point
(5.2).

## 6. Two scales eliminate every transferring three-state shell

The preceding example belongs to a general characteristic-zero argument
which does not use the factorial formula.  Let \(M_1,M_2,M_3\) be three
distinct monomials in one projected return fibre.  Translate the \(B\)-colour
and suppose that the scale-one and scale-two rows are

\[
 \sum_{j=1}^3a_jM_j(t)=0,
 \qquad
 \sum_{j=1}^3b_jM_j(t)^2=0,
 \qquad a_jb_j\ne0.
\tag{6.1}
\]

The coefficients may vary arbitrarily between the two scales.  In
particular they may contain factorial weights, carry residues, and either
fixed or character-power phases.

> **Theorem 6.1 (two-scale three-state elimination).**  If (6.1) holds on
> the complete one-colour Taylor orbit at a coefficient-torus point, then
> the three functions \(M_j(t)\) are pairwise proportional.  Consequently
> the three states have equal \(B\)-level and, because they lie in one return
> fibre, equal \(R\)-level.  Hence any three-state shell with a nonzero
> colour-level transfer has unit scale-one-plus-scale-two orbit ideal after
> saturation by all coefficient channels.

### Proof

Put \(A_j=a_jM_j(t)\) and
\(\lambda_j=b_j/a_j^2\ne0\).  The first equation gives
\(A_3=-A_1-A_2\).  Since every channel is nonzero, \(A_1\ne0\) in the
rational function field.  With \(z=A_2/A_1\), the second equation becomes

\[
 (\lambda_1+\lambda_3)
 +2\lambda_3z
 +(\lambda_2+\lambda_3)z^2=0.
\tag{6.2}
\]

This is a nonzero polynomial over the constant field: its linear
coefficient is \(2\lambda_3\ne0\).  Thus the rational function \(z(t)\) is
algebraic over the constants and is itself constant.  The first equation
then makes \(A_3/A_1\) constant as well.

If the common \(B\)-selection count is \(q_B\), then

\[
 \deg_t M_j(t)=sq_B-\ell_B(M_j).
\tag{6.3}
\]

Proportional nonzero polynomials have equal degree, so all three
\(B\)-levels agree.  Their total levels agree by the return equation, hence
their \(R\)-levels agree too.  The Nullstellensatz gives the asserted
saturated unit ideal whenever two of the levels differ. \(\square\)

If some \(a_j\) or \(b_j\) is zero, the scale has already lost support and
reduces to a one- or two-state shell.  If all three levels agree, every pair
difference splits into an \(R\)-only and a \(B\)-only return relation.  This
is the zero-transfer product block, not a genuinely mixed primitive-transfer
obstruction.  It is the only three-state architecture not removed by the
abstract two-scale argument.

The exact Singular census gives a useful finite audit of that residual
block.  For span two and colour counts at most two, there are 564 signed
triples modulo common sign.  Twenty-four survive the scale-one orbit ideal,
and all die at scale two under both phase laws.  Increasing both colour
counts to three gives 11,988 signed triples and 416 scale-one survivors:
160 transfer level and 256 are zero-transfer.  Every character-power case
and 415 fixed-sign cases die at scale two.  The sole fixed-sign survivor is

\[
 F_N=X^N-\binom{2N}{N}Y^N-
       \frac{(3N)!}{(N!)^3}Z^N,
\quad
 \begin{cases}
 X=R_1^3B_1^2,\\
 Y=R_1^3B_0B_2,\\
 Z=R_0R_1R_2B_1^2.
 \end{cases}
\tag{6.4}
\]

The complete \(B\)-orbit and the first two rows force

\[
 \frac{Y}{X}=\frac14,
 \qquad
 \frac{Z}{X}=\frac1{12}.
\tag{6.5}
\]

Indeed a constant \(B_0(t)B_2/B_1(t)^2\) has value \(1/4\) from its
leading coefficients.  At scale three the normalized fixed-sign row is

\[
 1-20\left(\frac14\right)^3
   -1680\left(\frac1{12}\right)^3
 =-\frac{41}{144}\ne0.
\tag{6.6}
\]

Thus its three-scale ideal is the unit ideal.  At span three with both
counts at most two, all 240 scale-one survivors among 8,408 signed triples
die already at scale two.  These are bounded regressions; Theorem 6.1 is the
all-span statement.

## 7. Exact primitive census and the historical route frontier

The accompanying script constructs the Graver basis with Normaliz, retains
the exact scaled-factorial collisions, generates the finite translation
closures, and tests their coefficient-torus ideals over \(\mathbb Q\) with
Singular.  Its span-four full-envelope run gives

\[
\begin{array}{c|r}
\text{quantity}&\text{count}\\ \hline
\text{raw Graver elements}&426\\
\text{normalized mixed support-}\geq5\text{ packets}&65\\
\text{factorial-partition obstructions}&48\\
\text{exact scalar-factorial survivors}&17\\
\text{one-colour torus survivors}&0\\
\text{independent two-colour torus survivors}&0\\
\text{diagonal torus survivors}&2
\end{array}
\tag{7.1}
\]

The two diagonal survivors are (4.1).  The one-colour regression through
span five has 2,225 raw Graver elements, 404 normalized mixed candidates,
279 factorial obstructions, 125 exact factorial survivors, and no torus
survivor in either colour.  Theorem 3.1 proves the zero-survivor conclusion
for every span; the census is an independent exact regression and records
the harmless diagonal exception.

This changes the counterexample search.  A projected candidate can no longer
be one primitive semigroup relation with indistinguishable prime-power
signatures.  It must instead satisfy all of the following.

1. At least several inequivalent primitive orbit rows cancel linearly before
   any one of them is inherited.
2. The cancellation survives both Cartesian translation directions over the
   same high-digit quotient.
3. It is not a conformal sum of pure-colour beta/centered blocks and does not
   lie on coefficient-support loss.
4. Its scale-dependent carry coefficients satisfy the exact Hall factorial
   formulas, rather than arbitrary linear syzygy coefficients.
5. A fixed multiplier leaves a nonzero tail; otherwise it is terminal rather
   than a GVC counterexample.

The next finite object is therefore the **translation-observability
matroid** of one bounded affine correction shell.  Its columns are the full
mixed derivative orbits of the low-correction states, and its circuits are
minimal linear orbit syzygies.  Section 7 proves that distinct high affine
rays cannot keep cancelling around such a circuit at all scales.  The global
problem is to obtain one bounded correction shell uniformly from the
positive-density Hall face.  An empty torus saturation proves promotion for
that shell.  A proper saturation, together with a surviving shifted-power
row and a lift to one actual binary pair, is the first route that can produce
a genuine counterexample.

Separating primitive binomials is therefore necessary and is now proved, but
it is not by itself sufficient: toric Graver decomposition is multiplicative,
whereas the unresolved Hall cancellation is linear.

Theorem 6.1 moves the first finite-prefix obstruction from three to four
states, after zero-transfer product blocks are separated.  The all-scale
condition is much more rigid than any fixed conic/cubic test.  The following
argument is the relevant earlier theorem.

Let \(K\) be a characteristic-zero constant field, let \(t\) be one Taylor
parameter, and fix nonzero rational functions

\[
 f_j,h_j\in\overline K(t)^\times,qquad 1\leq j\leq r.
\tag{7.2}
\]

The affine-ray term \(h_jf_j^N\) includes every fixed exponent profile

\[
 \prod_\nu X_\nu(t)^{Na_{j\nu}+b_{j\nu}}
 =\left(\prod_\nu X_\nu(t)^{b_{j\nu}}\right)
  \left(\prod_\nu X_\nu(t)^{a_{j\nu}}\right)^N.
\tag{7.3}
\]

Partition the indices by \(i\sim j\) when \(f_i/f_j\in\overline K^\times\).
For a class \(C\), choose a representative \(f_C\), absorb the constant
ratios into the scalar coefficients, and put

\[
 H_C=\operatorname {span}_{\overline K}\{h_j:j\in C\}.
\tag{7.4}
\]

> **Theorem 7.1 (arbitrary-coefficient affine-ray splitting).**  There is an
> integer \(N_0\), depending only on the fixed functions in (7.2), such that
> for every \(N\geq N_0\) the subspaces
> \[
>   f_C^N H_C\subset\overline K(t)
> \tag{7.5}
> \]
> are in direct sum.  Consequently, for completely arbitrary scalar rows
> \(c_j(N)\in\overline K\), including zeros and nonexponential carry weights,
> \[
>   \sum_{j=1}^r c_j(N)h_j(t)f_j(t)^N=0
> \tag{7.6}
> \]
> splits, for every sufficiently large \(N\), into one identity inside each
> proportional-base class.  The active subset may change with \(N\) inside
> this fixed finite affine-ray universe.

### Proof

Choose a basis \(h_{C,1},\ldots,h_{C,d_C}\) of every \(H_C\), and form the
Wronskian of all functions \(h_{C,a}f_C^N\).  Put

\[
 u_C=\frac{f_C'}{f_C},\qquad
 E=\sum_{C<D}d_Cd_D.
\tag{7.7}
\]

After factoring \(f_C^N\) from each of the \(d_C\) columns, the normalized
Wronskian is a polynomial in \(N\) over \(\overline K(t)\).  Its coefficient
of \(N^E\), up to one nonzero sign, is the confluent-Vandermonde product

\[
 \prod_C W(h_{C,1},\ldots,h_{C,d_C})
 \prod_{C<D}(u_D-u_C)^{d_Cd_D}.
\tag{7.8}
\]

For completeness, write

\[
 f_C^{-N}D^k(hf_C^N)=(D+Nu_C)^kh.
\tag{7.9}
\]

In the determinant, the top nonzero \(N\)-degree uses derivative orders
\(0,\ldots,d_C-1\) inside the \(C\)-block.  Their block determinants are the
Wronskians of the \(h_{C,a}\); between two blocks the ordinary Vandermonde
factor \(u_D-u_C\) occurs once for every pair of columns.  Terms in which a
derivative hits some \(u_C\) have smaller \(N\)-degree.  This gives (7.8).

Each first product factor is nonzero because linearly independent rational
functions have nonzero Wronskian in characteristic zero.  Each second product
factor is nonzero because

\[
 u_C=u_D\quad\Longleftrightarrow\quad(f_C/f_D)'=0
 \quad\Longleftrightarrow\quad f_C/f_D\in\overline K^\times.
\tag{7.10}
\]

Thus the normalized Wronskian is a nonzero polynomial in \(N\).  It vanishes
at only finitely many constant values of \(N\).  Outside those values the
displayed functions are linearly independent over the constant field, which
is precisely the direct sum (7.5).  Equation (7.6) then splits by classes.
\(\square\)

The proof is unaffected by a finite carry automaton or a congruence condition
on \(N\): restrict to one state or residue class and replace \(N=Lk+r\).
Then \(f_j^N=f_j^r(f_j^L)^k\), which has the same form.  Several Cartesian
translation variables can be restricted to a degree-separating
one-parameter curve.  To see this, choose finitely many regular points whose
evaluation matrices witness the required linear independence and whose
values witness every nonconstant ratio, then interpolate one polynomial
curve through those points.  Restriction to that curve preserves all the
needed nonvanishing statements simultaneously.

In the uncorrected projected model \(h_j=1\).  Theorem 7.1 says that every
all-scale row with a fixed finite state universe splits by proportionality of
the complete orbit functions, no matter how its factorial, carry, unit, or
torsion coefficient varies with \(N\).  Proportional Taylor polynomials have
equal degree, so every resulting class has zero marked-side level transfer.
The return equation gives zero transfer on the other side as well.  This
closes the cross-ray part of the formerly stated “nonexponential
affine-carry coefficient” gap.  For a general affine profile, however,
several distinct low corrections \(h_j\) can lie over the same \(f_C\).
Theorem 7.1 exposes their common high quotient but does not split the residual
identity inside \(H_C\).  That residual is the finite
translation-observability matroid on which the actual carry rows,
adjacent-channel identities, and marked-side data must still be imposed.

For fixed signed factorial rows there is a sharper description.  Let

\[
 F_N(t)=\sum_{j=1}^r\epsilon_jW_{\alpha_j}(N)M_j(t)^N,
 \qquad \epsilon_j\in\{1,-1\},
\tag{7.11}
\]

where the states and signs do not vary with \(N\), and
\(W_{\alpha_j}\) is the exact scalar factorial ray.

> **Theorem 7.2 (fixed-template all-scale pairing).**  If
> \(F_N(t)=0\) for every \(N\geq1\), then the terms split by equality of
> their complete factorial vectors, and within each vector class the
> positive and negative multisets of orbit functions \(M_j(t)\) are equal.
> In particular a four-state fixed-sign all-scale shell is a sum of at most
> two pairwise orbit identities; it cannot be a genuinely inequivalent
> four-state linear syzygy.  If instead the only coefficients are character
> powers \(\chi_j^N\), with no fixed signed multiplicities, no nonzero
> coefficient-torus shell can vanish at every scale.

### Proof

Factorial trace independence (FTI1), applied over the rational function field
of the translation parameters, separates the distinct canonical factorial
vectors: the coefficients \(\epsilon_jM_j(t)^N\) are exponential-rational
in \(N\).  Hence in each same-vector class \(C\),

\[
 p_N=\sum_{j\in C}\epsilon_jM_j(t)^N=0
 \qquad(N\geq1).
\tag{7.12}
\]

The rational generating function

\[
 \sum_{N\geq1}p_Nz^N
 =\sum_{j\in C}\epsilon_j
   \frac{M_j(t)z}{1-M_j(t)z}
\tag{7.13}
\]

vanishes.  Partial fractions show that the total signed multiplicity at each
nonzero value of \(M_j(t)\) is zero.  This is exactly the asserted multiset
pairing.  For character powers replace \(M_j\) by \(\chi_jM_j\); every
partial-fraction multiplicity is then positive, which is impossible on the
coefficient torus. \(\square\)

Theorem 7.2 uses factorial trace independence to identify the exact signed
multiplicities.  Theorem 7.1 is weaker inside one proportional-base class but
strictly broader across classes: it allows arbitrary scalar \(N\)-dependence
and changing support.  Inside one class it leaves the fixed correction space
\(H_C\); before common-quotient promotion the affine-ray universe may also
grow without bound with \(N\).

There is no conflict with the isoperiodic counterexample
\(\operatorname {CT}(z+z^{-1})^N=
\operatorname {CT}(z^2+z^{-2})^N\) in the translation-tangent note.
After coefficient extraction, each side is a sum over a number of Laurent
monomials growing with \(N\); it is not a fixed finite sum of the form (7.6).
That example is a model of exactly the unbounded-support escape left here.

The exact four-state audits illustrate the distinction.  At span two with
both counts at most two, 40 of 928 signed quartets survive scale one and all
die at scale two.  Increasing both counts to three gives 52,416 signed
quartets and 928 scale-one survivors.  Character-power phases kill all 928
at scale two.  Fixed external signs leave 60 through scale three.  All 60
are exact all-scale identities, but every one has the pairing in Theorem
7.2: two opposite-sign pairs have equal factorial rays and equal complete
one-colour orbit monomials.  Only 12 pair loci survive translation of the
other colour, and only six survive the two independent Taylor towers.  For
all six, exact ideal comparison gives

\[
 \left(R_1^2-4R_0R_2,\ B_1^2-4B_0B_2\right)
\tag{7.14}
\]

after coefficient-torus saturation.  This is the product of the two
quadratic Veronese loci, hence the already-safe homogeneous/pure-power block.
The 60 collisions are therefore useful obstructions to a naive claim that
every fixed-sign finite-prefix ideal is a unit ideal, but none is an
inequivalent all-scale Hall circuit or a GVC counterexample.

## 8. Four states with arbitrary coefficients need only three scales

There is also a sharp finite-prefix result which is independent of Theorem
7.1.  Let \(M_1,\ldots,M_4\) be nonzero rational orbit functions over a
characteristic-zero constant field and suppose

\[
 \sum_{i=1}^4a_iM_i=0,\qquad
 \sum_{i=1}^4b_iM_i^2=0,\qquad
 \sum_{i=1}^4c_iM_i^3=0,
 \qquad a_ib_ic_i\ne0.
\tag{8.1}
\]

Put

\[
 A_i=a_iM_i,\qquad
 \lambda_i=\frac{b_i}{a_i^2},\qquad
 \mu_i=\frac{c_i}{a_i^3}.
\tag{8.2}
\]

> **Theorem 8.1 (arbitrary-coefficient four-state pairing).**  Under (8.1),
> either every projective ratio \(M_i/M_j\) is constant, or, after relabelling,
> the three rows split into two pair identities:
> \[
> \begin{aligned}
> A_2&=-A_1,&\lambda_2&=-\lambda_1,&\mu_2&=\mu_1,\\
> A_4&=-A_3,&\lambda_4&=-\lambda_3,&\mu_4&=\mu_3.
> \end{aligned}
> \tag{8.3}
> \]
> Consequently every positive-dimensional four-state ratio component is
> already a sum of two proportional orbit blocks.  In a return fibre each
> block has zero level transfer in both colours.  Thus no genuinely
> inequivalent transferring four-state shell survives scales one, two, and
> three, even when all twelve row coefficients are unrelated.

### Proof

Work over the algebraic closure of the constant field.  Eliminate
\(A_4=-x-y-z\), where \(x=A_1,y=A_2,z=A_3\).  The last two rows become

\[
\begin{aligned}
 Q={}&\lambda_1x^2+\lambda_2y^2+\lambda_3z^2
       +\lambda_4(x+y+z)^2,\\
 C={}&\mu_1x^3+\mu_2y^3+\mu_3z^3
       -\mu_4(x+y+z)^3.
\end{aligned}
\tag{8.4}
\]

If the ratio image is not a point, its closure is a curve in
\(\mathbb P^2\) contained in \(Q=C=0\), so the conic and cubic share a curve
component.  The restriction of the nonsingular diagonal quadratic form
\(\sum_i\lambda_iA_i^2\) to \(\sum_iA_i=0\) has radical dimension at most
one.  Hence \(Q\) has rank three or two.

First suppose \(Q\) has rank three.  It is irreducible, so \(Q\mid C\).
Write

\[
 C=Q(\alpha x+\beta y+\gamma z),\qquad
 L_i=\lambda_i+\lambda_4\quad(1\leq i\leq3).
\tag{8.5}
\]

Comparing the six mixed-square coefficients in pairs and then the \(xyz\)
coefficient gives

\[
 L_1(\beta-\gamma)=
 L_2(\alpha-\gamma)=
 L_3(\alpha-\beta)=0,\qquad
 \lambda_4(\alpha+\beta+\gamma)=-3\mu_4.
\tag{8.6}
\]

These equations have no rank-three solution.  If none of the \(L_i\)
vanishes, then \(\alpha=\beta=\gamma\); comparing any mixed square with
the \(xyz\) equation forces the corresponding \(\lambda_i=0\).  If exactly
one \(L_i\) vanishes, the three linear coefficients are again equal, while
the mixed square using that \(L_i\) gives \(2\lambda_4\alpha=-3\mu_4\) and
the \(xyz\) equation gives \(\lambda_4\alpha=-\mu_4\).  If exactly two
vanish, say \(L_1=L_2=0\), the same coefficients give

\[
 (\lambda_1,\lambda_2,\lambda_3,\lambda_4)
 =(-d,-d,d,d),\qquad d=\lambda_4,
\tag{8.7}
\]

whose restricted determinant is zero.  If all three vanish, the mixed
squares give
\(2\lambda_4\alpha=2\lambda_4\beta=2\lambda_4\gamma=-3\mu_4\),
contradicting the \(xyz\) equation.  Thus the smooth-conic branch is empty.

It remains that \(Q\) has rank two.  Its determinant is

\[
 \left(\lambda_1\lambda_2\lambda_3\lambda_4\right)
 \left(\sum_{i=1}^4\lambda_i^{-1}\right)=0,
\tag{8.8}
\]

or, in expanded form,
\(\lambda_1\lambda_2\lambda_3+
\lambda_4(\lambda_1\lambda_2+\lambda_1\lambda_3+
\lambda_2\lambda_3)=0\).
Put \(r_i=\lambda_i^{-1}\).  Then \(\sum_i r_i=0\), and \(r\) is the
radical point of the conic on the hyperplane \(\sum_iA_i=0\).  The conic is
two distinct lines through \(r\).  Parametrize the common line of \(Q\) and
\(C\) by

\[
 A_i(s)=r_i+s u_i=r_i(1+t_is),\qquad t_i=\lambda_i u_i.
\tag{8.9}
\]

The cubic identity on that line is

\[
 \sum_{i=1}^4\mu_i r_i^3(1+t_is)^3=0.
\tag{8.10}
\]

For distinct values \(\tau\), the at most four polynomials
\((1+\tau s)^3\) are linearly independent: their coefficient matrix is a
Vandermonde matrix with nonzero binomial row factors.  Therefore the nonzero
weights \(\mu_ir_i^3\) must sum to zero separately in every equal-\(t_i\)
class.  No class is a singleton.  The \(t_i\) are not all equal, since then
\(u\) would be proportional to \(r\).  Four indices must consequently split
as two pairs, say

\[
 t_1=t_2=a,\qquad t_3=t_4=b,\qquad a\ne b.
\tag{8.11}
\]

Both \(r\) and \(u\) lie in \(\sum_iA_i=0\).  Hence

\[
 (r_1+r_2)+(r_3+r_4)=0,\qquad
 a(r_1+r_2)+b(r_3+r_4)=0.
\tag{8.12}
\]

Since \(a\ne b\), each pair sum is zero.  Thus
\(\lambda_2=-\lambda_1\) and \(\lambda_4=-\lambda_3\).  The two grouped
coefficients in (8.10) then give
\(\mu_2=\mu_1\) and \(\mu_4=\mu_3\).  Along the line,
\(A_2=-A_1\) and \(A_4=-A_3\), which is (8.3).  Conversely, (8.3) plainly
makes all three rows vanish pairwise.  This proves the theorem.
\(\square\)

The determinant formula in (8.8) is deliberately written in both compact
and expanded form: the compact expression is
\[
 (\lambda_1\lambda_2\lambda_3\lambda_4)
 \left(\sum_i\lambda_i^{-1}\right).
\]
A zero row coefficient is support loss and reduces to the three-state
theorem.  Thus Theorem 8.1 closes the entire arbitrary-coefficient
four-state finite-prefix branch, rather than merely the bounded fixed-sign
census.

The power form in (8.1) is essential.  A genuinely affine-offset term
\(h_i(t)f_i(t)^N\), with several low corrections over one common \(f_i\),
reduces instead to a linear identity among the \(h_i\) and belongs to the
residual common-base circuit described after Theorem 7.1.

## 9. The all-scale circuit ideal

The next computation has an exact theorem/counterexample dichotomy.  Fix one
stabilized affine correction alphabet and one common high-digit quotient.
For each scale \(N\), let

\[
 F_N(R,B)=\sum_{\sigma\in\mathcal C_N}
 \epsilon_\sigma(N)W_\sigma(N)R^{u_\sigma(N)}B^{v_\sigma(N)}
\tag{9.1}
\]

be the actual Hall shell: \(W_\sigma\) is the exact factorial weight and
\(\epsilon_\sigma\) is the carry/character coefficient, not a free
parameter.  With the two commuting Cartesian derivations, define

\[
 I_K=
 \left(
 D_1^aD_2^bF_N:
 1\leq N\leq K,\ a,b\geq0
 \right):\left(\prod_\nu R_\nu B_\nu\right)^\infty.
\tag{9.2}
\]

Each derivative list is finite.  The chain \(I_1\subseteq I_2\subseteq
\cdots\) lies in one Noetherian Laurent-coordinate ring, so the ideal
generated by all scales has a finite certificate.

There is an exact elimination once Theorem 7.1 has exposed one common base.
Let \(A\) be the coefficient Laurent ring, let
\(\mathcal D(J)\) denote the ideal generated by all iterates of the two
commuting Cartesian derivations on a set \(J\), and suppose

\[
 F_N=f^NH_N,\qquad
 H_N=\sum_{j=1}^r c_j(N)h_j,\qquad f\in A^\times.
\tag{9.3}
\]

The coefficients \(c_j(N)\) are arbitrary constants.  Put

\[
 C=\operatorname {span}_K
 \{c(N)=(c_1(N),\ldots,c_r(N)):N\geq1\}\subseteq K^r,
\quad
 T(c)=\sum_jc_jh_j.
\tag{9.4}
\]

> **Theorem 9.1 (common-base all-scale ideal collapse).**  In the Laurent
> coefficient ring,
> \[
>  \mathcal D\bigl(\{F_N:N\geq1\}\bigr)
>  =
>  \mathcal D\bigl(T(C)\bigr).
> \tag{9.5}
> \]
> In particular at most \(\dim C\leq r\) scale rows generate the complete
> all-scale Cartesian ideal.  If \(C=K^r\) and every \(h_j\) is a Laurent
> monomial, the ideal is the unit ideal on the coefficient torus.

### Proof

Fix \(N\).  Since \(f\) is a unit, \(H_N=f^{-N}F_N\) belongs to the ideal
generated by \(F_N\).  Conversely \(F_N=f^NH_N\).  The multivariate Leibniz
rule expresses every \(D_1^aD_2^bF_N\) as \(f^N D_1^aD_2^bH_N\) plus terms
in lower total derivatives of \(H_N\).  Induction on \(a+b\), multiplying by
\(f^{-N}\), gives

\[
 \mathcal D(F_N)=\mathcal D(H_N).
\tag{9.6}
\]

Taking the sum over \(N\) and using the constant-field linearity of the
derivations replaces the rows \(H_N=T(c(N))\) by any basis of their
coefficient span \(C\), proving (9.5).  If \(C=K^r\), every \(h_j=T(e_j)\)
lies in the ideal; a Laurent monomial is a unit.  \(\square\)

> **Corollary 9.2 (projective rigidity on a primitive orbit circuit).**
> Evaluate the correction functions \(h_j\) on a coefficient-torus
> Cartesian orbit.  If they form a matroid circuit, their relation kernel is
> one-dimensional.  A common-base shell can vanish at every scale only if
> all nonzero carry/factorial coefficient vectors \(c(N)\) are proportional
> to that one circuit vector.  Two nonproportional rows exclude the
> full-rank circuit stratum; any surviving coefficient-torus point lies on a
> further orbit-rank drop.

Indeed (9.5) puts \(C\) inside the relation kernel.  A circuit kernel has
dimension one.  Equivalently, after also saturating by one nonzero circuit
minor, two nonproportional rows give the unit ideal.  This is the exact role
for prime-power tomography in a bounded packet: it need only disprove
projective constancy of the actual coefficient vectors.  It need not compare
all packet sums or guess an infinite stabilization cutoff.  The remaining
rank-drop locus must be sent to the already-safe Veronese/beta blocks or
analyzed as a smaller orbit matroid.

The discriminant identity (5.3) is the calibration.  After a common-base
factor is inserted, every scale has coefficient vector proportional to
\((1,-4,-1)\), so \(C\) has rank one.  Its orbit-rank locus is the quadratic
Veronese block, already safe.  A genuine bounded obstruction must reproduce
that projective constancy with the exact Hall carry weights while landing on
a rank-drop component not covered by the terminal blocks and retaining a
shifted-power tail.

There is no longer an unclassified factorial mechanism inside that
projective-constancy test.  For a fixed affine correction state, write its
eventual scalar coefficient as

\[
 c_j(N)=\kappa_j\eta_j(N)
 \prod_{a,b}((aN+b)!)^{e_{j,a,b}},
 \qquad \kappa_j\ne0,\quad \eta_j(N)\in S_j\subset K^\times,
\tag{9.7}
\]

with finitely many integer offsets and a finite nonzero output alphabet
\(S_j\).  This includes signs, torsion-character values, and the output of a
fixed finite carry automaton.  The exact successor divisor is

\[
 \Delta_j=\sum_{a,b}e_{j,a,b}\sum_{u=1}^{a}
 \left[\frac{b+u}{a}\right]\in\mathbb Z[\mathbb Q].
\tag{9.8}
\]

Theorem 5.2 of
[FACTORIAL_TRACE_INDEPENDENCE.md](FACTORIAL_TRACE_INDEPENDENCE.md) proves
that equality of two such exact divisors is much stronger than equality of
their gamma-orbit signatures: their factorial ratio is already constant,
not a nontrivial exponential times a rational function.  Moreover every
such equality is generated by

\[
 (aN+k)!(cN+l-1)!
 =\frac ac\,(aN+k-1)!(cN+l)!,
 \qquad k/a=l/c.
\tag{9.9}
\]

> **Corollary 9.3 (affine-factorial projective test).**  Suppose every
> nonzero entry of a common-base coefficient vector has the form (9.7).
> The rows \(c(N)\) have rank one for all sufficiently large \(N\) only if
> all active \(\Delta_j\) agree, every relative finite-state law
> \(\eta_i(N)/\eta_j(N)\) is constant, and every pair of factorial profiles
> differs by a finite sum of the boundary transfers (9.9).  Conversely these
> conditions, with the resulting constants absorbed into \(\kappa_j\), give
> a projectively constant row.

Indeed, rank one makes each ratio \(c_i(N)/c_j(N)\) constant.  Its
factorial part is therefore a constant times the finite-valued sequence
\(\eta_j(N)/\eta_i(N)\).  Corollary 5.3 of the factorial-trace note says that
an integer-affine factorial ratio with finite image is constant.  Theorem
5.2 then gives equality of the exact divisors and the boundary presentation
(9.9), after which the relative finite-state law is constant.  The converse
is immediate.

The statement has a marked version: attach the operator, polynomial,
radial, or individual channel label to every point in (9.8).  The proof is a
direct sum over labels, so boundary transfer is allowed only between factors
with the same retained label.  Thus operator-side partitions,
polynomial-side partitions, and low digits see every transfer which moves a
boundary factor between different marked channels.  After those data are
retained, a bounded projectively constant obstruction cannot be blamed on a
mysterious affine-factorial identity.  It must come from an allowed
same-label redistribution together with a genuinely nongeometric carry-state
law having infinite output or additive packet sums, or from the further
orbit-rank locus in Corollary 9.2.

Additive packet coefficients require one more distinction.  Boundary
transfers do not classify addition: for example,
\[
 \frac{N!}{(N-1)!}+\frac{(N+1)!}{N!}
 -\frac{(2N+1)!}{(2N)!}=N+(N+1)-(2N+1)=0.
\tag{9.10}
\]
This is not a new obstruction when the carry law is periodic; it reduces to
finite rational-function elimination.

Fix a period \(L\).  Suppose that on every residue \(N=Lk+r\), each
coefficient is a fixed finite sum
\[
 c_j(Lk+r)=\sum_{\alpha}
 \gamma_{j,r,\alpha}
 \prod_{a,b}((aLk+ar+b)!)^{e_{j,r,\alpha,a,b}},
 \qquad \gamma_{j,r,\alpha}\in K.
\tag{9.11}
\]
For each summand, collect its signed slope vector \(v\) after the
substitution and use integer factorial shifts to write it uniquely as
\[
 R_{j,r,\alpha}(k)\Phi_v(k),
 \qquad R_{j,r,\alpha}(k)\in K(k).
\]
After summing equal vectors, write the coefficient of \(\Phi_v\) as
\(R_{j,r,v}(k)\).

> **Theorem 9.4 (periodic additive coefficient elimination).**  Let
> \(u=(u_1,\ldots,u_s)\in(K^\times)^s\).  The rows (9.11) lie on the fixed
> projective line \(Ku\) for every sufficiently large \(N\) if and only if,
> for every residue \(r\), vector \(v\), and pair \(i,j\),
> \[
>  u_iR_{j,r,v}(k)-u_jR_{i,r,v}(k)=0
>  \quad\text{in }K(k).
> \tag{9.12}
> \]
> Thus the all-scale projective test for every eventually periodic finite
> additive Hall alphabet is a finite list of exact polynomial identities
> after clearing denominators.

For necessity, projective constancy gives
\(u_ic_j(Lk+r)-u_jc_i(Lk+r)=0\).  Theorem 2.1 of the factorial-trace note
separates its distinct canonical signed slope vectors over rational
coefficients, giving (9.12).  Conversely (9.12) makes
\(c_j(Lk+r)/u_j\) independent of \(j\) on each residue, so every row lies on
\(Ku\).  There are only finitely many residues, vectors, and rational
identities.  \(\square\)

The same proof handles a finite union of periodic states after taking a
common period.  It covers all fixed torsion-character laws and every carry
output already known to depend eventually only on \(N\bmod L\).  It does
not cover a genuinely nonperiodic automatic carry language, nor a packet
whose affine profile alphabet grows with \(N\).  Those are now the two ways
scale dependence can evade exact finite rational elimination.

There is a precise two-prime interface.  Let \(\theta(N)\) be the finite
symbol which records the complete additive formula used at scale \(N\):
active summands, affine profiles, and finite coefficients.  Recall
[Cobham's theorem](https://doi.org/10.1017/CBO9780511546563.012): a
finite-alphabet sequence which is automatic in two multiplicatively
independent bases is ultimately periodic.

> **Corollary 9.5 (two-base automatic collapse).**  If one
> prime-independent coefficient-state sequence \(\theta(N)\) is both
> \(p\)-automatic and \(q\)-automatic for multiplicatively independent
> integers \(p,q\ge2\), then its bounded common-base projective test is the
> finite rational elimination of Theorem 9.4.

Cobham makes \(\theta\) ultimately periodic, so Theorem 9.4 applies after
discarding its finite preperiod.  \(\square\)

The hypothesis is the exact missing compatibility, not a consequence of
one-prime tomography.  The Hall construction currently supplies
prime-dependent carry automata and often probes only sparse exponents
\(N=qp^e\).  To invoke Corollary 9.5 one must prove that two such automata
encode the same intrinsic full-scale state sequence, or at least a common
automatic refinement.  A pair of unrelated automata which merely vanish on
their own prime-power subsequences is outside Cobham's theorem.  Conversely,
failure of this compatibility gives a concrete obstruction to extract: two
prime-dependent minimal state sequences with no common finite refinement.

The literal complete carry-position list is not even a finite-alphabet
output: its length and its marked positions grow with \(N\).  The following
argument does not rely on that type mismatch.  It passes to a one-bit
finite-alphabet quotient and proves that this quotient already forbids the
desired common refinement.

That obstruction already occurs in the smallest Kummer profile.  For an odd
prime \(p\), put
\[
 \chi_p(N)=
 \mathbf 1_{\{v_p\binom{2N}{N}>0\}}.
\tag{9.13}
\]
Equivalently, \(\chi_p(N)\) records whether adding \(N+N\) in base \(p\)
has at least one carry.  This is a quotient of the complete carry-position
record for the two-part partition \((N,N)\).

> **Theorem 9.6 (Cobham carry-refinement obstruction).**  For every odd
> prime \(p\):
>
> 1. \(\chi_p\) is a two-state \(p\)-automatic sequence but is not
>    ultimately periodic;
> 2. for every \(q,e\geq 0\),
>    \[
>     \chi_p(qp^e)=\chi_p(q);
>    \tag{9.14}
>    \]
> 3. if \(q\geq2\) is multiplicatively independent of \(p\), there is no
>    finite-alphabet sequence \(\theta\), automatic in both bases \(p\)
>    and \(q\), from which \(\chi_p\) is obtained by a letter map.
>
> In particular, no common finite two-base automatic state can retain even
> the “some carry occurred” quotient of a complete prime-specific Kummer
> record, and therefore none can refine the complete records at two
> unrelated primes.

Write \(h=(p-1)/2\).  There is no carry in \(N+N\) exactly when every
base-\(p\) digit of \(N\) is at most \(h\).  Hence a two-state automaton
starts in “no high digit seen” and changes permanently to “high digit seen”
on a digit \(>h\).  Its output is \(\chi_p\).

To prove nonperiodicity, let \(T>0\) be any proposed eventual period and
write \(T=p^a u\) with \(p\nmid u\).  Choose
\(m\in\{1,\ldots,p-1\}\) such that
\[
 mu\equiv-1\pmod p.
\tag{9.15}
\]
The digit in position \(a\) of \(mT\) is then \(p-1\).  Choose \(k\) so
large that \(p^k\) is beyond the proposed preperiod and \(p^k>mT\).  The
number \(p^k\) has only the digit \(1\leq h\), so
\(\chi_p(p^k)=0\).  The base-\(p\) expansion of \(p^k+mT\) is the disjoint
high digit \(1\) followed by the digits of \(mT\); it contains the digit
\(p-1\), so \(\chi_p(p^k+mT)=1\).  These two indices differ by a multiple
of \(T\), contradicting eventual periodicity.  Equation (9.14) is just a
digit shift.  Finally, if \(\theta\) existed, Cobham would make it
ultimately periodic, and its letter image \(\chi_p\) would be ultimately
periodic as well.  This contradiction proves the theorem.  \(\square\)

Theorem 9.6 closes the proposed Cobham promotion in the negative.  It is
not merely a failure to match labels chosen by two implementations: even
the one-bit quotient “some Kummer carry occurred” cannot be retained by a
common two-base automatic state.  A prime-independent quotient may still
satisfy Corollary 9.5, but it must forget this genuine separating datum.
Moreover (9.14) explains why tests restricted to \(qp^e\) do not detect the
problem: each such ray is stationary although the full-scale sequence is
nonperiodic.

This is an obstruction to the proof route, not a binary GVC counterexample.
It supplies neither an all-scale cancellation nor a nonzero shifted-power
tail.  Cobham can eliminate a common prime-independent quotient after that
quotient has been proved by other means; it cannot construct a common Hall
state from the prime-specific carry automata.

Passing to a synchronized multidimensional carry language does not preserve
the hoped-for escape.  Automatic relations are closed under coordinate
projection and finite output maps.  Thus any relation recognizable in both
bases from which “there exists a carry position” can be recovered would make
\(\chi_p\) automatic in both bases, contradicting Theorem 9.6.  An
infinite-output weighted construction which does not admit this finite
projection is outside the theorem, but it has also discarded the proposed
carry separator and requires a separate regular-sequence argument.

- If some \(I_K=(1)\), the displayed Nullstellensatz identity is a rigorous
  finite obstruction to that affine Hall template.
- If the all-scale ideal is proper, an algebraic coefficient-torus point
  survives every pure translation row.  One then computes the shifted-power
  or fixed-multiplier row.  A proved nonzero infinite tail, together with a
  lift from the packet coordinates to one actual binary \((\Lambda,P)\), is
  a genuine GVC(2) counterexample candidate.  The Hall-envelope theorem shows
  that no lift satisfying all of these conditions can exist.

Before a common base has been exposed, equality of two successive computed
ideals is not an all-scale certificate: a later row or a new affine profile
may enlarge the ideal.  A hypergeometric recurrence or carry automaton is
still needed to control the moving profile family.  After common-base
exposure, Theorem 9.1 replaces that infinite calculation by the exact
coefficient-span rank.

Theorems 7.1 and 8.1 change this target.  An all-scale relation supported on
any fixed finite affine-ray universe already splits by proportional high
orbit rays for arbitrary carry coefficients, and the uncorrected four-state
power branch splits after three rows.  The all-scale ideal (9.2) therefore
needs to be computed only on a **common-base low-correction circuit**.

There are now two precise unresolved inputs.

1. **Bounded isobase circuit.**  Theorem 9.1 reduces this case to the finite
   coefficient span \(C\), and Corollary 9.3 classifies its affine-factorial
   rank-one part by marked boundary transfers, while Theorem 9.4 eliminates
   every eventually periodic finite additive law.  Theorem 9.6 rules out a
   generic common automatic refinement which retains the actual carry state.
   Each primitive orbit circuit must therefore be tested directly with its
   prime-specific full-scale Hall coefficients, or be reduced to a proved
   prime-independent quotient which forgets no coefficient-rank separator.
   Infinite-output rows require a separate regular-sequence argument.
   Nonproportional rows force orbit-rank drop; a projectively constant
   boundary-transfer survivor must be classified as an already-safe packet
   or tested for a shifted-power tail.
2. **Unbounded moving face.**  On an unexposed positive-density Hall face,
   the number of affine exponent profiles grows with \(N\).  One must extract
   a fixed finite quotient class before Theorem 7.1 can be applied.

In the second case, normalize the exponent profiles in the compact return
polytope and take an accumulation face.  Prime-power tomography should be
applied to that face, not to individual fixed packets: its task is to force a
bounded affine correction alphabet over one common high quotient.  Failure
supplies a concrete obstruction—a sequence of minimal orbit syzygies whose
support cardinality tends to infinity and whose normalized profiles remain
positive-dimensional.

Thus the infinite all-scale ideal problem is closed for every bounded
common-base circuit: only a finite coefficient-rank calculation remains.
Its factorial subcalculation is the exact boundary signature (9.8).  The
periodic additive subcalculation is the rational identity list (9.12).  The
Cobham branch is now exact but conditional, and Theorem 9.6 proves that the
condition cannot generically be obtained while retaining Kummer carries.
The theorem left in this route must directly prove prime-specific nonperiodic
carry-state projective nonconstancy or safety for every actual primitive
correction circuit, and promote any growing positive-density face to such a
bounded circuit.  A projectively constant unsafe nonperiodic
boundary-transfer circuit or an unbounded sequence, plus a nonzero
shifted-power tail and a lift to one actual \((\Lambda,P)\), would have been
the remaining route to a counterexample.  GVC2ENV rules out an actual binary
counterexample independently of that classification.

## 10. Reproduction

Run

```bash
python3 scripts/research_binary_gvc_translation_observability.py \
  --radial-degree 4 \
  --output artifacts/generated-results/binary_gvc_translation_observability_span4.json

python3 scripts/research_binary_gvc_translation_observability.py \
  --radial-degree 5 \
  --modes operator,polynomial \
  --output artifacts/generated-results/binary_gvc_translation_observability_span5_one_colour.json

.venv/bin/python scripts/research_binary_gvc_all_scale_orbit_circuits.py \
  --structural-certificate --maximum-wronskian-rank 7 \
  --maximum-affine-slope 5 --maximum-affine-offset 3 \
  --output artifacts/generated-results/binary_gvc_all_scale_structural_certificate.json

python3 scripts/verify_binary_gvc_cobham_carry_obstruction.py \
  --primes 3,5,7,11,13 --maximum-period 512 \
  --output artifacts/generated-results/binary_gvc_cobham_carry_obstruction.json

python3 scripts/research_binary_gvc_all_scale_orbit_circuits.py \
  --span 2 --state-count 3 \
  --maximum-operator-count 3 --maximum-polynomial-count 3 \
  --maximum-scale 3 \
  --output artifacts/generated-results/binary_gvc_all_scale_orbit_circuits_span2_counts3.json

python3 scripts/research_binary_gvc_all_scale_orbit_circuits.py \
  --span 3 --state-count 3 \
  --maximum-operator-count 2 --maximum-polynomial-count 2 \
  --maximum-scale 3 \
  --output artifacts/generated-results/binary_gvc_all_scale_orbit_circuits_span3_counts2.json

python3 scripts/research_binary_gvc_all_scale_orbit_circuits.py \
  --span 2 --state-count 4 \
  --maximum-operator-count 2 --maximum-polynomial-count 2 \
  --maximum-scale 4 \
  --output artifacts/generated-results/binary_gvc_all_scale_orbit_circuits_span2_four_state_counts2.json

python3 scripts/research_binary_gvc_all_scale_orbit_circuits.py \
  --span 2 --state-count 4 \
  --maximum-operator-count 3 --maximum-polynomial-count 3 \
  --maximum-scale 3 --require-factorial-pairing \
  --output artifacts/generated-results/binary_gvc_all_scale_factorial_pair_circuits_span2_counts3.json
```

The structural command checks all 44 integer block partitions through
Wronskian rank seven, exact full-rank and projectively constant common-base
coefficient spans, the saturated smooth-conic exclusion, and the pair-block
converse.  It also computes all 106 Graver moves for slopes at most five and
offsets from -3 through 3 and certifies every one by the boundary
presentation of Theorem 5.2 in
[FACTORIAL_TRACE_INDEPENDENCE.md](FACTORIAL_TRACE_INDEPENDENCE.md).  An exact
rational-function example replays Theorem 9.4 and detects a perturbed
projective row.  The rank and affine-universe bounds are regressions for the
unbounded proofs, not bounded extrapolations.  Every output explicitly
labels a proper orbit ideal as a projected promotion obstruction, not a Hall
packet or GVC counterexample.

The Cobham command independently checks Kummer's carry count against
Legendre valuations through \(N=10{,}000\), replays the two-state digit
automata for \(p=3,5,7,11,13\), constructs exact contradictions for every
proposed period through \(512\) beyond index \(10{,}000\), and verifies the
stationary identity (9.14).  The bounded loops audit the formulas; the
arbitrary-\(T\) construction in Theorem 9.6 is the proof.
