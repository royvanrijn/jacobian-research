# Dickson Galois realizations with controlled local actions

The small-trace construction gives a quantitative family of locally cyclic
Galois extensions with an exact discriminant, as well as its ordinary class-group
quotient. Over the fixed real cyclotomic layer these are tame dihedral
extensions, pairwise linearly disjoint when their quadratic mirrors differ.
The same formulas identify a limitation for root-discriminant optimization:
the cyclotomic contribution grows with the level.
An explicit second coefficient choice realizes every nonabelian nonsplit
local subgroup specified in IG3, with the same degree-\(m+1\) discriminant
form and the same counting exponents.

The deductions below use [SD1–SD2](CYCLOTOMIC_MIRROR_TORSION_FROM_SMALL_DICKSON_TRACE.md).
They do not supply a new proof of those inputs or a historical-priority claim.

## The exact Galois counting consequence

Fix \(m=3^n\), \(n\ge2\), and the scale and local slice of SD1, including
its fixed denominator noncube test. Write
\[
 D=\operatorname{disc}K,\quad D_0=D/3,\quad
 M=\mathbb Q(\sqrt{-D_0}),\quad R=\mathbb Q(\zeta_m)^+,
 \quad B=MR.
\]
Here \(-D_0\) is the fundamental discriminant of \(M\), and
\([R:\mathbb Q]=m/3\). Let \(W\) be the normal closure of the degree-\(m\)
Dickson root field in SD2. A locally cyclic extension means a Galois extension
whose decomposition group at every place is cyclic.

**Theorem IG1.** The extension \(W/\mathbb Q\) is locally cyclic, has group
\(G_m=C_m\rtimes(\mathbb Z/m)^\times\), and satisfies
\[
 \begin{split}
 [W:\mathbb Q]&=2m^2/3,\\
 |\operatorname{disc}W|
   &=3^{(2n-1)m^2/3-m}|D_0|^{m^2/3},\\
 \operatorname{rd}(W)&=\operatorname{rd}(B)
   =3^{n-1-3/(2m)}\sqrt{|D|}.
 \end{split}                                                    \tag{1}
\]
The only ramified rational primes are \(3\) and the divisors of \(D_0\).
At three the decomposition and inertia groups coincide and are cyclic of
order \(m/3\). At each prime dividing \(D_0\) inertia has order two and
the decomposition group is contained in a conjugate of the multiplier group.
For each fixed sign of \(D\), the distinct such Galois fields of absolute
discriminant at most \(Y\) number at least
\[
                 \gg_m \frac{Y^{6/(m^2(m+1))}}{(\log Y)^2}.
                                                               \tag{2}
\]
They are totally real for \(D<0\), and totally imaginary for \(D>0\).
The corresponding positive-rank elliptic twist of SD2 still has conductor
\(12D^2\). On the genuine-prime subfamily of SD3 its moving-prime
Frobenius and ideal-order conclusions hold, with the logarithmic loss in
(2) replaced by a fixed power as in that corollary.

**Proof.** SD2 gives the group, the degree and an everywhere unramified
cyclic extension \(W/B\) of degree \(m\), completely split above three.
The discriminant tower formula therefore gives
\[
 |\operatorname{disc}W|=|\operatorname{disc}B|^m,
 \qquad \operatorname{rd}(W)=\operatorname{rd}(B).
\]
Substitute SD2's formula
\( |\operatorname{disc}B|=3^{(2n-1)m/3-1}|D_0|^{m/3}\)
and \([B:\mathbb Q]=2m/3\). This proves (1).

Three splits in \(M\). Since \(R/\mathbb Q\) is cyclic and totally
ramified at three, and \(W/B\) is completely split there, a completion
of \(W\) at three is exactly \(R_3\). This proves the stated local group.
At a prime \(p\mid D_0\), which is odd and different from three,
SD2's inertia calculation gives an affine involution. Conjugate it to
\(x\mapsto-x\) on \(\mathbb Z/m\). An affine map \(x\mapsto ax+b\)
commutes with this involution precisely when \(2b=0\); since \(m\) is
odd, this says \(b=0\). The normalizer of its order-two group is its
centralizer, hence is the cyclic multiplier group. Inertia is normal
in the decomposition group, so the latter lies in this normalizer and
is cyclic. At every other finite prime the extension is unramified,
so its decomposition group is cyclic. The archimedean groups have order
at most two. This proves local cyclicity and the ramification assertions.

For \(D>0\), the mirror \(M\) is imaginary, so the Galois field \(W\)
is totally imaginary. For \(D<0\), \(B\) is totally real and the odd
Galois extension \(W/B\) cannot acquire a complex completion. Thus
\(W\) is totally real.

On either sign SD2 supplies
\(\gg_m X^{2/(m+1)}/(\log X)^2\) distinct \(D\) with \(|D|\le X\).
Equation (1) writes \(|\operatorname{disc}W|=c_m|D|^{m^2/3}\)
for a positive constant \(c_m\). Different absolute values of \(D\)
give different discriminants, hence different fields. Substitute
\(X=(Y/c_m)^{3/m^2}\) to obtain (2). The remaining assertions retain
the same parameters and follow directly from SD2–SD3. \(\square\)

For example, at \(m=9\), these are degree-54 Galois fields with
\[
 |\operatorname{disc}W|=3^{72}|D/3|^{27},\qquad
 \operatorname{rd}(W)=3^{5/6}\sqrt{|D|},\qquad
 \#\{W:|\operatorname{disc}W|\le Y\}\gg Y^{1/135}/(\log Y)^2.
\]
The discriminant here belongs to the normal closure; the degree-nine root
field has the different discriminant recorded in SD2.

## The fixed base and the optimization limit

**Proposition IG2.** Over \(R\), every field in IG1 is a tame locally
cyclic extension with dihedral group \(C_m\rtimes C_2\) of order \(2m\),
completely split above three, and relative discriminant ideal
\[
                       \mathfrak d_{W/R}=(D_0\mathcal O_R)^m.
                                                               \tag{3}
\]
Fields with distinct quadratic mirrors are pairwise linearly disjoint over
\(R\). They all contain the nontrivial field \(R\), so are not pairwise
linearly disjoint over \(\mathbb Q\). No sequence of IG1 fields at levels
\(n\to\infty\) has bounded root discriminant.

**Proof.** The subgroup fixing \(R\) consists of the translations and
the two multipliers \(\{1,-1\}\). Its group is dihedral, and its unique
quadratic subfield over \(R\) is \(B\). The quadratic extension \(B/R\)
has relative discriminant \(D_0\mathcal O_R\): the discriminants of
\(M/\mathbb Q\) and \(R/\mathbb Q\) are coprime. Apply the relative
discriminant tower formula to the unramified extension \(W/B\) to get
(3). Its ramification indices are two, at primes of odd residue
characteristic. Local cyclicity passes to a base change from \(\mathbb Q\)
to an intermediate field, and the calculation at three in IG1 gives
complete splitting over \(R\).

Every nontrivial quotient of a dihedral group of order twice an odd
prime power has a quadratic quotient. Indeed any proper normal subgroup
is contained in the rotation group: a normal subgroup containing a
reflection contains all rotations, because two is invertible modulo
\(m\). A nontrivial intersection of two \(W\)'s over \(R\) would
therefore contain the distinguished quadratic subfield of each, forcing
their \(B\)'s to coincide. Distinct quadratic fields \(M\) remain
distinct after adjoining the odd-degree field \(R\), so this is
impossible. This proves pairwise linear disjointness over \(R\).

Finally (1) gives
\(\operatorname{rd}(W)\ge3^{n-1-3/(2\cdot3^n)}\), which tends to
infinity. \(\square\)

Thus the relative unramified extension adds degree at no root-discriminant
cost; increasing the cyclotomic level does have a cost. At fixed level,
optimization in this family amounts to minimizing the actual quadratic
discriminant subject to SD1. Neither a smallest example nor a root-discriminant
record is asserted here. Requiring all ramification to lie in a fixed finite
rational prime set also leaves only finitely many fields at fixed level:
there are finitely many possible quadratic mirrors with that support, and
each fixed \(B\) has only finitely many unramified cyclic degree-\(m\)
extensions by its finite ordinary class group.

## Primary-literature comparison

[Kim–König, Theorem 5.1](https://arxiv.org/pdf/1911.08742v3), PDF p.13,
states the Shafarevich consequence that every finite solvable group over a
number field has infinitely many linearly disjoint tame realizations whose
decomposition groups at ramified primes are cyclic and equal inertia.
The theorem statement and the following split-embedding discussion were
inspected. Thus qualitative locally cyclic realization of \(G_m\) is
already known. IG1 instead supplies the displayed discriminant/count and
its marked quadratic/elliptic data; it is wild at three and shares \(R\).
This comparison does not assert historical novelty for that package.

[Alberts, Corollaries 1.5 and 1.7](https://arxiv.org/pdf/1907.06289v3),
PDF pp.8–9, gives counting bounds for a transitive group with abelian normal
subgroup after fixing a solvable embedding problem, and lower bounds for
all fields once a realization exists. The exact statements were inspected.
They already apply to broad counting questions for these affine groups.
IG1 counts a restricted family with an unramified relative translation
extension and additional arithmetic data. No improvement of the unrestricted
counting bounds or proof of a Malle asymptotic is claimed.

[Dèbes, Theorem 1.1](https://arxiv.org/pdf/1404.4074), PDF p.3, gives
quantitative Galois specialization with unramified Frobenius prescriptions
at suitably large primes, conditional on a regular realization over
\(\mathbb Q\). Its exact hypotheses and discriminant exponent were
inspected. The radical descent here contains \(R\): simultaneous
inversion fixes the real cyclotomic constants. At \(n\ge2\) this is
not a \(\mathbb Q\)-regular cover, so that theorem cannot simply be
applied to this particular cover. No assertion that the abstract group
has no other regular realization is intended.

For IG3, Kim–König's Proposition 4.2 and its proof, PDF pp.8–9, are a close
methodological comparison: retain local completions by approximation and
keep new inertia involutive. Its regularity and locally abelian hypotheses
do not cover the nonabelian local group constructed here. Alberts's
Theorem 5.3, PDF pp.44–45, supplies counting with local conditions for a
fixed abelian-kernel embedding problem. The unramified translation condition
outside three in IG3 is a restriction at all other primes, while its
quadratic mirror varies. No direct application of those counting theorems
to this entire restricted family has been established here.

These are comparisons, not additional proof inputs to IG1–IG3. The
quadratic specialization input to IG3 is the precise Kulkarni–Levin theorem
cited in its proof. No historical-priority claim is made.

## Recovering the local free depth from the trace

Put \(F=\mathbb Q_3\), \(K_0=F(\pi)\), \(\pi^2=-3\),
\(\zeta=(-1+\pi)/2\) and \(\eta=(1+\pi^3)/(1-\pi^3)\), as in LD2.
Normalize \(v_{K_0}(\pi)=1\), so its restriction to \(F\) is \(2v_3\).

**Lemma IG3a.** Let \(\gamma=\zeta^j\eta^b\), where \(j\ne0\bmod3\)
and \(0\ne b\in\mathbb Z_3\). For \(T=\gamma+\gamma^{-1}\),
\[
             v_3(T-2)=1,\qquad v_3(T+1)=v_3(b)+2.                 \tag{4}
\]
Conversely, if \(T\in\mathbb Q_3\), \(T^2-4\) has squareclass \(-3\),
\(v_3(T-2)=1\), and \(T\ne-1\), either norm-one root of
\(X^2-TX+1\) has nonzero torsion component and free exponent of valuation
\(v_3(T+1)-2\).

**Proof.** The norm-one filtration in LD2 gives
\(v_{K_0}(\eta^b-1)=3+2v_3(b)\). Write \(u=\eta^b\).
The torsion component makes \(v_{K_0}(\gamma-1)=1\), and
\((\gamma-1)^2/\gamma=T-2\) gives the first equality. Since
\(\zeta^j+\zeta^{-j}=-1\),
\[
        T+1=(u-1)(\zeta^j-\zeta^{-j}/u).
\]
The second factor has valuation one: its summand
\(\zeta^j-\zeta^{-j}\) has valuation one and its correction from \(u\)
has valuation at least three. The product therefore has valuation
\(4+2v_3(b)\), proving (4).

For the converse the two roots are units reducing to one, and the first
identity forces depth one. The LD2 decomposition then has \(j\ne0\),
since its free subgroup starts at depth three. The condition \(T\ne-1\)
excludes a pure torsion root. Apply the equality just proved. \(\square\)

## Global realization of the nonsplit local action

Fix \(n\ge3\), \(1\le s\le n-2\), and put
\[
 m=3^n,\quad r=n-s,\quad A=2\cdot3^{s+2}-3,\quad
 T=2+A t^m,\quad H(t)=A t(4+A t^m).
                                                                    \tag{5}
\]
Use rational parameters on the fixed local slice
\[
             t\equiv1\pmod3,\qquad v_p(t)=0\quad(p\mid2A).          \tag{6}
\]
For the positive sign of \(H\) take \(t>0\); for its negative sign take
\(-(4/A)^{1/m}<t<0\). Restrict to quadratic
\(K=\mathbb Q(\sqrt{H(t)})\ne\mathbb Q(\zeta_3)\), and set
\[
 \gamma=(T+\sqrt{T^2-4})/2,\quad D=\operatorname{disc}K,\quad
 D_0=D/3,\quad M=\mathbb Q(\sqrt{-D_0}),\quad
 B=M\mathbb Q(\zeta_m)^+.
\]

**Theorem IG3.** The degree-\(m\) root field \(L\) of
\(\mathcal D_m(X,1)-T\) has normal closure \(W/\mathbb Q\) with full
affine group \(G_m=C_m\rtimes(\mathbb Z/m)^\times\). At three its
decomposition group, in this affine action, is conjugate to the LD3 subgroup
\[
 \Gamma_{n,s}=\left\{
 i\longmapsto ai+(\widetilde a-1)/3+3^s h\pmod {3^n}:
 a\equiv1\pmod3,\ h\bmod3^r\right\}.                             \tag{7}
\]
Here \(a\bmod3^n\) and \(\widetilde a\) is any lift congruent to one
modulo three. Changing the lift does not change the set, since
\(s\le n-2\). This is a nonabelian nonsplit extension
\[
  1\longrightarrow C_{3^r}\longrightarrow\Gamma_{n,s}
    \longrightarrow C_{3^{n-1}}\longrightarrow1,
       \qquad |\Gamma_{n,s}|=3^{2n-s-1}.                           \tag{8}
\]
At every rational prime different from three, inertia in \(W/\mathbb Q\)
is trivial or an affine involution; the decomposition groups there are
cyclic. In particular two is unramified. The extension \(W/B\) is cyclic
of order \(m\) and unramified away from three, with local cyclic degree
\(3^r\) at each of the two primes of \(B\) above three.

Both absolute discriminants are exact:
\[
 \begin{split}
 |\operatorname{disc}L|
      &=3^{\delta_n}|D_0|^{(m-1)/2},
       &\delta_n&=((2n+1)m-1)/2,\\
 |\operatorname{disc}W|
      &=3^{E_{n,s}}|D_0|^{m^2/3},
       &E_{n,s}&=\frac{(2n+1)m^2}{3}
           +\frac{3^{2n-s-2}-21\cdot3^s}{4}.
 \end{split}                                                       \tag{9}
\]
On each sign the distinct objects satisfy
\[
 \begin{array}{c|c}
 \text{bound}&\text{lower count}\\\hline
 |D|\le X&\gg_{n,s} X^{2/(m+1)}/(\log X)^2\\
 |\operatorname{disc}L|\le Y&\gg_{n,s} Y^{4/(m^2-1)}/(\log Y)^2\\
 |\operatorname{disc}W|\le Y&\gg_{n,s}Y^{6/(m^2(m+1))}/(\log Y)^2.
 \end{array}                                                       \tag{10}
\]
The fields \(W\) are totally real when \(D<0\) and totally imaginary
when \(D>0\). All counts are for fixed \(n,s\), with constants allowed
to depend on both.

**Proof.** First observe that \(A\) is a positive odd integer with
\(v_3(A)=1\) and \(A/3\equiv-1\bmod3\). By (6) and the binomial
valuation formula,
\[
 v_3(t^m-1)\ge n+1,\qquad
 T+1=2\cdot3^{s+2}+A(t^m-1).
\]
The second term has valuation at least \(n+2>s+2\). Consequently
\[
 v_3(T-2)=1,\quad v_3(T+1)=s+2,\quad
                   (T^2-4)/(-3)\equiv1\pmod3.                   \tag{11}
\]
Thus \(K_3=K_0\). Lemma IG3a writes either local branch as
\(\zeta^j\eta^b\), with \(j\ne0\) and \(v_3(b)=s\).
Inverting the branch if needed gives \(j=1\). LD3 now supplies precisely
(7), the local degree and conductor spectrum. The free exponent need not
equal \(3^s\); its valuation is the datum required by LD3. Multiplication
of the root labels by \(-1\) conjugates the two torsion orientations.

The group in (7) contains \(\sigma:i\mapsto4i+1\) and
\(\tau:i\mapsto i+3^s\), with
\(\sigma\tau\sigma^{-1}=\tau^4\ne\tau\), since \(r\ge2\).
Every lift of the quotient generator \(4\bmod3^n\) has translation
part \(b\equiv1\bmod3^s\). Its \(3^{n-1}\)-st power translates by
\(b(4^{3^{n-1}}-1)/3\), whose valuation is \(n-1\), so is nonzero
modulo \(3^n\). No such lift has order \(3^{n-1}\), proving nonsplitting.

At every prime dividing \(A\), the trace is integral by (6), and its
norm-one roots are units. At a prime \(p\nmid A\), an integral \(t\)
again gives unit roots; if \(t\) has a pole, the trace has valuation
\(m v_p(t)<0\). The two roots have the distinct integral Newton slopes
\(m v_p(t),-m v_p(t)\); the quadratic algebra splits there, and these
valuations are divisible by \(m\). Hence \(\gamma\) is \(m\)-virtual
at all finite primes. At two, \(T\) is odd and \(T^2-4\equiv5\bmod8\),
so \(K_2/\mathbb Q_2\) is unramified quadratic. In particular \(D\)
is odd, \(v_3(D)=1\), and \(-D_0\) is a fundamental discriminant
coprime to three; (11) makes three split in \(M\).

The local torsion component \(j\ne0\) shows that \(\gamma\) is not
a cube in \(K_3\), hence not in \(K\). As \(K\ne\mathbb Q(\zeta_3)\),
it is disjoint from \(\mathbb Q(\zeta_m)\), whose unique quadratic
subfield is \(\mathbb Q(\zeta_3)\). A cubic radical of \(\gamma\)
over \(K\) is nonnormal, since \(K\) lacks \(\zeta_3\). It cannot
lie in the abelian extension \(K(\zeta_m)/K\). Thus \(\gamma\) remains
noncubic there, and adjoining \(\theta=\gamma^{1/m}\) has degree \(m\).
The simultaneous involution \(\theta\mapsto\theta^{-1}\),
\(\zeta_m\mapsto\zeta_m^{-1}\), nontrivial on \(K\), is central
and fixes the \(m\) Dickson roots. Its fixed field is \(W\), with
full affine action and translation-fixed field \(B\), exactly as in
the Kummer descent of SD2. This establishes the global group, rather
than inferring it from local irreducibility.

Away from three, virtuality makes the radical over \(K(\zeta_m)\)
unramified. Cyclotomic ramification is absent there. Therefore rational
inertia is trivial or a reflection, and the relative translation extension
\(W/B\) is unramified. Two is unramified already in the quadratic and
cyclotomic factors. The affine normalizer calculation in IG1 proves
cyclicity of the other decomposition groups away from three. At three
the normal closure of the polynomial over \(\mathbb Q_3\) is its
completion in \(W\); LD3 identifies its translation subgroup over the
local real cyclotomic field with \(C_{3^r}\). This gives (8) and the
stated relative local degree at both primes of \(B\).

For completeness, the root-field different can also be calculated
directly. The polynomial translated by \(X=Z-1\) is Eisenstein at three:
\(\mathcal D_m(X,1)\equiv X^m\bmod3\),
\(\mathcal D_m(-1,1)=2\), and \(v_3(2-T)=1\).
At a root, differentiation of the Dickson identity gives
\[
 (\mathcal D'_m(X,1))^2=m^2(T^2-4)/(X^2-4).
\]
In the degree-\(m\) completion, \(Z\) is a uniformizer,
\(v(X^2-4)=v((Z+1)(Z-3))=1\), and \(v(T^2-4)=m\).
The different exponent is therefore \(nm+(m-1)/2=\delta_n\).
Each prime dividing \(D_0\) has a tame reflection with one fixed root,
giving root-field discriminant exponent \((m-1)/2\). This proves
the first formula in (9).

To calculate the normal-closure discriminant, use the local character
table from LD3. A character of exact order \(3^k\) of the local
\(C_{3^r}\) has conductor exponent
\[
 a_1=0,\qquad a_k=(3^{k-1}+1)/2\ (2\le k<r),\qquad
                         a_r=(3^n+1)/2.                         \tag{12}
\]
Restriction of global \(C_m\) characters to this local subgroup is
\(3^s\)-to-one, including the trivial character. Thus at each prime
of \(B\) over three the
[conductor–discriminant formula](sources/milne-cft403.pdf), Milne V,
Theorem 3.27 on printed p.165, gives exponent
\[
 C_{n,s}=3^s\left(
        3^{r-1}(m+1)+\sum_{k=2}^{r-1}3^{k-1}(3^{k-1}+1)\right).
                                                                    \tag{13}
\]
The empty sum is zero. Both base primes have norm three. The cyclotomic
calculation in LD2 gives discriminant
\(3^{((2n-1)m/3-1)/2}\) for \(\mathbb Q(\zeta_m)^+\), which ramifies
only at three. Its discriminant is coprime to that of \(M\), so their
discriminant product gives
\[
 |\operatorname{disc}B|=3^{(2n-1)m/3-1}|D_0|^{m/3}.
\]
The discriminant tower formula now yields
\[
 |\operatorname{disc}W|=|\operatorname{disc}B|^m3^{2C_{n,s}}.
\]
Summing the two geometric series in (13) gives the expression
\(E_{n,s}=m((2n-1)m/3-1)+2C_{n,s}\) in (9). The factor \(3^s\)
in (13) accounts for the smaller local translation subgroup and cannot
be omitted. For instance, the root discriminant of \(W\) is
\(3^{3E_{n,s}/(2m^2)}\sqrt{|D_0|}\).

Finally, \(H\) is squarefree of degree \(m+1\): its root zero is
simple, and its other roots are the simple nonzero roots of
\(4+A t^m\). The finite local conditions (6) are nonempty open sets
and can be combined with either specified real interval by weak
approximation. Choose a rational center in the product of these opens
and apply
[Kulkarni–Levin, Theorem 1.2](sources/kulkarni-levin-hilbert-irreducibility-2111.15582.pdf)
to the degree-two map given by the parameter minus that center on
\(y^2=H(t)\). The actual theorem permits a fixed finite set of places,
a fixed small neighborhood at each, and a fixed thin exclusion; its
Theorem 1.2 and proof on printed pp.2 and 6–7 were inspected. Translation
changes heights by a fixed multiplicative constant. Hence parameter
height at most \(U\) gives \(\gg_{n,s}U^2/(\log U)^2\) distinct
quadratic fields satisfying the local conditions. Discarding
\(\mathbb Q(\zeta_3)\) removes at most one distinct field.

For \(t=a/b\) in reduced form, its quadratic field is represented by
the integer
\[
                         A a(A a^m+4b^m),
\]
since \(m+1\) is even. Its absolute discriminant is
\(O_{n,s}(U^{m+1})\). This proves the first count in (10). On either
sign, distinct absolute values of \(D\) distinguish both the root
fields and their normal closures by (9). Substituting the exact
discriminant powers gives the other counts. The signatures follow
from the sign of the quadratic mirror and the odd Galois degree
\([W:B]=m\), as in IG1. \(\square\)

For the first nonabelian case \(n=3,s=1\), one has \(A=51\),
\(|\Gamma_{3,1}|=81\), \([W:\mathbb Q]=486\), and
\[
 |\operatorname{disc}L|=3^{94}|D/3|^{13},\qquad
 |\operatorname{disc}W|=3^{1692}|D/3|^{243}.
\]
The parameters \(t=1\) and \(t=-1/5\) give respectively
\(D=2805\) and \(D=-1519918441772460934899\). Each sign has
\(\gg Y^{1/3402}/(\log Y)^2\) distinct degree-486 normal closures.
These discriminant examples do not assert independent Galois-group
computation in degree 486; that group identification uses the proof above.

## Evidence scope

IG1–IG2 are written corollaries using the current SD2 theorem. Their new
steps are the affine normalizer calculation, discriminant conversion and
dihedral intersection argument. They require no new finite search.
The [SD checker](verify_cyclotomic_mirror_small_dickson_trace.py) checks
the input identities and representative field discriminants; it does not
prove the unrestricted counting theorem or the new group arguments.

The [IG3 checker](verify_dickson_prescribed_nonsplit_realizations.py) was
internally executed from the repository root with
`sage -python research/arithmetic/verify_dickson_prescribed_nonsplit_realizations.py`,
under SageMath 10.9 and PARI 2.17.3. It passed both torsion orientations
and free-unit trace fixtures of depths zero through three, parameter and
Eisenstein fixtures at \((n,s)=(3,1),(4,1),(4,2),(5,3)\), affine groups
and character restriction counts for the first three of those pairs, two
degree-27 root-field discriminants on opposite signs, and the degree-81
root-field discriminant at \((4,2),t=1\). Their discriminant exponents
at three are respectively \(94,94,364\).

One additional degree-81 compositum of the raw local representative
\(\mathcal D_{27}(X,1)-683/343\) with the real 27th cyclotomic layer
has discriminant \(3^{282}7^{72}\), agreeing with the conductor conversion
\(1692=(486/81)\cdot282\). Its quadratic field is
\(\mathbb Q(\zeta_3)\), so this is a control for the local normal
closure, not a parameter claimed to belong to the global family.
The all-level local-action identification, all-prime arguments and
asymptotic counts depend on the written proof and the named imported
theorems, not on these finite checks.
No independent replay, formal verification or external review is asserted.
