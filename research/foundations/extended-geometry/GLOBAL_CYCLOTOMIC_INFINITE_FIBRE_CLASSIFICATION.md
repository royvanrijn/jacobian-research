# A global cyclotomic classification of infinite odd-height fibres

Fix an odd positive integer \(a\) and a characteristic-zero field \(K\).
For a rational function \(u\in K(x)\) with \(u(0)=0,\ u'(0)=1\), its fibre
consists of the polynomial pairs \(C(0)=-D(0)=1\) for which
\[
 p=t^a(C+t^2D)/x,\qquad t^2+2xy=1,
\]
has all positive pure moments zero and selected endpoint coordinate
exactly \(u\). No bound on the degrees of \(C,D\) is imposed.

This note characterizes **all rational coordinates with infinite fibres**.
It removes the assumption that \(u/(xu')\) is polynomial from the
[Davenport--Zannier rigidity theorem](DAVENPORT_ZANNIER_COORDINATE_FIBRES_AND_SHABAT_RIGIDITY.md).
The answer is a polynomial composition with explicit ramification bounds
and a finite cyclotomic compatibility test. In particular, neither a
squarefree jet polynomial nor polynomial reciprocal units alone suffice.

The conclusion is a parameterization by rational maps and finite
algebraic conditions. It is not a list by coordinate degree, a
classification of arbitrary sphere phases, or a classification up to
all automorphisms of the sphere. Finite fibres without a polynomial
reciprocal unit remain governed by the retained fixed-coordinate theorem.
No historical priority, formal verification, or external review is claimed.

## 1. Normal form and the finite geometric data

First work over an algebraically closed characteristic-zero field.
Choose
\[
 R\in K[x],\quad R(0)=1,\quad R'(0)\ne0,\quad
 R\text{ nonconstant and not a proper polynomial power},
 \tag{G1}
\]
and
\[
 U\in K(z),\quad U(1/z)=-U(z),\quad U(1)=0,\quad
 U'(1)R'(0)=1.
 \tag{G2}
\]
Thus \(U\) has a simple zero at \(1\). Set \(u=U(R)\).
For every finite \(\zeta\), define
\[
 e_\zeta=\max_{R(\beta)=\zeta}
                  \operatorname{ord}_\beta(R-\zeta).
 \tag{G3}
\]
Every such fibre is nonempty and finite. The maximum, not an average
or a generic multiplicity, is required.

Define the finite **active set** \(T=T_{a,R,U}\) as follows.

| Value of \(U\) at \(\zeta\) | Local data | Include \(\zeta\) in \(T\) when |
|---|---|---|
| \(U(\zeta)=0\) | zero order \(s_\zeta\ge1\) | \(e_\zeta s_\zeta>a\) |
| \(U(\zeta)\ne0,\infty\) | \(f_\zeta=\operatorname{ord}_\zeta(U-U(\zeta))\ge1\) | \(e_\zeta f_\zeta>1\) |
| \(U(\zeta)=\infty\) | pole order \(m_\zeta\) | never |

Only finite zeros and critical points of \(U\), and finite critical
values of \(R\), can contribute. No condition at the source point
\(x=\infty\) is added.

The geometric conditions are
\[
 T\subset\mu_\infty,
 \tag{G4}
\]
where \(\mu_\infty\) denotes the roots of unity, and
\[
 \begin{array}{ll}
 e_\zeta(s_\zeta-a)\le a,
    &\zeta\in T,\quad U(\zeta)=0,\\[2pt]
 e_\zeta(f_\zeta-1)\le1,
    &\zeta\in T,\quad U(\zeta)\ne0,\infty.
 \end{array}
 \tag{G5}
\]
In particular, an active value \(0\) violates (G4).

Here is the remaining finite data. For each finite pole value \(\alpha\)
of \(U\), put
\[
 q_\alpha=\left\lfloor\frac{m_\alpha}{a+2}
                                  +\frac1{e_\alpha}\right\rfloor.
 \tag{G6}
\]
Choose one representative from each inversion pair
\(\{\alpha,\alpha^{-1}\}\) of finite nonzero poles, excluding the fixed
values \(1,-1\). The value \(1\) is not a pole. Put
\[
 h_\alpha(z)=\frac{z-\alpha}{1-\alpha z},\qquad
 H_{\mathbf n}(z)=\prod_\alpha h_\alpha(z)^{n_\alpha},
 \qquad -q_\alpha\le n_\alpha\le q_{\alpha^{-1}}.
 \tag{G7}
\]
The empty product is \(1\). Also set
\[
 k_{\min}=
 \begin{cases}
 -q_0,&U\text{ has a pole at }0,\\
 0,&U\text{ is regular at }0.
 \end{cases}
 \tag{G8}
\]
There are exactly
\[
 B_{a,R,U}=\prod_\alpha
                   (q_\alpha+q_{\alpha^{-1}}+1)
 \tag{G9}
\]
twist vectors, with empty product \(1\).

If (G4) holds, let
\[
 d=\operatorname{lcm}_{\zeta\in T}\operatorname{ord}(\zeta),
 \qquad d=1\text{ when }T=\varnothing.
 \tag{G10}
\]
For a twist vector, test the following finitely many equalities for one
integer \(k_0\in\{0,\ldots,d-1\}\):
\[
                  H_{\mathbf n}(\zeta)\zeta^{k_0}=-1
                       \quad\text{for every }\zeta\in T.
 \tag{G11}
\]
All values are defined and nonzero: \(T\) avoids poles of \(U\), and
(G4) excludes \(0\). Checking only that each value in (G11) is a root of
unity is insufficient; the **same** exponent must work at every value.

## 2. The global theorem

**Theorem GlobalCyclotomicInfiniteFibre.** A normalized rational
coordinate \(u\) has an infinite odd-height-\(a\) fibre if and only if
it is obtained from (G1)--(G2), satisfies (G4)--(G5), and at least one
twist vector satisfies (G11).

For any coordinate satisfying these conditions, \(R\) is its unique
normalized polynomial reciprocal generator, and \(U\) is uniquely
determined by \(u=U(R)\). The entire fibre, not merely a subfamily, is
the disjoint union of the following infinite sequences. For each
successful twist let \(k_*\) be the least integer at least \(k_{\min}\)
congruent to its solution \(k_0\) modulo \(d\). For \(j\ge0\), put
\[
 \eta_j=H_{\mathbf n}(R)R^{k_*+dj},\qquad
 b=\frac{u}{xu'},\qquad X_j=\frac b2(1+\eta_j),
 \tag{G12}
\]
\[
                 C_j=\frac{xX_j^a}{u},\qquad
                 D_j=-\frac{xX_j^{a+2}}u .
 \tag{G13}
\]
These are normalized polynomial pairs. Distinct displayed data give
distinct pairs. Each twist has at most one successful residue class
modulo \(d\). There are no additional isolated pairs in an infinite
fibre.

For any normalized coordinate that has a polynomial reciprocal
generator, the same construction gives (G1)--(G3), (G6)--(G9).
If (G4) or (G5) fails, its entire fibre has at most \(B_{a,R,U}\)
elements: at most one per twist. If the geometric conditions hold but
every test (G11) fails, the fibre is empty. Thus the theorem also gives
an explicit cardinality bound on the finite side when \(R\) exists.
It does not list the possible isolated exponents when (G4)--(G5) fail;
the retained fixed-coordinate and effective theorems address that task.

**Base fields.** Over an arbitrary characteristic-zero \(K\), the same
theorem holds with \(R\in K[x]\), \(U\in K(z)\), and the additional requirement
\[
                         H_{\mathbf n}\in K(z).
 \tag{G14}
\]
The ramification and torsion tests are geometric, over \(\overline K\).
For the geometric fibre omit (G14). A geometric successful twist is
not automatically asserted to descend to \(K\). The canonical \(R\)
itself always descends, by its uniqueness.

## 3. The retained interface and why the quotient is canonical

The [fixed-coordinate classification](ODD_START_FIXED_COORDINATE_FIBRE_CLASSIFICATION.md)
gives a bijection between pairs and normalized rational reciprocal units
\(\eta\). On the marked component of \(u(x)+u(y)=0\) these satisfy
\[
 \eta(0)=1,\qquad \eta(x)\eta(y)=1.
 \tag{G15}
\]
Their finite poles lie at poles of \(u\), with order at most
\(1+\lfloor m/(a+2)\rfloor\) above a pole of order \(m\). Their remaining
condition is \(\eta\equiv-1\pmod J\), where
\[
 \operatorname{ord}_\beta J=
 \begin{cases}
 \lceil\operatorname{ord}_\beta(u)/a\rceil-1,
       &\beta\ne0,\ u(\beta)=0,\\
 \operatorname{ord}_\beta u',
       &u(\beta)\ne0,\infty,\\
 0,&\text{otherwise}.
 \end{cases}
 \tag{G16}
\]
Only positive orders contribute. Formulae (G12)--(G13) recover the pair.
This is an all-order endpoint statement, not an inference from a finite
moment experiment.

That theorem also proves that an infinite fibre has a nonconstant
polynomial reciprocal unit. All such units are powers of one normalized
polynomial \(R_0\), every rational reciprocal unit belongs to \(K(R_0)\),
and \(u=U_0(R_0)\) with \(U_0(1/z)=-U_0(z)\). These assertions use the
marked component and the intermediate-field involution; they do not
assume rational return on the original \(x\)-line.

The canonical \(R_0\) cannot be a proper polynomial power. If
\(R_0=S^m,\ m>1\), normalize \(S(0)=1\). On the marked integral curve,
\((S(x)S(y))^m=1\), so \(S(x)S(y)\) is a constant root of unity.
Its marked value is \(1\); hence \(S\) is a smaller polynomial
reciprocal unit, a contradiction.

Conversely, (G1)--(G2) make \(R\) a reciprocal unit. Indeed the marked
return \(j(x)\) satisfies \(u(j(x))=-u(x)\). Both \(R(j(x))\) and
\(1/R(x)\) have constant term \(1\), and their images under \(U\)
equal \(-u(x)\). Since \(U'(1)\ne0\), local inversion gives
\(R(j(x))=1/R(x)\). The retained cyclicity theorem gives \(R=R_0^m\).
The nonpower condition forces \(m=1\). Over an algebraic closure it is
equivalent to the gcd of the finite root multiplicities of \(R\)
being one. Normalization makes any polynomial root descend whenever
\(R\) does.

The nonpower condition is necessary for the stated canonical
parameterization and exhaustion. For example,
\[
 u=\frac2m\frac{(1+x)^m-1}{(1+x)^m+1},\qquad m>1,
 \tag{G17}
\]
admits the displayed polynomial quotient \((1+x)^m\), but its
canonical reciprocal generator is \(1+x\). Its rational reciprocal
units are not all rational functions of \((1+x)^m\).

## 4. Ramification replaces every full jet

Let \(\beta\) be finite, \(R(\beta)=\zeta\), and
\(e=\operatorname{ord}_\beta(R-\zeta)\).

If \(U\) has a zero of order \(s\) at \(\zeta\), then \(u\) has zero
order \(es\). Away from the marked point, the jet order is
\[
                         q=\lceil es/a\rceil-1.
 \tag{G18}
\]
It is positive somewhere above \(\zeta\) exactly when
\(e_\zeta s>a\). The marked point causes no exception to this test:
there \(s=e=1\), so (G18) already gives zero. For every \(e\),
\[
 q\le e
 \quad\Longleftrightarrow\quad
 es\le a(e+1)
 \quad\Longleftrightarrow\quad
 e(s-a)\le a.
 \tag{G19}
\]
If \(s\le a\), the last inequality holds for every \(e\). If \(s>a\),
it is strongest at the maximum \(e_\zeta\). Thus the first line of
(G5) is equivalent to all these jet inequalities.

If \(U(\zeta)\ne0,\infty\) and its local degree is \(f\), then
\(u-u(\beta)\) has order \(ef\). Characteristic zero gives
\[
                q=\operatorname{ord}_\beta u'=ef-1.
 \tag{G20}
\]
The positive jets occur exactly at \(ef>1\), and
\[
                  q\le e\quad\Longleftrightarrow\quad e(f-1)\le1.
 \tag{G21}
\]
The maximum \(e_\zeta\) again gives exactly the strongest condition.
Poles of \(U\) give no roots of \(J\).

The retained elementary order test says that the class of \(R\) has
finite multiplicative order modulo \(J\) precisely when, at every
root \(\beta\) of jet order \(q\),
\[
 R(\beta)\in\mu_\infty,\qquad
             \operatorname{ord}_\beta(R-R(\beta))\ge q.
 \tag{G22}
\]
For completeness, if \(\zeta^d=1\), the first nonzero local term of
\(R^d-1\) is \(d\zeta^{d-1}(R-\zeta)\). Its coefficient is nonzero
in characteristic zero. Therefore its order is exactly \(e\), not
an arbitrarily increasing function of \(d\). Conversely a finite-order
class has a root-of-unity residue and this same order inequality.

Equations (G18)--(G21) identify (G22) **exactly** with (G4)--(G5).
The resulting order is (G10). Furthermore, for any rational \(H\)
regular at \(\zeta\), the difference \(H(R)-H(\zeta)\) vanishes to
order at least \(e\). Once (G4)--(G5) hold, every full-jet condition
on \(H(R)R^k+1\) is consequently equivalent to its value at \(\zeta\).
No squarefree replacement of \(J\) has been made.

Two immediate geometric restrictions are worth spelling out.

**Corollary QuotientRamificationRestriction.** Under (G4)--(G5), every
finite critical point of \(U\) with finite nonzero value is a root of unity, is simply
ramified (\(f=2\)), and has only unramified preimages under \(R\).
If \(U\) is unramified at a point \(\zeta\) with \(U(\zeta)\) finite
and nonzero, \(R\) may
have arbitrary ramification there, provided an active value is a
root of unity and (G11) holds.

**Marked and zero values.** If \(e_1>a\), the entire fibre is empty:
\(1\in T\), but every normalized reciprocal unit has value \(1\)
at \(z=1\), whereas (G11) requires \(-1\). At \(z=0\), infinitude
requires the value to be inactive. If \(U(0)=0\) of order \(s\), this
means \(e_0s\le a\). If \(U(0)\ne0,\infty\), it means \(e_0=f_0=1\).
If \(U\) has a pole at \(0\), there is no jet condition there and
only the exponent lower bound (G8) remains.

## 5. Pole bounds and the simultaneous cyclotomic test

Every rational reciprocal unit has the form \(V(R)\), with
\(V(1)=1,\ V(1/z)=V(z)^{-1}\). Its divisor is anti-invariant under
inversion. Its orders at \(1,-1\) are zero; the orders at \(0,\infty\)
are opposite. A finite nonzero pole of \(V\) must be a pole of \(U\),
because every finite \(z\)-value has a finite preimage under \(R\).
The paired zero must likewise occur at the inverse pole value of \(U\).
These facts give the unique factorization
\[
                          V(z)=z^kH_{\mathbf n}(z).
 \tag{G23}
\]
Each factor in (G7) equals \(1\) at \(1\), so normalization removes
the remaining scalar.

At a pole value \(\alpha\), a pole of \(V\) of order \(h\) pulls back
to order \(eh\), while \(u\) has pole order \(em_\alpha\). Its exact
allowed bound is
\[
 eh\le1+\left\lfloor\frac{em_\alpha}{a+2}\right\rfloor
 \quad\Longleftrightarrow\quad
 h\le\frac{m_\alpha}{a+2}+\frac1e.
 \tag{G24}
\]
The maximum \(e_\alpha\) gives (G6). The sign of \(n_\alpha\)
specifies which member of the pair is a pole, giving (G7). At zero
only negative \(k\) create a finite pole, giving (G8); infinity
imposes no upper bound. This proves the complete twist box directly,
as in the retained effective classification, without importing its
Orbit Problem theorem.

The fixed values need no extra factors. Anti-invariance of \(U\)
implies that at \(-1\) it has either a zero or a pole of odd order:
in a local coordinate in which inversion is negation its Laurent
series is odd. At \(1\) its zero is already simple by normalization.
Reciprocal units have order zero at both fixed values, and
\[
                  V(-1)=(-1)^{\,k+\sum_\alpha n_\alpha}.
 \tag{G25}
\]
Thus an active value \(-1\) supplies a parity condition, not a new
pole factor. A pole of \(U\) at \(-1\) contributes no jet condition.

If (G4)--(G5) hold, Section 4 reduces the jet condition to (G11).
The evaluation tuple
\((\zeta)_{\zeta\in T}\) has exact order \(d\); its integer powers
form a cyclic subgroup of \(\prod_{\zeta\in T}K^\times\).
Equation (G11) asks whether
\((-H_{\mathbf n}(\zeta)^{-1})_{\zeta\in T}\) belongs to that subgroup.
A nonempty solution set is one residue class modulo \(d\).
The lower bound on \(k\) never removes the entire class, proving
infinitude and (G12)--(G13).

Conversely, an infinite fibre supplies its canonical \(R\) and the
finite box (G23). Some twist has infinitely many surviving exponents.
Two such exponents imply a positive power of \([R]\) equals \(1\)
modulo \(J\), so (G4)--(G5) hold. Section 4 then gives (G11).
This proves necessity.

If the geometric conditions fail, a twist cannot have two surviving
exponents. If \(R\) vanishes at a root of \(J\), then \(U\) is not a
pole there, \(k_{\min}=0\), and every positive exponent makes the
candidate vanish at that root instead of taking the value \(-1\).
Otherwise \([R]\) is a unit of infinite order, so the difference of
two successful exponents is impossible. This proves the bound (G9).
When the geometric conditions hold every successful twist is infinite,
proving the no-isolated-pairs assertion.

Over a nonclosed field, the canonical generator descends by the
retained uniqueness argument. The rational function \(U\) then
descends too. Multiplication by \(R^k\in K(x)^\times\) cannot repair
a twist that is not defined over \(K\). Since substitution into the
transcendental element \(R\) is injective, descent is exactly (G14).
This proves the asserted base-field version.

## 6. When there are no nontrivial pole twists

Suppose every interval in (G7) consists only of zero, so \(H=1\).
If \(T\) is empty, every exponent \(k\ge k_{\min}\) works. Otherwise
(G11) is solvable if and only if the orders of **all** values in \(T\)
have the same positive 2-adic valuation. When it is solvable,
\[
                             k\equiv d/2\pmod d.
 \tag{G26}
\]
Indeed a root of unity of order \(n\) has \(k\)-th power \(-1\)
exactly when \(n\) is even and \(k\equiv n/2\pmod n\). Such a \(k\)
has 2-adic valuation \(v_2(n)-1\), proving necessity. Conversely write
every order as \(2^s m\), with the same \(s\ge1\) and \(m\) odd.
Then \(d/2\) is congruent to \(n/2\) modulo every \(n\).
The simultaneous modulus is their lcm \(d\).

This simplification does not apply when pole twists are allowed.
Section 8 gives an infinite fibre whose active orders are \(2\) and
\(8\), with every successful branch requiring a nontrivial twist.

## 7. Arbitrarily thick jet polynomials in infinite fibres

Work over \(\mathbb Q(i)\). For any \(e\ge2\), put
\[
 R=i+(1-i)(1+x)^e,\qquad
 U(z)=\frac2{e(1-i)}\,\frac{z-1}{z+1},\qquad u=U(R).
 \tag{G27}
\]
The normalization (G1)--(G2) holds. The zeros of \(R\) are simple,
so \(R\) is not a proper power. Its only finite critical point is
\(x=-1\), with value \(i\) and ramification index \(e\).
All preimages of \(R=1\) and \(R=-1\) are simple. Since \(U\)
is Möbius, for every odd \(a\ge1\)
\[
                    J_{a,u}=(x+1)^{e-1},\qquad T=\{i\}.
 \tag{G28}
\]
There are no nontrivial pole twists: the sole pole of \(U\) is the
inversion-fixed point \(-1\). The entire fibre is
\[
                             \eta=R^{2+4j},\qquad j\ge0.
 \tag{G29}
\]
For \(e\ge3\), \(J\) is nonreduced, with arbitrarily large
multiplicity as \(e\) increases. Nevertheless \(R-i\) vanishes to
order \(e\), so every required jet is satisfied. Here \(b=u/(xu')\)
has a pole of order \(e-1\) at \(-1\); these examples genuinely lie
outside the polynomial-\(b\) theorem.

## 8. A family for which pole twists are essential

Work over \(K=\mathbb Q(\sqrt2)\). Let
\[
 \alpha=\sqrt2-1,\quad R=-1+2(1+x)^e,\quad e\ge2,
 \qquad
 U(z)=\frac{1-\sqrt2}{2e}
             \frac{z^2-1}{z^2-2\sqrt2 z+1},\quad u=U(R).
 \tag{G30}
\]
Again \(u(0)=0,u'(0)=1\), and \(R\) has simple zeros. The poles of
\(U\) are \(\alpha,\alpha^{-1}=\sqrt2+1\), both simple. Their
preimages under \(R\) are unramified, so \(q_\alpha=q_{\alpha^{-1}}=1\)
at every odd height. Also \(k_{\min}=0\). Thus the complete twist
list is \(H^{-1},1,H\), where
\[
                           H(z)=\frac{z-\alpha}{1-\alpha z}.
 \tag{G31}
\]
The ordinary critical points of \(U\) are
\(\zeta=(1+i)/\sqrt2\) and \(\zeta^{-1}\), both simply ramified
and both of order \(8\). The map \(R\) is unramified above them.
Its only finite critical value is \(-1\), where \(U\) has a simple
zero and the unique preimage has ramification index \(e\).
Consequently
\[
 T=
 \begin{cases}
 \{\zeta,\zeta^{-1},-1\},&a<e,\\
 \{\zeta,\zeta^{-1}\},&a\ge e.
 \end{cases}
 \tag{G32}
\]
All of (G4)--(G5) hold, including the zero-jet bound
\(\lceil e/a\rceil-1\le e\).

Direct evaluation gives \(H(\zeta)=i=\zeta^2\),
\(H(\zeta^{-1})=-i\), and \(H(-1)=-1\). Therefore the full fibre has
the two branches
\[
                  \eta=H(R)R^{2+8j},\qquad
                  \eta=H(R)^{-1}R^{6+8j},
                         \qquad j\ge0,
 \tag{G33}
\]
at every odd height. The additional untwisted branch
\[
                              \eta=R^{4+8j}
 \tag{G34}
\]
occurs exactly when \(a\ge e\). For \(a<e\), the active orders \(2\)
and \(8\) are incompatible for pure powers, but both nontrivial pole
twists solve the simultaneous condition. All three twists are
defined over \(K\), so no coefficient-field extension is needed for
the displayed pairs.

One exact cancellation makes the construction particularly visible:
\[
 1+H(z)z^2
       =\frac{(z+1)(z^2-\sqrt2 z+1)}{1-\alpha z}.
 \tag{G35}
\]
The quadratic factor cancels the two ordinary critical jets of \(U\);
the factor \(z+1\) cancels the full ramified zero jet after pullback
by \(R\). At the first branch's initial exponent \(2\), (G12) becomes
the polynomial
\[
 X=\frac{(R-1)/x\;(R-\alpha)(1+x)^{e+1}}
                            {2e\sqrt2\,\alpha}.
 \tag{G36}
\]
It has \(X(0)=1\). This is a rational pole-twist phenomenon, not an
unjustified assumption that every admissible reciprocal unit is
polynomial. At \(a=1\), its jet at \(x=-1\) has order \(e-1\),
again allowing arbitrary thickness.

## 9. Countercontrols and source scope

Several tempting weakenings are false.

* Root-of-unity values alone do not suffice. With \(R=1+x\),
  \(w=(z-1)/(z+1)\), and
  \[
             U(z)=2\left(w+\frac23w^3+\frac15w^5\right),
  \tag{G37}
  \]
  the ordinary critical points \(z=\pm i\) have local degree \(3\).
  Their jet order is \(2\), while \(R\) has local degree \(1\).
  Thus (G5) fails despite cyclotomic residues.
* Compatible individual root values do not guarantee a common
  exponent. The no-twist active orders \(2,4\) already give an
  inconsistent system in (G26).
* Ignoring the marked value \(1\) would incorrectly admit unmarked
  roots of \(R-1\) with multiplicity greater than \(a\).
* Ignoring \(0\) would permit a nonunit class of \(R\) modulo \(J\).
* Ignoring maximum ramification in (G5)--(G6) can miss the unique
  offending preimage of a quotient value.
* Using a noncanonical proper-power quotient can omit rational
  reciprocal units, as in (G17).

The essential imported result is the retained fixed-coordinate
reciprocal-unit theorem, including its exact local polynomiality and
cyclicity statements. Its primary field input is Milne, *Fields and
Galois Theory*, version 5.10, Theorem 9.19 and Lemma 9.24
(printed pages 117--118): Lüroth and the rational-map field degree.
The marked-component and endpoint hypotheses are retained unchanged.
The explicit pole-twist argument is reproved in Section 5; no general
Orbit Problem theorem, Skolem--Mahler--Lech theorem, polynomial Pell
equation, or Belyi classification is imported.

The closest earlier theorem in the repository is the polynomial-\(b\)
Davenport--Zannier/Shabat classification. It is a genuine subcase;
Sections 7--8 show why its rigidity cannot be extended by silently
assuming squarefree jets or omitting rational twists.

The accompanying exact checker verifies bounded ramification tables,
simultaneous exponent controls, construction and moment identities,
the two displayed families, and deliberately invalid variants.
It does not prove the unrestricted theorem by finite testing.
The written proof establishes necessity, sufficiency, exhaustion and
descent. Precise source pins and actual replay counts are retained in
the source ledger and receipt.

The author replay passed 320 zero-ramification and 320
ordinary-ramification cases, 2880 pole-bound cases, 672 simultaneous
untwisted order systems, and 240 quotient-ring residue cases.
It reconstructed 14 normalized polynomial pairs and checked their
84 raw moments for exponents \(1,\ldots,6\), including 42 nontrivial
even-exponent checks. Three thickness controls and nine deliberately
invalid variants also passed their stated acceptance/rejection tests.
The coefficients are exact in \(\mathbb Q\), \(\mathbb Q(i)\), and
\(\mathbb Q(\sqrt2)\); floating-point acceptance is not used.

Run from the repository root:

    python3 research/foundations/extended-geometry/verify_global_cyclotomic_infinite_fibres.py --repo-root . --output research/foundations/extended-geometry/global_cyclotomic_infinite_fibres_verification.json

The [source ledger](global_cyclotomic_infinite_fibres_source_ledger.json)
pins every imported proof/helper/source, and the
[receipt](global_cyclotomic_infinite_fibres_verification.json) binds the
current proof, checker and ledger. These checks are internal exact
controls, not formal or external verification.
