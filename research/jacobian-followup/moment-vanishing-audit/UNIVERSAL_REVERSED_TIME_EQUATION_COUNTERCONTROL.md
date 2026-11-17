# Countercontrol to a universal reversed time equation

16 November 2026. Internal written counterproof; no novelty, independent
replay, formal verification or external review is claimed. See the
[checker](verify_reversed_time_obstruction.py),
[receipt](reversed_time_obstruction_verification.json) and
[source audit](reversed_time_obstruction_source_audit.json).

The direct attempt to prove a scalar time equation for EVERY actual pure,
eventually safe reversed-apolar source fails. An explicit pure SAFE pair
with nested actual ideals and an excellent linear gap has NO such equation,
even allowing coefficients formal in the marker parameter. Consequently a
universal equation/finite-profile-recurrence producer for the UNMODIFIED
source is not a valid full-gap obligation.

This is a countercontrol to that proof route, NOT an eventual-but-nonlinear
pair. The remaining universal positive-gap implication stays OPEN. This
note records a failed universal approach, rather than another positive
proper-subclass or normal-form reduction.

## 1. An actual pure safe pair, with exact gap and nested ideals

Let K be any characteristic-zero field, and put

    phi(Z)=sum_(j>=2) Z^(j!), T=phi(partial_x), P=x.         (1)

Since ord(phi)=2, every positive pure output is0 and

    T^m(Qx^m)=0 whenever m>deg Q.                          (2)

Each factor differentiates at least twice, whereas deg(Qx^m)<=degQ+m.
At the boundary Q=x^m the unique derivative order2m gives

    T^m x^(2m)=(2m)! !=0.                                 (3)

Thus the exact least real marked degree cutoff is1. Also phi=Z^2 U(Z)
with U(0)=1. Formal-unit invariance gives the ACTUAL heights and deficits

    h_m=floor(m/2), b_m=ceil(m/2).                         (4)

The actual ideals Ann(P^m)=(Z^(m+1)) are nested. Every requirement of the
qualitative premise and positive-gap conclusion is already satisfied.

For a concise higher-dimensional transcription, take
P=x+y^D+w^D, D>=2, and let T act only in x. It has the same heights,
sharp cutoff1, and source formula below with x replaced by P. The phase
is concise: its three first derivatives1,Dy^(D-1),Dw^(D-1) are linearly
independent. Its ideals are nested because partial_x P=1. This extension
is only to show that ambient dimension/phase degree do not repair the
failed universal equation assertion; it is not a new safety programme.

## 2. The reversed source is explicit

For the proper two-index exponential source, exponential eigenfunctions give

    H(a,b,x)=sum_(k,m) T^k x^m a^k b^m/(k!m!)
            =exp(bx+a phi(b)).

This identity is coefficientwise finite. Under the pure reversed substitution
a=1/z,b=tz it gives the ACTUAL source from GVCREVERSEDHEIGHT20261116:

    Psi(t,z,x)=exp(tzx+phi(tz)/z).                         (5)

The exponent phi(tz)/z starts at t^2 z, so(5) is a proper formal series.
Every source-z coefficient is polynomial in time, of degree<=2 times its
z degree. Formula(4) gives the same support bound directly.

**GVCREVERSETIMEOBSTRUCTION20261116.** There is NO nonzero scalar ordinary
time equation for(5) over K(x)((z))(t). Equivalently there is none with
coefficients in K(x)[[z]][t] and finite time powers/derivative orders.
Hence no finite recurrence of common bounded m-degree over K(x)[[z]]
exists for its actual reversed profiles. The counterpair is already pure,
safe and nested, with b_m>=m/2; the obstruction is not lack of a linear gap.

We prove the nonexistence in full, rather than infer it from a failed search.

## 3. A semiclassical leading coefficient forces an algebraic derivative

First remove the time exponential exp(tzx). A time equation for Psi gauges
to one for G=exp(phi(tz)/z), with partial_t replaced by partial_t+zx.
The gauged operator is nonzero and has coefficients in K(x)[[z]][t].
The function G is independent of x. To descend the equation to K[[z]][t],
choose a coefficientwise K-linear functional on K(x) nonzero on one
specified nonzero coefficient of the operator. For example take a suitable
Laurent coefficient at x=0, as in the retained suspension descent proof.
Applying that functional coefficientwise in z commutes with the equation
because G is independent of x. This yields a NONZERO equation over
K[[z]][t]. No unjustified evaluation at a pole is used.

Now replace t by v/z. This is legitimate in K((z))[[t]], giving

    E(v,z)=exp(phi(v)/z) in K((z))[[v]],
    partial_t=z partial_v.

The transformed operator has the finite form

    L=sum_(j=0)^r A_j(v,z) partial_v^j,
    A_j in K((z))[v], L!=0.                               (6)

Each A_j has finite v degree and a common finite lower z order. Write its
lowest nonzero z term as z^(nu_j) a_j(v), a_j(v)!=0, omitting A_j=0.
Repeated differentiation of E gives the finite Bell-polynomial identity

    E^(-1) partial_v^j E
      =z^(-j)(phi'(v))^j + sum_(ell<j) z^(-ell) B_(j,ell)(v), (7)

with B_(j,ell) formal series in v. Let q=min_j(nu_j-j).
Divide LE=0 by E and extract its z^q coefficient. The result is

    sum_(j:nu_j-j=q) a_j(v)(phi'(v))^j=0.                 (8)

This is a NONZERO polynomial equation over K[v] in phi'(v): distinct j
are distinct powers of its indeterminate and each a_j is nonzero.
The extraction is valid in K[[v]]((z)), because the finite Bell expressions
in(7) have a common lower z bound; this ring injects into K((z))[[v]].
Thus a scalar equation would force phi' algebraic over K(v).

But the coefficient sequence of phi' has nonzero terms at j!-1 and
arbitrarily long internal zero blocks followed by more nonzero terms.
An algebraic series in characteristic zero is D-finite: differentiating
its algebraic equation keeps derivatives in one finite separable field.
A D-finite series has a polynomial-coefficient finite recurrence; beyond
the finitely many roots of its leading coefficient, a sufficiently long
zero block forces all later coefficients to be0. This contradicts those
factorial gaps. Hence phi' is not algebraic, contradicting(8).
This proves the stated nonexistence, over EVERY characteristic-zero field.

For P=x+y^D+w^D the gauge exp(tzP) is still a polynomial-time rank-one
gauge. The same coefficientwise descent from K(x,y,w) removes independent
output parameters, so the nonexistence persists. Under a passive linear
suspension the exact gauge exp(tzu) preserves and reflects time-equation
existence by GVCREVERSESUSP20261116; therefore it also fails on the lifted
normal form with an explicitly passive operator variable.

## 4. What fails, and what is preserved

The sufficient scalar theorem MVTIMEINDICIAL20261116 and actual-height
criterion GVCREVERSEDHEIGHT20261116 remain VALID. Their hypothesis cannot
be asserted universally from individual eventuality, purity, nested ideals,
monotone heights, or even an already explicit positive linear gap.

Scalar-time finiteness of this unmodified source is also NOT invariant
under formal units. The unit-equivalent symbol Z^2 has identical heights,
all output-zero predicates and marked cutoffs, but its reversed source is

    exp(tzx+zt^2),

with first-order equation(partial_t-z(x+2t))Psi=0.
Thus the countercontrol does NOT disprove possible finiteness AFTER a
suitable height-preserving unit/phase-quotient replacement. Such a replacement
would need its own justified universal existence proof; it cannot be silently
substituted for a time equation of every actual input. The known unit and
phase-quotient closed classes remain unchanged.

The direct positive-gap obligation is invariant under these units, whereas
exact profile equations are not. The next universal argument must therefore
control the actual HEIGHT gap, or produce an actual eventual-but-nonlinear
pair, rather than retry this false unmodified-source finiteness claim.
No such gap counterexample or full positive theorem is obtained here.

## 5. A second direct attempt: geometric Chevalley estimates do not yet apply

The primary landing page of Adamus--Bierstone--Milman,
[Uniform Linear Bound in Chevalley's Lemma](https://doi.org/10.4153/CJM-2008-032-x),
describes estimates for analytic mappings regular in Gabrielov's sense.
The article is2008; the publisher's online date2018 is not a new theorem date.
Only its primary abstract/scope was checked; its full proof is NOT imported.

The natural linear map from symbols to the phase-exponential module is

    Phi(r)=exp(-sP) r(partial_x) exp(sP).

It is not an algebra homomorphism to ordinary functions. Already P=x^2 gives

    Phi(X)=2sx, Phi(X^2)=2s+4s^2x^2 != Phi(X)^2.            (9)

Consequently no geometric composition-map estimate can be applied to this
map merely by calling it a formal morphism. A genuine module/filtered estimate
would need to be proved. This observation does NOT refute all possible
Chevalley-type arguments; it records the exact missing hypothesis in this
attempt. No finite Rees model or product estimate has been supplied.

Small controls check the reversed-source coefficient identity,exact heights
and sharp marker boundary,formal-unit comparison,finite Bell leading terms,
and(9). They are transcription checks; infinite lacunarity and universal
equation nonexistence are the written arguments above. No increasing degree,
prime,rank,recurrence-fitting or coefficient sweep is performed.

Both full exceptional MZ cases remain open at their ACTUAL invariant-period
unit/relation problem. The actual four-variable finite radial lift already
has conditional linear cutoffs; this formal-source obstruction neither
settles nor replaces its mixed-isotypic purity/safety obligation. Closure-first
and deferred external review/outreach persist.
