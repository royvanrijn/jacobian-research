# Oaku graph certificates and actual polynomial-operator decision

16 November 2026. Internal written application of classical algorithms,
with an implemented rational backend and exact finite Weyl witnesses.
Independent replay, formal verification and external review remain false.
No novelty is claimed for the polynomial/convergent rational-slope mechanism
or for the generalized Bernstein algorithms. See the
[credit and scope note](BERNSTEIN_CLASSICAL_CREDIT_AND_FORMAL_EXTENSIONS.md).

For EVERY rational polynomial phase of positive degree and EVERY nonzero
zero-constant polynomial symbol, in ANY number of variables, the
[generator](algebraic_bernstein_generator.py) constructs two valid Bernstein certificates.
The [whole-source cutoff proof](WHOLE_GAUSSIAN_SOURCE_HOM_AND_EFFECTIVE_CUTOFFS.md)
explains how a checked recurrence yields numerical bounds for original markers.
It retains the operators, scalar torsion factor, explicit graph generators,
left-ideal membership witnesses and uniform specialization indices. A
separate sparse Weyl calculation checks the identities against the ORIGINAL
input-derived generators. No actual-image kernel computation is needed.

Together with the selector this gives a general pointwise algorithm for
qualitative all-marker safety and numerical sufficient affine cutoffs for
finite rational inputs. With a supplied all-positive-purity premise it gives
a linear cutoff when safe. This is not unrestricted positive GVC, not a
classification of all Gaussian safe normal forms, and not an effective
generator for arbitrary nonconvergent formal coefficient oracles.

See the [checker](verify_algebraic_bernstein_generator.py),
[historical internal-run receipt](algebraic_bernstein_generator_verification.json)
and [distributed unary certificate](algebraic_bernstein_examples/unary_homogeneous/certificate.json).
The checker contains the literal example inputs and generates its other
input, certificate, pipeline and CAS files; those generated files are not
distributed.

## 1. The explicit Fourier exponential module maps to the actual image

For P in Q[x_1,...,x_n], take the algebraic Weyl algebra on coefficient
coordinates xi_1,...,xi_n,eta. Let u be the generator of the quotient by

    xi_i-eta*P_i(-partial_xi),   i=1,...,n,
    partial_eta-P(-partial_xi).                            (1)

This is the partial Fourier transform of the rank-one exponential
connection exp(eta*P). In particular it is holonomic. Its presentation
is explicit: the exponential connection is obtained from the polynomial
structure sheaf by the Weyl automorphism partial_x_i->partial_x_i-P_i*eta,
partial_eta->partial_eta-P, followed by partial Fourier transform.

The retained actual model has

    xi_i=partial_x_i+eta*P_i,
    partial_xi_i=-x_i,
    partial_eta=actual partial_eta+P

on E=Q[x][[eta]]. Every generator in (1) kills1, so the algebraic
module maps to the actual cyclic image, with u sent to1. Polynomial
identities in this ambient algebraic module therefore descend to the
actual image. An isomorphism or full kernel computation is unnecessary
for that direction. We use neither a larger module's polynomial-Hom
dimension nor an arbitrary solution of its annihilator as a substitute
for the actual sequence.

Write f=lambda(xi), a nonzero zero-constant POLYNOMIAL. The two coefficient
functions used below are f and g=eta. The extension to an arbitrary formal
f is the separate retained MN1991 formal proof, not the present algorithm.

## 2. What Algorithms 3.6-3.8 supply, including torsion

The primary source is [Oaku, 1609.04676v1](https://arxiv.org/abs/1609.04676),
complete PDF15-18: graph presentation Theorem3.3, the torsion warning in
Lemma3.4, localization Algorithm3.6, graph criterion Proposition3.7 and
b-function Algorithm3.8. PDF3-4 also supplies holonomic specialization.

Algorithm3.6 computes the exact localized-power annihilator by FIRST
removing f-torsion via localization. We do not silently apply an older
torsion-free algorithm to arbitrary input. Instead Proposition3.7 says
that a GRAPH-module multiple of the b-function suffices, even with
f-torsion. Its scalar can be larger than the minimal localized b-function.
Algorithm3.8's graph/elimination stage computes such multiples. That is
all our actual-power estimates require.

For any module M and coefficient f, the graph module has underlying
direct sum of copies M*delta_f^(j), j>=0, with

    partial_i(v*delta^j)=(partial_i v)*delta^j-f_i*v*delta^(j+1),
    t(v*delta^j)=f*v*delta^j-j*v*delta^(j-1),
    partial_t(v*delta^j)=v*delta^(j+1), S=-partial_t*t.    (2)

Its annihilator is obtained by replacing partial_i in (1) by
partial_i+f_i*partial_t and adding t-f, as in Theorem3.3. The same
construction with two auxiliary pairs gives the explicit DOUBLE graph
of (f,eta). These are exact algebraic presentations with finite generators.

## 3. An unlocalized specialization at EVERY nonnegative integer

For k>=0 define on the graph direct sum

    Phi_k(sum_j v_j*delta_f^j)
      =sum_(0<=j<=k) (-1)^j*(k)_j*f^(k-j)*v_j,           (3)

where (k)_j=k(k-1)...(k-j+1). Terms j>k are zero; no negative power
of f is introduced. Formula (2) and the product rule prove

    Phi_k is D_base-linear,  Phi_k(S*w)=k*Phi_k(w).       (4)

For example, the extra product-rule factor (k-j) in differentiating
f^(k-j) combines with (k)_j to give (k)_(j+1), exactly matching the
second term in (2). At j=k that extra term is0; at j>k both sides
are0. Also S*delta^j=-f*delta^(j+1)+j*delta^j, whose image is
k times the image of delta^j. Thus (4) holds at the boundary indices
as well. This is a D_base[S]-map with S specialized to k; it is NOT
a fixed-exponent homomorphism for the whole auxiliary Weyl pair.

Consequently a true graph identity

    B(S)*u*delta_f=Q(S)*f*u*delta_f                       (5)

gives in the ORIGINAL unlocalized module

    B(k)*f^k*u=Q(k)*f^(k+1)*u, EVERY k>=0.               (6)

This includes f-torsion and avoids any exponent shift. No inverse f,
unspecified torsion bound, or completed-marker premise remains in (6).
It then maps to v_k=f(xi)^k*1 in the actual E-model.

Let J be the maximum partial_eta exponent in left PBW form of Q.
All coefficients preserve eta-order, partial_xi=-x preserves it,
and partial_eta lowers it by at most1. Choose k0>=1 larger than
every nonnegative integer root of B. For k>=k0, (6) gives BOTH
nonvanishing propagation v_k!=0=>v_(k+1)!=0 and the bound
0<=tau_(k+1)-tau_k<=J. These are the complete ordinary inputs
to the finite selector; they are no longer assumed metadata.

## 4. The generic first exponent and its scalar torsion witness

The second equation is constructed over Q(S), with S treated as a
CENTRAL coefficient parameter. To justify this, put

    A_f=D_base[S]*u*delta_f subset G_f(M).

The graph-to-localized-power map sends A_f to D_base[S]*u*f^S.
Its kernel becomes zero after tensoring with Q(S). Indeed a kernel
element has finitely many graph coefficients killed by powers of f
(Lemma3.4). Some power t^L therefore kills that element. Applying
partial_t^L gives

    partial_t^L*t^L=(-1)^L*S(S-1)...(S-L+1),             (7)

so a nonzero polynomial in S kills it. The generic first-exponent
module is consequently the usual generic-power module. It is
holonomic over Q(S): apply Oaku Theorem1.3 with the generic exponent
regarded as a constant in that ground field and pass to its actual
cyclic image. Zero generic modules are permitted here.

Apply the graph b-function construction to this generic module and
the coefficient eta, with second exponent Z. It gives nonzero
b(S,Z) in Q(S)[Z] and finite Q_2(S,Z) with rational S coefficients.
Clear their coefficient denominators. The resulting relation may
still hold only after scalar S localization. A NONZERO scalar
polynomial h(S) kills the remaining finite relation, giving on the
ORIGINAL double graph

    h(S)*[b(S,Z)-Q_2(S,Z)*eta]*u*delta_f*delta_eta=0.      (8)

Existence of h follows directly from the definition of central
localization. The inclusion A_f into G_f(M) stays injective under
the second graph construction, whose underlying module operation
is a direct sum of copies. Thus (8) is a relation in the explicit
full double graph, not only in an unknown actual image.

The implementation finds h by reducing successive
S^i*[b-Q_2*eta] modulo the explicit double-graph ideal and finding
the first rational linear dependence of their normal forms. The
known scalar torsion makes this terminate. It then computes LEFT
membership witnesses against the ORIGINAL graph generators.
The transcript retains both coefficient denominators and h; neither
is confused with a shift of the actual power index.

## 5. The polynomial constraint and specialization thresholds

Apply Phi_k for the first graph in (8), for any k>=0. The second
graph remains. Since the actual module embeds into E, eta is
injective and its localization embeds into E[eta^-1]. On this
second graph the usual falling-factorial specialization is defined
for EVERY integer j after eta localization. Therefore (8) gives

    h(k)*b(k,j)*eta^j*v_k
      =h(k)*Q_2(k,j)*eta^(j+1)*v_k.                      (9)

For finite tau_k choose j=-tau_k-1. The right operates on
eta^-tau_k*v_k in E and stays in E. The left has a possible simple
eta-pole with nonzero leading polynomial coefficient. Its scalar
must vanish. Thus

    p(k,tau_k)=0, p(X,Y)=h(X)*b(X,-Y-1)!=0               (10)

at EVERY nonnegative k with finite tau_k. Vertical factors h(X)
may make some indices uninformative; the selector already moves
beyond the roots of the leading Y coefficient. We do not divide
them out or falsely infer minimality of b.

The retained thresholds are explicit:

* ordinary graph relation and order constraint: valid from k=0;
* nonvanishing/step bound: k0=1+max({0} union nonnegative integer roots(B));
* later selector index K: computed from p,k0,J by root-free collars.

Thus all hypotheses of the implemented selector are generated from
finite original input and exact identities. Its finite original
T^K(P^j) evaluations then determine the actual slope/nilalternative
and produce numerical sufficient affine or pure-linear cutoffs.

The complete --decide operation also handles zero phase, zero operator,
nonzero symbol constant and constant phase directly. Zero phase/operator
gives identically zero positive outputs. A nonzero operator constant makes
T invertible on every bounded polynomial-degree space, so a nonzero phase
is unsafe already at marker1. A nonzero constant phase and positive-order
T are safe for m>deg(Q). These elementary cases require no graph search.

## 6. Complete finite certificate verification

The backend uses the installed Singular graph/elimination machinery,
with a fixed random seed. Older SDLoc is used as GRAPH-annihilator
elimination, not as an unsupported assertion of exact f-saturated
localized annihilators. Returned identities are checked independently
of that interpretation against the explicit double-graph generators.

In canonical coordinate-then-derivative PBW order the sparse verifier
uses the exact formula

    (x^a D^b)(x^c D^d)
      =sum_(l<=min(b,c)) binom(b,l)*(c)_l
                         x^(a+c-l) D^(b+d-l).            (11)

It reconstructs the graph generators from P and lambda and checks
both exported identities as literal sums of LEFT witness multiples.
For the ordinary identity it first lifts partial_eta to
partial_eta+partial_(t_eta), as required when adding the second
graph. This lift is essential and is checked in the recorded data.
No claimed normal-form-zero flag is trusted in place of the
coefficient-by-coefficient rational identity (11).

The verifier also recomputes the scalar roots, J, h*b specialization
and selector summary. Altering an operator, witness, graph generator
or summary therefore fails. These finite algebraic certificates are
stronger evidence than a bounded moment prefix; the universal
Fourier, graph and specialization proofs remain written arguments,
not a formal theorem-prover development or external review.

With no time cap, the stated exact algorithms terminate for finite
polynomial inputs. The CLI defaults to a per-stage time cap and
returns UNKNOWN on timeout. No classification may be inferred from
a timeout, missing witness or failed exact check. There is no
polynomial complexity or practical all-degree performance claim.

## 7. Example and executable use

For P=x^2 and lambda=xi^3 the generated ordinary polynomial is

    B(S)=(S+1)(3S+1)(3S+2)(3S+4), J=2, k0=1.

The second polynomial is

    b(S,Z)=(2Z+3S+2)(2Z+3S+3), h(S)=1.

Equation(10) becomes (3k-2tau_k)(3k-2tau_k+1)=0,
and the original-power query selects rho=3/2. The complete operators
and graph witnesses are retained in the linked example, not omitted
after reading off these polynomials. In the inactive-direction case
P=x,lambda=zeta, the backend instead finds h(S)=S. This genuinely
tests scalar torsion and the eventual-infinity branch.

Input JSON contains dimension, phase_terms and symbol_terms. Each
term is [exponent_vector,"rational coefficient"]. From the repository root,
supply `input.json`; it is not a bundled demonstration file. For the
example above its contents are
`{"dimension":1,"phase_terms":[[[2],"1"]],"symbol_terms":[[[3],"1"]]}`.
The commands generate their output files.

    .venv/bin/python research/moments/vanishing/algebraic_bernstein_generator.py input.json certificate.json
    .venv/bin/python research/moments/vanishing/algebraic_bernstein_generator.py input.json result.json --decide
    .venv/bin/python research/moments/vanishing/algebraic_bernstein_generator.py certificate.json verified.json --verify

The output includes both full operators, explicit graph generators,
both lists of left witnesses, cleared generic denominators, h(S),
specialization thresholds and the finite selector summary. The --decide
operation runs generation, direct verification, the selector's original-
power evaluations and cutoff extraction. The --verify operation replays
only rational Weyl arithmetic, without calling Singular. The checker
exercises the pipeline and CLI. The backend itself does not assume or
decide all-positive purity; --pure supplies that premise when a linear
cutoff is requested. A time limit of0 disables the computational cap.

## 8. Attribution and scope

Sabbah and Maisonobe receive credit for the classical analytic/convergent
Bernstein rational-slope mechanism, and Oaku for the generalized algorithms.
This polynomial implementation applies those inputs; it does not establish
novelty or independent verification.

The construction extends mathematically over effectively presented finitely
generated characteristic-zero fields with the required exact operations.
The executable backend here is over Q. No algorithm recognizing polynomial
representatives for arbitrary formal input, all-ring rational-slope
generation, least cutoff or complete Gaussian normal-form classification
is asserted.
