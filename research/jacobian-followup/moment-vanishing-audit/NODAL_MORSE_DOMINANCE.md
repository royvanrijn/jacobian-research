# Completing the nodal critical-line class and excluding the finite Morse locus

14 November 2026. Internal written analytic proofs and exact auxiliary
algebra. The checker does not certify the analysis. No external review,
independent specialist replay, formal verification or priority is claimed.

The subsequent [full formal-frame theorem](NODAL_FORMAL_FRAME.md)
closes every irreducible nodal-leading cubic, including all finite
degeneracies, coincident critical values and quartic resonances. Earlier
method boundaries below describe the scope of this retained argument.

## Results

Keep the [nodal affine normal form and full Gaussian exponent](NODAL_GAUSSIAN_CRITICAL_PROFILES.md):

    P=T*(uv+h)+a*u^3+d*v^3+b*u^2+c*v^2+e*u+f*v+k,
    a*d!=0, m=Q_TT!=0, det Hess Q!=0.                  (1)

All metric cross terms and Gaussian means are retained.

**Complete critical-line exclusion (G3NODALCENTERALL20261114).** Every
member with h=e=f=0 has a nonzero positive pure Gaussian moment, for all
k,b,c and affine Q. This removes the strict-triangle restriction from the
[preceding center-line theorem](NODAL_CENTER_LINE_DECAY.md). In particular
every irreducible nodal-leading ternary cubic with a positive-dimensional
affine critical set is pure-nonvanishing: the remaining nonisotropic-node
chart is (1), and the isotropic/coercive alternatives were already excluded.

**Finite Morse exclusion (G3NODALMORSE20261114).** Suppose the affine
critical set of P is finite, all its points are Morse, and their P-values
are pairwise distinct. Then some positive pure moment is nonzero, subject
to one additional condition on the quartic escape profile:

    h=f=0, e*c!=0: gamma=d*e*m/(2*c)-l_T !=0.          (2)

Here l_T is the linear T coefficient of Q. If e=0,f!=0, first exchange
u,v and all corresponding coefficients. On the other three escape profiles
there is no additional gamma condition. These conditions define a nonempty
open subset of each finite-critical profile; they do not cover its
degenerate, coincident-value or quartic-resonance boundaries.

The subsequent [resonant-pair theorem](NODAL_QUARTIC_RESONANCE.md)
also excludes gamma=0 when Xi!=0, with Xi given there. Thus the remaining
quartic resonance on this Morse/distinct-value locus satisfies both
gamma=0 and Xi=0. The earlier theorem below retains its original scope.

Neither theorem gives a numerical uniform pure prefix. General cubics
with cuspidal/reducible leading curve and the remaining nodal exceptional
strata are not closed by this note.

## 1. Local Morse periods and the actual global cut

The previous center-line proof supplies a flat cut basis of total rank7
and parameter-independent normalized tail bounds on h=e=f=0. If p is an
isolated Morse critical point of P with value alpha, it has one nearby
S_t=tP-Q critical point p(t), with

    p(t)=p+t^(-1)*(Hess P(p))^(-1)*grad Q(p)+O(t^(-2)),
    c(t)=S_t(p(t))=alpha*t-Q(p)
          -(1/(2t))*grad Q(p)^T*(Hess P(p))^(-1)*grad Q(p)
          +O(t^(-2)).                                (3)

The holomorphic Morse lemma with parameter epsilon=1/t applies to
P-epsilon*Q near p. The local thimble integral of the constant volume form
has nonzero leading term

    J_p(t)=exp(c(t))*t^(-3/2)*(C_p+o(1)), C_p!=0.      (4)

The determinant is nonzero, and the volume amplitude at p is1. The
orientation and square-root branch change C_p, not its nonvanishing.
The local result is [Pham, PartII Section1.1, formulas(1.5)'--(1.7),
printed page26](https://numdam.org/item/AST_1985__130__11_0.pdf).

We need (4) for the full cut period, not just a local integral. Truncate
its normalized phase cut at L=K*log|t|. For large t this corresponds to a
P-Q/t value disk of radius O(log|t|/|t|), inside the local Morse chart.
The cut's fiber class is the local vanishing generator, of rank one.
By the Gelfand--Leray slicing formula, its integral over this initial cut
equals that computed with the local representative: fiberwise integrals
of omega/dS_t depend only on that class. A spatially wandering representative
does not change the sliced period. The uniform tail bound makes the rest
o(exp(c(t))*|t|^(-3/2)) when K is sufficiently large. Thus (4) holds for
the actual full flat cut cycle. No assertion of formal exponential splitting
is used in this matching.

## 2. Two Morse points close the entire strict triangle

On h=e=f=0 the nonconstant critical values of P are

    alpha0=k, alpha1=k+4*b^3/(27*a^2),
    alpha2=k+4*c^3/(27*d^2).                          (5)

If0 is outside or on the boundary of their convex hull, the preceding
half-plane exclusion applies. Otherwise they are distinct vertices of a
nondegenerate triangle strictly containing0. In particular b,c!=0.
The two off-axis critical points are

    p1=(-2b/(3a),0,0), p2=(0,-2c/(3d),0),
    det Hess P(p1)=8*b^3/(9*a^2)!=0,
    det Hess P(p2)=8*c^3/(9*d^2)!=0.                  (6)

Their local Milnor numbers are one each. The other five dimensions of
phase-critical homology tend to the critical line, with critical actions
k*t+O(1). The preceding cut bounds show that every one of these five
scalar periods is O(exp(Re(k*t))) on each fixed ray.

There is a ray t=xi*R, R>0, such that, after interchanging the off-axis labels,

    Re(xi*alpha1)>Re(xi*alpha2)>0>Re(xi*k).          (7)

Here is an explicit construction, so the triangle is not silently replaced
by a half-plane. Write0=lambda0*k+lambda1*alpha1+lambda2*alpha2 with all
lambda_i>0 and sum1. Put delta=alpha2-alpha1 and
cross=Im(conj(alpha1)*alpha2)!=0. For sigma=sign(cross), let

    n=sigma*i*conj(delta), q=abs(cross)>0.

Then Re(n*alpha1)=Re(n*alpha2)=q and
Re(n*k)=-(lambda1+lambda2)*q/lambda0<0. Perturb n by a sufficiently small
positive epsilon*conj(delta). Both positive signs and the negative sign
persist, while the two off-axis real parts differ by epsilon*abs(delta)^2.
Normalize the direction if desired; its length is immaterial.

Choose independent decaying cut directions as in the preceding theorem.
The full rank7 basis consists of the two off-axis rank1 generators and
five line generators. If a constant scalar period were expressed in this
basis, divide by the larger growing expression in (4). The first coefficient
must be zero because C_p!=0 and every other term is exponentially smaller.
Repeat for the second off-axis coefficient. The remaining five terms decay
exponentially by (7), contradicting a nonzero constant. The retained
formal-seed/actual-period bridge now excludes all-zero pure moments.

This proves the complete critical-line theorem. The former method-boundary
example has

    (k,alpha1,alpha2)=(-1-i,2-i,-1+2i),
    xi=(27-33i)/10,
    (Re(xi*k),Re(xi*alpha1),Re(xi*alpha2))=(-6,21/10,39/10).

It is now a positive certificate for ray separation, not a surviving
triangle obstruction and not a pure-null witness.

## 3. Polynomial carriers for the finite-critical profiles

For each fixed nonzero t, S_t has a complex linear gradient lower bound
outside a compact set. Its cubic leading form has one ordinary node and
its quadratic value there is -m/2!=0, so the retained
[coercive-gradient criterion](COERCIVE_GRADIENT_MOMENTS.md) applies.
Unlike the center-line case, the constants and radius cannot be fixed
independently of t: the known critical branches escape.

Nevertheless, along any fixed ray t=xi*R they can be bounded by powers of R.
Indeed, consider the set of triples (R,rho,delta), with rho>=1 and0<delta<=1,
satisfying

    for all x, ||x||>=rho implies
      ||grad S_(xi*R)(x)||>=delta*||x||.             (8)

It is semialgebraic by real quantifier elimination and has a nonempty fiber
over every sufficiently large R. A cylindrical cell decomposition chooses
rho(R),delta(R) semialgebraically on an unbounded interval. Their graphs
are algebraic branches at infinity, so Cauchy root bounds give

    rho(R)+delta(R)^(-1)<=C*R^B.                   (9)

No explicit B is claimed. This is a one-parameter existence argument, not
an executed many-variable quantifier-elimination campaign.

Gradient transport from local critical neighborhoods now stays in

    ||x||^2<=C*R^D*(1+s).                         (10)

All critical points start inside rho(R); after enlarging it by one,
d||x||^2/ds<=2/delta(R) outside that ball proves (10).
Use the compactified constrained cut pairs and compatible Hardt triviality
from the center-line proof with this polynomial carrier. The constraint
still forbids spatial infinity at a finite cut length. It gives a jointly
semialgebraic flat basis with fixed multiplicities on an unbounded R-interval.
The same coordinate-minor area bound now yields

    |exp(-c_i(t))*tail_(s>=L)|
      <=C*R^A*(1+L)^M*exp(-L/2).                 (11)

Thus L=K*log R makes every normalized tail smaller than any assigned inverse
power. The proof establishes a joint parameter/cut-length bound, not only
fixed-t convergence. All seven generators are retained, because their
finite and escaping multiplicities exhaust the tame rank7 homology.

## 4. The seven nonzero leading amplitudes

If the finite P-critical points are Morse, each finite cluster has rank1
and (4) applies, using (11) for the global matching. Its critical action
is alpha*t+O(1), and its scalar prefactor is t^(-3/2) times a nonzero constant.

The [four escape profiles](NODAL_GAUSSIAN_CRITICAL_PROFILES.md) give n=1,3,4,5
escaping points. Their full-cut scalar periods have

    J_i(t)=exp(c_i(t))*t^(-1-1/n)*(C_i+o(1)), C_i!=0. (12)

For n=j+2 in{3,4,5}, write r=t^(-1/n), u=r^j*U, v=r*V, T=K/r.
After removing k*t the leading phase is r^(-2) times

    K*U*V+e*U+g_j*V^j-m*K^2/2.

Its critical Hessian determinant is n*m*K^2!=0. The volume factor r^(n-2)
and the three-variable Morse scale r^3 give t^(-1-1/n). The normalized
cut of logarithmic length stays in this chart, since its coordinate radius
is O(r*sqrt(log R)). Higher terms are uniformly negligible there. Equation(11)
then transfers the local Gaussian leading term to the full cut.

For n=1 use u=U/t, v=V/t, T=(h/m)*t+W. After removing the divergent scalar
action, the limiting phase is a nondegenerate quadratic in U,V,W:

    (h/m)*UV+(e-D*h/m)*U+(f-E*h/m)*V
       -m*W^2/2-l_T*W-q0.

Its volume factor is t^(-2), and all cubic corrections tend to zero on the
logarithmic Gaussian region. Again its Gaussian leading coefficient is
nonzero, proving (12). This check avoids assuming that an original-coordinate
Hessian determinant alone controls all higher terms.

## 5. Distinct positive actions and the quartic exception

For n=1 the escape action has a nonzero t^2 term, hence differs from every
finite action. For n=3,5 the escape actions have the common k*t plus
nonzero multiples of K^2*t^(2/n), with K running through n distinct roots.
Because n is odd, their K^2 are distinct. Finite critical values were
assumed pairwise distinct, and a finite action cannot equal an escape
action's positive part, even if its linear coefficient equals k.

For n=4, opposite K have the same t^(1/2) coefficient. One more term is needed:

    c_K(t)=k*t-m*K^2*t^(1/2)
              +(d*e*m/(2*c)-l_T)*K*t^(1/4)+O(1),
    K^4=-2*c*e^2/m.                               (13)

The coefficient follows from the first correction d*V^3-l_T*K in the
scaled phase, using V=-e/K. Corrections to the critical coordinates do not
alter it because the leading phase is stationary. Condition(2) distinguishes
opposite roots and hence every positive-action polynomial. If gamma=0,
the positive parts of opposite-root actions coincide; this case is retained
as an unresolved boundary, not treated as generically distinct.

Choose a t-ray avoiding the finitely many directions on which the first
nonzero positive-power term of an action, or of a difference of two distinct
actions, has zero real part. On this ray, all growing exponentials have
strict asymptotic growth order. Repeatedly divide a putative constant
period by the largest one in (4)/(12). Its nonzero leading scalar amplitude
forces that coefficient to vanish. After all growing contributions are
removed, every negative-action contribution decays exponentially. A bounded
action can only come from the at-most-one finite critical point with
alpha=0; its prefactor t^(-3/2) also decays. This contradicts a nonzero constant.
The seed-submodule bridge proves the finite Morse theorem.

This is successive dominance for an actual flat cut basis with verified
nonzero leading amplitudes. It is not a claim that arbitrary multiple
critical blocks have already been formally diagonalized.

## Evidence, sources and the new frontier

The [checker](verify_nodal_morse_dominance.py) and
[receipt](nodal_morse_dominance_verification.json) verify the universal
triangle identities, exact separating ray, both off-axis Hessians,
first critical displacement/action correction, the quartic term(13), and
squarefree finite critical-value polynomials at one point in each profile.
For h!=0, the example

    P=T*(uv+1)+u^3+v^3+u^2+2v^2+u+3v

has a squarefree degree6 critical polynomial and a squarefree degree6
critical-value resultant. For h=0, the examples (e,f,c)=(1,2,0),(1,0,1),
(1,0,0), with a=d=m=1,b=l_T=k=0, have finite value polynomials

    (27s^2+4)*(27s^2+32), (27s^2+4)*(27s-4), 27s^2+4.

These demonstrate nonempty open loci; they are not substituted for the
all-parameter written proofs. The initial checker used an inferred coefficient
field missing d; the rejected input and failure receipt are retained, and
the corrected common-field reduction passes.

The [source audit](nodal_morse_dominance_source_audit.json) records the new
Pham Morse passages and the retained cut estimates. It also distinguishes
the indexed abstract of Forjan's June2026 SSRN paper from its inaccessible
full text, and records why general formal direct-image claims were not
imported without their infinity hypotheses. No new formal stationary-phase
theorem is needed for either exclusion here.

Among irreducible nodal-leading cubics, any remaining possible pure-null
case must have a finite affine critical set and lie on at least one of:
a degenerate P-critical-point locus, a coincident P-critical-value locus,
or the quartic gamma=0 resonance locus (including the swapped version).
The center-line triangle is closed. The remaining finite exceptions need
grouped local analysis or a justified formal comparison; no result is
inferred from continuity across those discriminants.

General nonhomogeneous cuspidal/reducible-leading cubics, numerical uniform
pure gates, coefficient schemes and lower-power GVC remain open. All analytic
steps here remain internal written proofs with the named imports; the
symbolic checker explicitly does not certify them.
