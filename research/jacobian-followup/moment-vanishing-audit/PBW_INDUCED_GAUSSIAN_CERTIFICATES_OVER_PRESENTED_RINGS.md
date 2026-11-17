# PBW-induced Gaussian certificates over presented coefficient rings

16 November2026. Internal universal proof and implemented exact backend.
Independent replay, formal verification and external review remain false.
Priority remains unassessed. This closes the finite-primary certificate
implementation gate in the sealed [first-marker and real-region theorem](GAUSSIAN_FIRST_MARKER_SLOPES_AND_SHARP_REAL_CUTOFF_POLYHEDRA.md).
It does not classify every safe phase or actual Gaussian connection.
See the [producer](gaussian_primary_real_cutoff_backend.py),
[checker](verify_gaussian_primary_real_cutoffs.py),
[receipt](gaussian_primary_real_cutoff_verification.json) and
[source audit](gaussian_primary_real_cutoff_source_audit.json).

## 1. A faithful finite-algebra source without splitting fields

Let F be a characteristic-zero field and T=F[Y_1,...,Y_q]/I a nonzero
finite-dimensional commutative algebra. Neither reducedness nor a split
residue field is assumed. Form the LEFT Weyl quotient

    N_T=D_Y/(D_Y I)=D_Y tensor_(F[Y]) T.

Normal ordering with derivatives on the LEFT makes D_Y a free RIGHT
F[Y]-module with basis D_Y^alpha. Consequently

    N_T = direct_sum_(alpha>=0) D_Y^alpha tensor T,
    T --> N_T, a --> 1 tensor a, is injective.              (1)

This assertion uses the right-module PBW basis, not the false assertion
that I kills every vector of N_T. In the explicit derivative-polynomial
model F[zeta] tensor T, the Weyl actions are

    D_Yi = multiplication by zeta_i,
    Y_i  = M_Yi - partial_(zeta_i),                         (2)

where M_Yi is the exact multiplication table of Y_i on T. Thus
[D_Yi,Y_j]=delta_ij. Every polynomial in I kills the seed 1 tensor1;
all ordinary nilpotents survive (1). This is an infinite-dimensional
Weyl module, not an impossible finite-dimensional Weyl representation
and not a purported derivation on T.

For the order filtration, derivative-left normal ordering identifies
gr(D_Y I) with I*F[Y,zeta]. Hence

    Char(N_T)=V(I) times A^q_zeta,
    dim Char(N_T)=q.

For q=0 this says N_T=F. Thus N_T is holonomic. The PBW and characteristic-
variety conventions are classical; see [Sattelberger--Sturmfels,
1910.01395v1](https://arxiv.org/html/1910.01395v1), Section1, normal
monomials, Definition1.8, Theorem1.10 and Definition1.11. Equation(1)
and this short induced-module argument supply the explicit presentation
used here; no novelty is claimed for PBW induction or holonomic point
modules. In particular no splitting field, Macaulay inverse-system basis,
or common compositum is needed by this implementation.

## 2. Explicit Gaussian graph and both scalar certificates

Take a polynomial lift P(Y,X) of a T-valued phase, and fixed nonsingular
rational symmetric Gaussian precision H. Tensor N_T with the Gaussian
exponential module in X. Add the graph of P, then Fourier transform ONLY
X by X_i->D_ui, D_Xi->-u_i. The ORIGINAL input-derived generators are

    G_i=-u_i+sum_j H_ij D_uj+P_Xi(Y,D_u)D_v, i=1,...,n,
    I_j(Y), for every retained localized primary relation,
    G_0=v-P(Y,D_u),
    S=-D_v v, E=sum_i u_i D_ui.                            (3)

There is no omitted D_Y relation: the induced module has only I_j(Y)
as initial coefficient relations. All Y,D_Y have u-order weight0.
This is holonomic by tensor product, graph direct image and Fourier
transform. The exact unlocalized falling-factorial specialization from
the [Oaku packet](OAKU_GRAPH_CERTIFICATES_AND_ACTUAL_POLYNOMIAL_DECISION.md)
maps its seed at EVERY m>=0 to the original marker-generating series
V_m(u) with coefficients embedded by (1). Its u-order is exactly the
least original detecting-marker degree in T, including nilpotents.

The earlier Oaku and global weighted b-function arguments now apply
directly to (3), over F and generically over F(S). They give

    B(S)-Q(S)P(Y,D_u) = sum LEFT witnesses*generators,
    p(S,E)-R(S)       = sum LEFT witnesses*generators,       (4)

with nonzero scalar polynomials B in F[S], p in F[S,Z]. Every PBW
monomial of R has |u|-|D_u|>=1. Coefficient denominators and finite
central S-torsion are absorbed into the EXACT unlocalized identities;
no scalar factor is silently divided away and no minimal localized
b-function is asserted. Retain both complete operators and all LEFT
witnesses. Both identities specialize from m=0. Beyond the nonnegative
integer roots of B,

    -deg_X(P) <= d_(m+1)-d_m <= max total D_u degree(Q).

The Euler initial term gives p(m,d_m)=0 for every finite d_m, including
the exceptional indices at which p(m,Z) may vanish identically. No
actual-image kernel, equality of arbitrary solution spaces, or regularity
hypothesis enters this direction. Oaku supplies the ordinary graph
multiple with torsion; Andres et al.,1005.3257v1 Section5 supplies the
global weighted Euler b-function. These classical imports and the
separate formal/convergent credit boundaries stay as previously sealed.

## 3. General finite-input algorithm and scalar descent

Input is the retained family schema: R=QQ[c_0,...,c_(r-1)]/J, a finite
polynomial phase in R[X], and fixed rational precision H. The zero ring
is handled directly. The existing ordinary producer computes the FULL
primary decomposition J=intersection Q_j, including embedded components,
checks the intersection and prime radicals in Macaulay2, and chooses
maximal independent parameter sets U_j. It supplies faithful embeddings

    R --> product_j T_j,
    T_j=QQ(U_j)[Y_j]/Q_j, dim_(QQ(U_j)) T_j finite.

Primary decomposition remains explicitly TRUSTED CAS. The exact finite
multiplication tables, standard monomials, original relations, independent
sets and CAS transcripts are retained. Auxiliary ray parameters used by
the old table implementation are stripped from the new graph coefficient
field; no table or phase coefficient depends on them.

For each component of length r_j, test P^r_j=0 in its finite multiplication
table. This is equivalent to P being nilpotent: the radical of a finite
algebra of length r_j has r_j-th power zero. In that case all relevant
powers are finite, and exact heat orders give the region directly.
Otherwise use (3),(4). The backend implements exhaustive increasing
finite boxes of scalar polynomials, base operators and LEFT graph
multipliers over QQ(U_j). Each box is exact sparse linear algebra.
Existence of (4) proves uncapped termination in EVERY dimension/degree.
A cap returns UNKNOWN with retained completed components. There is no
complexity or practical universal-performance claim; budget checks are
cooperative outside the separately timed primary CAS invocation.

To reuse the rational collar arithmetic, clear every QQ(U_j) coefficient
denominator in each scalar constraint and expand the numerator in U_j.
Keep one NONZERO QQ component, together with the denominator, full cleared
numerator and selected parameter monomial. Vanishing in QQ(U_j) at
integer (m,d_m) forces every such component to vanish. Excluding the
integer roots of a chosen ordinary component is conservative and valid.
The resulting rational constraints and bounds feed the sealed arithmetic
on a_m=d_m+D*m; unshifted orders need not be monotone. Exact original heat
queries select rho and the finite exact late parallel-line list.

The sharp component regions combine by INTERSECTION. At every power the
full-ring order is the MINIMUM of component orders. Thus the global
asymptotic edge is the maximum of component edges; any unsafe component
makes the full ring unsafe; nilpotent components retain finite intercept
inequalities. At each retained supporting power the implementation
recomputes the full original-ring heat order, preserving the sharp region.
The E=0 slice is reported with its attained rational minimum when nonempty.

This closes the executable backend for finite PRESENTATIONS. Over an
arbitrary Q-algebra, the finite coefficient subalgebra gives the existence
theorem, but an unknown coefficient kernel is not an effective input.
No arbitrary nonconvergent-formal coefficient-oracle algorithm is implied.

## 4. Two exact optimizations with universal proofs

For a certified positive law d_m=rho*m+c_i at every m>=K, and C>=1/rho,
the deficit m-C*d_m is nonincreasing along each candidate line. Check all
early powers and the FIRST ADMISSIBLE integer point on EVERY candidate
line (integer d>=0). For rho=a/b reduced, c_i=c/b, solve

    a*m+c=0 mod b,
    m>=max(1,K,ceil(-c_i/rho)).

If all these bounds are <=E, the vertex is proved valid. A failed envelope
is INCONCLUSIVE because its line might never actually occur. The backend
then calls the sealed exact rational cutoff oracle, retaining an ORIGINAL
nonzero moment on failure. Finite vertex refinement and strict decreasing
first-index descent on finitely many lines terminate as already proved.
No periodicity or recurrence oracle is needed.

For a zero-constant spatial-linear phase P=l(Y).X, suppose some l_i is a
unit in T, with inverse q_i represented polynomially, and variance
v=l^T H^-1 l is nilpotent of index ell. Retain ordinary relation witnesses

    q_i*l_i-1=sum a_j I_j,  v^ell=sum b_j I_j.

Writing f=P(Y,D_u), T_i=-u_i+sum H_ij D_uj gives the exact ordinary identity

    S+1-q_i*T_i*f
      =-q_i*f*G_i-D_v*G_0+(q_i*l_i-1)*(f*D_v-1).        (5)

Let q(u)=u^T H^-1 u, A=E-q(u)-S, V=v*D_v^2. Then

    A+V=sum_i ((H^-1*u)_i+(H^-1*l)_i D_v)G_i+D_v G_0,
    [A,V]=-2V.

Induction gives F_ell(A)=product_(k=0..ell-1)(A+2k) modulo the original
graph ideal as (-1)^ell V^ell. Its latter term is killed by the b_j
relations. The code retains the full recursive LEFT witnesses. Set

    p(S,E)=product_(k=0..ell-1)(E-S+2k),
    R=p(S,E)-F_ell(E-q(u)-S).

Every monomial of R raises u-order by at least2: replacing E-S by
E-S-q(u) inserts a raising factor, while E preserves its degree. Thus
this is an exact shortcut for (4), for arbitrary length, with a general
search fallback outside its hypotheses. These identities are verified
by literal sparse Weyl multiplication; ordinary table membership alone
is never substituted for a graph witness.

## 5. Focused controls and remaining frontier

Six end-to-end controls retain all original presentations and transcripts:
free aZ over QQ[a]; Z+eW over dual numbers; Z+eW over
QQ[a,e]/(e^2,ae) with BOTH primary components; yZ over QQ[y]/(y^2-2)
without splitting the residue field; nilpotent eZ; and nonlinear aZ^2
through the general QQ(a) graph search. The checker also exercises a
zero ring, a zero phase and an unsafe real phase, scalar projection,
fractional line congruences, Weyl/induced-module actions and corruptions.
The exact proof of generality is Sections1--4, not these finite controls.

For Z+eW with E[ZW]=1 and e^2=0,

    heat(P^m)=Z^m+m e Z^(m-1)W+m(m-1)e Z^(m-2),
    d_1=1, d_m=m-2 for every m>=2.

The sharp region is C>=1,E>=2. The embedded component produces precisely
the same E>=2 obstruction; deleting it would return the WRONG E>=0.
For aZ the region is C>=1,E>=0 on the entire presented ring, including
a=0; for aZ^2 it is C>=1/2,E>=0. Nilpotent eZ has C+E>=1 in the
nonnegative quadrant, with asymptotic coefficient0 and pure minimum1.

Sabbah/Maisonobe receive the convergent rational-slope credit, Oaku the
generalized graph algorithms, and the cited sources their PBW/weighted
b-function machinery. The arbitrary nonconvergent formal extension and
all-ring descent remain separately scoped, with priority unassessed.
All historical Oaku operators, torsion factors, thresholds and sealed
receipts are preserved. Old capped cusp/cone/rational-ideal UNKNOWN jobs
are not restarted. The remaining broad goal is exhaustive safe/actual-
connection structure and regular transcendental realization outside the
known architectures. Closing this backend does not mark that goal complete.
