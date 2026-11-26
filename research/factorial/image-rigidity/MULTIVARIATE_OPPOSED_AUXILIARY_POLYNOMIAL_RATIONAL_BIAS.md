# Auxiliary polynomials at opposed rational-bias factorial contacts

Fix an ordinary two-monomial core, adjacent orders \(n,N=n+1\)
with \(n\ge2\), fixed positive gaps \(d,e\), and fixed signed
gap \(h=dr_+-er_-\). Append \(r_+\) copies of \((L,L+d)\)
and \(r_-\) copies of \((Y+e,Y)\). The heights and counts are
ordinary integers. Use the exact normalization and finite local
branches from the
[opposed diffuse branch theorem](MULTIVARIATE_OPPOSED_DIFFUSE_PUISEUX_HEIGHT_BRANCHES.md).
Fix a positive physical Gaussian contact
\((T_*,\lambda_*)\) with \(\lambda_*/T_*\in\mathbb Q\),
allowing a vanishing limiting curvature share.

**Theorem.** Any local branch with infinitely many ordinary
collisions approaching the contact is contained in a rational
algebraic curve in the two height coordinates. In particular, a
branch with a transcendental Puiseux coefficient of the
height ratio in the reciprocal active height contains only
finitely many sufficiently high ordinary collision profiles.
This excludes all branches continuing the specified root pair
in the [stationary rational-bias criterion](MULTIVARIATE_OPPOSED_DIFFUSE_PUISEUX_HEIGHT_BRANCHES.md#a-stationary-contact-of-rational-bias).
The stationary exclusion also holds on the boundary under
the same simple-root and nonzero forcing assumptions, as
shown below.
The rational-bias hypothesis supplies the one-parameter
coefficient field used in the dimension estimate below. No simplicity, transversality, or affirmative
test for a particular valuation is required for the algebraic
height-curve conclusion.

The proof uses an auxiliary polynomial in the two heights and
the common count parameter. The number of its coefficients grows
cubically with its degree. The possible Taylor coefficients of
its substitution lie in a function field of one variable with
linear pole growth, so the number of rational conditions to
impose a high order of vanishing grows quadratically. Integer
evaluation then forces an algebraic relation. Ax's differential
theorem transfers that relation to the heights alone.

## One transcendental coefficient parameter

As in the branch theorem, write \(r_\pm=c_\pm r+s_\pm\),
\(u=1/r\), \(p=r/L\), and \(q=r/Y\), with fixed integer
\(c_\pm,s_\pm\). Put \(\theta=T_*>0\). Write the limiting
\(p_*,q_*\) as \(\kappa_+\theta,\kappa_-\theta\) with rational
\(\kappa_\pm\). At least one \(\kappa_\pm\) is positive.
The clock \(w_*=\exp(-\theta/(2nN))\) is algebraic and
different from one. In particular \(\theta\) is transcendental
by Hermite–Lindemann, and \(\exp(a\theta)\) is algebraic for
every rational \(a\).

Let \(\Theta\) be an indeterminate. In the coefficient
expansion of the normalized parents, regard \(\Theta\) as
the symbol for \(\theta\) except in the fixed exponentials
of rational multiples of \(\theta\). Those are fixed algebraic
constants. Center the variables at \(p_*,q_*\), and denote
the increments by \(P,Q\). Each summand in a coefficient
logarithm is of the form
\[
     (c/u+s)\log(1+\ell u(\kappa\theta+P)),
\]
with rational \(c,s,\ell,\kappa\), or with \(Q\) in place of
\(P\). Its constant term in \(u\) is
\(c\ell(\kappa\theta+P)\). For \(a\ge1\), its coefficient of
\(u^a\) is polynomial in \(\Theta,P\) of \(\Theta\)-degree
at most \(a+1\le2a\). Thus, after exponentiation, every
coefficient of \(u^a P^b Q^c\) is a polynomial in \(\Theta\)
of degree at most \(2a\) over a fixed number field \(K\).
Indeed the exponential of the constant in \(u,P,Q\) is a
fixed algebraic number, and multiplication adds the powers
of \(u\) and the \(\Theta\)-degree bounds. The same is true
for the reciprocal of any such series with nonzero algebraic
constant term. This suffices for the off-diagonal quotient
resultant and the physical resultant. In passing to centered
logarithmic clock variables \(W=\log(w/w_*)\) and
\(V=\log(v/v_*)\), the change from \(P,Q\) is an invertible
rational linear map.

The degree bound can be expressed by a lift to
\[
                            K[[U_0,U_1,U_2,W,V]]
\]
with specialization \(U_i=u\Theta^i\). Each monomial
\(u^a\Theta^b\) with \(0\le b\le2a\) has a lift of total
\(U\)-degree \(a\): write \(b=2j+\epsilon\) with
\(\epsilon\in\{0,1\}\) and take
\(U_0^{a-j-\epsilon} U_1^\epsilon U_2^j\).
A formal power series lift therefore exists for both resultants.

Apply formal Weierstrass preparation to the lift of the angular
resultant in \(W\). By uniqueness, at \(U_0=U_1=U_2=0\)
its monic factor is \(W^a\), where \(a\) is the order of
vanishing of the central angular eliminant. Divide the lift
of the physical resultant by that monic polynomial in \(W\)
and take the norm of its remainder in the monic quotient.
Here the norm is the determinant of the matrix for multiplication
on the free quotient module. Upon setting all \(U_i\) to zero,
that matrix is triangular with diagonal equal to the physical
resultant at \(W=0\). The norm is thus regular in \(V\),
by the central isolated-fiber argument in the branch theorem.
A second formal preparation in \(V\) yields a monic
polynomial in \(V\) with coefficients in \(K[[U_0,U_1,U_2]]\).
Specializing \(U_i=u\Theta^i\) preserves formal preparation
and the norm identities. At \(\Theta=\theta\) these formal
preparations are the Taylor series of the analytic
preparations, by uniqueness. This keeps the field and
the coefficient-degree bound through both eliminations.

## A pole bound for formal root branches

Use the rational function field \(K(\Theta)\). Consider a formal
Laurent series in a parameter \(s\), with coefficients in a
finite extension \(E\) of \(K(\Theta)\). Say its coefficients
have linear pole growth if their poles lie in a fixed finite
set of places of \(E\), and at each place the pole order of
the coefficient of \(s^j\) is at most \(C(j+1)\) for \(j\ge0\).
Finitely many terms of negative order may be priced separately.
One can equivalently bound the pole divisor of the coefficient
by \((j+1)\mathcal D\) for one fixed effective divisor
\(\mathcal D\), after increasing \(\mathcal D\) where necessary.

Such series form a subfield of \(E((s))\). For multiplication,
coefficient indices add. For inversion, factor out the leading
term and use the geometric series for a series of positive
order. At order \(j\) a product in this last expansion has at
most \(j\) factors and the sum of their positive indices is \(j\).
Inverting the leading coefficient adds only the finite set
of poles of its inverse. The class is also closed under finite
ramification of \(s\).

For a monic polynomial in \(Z\) whose coefficients have no
negative powers of \(s\) and have linear pole growth, every
formal Puiseux root has linear pole growth after ramification
over a finite extension of \(E\). To see this, take the
squarefree part over \(E((s))\) by the polynomial Euclidean
algorithm, using the preceding subfield property. Its
coefficients still have no negative powers of \(s\), since
the original polynomial is monic and integral for the
\(s\)-valuation, hence all its roots are integral.
Formal Newton–Puiseux supplies a root over
\(\overline E[[s^{1/M}]]\), with \(M\) finite. Work after
ramification. If the derivative of the squarefree polynomial
at the root has order \(b\), truncate the root through order
\(\ell-1\) for any \(\ell>b\), writing the root as
\(Z_0(s)+s^\ell Z_1(s)\). The truncation has coefficients in
a finite extension \(E'/E\). On substitution divide the
polynomial by \(s^{b+\ell}\). Its constant and linear
coefficients are integral series and the linear one has a
nonzero constant term. The first assertion for the constant
coefficient follows by Taylor expansion at the actual root.
Higher powers of \(Z_1\) have coefficients divisible by \(s\)
because \(2\ell>b+\ell\). The new equation therefore has a
simple root at \(s=0\), with constant term in \(E'\).

Subtract that constant root, divide by the constant
derivative, and write the new series \(\zeta\) as
\(\zeta=\sum g_{a,b}s^a\zeta^b\), with the terms \((0,0)\)
and \((0,1)\) absent. The pole order of \(g_{a,b}\) at any
fixed place is at most \(C(a+1)\), with a common finite
support of poles. Inductively bound the pole order of
\([s^j]\zeta\) by \(C'(2j-1)\) for \(j\ge1\), taking
\(C'\ge2C\). In a term using \(b\) coefficients of
positive indices with sum \(j-a\), the bound is
\[
                       C(a+1)+C'(2j-2a-b)
                       \ \le\ C'(2j-1).
\]
For \(a=0\) use \(b\ge2\), and for \(a\ge1\) use
\(C(a+1)\le C'(2a-1)\). This proves the claim.

Apply the claim first to the monic polynomial in \(V\).
The root has zero constant term in the centered coordinate.
Substitution of this root into the monic polynomial in \(W\)
preserves linear pole growth. Indeed a term \(u^a V^b\)
that contributes to a fixed ramified order \(j\) has
\(Ma+b\le j\). The \(b\) positive coefficient indices
from \(V^b\) have sum at most \(j\), and the original
coefficient in \(\Theta\) has degree at most \(2a\).
Apply the claim again to the resulting \(W\) polynomial,
allowing further ramification and a finite field extension.
Every actual complex branch over \(\Theta=\theta\) is
obtained in this manner. An embedding of \(\overline{K(\Theta)}\)
in \(\mathbb C\) sending \(\Theta\) to the transcendental
\(\theta\) exists, and the Puiseux factorizations then
specialize to all formal roots of the actual monic polynomials.
The analytic branch theorem supplies convergence of the
actual expansions.

The functions \(p,q\) are linear in \(W,V\) and their fixed
centers, which are rational multiples of \(\Theta\).
Hence \(L=(up)^{-1}\), \(Y=(uq)^{-1}\), and \(r=u^{-1}\)
all have linear pole growth, with finitely many negative
orders in the ramified parameter. A coordinate \(p\) or \(q\)
identically zero cannot give a finite ordinary height and is
discarded. This applies at a positive contact with an inactive
curvature share as well as with both shares positive.

## An integer auxiliary polynomial

Fix one branch and a finite coefficient field \(E\) as above.
It is a function field in \(\Theta\) over a number field.
Its field of constants \(K'\) is a finite extension of \(K\).
On the smooth projective curve with function field \(E/K'\),
the space of functions with poles bounded by an effective
divisor \(D\) has dimension at most \(1+\deg D\) over \(K'\).
This follows by mapping to the principal parts at the support
of \(D\). The dimension of their space is \(\deg D\), and
the kernel consists of constants. Consequently the coefficients
with poles bounded by \(j\mathcal D\) belong to a rational
vector space of dimension \(O(j+1)\).

Let \(B\) be a positive integer at least the largest pole
order in the ramified parameter of \(L,Y,r\).
Choose a rational polynomial \(F(L,Y,r)\) of degree at most
\(D\) and evaluate it on the branch, multiplied by \(s^{BD}\).
A monomial of degree \(b\le D\) contributes
\(s^{B(D-b)}\) times a product of \(b\) analytic series,
each with coefficient pole divisor at most \((j+1)\mathcal D\)
at index \(j\). The coefficient at index \(j\) of the
result thus belongs to the space with poles bounded
by \((j+D)\mathcal D\). Making its first \(A\) coefficients
zero imposes at most \(O(A(A+D+1))\) rational linear
conditions. Take \(A=(B+1)D\). The number \(\binom{D+3}{3}\)
of unknown rational coefficients exceeds this quadratic
bound for sufficiently large \(D\). A nonzero \(F\) exists.

After specialization at \(\Theta=\theta\), the convergent
branch then satisfies \(F(L,Y,r)=O(s^D)\).
At ordinary integer samples its values lie in a fixed rational
lattice, so they are zero for all sufficiently small \(s>0\).
If these samples accumulate, the identity theorem applied
after multiplication by \(s^{BD}\) shows that \(F\) vanishes
identically on the branch.

## Removing the count coordinate

Assume ordinary collisions accumulate on the branch.
The common-factor minors vanish identically along the branch,
since they are holomorphic after the physical normalization
and vanish at all the samples. Thus the analytic family
has a common factor of degree at least two. If its heights
are algebraically dependent over \(\mathbb C\), they satisfy
a nonzero polynomial of some finite degree with complex
coefficients. Evaluate that degree's monomials at the
integer height pairs. The resulting rational row vectors
have deficient rank over \(\mathbb C\), hence also over
\(\mathbb Q\). A nonzero rational polynomial vanishes at
all integer samples and therefore on the branch.

It remains to exclude heights algebraically independent over
\(\mathbb C\). The auxiliary identity then implies \(r\)
is algebraic over \(\mathbb C(L,Y)\). For this case work
at a regular point of the branch with small positive parameter.
Strip and normalize the constants of the two parents. For
each nonconstant coefficient let \(H_{j,k}\) be the logarithm
of its repeated-block rational multiplier, namely
\[
 H_{j,k}=c_+\log (jL+1)_{kd}
       +c_-\log\!\left((jY+1)_{(j-k)e}/(jY+1)_{je}\right).
\]
The remaining core and offset factors are nonzero rational
functions of \(L,Y\). The \(s_0=n+N\) logarithms \(H_{j,k}\)
are independent over \(\mathbb Q\) modulo constants.
Indeed algebraic independence of the heights allows one to
read the divisors in \(\mathbb C(L,Y)\). A factor at
nonintegral offset \(t/j\) in height \(L\) or \(Y\) occurs
only at its own order. A gap of at least two provides such
a factor in each block of its factorial ramp. Its indicator
on coefficient indices is the corresponding suffix, with
a negative sign and reversed block order for the decreasing
row. Those suffixes span that order's block. If both gaps
are one, the increasing row provides suffixes starting
at indices \(1,\ldots,j-1\), and the decreasing row
provides suffixes starting at \(2,\ldots,j\). Distinct heights
give distinct factor divisors, so the combined columns again
have full rank in both orders.

Write \(H\) for this list and \(G=rH\). The exponential
tuple \(\exp G\) lies in a locus of codimension at least two
over \(\mathbb C(L,Y)\). In fact, after fixing the lower
polynomial, there are finitely many monic quadratic divisors
over an algebraic closure, counted with the possible
multiplicities. Modulo each such divisor, the classes of
\(z,z^2\) are linearly independent because its constant term
is nonzero. Thus two independent linear equations are
imposed on the upper coefficients. Nonzero rational
coefficient prefactors preserve this codimension.

For completeness, let \(\lambda\) be the dimension over
\(\mathbb Q\) of the relations modulo constants among
\((H,G)\). Their projection to the \(G\) coefficients is
injective. For a basis of relations write \((A+rB)H=c\),
where \(B\) has row rank \(\lambda\). The corresponding
maximal minor of \(A+rB\) is a polynomial with nonzero
leading coefficient in the nonconstant \(r\). Hence
\(\operatorname{trdeg}_{\mathbb C(L,Y,r)}
   \mathbb C(L,Y,r,H)\le s_0-\lambda\).
Since \(\exp H\) is rational in \(L,Y\) and \(r\) is
algebraic over \(\mathbb C(L,Y)\), the upper bound on
the transcendence degree of \((H,G,\exp H,\exp G)\)
over \(\mathbb C\) is
\[
                         2+s_0-\lambda+s_0-2
                         =2s_0-\lambda .
\]
Ax's one-derivation theorem gives at least \(2s_0-\lambda+1\)
by choosing a basis of the logarithms modulo constants.
The exponential map on the local analytic germs satisfies
the required logarithmic-derivative identities, and the
derivation has constants \(\mathbb C\).
[Ax, On Schanuel's conjectures, (SD)](https://annals.math.princeton.edu/1971/93-2/p03)

This contradiction proves the height-curve conclusion.
## The stationary criterion on the boundary

By formal Newton–Puiseux over \(\overline{\mathbb Q}\), a
rational algebraic height curve has algebraic coefficients
in the expansion of either height ratio in rational powers
of the reciprocal of the active height. A transcendental
coefficient in that parameter therefore excludes accumulation.

At a boundary contact exchange the columns and height labels
if needed so \(L\) has the positive limiting share. Write
\(\sigma=\lambda_*/T_*\). Under the hypotheses of the
stationary test in the [branch theorem](MULTIVARIATE_OPPOSED_DIFFUSE_PUISEUX_HEIGHT_BRANCHES.md#a-stationary-contact-of-rational-bias),
the first change in \(T\) is an algebraic nonzero multiple of
\(T_*^{2/k}u^{1/k}\), and
\(\lambda=\lambda_0(T)+O(u)\). Use \(L/Y\) instead of \(Y/L\).
In the midpoint curvature coordinates it equals
\((T-d\lambda)/(T+e\lambda)+O(u)\), with the denominator
nonzero at the contact since the \(L\) share is positive.
The derivative of that fraction along \(\lambda=\lambda_0(T)\) is
\[
 -\frac{(d+e)(\lambda_0'(T_*)-\sigma)}
              {T_*(1+e\sigma)^2}.
\]
It is a nonzero algebraic multiple of \(1/T_*\). As
\(1/L\) is asymptotic to a positive rational multiple of
\(uT_*\), the first coefficient of \(L/Y\) in powers of
\(1/L\), at exponent \(1/k\), is a nonzero algebraic multiple
of \(T_*^{1/k-1}\). This coefficient is transcendental since
\(k\ge2\) and \(T_*\) is transcendental. The boundary branches
meeting the stationary assumptions are thus excluded.

This auxiliary-polynomial argument leaves branches on
algebraic height curves for further analysis. Its
coefficient-field estimate is stated at rational bias.
