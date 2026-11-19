# The least unsafe prime degree across all two-height starts

19 November 2026. Work over a characteristic-zero
field with regular phases p=t^a(C+t²D)/x for ANY integer a>=0.
Degree means least sphere degree, equivalently radial-primitive balanced degree
of the original two-pair lift; balanced degree13 is ordinary totaldegree26.
There are exactly two displayed height terms. Purity means all positive
angular moments vanish; safety includes every fixed sphere marker and
every fixed original two-pair polynomial marker.

**Theorem TwoHeightMinimumUnsafePrimeThirteen.** Every pure phase of
prime least degree below13 in this entire class is safe. The least
prime unsafe degree is exactly13, supplied by the already retained
primitive, non-proper-power a1 hyperbolic witness. No unrestricted
prime-degree statement about arbitrary Image polynomials or additional
independent height profiles is asserted.

## 1. Even starts have even nontrivial pure degree

Normalize a unit center by C(0)=-D(0)=1. The reviewed all-start sharp
degree theorem gives rational X=P/J, gcd(P,J)=1, P0=J0=1, and

    C=HJ², D=-HP², u=xP^a/(HJ^(a+2)),
    omega=P^(a-1)/(HJ^(a+1)) dx.

For even a, all endpoint factors are nonzero, so purity gives Phi1 and
omega=du. Thus its normalized rational primitive R is exactlyu.
At a J-root of multiplicityj>0, u has pole orderh+(a+2)j whereas R
has strictly smaller orderh+(a+1)j-1. This contradicts R=u, so J1.
Write pdeg=degP,h=degH. A nontrivial unit-center phase has pdeg>=1
by the compensating-root conclusion of the sharp theorem.

The growth exponent of u at infinity is delta=1+a*pdeg-h. If delta
were nonzero, u' would have growth delta-1, whereas omega/dx has
growth(a-1)*pdeg-h=delta-pdeg-1. Equality u'=omega/dx would force
pdeg0, a contradiction. Hence delta0, so

    h=1+a*pdeg, degD=1+(a+2)*pdeg,
    degree p=(a+2)(pdeg+1).

This is even for even a. It includes a0; no odd-start rationality or
infinity hypothesis is silently applied there. The trivial unit phase
2t^a y is safe. The zero-center pure branch is safe by the same
all-start theorem. Consequently even starts contribute no unsafe odd
prime degrees. At prime degree2 the sharp lower bound excludes every
nontrivial unit phase, so that degree is safe as well.

## 2. Only one new odd-start case remains below13

For odd a>=1, every nontrivial pure unit-center phase has degree at
least2a+4. Among odd primes below13 the only candidate pairs(a,d) are

    (1,7), (1,11), (3,11).

The first two have no pure unit-center phase by the retained full
odd-cubic degree7 and degree11 theorems. Odd a>=5 already require
degree at least14. Thus only a3,d11 needs a further argument.

## 3. Starting height3, exact degree11

Put h=degH,pdeg=degP,j=degJ. The exact degree formula gives

    h+2pdeg<=7, h+2j<=9.

The sharp theorem supplies a P-root alpha with multiplicity at least1
and H-multiplicity at least3*v_alpha(P)+1>=4. Hence h>=4 and
pdeg>=1. The first bound forces pdeg1 and4<=h<=5. Normalize the
linear numerator as U=1+lambda*x, lambda!=0. Then

    P=U, H=U^4 H0, degH0<=1, degJ<=2.

The exact degree equals max(h+2j+2,h+6)=h+6, so degree11 forces
h5 and degH0=1. The auxiliary a1 phase with coefficients

    Ctilde=C/X²=U²H0 J4,
    Dtilde=D/X²=-U4H0 J²

is polynomial, regular and unit-centered. It has the SAME inverse
coordinate u and endpoint Jacobian Phi: its omega is dx/Ctilde,
which equals X²dx/C. Thus purity is equivalent for the two phases;
odd moments vanish by parity, and their even endpoint factors are
nonzero. This transformation uses P²|H, already proved by the local
multiplicity bound. It is not asserted for arbitrary coefficients.

For j0,1,2 the auxiliary degree pairs are respectively

    (degCtilde,degDtilde)=(3,5),(7,7),(11,9).

The first and last are excluded by the degree7 and degree11 unit-center
theorems. It remains to exclude a pure unit-center a1 phase with both
coefficient degrees7.

## 4. Equal coefficient degrees7 are impossible

Apply the canonical polynomial-descent theorem to that a1 phase.
The quotient degrees are equal, cbar=ebar=1+6/g, with g dividing6.
Thus cbar is7,4,3 or2. The degree2-center case has least phase degree4
and is excluded by the retained low-degree unit-center classification.

For cbar>=3 its ubar infinity zero order is cbar-1>=2. The genuine
quotient return j must fix infinity. Indeed

    j*(dz/Cbar)=dz/Cbar-2dubar.

Both terms on the right have positive order at infinity. If j mapped
infinity to a finite point, dz/Cbar would have a finite zero, impossible
for a reciprocal polynomial differential. Hence j=-z. Its Cbar roots
are paired under negation: repeated roots are detected by primitive
part poles, simple roots by nonzero logarithmic residues. Their count
Hbar is even and at least4, after the two-root exclusion.

At a paired root with multiplicities m<=n and common ubar pole orderk,
one has k<=n-1 and2k<=3m. Summing gives

    cbar>=(5/6)degubar+Hbar/2.

The infinity zero of ordercbar-1 and the simple base zero force
degubar>=cbar. Thus cbar>=3Hbar>=12, contradicting cbar<=7.
This excludes the middle auxiliary case, and hence a3,d11.

All prime degrees below13 are now covered. The retained hyperbolic
family at a1,n2 has exact balanced degree13 and an all-order fixed-marker
witness of unsafe behavior, radial primitivity and no scalar proper
power over any characteristic-zero extension. Therefore13 is the
exact least prime unsafe degree for the full two-height class.

The proof uses the reviewed ALL_START_HEIGHT_SHARP_UNSAFE_DEGREE.md,
the retained odd-cubic low-degree/degree7/degree11 theorems, canonical
polynomial descent, and HYPERBOLIC_ODD_HEIGHT_IMAGE_WITNESSES.md for
the existing upper bound. Degree19 is not a dependency. The attached
stdlib checker verifies the exact small degree-case arithmetic; it
does not replace these all-order sources or infer purity from a prefix.

The [checker](verify_two_height_minimum_unsafe_prime_thirteen.py),
[source ledger](two_height_minimum_unsafe_prime_thirteen_source_ledger.json)
and [receipt](two_height_minimum_unsafe_prime_thirteen_verification.json)
retain the exact degree reductions and their proof inputs.
