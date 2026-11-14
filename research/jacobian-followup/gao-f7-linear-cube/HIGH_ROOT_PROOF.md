# Gao F7: the sharp root-multiplicity bound five

14 November 2026. Registered claim **GAOF7MULT1**.

For exactly Gao's F7 over C, every root of its monic degree-twelve
eliminant has multiplicity at most five. The bound is sharp: Y=0 gives
R=t^5*(t^7-2t^4+2t-1), with seven other simple roots.

Together with [GAOF7CUBE1](../gao-f7-three-point-certificate/PROOF.md),
this bound left four lists:(5,4),(5,2,2),(4,3,2),(3,2,2,2). The subsequent
[GAOF7FIVE1 proof](../gao-f7-quintuple/PROOF.md) classifies all118 normalized
fivefold targets, each with seven-point fibres. Only(4,3,2) and(3,2,2,2)
remain **OPEN**. This is a bound on roots of the eliminant on c!=0, not a new
ramification locus in the etale F7 source or a proof of JC2.

The [separate verifier](verify_high_root.py) passed
[2426 checks](high-root-receipt.json) in28.94seconds. Its
[certificate](high-root-certificate.json) is738473bytes. It rederives
root-incidence coverage, checks exact polynomial determinants using two
arithmetic engines, and replays both exceptional-case rational unit
identities. It runs no Groebner basis or multivariate resultant algorithm.
The written coverage and interpolation arguments have no independent
end-to-end review or formal verification; novelty priority is unestablished.

## All multiple roots are covered without localization

Use the eliminant from [Gao2608.00222v2](https://arxiv.org/html/2608.00222v2):

    b(t)=Y4-3Y2*t^2+2Y1*t^3-t^4,
    a(t)=Y3-2Y2*t+Y1*t^2-t^3-t^2*b(t),
    R(t;Y)=a(t)^2+a(t)+t*b(t)+t^3+Y1*t-Y2.

Suppose s is a multiple root and write a=a(s), b=b(s). Put

    G2=a^2+a+s*b+s^3,
    G3=2s*G2+a+s^3+s^2*b,
    G4=3s^2*G2+b+s^4,
    X1=-det(d(G2,G3,G4)/d(s,a,b)).

The root equation and definitions force, successively,

    Y2=G2+s*Y1, Y3=G3+s^2*Y1, Y4=G4+s^3*Y1.

These use only constant denominators1. For arbitrary xi, substituting
Y=(xi,G2+s*xi,G3+s^2*xi,G4+s^3*xi) gives

    a(s;Y)=a, b(s;Y)=b, R(s;Y)=0,
    R'(s;Y)=xi-X1.

The derivative is taken with Y fixed before substitution. Hence any
multiple root has Y=X=(X1,G2+s*X1,G3+s^2*X1,G4+s^3*X1), at arbitrary
s,a,b, including all exceptional values. The checker proves every identity
above directly from R. This also supplies the converse coverage needed
here, rather than merely checking that a parameterized subset has roots.

A root of multiplicity at least six would satisfy R^(j)(s;X)=0 for
j=2,3,4,5. Let f2,f3,f4,f5 be the retained primitive integer polynomials
in [sixfold-input.json](../gao-f7-fibres/sixfold-input.json). The checker
independently reconstructs each derivative, divided by j!, and verifies
that it is a nonzero rational multiple of f_j. Thus all four f_j would
vanish. The old30-second direct Singular attempt on these three-variable
equations was inconclusive and is retained; it is not used as evidence.

## Eliminate the common quadratic term

Every f_j is quadratic in b. Write

    f2=A*b^2+B*b+C, A=4s^2*(3s^3-2)^2,
    A_j=[b^2]f_j, j=3,4,5.

For D=s*(3s^3-2)!=0, form three linear polynomials L_i=P_i*b+Q_i:

    L1=(A*f3-A_3*f2)/(4s*(3s^3-2)),
    L2=(A*f4-A_4*f2)/4,
    L3=(A*f5-A_5*f2)/(4s^2).

The numerators are exactly divisible as polynomials. Their integer
coefficients are in [high-root-bivariate-input.json](high-root-bivariate-input.json).
At a common zero of the f_j with D!=0, all L_i vanish. No P_i is divided
out, so its zero locus and all b-leading-coefficient degeneracies remain.

The following expressions are then necessary equations in a,s:

    C_ij=P_i*Q_j-P_j*Q_i, 1<=i<j<=3,
    E_i=A*Q_i^2-B*P_i*Q_i+C*P_i^2.

Indeed C_ij=P_i*L_j-P_j*L_i, and exactly

    E_i=P_i^2*f2+L_i*(A*Q_i-A*P_i*b-B*P_i).

Divide these expressions only by nonzero rational constants and powers
of s and3s^3-2. Denote the resulting primitive equations, in stored order,
by e0,e1,e2,e3,e4,e5, where the first three are C12,C13,C23 and the last
three are E1,E2,E3. All these divisions are valid on D!=0. The checker
replays the displayed identities and every exact normalization. Only
e0,e1,e2,e3 are needed below.

## Three univariate resultants, verified without trusting their producer

Let the resultants eliminate a:

    U1(s)=Res_a(e0,e1),
    U2(s)=Res_a(e0,e3),
    U3(s)=Res_a(e0,e2).

If all e_i vanish at a finite point (a,s), then every U_k(s)=0. This
follows directly from the Sylvester matrix: the nonzero vector of powers
of a annihilates its coefficient rows. This remains true when a leading
coefficient vanishes at that s; no degree-drop locus is omitted.

The certificate records integer coefficients of U1,U2,U3, of degrees
309,341,345. The verifier proves these are the exact determinants:

| Pair | Sylvester size | Degree bound in s | Distinct integer evaluations |
|---|---:|---:|---:|
| e0,e1 |14|364|365|
| e0,e3 |15|405|406|
| e0,e2 |14|385|386|

For a-degree m,n and s-degree d_f,d_g, each Sylvester determinant has
s-degree at most n*d_f+m*d_g. The verifier checks that the stored polynomial
also satisfies this bound, then evaluates both sides at bound+1 distinct
integers, using the original nominal a-degrees even at specialization
points where they drop. At every point FLINT integer determinants and a
separate plain-Python fraction-free Bareiss determinant both equal the
stored polynomial value. Equality at more points than the degree bound
proves the full characteristic-zero polynomial identity. These are exact
integer evaluations, not modular spot checks or heuristic interpolation.

The checker then verifies literal coefficient identities

    U1=c1*s^5*(3s^3-2)^8*p1,
    U2=c2*s^14*(3s^3-2)^9*p2,
    U3=c3*s^2*p3,

where c_i are nonzero rational constants and p_i are primitive integer
polynomials of degrees280,300,343. Thus D!=0 and all U_k(s)=0 would force
p1(s)=p2(s)=p3(s)=0.

Form the full integer Sylvester matrix of p1 and p2+p3, using their actual
characteristic-zero degrees280 and343. Its size is623x623, and

    det M = 919631120 mod1000000007.

The prime is checked by complete trial division. FLINT and independent
plain-Python modular Gaussian elimination agree on this nonzero residue.
Therefore the integer determinant is nonzero, and p1 and p2+p3 have no
common complex root. In particular p1,p2,p3 have no common root. This
contradicts a sixfold root on D!=0.

This is a conclusive integer-resultant argument. In contrast, a finite-field
unit-ideal report can lose a characteristic-zero component when leading
coefficients drop. Here the full actual-degree integer Sylvester matrix is
formed before reduction. The nonzero determinant certifies the relevant
characteristic-zero resultant itself.

## Both exceptional loci are excluded by exact unit identities

For s=0, specialize the four f_j and clear only rational constants. The
[primitive equations](high-root-szero-integral-input.json) have an explicit
rational unit identity with34 multiplier terms.

For3s^3=2, take each f_j's polynomial remainder modulo3s^3-2, clear only
rational constants, and adjoin3s^3-2 itself. The
[five equations](high-root-cubic-exception-integral-input.json) have an
explicit rational unit identity with355 multiplier terms. This covers all
three complex values of s at once; none is approximated or selected away.

Both identities are in the certificate. The verifier reconstructs the
specializations and remainders from the freshly checked f_j. For each
identity it clears a positive common coefficient denominator d and checks
sum_i Lambda_i*f_i=d as an integer polynomial in both FLINT and a separate
plain-Python sparse integer convolution. Hence neither exceptional locus
contains a common zero. Together with D!=0, this exhausts every s in C
and proves there is no root of multiplicity at least six.

## Sharpness and consequences

At the already known rational target Y=(0,0,0,0), direct substitution gives

    R=t^5*H7, H7=t^7-2t^4+2t-1.

H7(0)=-1 and disc(H7)=713!=0, so zero has multiplicity exactly five and
the other seven roots are simple. The verifier checks the full identity,
the discriminant and the squarefree gcd. Thus the global maximum root
multiplicity is exactly five. The initial finite-field seed search found
the reduction of this same known rational witness immediately; no Hensel
argument or new seven-point-fibre priority claim is needed.

The earlier exhaustive three-point list was9,7+2,6+3,5+4,5+2+2,4+3+2,
3+3+3,3+2+2+2. The maximum-five bound excludes9,7+2,6+3. GAOF7CUBE1 also
excludes3+3+3 and retains its independent bound sum floor(m_r/3)<=2.
This bound left four lists:5+4,5+2+2,4+3+2,3+2+2+2. The later
[complete fivefold locus](../gao-f7-quintuple/PROOF.md) excludes both
fivefold lists, leaving only4+3+2 and3+2+2+2. By GAOF7BND1 a three-point fibre cannot
lie on c=0. Global three-point existence, full F7 strata, F5 and JC2 remain
open.

## Reproduction and retained failures

With [the pinned dependencies](requirements-high-root.txt), from the
repository root on Linux:

```sh
timeout 180s .venv/bin/python research/jacobian-followup/gao-f7-linear-cube/verify_high_root.py
```

The replay uses Python3.12.3, SymPy1.14.0 and python-flint0.9.0. Its receipt
binds checker, coefficient input and certificate by SHA256. The derivation
script regenerates byte-identical coefficient inputs; the resultant producer
and identity exporter are untrusted by the verifier.

The previous four-variable seven-two system, its0.065second modular unit
output and20second rational-lift timeout remain retained. The new bivariate
system's modular basis finished in0.016seconds, while its direct rational
lift also hit20seconds. These reports do not support the theorem. Instead,
three exact univariate resultants completed in about0.03seconds each.
An extra redundant resultant and all generation checkpoints are retained.
One exceptional-case Singular input initially used rational coefficients
in a syntax rejected by the backend despite exit code0; its input and log
are preserved. Clearing denominators and using explicit output syntax
corrected that attempt. The separate checker accepts only literal identities.
Numerical generation had20second,1536MiB address-space and10MiB output caps;
the resultant producer had a120second time and2GiB address-space cap. The
final verifier completed within its180second time cap.

Source-version checks on14November still list Gao v2 of4November,
de Bondt2610.24889v1 of17October and Groebner.jl2607.06372v3 of1August.
The research direction and original eliminant come from the retained Gao
work; no new external theorem is imported into this exclusion. No formal
proof, independent end-to-end or external review, author endorsement or
novelty priority is claimed. MATH_STATUS.json is the status authority.
