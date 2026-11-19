# Odd-cubic degree-eleven safety and the sharp prime threshold thirteen

Date: 19 November 2026. The statements concern exactly regular winding-one odd-cubic-height phases p=t(C+t²D)/x over characteristic-zero fields. No unrestricted prime-degree assertion is intended.

## 1. Statement and retained inputs

**Theorem OddCubicDegreeElevenSafety.** Every all-pure phase of exact least sphere degree11 in this class is eventually safe for every fixed sphere and original polynomial marker. In fact there is no pure unit-center phase of exact degree11.

**Corollary SharpOddCubicPrimeThresholdThirteen.** Within this entire class, the least PRIME unsafe degree is13. Retained degree<=5 and odd-degree<=7 safety, together with the theorem, give the lower bound. The ALREADY RETAINED explicit prime-degree13 witness in HYPERBOLIC_ODD_HEIGHT_IMAGE_WITNESSES.md gives the upper bound. No new witness or global minimum prime degree is claimed. The retained degree9 witness remains unsafe and is compatible with this prime statement.

The all-degree zero-center classification in ODD_START_HEIGHT_RATIONALITY_AND_QUINTIC_SAFETY.md makes that branch safe. Normalize a unit-center phase by C(0)=-D(0)=1. Put c=degC,e=degD, h=number of distinct geometric roots of C. Retained rational-u descent and exact/log separation give

    u=x*sqrt(-D/C)/C in K(x), u=x+O(x²),
    u²=-Dx²/C³,
    dx/C=dR+Omega,
    R(v(u))-R(v(-u))=2u.                              (J1)

R is the rational primitive part, normalized R(0)=0. At a C-root of multiplicity m, it has a pole of exact order m-1 if m>=2; it has no other pole, including infinity. Thus

    D_R=deg_x R=c-h.                                  (J2)

Squarefree C is impossible. One-root centers have pure degrees3c or5c-2; two-root centers are impossible. These retained exclusions remove h<=2 and squarefree cases whenever they arise below.

## 2. The genuine stable field and its ramification filter

Set T=1/(C u'), and work over the algebraic closure k. Formal purity gives T(v(u))+T(v(-u))=2. Consequently

    E=k(u,R,T) subset k(x)

has a genuine involution sigma induced by formal replacement u by-u:

    sigma(u)=-u, sigma(R)=R-2u, sigma(T)=2-T.            (J3)

This does not assert an extension of sigma to k(x). Let g=[k(x):E]. By Luroth's theorem E=k(w), so it has a rational projective line model and sigma is a Mobius involution there. The field tower gives

    g divides D_R and deg_x u,
    deg_x(R-2u)=D_R,
    deg_x u<=2D_R.                                   (J4)

Indeed sigma preserves [E:k(R)], so the degree equality follows after multiplication by g. Equality in the last bound forces R and R_sigma=R-2u to have disjoint pole supports.

The differential omega=dx/C descends to E: omega=T du. Let pi:P1_x->P1_w be the degree-g map. Suppose u has zero order s>0 at x=infinity, and let b be the ramification index of pi there. If k0 is the zero order of u as a function on P1_w at w_infinity=pi(infinity), then

    s=b*k0.

Moreover ord_infinity(omega)=c-2. For a local differential eta on P1_w, pullback multiplies its order plus1 by b, then subtracts1. This follows directly by substituting a local parameter z=t^b times a unit; characteristic zero makes the leading derivative coefficient nonzero. Thus

    b divides c-1, and b divides s.                   (J5)

At x=0 the simple zero of u makes pi unramified, with a simple zero of u at w0=pi(0); formal sigma fixes w0.

In all applications below we have

    s>=D_R,   2D_R<2s+1.                              (J6)

These inequalities force w_infinity distinct from w0. Otherwise k0=1 and b=s, so the fiber over w0 contains x0 of multiplicity1 and infinity of multiplicitys, giving g>=s+1>D_R, contrary to(J4).

They also force sigma to FIX w_infinity. If it did not, u on P1_w would have zeros of order k0 at w_infinity and sigma(w_infinity), and a simple zero at w0. Hence

    deg_x u>=g(2k0+1)>=2s+g>=2s+1,

contradicting(J4),(J6). At a fixed point of a nontrivial Mobius involution the local derivative is-1. Since sigma(u)=-u, the positive zero order k0 there must be ODD. Counting its zero together with the distinct simple zero w0 now gives the reusable filter

    b|gcd(s,c-1), b<=g, k0=s/b odd,
    deg_x u>=g(k0+1),   deg_x u<=2D_R.                (J7)

No use has been made of a rational return on x.

The logarithmic differential Omega also descends to E, because omega and dR do. It is sigma-invariant:

    sigma(omega)=(2-T)d(-u)=omega-2du,
    sigma(dR)=dR-2du, hence sigma(Omega)=Omega.          (J8)

This explicit enlargement of the field avoids assuming that Omega descends to the smaller field k(u,R).

One additional value restriction shortens the lowest case. When u vanishes at infinity, its derivative has growth order -s-1, so

    growth_order_infinity T=s+1-c=(c-e)/2.             (J8a)

If c<e, then T(infinity)=0. But the already proved sigma-fixed point w_infinity would require T(w_infinity)=sigma(T)(w_infinity)=2-T(w_infinity), hence T(w_infinity)=1. Thus c<e is impossible under the common fixed-infinity hypotheses. The other case arguments below are retained explicitly and do not require this extra filter.

## 3. Pole-capacity facts, including rational return

At a simple C-root r, polynomial D and(J1) imply any pole of u has order at most1: if that order is k, then2k<=3. If deg_x u=2D_R, disjointness of R,R_sigma confines all poles of R_sigma to simple C-roots (there is no infinity pole in the present u-zero cases). Therefore

    D_R<=number of simple C-roots.                    (J9)

When g=1, sigma is a genuine automorphism j of k(x). It fixes x0 and infinity, so j(x)=-x. Thus u is odd in x. The set of distinct C-roots is invariant under negation: a repeated root maps to a pole of R_sigma, which cannot lie outside C's roots; a simple root has a nonzero logarithmic residue and maps to another Omega-pole. Since C(0)=1, all roots occur in distinct pairs r,-r, and h is even.

For each such pair let m<=n be its C-root multiplicities, and let k>=0 be the common pole order of u there (k0 means neither point is a pole). The pole orders of R and R_sigma are m-1,n-1 in opposite order, so

    k<=n-1.

Polynomial D at the smaller root gives2k<=3m. Consequently

    m+n>=5k/3+1,
    c>= (5/6)deg_x u+h/2.                             (J10)

This is a general capacity bound; it remains valid if leading pole terms cancel.

A maximal-degree version will be useful once. If g=1 and deg_x u=2D_R, R,R_sigma are disjoint. Each repeated C-root therefore has a SIMPLE opposite root. Its pole order m-1 transfers to that simple point and must be at most1, so m=2. Each negation pair contains at most one double root, and hence

    c<=3h/2 in this maximal case.                     (J11)

## 4. Exact degree11 case split

The exact least degree is max(c,e+2). Rational sqrt(-D/C) forces c,e to have the same parity. For degree11 the only possibilities are

    c=1,3,5,7,9 with e=9,
    or c=11 with e=1,3,5,7,9.

The c1 case is squarefree. For c3, all root patterns are one-root, two-root or squarefree; the one-root degrees are9 or13, never11. Thus both are excluded. In every remaining case h>=3 and

    s=ord_zero_infinity u=(3c-e-2)/2,
    D_R=c-h<=c-3.

This gives(c,s,D_R upper bound)=(5,2,2),(7,5,4),(9,8,6), or c11 with s=15,14,13,12,11 and D_R<=8. Each satisfies(J6).

### 4.1 c5,e9

Here s2,D_R<=2 satisfy(J6), so w_infinity is sigma-fixed. But(J8a) gives T(infinity)=0 because c5<e9. The fixed point would require T=1, a contradiction. This eliminates every root-multiplicity pattern without a quotient-log subcase.

### 4.2 c7,e9

Here s5 and gcd(s,c-1)=1 give b1,k0=5. Filter(J7) forces g1, since6g<=2D_R<=8. Negation symmetry makes h even; h>=6 gives2D_R<=2<6, so h4,D_R3 and deg_x u=6=2D_R. Bound(J11) would give c<=6, contrary to c7.

### 4.3 c9,e9

Here s8 and k0=s/b odd force b=8. But b<=g<=D_R<=6. Contradiction.

### 4.4 c11

Here b|gcd(s,10), k0=s/b odd, D_R<=8. All five e-cases are explicit:

* e1,s15: b1 or5. If b5, g>=5 and deg_x u>=4g>=20>16. Thus b1, and16g<=16 forces g1.
* e3,s14: b2,k0=7. Then8g<=16, so g2,D_R8,h3 and deg_x u16=2D_R. Bound(J9) would require eight simple C-roots, although h3 permits at most two. Contradiction.
* e5,s13: b1,k0=13, and14g<=16 forces g1.
* e7,s12: b1 or2 makes k0 even, impossible.
* e9,s11: b1,k0=11, and12g<=16 forces g1.

In every remaining g1 case the root count is even and at least4, while deg_x u>=12. Capacity bound(J10) gives

    c>= (5/6)*12+4/2=12,

contradicting c11. This exhausts degree11 and proves the theorem.

## 5. Sharp prime threshold and verification scope

The zero-center branch is safe for every degree. Retained odd-cubic degree<=5 safety handles primes2,3,5; retained degree7 safety handles7; the theorem handles11. The existing primitive all-order F13 phase supplies a fixed-marker unsafe example at prime degree13. Therefore13 is the exact minimum PRIME unsafe degree within this full odd-cubic-height class. The statements do not cover the three-profile phase t(C+t²B+t4D)/x, arbitrary septic/quintic phases, or unrestricted Image polynomials.

All-marker safety imports the repaired equivariance/output-degree transfer in the retained odd-start package. No assertion that E2 kills arbitrary nonzero weights is used. The written field, local-order and pole arguments supply universal quantifiers; finite case or original-contraction controls are separate.

The retained Milne Fields and Galois Theory v5.10 Theorem9.19 (printedp117) is used exactly for E=k(w). The local ramification/differential formula is proved directly in Section2. No extension of sigma to k(x) is used except in the explicitly proved g1 branches.

The [source ledger](odd_cubic_degree_eleven_safety_source_ledger.json) pins the seven retained inputs. The [checker](verify_odd_cubic_degree_eleven_safety.py) and [receipt](odd_cubic_degree_eleven_safety_verification.json) bind exact finite controls to these proof and input bytes.

```bash
python3 research/foundations/extended-geometry/verify_odd_cubic_degree_eleven_safety.py \
  --write-receipt research/foundations/extended-geometry/odd_cubic_degree_eleven_safety_verification.json
```
