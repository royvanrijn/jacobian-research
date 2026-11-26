# No-carry prime degree gates without a height restriction

This note concerns ordinary coefficient heights, not integrality of every Hurwitz coefficient. All supporting-line and factor-degree arguments are given below.

For q>=3 set

    R_h(A)=Lambda((Ax^2+x-q)^h)/h!,  Lambda(x^r)=r!,
    H_q(X)=X^q R_q(1/X)/q!=sum_(j=0)^q a_j X^j/j!.

The exact native coefficient formula is

    a_j=sum_(k=0)^j binom(j,k)(-q)^k
                  (2q-j-k)!/[q!(q-j)!].                 (1)

Thus the individual ordinary-coefficient summand at index j is

    T_jk=(-q)^k (2q-j-k)!/
                   [q!(q-j)! k!(j-k)!],  0<=k<=j<=q.    (2)

In particular T_j0=binom(2q-j,q)/j!. The two occurrences of q! in H_q's definition are intentional: R_q already divides the moment by q!.

**Theorem 1 (unrestricted no-carry ordinary polygon).** Let p be an odd prime dividing q and t=v_p(q)>=1. Suppose every base-p digit of q is at most (p-1)/2. Then H_q has exactly the same p-adic Newton polygon as its native C=0 polynomial

    sum_(j=0)^q binom(2q-j,q) X^j/j!.

It is also the polygon of sum_(j=0)^q X^j/j!. Every rational factor of R_q has degree divisible by p^t. No bound on floor(log_p(2q)) is required.

The key strict supporting-line theorem below does not need the no-carry assumption. That assumption is used only to make all native prefix coefficients units, ensuring that the supporting vertices are actually attained.

## 1. All factorial-prefix supporting lines

Write q in base p. Starting at H=0, add its nonzero digits in decreasing order of position. If the next nonzero digit is d in position s, the successive prefixes are H and H+d p^s. Because p^t divides q, every such s satisfies s>=t>=1. Put P_s=p^s. The line through their factorial heights is

    ell_s(j)=-v_p(H!)-lambda_s(j-H),
    lambda_s=(p^s-1)/[(p-1)p^s]=sum_(e=1)^s p^(-e).     (3)

It is a supporting line for all the points (j,-v_p(j!)), 0<=j<=q. Here is a proof over the entire index range. Let s_p denote digit sum, write H=Q p^s with p|Q and j=v p^s+w, 0<=w<p^s. One has v<=Q+d. The integer function s_p(v)-v is nonincreasing, since its increment at v is -(p-1)v_p(v+1). It is constant from Q through Q+d because d<p. Also s_p(w)-w/p^s>=0. Therefore

    s_p(j)-j/p^s >= s_p(Q)-Q=s_p(H)-H/p^s.             (4)

Legendre's formula turns (4) into -v_p(j!)>=ell_s(j). The two prefix endpoints attain equality. Their widths are d p^s and slopes are -lambda_s. As successive s decrease, these slopes increase. These are consequently all edges of the factorial lower polygon, including endpoints 0 and q.

For subsequent calculations an exact alternative form of the same line is useful:

    ell_s(j)=-sum_(e<=s) j/p^e
             -sum_(e>s) floor(q/p^e).                  (5)

Indeed H is divisible by p^(s+1), while q-H<p^(s+1), so floor(H/p^e)=floor(q/p^e) for e>s.

## 2. One exact factorial level

Fix 1<=k<=j<=q. For each e>=1 set P=p^e and write

    q=Q P+u, j=J P+v, k=K P+w, 0<=u,v,w<P.

The contribution of this level to the factorial part of v_p(T_jk) is exactly

    -J-K+T,
    T=floor((2u-v-w)/P)+1_(v>u)+1_(w>v).               (6)

This follows by substituting the five factorial floors in (2); the term kt from (-q)^k has not yet been allocated.

Two elementary facts hold for every u,v,w in [0,P):

    T>=-1; if T=-1, then v+w>P.                        (7)

The floor in (6) is at least -2. If it equals -2, then v+w>2u+P, which forces v>u and supplies at least one indicator. If it equals -1, it cannot have both indicators zero: v<=u and w<=v would give 2u-v-w>=0. Thus a negative T can only come from the first case, proving both assertions in (7).

For e<=s, subtracting that level of (5) from (6) gives

    -K+T+v/P >= -k/P.                                 (8)

If T>=0 this is immediate. If T=-1, inequality (7) gives T+(v+w)/P>0, which also proves (8).

For e>s the defect instead is

    Q-J-K+T >= -floor(k/P).                            (9)

If J<Q, the integer Q-J>=1 absorbs T>=-1. If J=Q, then j<=q implies v<=u. In this case T>=0: if w<=v, the floor in (6) is nonnegative; if w>v, it is at least -1 and the second indicator cancels it. This proves (9). In particular every level e>t, regardless of whether it is above or below s, contributes defect at least -k/p^e.

## 3. The low levels pay the entire tail, strictly

For e<=t, one has u=0 and e<=s. Allocate +k of the kt valuation from (-q)^k to each of these t levels. By (6)--(7), its resulting defect is at least

    k-floor(k/P)-1+v/P >= k-floor(k/P)-1.               (10)

For k>=2 this is nonnegative at every such level, since P>=3. Keep the surplus at e=1 and discard the other nonnegative low-level defects. By (8)--(9), the entire tail e>t is bounded below by

    -sum_(e>t) k/p^e=-k/[p^t(p-1)].                    (11)

Consequently, for k>=2,

    v_p(T_jk)-ell_s(j)
      >= k-floor(k/p)-1-k/[p^t(p-1)] >0.               (12)

The last strict inequality is explicit at all boundaries. For k=2, floor(2/p)=0, so the lower bound is at least 1-2/6=2/3. For k>=3, it is at least k-k/3-1-k/6=k/2-1>0.

For k=1, at every e<=t one has w=1 and T=0 exactly. If v=0, the negative floor -1 is canceled by 1_(w>v). If v>=1, the floor is -1 and is canceled by 1_(v>0), including v=P-1. The allocated defect is therefore 1+v/P>=1. Keeping one such level and using the same tail gives

    v_p(T_j1)-ell_s(j)
       >=1-1/[p^t(p-1)]>=5/6>0.                       (13)

This proof includes k=j and j=q: all factorial arguments in (2) remain nonnegative, (6) still holds with zero remainders as appropriate, and inequalities (8)--(13) use only 1<=k<=j<=q. Index j=0 has no positive-shift term.

We have proved the following stronger local statement without any digit or height assumption: for every odd p|q, every factorial-prefix supporting line (3), and EVERY positive-shift term in (2),

    v_p(T_jk)>ell_s(j),  1<=k<=j<=q.                   (14)

The estimate does not say v_p(a_j)>=0. A supporting line can lie strictly below -v_p(j!), allowing some Hurwitz coefficients a_j to have negative valuation.

## 4. Native units attain the polygon under no carry

The native ordinary term T_j0 is binom(2q-j,q)/j!, whose valuation is at least -v_p(j!), and hence at least every line ell_s(j). By (14), adding all positive-shift terms can never put an actual coefficient below any such line.

At a high-digit prefix j=H, adding q and q-H in base p has no carry under the digit hypothesis: only the suffix digits are doubled and each digit of q is at most (p-1)/2. The factorial digit formula gives

    v_p binom(2q-H,q)=0.

Thus the native coefficient at every prefix has height exactly -v_p(H!). Every positive-shift term has strictly larger valuation there by (14). Its sum cannot cancel the native term. All factorial-prefix vertices are retained, while no other point can lie below their supporting lines. This proves equality of the entire lower polygons in Theorem 1, with nonzero endpoints and full degree. The same reasoning also identifies the native C=0 polygon with this factorial polygon.

A root of H_q has valuation lambda_s for one of these faces: otherwise its vanishing sum has a unique least-valuation term. Its reduced denominator is p^s, since

    lambda_s=(1+p+...+p^(s-1))/p^s

has numerator prime to p. All s>=t. The conjugates of a monic irreducible factor over Q_p have one common valuation; the valuation of its nonzero constant coefficient is an integer. Hence its degree is divisible by p^s, and in particular p^t. Summing local degrees proves the rational factor-degree assertion. Reciprocity preserves degrees because the endpoint coefficients are nonzero. This proof uses no Coleman irreducibility theorem.

## 5. Exact resonance consequence and boundaries

The corrected factorial recurrence at its boundary-zero index gives

    R_q=(1+4qA)R_(q-2)+2A(2q-1)R_(q-1).                (15)

Thus every rational common factor of the two preceding resonant moments divides R_q. At EVERY odd p|q whose base-p digits have no doubling carry, its degree is divisible by p^v_p(q), with no height test. These distinct prime-power divisors may be multiplied. They may also be multiplied by the separately proved power-of-two divisor D(q), because the odd part is coprime to D(q). The resulting degree must still be <=q-2. This is an enlarged sufficient exclusion criterion; it does not settle all odd resonances.

**Corollary 2 (stronger partial dense-resonance exclusion).** Let K have characteristic zero and f=ax^2+bx+c with abc!=0. For odd q with at least four set bits in q-1, let D(q) be 2 to the power of the fourth-lowest set-bit position. Put

    J_0(q)=product_(p|q, p odd, all base-p digits of q<=(p-1)/2)
                      p^v_p(q).

If three consecutive factorial moments of f vanish beginning at n>=1, then

    q=-c/b=n+2 is odd,  s_2(q-1)>=4,
    D(q) J_0(q) <= q-2.                                (16)

No height bound is part of the definition of J_0. Hence the necessary set in (16) is contained in the previously proved height-guarded residual set. This is a necessary condition for a counterexample, not closure of every remaining odd q.

Here the exact dense resonance reduction, even-q exclusion, first-three binary-cluster exclusion and binary common-factor degree D(q) are the separately proved local statements in [the dense quadratic resonance note](UNIVARIATE_QUADRATIC_STRONG_FACTORIAL_RESONANCE.md), Sections1--5. Given a putative triple, these imply odd q, at least four set bits, and a nonconstant rational common factor G of R_(q-2),R_(q-1) with D(q)|deg G. Nonzero rational moment polynomials make any common value-zero algebraic over Q; taking its minimal polynomial produces G even when the coefficient field is not a subfield of C. Equation(15) puts G in R_q. Theorem1 at all no-carry odd prime divisors gives J_0(q)|deg G. Distinct prime powers and the power of two D(q) are coprime. Their product divides deg G<=q-2, proving(16). The imported binary/even statements are used only here; they are not premises of Theorem1.

The former counterfamily p=5,q=5^h+10,j=5^h-1,h>=19 remains a valid counterexample to no-carry-only Hurwitz integrality: its a_j has valuation 18-h<0. It does not contradict (14), whose ordinary supporting-line baseline is lower. No statement of integrality is inferred here. In particular this proof corrects neither that counterexample nor the previous height-guarded theorem; it supplies a different, stronger ordinary-polygon route.

No full dense quadratic Strong Factorial theorem, full SFC(2), SFC(3), worldwide priority or formal verification is asserted. The all-range supporting-line proof is (3)--(14).

## 6. Supplementary exact controls

The [portable standard-library checker](verify_univariate_quadratic_no_height_newton_polygon.py) and [retained finite receipt](univariate_quadratic_no_height_newton_polygon_verification.json) check the floor identities and strict term bounds on explicitly finite ranges, plus selected huge-family boundary terms using integer floor divisions only. The checker does not expand enormous factorials, read external inputs, or infer a universal theorem from a prefix. From the directory containing this note, run a fresh replay with ordinary Python:

```text
python verify_univariate_quadratic_no_height_newton_polygon.py --output fresh_no_height_verification.json
```

An existing output file is never overwritten. These controls are supplementary evidence; Theorem1 is the written all-range argument above. The fixed prime/carry polygon certificates and the earlier long prime-power proof are not premises.
