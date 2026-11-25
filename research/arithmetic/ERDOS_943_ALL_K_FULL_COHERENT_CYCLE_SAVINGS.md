# Explicit pointwise bounds for sums of two \(k\)-full integers

For every fixed fullness order \(k\ge3\), this proof gives the uniform prescribed-sum exponent \(33/(33k+35)\). Its difference from \(1/(k+1)\) is \(2/((k+1)(33k+35))\). The primitive positive-coefficient estimate is uniform, and an exact common-divisor normalization gives the result for all positive targets, including noncoprime summands.

The proof uses actual selected-edge graphs, original and complementary integer vertex potentials, and unit \(k\)-th-power congruences to recover missing roots. Eleven parameterized nonnegative identities, with their mirrors, bound every real exponent box. The older aggregate prescribed-sum theorem and any qualitative compactness saving are not premises. The three precise classical inputs are stated in Section 11.

## 1. Statement and canonical factors


A positive integer is k-full if each positive prime exponent is at least k. Let r_k(N) count ordered positive k-full value pairs with sum N.

**Theorem 1.** For every fixed integer \(k\ge3\), every \(\epsilon>0\), and every positive integer \(N\),
\[
                  r_k(N)\ll_{k,\epsilon}
                      N^{33/(33k+35)+\epsilon}.             \tag{1}
\]
The same exponent holds uniformly for every fixed positive coefficient pair in the primitive canonical equation (3). The saving below \(1/(k+1)\) is explicitly
\[
 \delta_k=\frac1{k+1}-\frac{33}{33k+35}
          =\frac{2}{(k+1)(33k+35)}>0.
\]
Uniformity of the implied constant in \(k\) is not asserted. For \(k=3\), the exponent is \(33/134\); that fixed numerical exponent is not asserted for all \(k\).

Every k-full integer has a unique factorization
\[
               m=a_k^k\prod_{j=k+1}^{2k-1}a_j^j,             \tag{2}
\]
where the coordinates a_j for j>k are squarefree and pairwise coprime. The unrestricted a_k may share primes with those coordinates. At a prime exponent e>=k write e=qk+r,0<=r<k. If r=0, put p^q in a_k. If r>0, put p once in a_(k+r) and p^(q-1) in a_k. This proves existence and uniqueness, including m=1. Replacing ANY one coordinate by a divisor of itself preserves all canonical conditions.

Consider positive canonical tuples satisfying
\[
 \kappa_1\prod_{j=k}^{2k-1}a_j^j+
 \kappa_2\prod_{j=k}^{2k-1}c_j^j=n,\qquad
 \gcd\left(\kappa_1\prod a_j^j,\kappa_2\prod c_j^j\right)=1. \tag{3}
\]
The positive integer coefficient pair is fixed during a count, and all implied constants are uniform in it. A nonempty problem has kappa_i<n. Every factor of either full summand is a unit modulo n. Coordinates on opposite sides of ONE representation are coprime. No corresponding cross-coprimality between unrelated representations is assumed.

Split all coordinates into ordinary dyadic intervals [A_j,2A_j),[C_j,2C_j), with lower endpoints at least1. Their lower exponents are nonnegative and satisfy
\[
       \sum_{j=k}^{2k-1}j\log_n A_j\le1,\qquad
       \sum_{j=k}^{2k-1}j\log_n C_j\le1.                     \tag{4}
\]
The fixed factors2 and the number O_k(log(2n)^(2k)) of boxes cost only n^epsilon. Hence it suffices to prove the displayed explicit exponent uniformly over all these boxes. Empty boxes and n=1 cause no difficulty.

## 2. Fixed pairs, omitted-coordinate codegrees and all heights

For every fixed degree e>=3, primitive positive solutions of P u^e+Q v^e=n, with P,Q nonzero integers, number O_(e,epsilon)(n^epsilon) uniformly in the coefficients. Indeed F(X,Y)=P X^e+QY^e has e distinct projective linear factors over the complex numbers. Evertse's stated theorem bounds its primitive Thue--Mahler solutions with primes dividing n by a constant depending on e raised to1+omega(n). The fixed value equation is a subset, and such a factor is n^epsilon for fixed e. Full summand coprimality in (3) supplies gcd(u,v)=1 whenever coordinates are omitted. We use only e in the fixed range[k,2k-1].

For an omitted exponent e, form a simple bipartite graph on the tuples of the other canonical coordinates. Keep one selected actual representation on each nonempty vertex pair. Every edge has at most n^epsilon labels by the fixed-pair result. Let P_1,P_2 include kappa_1 and the products at two distinct left vertices, and let Q include kappa_2 and the product at a shared right vertex. Selected edges satisfy
\[
                      P_i u_i^e+Q v_i^e=n\quad(i=1,2).
\]
Write u_i=g*alpha_i with gcd(alpha_1,alpha_2)=1, and v_i=h*gamma_i with gcd(gamma_1,gamma_2)=1. Since each actual pair(u_i,v_i) is coprime, the fraction
\[
               \frac uv=\frac{\alpha_1\gamma_2}{\alpha_2\gamma_1}
\]
is reduced: possible common factors would lie either in the two alpha coordinates, in the two gamma coordinates, or in a SAME-EDGE alpha/gamma pair. All four factors are n-units. Elimination gives
\[
 (u/v)^e\equiv P_2/P_1\pmod n,\qquad
        g^e(P_1u^e-P_2v^e)=n(\gamma_2^e-\gamma_1^e).        \tag{5}
\]
There are at most(2e)^omega(n) unit root classes in the congruence. For odd prime powers the unit group is cyclic and the e-th power kernel has at most e elements; for powers of2 it has at most2e. This includes primes dividing e. Each nonempty fiber is a coset of that kernel.

If both coordinates of primitive positive fractions in a single residue class are at most H, two distinct fractions have determinant a nonzero multiple of n. Their vectors have norm at most sqrt(2)H, so their angles are separated by at least n/(2H^2). The positive quadrant has angle length pi/2. Therefore the number of fractions in that class is at most1+pi*H^2/n. If the proposed separation exceeds the quadrant, there simply cannot be two fractions. This proves the estimate at EVERY height, including the height boundary.

For a fixed reduced fraction u/v, there are divisor-many ways to split u=alpha_1*gamma_2 and v=alpha_2*gamma_1. In an actual collision the coefficient P_1u^e-P_2v^e in (5) cannot vanish: otherwise gamma_1=gamma_2=1, so v_1=v_2, and the two complete first summands are equal. After removing known kappa_1, uniqueness of (2) makes the left vertices equal, a contradiction. Thus (5) determines the positive g uniquely, or rejects the tuple. Then
\[
             Qh^e=\frac{n-P_1(g\alpha_1)^e}{\gamma_1^e}
\]
is determined. Divide by known kappa_2 and use (2) to recover the entire right vertex and h uniquely. This decomposition is valid because replacing the omitted coordinate v_1 by its divisor h preserves canonical conditions, even when e>k is a squarefree coordinate. Candidates failing an integrality or positivity condition are rejected. All divisor arguments have polynomial height in n for fixed k.

Consequently, if the two omitted roots have dyadic scales n^u,n^v, the maximum codegree on either side is
\[
                \Delta\ll_{k,\epsilon}n^\epsilon
                                    (1+n^{2(u+v)-1}).       \tag{6}
\]
For part caps n^p,n^q, double counting neighbor pairs and Cauchy--Schwarz give E<<n^p+n^q+sqrt(Delta)*min(n^(p+q/2),n^(q+p/2)). Multiplying by the subpower selected-edge label multiplicity preserves the same exponent with epsilon. Thus the all-height omitted-coordinate exponent is
\[
 \max\{p,q,\min(p+q/2,q+p/2)+(u+v-1/2)^+\}.                \tag{7}
\]
This is an actual codegree calculation; it does not identify free interval radicands with canonical kernel monomials.

## 3. The five old box bounds

For one side put
\[
 x=\log_n A_k,\quad y=\log_n A_{k+1},\quad
 z=\sum_{j=k+2}^{2k-1}\log_n A_j,
\]
and use x',y',z' for the other side. Here k>=3 ensures the sum defining z contains at least the exponent k+2. Set
\[
 p=y+z,\ p'=y'+z',\quad
 r=x+z,\ r'=x'+z',\quad S=x+x'.                            \tag{8}
\]
The actual resources imply the weaker aggregate constraints
\[
                kx+(k+1)y+(k+2)z\le1                       \tag{9}
\]
and its primed twin. Coordinates with exponents larger than k+2 only strengthen these aggregate constraints.

Fixing a complete side of a representation fixes the complementary full summand; division by known kappa and (2) recovers its tuple uniquely. This gives L=x+y+z=x+p and L'=x'+p'. Formula(7), omitting k and k+1 respectively, gives
\[
 G_k=\max\{p,p',\min(p+p'/2,p'+p/2)+(S-1/2)^+\},
\]
\[
 G_{k+1}=\max\{r,r',\min(r+r'/2,r'+r/2)\}.                 \tag{10}
\]
The second height loss is zero because y+y'<=2/(k+1)<=1/2. Dyadic boundary constants in (6) do not change this exponent.

There is also the coefficient-uniform mixed-curve bound
\[
 J=z+z'+\min\left\{x+y'+\frac{\max(y,x')}{k+1},
                         y+x'+\frac{\max(x,y')}{k+1}\right\}. \tag{11}
\]
To prove it, fix all coordinates of exponents at least k+2, at cost n^(z+z'). Fix the coordinate of exponent k on the left and of exponent k+1 on the right. The remaining two unknowns satisfy an affine curve P U^(k+1)+Q V^k=n with P,Q positive. This curve is absolutely irreducible: as a polynomial in V over C[U], it is Eisenstein at any simple root of n-P U^(k+1). Such a root exists and is simple since n,P are nonzero in characteristic0. The degree is k+1. Bombieri--Pila therefore bounds its integer points in the enclosing square by O_(k,epsilon)(max(n^y,n^x')^(1/(k+1)+epsilon)), uniformly in all coefficients. For heights below the theorem's fixed degree-dependent threshold, the finite number of lattice points is absorbed into a constant depending on k. This proves the first term; reversing the fixed cross pair gives the second.

Hence the exponent
\[
                           F_{\rm old}=\min(L,L',G_k,G_{k+1},J) \tag{12}
\]
is valid in every box, uniformly in the fixed positive coefficients. This has been rederived here rather than imported as an old aggregate representation theorem.

## 4. General original and complementary k-th-power walks

Use the graph omitting k. Its vertex values are canonical kernels K=prod_(j>k)b_j^j, with squarefree pairwise coprime coordinates. Edges have selected primitive equations kappa_1*K_i*a_ij^k+kappa_2*L_j*c_ij^k=n. All roots and kernel factors are units modulo n.

An ordered simple six-cycle has positions11,22,33,12,23,31. Define the two alternating six-root products U=a11*a22*a33*c12*c23*c31 and V=a12*a23*a31*c11*c22*c33. A coherent cycle has U=V. Exactly as for rational edge ratios, it has unique globally primitive positive labels P_i,Q_j, with
\[
 a_{ij}=P_i/g_{ij},\quad c_{ij}=Q_j/g_{ij},\quad
                         g_{ij}=\gcd(P_i,Q_j).              \tag{13}
\]
Five tree edges give labels at most the product of five upper root caps, and in particular at most n^5 after harmless clipping of caps. Set X_i=kappa_1*K_i*P_i^k, Y_j=kappa_2*L_j*Q_j^k. Then X_i+Y_j=n*g_ij^k.

Fix an actual anchor edge and its prescribed original intrinsic gcd. This fixes its two labels and potentials. A next edge gcd divides a known label. The unknown potential is n*g^k minus the known one. Reject nonpositive or nondivisible values; divide by the known coefficient and use the UNIQUE full-k decomposition(2) to recover the next label and kernel coordinates. This is legitimate with the label as the unrestricted coordinate of exponent k. Four alternating steps visit all four other vertices. At every step prune label heights before making another divisor choice; all dividends are at most n^5 and positive potentials at most n^(1+5k). Final gcd, closing-edge, primitivity, canonical, cap, selected-edge and distinct-vertex checks only discard candidates. Thus an anchor plus original gcd gives n^epsilon cycles.

Let M=lcm(P_i,Q_j). At any prime, global primitivity makes the minimum of the six label valuations0; their total variation around the cycle is at least twice their maximum. Positivity and U=V give
\[
 v_p(U)=\tfrac12\sum_{vw}|v_p(z_v)-v_p(z_w)|\ge v_p(M),
\]
so M divides U. The complementary labels barP_i=M/P_i,barQ_j=M/Q_j are globally primitive, and with bar g_ij=gcd(barP_i,barQ_j),
\[
 a_{ij}=\bar Q_j/\bar g_{ij},\quad
 c_{ij}=\bar P_i/\bar g_{ij},\quad
                  g_{ij}\bar g_{ij}=M/(a_{ij}c_{ij}).       \tag{14}
\]
All complementary labels are also bounded by five exchanged root factors, hence n^5.

An anchor plus its prescribed complementary gcd fixes barP_1,barQ_1. At the next edge, a divisor of known barQ fixes a. FULL summand coprimality ensures gcd(a,kappa_2*L)=gcd(a,n)=1. The missing root obeys
\[
                     c^k\equiv n(\kappa_2L)^{-1}\pmod{a^k}. \tag{15}
\]
There are at most(2k)^omega(a) unit residue classes, and at most1+R_2/a^k representatives per class up to R_2. Recover K=(n-kappa_2*L*c^k)/(kappa_1*a^k), retaining only positive integral kernels whose prime exponents are among k+1,...,2k-1; their squarefree canonical coordinates are unique. Set the new complementary label to bar g*c. There is no free multiplier left. The alternating fixed-c step uses a^k congruent to n/(kappa_1*K) modulo c^k and is identical. Four steps contain two of each type. The same mandatory height pruning and all actual closing/selection/gcd checks apply, giving at most
\[
 n^\epsilon\mathcal L,\qquad
 \mathcal L=(1+R_2/r_1^k)^2(1+R_1/r_2^k)^2                 \tag{16}
\]
cycles per anchor plus prescribed complementary gcd, where r_i,R_i are the two lower and upper root caps. This includes a fixed root equal to1: its single residue class contains the full interval, already charged by(16).

Since M<=U<=(R_1R_2)^3, (14) gives every g*bar g<=D=(R_1R_2)^3/(r_1r_2). Original minimum gcd at most Y contributes at most E*Y*n^epsilon cycles. If that minimum exceeds Y, all complementary gcds are at most D/Y, so the complementary walk contributes at most E*(D/Y)*mathcalL*n^epsilon. Taking Y=sqrt(D*mathcalL) gives
\[
 C_{6,\rm coh}\ll_{k,\epsilon}En^\epsilon\sqrt{D\mathcal L}.
\]
In a dyadic box r_1~R_1~n^x,r_2~R_2~n^x', this is
\[
 C_{6,\rm coh}\ll_{k,\epsilon}E n^{w+\epsilon},\quad
       w=S+(x'-kx)^++(x-kx')^+.                            \tag{17}
\]
Thresholds less than1 cause no issue since the corresponding integer-gcd sector is empty. This is a counted anchor-and-gcd argument, not a fixed-label completion bound.

## 5. Noncoherent products and the new cycle exponent

The actual edge equations imply U^k=V^k modulo n. If U!=V, write U=g*u,V=g*v with gcd(u,v)=1. All factors are n-units, so g can be canceled. Since the nonzero integer u^k-v^k is divisible by n, h=max(u,v)>=n^(1/k). At most(2k)^omega(n) unit root classes are possible for u/v. For h in one dyadic interval of scale H, the angular bound from Section2 counts O_k(n^epsilon(1+H^2/n)) primitive fractions. Restore the common g in at most R/H ways, where R=(R_1R_2)^3. Summing the geometric dyadic series for H>=n^(1/k) gives
\[
       \#\{U,V\le R:(UV,n)=1,\ U\ne V,\ U^k=V^k\bmod n\}
              \ll_{k,\epsilon}n^\epsilon(Rn^{-1/k}+R^2/n). \tag{18}
\]
Empty height intervals contribute zero. Each product pair has n^epsilon allocations into its two ordered sextuples of root factors, by the fixed-order divisor bound. Given these twelve roots, the six equations are linear in the six kernel values. Their determinant is, up to sign,
\[
                     \kappa_1^3\kappa_2^3(U^k-V^k),
\]
which is nonzero for positive distinct U,V. Hence the six kernel values are unique; if they have the required canonical kernel form, their coordinate tuples are unique too. Thus all actual noncoherent cycles in the box satisfy
\[
 C_{6,\rm noncoh}\ll_{k,\epsilon}n^\epsilon
                     \{n^{3S-1/k}+n^{6S-1}\}.              \tag{19}
\]

For a simple graph, the quoted Morris--Saxton theorem gives E<<v^(4/3)+v*C6^(1/6). Here v<<n^V with V=max(p,p'). Insert(17),(19), take sixth roots termwise and absorb the E^(1/6) self-term. Explicitly if E<=C(B+D*E^(1/6)), either E<=2CB or E<=(2CD)^(6/5). After epsilon reallocation, this gives the valid box exponent
\[
 C_k=\max\left\{\frac{4V}{3},\frac{6V+w}{5},
                    V+\frac S2-\frac1{6k},\ V+S-\frac16\right\}. \tag{20}
\]
Root-label multiplicities per selected edge remain subpower by Section2. This step does not insert an anticipated representation estimate for E into a cycle bound.

## 6. Transfer of the nine ordinary integer identities

For the numerical certificate only, rename the aggregate coordinates by
\[
 (a,b,c;d,e,f)=(x,y,z;x',y',z').
\]
Thus \(p=b+c\), \(q=e+f\), \(r=a+c\), \(s=d+f\),
\(A=a+d\), and \(V=\max(p,q)\); these are the earlier \(p,p',r,r',S,V\).
Let
\[
 R_k=1-ka-(k+1)b-(k+2)c\ge0,\qquad
 S_k=1-kd-(k+1)e-(k+2)f\ge0.
\]
Write \(J_k\) for the expression \(J\) in (11). When \(A\le1/2\), the two graph expressions in (10) have no positive-part height loss. In this notation, write \(w_k=w=A+(d-ka)^++(a-kd)^+\) for (17). On the cone \(a\le3d\), \(d\le3a\), it equals \(A\) for every \(k\ge3\). Outside this cone, only the five bounds in (12) will be used.


For a fixed nonnegative t dominated by the relevant bounds, write
\[
 H=L-t,\quad H'=L'-t,\quad
 M_1=2p+q-2t,\quad M_2=p+2q-2t,
\]
\[
 N_1=2r+s-2t,\quad N_2=r+2s-2t,\quad
 I_p=p-t,\quad I_r=r-t,\quad I_s=s-t.
\]
Define J_b,J_d,J_a,J_e with the denominator4 EXACTLY as follows:
\[
 J_b=c+f+a+e+b/4-t,\quad J_d=c+f+a+e+d/4-t,
\]
\[
 J_a=c+f+b+d+a/4-t,\quad J_e=c+f+b+d+e/4-t.                 \tag{21}
\]
Since k+1>=4 and all coordinates are nonnegative, J_k is at most the analogous J_3 expression. More specifically, J_k>=t makes both corresponding terms in (11) at least t; replacing1/(k+1) by1/4 can only increase each. Thus the selected J_b or J_d is nonnegative according as b>=d or d>=b, and the selected J_a or J_e is nonnegative according as a>=e or e>=a. No unselected slack is assumed nonnegative.

For the first two cycle terms with V=q, put T1=4q/3-t and T2=A+6q-5t. The following table specifies nine EXACT identities. Each row with entries K,D,alpha,beta,gamma,gamma',Q means
\[
 \begin{split}
 K-[D+K(k-3)]t={}&\alpha R_k+\beta S_k
       +[\gamma+\alpha(k-3)]H
       +[\gamma'+\beta(k-3)]H'+Q.                         \tag{22}
 \end{split}
\]
In every row alpha+beta=K. The meanings of all other slacks are as just defined.

| Case | K | D | alpha | beta | gamma | gamma' | Q |
|---|---:|---:|---:|---:|---:|---:|---|
| p>=t,r>=t |25|104|25|0|71|0|21c+29Ip+4Ir|
| p>=t,s>=t,b>=d |100|416|58|42|90|84|5b+116Ip+42Is+84Jb|
| p>=t,s>=t,d>=b |200|832|137|63|243|84|31b+274Ip+63Is+168Jd|
| p>=t,G(k+1) mixed |75|312|58|17|116|22|46e+34f+116Ip+29N1|
| Gk mixed,r>=t,a>=e |99|405|67|32|126|64|52M1+6M2+67Ir+32Ja|
| Gk mixed,r>=t,e>=a |11|45|8|3|16|5|6M1+8Ir+4Je|
| both mixed,a>=3d |61|256|15|46|0|158|26M1+8b+23N2+22(a-3d)|
| both mixed,V=q,T1>=0 |12|49|4|8|8|16|4M1+4N2+9T1|
| both mixed,V=q,T2>=0 |33|134|11|22|16|38|14M1+11N2+6T2|

Every identity follows directly from its k=3 integer identity by the exact substitutions
\[
              R_3=R_k+(k-3)(H+t),\qquad
              S_3=S_k+(k-3)(H'+t).                        \tag{23}
\]
Moving the extra K(k-3)t to the left proves (22). This is also an immediate direct polynomial verification. All coefficients multiplying H,H' in (22) are nonnegative for k>=3; all entries in the table are nonnegative, and each Q has the required signs in its stated case.

The identity for the last table row, which supplies the final exponent, reads explicitly
\[
 \begin{split}
 33-(33k+35)t={}&11R_k+22S_k+(11k-17)H+(22k-28)H'\\
               &+14M_1+11N_2+6(A+6q-5t).                  \tag{24}
 \end{split}
\]
The entire table has a mirror under(a,b,c)<->(d,e,f); this exchanges the corresponding H/M/N/J slacks and preserves the domain and bounds.

For any ordinary row, its conclusion is
\[
 t\le\frac{K}{D+K(k-3)}\le\frac{33}{33k+35}.               \tag{25}
\]
All denominators are positive for k>=3. The second inequality is equivalent, after cross multiplication and cancellation of the k terms, to134K<=33D. For the nine rows the differences33D-134K are respectively82,328,656,246,99,11,274,9,0. Thus this comparison is exact and requires no numerical optimization.

## 7. The two exceptional cycle identities

The remaining cycle slacks, for V=q, are
\[
              T_{3,k}=q+A/2-1/(6k)-t,\qquad
              T_4=q+A-1/6-t.
\]
Their resource coefficients do not sum to the old left-hand constant, because the slacks themselves have constant terms. The correct exact identities are
\[
 \begin{split}
 (27-6/k)-(27k+27)t={}&9R_k+18S_k+9(k-3)H+(18k-36)H'\\
                         &+18M_1+9N_2+36T_{3,k},           \tag{26}
 \end{split}
\]
\[
 \begin{split}
 5-(6k+4)t={}&3R_k+3S_k+3(k-3)H+3(k-3)H'\\
                         &+6M_1+N_1+N_2+6T_4.             \tag{27}
 \end{split}
\]
All coefficients are nonnegative for k>=3. For (26), substitute(23) in the old25-108t identity and use T_(3,3)=T_(3,k)-(k-3)/(18k); moving that constant gives exactly27-6/k. Formula(27) follows from the old5-22t identity, whose resource weights sum to6 rather than5.

In their respective cases the identities yield
\[
 t\le\frac{9k-2}{9k(k+1)}\le\frac{33}{33k+35},\qquad
 t\le\frac5{6k+4}\le\frac{33}{33k+35}.                    \tag{28}
\]
The numerator of the first comparison after positive cross multiplication is48k+70; the second is33k-43. Both are positive for k>=3. The mirror identities handle V=p.

## 8. Complete coverage of the real domain

First handle A>1/2. If k>=4 this is impossible since ka,kd<=1. Thus k=3, where the resources imply p+q<(2-3/2)/4=1/8. Fixing kernel vertices and using the coefficient-uniform primitive binary cubic gives that representation exponent. It is strictly below33/134. No high-height codegree assertion is needed for this dispatch.

Now A<=1/2. Let t be any nonnegative number below all the five old box exponents, and also below C_k when inside the old cone. Then Rk,Sk,H,H' are nonnegative.

If p>=t, use the first four table rows: G(k+1)>=t forces r>=t, or s>=t, or both N1,N2>=0. In the s case choose the actual max(b,d) in J. The mirrors treat q>=t. If neither individual applies, p,q<t and Gk>=t forces M1,M2>=0.

Next if r>=t use table rows5--6, selecting the actual max(a,e); mirror for s>=t. Otherwise G(k+1)>=t forces N1,N2>=0. Thus both graphs are mixed.

Outside the old cone, row7 or its mirror applies, including cone boundaries. Inside the cone w_k=A and C_k>=t. If V=q, at least one of T1,T2,T_(3,k),T4 is nonnegative. Use its corresponding identity, namely row8,row9,(26),or(27). All other slacks in that identity have already been shown nonnegative. Mirror if V=p.

This exhausts every possible branch and bounds t by33/(33k+35). Therefore that exponent bounds the actual box count. The proof is on the entire real polytope, not on a finite list of exponents, and the constants in its arithmetic inputs remain uniform for every fixed coefficient pair.

## 9. Attainment for the specified menu

Put D=33k+35 and
\[
              (a,b,c;d,e,f)=(12,7,14;9,13,11)/D.             \tag{29}
\]
Use only the canonical coordinates of exponents k,k+1,k+2; set all higher exponent coordinates to zero. Both full resources equal1, since
\[
 12k+7(k+1)+14(k+2)=9k+13(k+1)+11(k+2)=D.
\]
The point lies strictly inside the old cone and has A=21/D<1/2. The two sides and two graph terms all equal33/D, while
\[
                         J_k=(41+13/(k+1))/D>33/D.
\]
The first two cycle terms are32/D and33/D. The third is69/(2D)-1/(6k), smaller than33/D by(D-9k)/(6kD)>0. The fourth is45/D-1/6, smaller by(D-72)/(6D)>0. Hence the particular five-old-bound plus cycle menu has exact maximum33/D, attained at(29). This says nothing about additional omitted-coordinate graphs, rectangular supersaturation improvements or other methods, which may lower the maximum.

## 10. Restore noncoprime original value pairs

Sections 6--8 show that every primitive dyadic box in (3) has representation exponent at most \(33/(33k+35)\). Summing the logarithmic number of boxes and absorbing their subpower label multiplicities therefore gives
\[
 \#\{\text{primitive canonical solutions of (3)}\}
       \ll_{k,\epsilon}n^{33/(33k+35)+\epsilon},             \tag{30}
\]
uniformly in all fixed positive coefficient pairs. This uses the direct arithmetic arguments of Sections 1--5, rather than an older aggregate representation theorem.

Let \(m_1+m_2=N\) be positive \(k\)-full integers and \(h=\gcd(m_1,m_2)\). Then \(h\) is \(k\)-full, and the quotients \(m_1/h,m_2/h\) are coprime. From either quotient, extract exactly those prime powers whose positive exponents are below \(k\) into a coefficient \(\kappa_i\). Every such prime divides \(h\), and its exponent in \(h\) is at least \(k\), so the extracted prime power divides \(h\). Hence \(\kappa_i\mid h\). The remaining quotient is \(k\)-full and has the unique canonical form (2). We obtain the primitive equation (3) at \(n=N/h\).

There are at most \(\tau(N)^3\) choices of \(h,\kappa_1,\kappa_2\). For each choice, a normalized tuple recovers the original values uniquely; invalid reconstructions are discarded. Apply (30), use \(N/h\le N\) and the positivity of \(33/(33k+35)\), and absorb the divisor factors in \(N^\epsilon\). This proves Theorem 1 for every positive target without any coprimality restriction on the original pair. Extracted coefficients need not be squarefree; their exponents below \(k\) are retained exactly.

The coefficient theorem used here is primitive-full-summand uniformity. No additional theorem for arbitrary nonprimitive weighted sums is being asserted.

## 11. Inputs and comparison scope


The proof uses three precise classical counting inputs.

* [Evertse, 1985, Theorem 6](sources/evertse-uniform-thue-mahler-1985.pdf), printed pages 84--85: for an integer binary form of fixed degree \(e\ge3\) with at least three distinct complex projective roots, the number of primitive Thue--Mahler solutions supported on \(t\) specified primes is at most \(7^{e^3(2t+3)}\), independently of the coefficients. At a fixed nonzero value \(n\), this gives \(O_{e,\epsilon}(|n|^\epsilon)\). Reducibility is not an excluded case. The theorem is used only in the fixed degree range \(k,\ldots,2k-1\). Its supplied derivation imports the underlying \(S\)-unit theorem; no new whole-proof audit of that deep input is asserted.
* [Bombieri--Pila, 1989, Theorem 5](sources/bombieri-pila-arcs-ovals-1989.pdf), section 3, printed page 17: integer points on an absolutely irreducible affine plane curve of fixed degree \(d\ge2\), in a square of side \(H\), satisfy a coefficient-independent bound \(O_{d,\epsilon}(H^{1/d+\epsilon})\). Its stated large-\(H\) threshold and subpower factor depend only on \(d\). The bounded smaller range is absorbed into the fixed-degree constant. Section 3 proves absolute irreducibility of the particular degree-\(k+1\) curves before applying this theorem.
* [Morris--Saxton, Theorem 1.5(a)](sources/morris-saxton-1309.2927v3.pdf), printed page 4, with cycle parameter three: above \(C v^{4/3}\) edges, a simple graph on \(v\) vertices contains at least \(c(E/v)^6\) simple six-cycles. This supplies the graph inequality used in Section 5. No conjectural rectangular supersaturation theorem or additional balanced-codegree assertion is used.

The canonical full-\(k\) factorization and omitted-coordinate graph construction are antecedents in [Powerful sum graphs](POWERFUL_SUM_GRAPHS.md), which records the exponent \(1/(k+1)\) for \(k\ge3\). The mixed-curve mechanism also appears in the [mixed-power product argument](ERDOS_943_DUAL_GRAPHS_MIXED_POWER_PRODUCTS.md). Every codegree, reconstruction, mixed-curve application and normalization required here is rederived above; an older aggregate sum bound is not a premise of the saving proof.

For the same ordered positive \(k\)-full prescribed-sum quantity, [Blomer--Schöbel, Corollary 2](https://doi.org/10.7169/facm/2013.49.2.12) gives the uniform benchmark \(N^{2/(2k+1)+\epsilon}\). This is a comparison of exact theorem scopes, not a claim about worldwide priority or the absence of further literature. The bound proved here is \(N^{33/(33k+35)+\epsilon}\), with the explicit saving \(2/((k+1)(33k+35))\) below \(1/(k+1)\). No common positive saving independent of \(k\), optimality among all methods, or fixed full-\(k\) exponent \(33/134\) is claimed. Fixed-shift estimates, target averages, unproved critical-cycle hypotheses and ordinary-prime statements are not used.
