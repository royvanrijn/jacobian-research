# A Dirichlet joint law for the stripped Hopf factors and a quintic-chart obstruction

All fields K below have characteristic zero. Use
\[
 \mathcal A=K[x,y,t]/(t^2+2xy-1),\qquad U=1+x,
\]
\[
 \Lambda(f)=\operatorname{CT}_x\frac12\int_{-1}^1
       f\left(x,\frac{1-t^2}{2x},t\right)dt.
                                                        \tag{1}
\]
Every integral denotes rational polynomial integration. A phase p is called pure when Lambda(p^m)=0 for every positive integer m. The following statements concern this specified functional and subalgebra; the global primitive degree-five existence question remains open.

## 1. The exact joint law

Define the regular polynomials
\[
 B=2y-(2+x)t^2=\frac{1-U^2t^2}{x},\qquad
 C=2y-t^2=\frac{1-Ut^2}{x}.                            \tag{2}
\]
Put
\[
 A_0=-C,\qquad D_0=C-B=Ut^2,\qquad E_0=1+B=UC.
                                                        \tag{3}
\]
The equality E_0=UC uses precisely t^2+2xy=1, and A_0+D_0+E_0=1. For a scalar a and integer n>=0 write (a)_n=a(a+1)...(a+n-1), with (a)_0=1.

**Theorem 1 (joint Dirichlet law).** For every i,j,k>=0,
\[
 \boxed{\Lambda(A_0^iD_0^jE_0^k)
       =\frac{(1/2)_i(1/2)_j(1/2)_k}{(3/2)_{i+j+k}}.}
                                                        \tag{4}
\]
Consequently, as a rational polynomial functional, (A_0,D_0,E_0) has the Dirichlet(1/2,1/2,1/2) law, equivalently the joint law of the squared coordinates of a uniform real unit two-sphere. Thus the native (B,C) functional equals the law of (Y^2-1,-X^2), with X,Y orthogonal real sphere coordinates. The statement is equality of ALL polynomial moments; no regular square roots of A_0,D_0,E_0 on the original quadric are required.

Here is a complete algebraic proof, including the zero-exponent boundaries. Set I=i+k. From (2)--(3),
\[
 \Lambda(A_0^iD_0^jE_0^k)
   =(-1)^i[x^I]U^{j+k}\frac12\int_{-1}^1
                      t^{2j}(1-Ut^2)^I\,dt.            \tag{5}
\]
There is a unique formal square root V=sqrt(1+x) in K[[x]] with constant term 1. It has V-1=x/2+O(x^2). Define the genuine polynomial
\[
 J_{I,j}(V)=\frac12\int_{-V}^{V}z^{2j}(1-z^2)^I dz.
\]
Its derivative vanishes to order at least I at V=1, so
\[
             (V-1)^{I+1}\mid J_{I,j}(V)-J_{I,j}(1).
                                                        \tag{6}
\]
For I=0 this is simply the order-one statement that a polynomial minus its value vanishes at that point, and remains sufficient.

Substitute z=tV in polynomial integration. The factor U^{j+k} t^{2j}dt becomes U^{k-1/2}z^{2j}dz. Thus (5) becomes
\[
          (-1)^i[x^I]U^{k-1/2}J_{I,j}(\sqrt U).
\]
By (6), replacing J(sqrt U) by J(1) leaves this coefficient unchanged. The factor U^{k-1/2} is a formal unit, so it does not lower the endpoint error order I+1. Therefore (5) equals
\[
          (-1)^i\binom{k-1/2}{i+k}J_{I,j}(1).           \tag{7}
\]
Splitting the half-integer factors above and below zero gives exactly
\[
 (-1)^i\binom{k-1/2}{i+k}
       =\frac{(1/2)_i(1/2)_k}{(i+k)!}.
\]
Also elementary integration by parts, or the beta integral, gives
\[
 J_{I,j}(1)=\frac{(1/2)_j I!}{(3/2)_{I+j}}.
\]
These two identities prove (4). They include I=0 and all i=0 or k=0 cases. Every coefficient in the argument is rational, so it is valid directly over every characteristic-zero K; no ordered-field or embedding premise is used.

For completeness the squared-sphere law in Theorem 1 follows by an elementary polar or Gaussian radial integral: the normalized moments of three coordinate squares are exactly the ratio in (4). Alternatively (4) defines the corresponding rational Dirichlet functional on the simplex. It is a classical beta law, not an additional conjectural purity theorem.

A useful closed form is
\[
 \boxed{\Lambda(B^mC^n)=(-1)^{m+n}
   \frac{(1/2)_n(m+n)!}{n!(3/2)_{m+n}}\quad(m,n\ge0).}
                                                        \tag{8}
\]
Indeed B=-(A_0+D_0), C=-A_0. Expand their powers and use (4); the rising-factorial Vandermonde sum is (1/2)_n(n+1)_m. This proves (8) at ALL orders. In particular Lambda(B)=-2/3, Lambda(C)=-1/3, Lambda(B^2)=8/15, Lambda(C^2)=1/5 and Lambda(BC)=4/15.

## 2. The normalized five-term quintic chart

Consider
\[
       p=BC+\alpha C^2+\beta B+\gamma C+\delta,
       \qquad \alpha,\beta,\gamma,\delta\in K.          \tag{9}
\]
The coefficient of BC is normalized to one; a nonzero general coefficient can be divided out without changing purity.

**Theorem 2 (finite obstruction in this chart).** There are no parameters in K, or any characteristic-zero extension, for which
\[
                  \Lambda(p^m)=0\quad(1\le m\le5).    \tag{10}
\]
In particular this chart has no pure phase. The conclusion uses a complete exact finite necessary-equation certificate and does not extrapolate a vanishing prefix into an all-orders assertion.

To describe the certificate transparently put u=A_0, v=D_0. Then B=-(u+v), C=-u and
\[
 p=a u^2+uv+b u+c v+d,
 \quad a=1+\alpha,\ b=-\beta-\gamma,\ c=-\beta,\ d=\delta.
                                                        \tag{11}
\]
This is an invertible affine parameter change. The first moment vanishes exactly when
\[
                d=d_*=-a/5-1/15-(b+c)/3.              \tag{12}
\]
Let n_m(a,b,c) denote the Dirichlet expectation of (a u^2+uv+b u+c v+d_*)^m. They are explicit rational polynomials. For example, the complete multinomial definition is
\[
 n_m=\sum_{r+s+t+h+l=m}
   \frac{m!}{r!s!t!h!l!}a^rb^tc^h d_*^l
   \frac{(1/2)_{2r+s+t}(1/2)_{s+h}}
        {(3/2)_{2r+2s+t+h}}.                          \tag{13}
\]
There is no omitted limiting or moment hypothesis in this finite formula.

Define integer polynomials
\[
 G_2=1575n_2,\quad G_3=3378375n_3,\quad
 G_4=861485625n_4,\quad G_5=245523403125n_5.
                                                        \tag{14}
\]
For orientation the first one is
\[
\begin{split}
 G_2={}&112a^2+240ab-120ac+140b^2-140bc+140c^2\\
       &+8a+20b+20c+8.
\end{split}
\]
The exact certificate file
[hopf_dirichlet_quintic_certificate.json](hopf_dirichlet_quintic_certificate.json)
contains EVERY coefficient of all four generators and of integer multipliers T_2,T_3,T_4,T_5 for the identity
\[
                 T_2G_2+T_3G_3+T_4G_4+T_5G_5=N\ne0.    \tag{15}
\]
The data key `nonzero_rhs` records N as an explicit positive 457-digit integer. Multiplier degrees are respectively 9,8,7,6, with 168,119,80,50 terms. Thus (15) is an ordinary rational unit-ideal certificate after division by N. The JSON coefficient lists are part of this exact finite proof data; their meaning is ordinary polynomial multiplication in K[a,b,c].

The standard-library verifier
[verify_hopf_dirichlet_quintic_obstruction.py](verify_hopf_dirichlet_quintic_obstruction.py)
reconstructs n_1,...,n_5 independently with Fraction arithmetic, verifies each scaling (14), and independently derives the same five moment polynomials by literal native Hopf Laurent multiplication. It then multiplies all integer polynomials in (15) and verifies that their sum has ONLY the constant coefficient N, with N>0. It also rejects an altered multiplier. This proof data needs no computer-algebra implementation or Groebner-basis trust for replay.

There is also a finite native route to Theorem 2 that does not assume the all-orders joint law. Substitute the literal Laurent expressions (2) into (11), with u=-C, v=C-B and d=d_*. Expand the powers of orders one through five and apply
\[
 \Lambda(x^h t^r)=
 \begin{cases}1/(r+1),&h=0\text{ and }r\text{ is even},\\
 0,&\text{otherwise},\end{cases}
 \qquad h\in\mathbb Z,\ r\ge0.
\]
This independently yields zero at order one and exactly the four scaled generator polynomials in (14), as the full literal coefficient reconstruction in the verifier checks. Together with (15), those five finite identities already prove the obstruction, independently of Theorem 1. The all-orders joint law remains a separate conclusion of the written endpoint argument.

If (10) held, (12)--(14) would give G_2=G_3=G_4=G_5=0. Identity (15) would then give N=0 in K, impossible in characteristic zero. This proves Theorem 2. The large certificate was generated by exact rational algebra; the independently verified identity, not the method used to discover it, is the finite mathematical implication.

## 3. Degree, primitivity and the precise excluded sector

The highest sphere-filtration symbols of B and C are -xt^2 and -t^2, of degrees three and two. Distinct monomials B^iC^j have distinct top x,t exponent pairs (i,2i+2j). Since K[x,t] injects into the associated graded domain K[x,y,t]/(t^2+2xy), B and C are algebraically independent and a polynomial in them has exact weighted degree with weights (B:3,C:2).

Thus (9) is exactly the normalized chart of least degree-five polynomials in K[B,C]: its only weight-five monomial is BC, while the other allowed monomials are C^2,B,C,1. Theorem 2 therefore excludes ALL degree-five phases in this particular polynomial subalgebra. It does not exclude arbitrary primitive quintics on the sphere or in the original Image problem.

For the original two-pair polynomial ring K[xi_1,xi_2,z_1,z_2], define
\[
 R=\xi_1z_1+\xi_2z_2,\quad Z=\xi_1z_2,\quad
 W=2\xi_2z_1,\quad T=\xi_1z_1-\xi_2z_2,
 \qquad T^2+2ZW=R^2.
\]
Every bilinear has bidegree (1,1). The original operator is
\(\mathcal E_2(\xi^\alpha z^\beta)=\partial_z^\alpha z^\beta\).
For balanced homogeneous H of degree d its output is constant and
\[
 \mathcal E_2(H)=(d+1)!\Lambda\bigl(H|_{R=1,\,Z=x,\,W=y,\,T=t}\bigr).
\]
To verify this equality it suffices to take H=xi_1^a xi_2^(d-a) z_1^b z_2^(d-b).
Both sides vanish for a!=b by the angular weight. For a=b, division of H by R^d gives 2^(-d)(1+t)^a(1-t)^(d-a), whose normalized interval integral is a!(d-a)!/(d+1)! by elementary polynomial integration. This equals the original contraction after multiplying by (d+1)!.

With
\[
 D_3=2WR^2-(2R+Z)T^2,\qquad E_2=2WR-T^2,
\]
the corresponding degree-five homogeneous lift is
\[
 F_5=D_3E_2+\alpha R E_2^2+\beta R^2D_3
                    +\gamma R^3E_2+\delta R^5.        \tag{16}
\]
At (xi_1,xi_2,z_1,z_2)=(1,1,1,-1), one has R=0,D_3=4,E_2=-4, so F_5=-16 for EVERY parameter choice. Hence every member of this normalized chart is radially primitive and has least balanced degree five. Nevertheless (10) fails for at least one of its first five powers. The factorial angular identity converts this directly into a nonzero original two-pair contraction. A nonproper-power condition cannot evade an obstruction to purity itself.

## 4. Classical comparison and limitations

The spherical Dirichlet density and beta moments are classical; the specific joint law (4) for the stripped Hopf factors is proved locally above. It is not an appeal to a general classification of complex polynomial moments for arbitrary simplex densities.

The usual Gaussian theorem for polynomials of Gaussian degree at most two does not apply directly to the chart (11). Homogenizing (11) on a two-sphere gives a quartic, even polynomial in THREE Gaussian coordinates, with a rank-one squared-variable term; it is not a quadratic in Gaussian coordinates. Likewise Liu--Sun, *The factorial conjecture and images of locally nilpotent derivations*, Bull. Aust. Math. Soc.101(2020), Theorems2.6--2.7, addresses HOMOGENEOUS two-variable factorial and mixed Gamma(1),Gamma(1/2) statements. Those source statements were inspected as context, not imported here to infer the three-half-parameter, nonhomogeneous chart. The finite identity (15) avoids any such unproved transfer.

This result leaves the full primitive Image degree-five problem open, including phases outside K[B,C]. The all-orders joint law follows from the endpoint proof in Section 1. The finite obstruction follows from the exact integer identity (15), whose complete coefficients are retained in the certificate. Additional bounded checks corroborate those two arguments without replacing either one.
