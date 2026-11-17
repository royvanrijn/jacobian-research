# Higher nilpotence and the reciprocal `R21` admission frontier

> **Status and scope.** This note proves an all-degree
> Cayley-integrable family with \(N^2\ne0\), and it sharpens the polynomial
> admission gate for the reciprocal `R21` packet. The higher-index family is
> an explicit symplectic automorphism and is Moyal-flat, so it is a control,
> not a `DC_2` candidate. For `R21`, the polynomial cotangent graph and four
> stable unimodular denominator charts are proved. The affine-contact route
> is excluded, while eleven explicit tame polynomial shears trivialize the
> graph form through degree four. The factorization \(R=xS\) then excludes
> every finite completion in the whole fiber-preserving stabilization
> subgroup. On the necessary stable-mixed \(U_2\) chart, twenty-six further
> exact shears remove the graph defect through degree six and leave an
> explicit septic remainder. Euler homotopy then eliminates every subsequent
> homogeneous defect in the completed local ring. This proves formal Darboux
> triviality. A constant-Pfaffian dilation path has a polynomial vector field
> but a nonpolynomial time-one map. The exact \(b=0\) kernel has no
> polynomial slice, so a normalizer's target-\(b\) component cannot have a
> polynomial Poisson mate on its zero fiber. In particular it cannot be an
> elementary coordinate; this excludes the complete twenty-six-factor
> correction subgroup and every elementary \(b\)-moving repair. The condition
> is nevertheless attainable: an explicit all-degree tame family supplies
> no-slice coordinates beginning in degree three. The complete polynomial
> constant rings and generic time divisors then make the split signature
> intrinsic and exclude every member of this family. For the closest row,
> \(k=2\), an exact polynomial conjugacy exists on \(I\ne0\), but its
> Jacobian is \(-I/5184\); the failure is concentrated on one affine-
> modification divisor. A new degree-seven Bezout coordinate now crosses
> that divisor and passes the complete \(b=0\) two-form admission gate. Its
> first transverse form equation also has a polynomial, divergence-free
> solution: an apparent \(1/C\) pole cancels against a constant-ring
> invariant. The resulting normal vector is not locally nilpotent, so this
> does not yet give a finite four-dimensional completion. At second order,
> the canonical minimal jet leaves only three normal defects. Bezout solves
> their form equation, but its volume equation has a uniquely determined
> nonterminating series. The obstruction is affine-linear in every invariant
> shift; grading reduces the relevant freedom to \(I^2U K[S]\), and a
> nonzero leading coefficient excludes that entire module. Thus the displayed
> degree-seven coordinate is excluded at second order. The only compatible
> Hamiltonian shifts are ambient symplectic shears, and all other Bezout
> companions are triangular changes \(W\mapsto W-TF\), so the exclusion
> covers the full reciprocal Bezout ansatz. A direct
> compactification audit further shows that raw pole order, leading contact
> degeneration on a fixed projective boundary, the pole divisor of an
> unnormalized primitive, and discrepancy data from the top wedge cannot be
> polynomial-Darboux invariants. The exponent-two to exponent-three
> constant-ring change is exactly one ordinary Rees modification, so that
> algebra realizes the required crossing rather than obstructing it. The
> viable global replacement is the Hamiltonian Derksen algebra. Exact
> calculation excludes the Hamiltonians of \(P,Q,e\) from its locally
> nilpotent locus and retains \(R\). More strongly, factorial closure and a
> genus-one generic-fiber calculation prove that every locally nilpotent
> Hamiltonian in the full stable-independent ring \(K[x,y,z]\) lies in
> \(K[R]\). The \(e\)-degree filtration eliminates the remaining sector and
> proves that the complete Hamiltonian-LND locus is \(K[R]\). Thus `R21` is
> **not globally polynomial-Darboux admitted**. On the quantum side, the
> natural graph-normal centralizer in the ambient inverse-Jacobian \(A_3\)
> is exactly a transported \(A_2\), so its restriction is an isomorphism, not
> a `DC_2` counterexample. The same Hamiltonian-LND theorem excludes every
> strict filtered PBW Weyl frame. Any surviving bridge must therefore be
> essentially filtration-collapsing or non-PBW.

## 1. A regular-index-four Hessian--Cayley family

In Darboux coordinates

\[
 (q_1,q_2,p_1,p_2),\qquad \{p_i,q_j\}=\delta_{ij},
\]

take arbitrary polynomials \(g,h\in K[p_1]\), a scalar \(c\), and

\[
 \mathcal A=p_1q_2+\frac c2p_2^2+g(p_1)p_2+h(p_1).
\tag{1.1}
\]

For \(N=\Pi\operatorname{Hess}(\mathcal A)\), direct multiplication gives

\[
 N^2=
 \begin{pmatrix}
 0&0&0&c\\
 0&0&-c&0\\
 0&0&0&0\\
 0&0&0&0
 \end{pmatrix},\qquad
 N^3=
 \begin{pmatrix}
 0&0&c&0\\
 0&0&0&0\\
 0&0&0&0\\
 0&0&0&0
 \end{pmatrix},\qquad N^4=0.
\tag{1.2}
\]

Hence \(c\ne0\) is a genuine pointwise Jordan-index-four open. The Cayley
matrix

\[
 C=I+N+\frac12N^2+\frac14N^3
\]

has closed rows and integrates to

\[
\begin{aligned}
 Q_1={}&q_1-q_2-g'(p_1)p_2-h'(p_1)+\frac c4p_1+\frac c2p_2,\\
 Q_2={}&q_2-g(p_1)-\frac c2p_1-cp_2,\\
 P_1={}&p_1,\\
 P_2={}&p_1+p_2.
\end{aligned}
\tag{1.3}
\]

The checker verifies \(JF=C\), \(JF\Pi JF^{\mathsf T}=\Pi\), and
\(\det JF=1\). Equation (1.3) is triangular: \(p_1=P_1\),
\(p_2=P_2-P_1\), after which \(q_2\) and \(q_1\) are recovered successively.
It is therefore a polynomial symplectic automorphism for every \(g,h,c\).
All odd Moyal bidifferentials of orders at least three vanish as well. This
is the first exact \(N^2\ne0\) row in the Hessian optimization branch, but it
is not structurally hostile.

## 2. Why the triangular branch cannot make \(N^2\) vary

The same calculation gives an all-degree rigidity statement. Start with

\[
 \mathcal A=p_1q_2+p_2^2+f(p_1,p_2).
\tag{2.1}
\]

For the Cayley matrix, every row curl vanishes identically except

\[
 \frac14 f_{222}-\frac12 f_{122},\qquad -\frac12f_{222}.
\tag{2.2}
\]

Thus row closure forces \(f_{222}=f_{122}=0\). In characteristic zero,
\(f_{22}\) is constant, and \(f\) has the form

\[
 f(p_1,p_2)=\frac a2p_2^2+g(p_1)p_2+h(p_1).
\]

Consequently the nonzero \(N^2\) in this entire triangular regular-\([4]\)
packet is necessarily constant. This explains the sparse searches: adding
more monomials inside this chain cannot produce the desired nonlinear
higher-nilpotence symbol.

The useful next Hessian search must therefore change the nilpotent flag, not
only enlarge \(g\) or \(h\). In particular, the HC4 cross-link should import
a moving regular-\([4]\) flag or a length-three packet rather than another
fixed triangular chain.

## 3. Exact `R21` graph data

Work over

\[
 K=\mathbb Q[q]/(q^2-4q+6).
\]

Put

\[
 A=1+xy^2,\qquad h(A)=q+(4q-6)A,
\]

\[
 B=A^2z+y^3h(A),\qquad P=AB,\qquad Q=y+xB,
\]

and

\[
 R=\int_0^{x/A}\bigl(1-t(Q-Pt)^2\bigr)\,dt.
\tag{3.1}
\]

The cancellation congruence removes every denominator in (3.1). The exact
certificate rechecks

\[
 R\in K[x,y,z],\qquad
 \det\frac{\partial(P,Q,R)}{\partial(x,y,z)}=-1,
\tag{3.2}
\]

and the quartic inverse incidence

\[
 T-\frac{Q^2T^2}{2}+\frac{2PQT^3}{3}
 -\frac{P^2T^4}{4}-R=0,
\qquad
 \partial_T=1-T(Q-PT)^2.
\tag{3.3}
\]

On \(\mathbb A^3\times\mathbb A^1_e\), the cotangent graph form

\[
 \Omega_{21}=dP\wedge dQ+de\wedge dR
\tag{3.4}
\]

is polynomial and has determinant one. This proves the universal local
rank-two graph, but polynomial admission still asks for an automorphism
\(H\in\operatorname{Aut}(\mathbb A^4_K)\) with
\(H^*\Omega_{21}=\omega_{\mathrm{std}}\).

## 4. Stable reciprocal powers are not yet graph admission

The identity \(A=1+xy^2\) gives, for every \(n\ge1\),

\[
 A^n-xy^2(1+A+\cdots+A^{n-1})=1.
\tag{4.1}
\]

Therefore

\[
 U_n=
 \begin{pmatrix}
 A^n&x\\
 y^2(1+A+\cdots+A^{n-1})&1
 \end{pmatrix}\in\operatorname{SL}_2(K[x,y]).
\tag{4.2}
\]

Acting with \(U_n\) on the two stable variables gives a polynomial
automorphism with a displayed polynomial inverse. For \(n=2,3\), affine
translations make the first stable coordinates respectively

\[
 B+xe,\qquad P+xe.
\tag{4.3}
\]

Thus the isolated \(A^2\) and \(A^3\) reconstruction factors are not stable
coordinate obstructions. This is real progress over merely observing a
localized chart.

It is not, however, the graph--Darboux theorem. Pulling (3.4) through either
natural inverse chart in (4.3), in the displayed standard polarization,
leaves all six independent symplectic-pair
defects nonzero. The reciprocal powers can be cleared separately, while the
coupled \(P,Q,R\) graph form still carries mixed terms. Calling (4.2) an
`R21` admission proof would conflate stable denominator removal with a
polynomial symplectic trivialization.

## 5. The affine-contact shortcut is impossible

Let \(J\) be the Jacobian of \((P,Q,R)\). An affine-contact
symplectization with constant unit multiplier would in particular give a
nonzero constant skew form \(C\) on the target whose pullback

\[
 J^{\mathsf T}CJ
\tag{5.1}
\]

is constant on the source. Write the three coefficients of \(C\) as
\(a_i+qb_i\). After reducing by \(q^2-4q+6\), the positive-degree
coefficients of (5.1) give a \(262\times6\) rational linear system. Its
rank is six. Hence \(C=0\).

This excludes the cheapest contact lift exactly. It does not exclude a
nonlinear contact form, a nonconstant polarization, or a general polynomial
Darboux change.

## 6. Tame normalization through degree four

Put

\[
 (u_0,u_1,u_2,u_3)=(z,y,e,x)
\]

and express the graph map \(\mathcal K=(P,Q,e,R)\) in these variables. Its
tangent map at the origin is the identity. With

\[
 D=J_{\mathcal K}\Pi J_{\mathcal K}^{\mathsf T}-\Pi,
\]

the sole linear defect is \(D_{12}=u_0\). The elementary automorphism

\[
 u_1\longmapsto u_1-u_0u_3
\tag{6.1}
\]

has inverse \(u_1\mapsto u_1+u_0u_3\) and kills the complete defect through
degree two. Writing \(a=u_0,b=u_1,d=u_3\), the cubic upper-triangular
entries that remain are

\[
\begin{array}{c|cccccc}
ij&01&02&03&12&13&23\\ \hline
D^{(3)}_{ij}
&-3b^2d&(21-15q)ab^2&3bd^2&(5q-6)b^3&0&3b^2d.
\end{array}
\tag{6.2}
\]

A layer of ordinary coordinate shears cannot remove (6.2): if the degree
four corrections \(V_0,V_1\) omit respectively \(u_0,u_1\), their
linearized contribution to the paired entry is
\(-\partial_0V_0-\partial_1V_1=0\), while \(D^{(3)}_{01}\ne0\). Coupled
shears do work.

For a Hamiltonian \(H=\gamma b^r\ell(a,d)^n\), where \(\ell\) is linear, use
the exact polynomial shear

\[
 (a,d)\longmapsto(a+H_d,d-H_a).
\tag{6.3}
\]

The linear form \(\ell\) is invariant under (6.3), so changing \(H\) to
\(-H\) gives its polynomial inverse. Apply (6.3) for

\[
 -\frac14b^2(a+d)^3,\quad
 -\frac14b^2(a-d)^3,\quad
 \frac12b^2a^3,
\tag{6.4}
\]

and then apply the elementary shear

\[
 u_2\longmapsto u_2+(6-5q)ab^3.
\tag{6.5}
\]

These four factors kill (6.2). The checker kills the complete quartic
defect with five more shears of type (6.3), having Hamiltonians

\[
 -\frac1{45}b(a+d)^5,\quad
 \frac1{45}b(a-d)^5,\quad
 \frac1{90}b(a+2d)^5,\quad
 -\frac1{90}b(a-2d)^5,\quad
 -\frac23bd^5,
\tag{6.6}
\]

followed by

\[
 u_2\longmapsto u_2+(15q-20)a^2b^2d.
\tag{6.7}
\]

Thus (6.1) and the ten factors (6.4)--(6.7) are genuine polynomial
automorphisms, not a formal coordinate change, and the corrected graph is
symplectic through degree four. Its first remaining defect is degree five:

\[
\begin{array}{c|cccccc}
ij&01&02&03&12&13&23\\ \hline
D^{(5)}_{ij}
&-5a^2d^3
&(27-15q)a^3d^2+(330-111q)b^5
&0
&(45q-66)a^2bd^2
&-\frac52ad^4
&5a^2d^3.
\end{array}
\tag{6.8}
\]

This is a finite tame-jet theorem, not a global Darboux theorem.

## 7. Why this tame chain cannot close globally

All eleven corrections above preserve the stabilization fibration. More
generally, consider any source automorphism of the form

\[
 H(a,b,c,d)=\bigl(h(a,b,d),\lambda c+G(a,b,d)\bigr),
\tag{7.1}
\]

where \(h\in\operatorname{Aut}\mathbb A^3_K\) and
\(\lambda\in K^\times\). Since the first, second, and fourth graph outputs
are independent of \(c\), a symplectic \(\mathcal K\circ H\) would satisfy

\[
 d(P\circ h)\wedge d(Q\circ h)
 +d(\lambda c+G)\wedge d(R\circ h)
 =da\wedge db+dc\wedge dd.
\tag{7.2}
\]

Comparing the \(dc\)-coefficients forces

\[
 d(R\circ h)=\lambda^{-1}dd,
\]

so \(R\circ h=\lambda^{-1}d+\kappa\). In particular, \(R\) would be a
polynomial coordinate.

The exact cancellation formula instead has

\[
 R=xS(x,y,z),\qquad S(0,y,z)=1,qquad S\notin K.
\tag{7.3}
\]

Thus \(x\) and \(S\) are both nonunits and \(R\) is reducible. A polynomial
coordinate is irreducible, so (7.3) contradicts the consequence of (7.2).
Therefore no finite composition in the entire subgroup (7.1) can trivialize
the graph form. This explains why the exact jet corrections keep succeeding
without producing global admission.

The next constructive gate is consequently sharper: a source change must
mix the stable variable \(c\) into at least one base coordinate, or the
target polarization must change. The \(U_2,U_3\) charts from Section 4 do
exactly such stable mixing, so they are no longer optional denominator
conveniences; they are the smallest remaining admission charts.

## 8. The stable-mixed \(U_2\) chart survives

Take the exponent-two chart from Section 4:

\[
 Z=B+xe,\qquad E=y^2(1+A)z+e,
\tag{8.1}
\]

and tangent-normalize by

\[
 (u_0,u_1,u_2,u_3)=(Z,Y,E,X).
\]

The graph map again has identity tangent. Its complete linear defect is

\[
 D_{02}=-u_2,\qquad D_{03}=u_3,\qquad D_{12}=u_0.
\tag{8.2}
\]

Apply, in order, the elementary shears

\[
\begin{aligned}
u_0&\longmapsto u_0+u_2u_3,\\
u_2&\longmapsto u_2-\frac12u_0^2,\\
u_2&\longmapsto u_2+2u_0u_1^2,
\end{aligned}
\tag{8.3}
\]

followed by the three linear-Hamiltonian shears (6.4). Every factor has a
displayed polynomial inverse. The exact jet calculation kills all defects
through degree three. Writing again \(a=u_0,b=u_1,d=u_3\), the first
remaining upper-triangular entries are

\[
\begin{array}{c|cccccc}
ij&01&02&03&12&13&23\\ \hline
D^{(4)}_{ij}
&-abd^2
&(60-50q)b^4-6a^2bd
&\frac13ad^3
&7ab^2d
&-\frac13bd^3
&abd^2.
\end{array}
\tag{8.4}
\]

Unlike the chain in Section 6, (8.3) genuinely mixes the stable coordinate
into the base, so the reducibility theorem of Section 7 does not exclude a
finite continuation. In fact, (8.4) is removed by five more shears of type
(6.3), with Hamiltonians

\[
 \frac1{360}b(a+d)^5,\quad
 -\frac1{360}b(a-d)^5,\quad
 -\frac1{720}b(a+2d)^5,\quad
 \frac1{720}b(a-2d)^5,\quad
 \frac1{12}bd^5,
\tag{8.5}
\]

followed by

\[
 u_2\longmapsto u_2-\frac72a^2b^2d+(12-10q)b^5.
\tag{8.6}
\]

Thus twelve exact factors remove every defect through degree four. The
first remaining entries are the simpler quintic row

\[
\begin{array}{c|cccccc}
ij&01&02&03&12&13&23\\ \hline
D^{(5)}_{ij}
&-\frac13a^2d^3
&(12-10q)b^5+4ab^3d
&0
&a^2bd^2-b^4d
&-\frac16ad^4
&\frac13a^2d^3.
\end{array}
\tag{8.7}
\]

The quintic correction equation also has a sparse solution:

\[
 V^{(6)}=
 \left(
 -\frac19a^3d^3,\quad 0,\quad
 -\frac13a^3bd^2+ab^4d+\frac{6-5q}{3}b^6,\quad
 \frac1{12}a^2d^4
 \right).
\tag{8.8}
\]

Its \((a,d)\)-part is Hamiltonian with \(H=-a^3d^4/36\). The seven
linear-Hamiltonian factors have Hamiltonians

\[
\begin{gathered}
 \frac{13}{60480}(a+d)^7,\quad
 \frac{13}{60480}(a-d)^7,\quad
 -\frac1{15120}(a+2d)^7,\quad
 -\frac1{15120}(a-2d)^7,\\
 \frac1{181440}(a+3d)^7,\quad
 \frac1{181440}(a-3d)^7,\quad
 -\frac1{3240}a^7.
\end{gathered}
\tag{8.9}
\]

Their leading vector is the \((a,d)\)-part of (8.8). The remaining component
is the elementary shear

\[
 u_2\longmapsto u_2-\frac13a^3bd^2+ab^4d+\frac{6-5q}{3}b^6.
\tag{8.10}
\]

All eight new factors have exact polynomial inverses. After twenty factors
in total, every graph defect through degree five vanishes. The first
remaining row is sextic:

\[
\begin{array}{c|cccccc}
ij&01&02&03&12&13&23\\ \hline
D^{(6)}_{ij}
&-\frac34b^4(3a+d)^2
&\frac14b^4\!\left(9a^2+(200q-234)ad-27d^2\right)
&b^3\!\left((d+3a)^3-36a^3\right)
&\frac13a^3d^3+(12-10q)b^5d
&\frac94b^4(3a^2-6ad-d^2)
&\frac34b^4(3a+d)^2.
\end{array}
\tag{8.11}
\]

The sextic correction is also explicit. A sparse degree-seven vector is

\[
\begin{aligned}
V^{(7)}_0={}&-\frac94a^3b^4-\frac94a^2b^4d-\frac34ab^4d^2,\\
V^{(7)}_1={}&0,\\
V^{(7)}_2={}&-\frac1{12}a^4d^3-(12-10q)ab^5d-\frac{27}{20}b^5d^2,\\
V^{(7)}_3={}&-\frac94a^3b^4+\frac{27}{4}a^2b^4d
 +\frac94ab^4d^2+\frac14b^4d^3.
\end{aligned}
\tag{8.12}
\]

The \((a,d)\)-part is Hamiltonian with

\[
 H=b^4\left(
 \frac9{16}a^4-\frac94a^3d-\frac98a^2d^2-\frac14ad^3
 \right).
\tag{8.13}
\]

The checker decomposes (8.13) as

\[
 b^4\left[
 \frac{571}{384}a^4-\frac{49}{32}(a+d)^4
 +\frac{59}{64}(a+2d)^4-\frac{37}{96}(a+3d)^4
 +\frac9{128}(a+4d)^4
 \right],
\tag{8.14}
\]

so five exact Hamiltonian shears and one elementary \(u_2\)-shear integrate
(8.12). After twenty-six factors, every defect through degree six vanishes.
The next row is nonzero in degree seven:

\[
\begin{array}{c|cccccc}
ij&01&02&03&12&13&23\\ \hline
D^{(7)}_{ij}
&-2ab^3d^3
&(126-91q)b^6d+24a^2b^3d^2
&\frac32ab^2d^4
&-\frac{21}{2}ab^4d^2
&-\frac12b^3d^4
&2ab^3d^3.
\end{array}
\tag{8.15}
\]

This answers the structural question. The quintic and sextic corrections
both exist, and Section 9 proves that the calculation is an all-order formal
Darboux tower. Formal solvability at every fixed degree is therefore a
theorem, not evidence for a finite polynomial completion. The next useful
problem is not to kill (8.15) blindly. Section 10 proves a stronger global
restriction: the target-\(b\) component of any polynomial Darboux normalizer
can have no polynomial Poisson mate on its zero fiber. Hence it cannot be an
elementary coordinate, excluding both the entire correction subgroup used so
far and every elementary \(b\)-moving repair. A surviving normalizer needs a
genuinely non-elementary no-slice target coordinate, or a different target
polarization. Sections 11--13 construct and then globally exclude the first
infinite tame no-slice family, while identifying its exact dense-open match
and the affine-modification divisor where it fails. Until a coordinate with
the corrected modification or a universal obstruction is proved, the R21
census row remains a localized candidate.

## 9. The all-order formal recurrence

The tower has a closed general construction. Work in the completed local
ring \(K[[u_0,u_1,u_2,u_3]]\), let \(W=-\Pi\) be the standard form matrix,
and suppose the current Poisson defect first appears in homogeneous degree
\(m\):

\[
 D_m=\left[J\Pi J^{\mathsf T}-\Pi\right]_m.
\]

The corresponding form defect is

\[
 E_m=WD_mW.
\tag{9.1}
\]

It is a closed two-form whose coefficients have degree \(m\). Let

\[
 \mathcal E=\sum_{i=0}^3u_i\partial_{u_i}
\]

be the Euler field. Cartan's formula gives

\[
 d\!\left(\iota_{\mathcal E}E_m\right)
 =\mathcal L_{\mathcal E}E_m=(m+2)E_m.
\]

Therefore define

\[
 \alpha_m=\frac{\iota_{\mathcal E}E_m}{m+2},
 \qquad
 V_{m+1}=\Pi\alpha_m.
\tag{9.2}
\]

Then \(\iota_{V_{m+1}}\omega=-\alpha_m\), and hence

\[
 \mathcal L_{V_{m+1}}\omega=-E_m.
\tag{9.3}
\]

Precomposition by the formal change
\(\operatorname{id}+V_{m+1}+O(\mathfrak m^{m+2})\) kills the complete degree
\(m\) defect. Repeating (9.1)--(9.3) constructs an \(\mathfrak m\)-adically
convergent formal automorphism

\[
 \widehat H=\cdots\circ
 \left(\operatorname{id}+V_{m+1}\right)\circ\cdots
\]

with

\[
 \widehat H^*\Omega_{21}=\omega.
\tag{9.4}
\]

The checker verifies closedness and exact cancellation for the displayed
degree-six and degree-seven rows. Their radial correction-vector term counts
are respectively \((9,8,8,4)\) and \((4,2,4,1)\), with Euler denominators
eight and nine.

There is also a one-shot Moser presentation. Choose a formal primitive
\(d\beta=\Omega_{21}-\omega\), put

\[
 \Omega_s=(1-s)\omega+s\Omega_{21},
\]

and solve

\[
 \iota_{X_s}\Omega_s=-\beta.
\tag{9.5}
\]

The time-one formal flow of \(X_s\) is \(\widehat H\). In dimension four,
the Pfaffian of the path is quadratic in \(s\), and both endpoints have
Pfaffian one. Thus

\[
 \operatorname{pf}(\Omega_s)
 =1+s(1-s)\Delta
\tag{9.6}
\]

for a polynomial \(\Delta\) with zero constant term in tangent-identity
coordinates. Consequently \(\Omega_s^{-1}\), and hence \(X_s\), has a
canonical geometric-series expansion in \(s(1-s)\Delta\).

For the exact untruncated \(U_2\) pullback form, the checker finds

\[
 \Delta=-u_0u_3+O(\mathfrak m^3)\ne0.
\tag{9.7}
\]

After reduction modulo \(q^2-4q+6\), it has 124 terms and source degrees from
two through twenty-seven. Since the coefficient ring is a domain, the
canonical geometric series is genuinely infinite. Thus the straight-line
Moser construction does not terminate naively. In fact, for the standard
primitive \(\beta\), the Pfaffian-adjugate numerator of
\(-\Omega_s^{-1}\beta\) has degree at most one in \(s\), while the nonconstant
denominator in (9.6) has degree two. Since \(\beta\ne0\), at least one
component of the straight-line Moser field is genuinely nonpolynomial. This
does not exclude a different path or a special closed-form polynomial
time-one map for the rational flow.

Equations (9.2) and (9.5) eliminate the graph defect to every formal order.
They do **not** prove polynomial admission. A polynomial theorem now requires
one of two genuinely global outcomes:

1. resum the Moser series to a finite polynomial flow with a polynomial
   inverse, using the \(U_2\) reciprocal identities; or
2. construct an invariant under all stable-mixed polynomial changes which
   proves that the formal normalizer cannot be polynomial.

Thus further homogeneous elimination is algorithmic but no longer
decisive. The next research object is the formal normalizer itself—its
recurrence, denominator \(\Delta\), and possible rational or algebraic
closed form.

## 10. Constant-Pfaffian dilation and the target-\(b\) slice obstruction

Let \(K\) be the exact tangent-identity graph map on the \(U_2\) chart. There
is a reciprocal-compatible path which removes the straight-line Pfaffian
denominator:

\[
 K_s(u)=\frac{K(su)}s,
 \qquad K_0=\operatorname{id},\quad K_1=K.
\tag{10.1}
\]

Since \(J K_s(u)=J K(su)\),

\[
 \det J K_s=1,
 \qquad \Omega_s=K_s^*\omega,
 \qquad \operatorname{pf}(\Omega_s)=1
\tag{10.2}
\]

for every \(s\). The canonical trivializing field is therefore polynomial:

\[
 V_s=-(J K_s)^{-1}\partial_sK_s
     =-\operatorname{adj}(J K_s)\partial_sK_s.
\tag{10.3}
\]

The exact component term counts are \((301,142,303,149)\); their maximum
source degrees are \((54,40,53,43)\), and their maximum \(s\)-degrees are
\((52,38,51,41)\). Direct multiplication verifies

\[
 J K_s\,V_s+\partial_sK_s=0.
\tag{10.4}
\]

Consequently the formal flow \(H_s\) satisfies
\(K_s\circ H_s=\operatorname{id}\). This path still cannot give a polynomial
time-one map. Indeed, the incidence polynomial (3.3) is \(g(P,Q,T)-R\),
which is irreducible over \(K(P,Q,R)\): over \(K(P,Q)\) it is primitive and
linear in the independent variable \(R\). It has degree four in \(T\), and

\[
 y=Q-PT,\qquad
 x=\frac{T}{1-T(Q-PT)^2},\qquad
 z=\frac{P/A-y^3h(A)}{A^2},\qquad
 A=\frac1{1-T(Q-PT)^2}
\tag{10.5}
\]

reconstruct the complete source field. Thus \(K\) has generic degree four,
so \(H_1=K^{-1}\) is formal but not polynomial. The dilation path succeeds
in polynomializing the vector field, not its time-one flow.

There is also a path-independent obstruction to the subgroup used by the
twenty-six corrections. Write \((a,b,c,d)=(u_0,u_1,u_2,u_3)\), restrict the
exact form to \(b=0\), and define its primitive kernel derivation by

\[
 \iota_{\delta_0}(da\wedge dc\wedge dd)=\Omega_{21}|_{b=0}.
\tag{10.6}
\]

The checker gives

\[
\begin{aligned}
 \delta_0(a)={}&1-ad+cd^2-\frac13a^2d^3
 +\frac56acd^4-\frac12c^2d^5,\\
 \delta_0(c)={}&-a+cd,\\
 \delta_0(d)={}&\frac16d^4(a-cd).
\end{aligned}
\tag{10.7}
\]

The constant term in \(\delta_0(a)\) makes this kernel vector primitive. If a
polynomial Darboux normalizer preserved \(K[b]\), it would send
\(b\mapsto\lambda b+\mu\). Restriction from the fiber
\(b=-\mu/\lambda\) to the fiber \(b=0\) would conjugate \(\delta_0\), up to a
nonzero constant, to the kernel \(\partial_a\) of \(dc\wedge dd\). Hence
\(\delta_0\) would be locally nilpotent.

But \(d\mid\delta_0(d)\) and \(\delta_0(d)\ne0\). For a locally nilpotent
derivation on a domain, \(f\mid\delta(f)\) forces \(\delta(f)=0\), as follows
immediately by comparing locally-nilpotent degrees in
\(\delta(f)=f g\). This contradiction proves:

> **\(b\)-fiber obstruction.** No polynomial Darboux normalizer of the exact
> \(U_2\) form preserves the subring \(K[b]\).

The same kernel has the stronger property of admitting no polynomial slice.
Put

\[
 t=a-cd,\qquad
 I=d^4t^2-12d,\qquad
 J=d^3c-2,\qquad
 y=d^2t.
\tag{10.8}
\]

Exact differentiation gives

\[
 \delta_0(I)=0,\qquad
 \delta_0(J)=\frac12d^3tJ,\qquad
 \delta_0(y)=d^2.
\tag{10.9}
\]

Hence \(I=1,J=0\) is an invariant curve. Its coordinate ring is

\[
 K\!\left[y,\frac1{y^2-1}\right],
\qquad
 d=\frac{y^2-1}{12},
\tag{10.10}
\]

and the induced derivation is

\[
 \delta_0(y)=\frac{(y^2-1)^2}{144}.
\tag{10.11}
\]

If \(\delta_0(s)=1\) for a polynomial \(s\), restriction to this curve would
give a rational primitive of

\[
 \frac{144\,dy}{(y^2-1)^2}.
\tag{10.12}
\]

Its residues at \(y=1,-1\) are respectively \(-36,36\), whereas the
derivative of a rational function has zero residue at every pole. Therefore
\(\delta_0\) has no polynomial slice.

This strengthens the normalizer obstruction. Suppose \(H^*\Omega_{21}=\omega\)
and put \(f=H_1\), the component mapping to the target coordinate \(b\).
Restriction of \(H\) gives an isomorphism from \(f=0\) to the target fiber
\(b=0\), carrying the Hamiltonian kernel of
\(\omega|_{f=0}\) to \(\delta_0\), up to a nonzero constant. Consequently:

> **Target-\(b\) slice obstruction.** The Hamiltonian derivation of \(f\) on
> \(K[a,b,c,d]/(f)\) cannot have a polynomial slice. Equivalently, \(f\)
> cannot have a polynomial Poisson mate modulo \(f\).

In particular \(f\) cannot be an elementary coordinate
\(\lambda u_i+F(u_0,\ldots,\widehat{u_i},\ldots,u_3)\): its standard paired
coordinate is an immediate Hamiltonian slice. This excludes every elementary
\(b\)-moving shear as well as all twenty-six factors in Section 8. A positive
construction now requires a genuinely non-elementary target-\(b\) coordinate
with no Poisson mate, or a different target polarization. This is not yet an
obstruction to all stable-mixed polynomial automorphisms.

## 11. Sharp tame no-slice coordinates and the split R21 signature

The target-coordinate obstruction is sharp rather than vacuous. For every
integer \(k\geq1\), compose two elementary determinant-one maps by putting

\[
 G=b+ad,
 \qquad
 F_k=a+c^kG.
\tag{11.1}
\]

Then

\[
 (a,b,c,d)\longmapsto(F_k,G,c,d)
\tag{11.2}
\]

is a tame polynomial automorphism, with inverse

\[
 a=F_k-c^kG,
 \qquad
 b=G-(F_k-c^kG)d.
\tag{11.3}
\]

Let \(D_k\) be the Hamiltonian derivation of \(F_k\) on \(F_k=0\), written
in the fiber coordinates \((G,c,d)\). Exact bracket calculation gives

\[
\begin{aligned}
 D_k(G)&=kG^2c^{2k-1}-1,\\
 D_k(c)&=-Gc^{2k},\\
 D_k(d)&=-kGc^{k-1}.
\end{aligned}
\tag{11.4}
\]

Set

\[
 t=Gc^k,
 \qquad I=t^2-2c,
 \qquad J=1+c^kd.
\tag{11.5}
\]

Then

\[
 D_k(t)=-c^k,
 \qquad D_k(I)=0,
 \qquad D_k(J)=-ktc^{k-1}J.
\tag{11.6}
\]

On the invariant curve \(I=1,J=0\), one has

\[
 c=\frac{t^2-1}{2},
 \qquad
 K[I=1,J=0]=K\!\left[t,\frac1{t^2-1}\right],
\tag{11.7}
\]

and the induced time form is

\[
 -\frac{2^k\,dt}{(t^2-1)^k}.
\tag{11.8}
\]

Its residue at \(t=1\) is

\[
 (-1)^k2^{1-k}\binom{2k-2}{k-1},
\tag{11.9}
\]

and the residue at \(t=-1\) is its negative. Both are nonzero. Hence no
\(D_k\) has a polynomial slice. In particular

\[
 F_1=a+bc+acd
\tag{11.10}
\]

is a cubic tame coordinate with no polynomial Poisson mate modulo \(F_1\).
The necessary target-coordinate condition from Section 10 is therefore
already attainable in degree three.

There is a useful low-degree boundary on the other side. Every polynomial
coordinate \(f\) of degree at most two has a linear Poisson mate. Indeed,
write \(\nabla f=Hu+\ell\), with \(H\) constant symmetric. Since a coordinate
has nowhere-vanishing gradient, the affine linear system \(Hu+\ell=0\) is
inconsistent. Thus there is a constant vector \(v\) with
\(v^{\mathsf T}H=0\) and \(v^{\mathsf T}\ell=1\), so \(v(f)=1\). The
symplectic matrix is invertible, hence \(v\) is the Hamiltonian vector field
of a linear polynomial \(g\), and \(\{f,g\}=1\). Consequently any R21
target-\(b\) component has degree at least three.

The family (11.1) also identifies the next compatibility test. From (11.6),

\[
 \frac{D_k(J)}J
 =k\frac{D_k(c)}c,
\tag{11.11}
\]

while the target derivation in (10.9) satisfies

\[
 \frac{\delta_0(J)}J
 =3\frac{\delta_0(d)}d.
\tag{11.12}
\]

The target invariant curve has time-form pole order two, but its transverse
logarithmic weight is three. For (11.1), both values equal \(k\). Thus
\(k=2\) matches the tangential pole order and \(k=3\) matches the transverse
weight.

As a first exact probe of the missing intrinsic-divisor statement, the
checker solves

\[
 D(p)=m\frac{D(c)}c\,p
\tag{11.13}
\]

through total degree seven for \(m=1,2,3,4\), and performs the analogous
calculation with \(c,D\) replaced by \(d,\delta_0\). The eigenspace dimensions
for the \(k=2\) tame model are

\[
 (2,2,2,3),
\tag{11.14}
\]

while those for the R21 kernel are

\[
 (1,1,2,2).
\tag{11.15}
\]

The first extra irreducible semi-invariant visible in this range occurs at
weight two for the tame model, \(1+c^2d\), but at weight three for R21,
\(d^3c-2\). This pattern is not an all-degree prime classification. Indeed,
at degree eight the \(k=2\) model has the additional polynomial

\[
 P_8=(G^2c^3-2)(1+c^2d)+c,
\tag{11.16}
\]

which satisfies

\[
 D_2(P_8)=\frac{D_2(c)}cP_8.
\tag{11.17}
\]

It is primitive and linear in \(d\) over \(K[G,c]\), hence irreducible.
Thus the bounded weight pattern cannot by itself be promoted to an
obstruction. The complete ring of polynomial constants is the correct
invariant.

## 12. Polynomial constants and the Danielewski exponent

The complete constant ring can be computed uniformly for (11.1). Put

\[
\begin{aligned}
 M_k&=G^2c^{2k-1}-2,\\
 I_k&=cM_k=G^2c^{2k}-2c,\\
 J_k&=1+c^kd,\\
 S_k&=M_k^kJ_k,\qquad \sigma_k=(-2)^k,\\
 U_k&=\frac{S_k(S_k-\sigma_k)}{I_k^k}.
\end{aligned}
\tag{12.1}
\]

The last expression is polynomial. Indeed, \(S_k\) contains \(M_k^k\), and

\[
 S_k-\sigma_k
 =(M_k^k-\sigma_k)+M_k^kc^kd
\tag{12.2}
\]

is divisible by \(c^k\), since \(M_k+2\) is divisible by
\(c^{2k-1}\). Direct differentiation gives

\[
 D_k(I_k)=D_k(S_k)=D_k(U_k)=0
\tag{12.3}
\]

and

\[
 I_k^kU_k=S_k(S_k-\sigma_k).
\tag{12.4}
\]

These are all the polynomial constants. To see this, localize at \(c\) and
use

\[
 t=Gc^k,\qquad L=\frac{J_k}{c^k}.
\tag{12.5}
\]

Then

\[
 A_c=K\!\left[I_k,L,t,\frac1{t^2-I_k}\right],
\qquad
 D_k=-2^{-k}(t^2-I_k)^k\partial_t,
\tag{12.6}
\]

so \(\ker(D_k|_{A_c})=K[I_k,L]\). After localizing \(I_k=cM_k\), both
\(c\) and \(M_k\) are units and \(S_k=I_k^kL\). Hence every global
constant lies in \(K[I_k^{\pm1},S_k]\).

The only possible poles are along \(c=0\) and \(M_k=0\). Along those two
divisors, \(S_k-\sigma_k\) and \(S_k\), respectively, vanish to order \(k\).
For \(I_k^{-r}P(I_k,S_k)\), divide \(P\) by
\(S_k(S_k-\sigma_k)\). Regularity on the two branches forces the remainder
to be divisible by \(I_k^r\); the quotient term replaces
\(S_k(S_k-\sigma_k)/I_k^k\) by \(U_k\). Induction on \(r\) proves

\[
 \ker D_k
 \cong
 \frac{K[I_k,S_k,U_k]}
 {\,I_k^kU_k-S_k(S_k-\sigma_k)\,}.
\tag{12.7}
\]

The R21 derivation has the same presentation with exponent three. Retain
\(t=a-cd\) and put

\[
\begin{aligned}
 N&=d^3t^2-12,&
 I_0&=dN,\\
 J_0&=d^3c-2,&
 S_0&=N^3J_0,\\
 \sigma_0&=3456,&
 U_0&=\frac{S_0(S_0-3456)}{I_0^3}.
\end{aligned}
\tag{12.8}
\]

Here \(U_0\) is polynomial, and the exact checker verifies

\[
 \delta_0(I_0)=\delta_0(S_0)=\delta_0(U_0)=0,
\qquad
 I_0^3U_0=S_0(S_0-3456).
\tag{12.9}
\]

With \(y=d^2t\) and \(L_0=J_0/d^3\), localization gives

\[
 A_d=K\!\left[I_0,L_0,y,\frac1{y^2-I_0}\right],
\qquad
 \delta_0=\frac{(y^2-I_0)^2}{144}\partial_y.
\tag{12.10}
\]

The same two-branch intersection argument therefore proves

\[
 \ker\delta_0
 \cong
 \frac{K[I_0,S_0,U_0]}
 {\,I_0^3U_0-S_0(S_0-3456)\,}.
\tag{12.11}
\]

We use the standard elementary invariant of the Danielewski rings

\[
 R_{n,\sigma}
 =K[X,Z,Y]/(X^nY-Z(Z-\sigma)),
\qquad \sigma\ne0.
\tag{12.12}
\]

For \(n\ge2\), their Makar--Limanov invariant is \(K[X]\). Every nonzero
locally nilpotent derivation is a \(K[X]\)-multiple of

\[
 \partial(X)=0,\qquad
 \partial(Z)=X^n,\qquad
 \partial(Y)=2Z-\sigma,
\tag{12.13}
\]

and the resulting canonical plinth ideal in \(K[X]\) is \((X^n)\). A
filtration by the \(X\)-adic order proves the first assertion; after
localizing \(X\), a locally nilpotent derivation is a multiple of
\(\partial_Z\), and preservation of \(Y=Z(Z-\sigma)/X^n\) forces that
multiple to be divisible by \(X^n\). Thus \(n\) is an isomorphism invariant.
For \(n=1\), the Makar--Limanov invariant is trivial, so it is also separated
from every \(n\ge2\) row.

Any polynomial conjugacy of fiber derivations induces an isomorphism of
their polynomial constant rings. Equations (12.7) and (12.11) therefore
force \(k=3\). But the generic time differential (11.8) has divisor

\[
 (2k-2)[\infty]-k[p_+]-k[p_-],
\tag{12.14}
\]

where \(p_\pm\) are the two punctures \(t^2=I_k\). The R21 time differential
in (12.10) has pole order two at its punctures, so conjugacy of the generic
fibers forces \(k=2\). The two requirements are incompatible:

> **Tame-family exclusion.** No coordinate
> \(F_k=a+c^k(b+ad)\), for any \(k\ge1\), can be the target-\(b\) component
> of a polynomial R21 Darboux normalizer.

This is the first global exclusion of an infinite genuinely non-elementary
no-slice coordinate family.

## 13. The exact dense-open \(k=2\) conjugacy

The incompatible exponents do not mean the \(k=2\) row was a poor guess.
It is exactly correct away from the invariant divisor \(I_0=0\). Define

\[
 G=\frac{t}{36},\qquad
 C=6d,\qquad
 W=-\frac1{36d^2}+\frac{I_0L_0}{864}.
\tag{13.1}
\]

Although written rationally, \(W\) is polynomial:

\[
 W=\frac d{864}\left(
 a^2cd^3-2a^2-2ac^2d^4+4acd+c^3d^5-2c^2d^2-12c
 \right).
\tag{13.2}
\]

Writing \(D_2\) for (11.4) at \(k=2\), exact differentiation gives

\[
 \delta_0(G,C,W)
 =-\frac1{36}D_2(G,C,W).
\tag{13.3}
\]

Moreover,

\[
 I_2(G,C)=I_0,
\qquad
 L_2(G,C,W)=\frac{I_0L_0}{864},
\tag{13.4}
\]

and

\[
 \det\frac{\partial(G,C,W)}{\partial(a,c,d)}
 =-\frac{I_0}{5184}.
\tag{13.5}
\]

After inverting \(I_0\), equations (13.1) have a polynomial inverse in the
localized rings: recover \(d=C/6\), \(t=36G\),
\(L_0=864L_2/I_0\), then \(c=L_0+2/d^3\) and \(a=t+cd\). Thus (13.1) is an
exact derivation conjugacy on \(I_0\ne0\). The Jacobian factor in (13.5) and
the constant-ring exponents prove that it cannot cross \(I_0=0\) as a
polynomial automorphism.

This identifies the obstruction geometrically. The \(k=2\) coordinate
matches the reciprocal dynamics on the dense open, but performs the wrong
affine modification over the two components of \(I_0=0\). A surviving
coordinate must retain pole order two on the generic leaf while changing
the Danielewski exponent of its polynomial constant ring from two to three.

## 14. A degree-seven coordinate crosses the fiber gate

The split signature can be realized globally. Start with standard
symplectic coordinates \((G,P,C,Q)\), put \(Q'=Q+G^2/2\), and define

\[
 A=\frac16GC^4,\qquad
 B=\frac13G^2C^3-1,\qquad
 R=\frac23G^3C^2,\qquad
 S=-(B+2).
\tag{14.1}
\]

The key identity is the two-term Bezout relation

\[
 RA+SB=1.
\tag{14.2}
\]

Consequently

\[
 F=AQ'-BP,\qquad W=RP+SQ'
\tag{14.3}
\]

complete to a polynomial coordinate system \((F,G,C,W)\). More precisely,

\[
 \det\frac{\partial(F,G,C,W)}{\partial(G,P,C,Q)}=1,
\tag{14.4}
\]

and the inverse is

\[
 P=AW+(B+2)F,\qquad
 Q=BW+RF-\frac12G^2.
\tag{14.5}
\]

The coordinate \(F\) has total degree seven. Its Hamiltonian row restricts
on \(F=0\) to

\[
 \{F,G\}=1-\frac13G^2C^3=-B,
 \qquad
 \{F,C\}=\frac16GC^4=A,
 \qquad
 \{F,W\}=-G.
\tag{14.6}
\]

Now use the target-fiber map

\[
 (a,c,d)=(G+CW,W,C),\qquad t=a-cd=G.
\tag{14.7}
\]

Its Jacobian is \(-1\). Substitution into the R21 kernel in Section 10 gives
exactly (14.6). More strongly, the restricted standard source form in these
coordinates is

\[
 \beta
 =G\,dG\wedge dC
  +A\,dG\wedge dW
  +B\,dC\wedge dW,
\tag{14.8}
\]

and the checker verifies

\[
 (G+CW,W,C)^*(\Omega_{21}|_{b=0})=\beta.
\tag{14.9}
\]

Thus (14.3) passes the complete hypersurface admission gate, not only the
constant-ring or characteristic-derivation tests. In particular, the
affine-modification mismatch of the \(k=2\) family was a defect of that
family, not an obstruction intrinsic to R21.

There is also a polynomial first transverse jet. Write

\[
 \mathscr D(G)=B,\qquad
 \mathscr D(C)=-A,\qquad
 \mathscr D(W)=G
\tag{14.10}
\]

and retain the rational invariants

\[
 I=C^4G^2-12C,\qquad L=W-\frac2{C^3}.
\tag{14.11}
\]

The normal volume equation is

\[
 \mathscr D(s)=-\frac13C^2G.
\tag{14.12}
\]

Its obvious solution \(-2/C\) is not polynomial. The constant-ring geometry
supplies exactly the missing cancellation:

\[
\begin{aligned}
 s
 &=-\frac2C-\frac{I^2L}{144}\\
 &=-\frac{C^2}{144}
 \left(C^6G^4W-2C^3G^4-24C^3G^2W+48G^2+144W\right).
\end{aligned}
\tag{14.13}
\]

This is polynomial and satisfies (14.12). Define

\[
\begin{aligned}
 n_d={}&\frac13GC^3+\frac16C^4Gs,\\
 n_c={}&-Gs,\\
 n_a={}&\frac13C^3GW-\frac23C^2G^2\\
 &\quad+\frac{s}{6}
 \left(C^4GW-2C^3G^2-6CG+6\right).
\end{aligned}
\tag{14.14}
\]

For the first jet

\[
 a=G+CW+Fn_a,\qquad b=F,\qquad
 c=W+Fn_c,\qquad d=C+Fn_d,
\tag{14.15}
\]

the pullback of the full four-dimensional R21 form agrees with the source
form at \(F=0\), including all normal--tangent coefficients. Moreover, in
the source fiber coordinates the normal vector is

\[
 V=(n_a-Wn_d-Cn_c)\partial_G+n_d\partial_C+n_c\partial_W,
\tag{14.16}
\]

and direct differentiation gives \(\operatorname{div}V=0\). This is the
first-order condition for (14.15) to arise from a polynomial
volume-preserving family of fiber automorphisms.

The canonical shortcut does not terminate: \(C\mid V(C)=n_d\) and
\(V(C)\ne0\), so the standard divisibility lemma for locally nilpotent
derivations proves that \(V\) is not locally nilpotent. Therefore
\(\exp(FV)\) is not a certified polynomial automorphism. Nevertheless its
second formal jet is polynomial, so the next obstruction can be computed
exactly.

## 15. The canonical second jet has a volume obstruction

Use the truncated formal flow of \(14.16\) through order \(F^2\). The
coefficient of \(F\) in the resulting form defect has no tangent--tangent
part. Its only entries are

\[
 dF\wedge
 \left(\rho_G\,dG+\rho_C\,dC+\rho_W\,dW\right),
\tag{15.1}
\]

where the three coefficients are independent of \(q\), have respectively
9, 11, and 7 terms, and have total degrees 22, 22, and 21. Put

\[
 \ell=(-S,R,0).
\]

The Bezout identity (14.2) says that \(\ell\) pairs to a unit with the
primitive kernel vector of \(14.8\). Since

\[
 \rho_GB-\rho_CA+\rho_WG=0,
\tag{15.2}
\]

the form equation has the polynomial solution

\[
 Z=\frac12
 \left(-\rho_WR,-\rho_WS,\rho_GR+\rho_CS\right).
\tag{15.3}
\]

Its component degrees are 26, 26, and 22. Every other solution differs from
\(15.3\) by \(\lambda\mathscr D\). The remaining condition that the second
fiber jet preserve volume is therefore

\[
 \mathscr D(\lambda)=-\operatorname{div}Z.
\tag{15.4}
\]

This equation has no polynomial solution for the canonical first jet
\(14.13\). Indeed, set

\[
 N=C^3G^2,\qquad J=C^3W.
\]

The exact divergence is

\[
 \operatorname{div}Z=\frac{C}{15552}H(N,J),
\tag{15.5}
\]

where

\[
\begin{aligned}
 H={}&JN^4-2N^4-42JN^3-24N^3+594JN^2-432N^2\\
    &-1728JN-2592N-12960J+46656.
\end{aligned}
\tag{15.6}
\]

Give the polynomial ring the grading

\[
 \operatorname{wt}(C)=2,\qquad
 \operatorname{wt}(G)=-3,\qquad
 \operatorname{wt}(W)=-6.
\tag{15.7}
\]

The derivation \(\mathscr D\) has weight three, while \(15.5\) has weight
two. Hence only the weight-minus-one part of \(\lambda\) can contribute.
Every monomial of that weight is \(GCN^iJ^j\), so necessarily

\[
 \lambda=GC\,p(N,J).
\tag{15.8}
\]

Direct differentiation gives

\[
 \mathscr D(GCp)=\frac C6\mathcal L(p),
\]

with

\[
 \mathcal L(p)
 =(N-6)p+N(N-12)p_N+3N(2-J)p_J.
\tag{15.9}
\]

Thus \(15.4\) is equivalent to

\[
 \mathcal L(p)=-\frac{H}{2592}.
\tag{15.10}
\]

Write \(p=\sum p_n(J)N^n\). Coefficient comparison uniquely determines the
\(p_n\). At degree five it gives the nonzero constant

\[
 p_5=-\frac{17}{5388768}.
\tag{15.11}
\]

Because the forcing in \(15.10\) has \(N\)-degree four, every later
coefficient is the nonzero constant

\[
 p_n=\frac{n}{6(2n+1)}p_{n-1}\qquad(n\ge6).
\tag{15.12}
\]

The unique series never terminates, proving the claim.

This is the first genuine obstruction after fiber admission. All polynomial
solutions of \(14.12\) are

\[
 s=s_{\min}+h,\qquad h\in\ker\mathscr D.
\tag{15.13}
\]

The invariant freedom can also be eliminated. First, an exact calculation
with an arbitrary polynomial \(h\), before imposing
\(\mathscr D(h)=0\), shows that the quadratic variation of
\(\operatorname{div}Z\) vanishes. Thus the second-order obstruction is
affine-linear in \(h\).

Under the grading (15.7),

\[
 \operatorname{wt}(I)=2,\qquad
 \operatorname{wt}(S)=0,\qquad
 \operatorname{wt}(U)=-6.
\tag{15.14}
\]

A weight-\(w\) invariant shift changes the divergence in weight \(w+4\).
Since (15.5) has weight two, only the weight-minus-two part of the constant
ring can cancel it. The presentation (12.11) gives

\[
 (\ker\mathscr D)_{-2}=I^2U\,K[S].
\tag{15.15}
\]

Indeed, a monomial \(I^aS^bU^c\) has weight \(2a-6c=-2\), hence
\(a=3c-1\); using \(I^3U=S(S-3456)\) reduces it to \(I^2U\) times a
polynomial in \(S\).

It remains to exclude this one free polynomial. Let
\(K(S)=k_mS^m+\cdots\) with \(k_m\ne0\), and put \(r=m+2\). The top
\(J^r\) coefficient contributed to \(p_5\) by \(I^2US^m\) is

\[
 c_r=
 \frac{
 (-1728)^r(r-1)(3r-4)(3r-2)(3r-1)
 (9r^2-24r-5)}
 {716636160}.
\tag{15.16}
\]

This is nonzero for every integer \(r\ge2\). All linear factors are
nonzero there; the quadratic factor equals \(-17\) at \(r=2\), equals
four at \(r=3\), and is strictly increasing thereafter. Consequently the
highest \(J\)-degree in \(p_5\) forces \(k_m=0\), a contradiction. Descending
on \(\deg K\) excludes every nonzero \(K\), while \(K=0\) leaves the
baseline value (15.11).

> **Second-order exclusion.** No polynomial first-jet choice
> \(s=s_{\min}+h\), \(h\in\ker\mathscr D\), admits a polynomial
> volume-compatible second jet. Therefore no polynomial R21 Darboux
> normalizer can have the degree-seven coordinate (14.3) as its target-\(b\)
> component.

This does not exclude R21 admission through a different target-\(b\)
coordinate or a different polarization. It does supply a reusable
candidate-screening invariant: after the fiber and first-jet gates, compute
the graded class of \(\operatorname{div}Z\) in
\(\operatorname{coker}\mathscr D\). A nonzero leading class excludes the
whole invariant family at second order without climbing an unbounded jet
ladder.

## 16. The complete Bezout ansatz and the general second-order gate

There are two apparent freedoms in (14.3), but neither produces a new
candidate. First replace \(Q'=Q+G^2/2\) by

\[
 Q_H=Q+H(G,C)
\]

and retain

\[
 F_H=AQ_H-BP,\qquad W_H=RP+SQ_H.
\tag{16.1}
\]

The Bezout identity still makes \((F_H,G,C,W_H)\) a determinant-one
polynomial coordinate system, with

\[
 P=AW_H+(B+2)F_H,\qquad
 Q=BW_H+RF_H-H(G,C).
\tag{16.2}
\]

An exact bracket calculation gives the particularly rigid formula

\[
 \{F_H,W_H\}|_{F_H=0}=-H_G.
\tag{16.3}
\]

The R21 fiber row (14.6) therefore forces

\[
 H_G=G,\qquad H=\frac12G^2+k(C).
\tag{16.4}
\]

But \(k(C)\) is only the ambient symplectic shear
\(Q\mapsto Q+k(C)\). Direct substitution in (16.2) shows that the complete
source two-form in coordinates \((F_H,G,C,W_H)\), including its
\(F_H\)-dependence, is identical to the form for \(k=0\). Hence the
second-order obstruction of Section 15 is unchanged.

Second, every Bezout solution of

\[
 R_TA+S_TB=1
\]

has the form

\[
 R_T=R+BT,\qquad S_T=S-AT
\tag{16.5}
\]

for a polynomial \(T(G,C)\). Its companion coordinate is merely

\[
 W_T=R_TP+S_TQ'=W-TF.
\tag{16.6}
\]

Thus it is a triangular change of companion coordinate preserving the same
target-\(b\) component \(F\). We obtain the stronger conclusion:

> **Bezout-family exclusion.** Every reciprocal-compatible coordinate in
> the ansatz (16.1), with every polynomial \(H\) satisfying the fiber gate
> and every Bezout companion (16.5), is excluded by the all-invariant
> second-order obstruction of Section 15.

The same conclusion holds for independent translations of both momenta.
Take arbitrary \(K,L\in K[G,C]\) and define

\[
 F=-BP+AQ+K(G,C),\qquad
 W=RP+SQ+L(G,C).
\tag{16.7}
\]

The inverse is

\[
\begin{aligned}
 P&=A(W-L)+(B+2)(F-K),\\
 Q&=B(W-L)+R(F-K).
\end{aligned}
\tag{16.8}
\]

Let \(\Omega_{K,L}\) be the standard source form in coordinates
\((F,G,C,W)\), and let \(\Omega_0\) denote the canonical form from
(14.5). Exact differentiation gives

\[
 \Omega_{K,L}-\Omega_0=-\frac{\mathcal E(K,L)}6\,dG\wedge dC,
\tag{16.9}
\]

where

\[
\begin{aligned}
\mathcal E(K,L)={}&
 C^4G L_C+2C^3G^2K_C-2C^3G^2L_G-4C^2G^3K_G\\
 &-6C^2G^2K+6G+6K_C+6L_G.
\end{aligned}
\tag{16.10}
\]

There are no other differing coefficients, including away from \(F=0\).
But the restricted fiber form agrees with R21 if and only if
\(\mathcal E(K,L)=0\). Hence fiber admission itself forces

\[
 \Omega_{K,L}=\Omega_0
\tag{16.11}
\]

as complete four-dimensional polynomial forms. Section 15 then excludes
every solution of the PDE (16.10), without needing to classify that
solution space.

> **Affine-momentum exclusion.** No R21 Darboux normalizer whose
> target-\(b\) coordinate and companion are affine-linear in the symplectic
> momenta \((P,Q)\), with coefficient row \((-B,A)\), can arise from any
> polynomial base translations \(K(G,C),L(G,C)\).

The mechanism is not special to the displayed formulas. Let
\(A=K[x_1,x_2,x_3]\), let \(\mu\) be its standard volume form, and suppose
a fiber-admission problem has restricted closed form

\[
 \beta=\iota_D\mu
\]

with primitive divergence-free characteristic derivation \(D\). After a
volume-compatible first jet, let the next normal form residual be a
one-form \(\rho\). Compatibility gives \(\rho(D)=0\). If

\[
 \iota_Z\beta=\rho
\tag{16.12}
\]

has one polynomial solution, every solution is \(Z+\lambda D\). Since
\(\operatorname{div}D=0\),

\[
 \operatorname{div}(Z+\lambda D)
 =\operatorname{div}Z+D(\lambda).
\tag{16.13}
\]

Therefore

\[
 \mathfrak o_2=[\operatorname{div}Z]\in A/D(A)
\tag{16.14}
\]

is independent of the chosen solution of the form equation. Its vanishing
is necessary and sufficient for a polynomial volume-compatible second
correction, once (16.12) is solvable. The remaining first-jet freedom acts
through \(\ker D\); grading \(\mathfrak o_2\) before solving is what reduced
the R21 calculation to (15.15)--(15.16).

This gives the symbol optimizer a structural score rather than another
bounded jet count: prefer fiber-admitted coordinates for which the leading
class \(\mathfrak o_2\) vanishes, and reject a whole invariant family as soon
as a graded component survives.

## 17. Fixed projective boundary data are not orbit invariants

There is a useful compactification calculation, but its first conclusion is
negative. Compactify the displayed affine coordinates
\((x,y,z,e)\) by \(\mathbb P^4\), with hyperplane at infinity
\(H_\infty\). If a closed polynomial two-form \(\eta\) has highest
coefficient degree \(m\), then its generic pole order along \(H_\infty\) is

\[
 m+3.                                                     \tag{17.1}
\]

Indeed, in a transverse coordinate \(t\), a degree-\(m\) coefficient
contributes at order at most \(t^{-m-3}\). Its top coefficient is the
radial contraction \(\iota_E\eta_m\) of the highest homogeneous part.
Since \(d\eta_m=0\), Cartan's formula gives

\[
 d(\iota_E\eta_m)=(m+2)\eta_m,
\tag{17.2}
\]

so this contraction cannot vanish when \(\eta_m\ne0\).

For \(\Omega_{21}\), in the order \((x,y,z,e)\), the six coefficient
degrees are

\[
\begin{array}{c|rrrrrr}
ij&01&02&03&12&13&23\\ \hline
\deg \Omega_{ij}&16&13&11&16&11&11.
\end{array}                                                \tag{17.3}
\]

The degree-sixteen part is

\[
 (\Omega_{21})_{16}
 =2x^5y^9z\,dy\wedge(3z\,dx+x\,dz).                       \tag{17.4}
\]

Thus this \(\mathbb P^4\) model has generic pole order nineteen. Its
highest homogeneous two-form has rank two, and its radial one-form satisfies
\(\alpha\wedge d\alpha=0\); the leading boundary distribution is therefore
integrable rather than contact. On the same compactification the standard
form has pole order three and its radial form gives the usual contact
structure on \(\mathbb P^3\).

Neither difference is invariant under polynomial changes of affine
coordinates. In standard Darboux coordinates \((a,b,c,d)\), for every
\(N\ge1\) the elementary automorphism

\[
 \phi_N(a,b,c,d)=(a,b+a^Nc,c,d)                            \tag{17.5}
\]

has Jacobian one and

\[
 \phi_N^*\omega_{\rm std}
 =\omega_{\rm std}+a^N\,da\wedge dc.                     \tag{17.6}
\]

This standard-orbit form has pole order \(N+3\) on the fixed \(\mathbb P^4\)
and a decomposable highest part. Taking \(N=16\) reproduces the coarse pole
order and leading-rank degeneration of (17.3)--(17.4) inside the standard
Darboux orbit itself.

> **Fixed-boundary no-go.** Pole order and leading contact degeneration on
> the coordinate \(\mathbb P^4\) compactification can exclude bounded-degree
> normalizers, but cannot obstruct the full
> \(\operatorname{Aut}(\mathbb A^4)\)-orbit.

The same warning applies more strongly to a primitive. There is already a
polynomial primitive

\[
 \lambda_{21}=P\,dQ+e\,dR,
 \qquad d\lambda_{21}=\Omega_{21}.                         \tag{17.7}
\]

Every other polynomial primitive is \(\lambda_{21}+df\), because polynomial
de Rham cohomology in degree one on affine space vanishes. Taking \(f\) of
arbitrarily high degree changes the projective pole profile without changing
the symplectic form. In particular, \(I=0\) is not an affine pole or residue
divisor of either \(\Omega_{21}\) or (17.7). It appears as the failure
divisor of the localized fiber conjugacy, which is different data. A useful
primitive invariant would first need a functorial minimization modulo exact
polynomial one-forms; no such minimization is proved here.

## 18. The top volume is blind, and the Rees crossing exists

The exact Pfaffian calculation gives

\[
 \frac12\Omega_{21}^{\wedge2}
 =dx\wedge dy\wedge dz\wedge de
 =\frac12\omega_{\rm std}^{\wedge2}.                      \tag{18.1}
\]

Consequently any compactification invariant which remembers only the
rational top form, its canonical divisor, or its discrepancies sees the same
volume section on both sides. This does not rule out a compactification
decorated by the complete meromorphic two-form, but it does rule out a
crepant-birational obstruction based only on \(K_X+D\) or on the top wedge.

The ordinary Rees algebra at the distinguished affine modification is also
explicit. Normalize the nonzero constant \(\sigma\), and write

\[
 R_n=K[I,S,U_n]/(I^nU_n-S(S-\sigma)).                      \tag{18.2}
\]

In \(R_2\), the center \(J=(I,U_2)\) is a complete intersection supported
at the two reduced points \(I=U_2=0\), \(S\in\{0,\sigma\}\). Hence its Rees
algebra is

\[
 \mathcal R_{R_2}(J)
 =R_2[X,Y]/(IY-U_2X).                                     \tag{18.3}
\]

On the \(X\ne0\) affine chart put \(V=Y/X=U_2/I\). Substitution gives

\[
 I^3V=S(S-\sigma),                                       \tag{18.4}
\]

so this chart is exactly \(R_3\). Thus the exponent-three `R21` constant
ring is obtained from the exponent-two tame ring by one standard affine
modification over \(I=0\).

> **Rees conclusion.** The unadorned Rees algebra does not obstruct the
> required crossing. It constructs it. The degree-seven coordinate of
> Section 14 realizes the resulting exponent-three fiber globally, and the
> failure occurs only when its symplectic normal neighborhood is extended to
> second order.

Accordingly the relevant object is a *decorated* Rees chart: (18.3) together
with the characteristic derivation, transverse volume, and the class
\(\mathfrak o_2\) of (16.14). Forgetting any of these last two layers loses
the known obstruction.

## 19. A global Hamiltonian Derksen invariant

There is an affine invariant which does not require choosing a projective
boundary. For a polynomial symplectic algebra \((A,\omega)\), define

\[
 \operatorname{HD}(A,\omega)
 =K[\,f\in A:X_f\text{ is locally nilpotent}\,],          \tag{19.1}
\]

where \(\iota_{X_f}\omega=df\). Also define the Hamiltonian
Makar--Limanov intersection

\[
 \operatorname{HML}(A,\omega)
 =\bigcap_{X_f\ {\rm LND}}\ker X_f.                       \tag{19.2}
\]

A polynomial symplectomorphism conjugates Hamiltonian derivations and
preserves local nilpotence. Therefore (19.1)--(19.2) are invariants of the
polynomial symplectic isomorphism class. For the standard form, the four
coordinate Hamiltonians are constant partial derivatives. They generate the
coordinate ring and their kernels intersect in \(K\), so

\[
 \operatorname{HD}(A,\omega_{\rm std})=A,
 \qquad
 \operatorname{HML}(A,\omega_{\rm std})=K.                \tag{19.3}
\]

Either a proper Hamiltonian Derksen algebra or a nontrivial Hamiltonian
Makar--Limanov intersection would therefore prove global non-admission.

The first exact `R21` screen is asymmetric. With the convention above,

\[
 X_R=\partial_e,                                          \tag{19.4}
\]

so every \(h(R)\) has locally nilpotent Hamiltonian
\(h'(R)\partial_e\). On the other hand the inverse-Jacobian calculation gives

\[
 x\mid X_P(x)\ne0,qquad
 x\mid X_Q(x)\ne0,qquad
 A\mid X_e(A)\ne0,qquad A=1+xy^2.                       \tag{19.5}
\]

For a locally nilpotent derivation on a domain,
\(g\mid D(g)\) forces \(D(g)=0\). Hence \(X_P,X_Q,X_e\) are not locally
nilpotent. The same calculation, together with the quartic inverse incidence,
proves the stronger target-subalgebra statement

\[
\begin{aligned}
 \{f\in K[P,Q,R]:X_f\text{ is LND}\}&=K[R],\\
 \{f\in K[e,R]:X_f\text{ is LND}\}&=K[R].                \tag{19.6}
\end{aligned}
\]

For the first row, \(x\mid X_f(x)\), so local nilpotence forces
\(X_f(x)=0\). If \(f_P,f_Q\) are not both zero, cancellation would require
\(-X_Q(x)/X_P(x)\in K(P,Q,R)\). In the quartic root coordinate \(T\), its
specialization at \(P=1,Q=0\) is

\[
 -\frac{3T(2T^3-5)}{4(T^3-4)}.                            \tag{19.7}
\]

But \(T^4-4T+4R\) is irreducible over \(K(R)\). If (19.7) equaled a base
scalar \(c(R)\), cross-multiplication and reduction by this quartic would
leave

\[
 -4cT^3-9T+24R+16c,                                      \tag{19.8}
\]

whose \(T\)-coefficient is nonzero. Thus \(f_P=f_Q=0\). The second row of
(19.6) follows directly because \(X_R(A)=0\) and every nonzero
\(e\)-derivative retains the last divisibility in (19.5).

This is not yet a Darboux obstruction. It proves only

\[
 K[R]\subseteq\operatorname{HD}(A,\Omega_{21}),           \tag{19.9}
\]

and excludes new Hamiltonian \(\mathbb G_a\)-directions in the full
three-coordinate target algebra \(K[P,Q,R]\) and in \(K[e,R]\). A
polynomial Darboux chart could use Hamiltonians far outside
\(K[P,Q,e,R]\), just as the degree-seven fiber coordinate lies outside the
elementary families.

There is nevertheless an all-source strengthening in the complete
\(e\)-independent sector. Put

\[
 C=K[x,y,z],\qquad R=x\Sigma,qquad \Sigma|_{x=0}=1.       \tag{19.10}
\]

If \(f\in C\), then \(X_f\) preserves \(C\) and \(X_f(R)=0\). If \(X_f\)
is locally nilpotent on \(C[e]\), its restriction \(D\) to \(C\) is locally
nilpotent. The kernel of an LND on a domain is factorially closed. Therefore
\(R=x\Sigma\in\ker D\) forces

\[
 D(x)=D(\Sigma)=0.                                       \tag{19.11}
\]

Extend scalars to \(K(x)\). A nonzero LND on the two-variable polynomial
ring \(K(x)[y,z]\), after algebraic closure, has kernel generated by a
coordinate. Hence \(\Sigma\in\ker D\) would make the normalization of every
generic irreducible component of \(\Sigma=c\) an affine line.

The exact generic fiber has different geometry. Since \(\Sigma-c\) is
quadratic in \(z\), its \(z\)-discriminant is

\[
 -\frac{x^3}{18}
 \left(6c(1+xy^2)^2-5xy^2-6\right).                     \tag{19.12}
\]

Writing \(u=xy^2\), the branch polynomial is

\[
 6cu^2+(12c-5)u+6(c-1),                                 \tag{19.13}
\]

whose discriminant is \(24c+25\) and whose constant term is
\(6(c-1)\). Over \(\overline{K(x,c)}\), (19.12) therefore has four distinct
roots in \(y\). The smooth projective normalization of the quadratic cover
has genus one. This contradicts the affine-line fibers forced by a nonzero
LND, so

\[
 D=0.                                                     \tag{19.14}
\]

It remains to recover \(f\). Since \(X_f\) vanishes on \(C\), it has the
form \(h\partial_e\), and contraction with \(\Omega_{21}\) gives

\[
 df=h\,dR.                                                \tag{19.15}
\]

The generic \(R\)-fiber is geometrically integral: the analogous branch
quadratic is

\[
 6ru^2+(12r-5x)u+6(r-x),                                 \tag{19.16}
\]

with discriminant \(x(24r+25x)\). Thus \(K(R)\) is relatively algebraically
closed in \(K(x,y,z)\), and (19.15) forces \(f\in K(R)\). Finally
\(R(x,0,0)=x\), so a rational function of \(R\) which is polynomial on
\(\mathbb A^3\) restricts to a polynomial in \(x\). Hence it belongs to
\(K[R]\).

> **Complete stable-independent classification.** For the exact `R21` form,
> \[
> \boxed{\{f\in K[x,y,z]:X_f\text{ is locally nilpotent}\}=K[R].}
> \tag{19.17}
> \]
> The forward inclusion is the genus-one rigidity argument above; the
> reverse inclusion is \(X_{h(R)}=h'(R)\partial_e\).

The positive-\(e\)-degree sector now descends to (19.17). Give
\(A=C[e]\) the filtration \(\deg_e(e)=1\), \(\deg_e(C)=0\). Since the form
is independent of \(e\), write, for \(g\in C\),

\[
 X_g=Y_g+H_g\partial_e,
 \qquad Z=X_e|_C,
 \qquad Z(R)=-1.                                        \tag{19.18}
\]

The Hamiltonian product rule gives the exact decomposition

\[
 X_{e^jg}
 =e^jY_g+e^jH_g\partial_e+j e^{j-1}gZ.                 \tag{19.19}
\]

The highest homogeneous component of a filtered locally nilpotent
derivation is locally nilpotent on the associated graded ring. This follows
directly by taking the highest \(e\)-degree term in every iterate. Let

\[
 f=\sum_{j=0}^m e^jf_j,
 \qquad f_m\ne0.                                        \tag{19.20}
\]

If \(m>0\) and \(Y_{f_m}\ne0\), the highest derivation component is
\(e^mY_{f_m}\). It fixes \(e\), so its local nilpotence would make
\(Y_{f_m}\) a nonzero LND of \(C\). Moreover \(Y_{f_m}(R)=0\). The proof of
(19.14) applies to every LND of \(C\) fixing \(R\), not merely to one already
known to be Hamiltonian. It excludes this possibility. Thus
\(Y_{f_m}=0\), and the differential argument (19.15)--(19.17) gives

\[
 f_m=h(R).                                               \tag{19.21}
\]

The derivation degree now drops to at most \(m-1\). Its homogeneous
degree-\(m-1\) component is

\[
 \delta=e^{m-1}\bigl(Y_{f_{m-1}}+m h(R)Z\bigr)
       +e^m h'(R)\partial_e.                             \tag{19.22}
\]

Again \(\delta\) is locally nilpotent. Since
\(\delta(e)=e^m h'(R)\) is divisible by \(e\), the LND divisibility lemma
forces \(h'(R)=0\). Hence \(h=c\in K\), and

\[
 \delta=e^{m-1}E,
 \qquad E=Y_{f_{m-1}}+mcZ.                              \tag{19.23}
\]

The variable \(e\) is fixed by (19.23), so local nilpotence of \(\delta\)
implies local nilpotence of \(E\) on \(C\). But

\[
 E(R)=-mc.                                               \tag{19.24}
\]

If \(c\ne0\), a scalar multiple of \(R\) is a slice for \(E\). The slice
theorem gives \(C=(\ker E)[R]\), so \((R)\) is prime. This contradicts the
nontrivial factorization \(R=x\Sigma\). Therefore \(c=0\), contradicting
\(f_m\ne0\). No \(m>0\) can occur.

Combining this descent with (19.17) gives the complete invariant:

> **Global Hamiltonian-LND classification.** For the exact `R21` form,
> \[
> \boxed{
> \{f\in K[x,y,z,e]:X_f\text{ is locally nilpotent}\}=K[R].
> }                                                       \tag{19.25}
> \]
> Consequently
> \[
> \operatorname{HD}(A,\Omega_{21})=K[R],
> \qquad
> \operatorname{HML}(A,\Omega_{21})=K[x,y,z].            \tag{19.26}
> \]

## 20. Global polynomial-Darboux obstruction

A polynomial symplectomorphism from \(\Omega_{21}\) to the standard form
would pull the four standard coordinate Hamiltonians back to locally
nilpotent Hamiltonians which generate \(A\). Equation (19.25) says that all
such Hamiltonians lie in the proper subalgebra \(K[R]\). Therefore

\[
 \boxed{
 (\mathbb A^4,\Omega_{21})
 \not\simeq_{\operatorname{Aut}(\mathbb A^4)}
 (\mathbb A^4,\omega_{\rm std}).
 }                                                       \tag{20.1}
\]

This settles polynomial Darboux admission for the reciprocal `R21` form. It
does **not** by itself prove \((DC_2)\): that conclusion still depends on a
separate quantization bridge which transfers the classical obstruction to
the proposed Weyl-algebra endomorphism. The compactification and Rees audit
remains useful diagnostically: it identifies why raw boundary profiles fail
and why the obstruction is instead detected by global additive-group
actions.

The proof uses two standard external facts about locally nilpotent
derivations: Rentschler's two-variable kernel theorem (after extending the
field \(K(x)\)) and the slice theorem. Every `R21`-specific polynomial,
discriminant, Hamiltonian field, factorization, and filtration coefficient
used above is included in the exact symbolic certificate.

## 21. Quantization bridge audit

The non-admission theorem removes `R21` from the standard pipeline which
first chooses a polynomial Darboux frame and then quantizes its four
coordinates in \(A_2\). There is nevertheless a canonical ambient
quantization. Let

\[
 C=K[x,y,z],
 \qquad
 \delta_i(P_j)=\delta_{ij},
 \qquad (P_1,P_2,P_3)=(P,Q,R),                          \tag{21.1}
\]

be the inverse-Jacobian derivations. Since the Jacobian determinant is a
unit, they are a polynomial frame of \(\operatorname{Der}_K(C)\), and

\[
 D(C)=C\langle\delta_P,\delta_Q,\delta_R\rangle\simeq A_3. \tag{21.2}
\]

The inverse-Jacobian endomorphism sends the target position generators to
\(P,Q,R\) and their target derivatives to (21.1). In the target \(A_3\),
take the graph-normal Weyl pair

\[
 a=\partial_P,
 \qquad b=P-\partial_Q,
 \qquad [a,b]=1.                                       \tag{21.3}
\]

Its transported pair is

\[
 a'=\delta_P,
 \qquad b'=P-\delta_Q.                                 \tag{21.4}
\]

The four operators

\[
 R,\qquad \delta_R,\qquad \delta_Q,\qquad Q-\delta_P  \tag{21.5}
\]

centralize (21.4). They form two commuting Weyl pairs:

\[
 [\delta_R,R]=1,
 \qquad
 [\delta_Q,Q-\delta_P]=1.                              \tag{21.6}
\]

These generators exhaust the centralizer. The geometric input is

\[
 \ker_C\delta_P\cap\ker_C\delta_Q=K[R].                \tag{21.7}
\]

Indeed an element of the source function field killed by both derivations is
algebraic over \(K(R)\). The geometric integrality of the generic \(R\)-fiber
from (19.16) makes \(K(R)\) relatively algebraically closed, and the section
\(R(x,0,0)=x\) cuts the polynomial intersection down to \(K[R]\).

For completeness, write an arbitrary differential operator uniquely as

\[
 T=\sum_\alpha c_\alpha
 \delta_P^{\alpha_P}\delta_Q^{\alpha_Q}\delta_R^{\alpha_R}.
                                                               \tag{21.8}
\]

Commutation with \(a'\) first gives \(\delta_P(c_\alpha)=0\). Comparing the
coefficient indexed by \(\beta\) in \([b',T]=0\) gives the triangular
recurrence

\[
 \delta_Q(c_\beta)
 =-(\beta_P+1)c_{\beta+e_P}.                            \tag{21.9}
\]

At maximal \(\delta_P\)-order, (21.7) puts the coefficient in \(K[R]\).
Descending in that order, (21.9) integrates only powers of the slice \(Q\),
with integration constants again in \(K[R]\). This is precisely the normal
ordering expansion in \(Q-\delta_P\). Therefore

\[
 \boxed{
 \operatorname{Cent}_{A_3}(\delta_P,P-\delta_Q)
 =K\langle R,\delta_R,\delta_Q,Q-\delta_P\rangle
 \simeq A_2.}                                          \tag{21.10}
\]

The right side is exactly the image of the target centralizer under the
inverse-Jacobian endomorphism. Consequently its restriction is an
isomorphism *onto* the transported copy (21.10). Composing it with any
identification of that copy with a fixed \(A_2\) gives an automorphism, not a
non-surjective endomorphism. The natural graph-centralizer reduction cannot
prove \((DC_2)\).

There is a second, filtration-level no-go. Let \(\mathcal B\) be a filtered
PBW quantization with

\[
 \operatorname{gr}\mathcal B\simeq(A,\Pi_{21})          \tag{21.11}
\]

as a Poisson algebra. Suppose four Weyl generators in \(\mathcal B\) form a
**strict filtered Weyl frame**, meaning that they satisfy the two-pair Weyl
relations and their principal symbols generate \(A\). Their inner
derivations are locally nilpotent. The highest-symbol component of a
filtered locally nilpotent derivation is locally nilpotent on the associated
graded ring, so every one of the four principal symbols has locally
nilpotent Hamiltonian field. Equation (19.25) puts all four in \(K[R]\),
which cannot generate \(A\). Hence

\[
 \boxed{
 \text{no strict filtered PBW Weyl frame quantizes }(A,\Pi_{21}).}
                                                               \tag{21.12}
\]

This is the exact scope of the bridge obstruction. It does not rule out an
abstract, unfiltered identification of a global `R21` quantization with
\(A_2\), because leading symbols of a filtration-collapsing generating set
need not generate the associated graded algebra. Such an identification
would have to be essentially non-PBW and cannot be the graph centralizer
(21.10). Constructing or excluding that exotic possibility is the remaining
`R21` route to \((DC_2)\).

## Reproduction

Run:

```bash
.venv/bin/python scripts/verify_dc2_higher_nilpotence_r21_frontier.py \
  --output \
  artifacts/generated-results/dc2_higher_nilpotence_r21_frontier.json
```

The 2026-08-11 repository audit replayed this command under Python 3.14.6
and SymPy 1.14.0. The extended generated JSON has SHA-256
`c1e95fbda0ce35ebe091cd17efd6615ecee39c0f2cb91313de5386aea5083aad`.
The status registry's checker-source field pins the checker itself at
SHA-256
`c92d33333a07def8714b60b3cd0eb591c3d0dd266059f205811629ee531a8028`.
The replay used the repository's `.python-version` and `requirements.txt`
locks; the generated result changed only by adding the audited boundary,
Rees, Hamiltonian-LND, and quantization-bridge certificate data.

The output records the exact higher-index matrices, the triangular inverse,
the two row-curl generators, the `R21` map and graph checks, the affine-contact
rank screen, the tame graph jet through degree four, and the four unimodular
stable charts, together with the fiber-preserving subgroup no-go theorem and
the stable-mixed \(U_2\) jet through degree six. It also checks closedness and
exact linearized cancellation for the Euler-homotopy recurrence on the
displayed degree-six and degree-seven defects, the exact polynomial field on
the constant-Pfaffian dilation path, and the target-\(b\) no-slice
obstruction.
It additionally certifies the degree-seven Bezout coordinate, its exact
R21 fiber two-form admission, and the polynomial divergence-free first
transverse jet (14.13)--(14.16). It also certifies the canonical second-order
Bezout correction and the all-degree recurrence (15.11)--(15.12) proving
that its volume equation has no polynomial solution. Finally, it verifies
affine-linearity in arbitrary invariant shifts and the leading coefficient
(15.16), excluding the complete \(I^2U K[S]\) freedom.
It also checks the general bracket formula (16.3), invariance under
\(Q\mapsto Q+k(C)\), and the complete Bezout-companion classification
(16.5)--(16.6), together with the affine-momentum identity
(16.9)--(16.11).
It finally computes the exact \(\mathbb P^4\) pole-degree profile and leading
rank-two form, verifies the arbitrary-pole elementary standard-orbit control,
checks equality of the top volume forms and the polynomial primitive, realizes
the exponent-two to exponent-three crossing as the \(I\)-chart of the Rees
algebra of \((I,U)\), and certifies the divisibility witnesses excluding the
Hamiltonians of \(P,Q,e\) from the locally nilpotent locus while retaining
\(X_R=\partial_e\). The quartic-root cancellation test strengthens this to
the exact Hamiltonian-LND locus \(K[R]\) inside \(K[P,Q,R]\), and separately
inside \(K[e,R]\). Finally, the quadratic-in-\(z\) discriminants (19.12) and
(19.16) certify the genus-one generic \(\Sigma\)-fiber and geometric
integrality of the generic \(R\)-fiber used to prove the complete
\(e\)-independent classification (19.17). The last certificate block records
the exact Hamiltonian decomposition (19.19), the two associated-graded
layers (19.21)--(19.24), the resulting complete locus (19.25), and the global
polynomial-Darboux obstruction (20.1). It also verifies the inverse-Jacobian
frame, graph-normal commutator, and two transported Weyl-pair relations used
in the centralizer exhaustion (21.1)--(21.10); the PBW exhaustion and strict
filtered no-go are the written algebraic arguments of Section 21.
