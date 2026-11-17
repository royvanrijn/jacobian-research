# Cyclic-coordinate nilpotence adapter and smaller global targets

This note concerns the known seventeen-dimensional Rodrigues local algebra. It does not assert any global membership in the moment ideal. Its purpose is to compute the least local exponents of global coordinate functions, so that a global support certificate can use lower-degree targets.

## Exact adapter

Let \(J=\begin{pmatrix}0&0&1\\0&1&0\\1&0&0\end{pmatrix}\), so \(q(x)=x^TJx\). Extract the degree-one, degree-two and degree-three homogeneous parts of the actual formal phase through parameter order six. Write these as
\[
H_1=\langle v,x\rangle,\qquad H_2=x^TJBx,\qquad H_3.
\]
Form the cyclic frame \(E=(v,Bv,B^2v)\) and its determinant \(D\). In the local algebra, \(D\) is a unit with residue \(1/8\). Set
\[
p=\tfrac12\operatorname{tr}(B^2),\quad r=\det B,\quad
s_j=\langle v,B^jv\rangle.
\]
If \(u=E^TJx\), then \(x=J\operatorname{adj}(E)^Tu/D\). Expanding the cubic after the numerator substitution gives tensor coordinates
\[
\overline U_{ijk}=D^3U_{ijk}.
\]
The script treats \(U\) as a symmetric tensor: an ordinary polynomial coefficient is divided by its multinomial multiplicity.

The replay checks all of the following in the actual seventeen-dimensional algebra:

- the stated full Gram formula and \(\det G=-D^2\);
- \(B^3=pB+rI\);
- \(E\operatorname{adj}(E)=DI\);
- all three harmonic trace equations \(\sum_{i,j}G_{ij}\overline U_{ijk}=0\);
- substitution back of the entire cubic, recovering \(D^3H_3\);
- exact normalization to \(D=1/8\).

For the last check, define \(\lambda=(1/(8D))^{1/6}\) by its finite binomial expansion in the nilpotent maximal ideal. Normalized coordinates are
\[
p'=\lambda^2p,\quad r'=\lambda^3r,\quad
s_j'=\lambda^{j+2}s_j,\quad
U_{ijk}'=\lambda^{-i-j-k-2}D^{-3}\overline U_{ijk}.
\]
All targets except the final one are therefore units times their unnormalized counterparts. For the final target,
\[
U_{222}'+128s_2'
=\lambda^{-8}D^{-3}(\overline U_{222}+2Ds_2),
\]
because \(\lambda^{12}=1/(64D^2)\). This gives a polynomial, unit-equivalent local target without choosing a numerical normalization root.

## Exact local exponent profile

The least positive exponents killing the fourteen normalized target functions are

| Target | Least exponent |
|---|---:|
| \(p\) | 7 |
| \(r\) | 5 |
| \(s_0\) | 7 |
| \(s_1\) | 5 |
| \(U_{000}\) | 3 |
| \(U_{001}\) | 3 |
| \(U_{002}\) | 4 |
| \(U_{011}\) | 4 |
| \(U_{012}\) | 5 |
| \(U_{022}\) | 7 |
| \(U_{111}\) | 5 |
| \(U_{112}\) | 7 |
| \(U_{122}\) | 3 |
| \(T=U_{222}+128s_2\) | 4 |

Every immediately preceding power is checked to be nonzero. The same profile holds at all three normalized Rodrigues points, by scalar symmetry.

In `global_target_nilindices.json`, the key `U222` denotes the target \(T\), not \(U_{222}\) alone. The file contains both normalized vectors and unit-equivalent vectors used during reconstruction.

## Simple adaptations

Exact multiplication gives the following smaller local exponents:
\[
\begin{array}{c|c}
4p+5s_0&3\\
35s_1-34r&2\\
U_{002}-6U_{011}&2\\
U_{111}+2U_{012}&2\\
4U_{022}-13U_{112}&3\\
T+\dfrac{896}{75}p^2&2.
\end{array}
\]
Again these are least exponents. They are checked using normalized vectors in the actual local algebra, not just weighted initial forms.

Assume the separately supplied global support lemma: the common zero set, on the first-six-moment cyclic scheme, of
\[
p,s_0,s_1,U_{000},U_{002},U_{022},T
\]
consists of the three normalized Rodrigues points. Under this lemma an equivalent seven-target list is
\[
\begin{aligned}
g_1&=4p+5s_0,&g_2&=s_0,&g_3&=s_1,\\
g_4&=U_{000},&g_5&=U_{002}-6U_{011},\\
g_6&=4U_{022}-13U_{112},&g_7&=T+\frac{896}{75}p^2.
\end{aligned}
\]
Indeed \(g_1=g_2=g_3=0\) forces \(p=s_0=s_1=0\). The Gram relation then gives \(s_2^3=1/64\), hence \(s_2\ne0\). The harmonic traces reduce to
\[
U_{011}=-2U_{002},\qquad U_{112}=-2U_{022}.
\]
Thus \(g_5=13U_{002}\), \(g_6=30U_{022}\), and \(g_7=T\) on this locus. The two lists have exactly the same common zero set.

Consequently, with that support lemma, the global gate is equivalently
\[
\dim_{\mathbb C}A_{14}=51
\quad\Longleftrightarrow\quad
g_1^3,g_2^7,g_3^5,g_4^3,g_5^2,g_6^3,g_7^2\in I_{14}.
\]
The forward implication uses the computed local exponents at the three points. The reverse implication uses the common-zero-set lemma, followed by the already known local lengths seventeen. **None of these seven global membership identities is established by the local computation.**

## A cubic relation replacing the seventh power

Put \(K_1=4p+5s_0\), \(K_6=4U_{022}-13U_{112}\), and retain the other adapted targets above. A bounded exact rank test on the five normalized local vectors
\[
s_0^3,\quad s_1^2,\quad s_2K_1,\quad K_1s_0^2,\quad K_6
\]
gives rank four and the unique relation normalized to have first coefficient one:
\[
\begin{aligned}
\Psi={}&s_0^3
+\frac{55782484748580}{23118017560861}s_1^2
-\frac{1009306881744}{399965701745}s_2K_1\\
&-\frac{34813910270351}{13198868157585}K_1s_0^2
-\frac{119786979}{79993140349}K_6.
\end{aligned}
\]
This is zero in the actual normalized local algebra. The primitive integer coefficient vector, in the displayed five-term order, is
\[
(3814472897542065,\ 9204109983515700,\ -9625759731192528,\ -10061220068131439,\ -5712042093615).
\]
All five terms have character six under the twelve-component grading below. The relation consequently holds at all three Rodrigues points.

An equivalent global certificate is therefore
\[
K_1^3,\ \Psi,\ s_1^5,\ U_{000}^3,\ (U_{002}-6U_{011})^2,\ K_6^3,\
\left(T+\frac{896}{75}p^2\right)^2\ \in I_{14}.
\]
For the reverse support implication, at any common closed point the other required powers give \(K_1=s_1=K_6=0\). Hence \(\Psi=s_0^3\), so \(s_0=0\), and the earlier support argument applies. The seven target polynomials have ordinary degrees \((3,3,5,3,2,3,4)\), with maximum five. Their characters are \((6,6,3,6,4,6,8)\).

This alternative replaces a seventh-power membership obligation by a cubic polynomial membership obligation. Its global membership remains unproved. `global_cubic_relation.py` reads the exported normalized vectors and checks the exact relation, including its primitive integer form.

## A twelve-component rational grading

The normalized cyclic model has a \(\mathbb Z/12\mathbb Z\) grading. Give coordinates weights
\[
\deg p=2,\quad\deg r=3,\quad\deg s_j=j+2,\quad
\deg U_{ijk}=-i-j-k-2\pmod {12}.
\]
All structural relations are homogeneous: the Gram determinant has weight twelve, hence zero modulo twelve, while the three trace equations have weights \(0,-1,-2\). Moment \(\mu_n\) has weight \(n\).

One geometric verification is to take \(\lambda^{12}=1\), set \(g=\lambda^6I\in O_3\), and act on phases by \(F(x)\mapsto\lambda F(gx)\). Then \(D\) is unchanged, because \(\lambda^6\det g=\lambda^{24}=1\), while \(h\) is invariant under \(g=\pm I\). The displayed coordinate weights and moment weights follow. This defines a cyclic order-twelve action; the grading and its linear-algebra decomposition are defined over \(\mathbb Q\), so no root-of-unity field extension is needed for membership matrices.

The seven adapted targets above are homogeneous of weights
\[
(2,2,3,10,8,6,4),
\]
and their required powers have weights
\[
(6,2,3,6,4,6,8)\pmod {12}.
\]
For a homogeneous target of weight \(w\), a witness multiplier on a generator of weight \(v\) need only retain monomials of weight \(w-v\). This is an exact reduction of any chosen finite coefficient ansatz; it neither bounds the required witness degree nor proves that an identity exists.

At a single Rodrigues point, the order-four stabilizer acts on transverse parameters as \((a,b)\mapsto(-ia,-b)\). This explains the observed support classes modulo four. It should not be confused with full weighted homogeneity in weights \((3,2)\).

## Replay

Keep the original `rodrigues_local_length17_certificate` folder beside the new scripts. Run

```sh
python global_target_nilindices.py > global_target_nilindices.log 2>&1
python global_adapted_targets.py > global_adapted_targets.log 2>&1
python global_cubic_relation.py > global_cubic_relation.log 2>&1
```

The first script reconstructs the full adapter and exports exact local vectors. The second independently reads those vectors and the retained multiplication matrices and checks the adapted targets. The third uses the same exports for the bounded five-vector relation test. All use rational arithmetic throughout.
