# Completed conic row factors: conic96, homogeneous96 and odd104

15 November 2026. Internal written refinement of the [conic endpoint proof](CONIC108_ODD126_GLOBAL130.md) and [weighted odd proof](WEIGHTED_ODD_PREFIX108.md). Retain the original Gaussian Wick functional, rational quotient operators, full bases, analytic exclusions and complete coefficient/frame covers.

| Scope | Sufficient pure prefix | Least-prefix interval |
|---|---:|---:|
| Every complex ternary polynomial of degree at most three | 130, unchanged | [20,130] |
| Odd ternary cubics | **104** | **[20,104]** |
| Homogeneous ternary cubics | **96** | **[14,96]** |
| Full sixteen-coefficient W-linear family | **96** | Not computed |
| Cubic conic exclusion | **96** | Not computed |

These are sufficient bounds. Leastness,20/14 sufficiency, expanded radical cofactors and independent analytic review remain open. Smooth-leading[20,94] and homogeneous-smooth[14,56] and both exact lower witnesses are unchanged. All final bounds are even, so parity supplies no one-step reduction.

The [checker](verify_conic_completed_row_factors.py), [certificate](conic_completed_row_factors_certificate.json), [receipt](conic_completed_row_factors_verification.json) and [source audit](conic_completed_row_factors_source_audit.json) verify auxiliary algebra. Universal all-degree module and path inductions, analytic nonzero zeroth coefficients and complete coefficient/frame coverage remain written inputs. There is no independent replay, external review, formal verification or expanded global radical certificate.

## 1. Factor the completed rows, keeping the actual seeds

The retained conic phase has the full structure

    W=A0+sum_(k=1)^p A_k*y^k+sum_(k=1)^q K^k*B_k*y^(-k),
    deg A0<=3, deg A_k,deg B_k<=3-k, deg K=d in{0,1,2}.

For q=3, B3 is a nonzero constant. Keep B=t*K_prime-2K/3,h=t^3*K and the complete logarithmic source

    S=K*(t*W_t-2W-t^2)+B*y*W_y,
    R_(j+1)=h*(R_j)_t-j*h_prime*R_j+t^2*B*y*(R_j)_y+S*R_j.

Its difference from the original connection is an exact logarithmic derivative with the original measure dy/y; no ordinary-volume divergence is added. Every derivative of h and the vertical field is retained.

The previously proved polynomial module requires the coefficient of y^n,n<0, to contain K^(-n). It holds for every unreduced actual seed. High cancellation at exponent n>=p divides by p*A_p and has witness exponent n-p>=0. Each negative output from this step contains the required K power because it arises from a negative phase coefficient. Low cancellation at n<-q divides by -q*K^q*B_q and has witness exponent n+q<0. Removing K^q from the input leaves K^(-n-q); after multiplication by any negative phase term the required power is again present. The t^2 Euler term also obeys this requirement. Thus each cancellation preserves the module, allowing inverses of the formal endpoint coefficients A_p,B_q.

Clear the already bounded A_p and B_q powers before specializing coefficients. This gives polynomial completed normal coefficients N_(k,j), with

    N_(k,j)=K^max(-k,0)*Nhat_(k,j),
    deg Nhat_(k,j)<=old_row_bound-d*max(-k,0).

This formal calculation precedes substitution of actual endpoints and K. Shared roots among these polynomials, and K(0)=0, do not add a coprimality open. Identically zero K is assigned to the already retained lower actual negative-order case. Nonzero endpoint constants are coefficient-field units.

Apply the division only after computing the actual seed columns, using the same selected rows and scalar denominators H_j. The row transformation is invertible over C(t), and its unit row factor is1. Each cofactor for a scalar relation on the selected rows has the same row factor; divide all of them by that factor. No connection matrix is transformed or differentiated. Thus no time-dependent basis derivative is missing. H_j and its valuation at least3j are unchanged. The original scalar relation and period exclusion therefore remain applicable to the same Wick series.

## 2. Completed row bounds for the retained endpoints

For p=q=3 the sharper retained old row bounds are3j+|k|,4j+max(-k,0),5j-k for d=0,1,2. After the completed-row division they become

    d0:3j+|k|; d1:4j; d2:5j-|k|.

On the six rows-3,...,2, every actual seed rank and row subset containing0 gives gates72,84,96. The last is5*21-9=96.

Apply the same subtraction d*max(-k,0) to each retained endpoint-degree bound. The full-rank numbers below also bound every smaller actual rank and every row choice containing0:

| Endpoint | d=0 | d=1 | d=2 |
|---|---:|---:|---:|
| p2,q3,deg A2=0 | 85 | 94 | 93 |
| p2,q3,deg A2=1 | 77 | 85 | 93 |
| p1,q3,deg A1=0 | 52 | 56 | 70 |
| p1,q3,deg A1=1 | 46 | 50 | 78 |
| p1,q3,deg A1=2 | 50 | 58 | 84 |
| p3,q2,deg B2=0 | 87 | 93 | 89 |
| p3,q2,deg B2=1 | 71 | 80 | 89 |

For p2,q3,d2, the old bound is7j-1-k and the new bound is7j-1-|k|. It uses only A2^(2j-1) normal clearing; the old padding by (K^3*B3)^max(0,2j-3) is unnecessary on the full polynomial module. For p3,q2 retain its positive Euler representative and B2^(2j-2) clearing; the completed-row lemma applies to its unchanged module as well. All source and denominator derivatives are retained in those representatives.

The regular positive-order-zero q3 case retains69; there is no need to improve it here. Every q<3 endpoint other than p3q2 and p2q2 retains a bound at most63. It remains to refine p2q2, whose older enlarged bound108 had prevented a uniform improvement.

## 3. The p2q2 module closes the remaining case

Use the same full structured phase, negative Euler representative and basis y^-2,y^-1,1,y. The source support is[-2,2]. A high cancellation cannot exit below-2, and a low cancellation cannot exit above1, so the high and low paths do not create new out-of-basis work for one another. At seed j the high and low depths are at most2j-1 and2j-2. Clear only

    E_j=A2^(2j-1)*B2^(2j-2), H_j=h^j*E_j.

The negative coefficient module supplies all K^2 factors needed by low divisions; no K denominator is cleared anew. Formal endpoint inversion and then clearing allow shared roots. After completed-row division, the following safe bounds suffice.

For constant K the source coefficients obey degree<=3-|n| and therefore seed coefficients obey3j-|n|. Each raw high or low step grows rational coefficient degree at infinity by at most twice its exponent displacement. Starting from |n|<=2j gives completed rational row degree<=5j+2|k| before endpoint clearing. Since both endpoints have degree at most1, E_j adds at most4j-3, so

    d0: deg Nhat_(k,j)<=9j-3+2|k|.

For linear K use weights(t,y)=(2,1). The phase weight is at most6 and seed cost8. A high whole step has cost at most6-2=4. A low whole step, after the module division by K^2, has cost at most6+2-4=4. Endpoint padding has cost at most2 and is covered by the same bound. Including both complete depths gives weighted completed numerator bound8j+4*(4j-3)=24j-12. After completed-row division,

    d1: deg Nhat_(k,j)<=floor((24j-12-k)/2)-max(-k,0).

For quadratic K use weights(1,1), phase bound3 and seed cost5. Both whole step costs are at most1, including endpoint padding. Thus

    d2: deg Nhat_(k,j)<=9j-3-k-2max(-k,0)=9j-3-|k|.

Every smaller rank and row subset containing0 gives gates86,93,74. The common H_j is unchanged by row division, remains polynomial, and has valuation at least3j, including shared endpoint/K/t roots.

The complete conic cover is now

    p3q3<=96; p2q3<=94; p1q3<=84; p0q3<=69;
    p3q2<=93; p2q2<=93; all other cases<=63.

Hence conic96 holds for cubic G, including arbitrary lower terms in this W-linear family. It is not an arbitrary-degree statement.

## 4. Reassemble the full coefficient and frame cover

Parent92 and rank-one95 give full W-linear max(96,92,95)=96. Homogeneous smooth56,cusp10,reducible10,nonisotropic-node72 and isotropic-node conic96 give homogeneous96.

The [weighted centered odd complements](WEIGHTED_ODD_PREFIX108.md) remain unchanged. Combine node100 with conic96 to give full odd node100; transverse98 with conic96 gives98; triangle78 with conic96 gives96. Odd cusp and tangent remain94 after parent92 and rank-one95-to94. The complete restricted cover is

| Leading type | Full odd sufficient prefix |
|---|---:|
| Degree at most two | 7 |
| Smooth cubic | 94 |
| Irreducible node | 100 |
| Irreducible cusp | 94 |
| Transverse line-smooth-conic | 98 |
| Tangent line-smooth-conic | 94 |
| Three independent lines | 96 |
| Three distinct concurrent lines | 104 |
| Repeated leading line(s) | 98 |

Its maximum104 gives odd104. Concurrent lines now control this odd bound; it is not a parity inference from even global130. The unrestricted cover remains quadratic7,smooth94,node129,cusp130,transverse127,tangent120,triangle118,concurrent104,repeated99, so global130 remains controlled by cusp130. No extra Gaussian frame hypothesis is imposed on global130.

## Evidence, failed probes and assurance

The checker verifies three generic p2q2 source strata and27 logarithmic gauge identities;162 completed-row identities in54 retained models;57 new exact seed identities in15 p2q2 and four fully structured p2q3 models; polynomial normal coefficients, scalar valuations and all actual-rank/row arithmetic. Shared endpoint/K/t roots and simple/double K zeros are included. The accepted run takes about17 seconds.

The [first failed input](conic_completed_row_factors_initial_input.py) and [receipt](conic_completed_row_factors_initial_failure.json) retained an arithmetic expectation95 for p2q3,d1; individual negative-row subtraction gives94, so this is a stronger value and does not change the final96 maximum. The [second failed input](conic_completed_row_factors_module_probe_input.py) and [receipt](conic_completed_row_factors_module_probe_failure.json), with [terminal log](conic_completed_row_factors_module_probe_failure.log), tried to apply the full module assertion to an older enlarged synthetic p2q3 model. That phase had only its y^-3 endpoint divisible by K^3: its y^-1,y^-2 coefficients were t^4+1,t^5+1. They lack K,K^2, so completed-row division is inapplicable. This does not invalidate that model's earlier coarser certificate. The new checker uses full structured p2q3 examples and preserves the counterexample to the weaker premise. Reopen this route if any actual Gaussian conic phase fails the full coefficient module, not merely an endpoint condition.

The universal formal-module, path and all-degree inductions and analytic exclusions remain independently unreviewed written proofs. Exact checks are auxiliary evidence, not formal certification of the complete analytic theorem. No fresh literature-coverage claim, external theorem or producer package is imported; the selected primary-source refresh remains retained. All earlier canonical proof bodies and receipts stay byte-identical. Least global[20,130],odd[20,104],homogeneous[14,96],20/14 sufficiency, expanded I130 radical cofactors, independent analytic review, exceptional coefficient schemes, global term minima and lower-power GVC remain open. Qualitative cubic coverage and both lower witnesses stay closed internally.
