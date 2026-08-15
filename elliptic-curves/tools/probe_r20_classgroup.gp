\\ R20 cubic-field / Brumer--Kramer probe. Scratch driver for GitHub Actions.
\\ The class-group computation is GRH-conditional until BNFCERTIFY_QUOTIENT=1.

default(parisizemax, 6000000000);
\p 100

x = 'x;
a4 = -4437412060110743641525245114305;
a6 = 3586842216822165612930264910099076801587288127;
E = ellinit([1,1,1,a4,a6]);

\\ If z=4*x_tors, this monic polynomial defines the cubic subfield of Q(E[2]).
f0 = x^3 + 5*x^2 - 70998592961771898264403921828872*x \
     + 229557901876618599227536954246340915301586440144;

\\ Exact bad-prime valuations from the pinned global minimal model.
bad = [
  [2,17], [3,2], [5,9], [7,6], [13,2], [17,4], [31,2], [79,2],
  [1049,1], [71889448247,1], [40200713707633,1],
  [491007790268548705232623905732119,1]
];

print("R20_PROBE_BEGIN");
print("PARI_VERSION=", version());
print("TWO_DIVISION_MONIC=", f0);
print("TWO_DIVISION_IRREDUCIBLE=", polisirreducible(f0));
print("TWO_DIVISION_POLY_DISC=", poldisc(f0));

fred = polredabs(f0);
print("REDUCED_POLYNOMIAL=", fred);
print("REDUCED_POLYNOMIAL_DISC=", poldisc(fred));

gettime();
nf = nfinit(fred);
print("NFINIT_MS=", gettime());
print("FIELD_POLYNOMIAL=", nf.pol);
print("FIELD_SIGNATURE=", nf.sign);
print("FIELD_DISCRIMINANT=", nf.disc);
print("FIELD_INDEX=", nf.index);
print("FIELD_INTEGRAL_BASIS=", nf.zk);

curve_disc = E.disc;
print("CURVE_MINIMAL_MODEL=", [1,1,1,a4,a6]);
print("CURVE_DISCRIMINANT=", curve_disc);
print("CURVE_DISCRIMINANT_SIGN=", sign(curve_disc));
print("CURVE_GLOBAL_REDUCTION=", ellglobalred(E));

phi_m = 0;
additive_term = 0;
for(i=1, #bad,
  p = bad[i][1];
  vd = bad[i][2];
  lr = elllocalred(E,p);
  split_count = #idealprimedec(nf,p);
  print("LOCAL|p=",p,"|vdisc=",vd,"|f=",lr[1],"|kod=",lr[2],"|tamagawa=",lr[4],"|cubic_primes=",split_count);
  if(lr[1] == 1 && vd % 2 == 0, phi_m++);
  if(lr[1] > 1, additive_term += split_count - 1);
);
archimedean_u = if(curve_disc > 0, 2, 1);
bk_n = phi_m + additive_term;
print("BK_U=",archimedean_u);
print("BK_PHI_M_COUNT=",phi_m);
print("BK_ADDITIVE_TERM=",additive_term);
print("BK_N=",bk_n);
print("BK_CLASS_2RANK_LOWER_FROM_RANK20=",20-archimedean_u-bk_n);

\\ Flag 0 is sufficient for the class-group structure and materially faster.
print("BNFINIT_BEGIN");
gettime();
bnf = bnfinit(fred,0);
print("BNFINIT_MS=", gettime());
print("CLASS_NUMBER=", bnf.no);
print("CLASS_GROUP_CYCLIC_FACTORS=", bnf.cyc);
class2 = sum(i=1,#bnf.cyc, if(bnf.cyc[i] % 2 == 0, 1, 0));
print("CLASS_GROUP_2RANK=", class2);
print("BK_SELMER_UPPER_FROM_CLASSGROUP=", class2+archimedean_u+bk_n);
print("BK_RESIDUAL_ABOVE_KNOWN_RANK20=", class2+archimedean_u+bk_n-20);

\\ bnfcertify(...,1) rigorously proves that the true class group is a quotient
\\ of the computed group.  If class2=13, this upper bound plus the independent
\\ Brumer--Kramer lower bound 13 proves the actual 2-rank is exactly 13.
print("BNFCERTIFY_QUOTIENT_BEGIN");
gettime();
certq = bnfcertify(bnf,1);
print("BNFCERTIFY_QUOTIENT_MS=", gettime());
print("BNFCERTIFY_QUOTIENT=", certq);
print("R20_PROBE_END");
quit;
