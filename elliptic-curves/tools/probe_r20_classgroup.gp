\\ R20 cubic-field / Brumer--Kramer probe. Scratch driver for GitHub Actions.
\\ No claim is made unless the emitted invariants are independently pinned.

default(parisizemax, 6000000000);
\p 100

x = 'x;
a4 = -4437412060110743641525245114305;
a6 = 3586842216822165612930264910099076801587288127;
E = ellinit([1,1,1,a4,a6]);

\\ If z=4*x_tors, this monic polynomial defines the cubic subfield of Q(E[2]).
f0 = x^3 + 5*x^2 - 70998592961771898264403921828872*x \
     + 229557901876618599227536954246340915301586440144;

print("R20_PROBE_BEGIN");
print("PARI_VERSION=", version());
print("TWO_DIVISION_MONIC=", f0);
print("TWO_DIVISION_IRREDUCIBLE=", polisirreducible(f0));
print("TWO_DIVISION_POLY_DISC=", poldisc(f0));
print("TWO_DIVISION_POLY_DISC_FACTOR=", factor(abs(poldisc(f0))));

fred = polredabs(f0);
print("REDUCED_POLYNOMIAL=", fred);
print("REDUCED_POLYNOMIAL_DISC=", poldisc(fred));

gettime();
nf = nfinit(fred);
print("NFINIT_MS=", gettime());
print("FIELD_POLYNOMIAL=", nf.pol);
print("FIELD_SIGNATURE=", nf.sign);
print("FIELD_DISCRIMINANT=", nf.disc);
print("FIELD_DISCRIMINANT_FACTOR=", factor(abs(nf.disc)));
print("FIELD_INDEX=", nf.index);
print("FIELD_INTEGRAL_BASIS=", nf.zk);
print("FIELD_GALOIS_GROUP=", polgalois(fred));

curve_disc = E.disc;
print("CURVE_MINIMAL_MODEL=", [1,1,1,a4,a6]);
print("CURVE_DISCRIMINANT=", curve_disc);
print("CURVE_DISCRIMINANT_SIGN=", sign(curve_disc));
curve_factor = factor(abs(curve_disc));
print("CURVE_DISCRIMINANT_FACTOR=", curve_factor);
print("CURVE_GLOBAL_REDUCTION=", ellglobalred(E));

phi_m = 0;
additive_term = 0;
for(i=1, matsize(curve_factor)[1],
  p = curve_factor[i,1];
  vd = curve_factor[i,2];
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

print("BNFINIT_BEGIN");
gettime();
bnf = bnfinit(fred,1);
print("BNFINIT_MS=", gettime());
print("CLASS_NUMBER=", bnf.no);
print("CLASS_GROUP_CYCLIC_FACTORS=", bnf.cyc);
class2 = sum(i=1,#bnf.cyc, if(bnf.cyc[i] % 2 == 0, 1, 0));
print("CLASS_GROUP_2RANK=", class2);
print("REGULATOR=", bnf.reg);
print("TORSION_UNITS=", bnf.tu);
print("FUNDAMENTAL_UNITS=", bnf.fu);

narrow = bnfnarrow(bnf);
print("NARROW_CLASS_NUMBER=", narrow[1]);
print("NARROW_CLASS_GROUP_CYCLIC_FACTORS=", narrow[2]);
narrow2 = sum(i=1,#narrow[2], if(narrow[2][i] % 2 == 0, 1, 0));
print("NARROW_CLASS_GROUP_2RANK=", narrow2);
print("BK_SELMER_UPPER_FROM_CLASSGROUP=", class2+archimedean_u+bk_n);
print("BK_RESIDUAL_ABOVE_KNOWN_RANK20=", class2+archimedean_u+bk_n-20);

print("BNFCERTIFY_BEGIN");
gettime();
cert = bnfcertify(bnf);
print("BNFCERTIFY_MS=", gettime());
print("BNFCERTIFY_RESULT=", cert);
print("R20_PROBE_END");
quit;
