default(realprecision,80);setrand(1);
P=[2, 3, 5, 7, 13, 17, 19, 47, 83, 7568348148323, 1023034923031265640083839109496335411267837093185764838343605884505996669629389738604483673];for(i=1,#P,if(!isprime(P[i]),error("nonprime hint")));addprimes(P);

nf=read("/home/royvanrijn/src/jacobian-research/research/artifacts/local/elliptic-curves/11952-rank-bound-audit-v1/nf.bin");setdebug("bnf",3);gettime();
print("STAGE|bnfinit");b=bnfinit(nf,0);print("MS_BNFINIT|",gettime());
print("COMPUTED_CYCLIC_FACTORS|",b.cyc);writebin("/home/royvanrijn/src/jacobian-research/research/artifacts/local/elliptic-curves/11952-rank-bound-audit-v1/uncertified_bnf.bin",b);
print("STAGE|bnfcertify_quotient");C=bnfcertify(b,1);
print("CLASS_QUOTIENT_CERTIFIED|",C);print("MS_CERTIFY|",gettime());
print("CLASS_2RANK_UPPER|",sum(i=1,#b.cyc,b.cyc[i]%2==0));
print("DONE|class_group");quit
