default(realprecision,80);default(parisizemax,8589934592);default(nbthreads,1);setrand(1);

main()={
  my(nf,b,g,C);
  nf=read("/home/royvanrijn/src/jacobian-research/research/artifacts/local/elliptic-curves/11952-rank-bound-audit-v1/nf.bin");
  if(nf.pol!=x^3+x^2-4060408041437129659386958883089829522172534180*x+69302996254333102097575929282499740768182420645699427883330979661600,error("cubic mismatch"));
  if(nf.disc!=30932024439246801322787810077880389223317162788400703970784702350579246728265890412316825281477068142364105,error("field discriminant mismatch"));
  C=nfcertify(nf);if(C!=[],error("maximal order not certified"));print("CERTIFY_ORDER|",C);
  print("STAGE|bnfinit");gettime();setdebug("bnf",1);
  b=bnfinit(nf,0);setdebug("bnf",0);print("BNFINIT_MS|",gettime());
  print("PROVISIONAL_CYCLIC_FACTORS|",b.cyc);
  g=sum(i=1,#b.cyc,b.cyc[i]%2==0);print("PROVISIONAL_2RANK|",g);
  writebin("/home/royvanrijn/src/jacobian-research/research/artifacts/local/elliptic-curves/11952-class-quotient-v1/provisional_bnf.bin",b);
  print("DONE|provisional");
};
iferr(main(),e,print("FAIL|",e);quit(1));quit(0)
