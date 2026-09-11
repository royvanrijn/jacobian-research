default(realprecision,80);default(parisizemax,8589934592);default(nbthreads,1);setrand(1);

main()={
  my(b,g,C);
  b=read("/home/royvanrijn/src/jacobian-research/research/artifacts/local/elliptic-curves/11952-class-quotient-v1/provisional_bnf.bin");
  if(b.nf.pol!=x^3+x^2-4060408041437129659386958883089829522172534180*x+69302996254333102097575929282499740768182420645699427883330979661600,error("cubic mismatch"));
  if(b.nf.disc!=30932024439246801322787810077880389223317162788400703970784702350579246728265890412316825281477068142364105,error("field discriminant mismatch"));
  print("PROVISIONAL_CYCLIC_FACTORS|",b.cyc);
  g=sum(i=1,#b.cyc,b.cyc[i]%2==0);if(g>16,error("NOT_SUFFICIENT"));
  print("STAGE|bnfcertify_flag1");gettime();C=bnfcertify(b,1);
  print("CLASSGROUP_QUOTIENT_CERTIFIED|",C);if(C!=1,error("quotient certification failed"));
  print("CERTIFY_MS|",gettime());print("G_UPPER|",g);print("DONE|quotient");
};
iferr(main(),e,print("FAIL|",e);quit(1));quit(0)
