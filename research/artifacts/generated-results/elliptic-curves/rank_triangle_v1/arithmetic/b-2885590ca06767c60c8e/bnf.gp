default(realprecision,80);setrand(1);
main()={my(nf,b);nf=read("/home/royvanrijn/src/jacobian-research/research/artifacts/generated-results/elliptic-curves/rank_triangle_v1/arithmetic/b-2885590ca06767c60c8e/nf.bin");b=bnfinit(nf,0);print("CYC|",b.cyc);print("DONE_BNF|1");};
iferr(main(),e,print("FAIL|",e);quit(1));quit(0)
