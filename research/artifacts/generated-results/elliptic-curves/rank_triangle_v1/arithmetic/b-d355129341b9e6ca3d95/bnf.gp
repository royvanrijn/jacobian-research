default(realprecision,80);setrand(1);
main()={my(nf,b);nf=read("/home/royvanrijn/src/jacobian-research/research/artifacts/generated-results/elliptic-curves/rank_triangle_v1/arithmetic/b-d355129341b9e6ca3d95/nf.bin");b=bnfinit(nf,0);print("CYC|",b.cyc);print("DONE_BNF|1");};
iferr(main(),e,print("FAIL|",e);quit(1));quit(0)
