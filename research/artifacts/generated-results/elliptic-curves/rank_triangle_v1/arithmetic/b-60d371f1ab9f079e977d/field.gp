
default(realprecision,80);setrand(1);
main()={
 my(E,M,v,D,f,nf,S,n,u,L,Q,p,c,fac);
 E=ellinit([0,0,0,-1349483544544960365470065304689450113428497854288507,19100725967956561417739236786166254841772557169223614917419535705136370935894]);v=0;M=ellminimalmodel(E,&v);
 print("MINIMAL_MODEL|",vector(5,i,M[i]));print("TRANSFORM|",v);
 print("DISCRIMINANT|",M.disc);print("C4|",M.c4);
 D=factor(abs(M.disc));for(i=1,matsize(D)[1],if(!isprime(D[i,1]),error("unproved factor")));
 S=vector(matsize(D)[1],i,D[i,1]);addprimes(S);
 for(i=1,#S,print("FACTOR|",[D[i,1],D[i,2]]));
 f=x^3+M.b2*x^2+8*M.b4*x+16*M.b6;
 print("CUBIC|",vector(4,i,polcoef(f,i-1)));
 if(!polisirreducible(f),error("rational2torsion"));
 nf=nfinit([f,concat([2],S)]);if(nfcertify(nf)!=[],error("uncertified order"));
 print("FIELD_DISC|",nf.disc);print("FIELD_INDEX|",nf.index);print("SIGNATURE|",nf.sign);
 n=0;for(i=1,#S,p=S[i];L=elllocalred(M,p);Q=idealprimedec(nf,p);c=if(L[1]==1,D[i,2]%2==0,#Q-1);n+=c;print("LOCAL|",[p,D[i,2],L[1],L[2],#Q,c]);print("SPLIT|",p,"|",vector(#Q,j,[Q[j].e,Q[j].f])));
 Q=idealprimedec(nf,2);print("LOCAL2_DIM|",#Q);
 print("REAL_LOCAL_DIM|",if(nf.sign[1]==3,1,0));
 print("BK_OFFSET|",n+if(M.disc>0,2,1));print("ROOT_NUMBER|",ellrootno(M));
 print("CONDUCTOR|",ellglobalred(M)[1]);
 writebin("/home/royvanrijn/src/jacobian-research/research/artifacts/generated-results/elliptic-curves/rank_triangle_v1/arithmetic/b-60d371f1ab9f079e977d/nf.bin",nf);print("DONE_FIELD|1");
};
iferr(main(),e,print("FAIL|",e);quit(1));quit(0)
