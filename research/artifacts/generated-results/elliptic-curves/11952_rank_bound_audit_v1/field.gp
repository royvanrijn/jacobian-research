default(realprecision,80);setrand(1);
P=[2, 3, 5, 7, 13, 17, 19, 47, 83, 7568348148323, 1023034923031265640083839109496335411267837093185764838343605884505996669629389738604483673];for(i=1,#P,if(!isprime(P[i]),error("nonprime hint")));addprimes(P);

E=ellinit([0,0,0,-84196621147240320617047979399750704971769668763392,206937637967498573533732205060490773473330035432339348502618816917963546624]);v=0;M=ellminimalmodel(E,&v);
print("MINIMAL_MODEL|",vector(5,i,M[i]));print("CHANGE_TO_MINIMAL|",v);
print("MINIMAL_DISCRIMINANT|",M.disc);print("C4|",M.c4);
D=factor(M.disc);if(prod(i=1,matsize(D)[1],D[i,1]^D[i,2])!=M.disc,error("factorization"));
for(i=1,matsize(D)[1],if(!isprime(D[i,1]),error("nonprime factor"));print("DISC_FACTOR|",[D[i,1],D[i,2]]));
if(M.a1!=0 || M.a3!=0,error("unsupported model"));
f=x^3+M.a2*x^2+M.a4*x+M.a6;
print("CUBIC_COEFFICIENTS_ASCENDING|",vector(4,i,polcoef(f,i-1)));
print("IRREDUCIBLE_MOD_23|",polisirreducible(Mod(1,23)*f));
print("STAGE|nfinit");gettime();nf=nfinit([f,P]);
print("MS_NFINIT|",gettime());print("CERTIFY_ORDER|",nfcertify(nf));
print("FIELD_DISCRIMINANT|",nf.disc);print("FIELD_SIGNATURE|",nf.sign);
print("FIELD_INDEX|",nf.index);print("FIELD_BASIS|",nf.zk);
n=0;
for(i=1,matsize(D)[1],p=D[i,1];L=elllocalred(M,p);Q=idealprimedec(nf,p);mult=(L[1]==1);add=(L[1]>1);contribution=if(mult,D[i,2]%2==0,if(add,#Q-1,0));n+=contribution;print("LOCAL|",[p,D[i,2],L[1],L[2],#Q,contribution]);print("SPLIT|",p,"|",vector(#Q,j,[Q[j].e,Q[j].f])));
u=if(M.disc>0,2,1);print("BOUND_OFFSET_U_PLUS_N|",u+n);
print("ROOT_NUMBER|",ellrootno(M));
writebin("/home/royvanrijn/src/jacobian-research/research/artifacts/local/elliptic-curves/11952-rank-bound-audit-v1/nf.bin",nf);print("DONE|field");quit
