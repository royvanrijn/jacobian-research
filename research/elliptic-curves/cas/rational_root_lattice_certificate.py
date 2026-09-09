"""Elementary rational-root certificates extracted from the fixed-cover replay.

No polynomial factorization. Callers must freeze their prime pool and limits.
"""
from sage.all import QQ, ZZ, matrix, vector

def evaluate(coeff,x,m):
    value=0
    for c in reversed(coeff):value=(value*x+int(c))%m
    return value
def gauss(a,b):
    a=vector(ZZ,a);b=vector(ZZ,b)
    while True:
        if b*b<a*a:a,b=b,a
        q=(2*(a*b)+a*a)//(2*(a*a))
        if q==0:return a,b
        b-=q*a
def check_reduced(a,b,M,r):
    assert abs(matrix(ZZ,[a,b]).det())==M
    assert all((x[0]-r*x[1])%M==0 for x in [a,b])
    assert a*a<=b*b and 2*abs(a*b)<=a*a
def root_certificate(f,pool):
    assert f[0]!=0 and f.leading_coefficient()!=0 and all(c in ZZ for c in f.list())
    coeff=list(map(int,f.list()));der=[i*c for i,c in enumerate(coeff)][1:]
    A=abs(ZZ(f[0]));B=abs(ZZ(f.leading_coefficient()));H=max(A,B)
    chosen=None;attempts=[]
    for p in pool:
        if coeff[-1]%p==0:attempts.append({'prime':p,'gate':'leading_coefficient_zero'});continue
        roots=[x for x in range(p) if evaluate(coeff,x,p)==0]
        simple=all(evaluate(der,x,p)!=0 for x in roots)
        attempts.append({'prime':p,'gate':'PASS' if simple else 'multiple_root','roots':roots})
        if simple:chosen=(p,roots);break
    assert chosen is not None,'No eligible prime in the unchanged pool; retain UNKNOWN.'
    p,roots=chosen;reports=[];rational=[]
    for root in roots:
        r=int(root);M=int(p)
        while M<=8*H*H:
            nextM=M*M
            r=(r-evaluate(coeff,r,nextM)*pow(evaluate(der,r,nextM),-1,nextM))%nextM
            M=nextM
        assert r%p==root and evaluate(coeff,r,M)==0
        a,b=gauss([M,0],[r,1]);check_reduced(a,b,M,r)
        record={'initial_root':root,'modulus':str(M),'lift':str(r),
                'reduced_basis':[[str(c) for c in v] for v in [a,b]]}
        if a*a>2*H*H:record['decision']='NO_BOUNDED_NUMERATOR_DENOMINATOR_VECTOR'
        else:
            assert a[1]!=0
            q=QQ(a[0])/a[1];is_root=bool(f(q)==0)
            record.update(decision='RATIONAL_ROOT' if is_root else 'ONLY_POSSIBLE_RATIO_IS_NOT_A_ROOT',candidate=str(q))
            if is_root:rational.append(str(q))
        reports.append(record)
    return {'prime':p,'prime_attempts':attempts,'all_roots_mod_p':roots,
            'numerator_bound':str(A),'denominator_bound':str(B),'common_bound':str(H),
            'lifts':reports,'rational_roots':rational}
