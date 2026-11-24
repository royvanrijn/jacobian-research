"""Independent certificate verifier. Standard library only; no project imports.
Reconstructs each row from the stated counting inequality, not solver matrices.
"""
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import json, math, sys, copy, gzip
R=range(1,7)
NAMES=tuple([f'{s}{i}' for s in 'abc' for i in (*R,0)])+('Fa','Fb','Fc','D')
# Coordinate i=0 denotes total alpha_s (the solver calls it s).
NAMES=tuple(n[:-1] if n.endswith('0') else n for n in NAMES)
POS={n:i for i,n in enumerate(NAMES)}
N=len(NAMES); D=POS['D']
def require(b,msg):
    if not b: raise ValueError(msg)
def var(n):
    a=[Q(0)]*N;a[POS[n]]=Q(1);return a
def add(*terms):
    out=[Q(0)]*N
    for c,v in terms:
        for i,x in enumerate(v):out[i]+=c*x
    return out
def coords(s,I):return add(*[(Q(1),var(f'{s}{i}')) for i in I])
def row(v,b=0):return (v,Q(b))
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def base(lam=Q(1)):
    L=add(*[(1,var(s)) for s in 'abc'])
    rr=[row(L,lam)]
    for s in 'abc':
        rr+=[row(add((1,coords(s,R)),(-1,var(s)))),row(add((7,var(s)),*[(i-7,var(f'{s}{i}')) for i in R]),1)]
    for x,y in (('a','b'),('a','c'),('b','c')):
        rr.append(row(add((1,var('D')),(-1,var(x)),(-1,var(y)))))
    for s in 'abc':
        for p in (2,3,5):rr.append(row(add((1,coords(s,[i for i in R if i%p==0])),(-1,var('F'+s)))))
    rr.append(row(add(*[(1,var('F'+s)) for s in 'abc'],(-4,L),(6,var('D')))))
    for x,y in (('b','a'),('c','b')):
        rr.append(row(add((1,var(x+'4')),(2,var(x+'5')),(-1,var(y+'4')),(-2,var(y+'5')))))
    return rr
def subset(I,empty=True):
    require(type(I) is list and (empty or I),'invalid subset')
    require(all(type(i) is int and i in R for i in I) and I==sorted(set(I)),'invalid subset indices')
    return I
def branches(n,mode):
    L=add(*[(1,var(s)) for s in 'abc'])
    if n['type']=='pairdev':
        require(set(n)=={'type','eps','children'},'pairdev keys')
        eps=Q(n.get('eps',Q(1,100)))
        pairs=[('a','b'),('a','c'),('b','c')]
        out=[]
        for i,j in pairs:
            out.append([row(add((1,coords(j,[3,4])),(-1,coords(i,[3,4]))),-eps)])
        close=[]
        for i,j in pairs:
            close.append(row(add((1,coords(i,[3,4])),(-1,coords(j,[3,4]))),eps))
            close.append(row(add((1,coords(j,[3,4])),(-1,coords(i,[3,4]))),eps))
        out.append(close)
        return out
    if n['type']=='geo':
        require(set(n)=={'type','I','J','K','children'},'geo keys')
        ss=[subset(n[k]) for k in ('I','J','K')]
        sig0=add(*[(1,coords(s,I)) for s,I in zip('abc',ss)])
        sig1=add(*[(i-1,var(s+str(i))) for s,I in zip('abc',ss) for i in I])
        return [[row(add((1,var('D')),(-1,L),(1,sig0)))],[row(add((1,var('D')),(-1,L),(-1,sig1)),-1)]]
    require(n['type']=='ss' and set(n)=={'type','A','pA','B','IB','kB','OC','children'},'ss keys')
    A,B=n['A'],n['B'];require(A in 'abc' and B in 'abc' and A!=B,'letters')
    C=next(s for s in 'abc' if s not in (A,B))
    pair=subset(n['pA'],False);I=subset(n['IB'],False)
    require(len(pair)==2 and (mode=='extended' or math.gcd(*pair)==1),'exponent pair')
    OC=subset(n['OC'])
    require((not OC or mode in ('value_lock','all','all_deg','extended')) and type(n['kB']) is int and n['kB'] in I,'OC/kB')
    p,q=(var(A+str(j)) for j in pair)
    fA=add((1,var(A)),(-1,coords(A,R if mode in ('liberated','combined','all','all_deg','extended') else pair)))
    fB=add((1,var(B)),(-1,coords(B,I)))
    w=add((7,var(B)),*[( -7 if i in I else i-7,var(B+str(i))) for i in R])
    rest=coords(A,[i for i in R if i not in pair]) if mode in ('liberated','combined','all','all_deg','extended') else [Q(0)]*N
    diag=row(add((1,var('D')),(-1,fA),(-1,fB),(-1,var(C))))
    out=[[diag]]
    fc=coords(C,OC)
    wc=add(*[(i,var(C+str(i))) for i in OC])
    value=add(*[(i,var(B+str(i))) for i in I],(-1,wc))
    legacy={(2,3):(2,3),(3,4):(3,4),(2,5):(5,),(3,5):(5,),(4,5):(5,)}
    for b in range(8):
        u,v,z=(b>>2)&1,(b>>1)&1,b&1
        kp=add((1,p),(1,q),(2,rest),(Q(u,2),add((2,p),(-1,w))),(Q(v,2),add((2,q),(-1,w))))
        lock=add((z,value)) if OC else add((1,coords(B,[i for i in I if i!=n['kB']])),(z,var(B+str(n['kB']))))
        def off(delta):return row(add((1,var('D')),(-1,fA),(-1,fB),(-Q(1,2),var(C)),(-Q(1,2),fc),(-Q(1,2),kp),(-Q(1,2),delta)))
        rr=[row(add(((-1 if u else 1),add((2,p),(-1,w))))),row(add(((-1 if v else 1),add((2,q),(-1,w))))),row(add(((-1 if z else 1),value if OC else var(B+str(n['kB']))))),off(lock)]
        for i in (() if OC else legacy.get(tuple(I),())):rr.append(off(add((2,var(B+str(i))))))
        if mode in ('all','all_deg','extended') and OC:rr.append(off(coords(B,I)))
        if mode in ('group23','combined','all','all_deg','extended'):
            for e in (range(2,7) if mode in ('all_deg','extended') else (2,3)):rr.append(off(add((2,coords(B,[i for i in I if i%e])))))
        out.append(rr)
    return out
def verify(cert):
    require(cert.get('schema')=='break213-v1','schema')
    mode=cert['mode'];require(mode in ('group23','liberated','combined','value_lock','all','all_deg','extended'),'mode')
    require(0<Q(cert['lambda'])<=2,'lambda')
    target=Q(cert['D_target']);require(target>0,'positive target')
    stats={'leaf':0,'geo':0,'ss':0,'pairdev':0};bounds=[];cache={}
    def visit(node,rr):
        kind=node['type'];require(kind in stats,'node type');stats[kind]+=1
        if kind=='leaf':
            require(set(node)=={'type','y','bound'},'leaf keys')
            require(type(node['y']) is dict and node['y'],'empty dual')
            acc=[Q(0)]*N;rhs=Q(0)
            for j,t in node['y'].items():
                require(j.isdigit() and str(int(j))==j and int(j)<len(rr),'row index')
                val=Q(t);require(val>0,'dual positivity')
                r,b=rr[int(j)];rhs+=val*b
                for i in range(N):acc[i]+=val*r[i]
            require(acc[D]==1 and all(v>=0 for i,v in enumerate(acc) if i!=D),'dual coefficients')
            require(rhs==Q(node['bound']) and rhs<target,'dual objective')
            bounds.append(rhs);return
        key=json.dumps({k:v for k,v in node.items() if k!='children'},sort_keys=True)
        if key not in cache: cache[key]=branches(node,mode)
        bb=cache[key]
        require(type(node['children']) is list and len(node['children'])==len(bb),'branch coverage')
        for child,extra in zip(node['children'],bb):visit(child,rr+extra)
    visit(cert['tree'],base(Q(cert['lambda'])))
    require(stats['leaf']==cert['leaves'],'leaf count')
    upper=max(bounds);require(upper==Q(cert['max_certified_D']),'maximum bound')
    return {'mode':mode,**stats,'upper':str(upper),'target':str(target),'gap_below_213_359':str(Q(213,359)-upper)}
def point_checks():
    # Separate direct formula: no solver or row generator used.
    X={'a':[30,0,43,43,0,0],'b':[30,0,43,43,0,0],'c':[0,43,43,29,0,0]}
    T={'a':120,'b':120,'c':119};results={}
    for mode in ('group23','liberated','combined'):
      best=Q(10)
      for A in 'abc':
       for B in 'abc':
        if A==B:continue
        C=next(s for s in 'abc' if s not in (A,B))
        for j,k in combinations(range(6),2):
         if math.gcd(j+1,k+1)!=1:continue
         p,q=Q(X[A][j],359),Q(X[A][k],359)
         fA=Q(T[A]-sum(X[A]) if mode in ('liberated','combined','all','all_deg','extended') else T[A]-X[A][j]-X[A][k],359)
         rest=Q(sum(X[A])-X[A][j]-X[A][k],359) if mode in ('liberated','combined','all','all_deg','extended') else Q(0)
         for m in range(1,64):
          I=[i for i in range(6) if m>>i&1];s=sum(X[B][i] for i in I)
          fB=Q(T[B]-s,359);w=Q(28+sum((i+1)*X[B][i] for i in range(6) if i not in I),359)
          kap=2*rest+p+q+Q(1,2)*(max(Q(0),2*p-w)+max(Q(0),2*q-w))
          dc=[Q(s,359)]
          if mode in ('group23','combined','all'):
           dc.extend(Q(2*sum(X[B][i] for i in I if (i+1)%d),359) for d in (2,3))
          old={(1,2):(1,2),(2,3):(2,3),(1,4):(4,),(2,4):(4,),(3,4):(4,)}
          dc.extend(Q(2*X[B][i],359) for i in old.get(tuple(I),()))
          rhs=fA+fB+max(Q(T[C],359),(Q(T[C],359)+kap+min(dc))/2)
          best=min(best,rhs)
      results[mode]=str(best)
    require(results['group23']=='201/359' and results['combined']=='184/359','point minima')
    # Independently check the whole-value-lock boundary, including its branch rows.
    value_rhs=Q(33+47+43,359)+Q(120-43+86,718)
    require(value_rhs==Q(409,718),'value-lock arithmetic')
    node={'type':'ss','A':'c','pA':[2,3],'B':'a','IB':[1,3],'kB':1,'OC':[4],'children':[]}
    vv=[Q(0)]*N
    for s in 'abc':
      vv[POS[s]]=Q(T[s],359)
      for i,x in enumerate(X[s],1):vv[POS[s+str(i)]]=Q(x,359)
    bb=branches(node,'value_lock')
    for value,expected in ((Q(213,359),False),(value_rhs,True)):
      vv[D]=value
      valid=any(all(dot(r,vv)<=b for r,b in br) for br in bb)
      require(valid==expected,'value-lock branch boundary')
    results['value_lock_selected']=str(value_rhs)
    return results
if __name__=='__main__':
    args=[a for a in sys.argv[1:] if not a.startswith('--')]
    path=Path(args[0]) if args else Path(__file__).with_name('cert_combined.json.gz')
    cert=json.loads(gzip.decompress(path.read_bytes()) if path.suffix=='.gz' else path.read_bytes());print(json.dumps(verify(cert),indent=2))
    print('independent point minima:',point_checks())
    if '--faults' in sys.argv:
      def leaf(n):return n if n['type']=='leaf' else leaf(n['children'][0])
      faults=[]
      x=copy.deepcopy(cert);l=leaf(x['tree']);key=next(iter(l['y']));l['y'][key]=str(Q(l['y'][key])+Q(1,1000));faults.append(x)
      x=copy.deepcopy(cert);x['tree']['children'].pop();faults.append(x)
      x=copy.deepcopy(cert);x['max_certified_D']='0';faults.append(x)
      for i,x in enumerate(faults,1):
       try:verify(x)
       except (ValueError,KeyError):print('fault',i,'REJECTED')
       else:raise ValueError('fault accepted')
    print('VERIFIED')
