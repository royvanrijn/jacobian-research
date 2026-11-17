"""Exact root/covector controls and a Laurent integer identity1; stdlib only."""
from pathlib import Path
from collections import defaultdict
from itertools import permutations
import hashlib,json

HERE=Path(__file__).resolve().parent
BASE=HERE
INPUT=HERE/'involutory_kch_obstruction_inputs.json'
DATA=json.loads(INPUT.read_text())

class Ring:
    def __init__(self,names):
        self.names=names.split();self.zero=(0,)*len(self.names)
    def c(self,n):return {self.zero:n} if n else {}
    def v(self,name,power=1):
        e=list(self.zero);e[self.names.index(name)]=power
        return {tuple(e):1}
    def add(self,*ps):
        out={}
        for p in ps:
            for e,a in p.items():
                out[e]=out.get(e,0)+a
                if not out[e]:del out[e]
        return out
    def neg(self,p):return {e:-a for e,a in p.items()}
    def sub(self,p,q):return self.add(p,self.neg(q))
    def mul(self,*ps):
        out=self.c(1)
        for p in ps:
            tmp={}
            for e,a in out.items():
                for f,b in p.items():
                    g=tuple(i+j for i,j in zip(e,f))
                    tmp[g]=tmp.get(g,0)+a*b
                    if not tmp[g]:del tmp[g]
            out=tmp
        return out
    def dump(self,p):return [[list(e),a] for e,a in sorted(p.items())]

def raw_first_loop():
    R=Ring('r R x X a b');r,rr,x,X,a,b=[R.v(n) for n in R.names]
    c=R.c;add,sub,mul=R.add,R.sub,R.mul
    pairing=[[c(2),r,x],[rr,c(2),a],[X,b,c(2)]]
    def dot(f,u):return add(*(mul(f[i],pairing[i][j],u[j]) for i in range(3) for j in range(3)))
    def difference(u,v,factor):return [sub(x,mul(factor,y)) for x,y in zip(u,v)]
    roots={k:[c(int(i==j)) for j in range(3)] for i,k in enumerate([1,2,6])}
    forms={k:[dict(a) for a in u] for k,u in roots.items()}
    for o,i,j in [(2,6,7),(7,2,3),(3,1,0)]:
        roots[j]=difference(roots[i],roots[o],dot(forms[o],roots[i]))
        forms[j]=difference(forms[i],forms[o],dot(forms[i],roots[o]))
        assert dot(forms[j],roots[j])==c(2)
    k=sub(c(1),mul(a,b))
    ss=add(mul(k,rr),mul(a,X));tt=add(mul(k,r),mul(b,x))
    assert roots[3]==[c(0),k,b] and forms[3]==[c(0),k,a]
    assert roots[0]==difference(roots[1],roots[3],ss)
    assert forms[0]==difference(forms[1],forms[3],tt)
    h=sub(c(2),mul(ss,tt))
    assert dot(forms[0],roots[1])==h and dot(forms[1],roots[0])==h
    closing_root=difference(roots[1],roots[0],h)
    closing_form=difference(forms[1],forms[0],h)
    assert closing_root==[sub(c(1),h),mul(h,ss,k),mul(h,ss,b)]
    assert closing_form==[sub(c(1),h),mul(h,tt,k),mul(h,tt,a)]
    return {'raw_first_loop_steps':3,'raw_directed_parameters':6,
            'closing_root_and_covector_formulas':True}

def paired_trace():
    R=Ring('x X y Y z Z');x,X,y,Y,z,Z=[R.v(n) for n in R.names]
    c=R.c;add,sub,mul,neg=R.add,R.sub,R.mul,R.neg
    pairing=[[c(2),x,y],[X,c(2),z],[Y,Z,c(2)]]
    def dot(f,u):return add(*(mul(f[i],pairing[i][j],u[j]) for i in range(3) for j in range(3)))
    def difference(u,v,a):return [sub(x,mul(a,y)) for x,y in zip(u,v)]
    U,V,W=[[c(int(i==j)) for j in range(3)] for i in range(3)]
    roots={0:U,4:[neg(a) for a in U],6:V,7:V,8:W}
    forms={k:[dict(a) for a in v] for k,v in roots.items()}
    for o,a,b in [(8,4,5),(5,8,9),(6,9,10),(10,0,12),(4,10,11)]:
        roots[b]=difference(roots[a],roots[o],dot(forms[o],roots[a]))
        forms[b]=difference(forms[a],forms[o],dot(forms[a],roots[o]))
        assert dot(forms[b],roots[b])==c(2)
    h=sub(c(1),mul(y,Y));v=add(mul(y,X),mul(z,h));C=add(mul(Y,x),mul(Z,h))
    w=sub(mul(add(c(2),h),y),mul(x,v));Q=sub(mul(add(c(2),h),Y),mul(X,C))
    expected_roots={5:[c(-1),c(0),Y],9:[y,c(0),h],10:[y,neg(v),h],
        11:[sub(y,w),neg(v),h],12:[sub(c(1),mul(Q,y)),mul(Q,v),neg(mul(Q,h))]}
    expected_forms={5:[c(-1),c(0),y],9:[Y,c(0),h],10:[Y,neg(C),h],
        11:[sub(Y,Q),neg(C),h],12:[sub(c(1),mul(w,Y)),mul(w,C),neg(mul(w,h))]}
    assert all(roots[i]==u for i,u in expected_roots.items())
    assert all(forms[i]==f for i,f in expected_forms.items())
    return {'paired_trace_steps':5,'directed_pairing_variables':6,
            'root_and_covector_formula_checks':10,'all_propagated_pairings':2}

def certificate():
    R=Ring('h y Y v C w Q p q r s alpha beta')
    h,y,Y,v,C,w,Q,p,q,r,s,alpha,beta=[R.v(n) for n in R.names]
    ai,bi=R.v('alpha',-1),R.v('beta',-1)
    add,sub,mul,neg,c=R.add,R.sub,R.mul,R.neg,R.c
    F=[sub(add(h,mul(y,Y)),c(1)),
       mul(p,sub(y,w)),add(c(1),mul(p,v)),sub(neg(mul(p,h)),beta),
       mul(q,sub(Y,Q)),add(c(1),mul(q,C)),sub(neg(mul(q,h)),bi),
       sub(c(-1),mul(r,sub(c(1),mul(Q,y)))),sub(neg(mul(r,Q,v)),alpha),
       sub(c(-1),mul(s,sub(c(1),mul(w,Y)))),sub(neg(mul(s,w,C)),ai)]
    # Each expression is accompanied by its complete list of ideal multipliers.
    def generator(i):return (F[i],[c(int(i==j)) for j in range(len(F))])
    def scale(a,e):return mul(a,e[0]),[mul(a,b) for b in e[1]]
    def total(*es):return add(*(e[0] for e in es)),[add(*(e[1][j] for e in es)) for j in range(len(F))]
    f=list(map(generator,range(len(F))))
    d=sub(w,y);D=sub(Q,Y)
    d_expr=total(scale(v,f[1]),scale(d,f[2]));assert d_expr[0]==d
    D_expr=total(scale(C,f[4]),scale(D,f[5]));assert D_expr[0]==D
    e=total(scale(neg(mul(bi,h)),f[2]),scale(neg(mul(bi,v)),f[3]))
    E=total(scale(neg(mul(beta,h)),f[5]),scale(neg(mul(beta,C)),f[6]))
    assert e[0]==sub(v,mul(bi,h)) and E[0]==sub(C,mul(beta,h))
    A=total(f[7],scale(neg(r),f[0]),scale(neg(mul(r,y)),D_expr))
    B=total(f[8],scale(neg(mul(Y,bi)),A),scale(mul(r,v),D_expr),scale(mul(r,Y),e))
    AA=total(f[9],scale(neg(s),f[0]),scale(neg(mul(s,Y)),d_expr))
    BB=total(f[10],scale(neg(mul(y,beta)),AA),scale(mul(s,C),d_expr),scale(mul(s,y),E))
    assert A[0]==sub(c(-1),mul(r,h))
    assert AA[0]==sub(c(-1),mul(s,h))
    assert B[0]==sub(mul(Y,bi),alpha)
    assert BB[0]==sub(mul(y,beta),ai)
    H=total(f[0],scale(neg(mul(beta,y)),B),scale(neg(alpha),BB));assert H[0]==h
    ONE=total(scale(c(-1),A),scale(neg(r),H));assert ONE[0]==c(1)
    # Independent reconstruction from the final expanded multipliers.
    assert add(*(mul(a,b) for a,b in zip(ONE[1],F)))==c(1)
    mod2=lambda p:{e:a%2 for e,a in p.items() if a%2}
    assert mod2(add(*(mul(mod2(a),mod2(b)) for a,b in zip(ONE[1],F))))==c(1)
    corrupt=[dict(a) for a in ONE[1]]
    exponent=next(iter(corrupt[0]));corrupt[0][exponent]+=1
    assert add(*(mul(a,b) for a,b in zip(corrupt,F)))!=c(1)
    out={'variables':R.names,'inverted_variables':['alpha','beta'],
         'equations':[R.dump(f) for f in F],'multipliers':[R.dump(f) for f in ONE[1]],
         'identity_constant':1,'coefficient_ring':'integer Laurent polynomials',
         'scope':'Necessary rank-one involutory root and independent-covector closure coefficients for13a373 over every field; these are KCH eigenvalue minusone only in characteristic not2. No common symmetric or Hermitian form is assumed. Noninvolutory meridian images and dependent seed covectors are outside this unit-certificate scope.'}
    assert out == DATA['unit_identity'], 'retained unit certificate differs'
    return {'integer_Laurent_identity':1,'equations':len(F),
            'multiplier_terms':sum(len(f) for f in ONE[1]),
            'identity_reduction_mod2':1,'corrupt_multiplier_rejected':True}



def paired_siblings():
    R=Ring('x X y Y z Z');x,X,y,Y,z,Z=[R.v(n) for n in R.names]
    c=R.c;add,sub,mul,neg=R.add,R.sub,R.mul,R.neg
    B=[[c(2),x,y],[X,c(2),z],[Y,Z,c(2)]]
    def dot(f,u):return add(*(mul(f[i],B[i][j],u[j]) for i in range(3) for j in range(3)))
    def vs(u,v,a):return [sub(x,mul(a,y)) for x,y in zip(u,v)]
    def cp(p):return [[dict(a) for a in v] for v in p]
    def conjugate(o,a):
        u,f=o;v,g=a
        return vs(v,u,dot(f,v)),vs(g,f,dot(g,u))
    U,V,W=[[c(int(i==j)) for j in range(3)] for i in range(3)]
    def pair(u):return (u,[dict(a) for a in u])
    F,G,J=U,V,W
    T=conjugate(pair(V),pair(W));A=conjugate(pair(U),T)
    t,tau=dot(F,T[0]),dot(T[1],U)
    assert t==sub(y,mul(z,x)) and tau==sub(Y,mul(Z,X))
    plans={375:[(6,9,10),(10,7,8),(4,10,11),(11,4,5),(8,11,12)],
           391:[(9,6,5),(6,9,10),(4,10,11),(5,11,12),(11,7,8)],
           392:[(12,6,5),(5,12,11),(4,11,10),(6,10,9),(10,7,8)]}
    reports=[]
    for name,plan in plans.items():
        roots={0:pair(U),4:pair([neg(a) for a in U]),6:pair(V),7:pair(V),
               12 if name==392 else 9:pair(W)}
        for o,a,b in plan:
            roots[b]=conjugate(roots[o],roots[a])
            assert dot(roots[b][1],roots[b][0])==c(2)
        if name in [375,391]:
            assert roots[10]==T and roots[11]==A
            assert roots[12][0][0]==neg(t) and roots[12][1][0]==neg(tau)
            if name==375:
                prod=sub(mul(t,tau),c(1))
                assert roots[5][0]==vs([mul(prod,a) for a in U],T[0],tau)
                assert roots[5][1]==vs([mul(prod,a) for a in F],T[1],t)
                assert roots[8]==conjugate(T,pair(V))
            else:
                D=conjugate(pair(W),pair(V))
                assert roots[5]==D and roots[12]==conjugate(D,A)
                assert roots[8]==conjugate(A,pair(V))
                assert conjugate(A,roots[8])==pair(V)
        else:
            h=sub(c(1),mul(z,Z))
            tt=add(mul(z,x),mul(h,y));TT=add(mul(Z,X),mul(h,Y))
            k=sub(mul(add(c(2),h),z),mul(tt,X))
            K=sub(mul(add(c(2),h),Z),mul(TT,x))
            er={5:[c(0),c(1),neg(Z)],11:[c(0),z,h],10:[neg(tt),z,h],
                9:[neg(tt),sub(z,k),h],8:[mul(K,tt),sub(c(1),mul(K,z)),neg(mul(K,h))]}
            ef={5:[c(0),c(1),neg(z)],11:[c(0),Z,h],10:[neg(TT),Z,h],
                9:[neg(TT),sub(Z,K),h],8:[mul(k,TT),sub(c(1),mul(k,Z)),neg(mul(k,h))]}
            assert all(roots[i]==(er[i],ef[i]) for i in er)
        reports.append({'name':'13a_'+str(name),'paired_steps':5,'root_and_covector_formulas':True})
    return reports

def local_folds():
    R=Ring('a A b B c C');a,aa,b,bb,c,cc=[R.v(n) for n in R.names]
    inv=R.v('a',-1);num=R.c;add,sub,mul,neg=R.add,R.sub,R.mul,R.neg
    U,V,W=[[num(int(i==j)) for j in range(3)] for i in range(3)]
    def pair(v):return (v,[dict(a) for a in v])
    def setup(matrix):
        def dot(f,u):return add(*(mul(f[i],matrix[i][j],u[j]) for i in range(3) for j in range(3)))
        def conjugate(o,p):
            u,f=o;v,g=p
            return ([sub(x,mul(dot(f,v),y)) for x,y in zip(v,u)],
                    [sub(x,mul(dot(g,u),y)) for x,y in zip(g,f)])
        return dot,conjugate
    raw=[[num(2),a,b],[aa,num(2),c],[bb,cc,num(2)]]
    dot,conjugate=setup(raw)
    for typ in ['I','II']:
        p2=conjugate(pair(U),pair(V)) if typ=='I' else conjugate(pair(W),pair(V))
        p3=conjugate(pair(W),p2) if typ=='I' else conjugate(pair(U),p2)
        closing=conjugate(p3,pair(U))
        t=sub(cc,mul(a,bb)) if typ=='I' else sub(a,mul(cc,b))
        T=sub(c,mul(aa,b)) if typ=='I' else sub(aa,mul(c,bb))
        H=sub(neg(aa),mul(T,bb)) if typ=='I' else neg(T)
        K=sub(neg(a),mul(t,b)) if typ=='I' else neg(t)
        assert dot(p3[1],U)==H and dot(U,p3[0])==K
        if typ=='I':
            assert closing[0]==[add(num(1),mul(a,H)),neg(H),mul(t,H)]
            assert closing[1]==[add(num(1),mul(aa,K)),neg(K),mul(T,K)]
        else:
            assert closing[0]==[add(num(1),mul(t,H)),neg(H),mul(cc,H)]
            assert closing[1]==[add(num(1),mul(T,K)),neg(K),mul(c,K)]
        special=([[num(2),a,b],[inv,num(2),mul(inv,b)],[bb,mul(a,bb),num(2)]] if typ=='I' else
                 [[num(2),a,b],[inv,num(2),num(0)],[bb,num(0),num(2)]])
        dd,conj=setup(special)
        q2=conj(pair(U),pair(V)) if typ=='I' else conj(pair(W),pair(V))
        q3=conj(pair(W),q2) if typ=='I' else conj(pair(U),q2)
        q4=conj(pair(V),q3);qj=conj(q2,pair(W))
        assert q4==([neg(mul(a,u)) for u in U],[neg(mul(inv,u)) for u in U])
        assert qj==pair(W)
        assert conj(q3,pair(U))==([mul(inv,u) for u in V],[mul(a,u) for u in V])
    # Alternate old-four local loop, with arbitrary six directed pairings.
    dot,conj=setup(raw)
    p0=conj(pair(V),pair(U));p4=conj(pair(U),pair(V));p2=conj(p0,pair(V))
    k=sub(num(1),mul(a,aa))
    assert p2==([a,k,num(0)],[aa,k,num(0)])
    p7=conj(p2,pair(W))
    r=dot(p2[1],W);s=dot(W,p2[0])
    assert p7==([neg(mul(r,a)),neg(mul(r,k)),num(1)],
                 [neg(mul(s,aa)),neg(mul(s,k)),num(1)])
    p=dot(p7[1],U);q=dot(U,p7[0]);closing=conj(p7,pair(U))
    assert closing==([add(num(1),mul(p,r,a)),mul(p,r,k),neg(p)],
                     [add(num(1),mul(q,s,aa)),mul(q,s,k),neg(q)])
    special=[[num(2),a,num(0)],[inv,num(2),c],[num(0),cc,num(2)]]
    dd,conj=setup(special)
    p0=conj(pair(V),pair(U));p4=conj(pair(U),pair(V));p2=conj(p0,pair(V));p7=conj(p2,pair(W))
    assert p4==([neg(mul(a,u)) for u in p0[0]],[neg(mul(inv,u)) for u in p0[1]])
    assert p7==pair(W)
    return {'A6_raw_local_types_checked':2,'paired_normalized_folds_checked':3,
            'independent_formal_pairings':6,'only_inverted_parameter':'a'}

def crossing_set(pd):
    parent={e:e for row in pd for e in row}
    def find(a):
        while parent[a]!=a:a=parent[a]
        return a
    for row in pd:parent[find(row[1])]=find(row[3])
    arcs=defaultdict(list)
    for e in parent:arcs[find(e)].append(e)
    arcs=sorted(sorted(v) for v in arcs.values())
    where={e:i for i,a in enumerate(arcs) for e in a}
    return {(where[b],)+tuple(sorted((where[a],where[c]))) for a,b,c,d in pd}

def incidences16():
    inputs=[BASE/'seven_sheet_cover_inputs.json',BASE/'a6_affine_inputs.json']
    cases={r['name']:r for p in inputs for r in json.loads(p.read_text())['cases']}
    key=lambda t:(t[0],)+tuple(sorted(t[1:]))
    local=[(2,6,7),(7,2,3),(3,1,0),(1,3,4),(0,1,2)]
    def reduced(name,local,fold):
        case=cases[name];cr=crossing_set(case['pd'])
        loc={key(t) for t in local};assert loc<=cr and len(loc)==5
        remaining={key(tuple(fold.get(i,i) for i in t)) for t in cr-loc}
        assert len(remaining)==8
        known=set(case['seed_arcs'])
        for o,a,b in case.get('seed_trace',case.get('trace')):
            assert o in known and a in known and b not in known
            assert key((o,a,b)) in cr;known.add(b)
        assert known==set(range(13))
        return remaining
    old={i:reduced('13a_'+str(i),local,{4:0,7:6}) for i in [373,375,391,392]}
    alt=[(3,1,0),(1,3,4),(0,3,2),(2,6,7),(7,1,2)]
    reports=[]
    for new,prior in [(3264,373),(3266,375),(3275,391),(3276,392)]:
        assert reduced('13a_'+str(new),alt,{4:0,7:6})==old[prior]
        reports.append({'name':'13a_'+str(new),'reduced_system':'13a_'+str(prior),'identity_relabeling':True})
    specs=[(304,3232,6,7,[0,5,6,8,9,10,11,12],[9,10,11,12,0,5,6,8],375),
           (362,3259,6,7,[0,5,6,8,9,10,11,12],[10,11,12,0,5,6,8,9],375),
           (623,3374,8,9,[0,5,6,7,8,10,11,12],[8,6,5,0,12,11,10,9],373),
           (629,3375,8,9,[0,5,6,7,8,10,11,12],[5,0,12,11,10,9,8,6],373)]
    for first,second,j,jj,labels,target,prior in specs:
        mapping=dict(zip(labels,target))
        for new,typ in [(first,'I'),(second,'II')]:
            loc=([(0,1,2),(j,2,3),(3,0,1),(1,3,4),(2,j,jj)] if typ=='I' else
                 [(j,1,2),(0,2,3),(3,0,1),(1,3,4),(2,j,jj)])
            reduced_set=reduced('13a_'+str(new),loc,{4:0,jj:j})
            assert {key(tuple(mapping[i] for i in t)) for t in reduced_set}==old[prior]
            reports.append({'name':'13a_'+str(new),'local_type':typ,'reduced_system':'13a_'+str(prior),'relabeling':mapping})
    return {'base_reduced_systems':4,'transfers':reports,
            'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}}

def basis_changes():
    paths=[BASE/'seven_sheet_cover_inputs.json',BASE/'a6_affine_inputs.json']
    cases={r['name']:r for p in paths for r in json.loads(p.read_text())['cases']}
    specs=[(373,[1,2,6,8],[]),(375,[1,2,6,9],[]),
      (391,[1,2,6,9],[(9,5,6)]),(392,[1,2,6,12],[(5,11,12),(12,5,6)]),
      (3264,[1,3,6,10],[(3,1,0),(10,0,12),(12,5,6)]),
      (3266,[1,3,6,9],[(3,1,0),(5,8,9),(9,0,12),(12,5,6)]),
      (3275,[1,3,6,9],[(9,5,6)]),(3276,[1,3,6,9],[(3,1,0),(9,0,12),(12,5,6)]),
      (304,[0,1,6,9],[(9,5,6)]),(3232,[0,1,6,9],[(9,5,6)]),
      (362,[0,1,6,11],[(11,5,6)]),(3259,[0,1,6,11],[(11,5,6)]),
      (623,[0,1,5,8],[]),(3374,[0,1,5,8],[]),
      (629,[0,1,5,8],[(5,7,8)]),(3375,[0,1,5,8],[(5,7,8)])]
    R=Ring(' '.join('a%d%d'%(i,j) for i in range(4) for j in range(4) if i!=j))
    num=R.c;add,sub,mul=R.add,R.sub,R.mul
    B=[[num(2) if i==j else R.v('a%d%d'%(i,j)) for j in range(4)] for i in range(4)]
    def dot(f,u):return add(*(mul(f[i],B[i][j],u[j]) for i in range(4) for j in range(4)))
    def det(matrix):
        terms=[]
        for p in permutations(range(4)):
            parity=sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))%2
            terms.append(mul(num(-1 if parity else 1),*(matrix[i][p[i]] for i in range(4))))
        return add(*terms)
    result=[]
    key=lambda t:(t[0],)+tuple(sorted(t[1:]))
    for index,new,plan in specs:
        name='13a_'+str(index);case=cases[name];old=case['seed_arcs'];cr=crossing_set(case['pd'])
        roots={a:[num(int(i==j)) for j in range(4)] for i,a in enumerate(old)}
        forms={a:[dict(x) for x in v] for a,v in roots.items()}
        for o,a,b in plan:
            assert key((o,a,b)) in cr and o in roots and a in roots and b not in roots
            root_factor=dot(forms[o],roots[a]);form_factor=dot(forms[a],roots[o])
            roots[b]=[sub(x,mul(root_factor,y)) for x,y in zip(roots[a],roots[o])]
            forms[b]=[sub(x,mul(form_factor,y)) for x,y in zip(forms[a],forms[o])]
            assert dot(forms[b],roots[b])==num(2)
        assert det([roots[a] for a in new])==num(1)
        assert det([forms[a] for a in new])==num(1)
        result.append({'name':name,'original_seeds':old,'changed_seeds':new,
                       'literal_pair_derivations':[list(t) for t in plan],
                       'root_change_determinant':1,'covector_change_determinant':1})
    return result


def main():
    assert len(DATA['names']) == len(set(DATA['names'])) == 16
    assert set(DATA['names']) == {r['name'] for r in DATA['basis_changes']}
    primary = {
      'sources/cornwell-kch-1310.7526v3.pdf':'991940273daed005b8af270233cceaa58c11d69168efb7006f3bfa9bd598bcbd',
      'sources/cornwell-character-cord-1509.04962v1.pdf':'1845d56c7be1f824b22aceee6d26a816e84fab05d0baac80c15bb2bff16c067b'}
    assert primary == DATA['primary_source_sha256']
    for name,pin in primary.items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest() == pin
    raw = raw_first_loop()
    trace = paired_trace()
    unit = certificate()
    siblings = paired_siblings()
    folds = local_folds()
    changes = basis_changes()
    diagrams = json.loads(json.dumps(incidences16()))
    assert changes == DATA['basis_changes']
    assert diagrams == DATA['diagram_incidence_checks']
    assert folds == DATA['raw_local_folds']
    files = ['INVOLUTORY_KCH_OBSTRUCTION.tex',Path(__file__).name,INPUT.name,
             'seven_sheet_cover_inputs.json','a6_affine_inputs.json'] + list(primary)
    report = {'status':'PASS','named_cases':16,'names':DATA['names'],
        'raw_first_loop':raw,'paired_373_trace':trace,'integer_unit_identity':unit,
        'other_paired_systems':siblings,'raw_local_folds':folds,
        'basis_changes':changes,'diagram_incidence_checks':diagrams,
        'representation_scope':DATA['representation_scope'],
        'augmentation_scope':DATA['augmentation_scope'],
        'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files},
        'assurance':'Internal exact integer/Laurent algebra and literal PD incidence checks. Four-seed group generation and named identities are bound to the retained diagram inputs. The augmentation-irrep correspondence is an explicit primary import; neither formal verification nor external specialist review is claimed.'}
    (HERE/'involutory_kch_obstruction_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'PASS','named_cases':16,'integer_unit_identity':1,
        'root_and_covector_basis_changes':16,'literal_incidence_transfers':12,
        'representation_fields':'all','KCH_eigenvalue_in_characteristic_not2':-1,
        'irreducible_dimension_upper':3},indent=2))

if __name__ == '__main__':
    main()
