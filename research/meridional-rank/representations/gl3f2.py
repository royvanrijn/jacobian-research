#!/usr/bin/env python3
"""Small exact binary matrix group, modules and tuple orbits."""
from itertools import product


def matmul(a,b):
    return tuple(xor(b[k] for k in range(len(a)) if row>>k&1) for row in a)
def xor(values):
    out=0
    for v in values:out^=v
    return out
def matvec(a,v):return sum(((r&v).bit_count()%2)<<i for i,r in enumerate(a))
def transpose(a):return tuple(sum(((r>>i)&1)<<j for j,r in enumerate(a)) for i in range(len(a)))


def rref(rows,columns):
    a=list(rows);pivots=[];k=0
    for j in range(columns):
        i=next((i for i in range(k,len(a)) if a[i]>>j&1),None)
        if i is None:continue
        a[k],a[i]=a[i],a[k]
        for i in range(len(a)):
            if i!=k and a[i]>>j&1:a[i]^=a[k]
        pivots.append(j);k+=1
        if k==len(a):break
    return a[:k],pivots
def kernel(rows,columns):
    a,pivots=rref(rows,columns);out=[]
    for j in range(columns):
        if j in pivots:continue
        out.append((1<<j)^xor(1<<p for p,row in zip(pivots,a) if row>>j&1))
    return out
def rank(rows,columns):return len(rref(rows,columns)[1])


def generated(values,table):
    seen={0};queue=[0]
    for a in queue:
        for b in values:
            c=table[a][b]
            if c not in seen:seen.add(c);queue.append(c)
    return frozenset(seen)


def group_model():
    unit=(1,2,4)
    matrices=sorted(tuple((x>>(3*i))&7 for i in range(3)) for x in range(512)
                    if rank([(x>>(3*i))&7 for i in range(3)],3)==3)
    matrices.remove(unit);matrices.insert(0,unit);assert len(matrices)==168
    index={a:i for i,a in enumerate(matrices)}
    table=[[index[matmul(a,b)] for b in matrices] for a in matrices]
    inverses=[next(j for j in range(168) if table[i][j]==0) for i in range(168)]
    orders=[]
    for i in range(168):
        x=i;n=1
        while x:x=table[x][i];n+=1
        orders.append(n)
    dual=[index[transpose(matrices[inverses[i]])] for i in range(168)]
    inner=[tuple(table[table[g][i]][inverses[g]] for i in range(168)) for g in range(168)]
    autos=sorted(set(inner+[tuple(a[dual[i]] for i in range(168)) for a in inner]))
    assert len(autos)==336 and len(set(dual))==168
    assert all(dual[table[i][j]]==table[dual[i]][dual[j]] for i in range(168) for j in range(168))
    pairs=[(a,b) for a in range(168) for b in range(168) if orders[a]==2 and orders[b]==3 and orders[table[a][b]]==7]
    assert len(pairs)==336 and len(generated(pairs[0],table))==168
    assert len({(a[pairs[0][0]],a[pairs[0][1]]) for a in autos})==336
    core=[[table[table[b][inverses[a]]][b] for a in range(168)] for b in range(168)]
    return matrices,table,inverses,orders,autos,dual,inner,core,pairs[0]


def all_subgroups(table):
    unit=frozenset([0]);gens={unit:()};queue=[unit]
    for h in queue:
        for g in range(len(table)):
            if g in h:continue
            values=gens[h]+(g,);larger=generated(values,table)
            if larger not in gens:gens[larger]=values;queue.append(larger)
    return sorted(gens,key=lambda h:(len(h),tuple(sorted(h))))


def tuple_orbits(arity,table,orders,autos):
    result=[]
    for first_position in range(arity-1):
        for order in [2,3,4,7]:
            first=next(i for i in range(1,168) if orders[i]==order)
            assert len({a[first] for a in autos})==sum(x==order for x in orders)
            stabilizer=[a for a in autos if a[first]==first]
            for tail in product(range(168),repeat=arity-first_position-1):
                values=(0,)*first_position+(first,)+tail
                if any(tuple(a[i] for i in values)<values for a in stabilizer):continue
                if len(generated(values,table))==168:result.append(values)
    return sorted(result)


def steinberg(matrices,inverses):
    # H1 of the connected bipartite Fano incidence graph has dimension8.
    flags=[(p,l) for p in range(1,8) for l in range(1,8) if (p&l).bit_count()%2==0]
    assert len(flags)==21
    ends=[(p-1,7+l-1) for p,l in flags];parent=list(range(14));tree=[]
    def root(i):
        while parent[i]!=i:i=parent[i]
        return i
    for e,(a,b) in enumerate(ends):
        if root(a)!=root(b):parent[root(a)]=root(b);tree.append(e)
    assert len(tree)==13
    chords=[e for e in range(21) if e not in tree];assert len(chords)==8
    adjacent=[[] for _ in range(14)]
    for e in tree:
        a,b=ends[e];adjacent[a].append((b,e));adjacent[b].append((a,e))
    cycles=[]
    for e in chords:
        a,b=ends[e];paths={a:0};queue=[a]
        for v in queue:
            for w,f in adjacent[v]:
                if w not in paths:paths[w]=paths[v]^(1<<f);queue.append(w)
        cycles.append(paths[b]^(1<<e))
    output=[]
    for g,a in enumerate(matrices):
        dual=transpose(matrices[inverses[g]])
        perm=[flags.index((matvec(a,p),matvec(dual,l))) for p,l in flags]
        image_cycles=[xor(1<<perm[e] for e in range(21) if c>>e&1) for c in cycles]
        m=tuple(sum(((c>>e)&1)<<j for j,c in enumerate(image_cycles)) for e in chords)
        for c,v in zip(image_cycles,[sum(((c>>e)&1)<<i for i,e in enumerate(chords)) for c in image_cycles]):
            assert xor(cycles[j] for j in range(8) if v>>j&1)==c
        output.append(m)
    assert output[0]==tuple(1<<i for i in range(8))
    return output,{'flags':flags,'tree':tree,'chords':chords,'cycles':cycles}


def cocycle_matrix(images,inverses,words):
    d=len(images[0]);unit=tuple(1<<i for i in range(d));out=[]
    for word in words:
        prefix=unit;block=[0]*d
        for letter in word:
            i=abs(letter)-1
            if letter<0:prefix=matmul(prefix,inverses[i])
            for j,row in enumerate(prefix):block[j]^=row<<(d*i)
            if letter>0:prefix=matmul(prefix,images[i])
        assert prefix==unit;out.extend(block)
    return out


def cocycles(images,inverses,words):
    d=len(images[0]);rows=cocycle_matrix(images,inverses,words)
    principal=[row^(1<<j) for a in images for j,row in enumerate(a)]
    assert all(xor(principal[j] for j in range(len(principal)) if row>>j&1)==0 for row in rows)
    basis=kernel(rows,d*len(images));inv=d-rank(principal,d)
    return {'dimension':d,'rows':rows,'columns':d*len(images),'kernel_basis':basis,
            'Z1_dimension':len(basis),'invariants':inv,'H1_dimension':len(basis)-d+inv,
            'sharp_cover_rank3':len(basis)>2*d}


def root_cocycles(images,words):
    roots=[]
    for a in images:
        diff=tuple(row^(1<<j) for j,row in enumerate(a));assert rank(diff,3)==1
        roots.append(next(sum(((r>>k)&1)<<j for j,r in enumerate(diff)) for k in range(3) if any(r>>k&1 for r in diff)))
    rows=[]
    for word in words:
        prefix=(1,2,4);block=[0]*3
        for letter in word:
            i=abs(letter)-1;v=matvec(prefix,roots[i])
            for j in range(3):block[j]^=((v>>j)&1)<<i
            prefix=matmul(prefix,images[i])
        assert prefix==(1,2,4);rows.extend(block)
    basis=kernel(rows,len(images))
    return {'roots':roots,'rows':rows,'columns':len(images),'kernel_basis':basis,
            'dimension':len(basis),'sharp_four':len(basis)>=4}
