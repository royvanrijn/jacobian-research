#!/usr/bin/env python3
"""Standalone exact positive rank3 certificates for8 named A6 knots.

Standard library only; no producer/canonical checker imports. Polynomial
convolution followed by monic-Phi reduction is independent of the producer's
cyclic-convolution implementation. --census re-reads the named census.
"""
import argparse
from copy import deepcopy
import hashlib,json,math
from itertools import permutations
from pathlib import Path

HERE=Path(__file__).resolve().parent
INPUT=HERE/'a6_rank_three_augmentation_inputs.json'
SUPPORT=HERE


class Field:
    def __init__(self,record):
        self.p=p=record['prime']
        assert p>2 and all(p%q for q in range(2,math.isqrt(p)+1))
        assert record['minimal_polynomial_coefficients']==[1]*p
        # Explicit Eisenstein check for Phi_p(X+1): its coefficients are
        # binomial(p,j+1), j=0,...,p-1. Constantp is not divisible byp^2,
        # all other nonleading coefficients are divisible byp.
        shifted=[math.comb(p,j+1) for j in range(p)]
        assert shifted[-1]==1 and shifted[0]==p and shifted[0]%(p*p)!=0
        assert all(a%p==0 for a in shifted[:-1])
        self.n=p-1;self.zero=(0,)*self.n;self.one=(1,)+(0,)*(self.n-1)
    def c(self,n):return (n,)+(0,)*(self.n-1)
    def add(self,*ps):return tuple(sum(a[i] for a in ps) for i in range(self.n))
    def neg(self,a):return tuple(-x for x in a)
    def sub(self,a,b):return self.add(a,self.neg(b))
    def mul(self,a,b):
        coefficients=[0]*(2*self.n-1)
        for i,x in enumerate(a):
            if x:
                for j,y in enumerate(b):
                    if y:coefficients[i+j]+=x*y
        for k in range(len(coefficients)-1,self.n-1,-1):
            value=coefficients[k]
            if value:
                coefficients[k]=0
                for j in range(self.n):coefficients[k-self.n+j]-=value
        return tuple(coefficients[:self.n])
    def vector(self,v):
        assert len(v)==self.n and all(isinstance(a,int) for a in v)
        return tuple(v)
    def determinant(self,A):
        value=self.zero;n=len(A)
        for p in permutations(range(n)):
            term=self.c((-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n)))
            for i in range(n):term=self.mul(term,A[i][p[i]])
            value=self.add(value,term)
        return value


def mm(F,A,B):
    return [[F.add(*(F.mul(A[i][k],B[k][j]) for k in range(len(B)))) for j in range(len(B[0]))]
            for i in range(len(A))]
def dot(F,B,u,v):
    return F.add(*(F.mul(F.mul(u[i],B[i][j]),v[j]) for i in range(3) for j in range(3)))


def recipe(F,record,roots,B):
    """Bind the manuscript's complete recurrence/sign recipe to the data."""
    name=record['name'];p=F.p
    groups=[(['13a_304','13a_3232'],6,[0,5,6,8,9,10,11,12],
             [17,22,15,34,0,30,1,6],[1,1,1,1,1,-1,-1,1]),
            (['13a_362','13a_3259'],6,[0,5,6,8,9,10,11,12],
             [22,15,34,0,30,1,6,17],[1,1,1,1,-1,-1,1,1]),
            (['13a_623','13a_3374'],8,[0,5,6,7,8,10,11,12],
             [1,1,2,0,5,1,6,3],[1,-1,-1,1,1,1,1,1]),
            (['13a_629','13a_3375'],8,[0,5,6,7,8,10,11,12],
             [2,0,5,1,6,3,1,1],[-1,1,1,1,1,1,1,-1])]
    group_index=next(i for i,row in enumerate(groups) if name in row[0])
    names,j,labels,colours,oldsigns=groups[group_index]
    finalsigns=[[1,-1,1,-1,1,-1,-1,1,1,1,-1,-1,1],
                [1,-1,1,-1,1,-1,-1,1,1,-1,-1,1,1],
                [1,-1,1,-1,1,1,1,-1,-1,1,1,1,1],
                [1,-1,1,-1,1,1,1,1,1,-1,-1,-1,1]][group_index]
    assert record['arc_signs']==finalsigns
    z=(0,1)+(0,)*(p-3);zi=F.one
    for _ in range(p-1):zi=F.mul(zi,z)
    q=F.add(z,zi)
    r=[[F.one,F.zero],[F.zero,F.one]]
    for _ in range(p-1):r.append([F.sub(F.mul(q,a),b) for a,b in zip(r[-1],r[-2])])
    assert r[p]==r[0]
    plane={i:[F.mul(F.c(sign),a) for a in r[colour]]
           for i,colour,sign in zip(labels,colours,oldsigns)}
    A,C=plane[0],plane[j]
    det=lambda u,v:F.sub(F.mul(u[0],v[1]),F.mul(u[1],v[0]))
    delta=det(A,C);assert delta!=F.zero
    unsigned=[[F.mul(F.c(sign),a) for a in u] for u,sign in zip(roots,finalsigns)]
    for i,v in plane.items():
        assert F.mul(delta,unsigned[i][0])==det(v,C)
        assert F.mul(delta,unsigned[i][1])==det(A,v)
        assert unsigned[i][2]==F.zero
    eA=[F.one,F.zero,F.zero];eC=[F.zero,F.one,F.zero];eE=[F.zero,F.zero,F.one]
    EmA=[F.c(-1),F.zero,F.one]
    assert unsigned[1]==eE and unsigned[3]==EmA and unsigned[4]==[F.c(-1),F.zero,F.zero]
    assert unsigned[2]==(EmA if record['local_type']=='I' else eE)
    assert unsigned[j+1]==eC
    b=F.add(F.mul(F.c(2),F.mul(A[0],C[0])),F.mul(q,F.mul(A[0],C[1])),
            F.mul(q,F.mul(A[1],C[0])),F.mul(F.c(2),F.mul(A[1],C[1])))
    c=b if record['local_type']=='I' else F.zero
    assert B==[[F.c(2),b,F.one],[b,F.c(2),c],[F.one,c,F.c(2)]]


def diagram(pd):
    ends={};adj={e:set() for row in pd for e in row}
    for i,row in enumerate(pd):
        assert len(row)==4 and len(set(row))==4
        for j,e in enumerate(row):ends.setdefault(e,[]).append(4*i+j)
        adj[row[1]].add(row[3]);adj[row[3]].add(row[1])
    assert all(len(x)==2 for x in ends.values()) and len(ends)==2*len(pd)
    alpha={a:b for pair in ends.values() for a,b in (pair,pair[::-1])}
    unseen=set(alpha);faces=0
    while unseen:
        start=x=min(unseen)
        while True:
            unseen.remove(x);y=alpha[x];x=4*(y//4)+(y+1)%4
            if x==start:break
        faces+=1
    assert faces==len(pd)+2
    x=0;visited=set()
    while x not in visited:
        visited.add(x);i,j=divmod(x,4);x=alpha[4*i+(j+2)%4]
    assert x==0 and len(visited)==2*len(pd)
    components=[];unused=set(adj)
    while unused:
        queue=[min(unused)];cc=set()
        for x in queue:
            if x in cc:continue
            cc.add(x);queue.extend(adj[x]-cc)
        components.append(sorted(cc));unused-=cc
    components.sort();assert len(components)==len(pd)
    arc={e:i for i,cc in enumerate(components) for e in cc}
    return [(arc[b],arc[a],arc[c]) for a,b,c,d in pd],alpha


def named_map(record):
    _,a=diagram(record['pd']);_,b=diagram(record['census_pd'])
    cert=record['census_map_witness'];f=cert['dart_images']
    orientation=cert['orientation'];swap=cert['over_under_swap']
    assert sorted(f)==list(range(len(a))) and orientation in (-1,1) and swap in (0,1)
    for x,y in enumerate(f):
        assert f[a[x]]==b[y]
        assert (x%2)^(y%2)==swap
        assert f[4*(x//4)+(x+1)%4]==4*(y//4)+(y+orientation)%4


def printed_crossings(record):
    """Reconstruct the manuscript's local5/remaining8 from literal PD."""
    rows=[(['13a_304','13a_3232'],6,
           [(4,8,9),(11,4,5),(5,11,12),(9,5,6),
            (7,9,10),(12,7,8),(8,10,11),(10,12,0)]),
          (['13a_362','13a_3259'],6,
           [(4,10,11),(8,4,5),(5,8,9),(11,5,6),
            (7,9,10),(12,7,8),(9,11,12),(10,12,0)]),
          (['13a_623','13a_3374'],8,
           [(4,6,7),(10,4,5),(5,11,12),(9,5,6),
            (6,12,0),(7,10,11),(11,7,8),(12,9,10)]),
          (['13a_629','13a_3375'],8,
           [(4,10,11),(11,4,5),(5,7,8),(9,5,6),
            (6,12,0),(10,6,7),(7,11,12),(12,9,10)])]
    names,j,remaining=next(row for row in rows if record['name'] in row[0])
    local=([(0,1,2),(j,2,3),(3,0,1),(1,3,4),(2,j,j+1)]
           if record['local_type']=='I' else
           [(j,1,2),(0,2,3),(3,0,1),(1,3,4),(2,j,j+1)])
    key=lambda row:(row[0],)+tuple(sorted(row[1:]))
    actual=diagram(record['pd'])[0]
    actual_keys=set(map(key,actual));local_keys=set(map(key,local))
    assert len(local_keys)==5 and local_keys<=actual_keys
    extracted=[row for row in actual if key(row) not in local_keys]
    assert len(extracted)==8 and set(map(key,extracted))==set(map(key,remaining))
    assert not local_keys&set(map(key,remaining))
    assert local_keys|set(map(key,remaining))==actual_keys


def verify(record):
    named_map(record);printed_crossings(record);F=Field(record['field'])
    roots=[[F.vector(a) for a in u] for u in record['normal_vectors']]
    B=[[F.vector(a) for a in row] for row in record['gram']]
    saved=[[F.vector(a) for a in row] for row in record['cord_matrix']]
    assert len(roots)==13 and all(len(u)==3 for u in roots)
    assert len(B)==3 and all(len(row)==3 for row in B)
    assert all(B[i][j]==B[j][i] for i in range(3) for j in range(3))
    recipe(F,record,roots,B)
    C=[[dot(F,B,u,v) for v in roots] for u in roots]
    assert C==saved and all(C[i][i]==F.c(2) for i in range(13))
    indices=record['nonzero_principal_minor']['indices']
    minor=F.determinant([[C[i][j] for j in indices] for i in indices])
    assert minor==F.vector(record['nonzero_principal_minor']['value']) and minor!=F.zero
    assert F.determinant(B)==minor
    root_det=F.determinant([roots[i] for i in indices]);assert root_det!=F.zero
    covectors=[[F.add(*(F.mul(u[k],B[k][j]) for k in range(3))) for j in range(3)] for u in roots]
    form_det=F.determinant([covectors[i] for i in indices]);assert form_det!=F.zero
    I=[[F.c(int(i==j)) for j in range(3)] for i in range(3)]
    matrices=[[[F.sub(I[i][j],F.mul(u[i],f[j])) for j in range(3)] for i in range(3)]
              for u,f in zip(roots,covectors)]
    assert all(mm(F,M,M)==I for M in matrices)
    crossings,_=diagram(record['pd'])
    for o,a,b in crossings:
        assert mm(F,mm(F,matrices[o],matrices[a]),matrices[o])==matrices[b]
        for l in range(13):assert F.add(C[l][a],C[l][b])==F.mul(C[l][o],C[o][a])
    # Each displacement is a nonzero rank1 outer product, sinceu_i andf_i
    # are nonzero and f_i(u_i)=2. Zero common covector kernel plus connected
    # nonzero pairing graph proves irreducibility by root propagation.
    assert all(any(a!=F.zero for a in u) for u in roots)
    assert all(any(a!=F.zero for a in f) for f in covectors)
    graph={i:{j for j in range(13) if j!=i and C[i][j]!=F.zero} for i in range(13)}
    assert graph[0]==set(range(1,13))
    queue=[0];seen={0}
    for i in queue:
        for j in graph[i]-seen:seen.add(j);queue.append(j)
    assert len(seen)==13
    return {'name':record['name'],'cyclotomic_prime':F.p,'field_irreducibility':'Eisenstein after X+1',
            'cord_rank':3,'root_span_dimension':3,'common_fixed_dimension':0,
            'nonzero_principal_minor_indices':indices,'minor':list(minor),
            'root_pairing_graph_connected':True,'irreducible_dimension':3,
            'pairing_graph_contains_star_at_zero':True,
            'printed_recurrence_and_sign_recipe_verified':True,
            'printed_local5_remaining8_partition_verified':True,
            'literal_reflection_matrix_crossings':13,'cord_scalar_equations':169,
            'named_planar_map_checked':True}


def census(data):
    import importlib.metadata as md
    import snappy_15_knots
    from spherogram import Link
    provenance=data['census_provenance'];root=Path(snappy_15_knots.__file__).parent
    assert md.version('snappy_15_knots')==provenance['version']
    assert md.version('spherogram')==provenance['spherogram_version']
    for pin in provenance['source_files']:
        raw=(root/pin['path']).read_bytes()
        assert len(raw)==pin['bytes'] and hashlib.sha256(raw).hexdigest()==pin['sha256']
    for record in data['cases']:
        link=Link('K'+record['name'].replace('_',''))
        assert [list(row) for row in link.PD_code()]==record['census_pd']
        assert json.loads(json.dumps(link.DT_code()))==record['census_DT']
    return len(data['cases'])


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--census',action='store_true');args=parser.parse_args()
    data=json.loads(INPUT.read_text());cases=data['cases']
    assert len(cases)==len(set(data['names']))==8 and data['names']==[c['name'] for c in cases]
    for file,sha in data['source_sha256'].items():
        assert hashlib.sha256((SUPPORT/file).read_bytes()).hexdigest()==sha
    original=json.loads((SUPPORT/'a6_affine_inputs.json').read_text())
    for record in cases:
        source=next(c for c in original['cases'] if c['name']==record['name'])
        assert all(record[k]==source[k] for k in ['pd','census_pd','census_DT','census_map_witness'])
    results=[verify(c) for c in cases]
    bad=deepcopy(cases[0]);bad['cord_matrix'][0][0][0]+=1
    try:verify(bad)
    except AssertionError:pass
    else:raise AssertionError('corrupted cord matrix accepted')
    bad=deepcopy(cases[0]);bad['field']['minimal_polynomial_coefficients'][0]+=1
    try:verify(bad)
    except AssertionError:pass
    else:raise AssertionError('wrong field polynomial accepted')
    receipt={'status':'PASS','results':results,'cases':8,
             'matrix_crossings_checked':104,'cord_scalar_equations_checked':1352,
             'optional_census_rereads':census(data) if args.census else None,
             'corruptions_rejected':['cord_matrix','field_polynomial'],
             'scope':data['scope'],
             'input_sha256':{INPUT.name:hashlib.sha256(INPUT.read_bytes()).hexdigest(),
                             Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                             'A6_RANK_THREE_AUGMENTATIONS.tex':hashlib.sha256((HERE/'A6_RANK_THREE_AUGMENTATIONS.tex').read_bytes()).hexdigest()}|data['source_sha256']}
    receipt_path=HERE/'a6_rank_three_augmentation_verification.json'
    if receipt_path.is_file():
        previous=json.loads(receipt_path.read_text())
        if 'internal_pdf_reviews' in previous:
            # These reviews retain their own source/PDF hashes. A later
            # source edit does not promote old-edition review to new assurance.
            receipt['internal_pdf_reviews']=previous['internal_pdf_reviews']
    receipt_path.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['status','cases','matrix_crossings_checked','cord_scalar_equations_checked','optional_census_rereads']},indent=2))


if __name__=='__main__':main()
