#!/usr/bin/env sage-python
"""The prescribed determinant4100 common core has incompatible two-primary glue.

Producer: exact frame coordinates and PARI short vectors.
Replay: independent exact LDL sphere enumeration and elementary color
refinement prove the core automorphism group is only +/-1.
"""
import argparse,hashlib,json,zipfile
from pathlib import Path
from sage.all import QQ,ZZ,matrix,vector,gcd,pari
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
MARK=ART/'x1092_class1_realization_marking_v1.json'
PACKETS=ART/'det1092_pruned_anchor_packets_v1.zip'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
OUTPUT=ART/'curve302_class1_prescribed_core_glue_obstruction_v1.json'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rows(M):return [[int(x)for x in row]for row in M]
def canonical(v):
    t=tuple(map(int,v))
    return t if next(x for x in t if x)>0 else tuple(-x for x in t)
def exact_short_vectors(G):
    """Complete rational-LDL recursion; no numerical square-root bounds."""
    n=G.nrows();L=matrix.identity(QQ,n);diagonal=[]
    for i in range(n):
        d=G[i,i]-sum(L[i,k]**2*diagonal[k]for k in range(i));assert d>0
        diagonal.append(d)
        for j in range(i+1,n):
            L[j,i]=(G[j,i]-sum(L[j,k]*L[i,k]*diagonal[k]for k in range(i)))/d
    assert L*matrix.diagonal(QQ,diagonal)*L.transpose()==G
    word=[ZZ(0)]*n;out=[]
    def visit(i,remaining):
        if i<0:
            if any(word):out.append(tuple(word))
            return
        shift=sum((L[j,i]*word[j]for j in range(i+1,n)),QQ(0))
        bound=remaining/diagonal[i]
        radius=ZZ(bound.floor()).isqrt()+1
        for value in range(int((-shift).floor()-radius-1),int((-shift).ceil()+radius+2)):
            rest=remaining-diagonal[i]*(value+shift)**2
            if rest>=0:word[i]=ZZ(value);visit(i-1,rest)
        word[i]=0
    visit(n-1,QQ(4));return out
def refine(V,G):
    H=V*G*V.transpose();N=V.nrows()
    absolute=[[abs(int(H[i,j]))for j in range(N)]for i in range(N)]
    colors=[0]*N;counts=[1]
    while len(set(colors))<N:
        signatures=[(colors[i],tuple(sorted((colors[j],absolute[i][j])for j in range(N)if j!=i)))for i in range(N)]
        distinct=sorted(set(signatures));lookup={x:i for i,x in enumerate(distinct)}
        new=[lookup[s]for s in signatures]
        assert len(set(new))>len(set(colors)),'refinement did not prove line rigidity'
        colors=new;counts.append(len(set(colors)))
    seen={0};todo=[0]
    while todo:
        i=todo.pop()
        for j in range(N):
            if j not in seen and absolute[i][j]:seen.add(j);todo.append(j)
    assert len(seen)==N and V.rank()==G.nrows()
    return counts
def compute(check=False):
    mark=read(MARK);core=mark['requested_shared_core_witness']
    C=matrix(ZZ,core['common_core_basis_in_niemeier']);H=matrix(ZZ,core['common_core_gram'])
    assert H.det()==core['common_core_determinant']==4100
    with zipfile.ZipFile(PACKETS)as z:packet=json.loads(z.read('anchor-16.json'))
    embeddings=[]
    for index in (25,33):
        e=packet['embeddings'][index];G=matrix(ZZ,e['gram']);B=matrix(ZZ,e['complement_basis_in_ambient'])
        K=matrix(ZZ,B.transpose().solve_right(C.transpose()).transpose())
        assert K*G*K.transpose()==H and K.row_module().saturation()==K.row_module()
        normal=(K*G).right_kernel().basis_matrix();assert normal.nrows()==1
        v=normal.row(0);T=K.stack(normal);index_value=abs(T.det());assert index_value==2050
        norm=v*G*v;div=gcd(G*v);assert norm==1119300 and div==546 and G.det()==1092
        coordinates=T.inverse()
        selected=next(i for i,row in enumerate(coordinates)if row[-1].denominator()==index_value)
        r=coordinates.row(selected)
        generator=[int((x*index_value)%index_value)for x in r[:-1]]
        coefficient=int((r[-1]*index_value)%index_value)
        assert gcd([index_value,*generator])==1 and gcd(coefficient,index_value)==1
        # Since K and the normal are primitive, the projection has order2050.
        # Verify every frame basis vector projects into this displayed cyclic group.
        inverse=pow(coefficient,-1,int(index_value))
        for row in coordinates:
            multiple=int(row[-1]*index_value)*inverse%index_value
            assert all((index_value*x-multiple*a)%index_value==0 for x,a in zip(row[:-1],generator))
        two=[x%2 for x in generator]
        embeddings.append({'embedding':index,'frame_gram':rows(G),'core_in_frame':rows(K),
                           'normal_in_frame':list(map(int,v)),'normal_norm':int(norm),'normal_divisibility':int(div),
                           'index':int(index_value),'generator_frame_row':selected,
                           'core_projection_generator_mod_2050':generator,
                           'normal_coefficient_mod_2050':coefficient,'unique_order_two_glue_bits':two})
    assert embeddings[0]['unique_order_two_glue_bits']!=embeddings[1]['unique_order_two_glue_bits']
    # Complete minimal-vector lines determine all core automorphisms.
    if check:
        stored=read(OUTPUT);U=matrix(ZZ,stored['lll_columns']);assert abs(U.det())==1
        reduced=U.transpose()*H*U
        short=exact_short_vectors(reduced)
        minimum=set()
        for row in short:
            v=U*vector(ZZ,row);assert v*H*v==4
            minimum.add(canonical(v))
        assert sorted(minimum)==[tuple(v)for v in stored['minimal_vector_lines']]
    else:
        U=matrix(ZZ,pari.qflllgram(H));assert abs(U.det())==1
        raw=pari.qfminim(H,4);V=matrix(ZZ,raw[2]).transpose()
        minimum={canonical(v)for v in V}
        assert int(raw[0])==2*len(minimum)
    V=matrix(ZZ,sorted(minimum));assert V.nrows()==516
    counts=refine(V,H)
    # The actual low-degree realization has a different literal intersection.
    G=matrix(ZZ,read(PARENT)['generic_height_gram']);w=vector(ZZ,mark['fibre_D'][2:])
    actual=(matrix(ZZ,1,17,list(G*w))).right_kernel().basis_matrix()
    actual_gram=actual*G*actual.transpose();assert actual_gram.det()==13104
    result={'schema':'curve302.class1-prescribed-core-glue-obstruction.v1',
            'bindings':{str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),MARK,PACKETS,PARENT]},
            'core_gram':rows(H),'embeddings':embeddings,'lll_columns':rows(U),
            'minimal_vector_lines':rows(V),'refinement_class_counts':counts,
            'minimal_vectors_span_rank':16,'nonorthogonality_graph_connected':True,
            'core_automorphisms':'exactly +I and -I',
            'existing_realization_common_core_gram':rows(actual_gram),
            'existing_realization_common_core_determinant':13104,
            'status':'PASS_PRESCRIBED_COMMON_CORE_HAS_INCOMPATIBLE_2_GLUE',
            'theorem':'No integral isometry U+(-frame25) -> U+(-frame33) identifies the two retained determinant4100 core embeddings, even allowing an arbitrary core isometry. Their cyclic order2050 projection subgroups have different order2 elements, while every core isometry is +/-I and fixes each such subgroup.',
            'boundary':'This excludes simultaneous marked geometric realization of this particular retained common-core witness. Class1 and class6 still have rational realizations; other overlaps are not excluded.',
            'limits':{'producer_seconds':60,'independent_replay_seconds':120,'core_rank':16,'sphere_norm_bound':4,'new_fibrations':0}}
    print('PASS core automorphisms +/-1; refinement',counts,'; incompatible order2 glue',flush=True)
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args()
    result=compute(args.check)
    if args.check:assert read(OUTPUT)==result
    else:
        with OUTPUT.open('x')as out:json.dump(result,out,indent=2,sort_keys=True);out.write('\n')
