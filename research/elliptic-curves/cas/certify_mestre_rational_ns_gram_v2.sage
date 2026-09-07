#!/usr/bin/env sage-python
"""Reconstruct six saturated rational divisor lattices; no fibration search."""
import argparse, json, hashlib, importlib.util
from importlib.machinery import SourceFileLoader
from pathlib import Path
from sage.all import QQ, ZZ, GF, PolynomialRing, EllipticCurve, matrix, vector, block_diagonal_matrix
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
INPUT=ART/'mestre_468_replay_bundle_v1.json'
OUT=ART/'mestre_rational_ns_gram_v2.json'
LOCAL=ROOT/'artifacts/local/elliptic-curves/mestre-rational-ns-gram-v2'
helper=ROOT/'elliptic-curves/cas/mestre_generic_height_audit.sage'
loader=SourceFileLoader('height_source',str(helper));spec=importlib.util.spec_from_loader(loader.name,loader);height_source=importlib.util.module_from_spec(spec);loader.exec_module(height_source)

def compute():
    data=json.loads(INPUT.read_text());rows=[]
    for row in data['rows']:
        h=row['generic_heights'];R=PolynomialRing(QQ,'T');K=R.fraction_field();A=R(h['curve_A_coefficients']);B=R(h['curve_B_coefficients']);E=EllipticCurve(K,[A,B])
        raw=[E([K(x),K(y)]) for x,y in h['covariant_points']];cloud=[raw[0]]
        for P in raw[1:]:
            halves=(P-raw[0]).division_points(2)
            assert len(halves)==1 and 2*halves[0]==P-raw[0];cloud.append(halves[0])
        points=[cloud[i] for i in h['seed_indices']];geometry=height_source.setup_height(A,B)
        profiles=[]
        H=matrix(QQ,h['seed_height_gram'])
        for i,P in enumerate(points):
            value,profile=height_source.height(P,A,B,geometry)
            assert value==H[i,i];profiles.append(profile)
        G=matrix(QQ,18,18);G[0,1]=G[1,0]=1;G[1,1]=-2
        for j in range(2,7):G[j,j]=-2
        G[4,5]=G[5,4]=G[5,6]=G[6,5]=1
        for i,p in enumerate(profiles):
            j=i+7;G[j,j]=-2;G[0,j]=G[j,0]=1
            G[1,j]=G[j,1]=p['zero_section_intersection']
            c=p['components']
            for k in (0,1):G[2+k,j]=G[j,2+k]=c[k]
            if c[2]:G[3+c[2],j]=G[j,3+c[2]]=1
            for k in range(i):
                q=profiles[k];d=q['components']
                correction=QQ(c[0]*d[0]+c[1]*d[1])/2+min(c[2],d[2])-QQ(c[2]*d[2])/4
                value=2+p['zero_section_intersection']+q['zero_section_intersection']-H[i,k]-correction
                assert value>=0 and value.denominator()==1
                G[j,7+k]=G[7+k,j]=value
        assert all(c.denominator()==1 for c in G.list());G=matrix(ZZ,G)
        assert G.det()==-468 and all(G[i,i]%2==0 for i in range(18))
        # Split off the displayed U explicitly; the remaining frame is negative definite.
        Uchange=matrix.identity(ZZ,18);Uchange[1,0]=1
        for i in range(7,18):
            Uchange[i,1]=-1;Uchange[i,0]=-profiles[i-7]['zero_section_intersection']-2
        split=Uchange*G*Uchange.transpose()
        assert abs(Uchange.det())==1 and split[:2,:2]==matrix(ZZ,[[0,1],[1,0]])
        assert not any(split[i,j] for i in range(2) for j in range(2,18))
        assert (-split[2:,2:]).is_positive_definite()
        elementary=[int(c) for c in G.elementary_divisors() if abs(c)!=1]
        kernel=matrix(GF(2),G).right_kernel();glues=[]
        for vv in kernel:
            v=vector(ZZ,[int(c) for c in vv])
            if not any(v) or (v*G*v-4)%8:continue
            assert all(c%2==0 for c in G*v)
            j=next(i for i,c in enumerate(v) if c)
            extension=block_diagonal_matrix(G,matrix(ZZ,[[-4]]))
            change=matrix.identity(QQ,19)
            change[j]=vector(QQ,list(v)+[1])/2
            geo=change*extension*change.transpose()
            assert abs(change.det())==QQ(1)/2 and geo.det()==468
            assert all(c.denominator()==1 for c in geo.list()) and all(geo[i,i]%2==0 for i in range(19))
            action=matrix.diagonal(QQ,[1]*18+[-1]);galois=change*action*change.inverse()
            assert all(c.denominator()==1 for c in galois.list()) and galois*galois==matrix.identity(QQ,19)
            assert galois*geo*galois.transpose()==geo
            glues.append({'fixed_numerator':[int(c) for c in v],'replace_basis_index':j,'geometric_Gram':[[int(c) for c in r] for r in geo.rows()], 'Galois_action':[[int(c) for c in r] for r in galois.rows()]})
        result={'outer_u':row['outer_u'],'basis':['F','O']+['I2_minus','I2_plus','I4_1','I4_2','I4_3']+['C'+str(i) for i in h['seed_indices']],
          'section_points':[[str(c) for c in P.xy()] for P in points],'section_profiles':profiles,
          'rational_NS_Gram':[[int(c) for c in r] for r in G.rows()], 'determinant':int(G.det()),'signature':[1,17],
          'U_split_change':[[int(c) for c in r] for r in Uchange.rows()], 'discriminant_group_invariants':elementary,
          'two_primary_kernel_dimension':int(kernel.dimension()),'eligible_index_two_geometric_glues':glues,
          'geometric_Gram_status':'UNIQUE_GLUE' if len(glues)==1 else 'UNRESOLVED_GLUE',
          'new_fibration':'NOT_CONSTRUCTED','rootless_frame_search':False}
        rows.append(result);print('u'+row['outer_u'],'NS det',G.det(),'disc group',elementary,'eligible glues',len(glues),flush=True)
    paths=[Path(__file__).resolve(),INPUT,helper,ART/'mestre_468_portable_replay_v2.json']
    return {'schema':'elliptic-curves.mestre-rational-ns-gram.v2','status':'PASS','rows':rows,
      'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
      'scope':'Exact full rational NS intersection matrices, using already proved rank and saturation. Complete finite index-two gluing possibilities with the known primitive anti class of square -4. No rootless-frame inspection, neighbour enumeration, transcendental marking or equation search.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();result=compute()
    if a.check:assert result==json.loads(OUT.read_text())
    else:
        if OUT.exists():raise FileExistsError('preserve NS Gram result')
        OUT.write_text(json.dumps(result,indent=2)+'\n')
    print('PASS6 RATIONAL NS MATRICES',flush=True)
