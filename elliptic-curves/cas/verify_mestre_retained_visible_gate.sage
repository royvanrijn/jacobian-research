#!/usr/bin/env sage-python
"""Check all retained visible ceilings and join the two exact frame replays."""
import argparse,json
from pathlib import Path
from sage.all import QQ,ZZ,matrix,vector
def main(path):
    b=json.loads(path.read_text());ns=b['ns'];G=matrix(QQ,ns['rational_NS_Gram']);inv=G.inverse();H=matrix(QQ,b['heights']['seed_height_gram'])
    std=matrix.identity(QQ,18).rows();F,O=std[:2]
    def section(w):
        if not any(w):return O
        comp=[sum(w[i]*ns['section_profiles'][i]['components'][j] for i in range(11))%n for j,n in enumerate((2,2,4))]
        corr=QQ(comp[0]+comp[1])/2+QQ(comp[2]*(4-comp[2]))/4;oo=(w*H*w-4+corr)/2
        intersections=[]
        for i,p in enumerate(ns['section_profiles']):
            other=p['components'];c=QQ(comp[0]*other[0]+comp[1]*other[1])/2+min(comp[2],other[2])-QQ(comp[2]*other[2])/4
            intersections.append(2+oo+p['zero_section_intersection']-(w*H)[i]-c)
        v=inv*vector(QQ,[1,oo,*comp[:2],*[int(comp[2]==i) for i in (1,2,3)],*intersections])
        assert all(c in ZZ for c in v) and v*G*v==-2
        return v,comp
    visible=[vector(QQ,v) for v in b['visible_curve_classes']]
    assert visible[:9]==[O]+std[2:7]+[F-std[2],F-std[3],F-sum(std[4:7])]
    for v in visible[9:]:assert section(vector(ZZ,v[7:]))[0]==v
    V=matrix(QQ,visible);VG=V*G
    exact={tuple(r['fibre_class']):r['generic_Q_MW_rank'] for r in b['prior_frames']+b['new_frames']}
    assert len(exact)==203
    bases={r['source_index']:vector(QQ,r['class']) for r in b['base_classes']}
    assert list(bases[27])==[3,2,0,0,-1,-2,-1,1,-1,0,0,0,0,0,0,0,0,0]
    assert list(bases[35])==[1,1,0,0,-1,-2,-1,1,0,0,0,1,-1,0,0,0,0,0]
    seen=set();max_bound=0
    for row in b['admission_rows']:
        base=bases[row['source_index']];w=vector(ZZ,row['subtract_section_word']);S,sc=section(w);assert base*G*S==2
        bc=base*G;comp=bc[2:7];full=[[2-comp[0],comp[0]],[2-comp[1],comp[1]],[2-sum(comp[2:]),*comp[2:]]]
        rotated=[[a[(j+c)%len(a)] for j in range(len(a))] for a,c in zip(full,sc)]
        cmp=[rotated[0][1],rotated[1][1],*rotated[2][1:]]
        cross=[]
        for j in range(11):
            ss=w+vector(ZZ,[int(i==j) for i in range(11)]);point=O if not any(ss) else section(ss)[0]
            cross.append(base*G*point)
        v=inv*vector(QQ,[2,2,*cmp,*cross]);assert list(v)==row['bisection_class'] and v*G*v==-2
        f=O+v;assert list(f)==row['fibre_class'] and f*G*f==0
        values=VG*f;assert min(values)>=0
        zero=[i for i,c in enumerate(values) if c==0];one=[i for i,c in enumerate(values) if c==1]
        assert zero==row['visible_vertical_indices'] and one==row['visible_section_indices'] and one
        ceiling=17-matrix(QQ,[visible[i] for i in zero]+[v]).rank();assert ceiling==row['visible_MW_upper_bound']
        final=exact[tuple(f)] if ceiling>11 else ceiling;assert final<=11;max_bound=max(max_bound,final)
        assert tuple(f) not in seen;seen.add(tuple(f))
    assert len(seen)==554
    print('PASS554 distinct fibre classes;203 exact-frame inputs;351 visible ceilings; all MW_Q<=11',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);a=p.parse_args();main(a.input)
