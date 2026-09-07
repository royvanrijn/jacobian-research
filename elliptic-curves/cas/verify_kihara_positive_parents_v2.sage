#!/usr/bin/env sage-python
"""Portable section-height and semistable surface checks on three parents."""
import argparse,json
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,prod,matrix
def main(path,directory,indices):
    data=json.loads(path.read_text());R=PolynomialRing(QQ,'T');T=R.gen();K=R.fraction_field();S=PolynomialRing(R,'x');x=S.gen();out=[]
    for idx,row in enumerate(data['rows']):
        if idx not in indices:continue
        v=QQ(row['ratio']);centres=[QQ(0),(2*v*v+v+2)**2,2*(v+1)**2*(2*v*v+v+1),4*v*v-v+4,v*(2*v-1)*(2*v*v+4*v+5),4*v**4+8*v**3+9*v*v-2*v+2]
        norm=centres[1];centres=[a/norm for a in centres];assert list(map(str,centres))==row['normalized_roots']
        product=prod((x-c)**2-T*T for c in centres);g={6:R(1)}
        for j in range(5,-1,-1):g[j]=(product[6+j]-sum(g[i]*g[k] for i in g for k in g if i+k==6+j))/2
        square=S([g[i] for i in range(7)]);quartic=S([R(a) for a in row['quartic_coefficients']]);assert square**2-product==T*T*quartic
        e,d,c,b,a=quartic.list();A=-27*(12*a*e-3*b*d+c*c);B=-27*(72*a*c*e+9*b*c*d-27*a*d*d-27*b*b*e-2*c**3)
        assert list(map(str,A.list()))==row['raw_A'] and list(map(str,B.list()))==row['raw_B']
        delta=-16*(4*A**3+27*B**2);assert (A.degree(),B.degree())==(8,12) and delta.degree() in (18,20) and delta.gcd(A)==1
        infinity_n=24-delta.degree();multiple=infinity_n
        factors=list(delta.squarefree_decomposition());assert all(m in (1,2) for f,m in factors)
        node=-3*B[12]/(2*A[8]);assert node**3+A[8]*node+B[12]==0 and (3*node).is_square()
        pts=[tuple(map(K,p)) for p in row['quartic_points']];assert all(y*y==quartic(xx) for xx,y in pts)
        x0,y0=pts[0];tr=quartic(x+x0);ee,dd,cc,bb,aa=map(K,tr.list());E=EllipticCurve(K,[A,B]);images=[]
        for xx,y in pts[1:]:
            z=xx-x0
            if not z:
                assert y in (y0,-y0)
                if y==y0:images.append(E(0));continue
                wx=dd*dd/(4*y0*y0)-cc;wy=-dd*wx/(2*y0)-y0*bb
            else:
                wx=(2*y0*(y+y0)+dd*z)/z**2;wy=(2*(wx*wx-4*y0*y0*aa)*z-2*dd*wx-4*y0*y0*bb)/(4*y0)
            images.append(E([9*wx+3*cc,27*wy]))
        wx=dd*dd/(4*y0*y0)-cc;wy=-dd*wx/(2*y0)-y0*bb;involution=E([9*wx+3*cc,27*wy]);assert 6*involution==sum(images[:11],E(0));images[0]=involution
        # The infinity multiplicity is even and kills every finite I2 group.
        # Compute heights after that multiple, without local node charts.
        four=[multiple*P for P in images]
        def scaled_height(P):
            if not P:return QQ(0)
            n,d=P[0].numerator(),P[0].denominator();return QQ(4+max(d.degree(),n.degree()-4))/(multiple*multiple)
        H=matrix(QQ,12)
        for i,P in enumerate(four):
            H[i,i]=scaled_height(P)
            for j in range(i):
                plus=scaled_height(P+four[j]);minus=scaled_height(P-four[j]);assert plus+minus==2*(H[i,i]+H[j,j]);H[i,j]=H[j,i]=(plus-H[i,i]-H[j,j])/2
        rank=H.rank();inds=list(H.pivots());minor=H.matrix_from_rows_and_columns(inds,inds);assert minor.is_positive_definite()
        if 'height_Gram' in row:assert H==matrix(QQ,row['height_Gram'])
        counts=[]
        for p in (131,239,251):
            try:
                field=GF(p);rr=PolynomialRing(field,'u');u=rr.gen();av=rr(A);bv=rr(B);dv=-16*(4*av**3+27*bv**2)
            except (ValueError,ZeroDivisionError):counts.append({'prime':p,'status':'NOT_GOOD_IN_FIXED_MODEL'});continue
            if (av.degree(),bv.degree(),dv.degree())!=(8,12,delta.degree()) or dv.gcd(av)!=1:
                counts.append({'prime':p,'status':'NOT_GOOD_IN_FIXED_MODEL'});continue
            if [(f.degree(),m) for f,m in dv.squarefree_decomposition()]!=[(f.degree(),m) for f,m in factors]:
                counts.append({'prime':p,'status':'NOT_GOOD_IN_FIXED_MODEL'});continue
            # Rationality of the finite I2 base places must survive reduction.
            doubled=next((f for f,m in factors if m==2),R(1));assert all(f.degree()==1 for f,m in doubled.factor())
            total=infinity_n*p
            for t in field:
                aa,bb=av(t),bv(t)
                if 4*aa**3+27*bb**2:total+=EllipticCurve(field,[aa,bb]).cardinality()
                else:
                    nod=-3*bb/(2*aa);split=(3*nod).is_square();mult=dv.valuation(u-t);assert mult in (1,2)
                    total+=p+(0 if split else 2)+(mult-1)*p
            counts.append({'prime':p,'status':'PASS','surface_point_count':int(total)})
        if 'surface_counts' in row:
            assert [(a['prime'],a.get('surface_point_count')) for a in counts]==[(a['prime'],a.get('surface_point_count')) for a in row['surface_counts']]
        result={'ratio':str(v),'status':'PASS','finite_fibres':[{'multiplicity':int(m),'total_degree':int(f.degree())} for f,m in factors],'infinity_fibre':'split I'+str(infinity_n),'height_multiplier':int(multiple),'generic_section_span_rank':int(rank),'independent_indices':inds,'independent_height_determinant':str(minor.det()),'height_Gram':[list(map(str,r)) for r in H.rows()],'sections':[[str(c) for c in P] for P in images],'surface_counts':counts,'full_generic_rank':'UNKNOWN','full_NS':'UNKNOWN'}
        with (directory/('parent'+str(idx)+'.json')).open('x') as f:json.dump(result,f,indent=2);f.write('\n')
        out.append(result);print('PASS positive parent',v,'section rank',rank,'det',minor.det(),flush=True)
    with (directory/'result.json').open('x') as f:json.dump({'status':'PASS','rows':out},f,indent=2);f.write('\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output-directory',type=Path,required=True);p.add_argument('--indices',default='0,1,2');a=p.parse_args();main(a.input,a.output_directory,{int(i) for i in a.indices.split(',')})
