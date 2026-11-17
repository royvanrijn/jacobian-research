#!/usr/bin/env python3
"""Compile bounded exact moments for the full normalized mixed WT/WZ2 chart.

No solver is invoked. A separate explicitly bounded Singular run can read the
export. The coefficient a remains an indeterminate, not a sampled value.
"""
from pathlib import Path
from math import factorial
import argparse,hashlib,json,time

G=[('r0',0,0),('u0',0,1),('v0',0,2),('s',0,3),('r1',1,0),('u1',1,1),('v1',1,2),('r2',2,0),('u2',2,1),('r3',3,0)]
VARIABLES=['a']+[x[0] for x in G]


def compositions(n,k):
    if k==1:
        yield(n,);return
    for i in range(n+1):
        for tail in compositions(n-i,k-1):yield(i,)+tail


def moments(m):
    out={}
    fact=[factorial(i) for i in range(3*m+1)]
    for n in range(m+1):
        k=m-n
        for counts in compositions(n,len(G)):
            zz=sum(c*g[1] for c,g in zip(counts,G))
            if k<zz or (k-zz)%2:continue
            ell=(k-zz)//2
            tt=sum(c*g[2] for c,g in zip(counts,G))
            denG=1
            for c in counts:denG*=fact[c]
            for j in range(k-ell+1):
                if (tt+j)%2:continue
                aw=k-ell-j
                den=denG*fact[ell]*fact[j]*fact[aw]
                assert fact[m]%den==0
                momT=fact[tt+j]//(2**((tt+j)//2)*fact[(tt+j)//2])
                value=(fact[m]//den)*fact[k]*momT
                if value:
                    mon=(aw,)+counts;out[mon]=out.get(mon,0)+value
    return out


def singular_expr(poly):
    terms=[]
    for powers,c in sorted(poly.items()):
        factors=[str(c)]+[name+(('^'+str(e)) if e!=1 else '') for name,e in zip(VARIABLES,powers) if e]
        terms.append('*'.join(factors))
    return '+'.join(terms) if terms else '0'


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cutoff',type=int,default=8);args=ap.parse_args()
    assert 4<=args.cutoff<=12
    start=time.monotonic();rows=[];polys=[]
    for m in range(1,args.cutoff+1):
        poly=moments(m);polys.append(poly);rows.append({'m':m,'terms':[[list(k),v] for k,v in sorted(poly.items())]})
    # Exact first-appearance coefficients for triangular elimination; these
    # declare all exceptional a values before a generic function-field solver.
    def derivative_coefficient(poly,name):
        j=VARIABLES.index(name);ans={}
        for mon,c in poly.items():
            if mon[j]:
                assert mon[j]==1
                key=list(mon);key[j]=0;ans[tuple(key)]=c
        return singular_expr(ans)
    coefficients=[derivative_coefficient(polys[m-1],name) for m,name in ((1,'r0'),(2,'r1'),(3,'r2'),(4,'r3'))]
    assert coefficients==['1','2*a','6+6*a^2','72*a+24*a^3'],coefficients
    rec={'chart':'P=W*(a+T+Z^2)+G(Z,T), deg G<=3','variables':VARIABLES,'cutoff':args.cutoff,'moments':rows,
         'triangular_coefficients':dict(zip(['r0','r1','r2','r3'],coefficients)),
         'generic_parameter_open':'a*(a^2+1)*(a^2+3)!=0',
         'excluded_parameter_strata':['a=0','a^2+1=0','a^2+3=0'],
         'coverage':'Exact full normalized chart moment polynomials; no solution or all-order conclusion is inferred.',
         'compiler_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    dest=Path(__file__).parent
    (dest/f'mixed_wt_parent_moments_{args.cutoff}.json').write_text(json.dumps(rec,separators=(',',':'))+'\n')
    lines=['ring rr=(0,a),(r3,r2,r1,r0,s,v1,u2,u1,v0,u0),dp;']
    lines += [f'poly M{i}='+singular_expr(poly)+';' for i,poly in enumerate(polys,1)]
    lines += ['poly R0=-(M1-r0);','poly N2=subst(M2,r0,R0);','poly R1=-(N2-2*a*r1)/(2*a);',
              'poly N3=subst(subst(M3,r0,R0),r1,R1);','poly R2=-(N3-6*(a^2+1)*r2)/(6*(a^2+1));',
              'poly N4=subst(subst(subst(M4,r0,R0),r1,R1),r2,R2);','poly R3=-(N4-24*a*(a^2+3)*r3)/(24*a*(a^2+3));']
    for m in range(5,args.cutoff+1):lines.append(f'poly J{m}=subst(subst(subst(subst(M{m},r0,R0),r1,R1),r2,R2),r3,R3);')
    lines += ['ideal IJ='+','.join('J'+str(m) for m in range(5,args.cutoff+1))+';',
              'ring r6=(0,a),(s,v1,u2,u1,v0,u0),dp;',
              'map phi=rr,0,0,0,0,s,v1,u2,u1,v0,u0;',
              'ideal I=phi(IJ);','print("COMPILED_REDUCED_INPUT");','size(I);',
              'option(redSB);','ideal B=std(I);','print("EXACT_QA_BASIS");','size(B);','if (size(B)<=30) { B; } else { print("BASIS_TOO_LARGE_FOR_LOG"); }',
              'print("QUOTIENT_DIMENSION");','dim(B);','quit;']
    (dest/f'mixed_wt_parent_{args.cutoff}.sing').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'cutoff':args.cutoff,'terms_per_moment':[len(p) for p in polys],'elapsed_seconds':time.monotonic()-start,'triangular_coefficients':rec['triangular_coefficients'],'open':rec['generic_parameter_open']},indent=2))


if __name__=='__main__':main()
