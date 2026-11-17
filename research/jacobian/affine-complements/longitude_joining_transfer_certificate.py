"""Exact Z[t] certificate for genuine cables (r,k)=(36t+1,36t-5).

For t>=1 both parameters exceed one. Each covered cable-space component
with outer lattice (e,0),(0,2), e=1 or3, has H1=Z^2 on (m_new,L_old),
with m_old=r*m_new and L_new=r*L_old. Every operation below is unimodular
over Z[t]; substitution proves ALL integer t, not just finite samples.
"""
import sympy as S
from full_peripheral_linking import Flows


def certificate(e):
    if e not in (1,3):raise ValueError('this certificate covers e=1 or3 only')
    t=S.Symbol('t',integer=True);r=36*t+1;k=36*t-5;d=2*e
    cyc=lambda j:tuple((i+j)%d for i in range(d))
    f=Flows([cyc(2 if e==3 else 0),cyc(e),cyc((2 if e==3 else 0)+e)])
    m,l,z=f.generators
    def powerpoly(a,exponent):
        small=int(S.expand(exponent).subs(t,0))%6
        q=S.expand((exponent-small)/6)
        assert all(c.is_Integer for c in S.Poly(q,t).all_coeffs())
        b=f.power(a,small);c=f.power(a,6);assert c[0]==f.ident
        return b[0],[[S.expand(v+q*w) for v,w in zip(b[1][s],c[1][b[0][s]])] for s in range(d)]
    # (1-6t)*r+(6t)*k=1; the longitude exponent uses minus the second coefficient.
    assert S.expand((1-6*t)*r+6*t*k)==1
    mn=f.product(powerpoly(m,1-6*t),powerpoly(l,-6*t),powerpoly(z,6*t))
    fibre=f.product(powerpoly(m,k),powerpoly(l,r))
    ln=f.product(powerpoly(mn,-r*k),fibre)
    rels=[f.product(m,l,f.inv(m),f.inv(l)),f.product(powerpoly(z,r),f.inv(fibre))]
    assert all(rel[0]==f.ident for rel in rels)
    assert mn[0]==m[0] and ln[0]==l[0]
    par=list(range(d));tree=[]
    def rt(i):
        while par[i]!=i:i=par[i]
        return i
    for g,a in enumerate(f.generators):
        for i in range(d):
            j=a[0][i]
            if rt(i)!=rt(j):par[rt(i)]=rt(j);tree.append(g*d+i)
    assert len(tree)==d-1
    keep=[i for i in range(3*d) if i not in tree]
    proj=lambda row:[S.expand(row[i]) for i in keep]
    rows=[proj(row) for rel in rels for row in rel[1]]
    ports=[proj(f.power(m,e)[1][0]),proj(f.power(l,2)[1][0]),
           proj(f.power(mn,e)[1][0]),proj(f.power(ln,2)[1][0])]
    original_rows=[row[:] for row in rows]
    original_ports=[row[:] for row in ports]
    ops=([('pivot',0,0,1),('row_subtract',2,1),('pivot',2,1,-1),
          ('column_add',0,1,S.Integer(1)),('column_add',1,0,6*t-1),('pivot',1,0,-1)] if e==1 else
         [('pivot',0,5,1),('pivot',1,1,-1),('pivot',1,1,1),('pivot',1,0,1),('pivot',0,0,1),
          ('row_subtract',3,1),('pivot',3,2,-1),('row_subtract',4,1),('pivot',4,3,-1),
          ('row_subtract',4,2),('pivot',4,2,-1),('row_subtract',2,1),('pivot',2,1,-1),
          ('row_subtract',2,1),('pivot',2,0,-1),('column_add',0,1,18*t),('pivot',1,1,1)])
    width=len(keep)
    for op in ops:
        if op[0]=='row_subtract':
            _,i,j=op;rows[i]=[S.expand(a-b) for a,b in zip(rows[i],rows[j])]
        elif op[0]=='column_add':
            _,i,j,c=op
            assert all(v.is_Integer for v in S.Poly(c,t).all_coeffs())
            for row in rows+ports:row[j]=S.expand(row[j]+c*row[i])
        else:
            _,i,j,unit=op;pivot=rows[i];assert pivot[j]==unit and abs(unit)==1
            def reduce(row):return [S.expand(row[h]-row[j]*unit*pivot[h]) for h in range(width) if h!=j]
            rows=[reduce(row) for h,row in enumerate(rows) if h!=i]
            ports=[reduce(row) for row in ports];width-=1
    assert width==2 and all(x==0 for row in rows for x in row)
    det=S.expand(S.Matrix([ports[2],ports[1]]).det());assert det in (-1,1)
    assert all(S.expand(a-r*b)==0 for a,b in zip(ports[0],ports[2]))
    assert all(S.expand(a-r*b)==0 for a,b in zip(ports[3],ports[1]))
    return dict(normal_index=e,longitudinal_degree=2,shear=0,winding=str(r),cable_parameter=str(k),
                coefficient_ring='Z[t]',parameter_range='all integers t>=1',
                original_rows=[[str(x) for x in row] for row in original_rows],
                original_ports=[[str(x) for x in row] for row in original_ports],
                unimodular_operations=[[str(x) for x in op] for op in ops],
                final_rows=[[str(x) for x in row] for row in rows],
                final_ports=[[str(x) for x in row] for row in ports],basis_determinant=int(det),
                H1_free_rank=2,H1_torsion=[],basis=['m_new','L_old'],
                all_parameter_port_identities=['m_old=r*m_new','L_new=r*L_old'],
                parameter_samples_used_as_proof=False)


if __name__=='__main__':
    import json
    print(json.dumps([certificate(1),certificate(3)],indent=2))
