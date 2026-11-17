"""Exact two-meridian Artin relation in every PSL2(11) A_k stratum."""
from collections import Counter
from psl11_curve_spectrum import closure
from local_image import fixed_sheets,subgroup_kind

def power(a,n,identity,table):
    out=identity
    while n:
        if n&1:out=table[out][a]
        a=table[a][a];n//=2
    return out

def alternating(a,b,m,identity,table):
    out=identity
    for j in range(m):
        out=table[out][a if j%2==0 else b]
    return out

def pair_orbits(identity,table,inverses,orders,perms):
    """Disjoint simultaneous-conjugacy orbits of all429^2 actual meridian pairs."""
    allowed={i for i,o in enumerate(orders) if o in [2,3,5]}
    remaining={(a,b) for a in allowed for b in allowed};rows=[]
    while remaining:
        a,b=min(remaining)
        orb={(table[table[g][a]][inverses[g]],table[table[g][b]][inverses[g]])
             for g in range(len(table))}
        if not orb<=remaining:raise AssertionError("pair orbits overlap")
        remaining-=orb
        H=closure((a,b),identity,table)
        marks=[len(fixed_sheets(H,p)) for p in perms]
        prod=table[a][b];order=orders[prod];period=2*order
        residues=[m for m in range(period)
                  if alternating(a,b,m,identity,table)==alternating(b,a,m,identity,table)]
        rows.append(dict(a=a,b=b,orbit_size=len(orb),orders=[orders[a],orders[b]],
            product_order=order,period=period,residues=residues,
            subgroup=sorted(H),kind=subgroup_kind(H,orders,marks),marks=marks))
    return rows

def satisfies(row,m):
    return m%row['period'] in row['residues']

def spectrum(rows,m,allowed_orders=(2,3,5)):
    answer=Counter()
    for row in rows:
        if set(row['orders'])<=set(allowed_orders) and satisfies(row,m):
            answer[row['kind']]+=row['orbit_size']
    return dict(sorted(answer.items()))

def full_group_orders(rows,m,allowed_orders=(2,3,5)):
    answer=Counter()
    for row in rows:
        if row['kind']=='G' and set(row['orders'])<=set(allowed_orders) and satisfies(row,m):
            answer[tuple(row['orders']+[row['product_order']])]+=row['orbit_size']
    return answer
