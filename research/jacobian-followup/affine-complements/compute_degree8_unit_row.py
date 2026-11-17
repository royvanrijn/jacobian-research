"""Reproduce PF2RULED28's degree-three unit row by exact rational linear algebra.
Mod1009 chooses independent Macaulay rows only. Exact FLINT Dixon solving
and the final full polynomial identity are the correctness gates.
"""
from pathlib import Path
from time import monotonic
import gzip,json
from flint import fmpq,fmpq_poly,fmpq_mat,nmod_mat
from verify_degree8_unit_row import HERE,CTX,primary_field,read_gz,encoded_polynomial,rational_polynomial


def monomials(maximum):
    for total in range(maximum+1):
        for a in range(total,-1,-1):
            for b in range(total-a,-1,-1):
                for c in range(total-a-b,-1,-1):yield (a,b,c,total-a-b-c)

def main():
    start=monotonic();primary=read_gz(HERE/'degree8_affine_reduction_certificates.json.gz');result=[]
    for old in primary['sectors']:
        k=old['multiplicities'][0];p8,source=primary_field(old);pcoeff=[fmpq(str(p8.nth(j))) for j in range(9)];pfield=fmpq_poly(pcoeff)
        generators=[]
        for poly in source:
            terms={}
            for ex,c in poly.to_dict().items():
                assert ex[0]<8
                terms.setdefault(tuple(ex[1:]),[fmpq(0) for _ in range(8)])[ex[0]]+=c
            generators.append({ex:tuple(c) for ex,c in terms.items()})
        prime=1009;root=990
        def residue(c):return int(c.numerator)%prime*pow(int(c.denominator)%prime,-1,prime)%prime
        def field_residue(c):return sum(residue(t)*pow(root,j,prime) for j,t in enumerate(c))%prime
        assert field_residue(pcoeff)==0
        mons=list(monomials(3));assert len(mons)==35
        columns=[{tuple(a+b for a,b in zip(ex,m)):c for ex,c in g.items()} for g in generators for m in mons]
        support=sorted(set().union(*(row.keys() for row in columns)),key=lambda ex:(sum(ex),ex));positions={ex:i for i,ex in enumerate(support)}
        modular=nmod_mat(len(support),len(columns),prime)
        for j,col in enumerate(columns):
            for ex,c in col.items():modular[positions[ex],j]=field_residue(c)
        rr,rank=modular.transpose().rref();assert rank==175
        rows=[next(j for j in range(rr.ncols()) if rr[i,j]) for i in range(rank)]
        assert len(set(rows))==175
        # Restriction of scalars expands each K coefficient to its8-by8
        # multiplication matrix on 1,y,...,y^7. No number-field elimination.
        size=175;Q=fmpq_mat(size*8,size*8);rhs=fmpq_mat(size*8,1);cache={}
        for i,row in enumerate(rows):
            ex=support[row];rhs[i*8,0]=int(ex==(0,0,0,0))
            for j,col in enumerate(columns):
                value=col.get(ex)
                if value is None:continue
                if value not in cache:
                    block=[]
                    for power in range(8):
                        product=fmpq_poly([fmpq(0)]*power+list(value))%pfield
                        block.append([product[a] for a in range(8)])
                    cache[value]=block
                block=cache[value]
                for b in range(8):
                    for a in range(8):
                        if block[b][a]:Q[i*8+a,j*8+b]=block[b][a]
        print(f'sector{k}: exact1400-by1400 rational Dixon solve',flush=True)
        solution=Q.solve(rhs,algorithm='dixon');assert Q*solution==rhs
        encode=lambda data:[[list(map(int,ex)),[str(c) for c in values]] for ex,values in data.items()]
        Cs=[[[list(ex),[str(solution[(i*35+j)*8+b,0]) for b in range(8)]] for j,ex in enumerate(mons)] for i in range(5)]
        Ds=[encode(g) for g in generators];p=rational_polynomial(list(map(str,pcoeff)))
        full=sum((encoded_polynomial(c)*encoded_polynomial(d) for c,d in zip(Cs,Ds)),CTX.constant(0))
        assert divmod(full-1,p)[1]==0
        result.append(dict(P8_coefficients=list(map(str,pcoeff)),variables=['q4','q5','q6','q7'],generators=Ds,
            coefficients=Cs,sector=k,multiplier_degree=3,
            discovery_method='degree3 Macaulay rows selected mod1009; restriction of scalars and exact FLINT Dixon; full characteristic-zero identity checked',
            macaulay_row_labels=[list(map(int,support[i])) for i in rows],Macaulay_field_matrix_size=175,rational_matrix_size=1400))
        print(f'sector{k}: full exact polynomial unit verified; elapsed{monotonic()-start:.3f}s',flush=True)
    out=HERE/'degree8_unit_row_certificates.json.gz'
    out.write_bytes(gzip.compress((json.dumps(dict(as_of='2026-11-14',sectors=result),indent=2)+'\n').encode(),mtime=0))
    print('wrote both exact unit-row certificates',flush=True)

if __name__=='__main__':main()
