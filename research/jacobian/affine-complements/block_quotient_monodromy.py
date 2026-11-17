"""Small exact PF122 controls; all-degree deductions are written in the note."""


def induced_blocks(p,blocks):
    where={x:i for i,B in enumerate(blocks) for x in B}
    result=[]
    for B in blocks:
        images={where[p[x]] for x in B}
        assert len(images)==1
        result.append(images.pop())
    return result


def support(p):
    return {i for i,j in enumerate(p) if i!=j}


def binary_rank(vectors):
    basis={}
    for z in vectors:
        while z:
            k=z.bit_length()-1
            if k in basis:z^=basis[k]
            else:
                basis[k]=z
                break
    return len(basis)
