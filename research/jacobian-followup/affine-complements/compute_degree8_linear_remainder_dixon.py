"""Exact rational Dixon stage for PF2RULED31; invoked by the Sage producer."""
from pathlib import Path
import gzip,json,sys
from flint import fmpq,fmpq_mat

def main():
    data=json.loads(gzip.decompress(Path(sys.argv[1]).read_bytes()))
    A=fmpq_mat([[fmpq(c) for c in row] for row in data['rows']])
    b=fmpq_mat([[fmpq(c)] for c in data['rhs']])
    assert A.nrows()==A.ncols()==140 and b.nrows()==140
    solution=A.solve(b,algorithm='dixon')
    assert A*solution==b
    Path(sys.argv[2]).write_text(json.dumps([str(solution[i,0]) for i in range(140)])+'\n')
    print('PASS exact140 rational Dixon solve and full matrix identity',flush=True)
if __name__=='__main__':main()
