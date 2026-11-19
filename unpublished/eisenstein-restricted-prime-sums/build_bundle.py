#!/usr/bin/env python3
"""Build the research or TeX-only submission archive and verify PDF rebuilds."""
from pathlib import Path
from hashlib import sha256
from zipfile import ZipFile, ZipInfo, ZIP_DEFLATED
import argparse,json,os,subprocess,tempfile

def digest(path):
    return sha256(path.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify-rebuild',action='store_true')
    parser.add_argument('--hexagon',action='store_true',
                        help='Build hexagon-source.zip containing TeX files only')
    args=parser.parse_args()
    root=Path(__file__).resolve().parent
    names=('main.tex','references.tex') if args.hexagon else (
        'main.tex','references.tex','README.md','metadata.json','verification.json','build_bundle.py')
    files=[root/name for name in names]
    files+=sorted((root/'sections').glob('*.tex'))
    if not args.hexagon:
        files+=sorted(path for path in (root/'certificates').rglob('*')
                      if path.is_file() and '__pycache__' not in path.parts
                      and path.suffix in ('.py','.java','.json','.pdf'))
    for path in files:
        if not path.is_file() or not path.resolve().is_relative_to(root):
            raise RuntimeError('Missing or external source: '+str(path))
    target=root/('hexagon-source.zip' if args.hexagon else 'source.zip')
    with ZipFile(target,'w') as archive:
        for path in sorted(files):
            info=ZipInfo(path.relative_to(root).as_posix(),(1980,1,1,0,0,0))
            info.compress_type=ZIP_DEFLATED
            info.external_attr=0o100644<<16
            archive.writestr(info,path.read_bytes())
    result={'archive':target.name,'members':len(files),'archive_sha256':digest(target)}
    if args.hexagon:
        assert all(path.suffix=='.tex' for path in files)
        result['contents']='TeX only; no certificate code, receipts or build scripts'
    if args.verify_rebuild:
        with tempfile.TemporaryDirectory(prefix='exact-norms-source-') as temp:
            dest=Path(temp)
            with ZipFile(target) as archive:
                for name in archive.namelist():
                    if not (dest/name).resolve().is_relative_to(dest):
                        raise RuntimeError('Unsafe archive path')
                archive.extractall(dest)
            environment=dict(os.environ,SOURCE_DATE_EPOCH='1795046400',FORCE_SOURCE_DATE='1')
            build=subprocess.run(['latexmk','-pdf','-jobname=paper','-interaction=nonstopmode','-halt-on-error','main.tex'],cwd=dest,env=environment,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
            if build.returncode:
                raise RuntimeError(build.stdout[-8000:])
            log=(dest/'paper.log').read_text(errors='replace')
            if any(x in log for x in ('undefined references','undefined citations','Overfull','LaTeX Warning:')):
                raise RuntimeError('Clean rebuild diagnostics: '+chr(10).join(line for line in log.splitlines() if any(t in line for t in ('Warning','undefined','Overfull'))))
            original=digest(root/'paper.pdf'); rebuilt=digest(dest/'paper.pdf')
            if original!=rebuilt:
                raise RuntimeError('PDF is not byte-identical: '+original+' != '+rebuilt)
            result['clean_source_rebuild']='passed; byte-identical PDF'
            result['pdf_sha256']=original
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    main()
