"""Build a self-hosted Pages artifact with hash-pinned WASM and viewer assets."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
import tarfile
from urllib.request import Request, urlopen
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def download(url, expected, cache, name):
    path = cache/name
    if not path.exists():
        print(f'Downloading {name}', flush=True)
        with urlopen(Request(url, headers={'User-Agent': 'PocketSeekR-site-builder'}), timeout=120) as response:
            data = response.read()
        if sha256(data).hexdigest() != expected:
            raise ValueError(f'Asset hash mismatch: {name}')
        path.write_bytes(data)
    if sha256(path.read_bytes()).hexdigest() != expected:
        raise ValueError(f'Cached asset hash mismatch: {name}')
    return path


def build(output, cache):
    manifest = json.loads((ROOT/'web/runtime-lock.json').read_text())
    output.mkdir(parents=True, exist_ok=True)
    cache.mkdir(parents=True, exist_ok=True)
    for name in ['index.html', 'style.css', 'app.js', 'worker.js', 'header-banner.svg', 'favicon.svg']:
        shutil.copy2(ROOT/'web'/name, output/name)
    (output/'.nojekyll').write_text('')
    runtime = output/'runtime'
    runtime.mkdir(exist_ok=True)
    pyodide = manifest['pyodide']
    archive = download(pyodide['archive_url'], pyodide['archive_sha256'], cache, 'pyodide-core.tar.bz2')
    with tarfile.open(archive) as tar:
        for name in ['pyodide.mjs', 'pyodide.asm.mjs', 'pyodide.asm.wasm', 'python_stdlib.zip', 'pyodide-lock.json']:
            (runtime/name).write_bytes(tar.extractfile('pyodide/'+name).read())
    lock = json.loads((runtime/'pyodide-lock.json').read_text())
    for name, package in pyodide['packages'].items():
        upstream = lock['packages'][name]
        for key in ['file_name', 'sha256', 'version', 'depends']:
            if upstream[key] != package[key]:
                raise ValueError(f'Package lock mismatch: {name}/{key}')
        if any(dep not in pyodide['packages'] for dep in package['depends']):
            raise ValueError(f'Unbundled dependency: {name}')
        path = download(pyodide['package_base_url']+package['file_name'], package['sha256'], cache, package['file_name'])
        shutil.copy2(path, runtime/package['file_name'])
    viewer = manifest['viewer']
    archive = download(viewer['archive_url'], viewer['archive_sha256'], cache, '3dmol.tgz')
    vendor = output/'vendor'
    vendor.mkdir(exist_ok=True)
    with tarfile.open(archive) as tar:
        for source, destination in [
            ('package/build/3Dmol-min.js', '3Dmol-min.js'),
            ('package/build/3Dmol-min.js.LICENSE.txt', '3Dmol-license.txt'),
            ('package/LICENSE', '3Dmol-LICENSE'),
        ]:
            (vendor/destination).write_bytes(tar.extractfile(source).read())
    with zipfile.ZipFile(output/'pocketseekr.zip', 'w', zipfile.ZIP_DEFLATED) as zip:
        for path in sorted((ROOT/'src/pocketseekr').rglob('*.py')):
            zip.write(path, path.relative_to(ROOT/'src'))
    # A real-structure RNA-only example, serialized by the same parser used in inference.
    import sys
    sys.path.insert(0, str(ROOT/'src'))
    from pocketseekr.browser import atoms_to_pdb
    from pocketseekr.utils.io import _extract_rna_atoms, parse_pdb_atoms
    examples = output/'examples'
    examples.mkdir(exist_ok=True)
    source = ROOT/'tests/data/1F1T.pdb'
    (examples/'1F1T_rna.pdb').write_text(atoms_to_pdb(_extract_rna_atoms(parse_pdb_atoms(source), source)))
    shutil.copy2(ROOT/'examples/synthetic_shell.pdb', examples/'synthetic_shell.pdb')
    shutil.copy2(ROOT/'docs/pocket_v5_algorithms.pdf', output/'algorithms.pdf')
    shutil.copy2(ROOT/'web/runtime-lock.json', output/'runtime-lock.json')
    print(f'Built {output}: {sum(p.stat().st_size for p in output.rglob("*") if p.is_file()) / 1024**2:.1f} MiB', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'outputs/site')
    parser.add_argument('--cache', type=Path, default=ROOT/'outputs/web-assets')
    args = parser.parse_args()
    build(args.output, args.cache)
