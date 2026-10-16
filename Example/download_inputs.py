"""Re-download the example inputs from input_manifest.json via SSH/scp."""
import hashlib
import json
from pathlib import Path
import subprocess


def checksum(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    root = Path(__file__).resolve().parent
    manifest = json.loads((root / 'input_manifest.json').read_text())
    print('Total input bytes:', sum(item['bytes'] for item in manifest['files']))
    for item in manifest['files']:
        destination = (root / item['local']).resolve()
        if not (destination.is_relative_to(root / 'data') or destination.is_relative_to(root.parent / 'Tcal')):
            raise ValueError('Manifest destination must be under Example/data or package Tcal')
        if destination.is_file() and destination.stat().st_size == item['bytes']:
            if checksum(destination) == item['sha256']:
                print('Already verified:', destination.name)
                continue
        destination.parent.mkdir(parents=True, exist_ok=True)
        partial = destination.with_name(destination.name + '.part')
        subprocess.run(['scp', '-o', 'BatchMode=yes',
                        manifest['host'] + ':' + item['remote'], str(partial)], check=True)
        if partial.stat().st_size != item['bytes'] or checksum(partial) != item['sha256']:
            raise RuntimeError('Verification failed: ' + str(partial))
        partial.replace(destination)
        print('Verified:', destination.name)


if __name__ == '__main__':
    main()
