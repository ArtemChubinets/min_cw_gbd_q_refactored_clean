#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Fetch the pinned ACM CALGO 994 release, apply the witness patch, and build."""

import argparse
import hashlib
from io import BytesIO
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile
from zipfile import ZipFile

URL = 'https://calgo.acm.org/994.zip'
SHA256 = 'f08357bcbd9b48abd42107a7619b3e07d95986a492f7ce7a77567f03e4c46d72'
HERE = Path(__file__).resolve().parent
# Keep CALGO's generation path and the original kernels it references. The
# upstream command-line program, datasets and configure system are not needed.
SOURCES = '''alg_basic alg_common alg_gray alg_optimized alg_saved
alg_saved_unrolled alg_stack computing_kernels generations luf2 utils'''.split()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prefix', type=Path,
                        default=Path.home() / '.local/share/qgbd-tools/acm994')
    parser.add_argument('--archive', type=Path,
                        help='use this local copy of 994.zip; still verify its SHA-256')
    args = parser.parse_args()
    prefix = args.prefix.expanduser().resolve()
    bin_dir = prefix / 'bin'
    license_dir = prefix / 'share/licenses/acm994-codeword'
    cc = shlex.split(os.environ.get('CC', 'cc'))
    cflags = shlex.split(os.environ.get('CFLAGS', '-O3'))
    with tempfile.TemporaryDirectory(prefix='acm994-build-') as directory:
        work = Path(directory)
        archive = args.archive
        if archive is None:
            archive = work / '994.zip'
            subprocess.run(['curl', '--proto', '=https', '-fsSL', '--retry', '2',
                            '--output', str(archive), URL], check=True)
        archive_data = archive.read_bytes()
        digest = hashlib.sha256(archive_data).hexdigest()
        if digest != SHA256:
            parser.exit(1, f'SHA-256 mismatch for {archive}: {digest}; expected {SHA256}\n')
        source_names = {name + '.c' for name in SOURCES}
        with ZipFile(BytesIO(archive_data)) as upstream:
            for member in upstream.namelist():
                parent, _, name = member.rpartition('/')
                if parent == '994/Src/src' and (name.endswith('.h') or name in source_names):
                    (work / name).write_bytes(upstream.read(member))
            license_text = upstream.read('994/Doc/LICENSE')
        subprocess.run(['patch', '--batch', '--forward', '--fuzz=0', '-p1',
                        '-i', str(HERE / 'codeword.patch')], cwd=work, check=True)
        # The witness patch is for CALGO's scalar saved-combination kernel.
        (work / 'config.h').write_text('/* Scalar CALGO 994 configuration. */\n')
        executable = work / 'acm994-codeword'
        subprocess.run(cc + cflags + ['-std=gnu11', '-fopenmp', '-I', str(work),
                        '-I', str(HERE), '-o', str(executable),
                        str(HERE / 'codeword_cli.c')]
                       + [s + '.c' for s in SOURCES] + ['-lm'], cwd=work, check=True)
        smoke = subprocess.run([str(executable)], input=b'1 1\n1\n', check=True,
                               stdout=subprocess.PIPE)
        if smoke.stdout != b'RESULT 1 1\n':
            parser.exit(1, 'built executable failed the protocol smoke test\n')
        bin_dir.mkdir(parents=True, exist_ok=True)
        license_dir.mkdir(parents=True, exist_ok=True)
        (license_dir / 'LICENSE').write_bytes(license_text)
        # Install atomically so a concurrent reader sees a complete executable.
        with tempfile.TemporaryDirectory(prefix='.install-', dir=bin_dir) as directory:
            staged = Path(directory) / executable.name
            shutil.copy2(executable, staged)
            staged.replace(bin_dir / executable.name)
    print(bin_dir / executable.name)


if __name__ == '__main__':
    main()
