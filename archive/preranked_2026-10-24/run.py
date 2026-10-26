#!/usr/bin/env python3
"""Reproduce scientific figures and tables from frozen inputs, from any working directory."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--rscript',default='Rscript',help='Rscript executable with the packages in R_packages.txt')
    ap.add_argument('--validate-only',action='store_true',help='Check frozen inputs and existing results without recomputing')
    args=ap.parse_args()
    root=Path(__file__).resolve().parent
    stages=[('validation',[sys.executable,str(root/'scripts/validate_release.py')])]
    if not args.validate_only:
        stages=[('analysis',[sys.executable,str(root/'scripts/analyze.py')]),
                ('cell_fate',[args.rscript,str(root/'scripts/cell_fate.R')]),
                ('figures',[sys.executable,str(root/'scripts/render_results.py')]),
                ('figure_single_degron',[sys.executable,str(root/'scripts/render_results.py'),'--single-degron'])]
        stages.append(('validation',[sys.executable,str(root/'scripts/validate_release.py')]))
    for name,command in stages:
        print(f'Running {name}...',flush=True)
        logfile=root/'qc'/f'{name}.log'
        with logfile.open('w') as stream:
            result=subprocess.run(command,cwd=root,stdout=stream,stderr=subprocess.STDOUT,
                                  env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
        if result.returncode:
            print(logfile.read_text(),file=sys.stderr)
            raise SystemExit(f'{name} failed; see {logfile}')
    print('Complete. Main figure: results/figures/figure_1_chip48_and_cell_fate_single_degron.pdf')

if __name__=='__main__':main()
