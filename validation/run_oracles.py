#!/usr/bin/env python3
"""Regenerate golden values from every oracle.

Runs each Python oracle script in validation/oracles/ (files ending in
`_oracle.py`) and, unless --python-only is given, each R oracle script
(files ending in `.R`). Scripts that need a package that is not installed
report it and are skipped, so a partial environment still produces every
golden file it can.

    python3 validation/run_oracles.py
    python3 validation/run_oracles.py --python-only
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def main() -> int:
    python_only = "--python-only" in sys.argv[1:]
    failed = []
    for script in sorted((HERE / "oracles").glob("*_oracle.py")):
        print(f"== {script.name}")
        if subprocess.run([sys.executable, str(script)], cwd=ROOT).returncode:
            failed.append(script.name)
    if not python_only:
        rscript = shutil.which("Rscript")
        for script in sorted((HERE / "oracles").glob("*.R")):
            print(f"== {script.name}")
            if rscript is None:
                print("   Rscript not found; skipped")
                continue
            if subprocess.run([rscript, str(script)], cwd=ROOT).returncode:
                failed.append(script.name)
    if failed:
        print("Oracle scripts failed: " + ", ".join(failed))
        return 1
    return 0


if __name__ == "__main__":
    if {"-h", "--help"} & set(sys.argv[1:]):
        print(__doc__.strip())
        sys.exit(0)
    sys.exit(main())
