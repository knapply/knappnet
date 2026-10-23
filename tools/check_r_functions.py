#!/usr/bin/env python3
"""Fails if two files in r/R define the same top-level function, since the
later definition would silently replace the earlier one.

    python3 tools/check_r_functions.py
"""

import re
import sys
from collections import defaultdict
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parent.parent
    defs = defaultdict(list)
    for path in sorted((root / "r" / "R").glob("*.R")):
        for n, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            m = re.match(r"^`?([A-Za-z.][\w.]*)`?\s*(<-|=)\s*function\b", line)
            if m:
                defs[m.group(1)].append(f"{path.relative_to(root)}:{n}")
    dups = {k: v for k, v in defs.items() if len(v) > 1}
    for name, where in sorted(dups.items()):
        print(f"{name} is defined more than once: {', '.join(where)}")
    if not dups:
        print(f"{len(defs)} R functions, none defined twice.")
    return 1 if dups else 0


if __name__ == "__main__":
    if {"-h", "--help"} & set(sys.argv[1:]):
        print(__doc__.strip())
        sys.exit(0)
    sys.exit(main())
