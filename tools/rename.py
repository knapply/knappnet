#!/usr/bin/env python3
"""Rename the project, its GitHub repository, or both.

The current name and repository are read from project.toml. Every text file
in the repository is rewritten, along with any file or directory names that
contain the old name.

    python3 tools/rename.py newname                    # package name only
    python3 tools/rename.py --repo owner/newrepo       # GitHub location only
    python3 tools/rename.py newname --repo owner/new   # both
    python3 tools/rename.py newname --dry-run          # show what would change
    python3 tools/rename.py --self-test                # verify renaming works

After renaming the repository on GitHub, also update the remote:
    git remote set-url origin https://github.com/owner/newrepo.git
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SKIP_DIRS = {".git", "target", ".cargo", "__pycache__", ".Rproj.user", "node_modules"}
# R package names: letters, digits, and dots; must start with a letter.
NAME_RE = re.compile(r"^[A-Za-z][A-Za-z0-9.]*[A-Za-z0-9]$")
REPO_RE = re.compile(r"^[A-Za-z0-9-]+/[A-Za-z0-9._-]+$")


def current(root: Path) -> tuple[str, str]:
    data = tomllib.loads((root / "project.toml").read_text())["project"]
    return data["name"], data["repository"]


def text_files(root: Path):
    for path in sorted(root.rglob("*")):
        if any(part in SKIP_DIRS for part in path.relative_to(root).parts):
            continue
        if not path.is_file():
            continue
        try:
            yield path, path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, PermissionError):
            continue


def rename(root: Path, new_name: str | None, new_repo: str | None, dry_run: bool) -> int:
    old_name, old_repo = current(root)
    subs: list[tuple[re.Pattern, str]] = []
    # The repository substitution runs first because it may contain the name.
    if new_repo and new_repo != old_repo:
        subs.append((re.compile(re.escape(old_repo)), new_repo))
    if new_name and new_name != old_name:
        subs.append((re.compile(rf"(?<![A-Za-z0-9]){re.escape(old_name)}(?![A-Za-z0-9])"),
                     new_name))
    if not subs:
        print("Nothing to change.")
        return 0

    changed = 0
    for path, text in text_files(root):
        new = text
        for pattern, repl in subs:
            new = pattern.sub(repl, new)
        if new != text:
            changed += 1
            print(f"  edit   {path.relative_to(root)}")
            if not dry_run:
                path.write_text(new, encoding="utf-8")

    if new_name and new_name != old_name:
        # Deepest paths first so parents are renamed after their children.
        for path in sorted(root.rglob(f"*{old_name}*"), key=lambda p: -len(p.parts)):
            if any(part in SKIP_DIRS for part in path.relative_to(root).parts):
                continue
            target = path.with_name(path.name.replace(old_name, new_name))
            print(f"  move   {path.relative_to(root)} -> {target.relative_to(root)}")
            if not dry_run:
                path.rename(target)
            changed += 1

    print(f"{'Would change' if dry_run else 'Changed'} {changed} item(s).")
    return changed


def self_test() -> int:
    """Renames a scratch copy and checks that no trace of the old name remains."""
    old_name, old_repo = current(ROOT)
    probe, probe_repo = "zzrenameprobe", "probe-owner/zzrenameprobe-repo"
    with tempfile.TemporaryDirectory() as tmp:
        copy = Path(tmp) / "repo"
        shutil.copytree(ROOT, copy, ignore=shutil.ignore_patterns(*SKIP_DIRS))
        rename(copy, probe, probe_repo, dry_run=False)
        leftovers = [str(p.relative_to(copy)) for p, t in text_files(copy)
                     if old_name in t or old_repo in t]
        leftovers += [str(p.relative_to(copy)) for p in copy.rglob(f"*{old_name}*")]
        required = {
            "r/DESCRIPTION": f"Package: {probe}",
            "r/NAMESPACE": f"useDynLib({probe},",
            "r/src/init.c": f"R_init_{probe}(",
            "r/tests/testthat.R": f'test_check("{probe}")',
            "project.toml": f'repository = "{probe_repo}"',
        }
        missing = [f for f, s in required.items() if s not in (copy / f).read_text()]
        sync = subprocess.run([sys.executable, str(copy / "tools/sync_r_sources.py"), "--check"],
                              capture_output=True, text=True)
        ok = not leftovers and not missing and sync.returncode == 0
        for f in leftovers:
            print(f"  old name still present in {f}")
        for f in missing:
            print(f"  expected the new name in {f}")
        if sync.returncode:
            print(sync.stdout)
    print("Rename self-test passed." if ok else "Rename self-test FAILED.")
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("name", nargs="?", help="new package name")
    ap.add_argument("--repo", help="new GitHub repository, as owner/name")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        return self_test()
    if not args.name and not args.repo:
        ap.error("give a new name, --repo, or both")
    if args.name and not NAME_RE.match(args.name):
        ap.error("package names must start with a letter and contain only letters, "
                 "digits, and dots")
    if args.repo and not REPO_RE.match(args.repo):
        ap.error("--repo must look like owner/name")
    rename(ROOT, args.name, args.repo, args.dry_run)
    if not args.dry_run:
        print("Next: review with `git diff`, then run the tests.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
