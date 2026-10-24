"""Builds the extension module: compiles the engine's C ABI as a static
library with Cargo, then links it into knappnet._core.

The Rust workspace is the parent directory of this one (override with
KNAPPNET_SOURCE). Cargo's target directory follows CARGO_TARGET_DIR as
usual, so the library is shared with the rest of the repository's builds.
The workspace has no third-party crates, so the build runs offline.
"""

import os
import subprocess
import sys
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("KNAPPNET_SOURCE", HERE.parent)).resolve()

# System libraries the Rust standard library needs when linked statically
# (as reported by `rustc --print native-static-libs`; a superset is harmless).
if sys.platform == "win32":
    SYSTEM_LIBS = ["ws2_32", "userenv", "bcrypt", "ntdll", "advapi32", "kernel32",
                   "synchronization"]
    EXTRA_LINK = []
elif sys.platform == "darwin":
    SYSTEM_LIBS = ["iconv"]
    EXTRA_LINK = ["-framework", "CoreFoundation"]
else:
    SYSTEM_LIBS = ["pthread", "dl", "m", "rt", "util"]
    EXTRA_LINK = []


def target_dir() -> Path:
    env = os.environ.get("CARGO_TARGET_DIR")
    return Path(env).resolve() if env else ROOT / "target"


class CargoBuildExt(build_ext):
    def run(self):
        cmd = ["cargo", "build", "--release", "--offline", "-p", "capi"]
        print(" ".join(cmd), flush=True)
        subprocess.run(cmd, cwd=ROOT, check=True)
        name = "capi.lib" if sys.platform == "win32" else "libcapi.a"
        lib = target_dir() / "release" / name
        if not lib.exists():
            raise RuntimeError(f"Cargo did not produce {lib}")
        for ext in self.extensions:
            ext.extra_objects = [str(lib)]
        super().run()


core = Extension(
    "knappnet._core",
    sources=["src/knappnet/_core.c"],
    include_dirs=[str(ROOT / "crates" / "capi" / "include")],
    libraries=SYSTEM_LIBS,
    extra_link_args=EXTRA_LINK,
    define_macros=[("Py_LIMITED_API", "0x030A0000")],
    py_limited_api=True,
)

setup(
    ext_modules=[core],
    cmdclass={"build_ext": CargoBuildExt},
    options={"bdist_wheel": {"py_limited_api": "cp310"}},
)
