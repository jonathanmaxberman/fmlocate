"""Where FMLocate keeps downloaded licence data, map tiles and saved measurements.

Override with FMLOCATE_CACHE / FMLOCATE_DATA, or --cache-dir on the command line.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path


def _default(kind: str) -> Path:
    home = Path.home()
    if sys.platform == "darwin":
        return home / "Library" / ("Caches" if kind == "cache" else "Application Support") / "fmlocate"
    xdg = os.environ.get("XDG_CACHE_HOME" if kind == "cache" else "XDG_DATA_HOME")
    return Path(xdg) / "fmlocate" if xdg else home / (".cache" if kind == "cache" else ".local/share") / "fmlocate"


CACHE_DIR = Path(os.environ.get("FMLOCATE_CACHE") or _default("cache"))
DATA_DIR = Path(os.environ.get("FMLOCATE_DATA") or _default("data"))


def set_cache_dir(path: str | Path) -> None:
    global CACHE_DIR
    CACHE_DIR = Path(path).expanduser()


def cache(*parts: str) -> Path:
    return CACHE_DIR.joinpath(*parts)


def measurements() -> Path:
    return DATA_DIR / "measurements"
