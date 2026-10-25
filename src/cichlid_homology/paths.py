"""Where outputs go: ~/share, so the home directory stays clean."""
from __future__ import annotations

import time
from pathlib import Path

from . import __version__


def default_out(stem: str, ext: str = "", base: str = "~/share", now: str | None = None) -> Path:
    """~/share/<stem>_<version>_<timestamp><ext>; nothing is created here."""
    stamp = now or time.strftime("%Y%m%d-%H%M%S")
    return Path(base).expanduser() / f"{stem}_{__version__}_{stamp}{ext}"
