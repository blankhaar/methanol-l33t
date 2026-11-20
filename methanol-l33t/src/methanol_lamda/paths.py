"""Portable project-path helpers used by commands and committed manifests."""

from __future__ import annotations

from pathlib import Path


def resolve_from_root(root: str | Path, value: str | Path) -> Path:
    """Resolve a user path, interpreting relative values from the project root."""
    project_root = Path(root).expanduser().resolve()
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = project_root / path
    return path.resolve()


def portable_path(root: str | Path, value: str | Path) -> str:
    """Return a POSIX project-relative path when *value* is inside *root*."""
    project_root = Path(root).expanduser().resolve()
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = project_root / path
    path = path.resolve()
    try:
        return path.relative_to(project_root).as_posix()
    except ValueError:
        return str(path)
