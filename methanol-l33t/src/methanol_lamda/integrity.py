"""Cryptographic bindings between normalized data and immutable raw artifacts."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Mapping


class RawBindingError(ValueError):
    """A completeness assertion is not bound to a pinned raw artifact."""


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _checksum_pins(root: Path) -> dict[str, str]:
    path = root / "data" / "raw" / "SHA256SUMS"
    if not path.exists():
        raise RawBindingError("data/raw/SHA256SUMS is missing")
    pins: dict[str, str] = {}
    for number, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue
        fields = line.split(maxsplit=1)
        if len(fields) != 2 or len(fields[0]) != 64:
            raise RawBindingError(f"malformed SHA256SUMS line {number}")
        relative = fields[1].lstrip("*")
        pins[relative] = fields[0].lower()
    return pins


def verify_raw_binding(root: str | Path, artifact: Mapping[str, object]) -> None:
    """Require one artifact under data/raw to match its declaration and pin."""
    root = Path(root).resolve()
    relative = str(artifact.get("path", ""))
    declared = str(artifact.get("sha256", "")).lower()
    if not relative or len(declared) != 64:
        raise RawBindingError("raw artifact needs path and a full SHA-256")
    target = (root / relative).resolve()
    raw_root = (root / "data" / "raw").resolve()
    try:
        target.relative_to(raw_root)
    except ValueError as exc:
        raise RawBindingError(f"raw artifact escapes data/raw: {relative}") from exc
    if not target.is_file():
        raise RawBindingError(f"raw artifact is absent: {relative}")
    actual = sha256_path(target)
    if actual != declared:
        raise RawBindingError(f"raw artifact checksum mismatch: {relative}")
    pinned = _checksum_pins(root).get(relative)
    if pinned != declared:
        raise RawBindingError(f"raw artifact is not pinned identically in SHA256SUMS: {relative}")


def verify_pinned_sha256(root: str | Path, declared: str) -> None:
    """Require a hash to resolve to at least one existing, identically pinned raw file."""
    root = Path(root).resolve()
    declared = declared.lower()
    if len(declared) != 64:
        raise RawBindingError("a full raw SHA-256 is required")
    matches = [path for path, digest in _checksum_pins(root).items() if digest == declared]
    if not matches:
        raise RawBindingError("raw SHA-256 is absent from SHA256SUMS")
    failures: list[str] = []
    for relative in matches:
        try:
            verify_raw_binding(root, {"path": relative, "sha256": declared})
        except RawBindingError as exc:
            failures.append(str(exc))
        else:
            return
    raise RawBindingError("; ".join(failures))


def verify_dataset_record(
    root: str | Path,
    record: Mapping[str, object],
    *,
    dataset_id: str,
    normalized_transition_count: int,
) -> None:
    """Verify an exact dataset identity, row count, and every raw binding."""
    if record.get("dataset_id") != dataset_id:
        raise RawBindingError(f"expected exact dataset id {dataset_id}")
    if record.get("selection_complete") is not True:
        raise RawBindingError(f"{dataset_id} is not asserted selection-complete")
    try:
        declared_count = int(record.get("normalized_transition_count", -1))
    except (TypeError, ValueError) as exc:
        raise RawBindingError(f"{dataset_id} has an invalid normalized row count") from exc
    if normalized_transition_count <= 0 or declared_count != normalized_transition_count:
        raise RawBindingError(
            f"{dataset_id} normalized count {normalized_transition_count} "
            f"does not match declaration {declared_count}"
        )
    artifacts = record.get("raw_artifacts")
    if not isinstance(artifacts, list) or not artifacts:
        raise RawBindingError(f"{dataset_id} has no raw-artifact bindings")
    for artifact in artifacts:
        if not isinstance(artifact, dict):
            raise RawBindingError(f"{dataset_id} has a malformed raw-artifact binding")
        verify_raw_binding(root, artifact)
