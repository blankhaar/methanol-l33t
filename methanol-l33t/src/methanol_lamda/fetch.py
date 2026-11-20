"""Immutable, allow-listed source capture from ``config/sources.json``."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import tempfile
import urllib.parse
import urllib.request


class FetchError(RuntimeError):
    pass


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _capture_status(prefix: str, actual: str, expected: str | None) -> str:
    """Describe capture integrity without treating an unpinned file as verified."""

    if not expected:
        return f"{prefix}_unpinned"
    return f"{prefix}_verified" if expected.lower() == actual else "checksum_mismatch"


def _require_allowed_https(url: str, allowed_hosts: set[str], source_id: str) -> None:
    parsed = urllib.parse.urlparse(url)
    host = parsed.hostname.lower() if parsed.hostname else None
    if parsed.scheme.lower() != "https" or host not in allowed_hosts:
        raise FetchError(
            f"source {source_id} resolved outside the HTTPS allow list: {url}"
        )


def _response_metadata(response: object, requested_url: str) -> dict[str, object]:
    geturl = getattr(response, "geturl", None)
    final_url = geturl() if callable(geturl) else requested_url
    metadata: dict[str, object] = {
        "final_url": final_url,
        "retrieved_utc": _utc_now(),
    }
    status = getattr(response, "status", None)
    if status is None:
        getcode = getattr(response, "getcode", None)
        status = getcode() if callable(getcode) else None
    if status is not None:
        metadata["http_status"] = status
    headers = getattr(response, "headers", None)
    content_type = headers.get("Content-Type") if headers is not None else None
    if content_type:
        metadata["content_type"] = content_type
    return metadata


def _stage_response(response: object, target: Path) -> tuple[Path, str, int]:
    """Stream a response into a durable sibling temporary file."""

    digest = hashlib.sha256()
    byte_count = 0
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=target.parent,
            prefix=f".{target.name}.",
            suffix=".part",
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            while True:
                block = response.read(1024 * 1024)  # type: ignore[attr-defined]
                if not block:
                    break
                handle.write(block)
                digest.update(block)
                byte_count += len(block)
            handle.flush()
            os.fsync(handle.fileno())
        return temporary_path, digest.hexdigest(), byte_count
    except Exception:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
        raise


def fetch_sources(root: str | Path, *, offline: bool = False) -> dict[str, object]:
    root = Path(root)
    config_path = root / "config" / "sources.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    ledger: dict[str, object] = {
        "generated_utc": _utc_now(), "sources": []
    }
    raw_root = root / "data" / "raw"
    raw_root.mkdir(parents=True, exist_ok=True)
    allowed_hosts = {str(host).lower() for host in config.get("allowed_hosts", [])}
    for item in config["sources"]:
        record: dict[str, object] = {
            "id": item["id"], "configured_status": item["status"], "artifacts": []
        }
        artifacts: list[dict[str, object]] = record["artifacts"]  # type: ignore[assignment]
        for local in item.get("local_paths", []):
            target = root / local["path"]
            if not target.exists():
                artifacts.append({"path": local["path"], "status": "missing_local_cache"})
                continue
            actual = sha256_path(target)
            expected = local.get("sha256")
            status = _capture_status("present", actual, expected)
            artifacts.append(
                {
                    "path": local["path"],
                    "status": status,
                    "sha256": actual,
                    "expected_sha256": expected,
                    "bytes": target.stat().st_size,
                }
            )
        for download in item.get("downloads", []):
            target = raw_root / download["target"]
            artifact: dict[str, object] = {"path": str(target.relative_to(root)), "url": download["url"]}
            if target.exists():
                actual = sha256_path(target)
                expected = download.get("sha256")
                artifact.update(
                    {
                        "status": _capture_status("present", actual, expected),
                        "sha256": actual,
                        "expected_sha256": expected,
                        "bytes": target.stat().st_size,
                    }
                )
                artifacts.append(artifact)
                continue
            if offline:
                artifact["status"] = "missing_offline"
                artifacts.append(artifact)
                continue
            _require_allowed_https(download["url"], allowed_hosts, item["id"])
            temporary_path: Path | None = None
            try:
                request = urllib.request.Request(download["url"], headers={"User-Agent": "methanol-lamda-builder/0.1"})
                with urllib.request.urlopen(request, timeout=60) as response:
                    metadata = _response_metadata(response, download["url"])
                    artifact.update(metadata)
                    try:
                        _require_allowed_https(str(metadata["final_url"]), allowed_hosts, item["id"])
                    except FetchError as exc:
                        artifact.update({"status": "redirect_rejected", "error": str(exc)})
                    else:
                        target.parent.mkdir(parents=True, exist_ok=True)
                        temporary_path, actual, byte_count = _stage_response(response, target)
                        expected = download.get("sha256")
                        status = _capture_status("downloaded", actual, expected)
                        artifact.update(
                            {
                                "status": status,
                                "sha256": actual,
                                "expected_sha256": expected,
                                "bytes": byte_count,
                            }
                        )
                        if status != "checksum_mismatch":
                            temporary_path.chmod(0o444)
                            os.replace(temporary_path, target)
                            temporary_path = None
            except Exception as exc:
                artifact.update({"status": "fetch_failed", "error": str(exc)})
            finally:
                if temporary_path is not None:
                    temporary_path.unlink(missing_ok=True)
            artifacts.append(artifact)
        if not artifacts:
            record["capture_status"] = "metadata_or_manual_only"
        elif any(artifact["status"] == "checksum_mismatch" for artifact in artifacts):
            record["capture_status"] = "checksum_mismatch"
        elif any(artifact["status"] == "redirect_rejected" for artifact in artifacts):
            record["capture_status"] = "redirect_rejected"
        elif all(
            artifact["status"] in {"present_verified", "downloaded_verified"}
            for artifact in artifacts
        ):
            record["capture_status"] = "all_configured_artifacts_verified"
        elif all(
            str(artifact["status"]).startswith(("present_", "downloaded_"))
            for artifact in artifacts
        ):
            record["capture_status"] = "all_configured_artifacts_present_with_unpinned"
        else:
            record["capture_status"] = "configured_artifacts_incomplete"
        ledger["sources"].append(record)  # type: ignore[union-attr]
    ledger_path = raw_root / "source-ledger.json"
    ledger_path.write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    mismatches = [
        f"{record['id']}:{artifact['path']}"
        for record in ledger["sources"]  # type: ignore[index]
        for artifact in record["artifacts"]
        if artifact["status"] == "checksum_mismatch"
    ]
    if mismatches:
        raise FetchError("raw-source checksum mismatch: " + ", ".join(mismatches))
    rejected_redirects = [
        f"{record['id']}:{artifact['path']}"
        for record in ledger["sources"]  # type: ignore[index]
        for artifact in record["artifacts"]
        if artifact["status"] == "redirect_rejected"
    ]
    if rejected_redirects:
        raise FetchError(
            "raw-source redirect outside allow list: " + ", ".join(rejected_redirects)
        )
    return ledger
