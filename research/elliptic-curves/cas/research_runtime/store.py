"""Content-addressed exact facts with locked discovery and atomic publication."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import asdict, is_dataclass
from fractions import Fraction
from hashlib import sha256
import fcntl
import json
import os
from pathlib import Path
import tempfile
from typing import Callable


def canonical(value):
    if is_dataclass(value):
        value = asdict(value)
    if isinstance(value, Fraction):
        return {"rational": [value.numerator, value.denominator]}
    if isinstance(value, dict):
        if any(not isinstance(key, str) for key in value):
            raise TypeError("cache mappings require string keys")
        return {key: canonical(item) for key, item in sorted(value.items())}
    if isinstance(value, (tuple, list)):
        return [canonical(item) for item in value]
    if value is None or type(value) in (str, int, bool):
        return value
    raise TypeError(f"not exact JSON data: {type(value).__name__}")


def encoded(value) -> bytes:
    return json.dumps(canonical(value), sort_keys=True, separators=(",", ":")).encode()


def digest(value) -> str:
    return sha256(encoded(value)).hexdigest()


def atomic_write(path: Path, data: bytes) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        directory_fd = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        Path(temporary).unlink(missing_ok=True)


def checkpoint(path: Path, value) -> None:
    # Timing measurements may be floating point; proof/cache identity may not.
    atomic_write(path, (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode())


class FactStore:
    """An immutable result per (namespace, exact inputs, algorithm version).

    A missing fact stays missing during replay. Discovery is always explicit.
    A lock spans the lookup and builder so competing workers do not duplicate
    arithmetic. Interrupted builders publish nothing. Corruption is an error,
    never a reason to silently regenerate a proof.
    """

    def __init__(self, root: Path):
        self.root = Path(root)
        self._accessed = set()

    @contextmanager
    def lock(self, key):
        path = self.root / "locks" / (digest(key) + ".lock")
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("a+b") as handle:
            fcntl.flock(handle, fcntl.LOCK_EX)
            yield

    def _key(self, namespace, inputs, version):
        if not namespace or not version:
            raise ValueError("facts require a namespace and algorithm version")
        return {"namespace": namespace, "inputs": canonical(inputs), "version": version}

    def get(self, namespace, inputs, *, version="1"):
        key = self._key(namespace, inputs, version)
        path = self.root / "index" / (digest(key) + ".json")
        if not path.exists():
            return None
        index = json.loads(path.read_bytes())
        blob_hash = index["sha256"]
        if len(blob_hash) != 64 or any(c not in "0123456789abcdef" for c in blob_hash):
            raise ValueError("invalid fact content address")
        raw = (self.root / "objects" / (blob_hash + ".json")).read_bytes()
        if sha256(raw).hexdigest() != blob_hash:
            raise ValueError("corrupt fact content")
        record = json.loads(raw)
        if record["key"] != key:
            raise ValueError("fact input/algorithm mismatch")
        self._accessed.add(blob_hash)
        return record["value"]

    def require(self, namespace, inputs, *, version="1"):
        result = self.get(namespace, inputs, version=version)
        if result is None:
            raise FileNotFoundError(f"missing retained {namespace} fact; run discovery explicitly")
        return result

    def discover(self, namespace, inputs, build: Callable, *, version="1", validate=None):
        key = self._key(namespace, inputs, version)
        with self.lock(key):
            result = self.get(namespace, inputs, version=version)
            if result is None:
                result = canonical(build())
                if result is None:
                    raise ValueError("a fact cannot have a null payload")
                if validate is not None:
                    validate(result)
                raw = encoded({"key": key, "value": result})
                blob_hash = sha256(raw).hexdigest()
                atomic_write(self.root / "objects" / (blob_hash + ".json"), raw)
                atomic_write(self.root / "index" / (digest(key) + ".json"), encoded({"sha256": blob_hash}))
                self._accessed.add(blob_hash)
            elif validate is not None:
                validate(result)
            return result

    def snapshot(self):
        """Portable exact facts touched by this run, alongside its proof witnesses.

        This preserves expensive setup for no-cache replay. Content hashes
        attest integrity only; mathematical replayers must still check witnesses.
        """
        facts = []
        for blob_hash in sorted(self._accessed):
            raw = (self.root / "objects" / (blob_hash + ".json")).read_bytes()
            if sha256(raw).hexdigest() != blob_hash:
                raise ValueError("corrupt snapshot fact")
            facts.append({"sha256": blob_hash, "record": json.loads(raw)})
        return {"schema": "elliptic-curves.arithmetic-facts.v1", "facts": facts}

    def import_snapshot(self, snapshot):
        """Restore retained facts without invoking any arithmetic builder."""
        if snapshot.get("schema") != "elliptic-curves.arithmetic-facts.v1":
            raise ValueError("unknown arithmetic snapshot schema")
        for item in snapshot["facts"]:
            record, blob_hash = item["record"], item["sha256"]
            raw = encoded(record)
            if sha256(raw).hexdigest() != blob_hash:
                raise ValueError("snapshot content hash mismatch")
            key = record["key"]
            if key != self._key(key["namespace"], key["inputs"], key["version"]):
                raise ValueError("invalid snapshot fact key")
            with self.lock(key):
                existing = self.get(key["namespace"], key["inputs"], version=key["version"])
                if existing is not None and existing != record["value"]:
                    raise ValueError("snapshot conflicts with an existing immutable fact")
                if existing is None:
                    atomic_write(self.root / "objects" / (blob_hash + ".json"), raw)
                    atomic_write(self.root / "index" / (digest(key) + ".json"), encoded({"sha256": blob_hash}))
                self._accessed.add(blob_hash)


def default_store() -> FactStore:
    path = os.environ.get("EC_ARITHMETIC_CACHE")
    return FactStore(Path(path) if path else Path(__file__).resolve().parents[3] / "artifacts" / "local" / "arithmetic-cache")


class FiniteFieldFacts:
    """Shared traces, quotient bases, fibres and Kummer maps at exact primes.

    Isomorphic presentations are shared only after an explicit transport; the
    j-invariant alone is deliberately not an identity. Extension degree and
    coordinate labels distinguish different local quotient spaces.
    """

    def __init__(self, store: FactStore):
        self.store = store

    def query(self, model, prime, kind, *, degree=1, labels=(), build=None, version="1"):
        if type(prime) is not int or prime < 2 or type(degree) is not int or degree < 1:
            raise ValueError("invalid finite-field parameters")
        inputs = {"model": model, "prime": prime, "degree": degree, "labels": labels}
        if build is None:
            return self.store.require("finite-field/" + kind, inputs, version=version)
        return self.store.discover("finite-field/" + kind, inputs, build, version=version)
