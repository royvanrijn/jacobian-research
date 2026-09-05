"""One process-group supervisor for bounded research workers (Linux/POSIX).

A separate pipe watchdog also kills the group if the supervisor is SIGKILLed.
Limits apply before exec; group RSS includes descendants. No shell is used.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from contextlib import nullcontext
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import resource
import signal
import subprocess
import sys
import tempfile
import time

from .store import checkpoint


class WorkerFailure(RuntimeError):
    def __init__(self, record):
        self.record = record
        super().__init__(f"research worker stopped: {record['outcome']}")


def preserve_previous(path):
    """Keep overwritten discovery logs/results by content hash."""
    path = Path(path)
    if not path.exists():
        return None
    value = sha256(path.read_bytes()).hexdigest()
    retained = path.parent / "previous" / value / path.name
    retained.parent.mkdir(parents=True, exist_ok=True)
    if retained.exists():
        if sha256(retained.read_bytes()).hexdigest() != value:
            raise ValueError("corrupt retained worker output")
        path.unlink()
    else:
        path.replace(retained)
    return retained


def capture(command, *, limits, log_path=None, cwd=None, input_text=None, check=True, env=None, separate_stderr=False):
    """Bounded capture for version probes and small CAS adapters.

    Raw logs survive both success and exceptions. Lifetime handling stays in
    run(); callers receive familiar CompletedProcess/TimeoutExpired objects.
    """
    if log_path is None:
        from .store import default_store, digest
        directory = default_store().root.parent / "runtime-workers" / digest(
            {"argv": list(map(str, command)), "invocation_ns": time.time_ns()})
        log_path = directory / "worker.log"
    stderr_path = Path(log_path).with_suffix(".stderr.log") if separate_stderr else None
    record = run(command, limits=limits, log_path=log_path, cwd=cwd, env=env, input_text=input_text, stderr_path=stderr_path,
                 checkpoint_path=Path(log_path).with_suffix(".supervisor.json"))
    output = Path(log_path).read_text(errors="replace")
    errors = stderr_path.read_text(errors="replace") if stderr_path else ""
    if record["outcome"] == "strict_wall_timeout":
        raise subprocess.TimeoutExpired(command, limits.wall_seconds, output=output, stderr=errors)
    if record["outcome"] not in ("completed", "backend_failure"):
        raise WorkerFailure(record)
    result = subprocess.CompletedProcess(command, record["returncode"], output, errors)
    result.supervision = record
    if check:
        result.check_returncode()
    return result


def log_summary(path, *, tail_bytes=100_000):
    """Bound memory while hashing a potentially large arithmetic log."""
    path=Path(path); h=sha256(); size=0
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(1<<20),b''):
            h.update(chunk);size+=len(chunk)
        stream.seek(max(0,size-tail_bytes));tail=stream.read().decode(errors='replace')
    return {"byte_count":size,"sha256":h.hexdigest(),"tail":tail,"tail_truncated":size>tail_bytes,"path":str(path)}


def capture_record(command, *, limits, log_path=None, cwd=None, input_text=None, env=None):
    """Structured capture including failures, with separate retained streams."""
    from .store import default_store, digest, atomic_write
    if log_path is None:
        directory = default_store().root.parent/"runtime-workers"/digest(
            {"argv":list(map(str,command)),"invocation_ns":time.time_ns()})
        log_path = directory/"worker.log"
    log_path = Path(log_path)
    if input_text is not None:
        source = log_path.with_suffix('.stdin')
        preserve_previous(source)
        atomic_write(source,input_text.encode())
    record = run(command,limits=limits,log_path=log_path,cwd=cwd,env=env,input_text=input_text,
        stderr_path=log_path.with_suffix('.stderr.log'),checkpoint_path=log_path.with_suffix('.supervisor.json'))
    record.update({"stdout":log_path.read_text(errors='replace'),
                   "stderr":log_path.with_suffix('.stderr.log').read_text(errors='replace')})
    return record


def captured_run(args, *, input=None, capture_output=False, text=False, encoding=None,
                 errors=None, universal_newlines=None, stdout=None, stderr=None,
                 check=False, timeout=None, cwd=None, env=None, shell=False,
                 limits=None, log_path=None):
    """Bounded subprocess.run adapter for legacy CAS bridges.

    Default bounds are 300 seconds / 1 GiB and can be explicitly configured
    with EC_WORKER_WALL_SECONDS / EC_WORKER_RSS_BYTES. Caller timeouts override
    the wall default. Separate stdout/stderr and native return types survive;
    ownership, limits and retry evidence live exclusively in this module.
    """
    if shell or isinstance(args, str):
        raise ValueError("research workers require an explicit argv, without a shell")
    if capture_output:
        if stdout is not None or stderr is not None:
            raise ValueError("capture_output conflicts with stdout/stderr")
        stdout = stderr = subprocess.PIPE
    if input is not None and not isinstance(input, (str, bytes)):
        raise TypeError("worker input must be text or bytes")
    codec = encoding or 'utf-8'
    input_text = input.decode(codec, errors or 'strict') if isinstance(input,bytes) else input
    if limits is None:
        limits = Limits(float(timeout if timeout is not None else os.environ.get('EC_WORKER_WALL_SECONDS',300)),
            int(os.environ.get('EC_WORKER_RSS_BYTES',1_073_741_824)),
            pari_stack_bytes=int(os.environ.get('EC_PARI_STACK_BYTES',256_000_000)))
    elif timeout is not None:
        from dataclasses import replace
        limits = replace(limits,wall_seconds=min(limits.wall_seconds,timeout))
    result = capture(args,limits=limits,log_path=log_path,cwd=cwd,env=env,input_text=input_text,
                     check=False,separate_stderr=stderr!=subprocess.STDOUT)
    want_text = text or universal_newlines or encoding is not None or errors is not None
    def deliver(value, target, fallback):
        converted = value if want_text else value.encode(codec,errors or 'strict')
        if target == subprocess.PIPE:
            return converted
        if target == subprocess.DEVNULL:
            return None
        if target is None:
            fallback.write(value);fallback.flush()
        elif hasattr(target,'write'):
            try:target.write(converted)
            except TypeError:target.write(value if isinstance(converted,bytes) else value.encode(codec))
            target.flush()
        elif isinstance(target,int) and target>=0:
            os.write(target,value.encode(codec,errors or 'strict'))
        elif target != subprocess.STDOUT:
            raise ValueError("unsupported worker stream target")
        return None
    output = deliver(result.stdout,stdout,sys.stdout)
    error = None if stderr == subprocess.STDOUT else deliver(result.stderr,stderr,sys.stderr)
    completed = subprocess.CompletedProcess(args,result.returncode,output,error)
    completed.supervision = result.supervision
    if check:completed.check_returncode()
    return completed


@dataclass(frozen=True)
class Limits:
    wall_seconds: float
    rss_bytes: int
    address_space_bytes: int | None = None
    pari_stack_bytes: int = 256_000_000
    terminate_grace_seconds: float = 1.0
    poll_seconds: float = 0.05
    cgroup: str | None = None

    def __post_init__(self):
        for name in ("wall_seconds", "rss_bytes", "pari_stack_bytes", "terminate_grace_seconds", "poll_seconds"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be positive and finite")
        if self.address_space_bytes is not None and (type(self.address_space_bytes) is not int or self.address_space_bytes <= 0):
            raise ValueError("invalid address-space limit")


def _start_token(pid):
    try:
        # comm can contain spaces or parentheses.
        return Path(f"/proc/{pid}/stat").read_text().rsplit(")", 1)[1].split()[19]
    except (FileNotFoundError, ProcessLookupError):
        return None


def _process_scope(pgid):
    """Owned process group and its descendants, including nested supervisors."""
    processes = {}
    for entry in Path("/proc").iterdir():
        if not entry.name.isdecimal():
            continue
        try:
            fields = (entry / "stat").read_text().rsplit(")", 1)[1].split()
            processes[int(entry.name)] = (int(fields[1]), int(fields[2]),
                int(fields[21])*os.sysconf("SC_PAGE_SIZE"), fields[19])
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            pass
    owned = {pid for pid, row in processes.items() if row[1] == pgid}
    owned.add(pgid)
    while True:
        children = {pid for pid, row in processes.items() if row[0] in owned}
        if children <= owned:
            break
        owned.update(children)
    return {pid: processes[pid] for pid in owned if pid in processes}


def _group_rss(pgid):
    return sum(row[2] for row in _process_scope(pgid).values())


def _signal_group(pgid, sig):
    try:
        os.killpg(pgid, sig)
    except ProcessLookupError:
        pass


def _stop_group(pgid, grace):
    # Kill descendants even when the original group leader has already exited.
    descendants = _process_scope(pgid)
    for pid, row in descendants.items():
        if row[1] != pgid and _start_token(pid) == row[3]:
            try:
                os.kill(pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
    _signal_group(pgid, signal.SIGTERM)
    deadline = time.monotonic() + grace
    while time.monotonic() < deadline:
        try:
            os.killpg(pgid, 0)
        except ProcessLookupError:
            if not any(_start_token(pid) == row[3] for pid, row in descendants.items() if row[1] != pgid):
                return
        time.sleep(min(0.025, max(0, deadline - time.monotonic())))
    _signal_group(pgid, signal.SIGKILL)
    for pid, row in descendants.items():
        if row[1] != pgid and _start_token(pid) == row[3]:
            try:
                os.kill(pid, signal.SIGKILL)
            except ProcessLookupError:
                pass


def _watchdog(read_fd, pgid, token, grace):
    # Parent holds the only write end. Explicit cleanup sends a byte; EOF
    # without that byte means interruption/crash and triggers orphan cleanup.
    try:
        disarmed = os.read(read_fd, 1)
        if not disarmed:
            current = _start_token(pgid)
            if current is None or current == token:
                _stop_group(pgid, grace)
    finally:
        os.close(read_fd)


def _launch(config_path, ready_fd):
    # Do not start the CAS before its watchdog exists. If the supervisor dies
    # in that startup window, EOF makes this launcher exit instead of orphaning
    # an unguarded worker.
    try:
        if os.read(ready_fd,1)!=b"1":raise SystemExit("supervisor disappeared before worker admission")
    finally:
        os.close(ready_fd)
    config = json.loads(Path(config_path).read_text())
    limits = Limits(**config["limits"])
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    if limits.address_space_bytes is not None:
        resource.setrlimit(resource.RLIMIT_AS, (limits.address_space_bytes, limits.address_space_bytes))
    if limits.cgroup:
        cgroup = Path(limits.cgroup)
        # The caller supplies a delegated, per-job cgroup. Refuse a hierarchy root.
        if cgroup.resolve() == Path("/sys/fs/cgroup"):
            raise ValueError("use a dedicated delegated cgroup")
        (cgroup / "memory.max").write_text(str(limits.rss_bytes))
        (cgroup / "cgroup.procs").write_text(str(os.getpid()))
    os.execvpe(config["command"][0], config["command"], os.environ)


def run(command, *, limits: Limits, log_path: Path, result_path=None,
        checkpoint_path=None, cwd=None, env=None, input_text=None, stderr_path=None):
    """Run one worker and return a structured outcome, preserving raw logs.

    A zero exit does not attest a mathematical result. If result_path is given,
    a fresh nonempty JSON result is required for the outcome ``completed``.
    Existing result files must be explicitly handled by the caller before run.
    """
    command = list(map(str, command))
    if not command or not command[0]:
        raise ValueError("an explicit argument vector is required")
    log_path = Path(log_path).resolve()
    result_path = Path(result_path).resolve() if result_path else None
    if result_path and result_path.exists():
        raise FileExistsError("preserve or explicitly retire the existing worker result")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    preserve_previous(log_path)
    if stderr_path:
        stderr_path = Path(stderr_path).resolve()
        if stderr_path == log_path:
            raise ValueError("separate stderr needs a different log path")
        stderr_path.parent.mkdir(parents=True, exist_ok=True)
        preserve_previous(stderr_path)
    if checkpoint_path:
        preserve_previous(checkpoint_path)
    if result_path:
        result_path.parent.mkdir(parents=True, exist_ok=True)
    child_env = dict(os.environ)
    child_env.update(env or {})
    cas = str(Path(__file__).resolve().parents[1])
    child_env["PYTHONPATH"] = cas + os.pathsep + child_env.get("PYTHONPATH", "")
    child_env["EC_PARI_STACK_BYTES"] = str(limits.pari_stack_bytes)
    started = time.monotonic()
    peak_rss = 0
    outcome = "running"
    process = watcher = None
    write_fd = read_fd = None
    ready_read = ready_write = None
    record = {"command": command, "limits": asdict(limits), "outcome": outcome}
    with tempfile.TemporaryDirectory(prefix="ec-worker-") as directory:
        config = Path(directory) / "launch.json"
        checkpoint(config, {"command": command, "limits": asdict(limits)})
        with log_path.open("w") as log, (stderr_path.open("w") if stderr_path else nullcontext(subprocess.STDOUT)) as error_log, tempfile.TemporaryFile(mode="w+") as stdin:
            if input_text is not None:
                stdin.write(input_text)
                stdin.seek(0)
            try:
                ready_read,ready_write=os.pipe()
                process = subprocess.Popen(
                    [sys.executable, "-m", "research_runtime.supervisor", "exec", str(config),str(ready_read)],
                    stdin=stdin, stdout=log, stderr=error_log,
                    start_new_session=True, cwd=cwd, env=child_env,pass_fds=(ready_read,),
                )
                os.close(ready_read);ready_read=None
                read_fd, write_fd = os.pipe()
                token = _start_token(process.pid)
                watcher = subprocess.Popen(
                    [sys.executable, "-m", "research_runtime.supervisor", "watch",
                     str(read_fd), str(process.pid), str(token), str(limits.terminate_grace_seconds)],
                    pass_fds=(read_fd,), stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
                    start_new_session=True, env=child_env,
                )
                os.close(read_fd)
                read_fd = None
                record.update({"pid": process.pid, "start_token": token})
                if checkpoint_path:
                    checkpoint(Path(checkpoint_path), record)
                os.write(ready_write,b"1");os.close(ready_write);ready_write=None
                while process.poll() is None:
                    peak_rss = max(peak_rss, _group_rss(process.pid))
                    if peak_rss > limits.rss_bytes:
                        outcome = "strict_rss_limit"
                        break
                    if time.monotonic() - started >= limits.wall_seconds:
                        outcome = "strict_wall_timeout"
                        break
                    time.sleep(limits.poll_seconds)
                if outcome == "running":
                    outcome = "completed" if process.returncode == 0 else "backend_failure"
            except BaseException:
                outcome = "interrupted"
                raise
            finally:
                if ready_read is not None:os.close(ready_read)
                if ready_write is not None:os.close(ready_write)
                if process is not None:
                    _stop_group(process.pid, limits.terminate_grace_seconds)
                    process.wait()
                if write_fd is not None:
                    try:
                        os.write(write_fd, b"1")
                    except BrokenPipeError:
                        pass
                    os.close(write_fd)
                if read_fd is not None:
                    os.close(read_fd)
                if watcher is not None:
                    watcher.wait(timeout=limits.terminate_grace_seconds + 5)
                record.update({"outcome": outcome, "returncode": process.returncode if process else None,
                               "wall_seconds": time.monotonic() - started,
                               "peak_observed_rss_bytes": peak_rss})
                if checkpoint_path:
                    checkpoint(Path(checkpoint_path), record)
    if outcome == "completed" and result_path:
        try:
            result = json.loads(result_path.read_text())
            if not isinstance(result, (dict, list)) or not result:
                raise ValueError("empty worker result")
        except (FileNotFoundError, ValueError):
            outcome = "missing_or_invalid_result"
    record.update({"outcome": outcome, "timeout_seconds": limits.wall_seconds,
                   "failure_reason":None if outcome=="completed" else outcome,
                   "termination_signal":-record["returncode"] if record["returncode"] is not None and record["returncode"]<0 else None,
                   "rss_limit_bytes": limits.rss_bytes, "log": str(log_path),
                   "log_sha256": sha256(log_path.read_bytes()).hexdigest(),
                   "worker_result": str(result_path) if result_path and result_path.exists() else None,
                   "worker_result_sha256": sha256(result_path.read_bytes()).hexdigest() if result_path and result_path.exists() else None})
    if stderr_path:
        record.update({"stderr_log":str(stderr_path),"stderr_sha256":sha256(stderr_path.read_bytes()).hexdigest()})
    if checkpoint_path:
        checkpoint(Path(checkpoint_path), record)
    return record


def supervise_source(sage_python, worker_source, payload, result_path, log_path, *, timeout, rss_limit_bytes):
    """Compatibility adapter for existing retained Python worker protocols."""
    result_path = Path(result_path).resolve()
    log_path = Path(log_path).resolve()
    # Inputs and source remain reviewable after completion or interruption.
    input_path = log_path.with_suffix(log_path.suffix + ".input.json")
    worker_path = log_path.with_suffix(log_path.suffix + ".worker.py")
    preserve_previous(input_path)
    preserve_previous(worker_path)
    checkpoint(input_path, payload)
    from .store import atomic_write
    atomic_write(worker_path, worker_source.replace("INPUT_PATH", repr(str(input_path)))
                 .replace("OUTPUT_PATH", repr(str(result_path))).encode())
    return run([sage_python, str(worker_path)], limits=Limits(timeout, rss_limit_bytes,
               pari_stack_bytes=int(payload.get("pari_stack_bytes", payload.get("stack_bytes", 256_000_000)))),
               log_path=log_path, result_path=result_path, checkpoint_path=log_path.with_suffix(".supervisor.json"))


if __name__ == "__main__":
    if sys.argv[1] == "exec":
        _launch(sys.argv[2],int(sys.argv[3]))
    elif sys.argv[1] == "watch":
        _watchdog(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], float(sys.argv[5]))
    else:
        raise SystemExit("unknown supervisor mode")
