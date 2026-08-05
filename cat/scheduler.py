"""
cat/scheduler.py - Cluster scheduler abstraction.
"""

from __future__ import annotations

import logging
import os
import re
import shlex
import shutil
import stat
import subprocess
import time
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence

logger = logging.getLogger(__name__)


def _submit_with_retries(
    cmd: Sequence[str],
    *,
    attempts: int = 4,
    base_delay_s: float = 5.0,
) -> subprocess.CompletedProcess:
    """Run a job-submission command, retrying transient scheduler failures.

    Busy SLURM/SGE controllers intermittently reject submissions (socket
    timeouts, "Batch job submission failed", momentary daemon unavailability).
    A single such hiccup should not abort a rule, so retry with exponential
    backoff before giving up. Deterministic errors (bad script, bad flags) will
    fail identically on every attempt and still surface after the last try.
    """
    last_exc: Optional[Exception] = None
    for attempt in range(1, attempts + 1):
        try:
            return subprocess.run(cmd, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as exc:
            last_exc = exc
            if attempt == attempts:
                break
            delay = base_delay_s * (2 ** (attempt - 1))
            logger.warning(
                f"submission '{cmd[0]}' failed (attempt {attempt}/{attempts}, "
                f"rc={exc.returncode}): {(exc.stderr or '').strip()[:200]}; "
                f"retrying in {delay:.0f}s"
            )
            time.sleep(delay)
    assert last_exc is not None
    raise last_exc


# ──────────────────────────────────────────────────────────────────────────────
# Result types
# ──────────────────────────────────────────────────────────────────────────────


@dataclass
class JobResult:
    """Outcome of a Scheduler.wait() call."""

    ok: bool
    completed: int = 0
    failed: int = 0
    total: int = 0
    detail: str = ""

    def __bool__(self) -> bool:
        return self.ok


# ──────────────────────────────────────────────────────────────────────────────
# Scheduler base class
# ──────────────────────────────────────────────────────────────────────────────


class Scheduler(ABC):
    """Abstract cluster scheduler.

    Subclasses implement header(), submit(), wait(), and cancel(). The base
    class provides the shared script_preamble() helper and the write_script()
    utility (which sets the executable bit and returns the path).
    """

    name: str = "abstract"

    # ── Script generation ─────────────────────────────────────────────────

    @abstractmethod
    def header(
        self,
        *,
        job_name: str,
        cpus: int,
        mem: str,
        walltime: str,
        log_out: str,
        log_err: str,
        partition: Optional[str] = None,
        queue: Optional[str] = None,
        exclude: Optional[str] = None,
        array: Optional[tuple[int, int]] = None,
        max_concurrent: Optional[int] = None,
        dependency: Optional[str] = None,
        module_load: Optional[str] = None,
        extra_directives: Optional[Sequence[str]] = None,
    ) -> str:
        """Build the shebang + scheduler-directive block for a job script.

        Parameters
        ----------
        job_name : str
            Human-readable name shown in qstat/squeue.
        cpus : int
            Cores per task (translates to SLURM --cpus-per-task or SGE -pe smp N).
        mem : str
            Memory request with size suffix (e.g. "128G"). The unit is preserved
            verbatim for SLURM and re-emitted under the configured SGE resource
            flag (h_vmem by default).
        walltime : str
            Wallclock limit "HH:MM:SS".
        log_out, log_err : str
            Stdout/stderr log paths. May contain scheduler placeholders such as
            "%A_%a" (SLURM) or "$JOB_ID.$TASK_ID" (SGE); callers that need
            portability should let the scheduler choose the placeholder
            (e.g. by calling array_log_paths()).
        partition : str, optional
            SLURM partition. Ignored by SGE / Local.
        queue : str, optional
            SGE queue. Ignored by SLURM / Local.
        exclude : str, optional
            Hostname exclude expression. SLURM accepts a comma list
            ("host1,host2"). SGE accepts its native form
            ("!host1&!host2"); this module converts comma lists to the
            SGE form automatically.
        array : (start, end), optional
            Task ID range. Always interpret as **1-based, inclusive**.
        max_concurrent : int, optional
            Cap on simultaneously running array tasks.
        dependency : str, optional
            Backend-specific dependency expression. For SLURM, anything
            valid in --dependency=; for SGE, a job ID or comma list (the
            scheduler emits -hold_jid). Use depends_on_job_id() for the
            common "after this finishes ok" case.
        module_load : str, optional
            ``module load NAME`` line inserted after the directives. Empty/None
            skips the line.
        extra_directives : sequence of str, optional
            Raw extra directive lines, already in the backend's native syntax.
        """

    def array_log_paths(self, directory: str | os.PathLike, stem: str) -> tuple[str, str]:
        """Return (stdout, stderr) paths with backend-specific array placeholders.

        SLURM expands ``%A`` / ``%a``; SGE expands ``$JOB_ID`` / ``$TASK_ID``.
        Passing SLURM tokens to SGE leaves a literal ``%A_%a`` filename and
        makes failures nearly undiagnosable.
        """
        d = str(directory).rstrip("/")
        if self.name == "sge":
            return (
                f"{d}/{stem}.$JOB_ID.$TASK_ID.out",
                f"{d}/{stem}.$JOB_ID.$TASK_ID.err",
            )
        return (f"{d}/{stem}_%A_%a.out", f"{d}/{stem}_%A_%a.err")

    def job_log_paths(self, directory: str | os.PathLike, stem: str) -> tuple[str, str]:
        """Return (stdout, stderr) paths for a non-array job."""
        d = str(directory).rstrip("/")
        if self.name == "sge":
            return (f"{d}/{stem}.$JOB_ID.out", f"{d}/{stem}.$JOB_ID.err")
        return (f"{d}/{stem}_%j.out", f"{d}/{stem}_%j.err")

    def script_preamble(
        self,
        *,
        conda_env: Optional[str] = None,
        set_strict: bool = True,
        extra_env: Optional[dict[str, str]] = None,
    ) -> str:
        """Boilerplate body lines: conda activation + ``set -euo pipefail``.

        Replaces the conda/set-strict snippets that used to live verbatim in
        align_transcripts_slurm.py, parent_gene_assignment_slurm.py, and
        consensus_parallel.py.
        """
        lines: list[str] = []
        if conda_env:
            # Activate before strict mode to avoid unbound-variable errors from
            # conda's init scripts. Prefer CONDA_EXE (propagated with -V) so
            # non-~/miniconda3 installs work on SGE compute nodes.
            lines.append(
                'if [ -n "${CONDA_EXE:-}" ] && '
                '[ -f "$(dirname "$(dirname "$CONDA_EXE")")/etc/profile.d/conda.sh" ]; then\n'
                '  # shellcheck source=/dev/null\n'
                '  source "$(dirname "$(dirname "$CONDA_EXE")")/etc/profile.d/conda.sh"\n'
                'elif [ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]; then\n'
                '  # shellcheck source=/dev/null\n'
                '  source "${HOME}/miniconda3/etc/profile.d/conda.sh"\n'
                'elif [ -f "${HOME}/mambaforge/etc/profile.d/conda.sh" ]; then\n'
                '  # shellcheck source=/dev/null\n'
                '  source "${HOME}/mambaforge/etc/profile.d/conda.sh"\n'
                'fi'
            )
            lines.append(f"conda activate {shlex.quote(conda_env)}")
            lines.append('export PATH="${CONDA_PREFIX}/bin:$PATH"')
        if set_strict:
            lines.append("set -euo pipefail")
        if extra_env:
            for k, v in extra_env.items():
                lines.append(f"export {k}={shlex.quote(v)}")
        return "\n".join(lines)

    def write_script(self, content: str, dest: str | os.PathLike) -> str:
        """Write *content* to *dest*, chmod 0755, return absolute path string."""
        path = Path(dest)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
        path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        return str(path.resolve())

    # ── Backend identity ──────────────────────────────────────────────────

    @abstractmethod
    def task_id_env(self) -> str:
        """Environment variable holding the array task ID inside the job."""

    @abstractmethod
    def array_index_base(self) -> int:
        """0 for SLURM (allows either; we treat as 1), 1 for SGE."""

    # ── Submission and waiting ────────────────────────────────────────────

    @abstractmethod
    def submit(self, script_path: str | os.PathLike) -> str:
        """Submit *script_path*; return the backend's job ID."""

    @abstractmethod
    def wait(
        self,
        job_id: str,
        *,
        num_tasks: Optional[int] = None,
        timeout_s: int = 12 * 3600,
        check_interval_s: int = 30,
        sentinel_dir: Optional[str | os.PathLike] = None,
    ) -> JobResult:
        """Block until *job_id* completes; return a JobResult.

        For array jobs (num_tasks > 1), the scheduler reports per-task progress.
        The SGE backend requires *sentinel_dir* (each task writes its exit code
        there) since qacct output is not portable. The SLURM backend ignores
        *sentinel_dir* and uses sacct.
        """

    @abstractmethod
    def cancel(self, job_id: str) -> None:
        """Cancel *job_id* (best-effort)."""

    @abstractmethod
    def job_present(self, job_id: str) -> bool:
        """True if the scheduler still has a record of *job_id* (queued/running).

        Used by call sites that want to drive their own polling loop (e.g.
        hints_db.py's per-script success/error marker inspection) instead of
        using ``wait()``. Always best-effort: transient scheduler errors
        return True so the caller keeps polling.
        """

    def verify_completed(self, job_id: str) -> JobResult:
        """After a job has drained from the queue, report whether it succeeded.

        Used by call sites that own their own poll loop (i.e. they already
        called ``job_present()`` in a loop until it returned False) and now
        need to know if the job finished cleanly.

        - SLURM: uses ``sacct`` to inspect each task's State / ExitCode.
        - SGE:   returns ``ok=True`` with a "cannot verify" detail, since
                 ``qacct`` output formats are not portable across SGE flavours.
                 Callers that need precise per-task verification on SGE should
                 wire sentinel files via :meth:`trap_sentinel` instead.
        - Local: returns ``ok=True`` (jobs are synchronous; failures would
                 have already raised by the time ``submit()`` returned).

        The default implementation here returns success; backends override
        when they have a richer post-mortem channel.
        """
        return JobResult(ok=True, detail="verify_completed not implemented for this backend")

    # ── Helpers used by callers ───────────────────────────────────────────

    def depends_on_job_id(self, job_id: str) -> str:
        """Return a dependency expression for "after *job_id* completes ok"."""
        return self._depends_on_job_id(job_id)

    @abstractmethod
    def _depends_on_job_id(self, job_id: str) -> str: ...

    def sentinel_lines(self, sentinel_dir: str | os.PathLike, task_id_var: Optional[str] = None) -> str:
        """Bash snippet that writes per-task success/failure sentinels.

        Inserted near the end of the job body. After this prelude:

            if some_command; then
                touch_sentinel ok
            else
                touch_sentinel "$?"
            fi

        Callers usually use the simpler emit_sentinel() helper instead.
        """
        var = task_id_var or self.task_id_env()
        return (
            f'_CAT_SENT_DIR={shlex.quote(str(sentinel_dir))}\n'
            f'mkdir -p "$_CAT_SENT_DIR"\n'
            f'_CAT_SENT_FILE="$_CAT_SENT_DIR/done.${{{var}:-0}}"\n'
            'touch_sentinel() {\n'
            '  echo "$1" > "$_CAT_SENT_FILE"\n'
            '}\n'
        )

    def trap_sentinel(self, sentinel_dir: str | os.PathLike, task_id_var: Optional[str] = None) -> str:
        """Bash snippet that emits a sentinel on EXIT, capturing the exit code.

        Use at the top of a job body so any path (success or failure) writes a
        sentinel file. Requires ``set -e`` to be either off or paired with an
        explicit final ``exit 0`` for the success case.
        """
        var = task_id_var or self.task_id_env()
        return (
            f'_CAT_SENT_DIR={shlex.quote(str(sentinel_dir))}\n'
            f'mkdir -p "$_CAT_SENT_DIR"\n'
            f'_CAT_SENT_FILE="$_CAT_SENT_DIR/done.${{{var}:-0}}"\n'
            'trap \'echo $? > "$_CAT_SENT_FILE"\' EXIT\n'
        )


# ──────────────────────────────────────────────────────────────────────────────
# SLURM
# ──────────────────────────────────────────────────────────────────────────────


# sacct states considered terminal failures (mirrors historical *_slurm.py).
_SLURM_FAILED_STATES = {"FAILED", "CANCELLED", "TIMEOUT", "NODE_FAIL", "OUT_OF_MEMORY", "PREEMPTED", "BOOT_FAIL"}
_SLURM_SUCCESS_STATES = {"COMPLETED"}


class SlurmScheduler(Scheduler):
    """SLURM backend. Matches the historical sbatch/sacct behaviour 1:1."""

    name = "slurm"

    def __init__(self, cfg: Optional[dict] = None):
        cfg = cfg or {}
        self.default_partition: Optional[str] = cfg.get("partition")
        self.default_exclude: str = cfg.get("exclude_nodes", "") or ""
        self.default_module_load: str = cfg.get("module_load", "") or ""

    # ── identity ─────────────────────────────────────────────────────────
    def task_id_env(self) -> str:
        return "SLURM_ARRAY_TASK_ID"

    def array_index_base(self) -> int:
        # SLURM accepts either 0- or 1-based; we standardise on 1 internally
        # so call sites are portable to SGE.
        return 1

    # ── header ────────────────────────────────────────────────────────────
    def header(
        self,
        *,
        job_name: str,
        cpus: int,
        mem: str,
        walltime: str,
        log_out: str,
        log_err: str,
        partition: Optional[str] = None,
        queue: Optional[str] = None,
        exclude: Optional[str] = None,
        array: Optional[tuple[int, int]] = None,
        max_concurrent: Optional[int] = None,
        dependency: Optional[str] = None,
        module_load: Optional[str] = None,
        extra_directives: Optional[Sequence[str]] = None,
    ) -> str:
        partition = partition or self.default_partition
        if exclude is None:
            exclude = self.default_exclude
        if module_load is None:
            module_load = self.default_module_load

        lines: list[str] = ["#!/bin/bash", f"#SBATCH --job-name={job_name}"]
        if partition:
            lines.append(f"#SBATCH --partition={partition}")
        if exclude:
            lines.append(f"#SBATCH --exclude={exclude}")
        lines.append("#SBATCH --nodes=1")
        lines.append("#SBATCH --ntasks=1")
        lines.append(f"#SBATCH --cpus-per-task={cpus}")
        lines.append(f"#SBATCH --mem={mem}")
        lines.append(f"#SBATCH --time={walltime}")
        if array:
            start, end = array
            spec = f"{start}-{end}"
            if max_concurrent and max_concurrent < (end - start + 1):
                spec += f"%{max_concurrent}"
            lines.append(f"#SBATCH --array={spec}")
        if dependency:
            lines.append(f"#SBATCH --dependency={dependency}")
        lines.append(f"#SBATCH --output={log_out}")
        lines.append(f"#SBATCH --error={log_err}")
        if extra_directives:
            for d in extra_directives:
                lines.append(d if d.startswith("#") else f"#SBATCH {d}")
        if module_load:
            lines.append("")
            lines.append(f"module load {module_load}")
        # Drop a leaked controller PYTHONPATH (wrong site-packages / py version).
        # Conda packages themselves do NOT live on PYTHONPATH — jobs must still
        # `conda activate` (see build_sbatch_header / script_preamble) so PATH
        # points at the env's python. unset alone is not enough on SGE.
        lines.append("")
        lines.append("unset PYTHONPATH")
        lines.append("")
        return "\n".join(lines) + "\n"

    def _depends_on_job_id(self, job_id: str) -> str:
        return f"afterok:{job_id}"

    # ── submit ────────────────────────────────────────────────────────────
    def submit(self, script_path: str | os.PathLike) -> str:
        script_path = str(script_path)
        logger.info(f"sbatch {script_path}")
        result = _submit_with_retries(["sbatch", script_path])
        # Output format: "Submitted batch job 12345" (or with --parsable: just "12345" or "12345;cluster")
        token = result.stdout.strip().split()[-1]
        job_id = token.split(";")[0]
        logger.info(f"Submitted SLURM job {job_id}")
        return job_id

    # ── wait ──────────────────────────────────────────────────────────────
    def wait(
        self,
        job_id: str,
        *,
        num_tasks: Optional[int] = None,
        timeout_s: int = 12 * 3600,
        check_interval_s: int = 30,
        sentinel_dir: Optional[str | os.PathLike] = None,
    ) -> JobResult:
        """Poll sacct for per-array-task state.

        Mirrors wait_for_slurm_job in cat/align_transcripts_slurm.py:
        - Filter out .batch and .extern step rows
        - Parent row (just JOBID) is fallback if no array rows seen
        - Terminal: all tasks in COMPLETED ∪ FAILED_STATES
        """
        logger.info(f"Waiting for SLURM job {job_id} (timeout: {timeout_s / 3600:.1f}h)")
        start = time.time()
        last_report = 0.0
        status_report_interval = 300.0
        # If the caller passed num_tasks at all, the job was submitted as an
        # array (``#SBATCH --array=1-N``) and sacct will emit per-task rows
        # like "JOBID_1", "JOBID_2", ... — never a bare "JOBID" parent row.
        # This is true even for N==1 (a one-task array), so we must not gate
        # the array branch on num_tasks > 1; otherwise the wait loop polls
        # forever waiting for a parent row that never appears.
        is_array = bool(num_tasks and num_tasks >= 1)

        while time.time() - start < timeout_s:
            try:
                proc = subprocess.run(
                    ["sacct", "-j", job_id, "-n", "-o", "JobID,State", "--parsable2"],
                    capture_output=True, text=True, check=True,
                )
            except subprocess.CalledProcessError as e:
                logger.warning(f"sacct error: {e}")
                time.sleep(check_interval_s)
                continue

            task_states: dict[str, str] = {}
            parent_state: Optional[str] = None
            for line in proc.stdout.strip().splitlines():
                parts = line.strip().split("|")
                if len(parts) != 2:
                    continue
                spec, state = parts
                if spec == job_id:
                    parent_state = state
                elif "_" in spec and not any(x in spec for x in (".batch", ".extern")):
                    task_states[spec] = state

            if is_array and task_states:
                completed = sum(1 for s in task_states.values() if s in _SLURM_SUCCESS_STATES)
                failed = sum(1 for s in task_states.values() if s in _SLURM_FAILED_STATES)
                running = sum(1 for s in task_states.values() if s not in _SLURM_SUCCESS_STATES and s not in _SLURM_FAILED_STATES)
                total = len(task_states)

                elapsed = time.time() - start
                if elapsed - last_report > status_report_interval:
                    logger.info(
                        f"Job {job_id} progress: {completed}/{total} done, "
                        f"{failed} failed, {running} running ({elapsed / 60:.1f} min)"
                    )
                    last_report = elapsed

                if completed + failed == total:
                    if failed > 0:
                        return JobResult(
                            ok=False, completed=completed, failed=failed, total=total,
                            detail=f"{failed}/{total} array tasks failed",
                        )
                    return JobResult(ok=True, completed=completed, total=total)
            elif parent_state is not None:
                if parent_state in _SLURM_SUCCESS_STATES:
                    return JobResult(ok=True, completed=1, total=1)
                if parent_state in _SLURM_FAILED_STATES:
                    return JobResult(ok=False, total=1, failed=1, detail=parent_state)
                # Pending / Running — keep polling

            time.sleep(check_interval_s)

        return JobResult(ok=False, detail=f"timeout after {timeout_s / 3600:.1f}h")

    def cancel(self, job_id: str) -> None:
        try:
            subprocess.run(["scancel", str(job_id)], check=False)
        except FileNotFoundError:
            logger.warning("scancel not available; cannot cancel job")

    def job_present(self, job_id: str) -> bool:
        try:
            res = subprocess.run(
                ["squeue", "-j", str(job_id), "-h", "-o", "%i"],
                capture_output=True, text=True,
            )
        except FileNotFoundError:
            return True  # squeue gone? Be conservative.
        # squeue exits 0 with empty stdout when the job is no longer queued.
        return bool(res.stdout.strip())

    def verify_completed(self, job_id: str) -> JobResult:
        """Inspect sacct for per-task State / ExitCode after the job drains.

        Mirrors the historical sacct-validation logic that lived inline in
        augustus_parallel.py / augustus_pb_parallel.py before the migration.
        Returns ok=False if any non-step task ended in a failed state OR with
        a non-zero exit code. Step rows (``.batch`` / ``.extern``) are
        ignored. DependencyNeverSatisfied is reported via the parent row's
        state being something other than COMPLETED.
        """
        try:
            res = subprocess.run(
                ["sacct", "-j", str(job_id), "--format=JobID,State,ExitCode",
                 "--noheader", "--parsable2"],
                capture_output=True, text=True, check=True,
            )
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            return JobResult(ok=True, detail=f"sacct unavailable: {e}; assuming success")

        failed: list[str] = []
        completed = 0
        total = 0
        for line in res.stdout.strip().splitlines():
            parts = line.split("|")
            if len(parts) < 3:
                continue
            spec, state, exit_code = parts[0], parts[1], parts[2]
            if ".batch" in spec or ".extern" in spec:
                continue
            total += 1
            if state in _SLURM_SUCCESS_STATES and exit_code == "0:0":
                completed += 1
            else:
                failed.append(f"{spec}(state={state},exit={exit_code})")

        if failed:
            detail = f"{len(failed)}/{total} task(s) failed: " + ", ".join(failed[:5])
            if len(failed) > 5:
                detail += f", ... +{len(failed) - 5} more"
            return JobResult(ok=False, completed=completed, failed=len(failed), total=total, detail=detail)
        return JobResult(ok=True, completed=completed, total=total)


# ──────────────────────────────────────────────────────────────────────────────
# SGE (site-agnostic)
# ──────────────────────────────────────────────────────────────────────────────


# Regex covers single-job ("Your job 12345 (...)"), array
# ("Your job-array 12345.1-N:1 (...)"), and -terse mode (just the ID).
_QSUB_JOB_ID_RE = re.compile(r"Your\s+(?:job|job-array)\s+(\d+)")


class SgeScheduler(Scheduler):
    """Site-agnostic Grid Engine backend.

    The expensive site-specific knobs (parallel environment name, memory
    resource flag, hostname-exclude syntax flavour) come from the config
    dict passed at construction.

    Polling uses filesystem sentinel files because qacct output formats vary
    between SGE forks. Each task writes its exit code into
    ``<sentinel_dir>/done.<TASK_ID>``; the wait loop checks for all N files
    plus a fallback that detects job-disappearance via ``qstat -j``.
    """

    name = "sge"

    def __init__(self, cfg: Optional[dict] = None):
        cfg = cfg or {}
        self.queue: Optional[str] = cfg.get("queue") or None
        # Common values: smp (UGE default), threaded, openmp, mpi.
        self.parallel_env: str = cfg.get("parallel_env", "smp")
        # h_vmem is most portable; some clusters require mem_free or s_vmem.
        self.memory_flag: str = cfg.get("memory_flag", "h_vmem")
        # When True (default), YAML/SLURM-style ``mem`` is treated as *total*
        # job memory and divided by ``cpus`` before writing #$ -l <flag>=...,
        # because h_vmem / mem_free are per-slot on most SGE sites. Set False
        # only if your YAML values are already per-slot.
        self.memory_per_slot: bool = bool(cfg.get("memory_per_slot", True))
        self.hostname_exclude: str = cfg.get("hostname_exclude", "") or ""
        self.default_module_load: str = cfg.get("module_load", "") or ""
        self.extra_qsub_flags: list[str] = list(cfg.get("extra_qsub_flags", []) or [])

    # ── identity ─────────────────────────────────────────────────────────
    def task_id_env(self) -> str:
        return "SGE_TASK_ID"

    def array_index_base(self) -> int:
        return 1

    # ── helpers ──────────────────────────────────────────────────────────
    def _to_sge_hostname_expr(self, raw: str) -> str:
        """Translate a SLURM-style comma list into SGE hostname-exclude syntax.

        - "host1,host2"            -> "!host1&!host2"
        - "!h1&!h2" (already SGE)  -> returned unchanged
        - "" / None                -> ""

        The "&" form is the most portable across UGE / OGE / SGE-classic.
        """
        raw = raw.strip()
        if not raw:
            return ""
        if "&" in raw or raw.startswith("!"):
            return raw  # already SGE form
        hosts = [h.strip() for h in raw.split(",") if h.strip()]
        return "&".join(f"!{h}" for h in hosts)

    @staticmethod
    def _mem_for_slots(mem: str, cpus: int) -> str:
        """Convert total job memory to a per-slot request (ceil division).

        ``mem`` uses the same strings as SLURM ``--mem`` (e.g. ``16G``, ``512M``).
        """
        import math

        cpus = max(1, int(cpus or 1))
        s = str(mem).strip().upper()
        if s.endswith("G"):
            total, unit = float(s[:-1]), "G"
        elif s.endswith("M"):
            total, unit = float(s[:-1]), "M"
        else:
            total, unit = float(s), "G"
        per = math.ceil(total / cpus)
        if unit == "G" and per < 1:
            # Tiny totals: request at least 1M per slot rather than 0G.
            return f"{max(1, math.ceil(total * 1024 / cpus))}M"
        return f"{max(1, int(per))}{unit}"

    # ── header ────────────────────────────────────────────────────────────
    def header(
        self,
        *,
        job_name: str,
        cpus: int,
        mem: str,
        walltime: str,
        log_out: str,
        log_err: str,
        partition: Optional[str] = None,  # accepted but unused (SLURM-only)
        queue: Optional[str] = None,
        exclude: Optional[str] = None,
        array: Optional[tuple[int, int]] = None,
        max_concurrent: Optional[int] = None,
        dependency: Optional[str] = None,
        module_load: Optional[str] = None,
        extra_directives: Optional[Sequence[str]] = None,
    ) -> str:
        queue = queue or self.queue
        if exclude is None:
            exclude = self.hostname_exclude
        if module_load is None:
            module_load = self.default_module_load

        # YAML / SLURM ``mem`` is total job memory. On most SGE sites h_vmem and
        # mem_free are per-slot, so divide by cpus unless memory_per_slot=False.
        mem_req = (
            self._mem_for_slots(mem, cpus) if self.memory_per_slot else mem
        )
        mem_expr = f"{self.memory_flag}={mem_req}"
        # walltime: h_rt is universally supported.
        walltime_expr = f"h_rt={walltime}"

        lines: list[str] = [
            "#!/bin/bash",
            f"#$ -N {job_name}",
            "#$ -cwd",                     # run in submit dir (matches `set -euo pipefail` expectations)
            "#$ -V",                       # propagate environment
            "#$ -j n",                     # keep stdout/stderr split
        ]
        if queue:
            lines.append(f"#$ -q {queue}")
        if cpus and cpus > 1:
            lines.append(f"#$ -pe {self.parallel_env} {cpus}")
        lines.append(f"#$ -l {mem_expr}")
        lines.append(f"#$ -l {walltime_expr}")
        host_expr = self._to_sge_hostname_expr(exclude or "")
        if host_expr:
            # Single-quoted to protect the '!' from history expansion.
            lines.append(f"#$ -l hostname='{host_expr}'")
        if array:
            start, end = array
            if start < 1 or end < 1:
                raise ValueError(f"SGE arrays must be 1-based, got {array!r}")
            lines.append(f"#$ -t {start}-{end}")
            if max_concurrent and max_concurrent < (end - start + 1):
                lines.append(f"#$ -tc {max_concurrent}")
        if dependency:
            # SGE accepts a comma list of job IDs or names after -hold_jid.
            lines.append(f"#$ -hold_jid {dependency}")
        lines.append(f"#$ -o {log_out}")
        lines.append(f"#$ -e {log_err}")
        for flag in self.extra_qsub_flags:
            lines.append(f"#$ {flag}")
        if extra_directives:
            for d in extra_directives:
                lines.append(d if d.startswith("#$") else f"#$ {d}")
        if module_load:
            lines.append("")
            lines.append(f"module load {module_load}")
        # See SlurmScheduler.header: drop leaked PYTHONPATH; jobs must still
        # conda-activate so PATH resolves to the env python (not system python3).
        lines.append("")
        lines.append("unset PYTHONPATH")
        lines.append("")
        return "\n".join(lines) + "\n"

    def _depends_on_job_id(self, job_id: str) -> str:
        return str(job_id)

    # ── submit ────────────────────────────────────────────────────────────
    def submit(self, script_path: str | os.PathLike) -> str:
        script_path = str(script_path)
        logger.info(f"qsub {script_path}")
        try:
            result = _submit_with_retries(["qsub", script_path])
        except FileNotFoundError as e:
            raise RuntimeError("qsub not found on PATH; is this an SGE host?") from e
        m = _QSUB_JOB_ID_RE.search(result.stdout)
        if not m:
            raise RuntimeError(
                f"Could not parse SGE job ID from qsub output:\nSTDOUT: {result.stdout!r}\nSTDERR: {result.stderr!r}"
            )
        job_id = m.group(1)
        logger.info(f"Submitted SGE job {job_id}")
        return job_id

    # ── wait ──────────────────────────────────────────────────────────────
    def _job_present(self, job_id: str) -> bool:
        """True if SGE still knows about *job_id* (queued or running)."""
        try:
            res = subprocess.run(
                ["qstat", "-j", str(job_id)], capture_output=True, text=True
            )
        except FileNotFoundError as e:
            raise RuntimeError("qstat not found on PATH; is this an SGE host?") from e
        # qstat -j returns 0 if the job exists, non-zero (with a "do not exist"
        # message on stderr) when it has finished. The exit code is what UGE,
        # OGE, and SGE-classic all share; the stderr wording does not.
        return res.returncode == 0

    def wait(
        self,
        job_id: str,
        *,
        num_tasks: Optional[int] = None,
        timeout_s: int = 12 * 3600,
        check_interval_s: int = 30,
        sentinel_dir: Optional[str | os.PathLike] = None,
    ) -> JobResult:
        if sentinel_dir is None:
            raise ValueError(
                "SgeScheduler.wait requires sentinel_dir. Have each task write "
                "its exit code into <sentinel_dir>/done.<SGE_TASK_ID> via the "
                "trap_sentinel() helper."
            )
        sentinel_dir = Path(sentinel_dir)
        total = num_tasks if (num_tasks and num_tasks > 0) else 1
        expected = [sentinel_dir / f"done.{i}" for i in range(1, total + 1)]
        # For single (non-array) jobs the task ID env var is not set, so we
        # also accept "done.0" (matches sentinel_lines default fallback).
        if total == 1:
            expected.append(sentinel_dir / "done.0")

        logger.info(f"Waiting for SGE job {job_id} ({total} task{'s' if total > 1 else ''}, timeout {timeout_s / 3600:.1f}h)")
        start = time.time()
        last_report = 0.0
        status_report_interval = 300.0

        while time.time() - start < timeout_s:
            present_sentinels = [p for p in expected if p.exists()]
            # Deduplicate for the single-task case where both done.0 and done.1 may match.
            done_ids = {p.name for p in present_sentinels}
            done = len(done_ids) if total == 1 else len([p for p in expected[:total] if p.exists()])
            if done >= total:
                # Read exit codes (treat unreadable / missing as 0 in the rare
                # case of a partially-written file racing with us).
                failed = 0
                for p in present_sentinels:
                    try:
                        code = p.read_text().strip() or "0"
                    except OSError:
                        code = "0"
                    if code != "0":
                        failed += 1
                if failed:
                    return JobResult(
                        ok=False, completed=total - failed, failed=failed, total=total,
                        detail=f"{failed}/{total} tasks exited non-zero",
                    )
                return JobResult(ok=True, completed=total, total=total)

            # Fallback: if all sentinels haven't appeared but the job is gone
            # from qstat, assume crash.
            if not self._job_present(job_id):
                # Give the filesystem a moment in case sentinels are racing.
                time.sleep(2)
                done_ids = {p.name for p in expected if p.exists()}
                done2 = len(done_ids) if total == 1 else sum(1 for p in expected[:total] if p.exists())
                if done2 >= total:
                    failed = 0
                    for p in expected:
                        if not p.exists():
                            continue
                        try:
                            code = p.read_text().strip() or "0"
                        except OSError:
                            code = "0"
                        if code != "0":
                            failed += 1
                    if failed:
                        return JobResult(
                            ok=False, completed=total - failed, failed=failed, total=total,
                            detail=f"{failed}/{total} tasks exited non-zero",
                        )
                    return JobResult(ok=True, completed=total, total=total)
                return JobResult(
                    ok=False, completed=done2, failed=total - done2, total=total,
                    detail=f"job {job_id} no longer in queue but only {done2}/{total} sentinels present",
                )

            elapsed = time.time() - start
            if elapsed - last_report > status_report_interval:
                logger.info(f"SGE job {job_id} progress: {done}/{total} done ({elapsed / 60:.1f} min)")
                last_report = elapsed
            time.sleep(check_interval_s)

        return JobResult(ok=False, detail=f"timeout after {timeout_s / 3600:.1f}h")

    def cancel(self, job_id: str) -> None:
        try:
            subprocess.run(["qdel", str(job_id)], check=False)
        except FileNotFoundError:
            logger.warning("qdel not available; cannot cancel job")

    def job_present(self, job_id: str) -> bool:
        return self._job_present(job_id)


# ──────────────────────────────────────────────────────────────────────────────
# Local (single-machine)
# ──────────────────────────────────────────────────────────────────────────────


class LocalScheduler(Scheduler):
    """Run job bodies in a local subprocess.

    Strips scheduler-specific directive lines and runs the body with
    ``bash -euo pipefail -c``. Array jobs run serially in a foreground loop
    with TASK_ID iterating from start to end.
    """

    name = "local"

    def __init__(self, cfg: Optional[dict] = None):
        # cfg is accepted for symmetry with the other backends.
        self._next_id = 1
        # Track running jobs for cancel(): currently a no-op since wait()
        # runs synchronously, but stub anchored for future async paths.
        self._known: dict[str, dict] = {}

    def task_id_env(self) -> str:
        return "CAT_LOCAL_TASK_ID"

    def array_index_base(self) -> int:
        return 1

    def header(
        self,
        *,
        job_name: str,
        cpus: int,
        mem: str,
        walltime: str,
        log_out: str,
        log_err: str,
        partition: Optional[str] = None,
        queue: Optional[str] = None,
        exclude: Optional[str] = None,
        array: Optional[tuple[int, int]] = None,
        max_concurrent: Optional[int] = None,
        dependency: Optional[str] = None,
        module_load: Optional[str] = None,
        extra_directives: Optional[Sequence[str]] = None,
    ) -> str:
        # Header is just the shebang + a comment block documenting the request;
        # nothing is honoured by Local but the comments help when scripts are
        # later promoted to a real cluster.
        lines: list[str] = ["#!/bin/bash"]
        lines.append(f"# job_name={job_name}")
        lines.append(f"# cpus={cpus} mem={mem} walltime={walltime}")
        if array:
            start, end = array
            lines.append(f"# array={start}-{end} max_concurrent={max_concurrent or 'unbounded'}")
        if module_load:
            lines.append("")
            lines.append(f"module load {module_load} 2>/dev/null || true")
        lines.append("")
        # Stash the array range and log paths so submit() can act on them.
        self._pending = {
            "array": array, "log_out": log_out, "log_err": log_err,
        }
        return "\n".join(lines) + "\n"

    def _depends_on_job_id(self, job_id: str) -> str:
        return job_id  # Local jobs already ran by the time submit() returns.

    def submit(self, script_path: str | os.PathLike) -> str:
        script_path = str(script_path)
        pending = getattr(self, "_pending", {}) or {}
        array = pending.get("array")
        log_out = pending.get("log_out")
        log_err = pending.get("log_err")
        job_id = f"local-{self._next_id}"
        self._next_id += 1
        body = self._strip_directives(Path(script_path).read_text())
        if array:
            start, end = array
            results: list[int] = []
            for i in range(start, end + 1):
                rc = self._run_one(body, env_extra={"CAT_LOCAL_TASK_ID": str(i)}, log_out=log_out, log_err=log_err)
                results.append(rc)
            self._known[job_id] = {"results": results, "total": end - start + 1}
        else:
            rc = self._run_one(body, env_extra={"CAT_LOCAL_TASK_ID": "0"}, log_out=log_out, log_err=log_err)
            self._known[job_id] = {"results": [rc], "total": 1}
        return job_id

    def wait(
        self,
        job_id: str,
        *,
        num_tasks: Optional[int] = None,
        timeout_s: int = 12 * 3600,
        check_interval_s: int = 30,
        sentinel_dir: Optional[str | os.PathLike] = None,
    ) -> JobResult:
        info = self._known.get(job_id)
        if info is None:
            return JobResult(ok=False, detail=f"unknown local job {job_id}")
        total = info["total"]
        failed = sum(1 for rc in info["results"] if rc != 0)
        return JobResult(
            ok=(failed == 0), completed=total - failed, failed=failed, total=total,
            detail="" if failed == 0 else f"{failed}/{total} tasks failed",
        )

    def cancel(self, job_id: str) -> None:
        # Local jobs are synchronous; nothing to cancel after submit() returns.
        return

    def job_present(self, job_id: str) -> bool:
        # Local jobs are synchronous; once submit() returns they're done.
        return False

    # ── private helpers ───────────────────────────────────────────────────
    @staticmethod
    def _strip_directives(text: str) -> str:
        out_lines: list[str] = []
        for line in text.splitlines():
            stripped = line.lstrip()
            if stripped.startswith("#SBATCH") or stripped.startswith("#$"):
                continue
            if line.startswith("#!"):
                continue  # we re-launch via bash -c, no need for shebang
            out_lines.append(line)
        return "\n".join(out_lines)

    @staticmethod
    def _run_one(body: str, *, env_extra: dict[str, str], log_out: Optional[str], log_err: Optional[str]) -> int:
        env = os.environ.copy()
        env.update(env_extra)
        stdout = open(log_out, "a") if log_out else None
        stderr = open(log_err, "a") if log_err else None
        try:
            return subprocess.call(["bash", "-c", body], env=env, stdout=stdout, stderr=stderr)
        finally:
            if stdout:
                stdout.close()
            if stderr:
                stderr.close()


# ──────────────────────────────────────────────────────────────────────────────
# Factory
# ──────────────────────────────────────────────────────────────────────────────


_VALID_MODES = ("slurm", "sge", "local")


def detect_execution_mode() -> str:
    """Best-effort auto-detection of the batch scheduler on this host.

    Only unambiguous signals are honoured so we never emit directives for the
    wrong scheduler:

    - ``sbatch`` on PATH                     -> ``slurm``
    - ``qsub`` on PATH *and* ``$SGE_ROOT`` set -> ``sge``   (``$SGE_ROOT``
      distinguishes Grid Engine from PBS/Torque, which also ship ``qsub``)
    - otherwise                              -> ``local``   (runs everywhere)

    Users on a scheduler that cannot be auto-detected (PBS/Torque, LSF, or an
    SGE install without ``$SGE_ROOT``) should set ``execution_mode`` explicitly.
    """
    if shutil.which("sbatch"):
        return "slurm"
    if shutil.which("qsub") and os.environ.get("SGE_ROOT"):
        return "sge"
    return "local"


def resolve_execution_mode(mode: Optional[str]) -> str:
    """Resolve a possibly-``auto``/empty execution_mode to a concrete backend."""
    mode = (mode or "auto").lower()
    if mode == "auto":
        detected = detect_execution_mode()
        logger.info(f"execution_mode=auto resolved to '{detected}'")
        return detected
    return mode


def get_scheduler(mode: str, config: Optional[dict] = None) -> Scheduler:
    """Construct the Scheduler for *mode*, reading site-specific config.

    *mode* may be ``auto`` (or empty), in which case the scheduler is detected
    from the environment via :func:`detect_execution_mode`.

    *config* is the full snakemake config dict (or any dict with a "cluster"
    sub-block). Both ``config['cluster']['slurm']`` / ``config['cluster']['sge']``
    and the legacy ``config['slurm']`` layout are accepted.
    """
    mode = resolve_execution_mode(mode)
    if mode not in _VALID_MODES:
        raise ValueError(f"execution_mode must be one of {_VALID_MODES} or 'auto', got {mode!r}")

    config = config or {}
    cluster_cfg = config.get("cluster", {}) if isinstance(config, dict) else {}

    if mode == "slurm":
        # New layout: cluster.slurm.{partition,exclude_nodes,module_load}
        # Legacy layout: slurm.{partition,exclude_nodes,module_load}
        sub = cluster_cfg.get("slurm") if isinstance(cluster_cfg, dict) else None
        if not sub:
            sub = config.get("slurm", {}) if isinstance(config, dict) else {}
        return SlurmScheduler(sub or {})
    if mode == "sge":
        sub = cluster_cfg.get("sge") if isinstance(cluster_cfg, dict) else None
        return SgeScheduler(sub or {})
    return LocalScheduler()
