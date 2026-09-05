"""Unit tests for the annotation worker-timeout fallback in
eggnogmapper.annotator.e7.annotate.AnnotationEngine.annotate_batch.

A sub-batch that exceeds the per-task deadline must be recomputed in-process
rather than aborting the whole batch. The in-process path is the same code the
pool workers execute, so the recomputed result is identical -- only slower.

Only the parallel dispatch/collect control flow is exercised here; the pool and
the annotation engine itself are mocked.
"""

import multiprocessing
import logging
import pytest
from unittest.mock import MagicMock

from eggnogmapper.annotator.e7.annotate import AnnotationEngine


class _FakeAsyncResult:
    """Stands in for multiprocessing.pool.AsyncResult."""

    def __init__(self, result=None, time_out=False):
        self._result = result
        self._time_out = time_out

    def get(self, timeout=None):
        if self._time_out:
            raise multiprocessing.TimeoutError()
        return self._result


def _engine_with_pool(async_results, inproc_calls):
    """Build an AnnotationEngine whose pool returns the given AsyncResults."""
    engine = AnnotationEngine.__new__(AnnotationEngine)

    pool = MagicMock()
    pool.apply_async.side_effect = list(async_results)

    def _inproc(seed_ids, **kwargs):
        inproc_calls.append(list(seed_ids))
        return {s: {"recomputed": True} for s in seed_ids}

    engine._annotate_batch_inproc = _inproc
    return engine, pool


def test_timed_out_subbatch_is_recomputed_in_process(caplog):
    """One slow sub-batch must not discard the seeds that already succeeded."""
    inproc_calls = []
    engine, pool = _engine_with_pool(
        [
            _FakeAsyncResult({1: {"ok": True}, 2: {"ok": True}}),
            _FakeAsyncResult(time_out=True),           # the pathological one
            _FakeAsyncResult({5: {"ok": True}, 6: {"ok": True}}),
        ],
        inproc_calls,
    )

    with caplog.at_level(logging.WARNING):
        merged = engine.annotate_batch(
            [1, 2, 3, 4, 5, 6], pool=pool, sub_batch_size=2
        )

    # every seed is present -- previously this raised RuntimeError and lost all six
    assert sorted(merged) == [1, 2, 3, 4, 5, 6]
    # exactly the timed-out sub-batch was recomputed, and with the right seeds
    assert inproc_calls == [[3, 4]]
    assert merged[3] == {"recomputed": True}
    assert "in-process" in caplog.text


def test_no_timeout_leaves_results_untouched(monkeypatch):
    """The happy path must not invoke the in-process fallback at all."""
    inproc_calls = []
    engine, pool = _engine_with_pool(
        [
            _FakeAsyncResult({1: {"ok": True}, 2: {"ok": True}}),
            _FakeAsyncResult({3: {"ok": True}, 4: {"ok": True}}),
        ],
        inproc_calls,
    )

    merged = engine.annotate_batch([1, 2, 3, 4], pool=pool, sub_batch_size=2)

    assert sorted(merged) == [1, 2, 3, 4]
    assert inproc_calls == []


def test_env_var_overrides_the_deadline(monkeypatch):
    """EGGNOG_ANNOT_TASK_TIMEOUT replaces the size-derived deadline."""
    seen = {}

    class _Recording(_FakeAsyncResult):
        def get(self, timeout=None):
            seen["timeout"] = timeout
            return {1: {"ok": True}, 2: {"ok": True}}

    inproc_calls = []
    engine, pool = _engine_with_pool([_Recording(), _Recording()], inproc_calls)
    monkeypatch.setenv("EGGNOG_ANNOT_TASK_TIMEOUT", "900")

    engine.annotate_batch([1, 2, 3], pool=pool, sub_batch_size=2)

    assert seen["timeout"] == 900


def test_invalid_env_var_falls_back_to_computed_deadline(monkeypatch, caplog):
    """A non-integer override warns and keeps max(120, sub_batch_size * 2)."""
    seen = {}

    class _Recording(_FakeAsyncResult):
        def get(self, timeout=None):
            seen["timeout"] = timeout
            return {1: {"ok": True}, 2: {"ok": True}}

    inproc_calls = []
    engine, pool = _engine_with_pool([_Recording(), _Recording()], inproc_calls)
    monkeypatch.setenv("EGGNOG_ANNOT_TASK_TIMEOUT", "not-a-number")

    # 3 seeds > sub_batch_size 2, so the parallel path (and the deadline) is used
    with caplog.at_level(logging.WARNING):
        engine.annotate_batch([1, 2, 3], pool=pool, sub_batch_size=2)

    assert seen["timeout"] == 120          # max(120, 2 * 2)
    assert "not an integer" in caplog.text
