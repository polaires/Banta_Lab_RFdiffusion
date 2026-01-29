"""
Inference Backend — Abstracts in-process vs HTTP execution.

Provides InProcessBackend (wraps inference_utils) and HTTPBackend
(wraps api_client) behind a common interface so pipeline modules
don't care which backend is in use.
"""

import inspect
import logging
from typing import Any, Callable, Dict, Set

logger = logging.getLogger(__name__)


def _accepted_params(func: Callable) -> Set[str]:
    """Get the set of parameter names a function accepts.

    Returns all param names. If the function accepts **kwargs,
    returns an empty set (meaning "accept anything").
    """
    try:
        sig = inspect.signature(func)
        has_var_keyword = any(
            p.kind == inspect.Parameter.VAR_KEYWORD
            for p in sig.parameters.values()
        )
        if has_var_keyword:
            return set()  # empty = accept everything
        return set(sig.parameters.keys())
    except (ValueError, TypeError):
        return set()


def _filter_params(params: Dict[str, Any], func: Callable) -> Dict[str, Any]:
    """Filter params dict to only keys accepted by func.

    Strips 'task' key unconditionally. If the function's accepted
    params can't be determined (has **kwargs), passes everything.
    """
    allowed = _accepted_params(func)
    filtered = {k: v for k, v in params.items() if k != "task"}
    if allowed:
        filtered = {k: v for k, v in filtered.items() if k in allowed}
    return filtered


class InProcessBackend:
    """Runs inference directly via Foundry Python modules (GPU required)."""

    def __init__(self):
        self._rfd3 = None
        self._mpnn = None
        self._rf3 = None

    def _ensure_imports(self):
        if self._rfd3 is None:
            from inference_utils import (
                run_rfd3_inference,
                run_mpnn_inference,
                run_rf3_inference,
            )
            self._rfd3 = run_rfd3_inference
            self._mpnn = run_mpnn_inference
            self._rf3 = run_rf3_inference

    def run_rfd3(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Run RFD3 backbone generation in-process."""
        self._ensure_imports()
        return self._rfd3(**_filter_params(params, self._rfd3))

    def run_mpnn(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Run MPNN sequence design in-process."""
        self._ensure_imports()
        return self._mpnn(**_filter_params(params, self._mpnn))

    def run_rf3(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Run RF3 structure prediction in-process."""
        self._ensure_imports()
        return self._rf3(**_filter_params(params, self._rf3))


class HTTPBackend:
    """Runs inference via HTTP API (Docker container must be running)."""

    def __init__(self, api_url: str = None, timeout: int = 600):
        from api_client import API_URL, DEFAULT_TIMEOUT
        self.api_url = api_url or API_URL
        self.timeout = timeout or DEFAULT_TIMEOUT

    def run_rfd3(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Run RFD3 via HTTP API."""
        from api_client import call_api
        result = call_api("rfd3", params, timeout=self.timeout)
        if result.get("status") == "failed":
            raise RuntimeError(f"RFD3 API error: {result.get('error')}")
        return result.get("result", {})

    def run_mpnn(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Run MPNN via HTTP API."""
        from api_client import call_api
        result = call_api("mpnn", params, timeout=self.timeout)
        if result.get("status") == "failed":
            raise RuntimeError(f"MPNN API error: {result.get('error')}")
        return result.get("result", {})

    def run_rf3(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Run RF3 via HTTP API."""
        from api_client import call_api
        result = call_api("rf3", params, timeout=self.timeout)
        if result.get("status") == "failed":
            raise RuntimeError(f"RF3 API error: {result.get('error')}")
        return result.get("result", {})


def create_backend(mode: str = "auto", **kwargs) -> Any:
    """
    Factory function to create the appropriate backend.

    Args:
        mode: "in_process", "http", or "auto" (auto-detects Foundry availability)
        **kwargs: Passed to the backend constructor

    Returns:
        InProcessBackend or HTTPBackend instance
    """
    if mode == "in_process":
        return InProcessBackend()
    elif mode == "http":
        return HTTPBackend(**kwargs)
    elif mode == "auto":
        try:
            from inference_utils import check_foundry_available
            if check_foundry_available():
                logger.info("Auto-detected Foundry — using InProcessBackend")
                return InProcessBackend()
        except ImportError:
            pass
        logger.info("Foundry not available — using HTTPBackend")
        return HTTPBackend(**kwargs)
    else:
        raise ValueError(f"Unknown backend mode: {mode!r}. Use 'in_process', 'http', or 'auto'.")
