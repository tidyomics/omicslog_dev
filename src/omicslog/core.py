from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any

import anndata as ad
import pandas as pd

LOG_KEY = "_omicslog"

def _safe_deepcopy_dict(d: dict) -> dict:
    result = {}
    for k, v in d.items():
        try:
            result[k] = deepcopy(v)
        except (TypeError, Exception):
            pass
    return result


def _timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _format_log_message(operation: str, message: str, ts: str | None = None) -> list[str]:
    stamp = ts or _timestamp()
    return [stamp, operation, message]


def _ensure_log_container(adata: ad.AnnData) -> pd.DataFrame:
    current = adata.uns.get(LOG_KEY)
    if not isinstance(current, list):
        adata.uns[LOG_KEY] = pd.DataFrame(columns=["Time","Operation","Message"])
    return adata.uns[LOG_KEY]


def _append_log_messages(adata: ad.AnnData, messages: list[str] | tuple[str, ...]) -> None:
    if not messages:
        return
    container = _ensure_log_container(adata)
    new_rows = pd.DataFrame(messages, columns=["Time", "Operation", "Message"])
    adata.uns[LOG_KEY] = pd.concat([container, new_rows], ignore_index=True)

def _inherit_and_append(
    parent: ad.AnnData,
    child: ad.AnnData,
    messages: list[str] | tuple[str, ...],
) -> None:
    child.uns[LOG_KEY] = _ensure_log_container(parent)
    _append_log_messages(child, messages)

def _parent_set(obj, attr: str, value) -> None:
    """Set an attribute via the first parent class that defines it.
    Works for both plain properties (.fset) and custom descriptors (.__set__).
    """
    for base in type(obj).__mro__[1:]:
        if attr in base.__dict__:
            base.__dict__[attr].__set__(obj, value)
            return
    object.__setattr__(obj, attr, value)

class _LoggingProxy:
    """
    Transparent proxy for dict-like AnnData components (layers, obsm, varm, ...).
    Intercepts __setitem__ and __delitem__ to log mutations automatically.
    """

    def __init__(self, wrapped, owner: "LoggedAnnDataStandalone", label: str):
        object.__setattr__(self, "_w", wrapped)
        object.__setattr__(self, "_owner", owner)
        object.__setattr__(self, "_label", label)

    def __setitem__(self, key: str, value) -> None:
        w = object.__getattribute__(self, "_w")
        owner = object.__getattribute__(self, "_owner")
        label = object.__getattribute__(self, "_label")
        verb = "updated" if key in w else "added"
        w[key] = value
        _append_log_messages(owner, [_format_log_message(label, f"'{key}' {verb}")])

    def __delitem__(self, key: str) -> None:
        w = object.__getattribute__(self, "_w")
        owner = object.__getattribute__(self, "_owner")
        label = object.__getattribute__(self, "_label")
        del w[key]
        _append_log_messages(owner, [_format_log_message(label, f"'{key}' removed")])

    def __getitem__(self, key):
        return object.__getattribute__(self, "_w")[key]

    def __getattr__(self, name):
        return getattr(object.__getattribute__(self, "_w"), name)

    def __contains__(self, key):
        return key in object.__getattribute__(self, "_w")

    def __iter__(self):
        return iter(object.__getattribute__(self, "_w"))

    def __len__(self):
        return len(object.__getattribute__(self, "_w"))

    def __repr__(self):
        return repr(object.__getattribute__(self, "_w"))


@dataclass
class AnnDataSnapshot:
    """Captures the full state of an AnnData object for diffing."""
    n_obs: int
    n_vars: int
    obs_cols: list[str] = field(default_factory=list)
    var_cols: list[str] = field(default_factory=list)
    layers: list[str] = field(default_factory=list)
    obsm: list[str] = field(default_factory=list)
    varm: list[str] = field(default_factory=list)
    obsp: list[str] = field(default_factory=list)
    varp: list[str] = field(default_factory=list)

    @classmethod
    def from_anndata(cls, adata: ad.AnnData) -> "AnnDataSnapshot":
        return cls(
            n_obs=adata.n_obs,
            n_vars=adata.n_vars,
            obs_cols=list(adata.obs.columns),
            var_cols=list(adata.var.columns),
            layers=list(adata.layers.keys()),
            obsm=list(adata.obsm.keys()),
            varm=list(adata.varm.keys()),
            obsp=list(adata.obsp.keys()),
            varp=list(adata.varp.keys()),
        )


def _diff_keys(
    pre: list[str],
    post: list[str],
    label: str,
    operation: str,
    ts: str,
) -> list[str]:
    msgs = []
    for k in sorted(set(pre) - set(post)):
        msgs.append(_format_log_message(operation, f"{label} removed: '{k}'", ts))
    for k in sorted(set(post) - set(pre)):
        msgs.append(_format_log_message(operation, f"{label} added: '{k}'", ts))
    return msgs


def _subset_messages(
    pre: AnnDataSnapshot,
    post: AnnDataSnapshot,
    operation: str = "subset",
) -> list[str]:
    msgs: list[str] = []
    ts = _timestamp()

    if pre.n_vars != post.n_vars:
        removed = pre.n_vars - post.n_vars
        pct = round((removed / pre.n_vars) * 100) if pre.n_vars else 0
        msgs.append(_format_log_message(
            operation,
            f"removed {removed} genes ({pct}%), {post.n_vars} genes remaining",
            ts,
        ))

    if pre.n_obs != post.n_obs:
        removed = pre.n_obs - post.n_obs
        pct = round((removed / pre.n_obs) * 100) if pre.n_obs else 0
        msgs.append(_format_log_message(
            operation,
            f"removed {removed} samples ({pct}%), {post.n_obs} samples remaining",
            ts,
        ))

    msgs += _diff_keys(pre.obs_cols, post.obs_cols, "obs column", operation, ts)
    msgs += _diff_keys(pre.var_cols, post.var_cols, "var column", operation, ts)
    msgs += _diff_keys(pre.layers,   post.layers,   "layer",      operation, ts)
    msgs += _diff_keys(pre.obsm,     post.obsm,     "obsm",       operation, ts)
    msgs += _diff_keys(pre.varm,     post.varm,     "varm",       operation, ts)
    msgs += _diff_keys(pre.obsp,     post.obsp,     "obsp",       operation, ts)
    msgs += _diff_keys(pre.varp,     post.varp,     "varp",       operation, ts)

    return msgs


class LoggedAnnDataStandalone(ad.AnnData):
    """Standalone subclass strategy with local logging helpers and message style."""

    def __init__(self, *args: Any, **kwargs: Any):
        super().__init__(*args, **kwargs)
        _ensure_log_container(self)

    @classmethod
    def _safe_component_copy(cls, value):
        return value.copy() if hasattr(value, "copy") else deepcopy(value)

    @classmethod
    def from_anndata(cls, adata: ad.AnnData) -> "LoggedAnnDataStandalone":
        if isinstance(adata, cls):
            _ensure_log_container(adata)
            return adata

        kwargs: dict[str, Any] = {
            "X": cls._safe_component_copy(adata.X) if adata.X is not None else None,
            "obs": adata.obs.copy(),
            "var": adata.var.copy(),
            "uns": _safe_deepcopy_dict(dict(adata.uns)),
            "obsm": {k: cls._safe_component_copy(v) for k, v in adata.obsm.items()},
            "varm": {k: cls._safe_component_copy(v) for k, v in adata.varm.items()},
            "layers": {k: cls._safe_component_copy(v) for k, v in adata.layers.items()},
            "obsp": {k: cls._safe_component_copy(v) for k, v in adata.obsp.items()},
            "varp": {k: cls._safe_component_copy(v) for k, v in adata.varp.items()},
        }
        
        if adata.raw is not None:
            kwargs["raw"] = {
                "X": cls._safe_component_copy(adata.raw.X),
                "var": adata.raw.var.copy(),
                "varm": {k: cls._safe_component_copy(v) for k, v in adata.raw.varm.items()},
            }

        logged = cls(**kwargs)
        _ensure_log_container(logged)
        return logged

    # --- proxied properties: each needs a getter AND a setter ---

    @property
    def layers(self):
        return _LoggingProxy(super().layers, self, "layers")

    @layers.setter
    def layers(self, value):
        _parent_set(self, "layers", value)

    @property
    def obsm(self):
        return _LoggingProxy(super().obsm, self, "obsm")

    @obsm.setter
    def obsm(self, value):
        _parent_set(self, "obsm", value)

    @property
    def varm(self):
        return _LoggingProxy(super().varm, self, "varm")

    @varm.setter
    def varm(self, value):
        _parent_set(self, "varm", value)

    @property
    def obsp(self):
        return _LoggingProxy(super().obsp, self, "obsp")

    @obsp.setter
    def obsp(self, value):
        _parent_set(self, "obsp", value)

    @property
    def varp(self):
        return _LoggingProxy(super().varp, self, "varp")

    @varp.setter
    def varp(self, value):
        _parent_set(self, "varp", value)

    @property
    def obs(self):
        return _LoggingProxy(super().obs, self, "obs")

    @obs.setter
    def obs(self, value):
        _parent_set(self, "obs", value)

    @property
    def var(self):
        return _LoggingProxy(super().var, self, "var")

    @var.setter
    def var(self, value):
        ad.AnnData.var.fset(self, value)

    # --- snapshot & subsetting ---

    def _snapshot(self) -> AnnDataSnapshot:
        return AnnDataSnapshot.from_anndata(self)

    def __getitem__(self, index):
        pre = self._snapshot()
        result = super().__getitem__(index)
        logged_result = self.from_anndata(result)
        msgs = _subset_messages(pre, logged_result._snapshot(), operation="subset")
        _inherit_and_append(self, logged_result, msgs)
        return logged_result

    def _inplace_subset(self, index):
        pre = self._snapshot()
        super()._inplace_subset(index)
        _append_log_messages(self, _subset_messages(pre, self._snapshot(), operation="subset"))

    def _operation_log_block(self) -> str:
        logs = self.uns.get(LOG_KEY, [])
        if isinstance(logs, pd.DataFrame):
            if logs.empty:
                return ""
            rows = logs.apply(lambda r: f"[{r['Time']}] {r['Operation']}: {r['Message']}", axis=1)
            return "\n\nOperation log:\n" + "\n".join(rows)
        if not logs:
            return ""
        return "\n\nOperation log:\n" + "\n".join(str(x) for x in logs)

    def __repr__(self) -> str:
        return super().__repr__() + self._operation_log_block()

    def __str__(self) -> str:
        return self.__repr__()

    def operation_log(self) -> list[str]:
        return list(self.uns.get(LOG_KEY, []))


def log_start(adata: ad.AnnData) -> LoggedAnnDataStandalone:
    return LoggedAnnDataStandalone.from_anndata(adata)