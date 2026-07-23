"""gff2mss: canonical INSDC GFF3 -> DDBJ MSS conversion."""

__all__ = ["convert", "emit_ann", "emit_fasta", "load_config", "load_common"]

# Lazy re-exports (PEP 562). Importing the gff2mss package itself must not pull
# in the optional 'ddbj-gff' dependency — that keeps `import gff2mss.cli` and
# `gff2mss --help` working without it, and lets the CLI emit a friendly error.
# Each name is resolved on first access, which is where ddbj-gff is truly needed.
_LAZY = {
    "load_common": ".config",
    "load_config": ".config",
    "convert": ".convert",
    "emit_ann": ".emit",
    "emit_fasta": ".emit",
}


def __getattr__(name):
    module = _LAZY.get(name)
    if module is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    import importlib
    return getattr(importlib.import_module(module, __name__), name)


def __dir__():
    return sorted(__all__)
