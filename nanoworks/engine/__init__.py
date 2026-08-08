"""DFT calculation backends used by Nanoworks."""

def normalize_engine_name(engine):
    """Return the canonical Nanoworks engine name."""
    return str(engine).strip().upper()
