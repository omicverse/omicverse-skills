from .registry import list_skills, load_skill_text, skill_root
from .notebook_index import (
    SKILL_NOTEBOOK_MAP,
    build_index,
    check_freshness,
    format_freshness_report,
    index_path,
    load_index,
)

__all__ = [
    "list_skills",
    "load_skill_text",
    "skill_root",
    # notebook-index API
    "SKILL_NOTEBOOK_MAP",
    "build_index",
    "check_freshness",
    "format_freshness_report",
    "index_path",
    "load_index",
]

__version__ = "0.3.1"
