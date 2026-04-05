from __future__ import annotations

from pathlib import Path
from typing import Dict, List


def skill_root() -> Path:
    return Path(__file__).resolve().parent / "skills"


def _parse_frontmatter(text: str) -> Dict[str, str]:
    payload: Dict[str, str] = {}
    lines = (text or "").splitlines()
    if len(lines) < 3 or lines[0].strip() != "---":
        return payload
    for line in lines[1:]:
        if line.strip() == "---":
            break
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        payload[key.strip().lower()] = value.strip().strip('"').strip("'")
    return payload


def list_skills() -> List[Dict[str, str]]:
    items: List[Dict[str, str]] = []
    root = skill_root()
    if not root.exists():
        return items

    for skill_file in sorted(root.glob("*/SKILL.md")):
        text = skill_file.read_text(encoding="utf-8", errors="ignore")
        meta = _parse_frontmatter(text)
        skill_dir = skill_file.parent
        entry = {
            "slug": skill_dir.name,
            "name": meta.get("name", skill_dir.name),
            "title": meta.get("title", ""),
            "description": meta.get("description", ""),
            "path": str(skill_file),
            "reference_path": str(skill_dir / "reference.md") if (skill_dir / "reference.md").exists() else "",
        }
        items.append(entry)
    return items


def load_skill_text(slug: str) -> str:
    skill_file = skill_root() / slug / "SKILL.md"
    if not skill_file.exists():
        raise FileNotFoundError(f"Skill not found: {slug}")
    return skill_file.read_text(encoding="utf-8", errors="ignore")
