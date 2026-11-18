"""Shared by the checks: recording results, and making sure that
neither sceptre nor pysceptre was loaded."""
import json
import sys
from pathlib import Path

OUT = Path("verify/out")
(OUT / "results").mkdir(parents=True, exist_ok=True)


class Results:
    def __init__(self, name):
        self.name, self.rows = name, []

    def add(self, check, ok, detail=""):
        """`ok` is True (pass), False (fail) or a string (an observation, not a verdict)."""
        self.rows.append({"check": check, "result": "pass" if ok is True else "fail" if ok is False else str(ok), "detail": str(detail)})

    def save(self):
        loaded = sorted({m.split(".")[0] for m in sys.modules} & {"pysceptre", "sceptre"})
        self.add("neither sceptre nor pysceptre was loaded", not loaded, ", ".join(loaded))
        (OUT / "results" / f"{self.name}-python.json").write_text(json.dumps(self.rows, indent=1) + "\n")
        failed = [r for r in self.rows if r["result"] == "fail"]
        print(f"{self.name} (Python): {len(self.rows) - len(failed)} of {len(self.rows)} checks pass")
        for r in failed:
            print("  FAIL", r["check"], "-", r["detail"])
