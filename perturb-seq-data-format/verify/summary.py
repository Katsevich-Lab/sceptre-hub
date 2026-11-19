"""Prints what every check of the verification suite found, and exits with an
error if any failed.

    python verify/summary.py
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

rows = defaultdict(dict)
for path in sorted(Path("verify/out/results").glob("*.json")):
    name, language = path.stem.rsplit("-", 1)
    rows[name][language] = json.loads(path.read_text())
cases = {c["name"]: c for c in json.loads(Path("conformance/cases.json").read_text())}
for language in ("python", "r"):
    results = json.loads(Path(f"conformance/results-{language}.json").read_text())
    rows["conformance"][language] = [
        {"check": name, "result": "pass" if sorted(r["rules"]) == cases[name]["expect"] and r["conforms"] == cases[name]["conforms"]
         and False not in [v for k, v in r.items() if k != "conforms"] else "fail", "detail": ""}
        for name, r in results.items()]

failed, observations = [], []
print(f"{'suite':22} {'language':8} {'checks':>6} {'pass':>5} {'fail':>5} {'observations':>13}")
for name in sorted(rows):
    for language in sorted(rows[name]):
        checks = rows[name][language]
        fails = [c for c in checks if c["result"] == "fail"]
        notes = [c for c in checks if c["result"] not in ("pass", "fail")]
        failed += [(name, language, c) for c in fails]
        observations += [(name, language, c) for c in notes]
        print(f"{name:22} {language:8} {len(checks):>6} {len(checks) - len(fails) - len(notes):>5} {len(fails):>5} {len(notes):>13}")
if observations:
    print("\nObservations (recorded, not verdicts):")
    seen = set()
    for name, language, c in observations:
        key = (c["result"], c["detail"])
        if key not in seen:
            seen.add(key)
            n = sum(1 for _, _, o in observations if (o["result"], o["detail"]) == key)
            print(f"  [{name}, {language}] {c['check']}: {c['result']}" + (f"  (and {n - 1} more like it)" if n > 1 else ""))
for name, language, c in failed:
    print(f"FAIL [{name}, {language}] {c['check']}: {c['detail']}")
print("\n" + ("ALL CHECKS PASS" if not failed else f"{len(failed)} CHECK(S) FAILED"))
sys.exit(1 if failed else 0)
