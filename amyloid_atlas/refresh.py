"""Rebuild amyloid_atlas.json from the live Amyloid Atlas page.

Run by .github/workflows/refresh-amyloid-atlas.yml on a weekly schedule (and
by hand: `python amyloid_atlas/refresh.py`). helicon's own
`helicon.lib.dataset.get_amyloid_atlas` used to scrape the page directly, but
that page sends no CORS header, so a page fetching it from a different origin
(the Shinylive/Pyodide build of `helicon webApps`) cannot load it at all.
`get_amyloid_atlas` now reads this file's raw.githubusercontent.com URL
instead, which does send one, and this script is what keeps the file current.

Needs ``helicon[html]`` or ``pandas[html]`` installed (html5lib, bs4) for
``pandas.read_html``, beyond plain ``helicon``.
"""

import sys
from pathlib import Path

import helicon
import pandas as pd
from helicon.lib.dataset import EMDB

OUT = Path(__file__).parent / "amyloid_atlas.json"

# a PDB ID the atlas carries that no longer resolves to it in EMDB's own table
REPLACED_PDB_IDS = {"7z40": "8ade"}


def build(url: str = "https://people.mbi.ucla.edu/sawaya/amyloidatlas") -> pd.DataFrame:
    df = pd.read_html(url, header=0, flavor="html5lib")[0]

    mask = df["PDB ID"].isin(REPLACED_PDB_IDS)
    df.loc[mask, "PDB ID"] = df.loc[mask, "PDB ID"].str.lower().map(REPLACED_PDB_IDS)
    df = df[df["Method"].str.lower() == "cryoem"].copy()

    emdb = EMDB()
    assert emdb.meta is not None, "Failed to get the list of EMDB entries"
    pdb2emd = {
        pdb_id: row["emd_id"]
        for _, row in emdb.meta.iterrows()
        for pdb_id in str(row["pdb"]).lower().split(",")
        if pdb_id
    }
    df["emd_id"] = df["PDB ID"].map(pdb2emd)

    df["sample"] = df["Protein"] + " - " + df["Fibril Origins"]
    # found by name rather than hard-coded: the page's own text for the "Å"
    # in this header is not always the same Unicode codepoint as what one
    # types, and a plain string literal here has silently stopped matching
    # before
    resolution_col = next(c for c in df.columns if "esol" in c)
    cols = {
        "emd_id": "emd_id",
        "resolution": resolution_col,
        "pdb_id": "PDB ID",
        "sample": "sample",
    }
    others = ["Residues Ordered", "Reference"]
    df = df[helicon.flatten([*cols.values(), others])]
    df = df.rename(columns={v: k for k, v in cols.items()})
    df = df.drop_duplicates(subset=["emd_id", "pdb_id"])
    df = df.dropna(subset=["emd_id"])
    return df.reset_index(drop=True)


def main() -> int:
    df = build()
    if len(df) < 100:
        # the page's own table failed to parse as expected (a layout change,
        # a transient error page, ...) sooner than a usable table shrinking
        # by three quarters in real life
        print(f"got only {len(df)} rows, expected several hundred; not writing", file=sys.stderr)
        return 1
    df.to_json(OUT, orient="records", indent=2)
    print(f"wrote {len(df)} rows to {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
