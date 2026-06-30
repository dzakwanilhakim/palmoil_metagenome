"""
epi2me_extract.py
Extract taxonomy abundance tables from EPI2ME wf-metagenomics / wf-16s HTML reports.

Data lives in HTML <table> tags with the rank name as the first column header.
Works for: phylum, class, order, family, genus, species

Usage:
    from epi2me_extract import extract_abundance_from_html
    df = extract_abundance_from_html("report.html", "genus")
"""

import re
import pandas as pd
from io import StringIO
from pathlib import Path


VALID_RANKS = ["phylum", "class", "order", "family", "genus", "species"]


def extract_abundance_from_html(html_path: str, rank: str = "genus") -> pd.DataFrame:
    """
    Extract an abundance table (taxa × samples) from an EPI2ME HTML report.

    Args:
        html_path : path to the .html report file
        rank      : taxonomic rank — one of phylum/class/order/family/genus/species

    Returns:
        pd.DataFrame  index = taxon names, columns = sample/barcode names, values = read counts
    """
    rank_lower = rank.lower()
    if rank_lower not in VALID_RANKS:
        raise ValueError(f"Invalid rank '{rank}'. Choose from: {VALID_RANKS}")

    html = Path(html_path).read_text(encoding="utf-8")

    # EPI2ME embeds data as plain HTML <table> elements inside DataTable widgets.
    # Each table has the rank name as its first <th scope="col"> header.
    table_ids = re.findall(r'id="(DataTable_[^"]+_inner)"', html)
    if not table_ids:
        raise ValueError(
            f"No DataTable elements found in {html_path}. "
            "Is this a valid EPI2ME wf-metagenomics or wf-16s report?"
        )

    for tid in table_ids:
        idx = html.find(f'id="{tid}"')
        if idx == -1:
            continue
        end_idx = html.find("</table>", idx)
        if end_idx == -1:
            continue
        table_html = html[idx : end_idx + 8]

        headers = re.findall(r'<th[^>]*scope="col"[^>]*>(.*?)</th>', table_html)
        if not headers or headers[0].strip().lower() != rank_lower:
            continue

        df = pd.read_html(StringIO(f"<table>{table_html}</table>"))[0]
        df = df.set_index(df.columns[0])
        df.index.name = rank_lower
        return df

    raise ValueError(
        f"No table with rank '{rank_lower}' found in {html_path}. "
        f"Available ranks in this file: {_detect_available_ranks(html)}"
    )


def _detect_available_ranks(html: str) -> list[str]:
    """Return which taxonomy ranks are present in the HTML."""
    found = []
    table_ids = re.findall(r'id="(DataTable_[^"]+_inner)"', html)
    for tid in table_ids:
        idx = html.find(f'id="{tid}"')
        if idx == -1:
            continue
        end_idx = html.find("</table>", idx)
        if end_idx == -1:
            continue
        table_html = html[idx : end_idx + 8]
        headers = re.findall(r'<th[^>]*scope="col"[^>]*>(.*?)</th>', table_html)
        if headers and headers[0].strip().lower() in VALID_RANKS:
            rank = headers[0].strip().lower()
            if rank not in found:
                found.append(rank)
    return found