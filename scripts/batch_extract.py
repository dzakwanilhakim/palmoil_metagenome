"""
batch_extract.py
Batch extract taxonomy abundance tables from all EPI2ME HTML reports in a directory.
Outputs one TSV per HTML file: [report_filename]_[rank].tsv

Usage:
    python batch_extract.py data/raw --rank genus
    python batch_extract.py data/raw --all-ranks
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from epi2me_extract import extract_abundance_from_html, VALID_RANKS


def batch_extract(html_dir: str, rank: str = "genus", all_ranks: bool = False):
    html_dir = Path(html_dir)
    html_files = sorted(html_dir.glob("*.html"))

    if not html_files:
        print(f"[!] No .html files found in {html_dir}")
        sys.exit(1)

    ranks = VALID_RANKS if all_ranks else [rank]

    print(f"[→] {len(html_files)} HTML files | ranks: {ranks}\n")

    errors = []

    for html_file in html_files:
        print(f"[FILE] {html_file.name}")
        for rank_name in ranks:
            out_path = html_dir / f"{html_file.stem}_{rank_name}.tsv"
            try:
                df = extract_abundance_from_html(str(html_file), rank_name)
                df.to_csv(out_path, sep="\t")
                print(f"  [✓] {rank_name}: {df.shape[0]} taxa × {df.shape[1]} samples → {out_path.name}")
            except Exception as e:
                msg = f"  [✗] {rank_name}: {e}"
                print(msg)
                errors.append(f"{html_file.name} / {rank_name}: {e}")
        print()

    print("─" * 50)
    if errors:
        print(f"[!] {len(errors)} error(s):")
        for err in errors:
            print(f"    {err}")
        sys.exit(1)
    else:
        print("[✓] All done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Batch extract EPI2ME abundance tables from all HTML files in a directory"
    )
    parser.add_argument("html_dir", help="Directory containing EPI2ME .html report files")
    parser.add_argument(
        "--rank", default="genus",
        choices=VALID_RANKS,
        help="Taxonomic rank to extract (default: genus)",
    )
    parser.add_argument(
        "--all-ranks", action="store_true",
        help="Extract all ranks: phylum, class, order, family, genus, species",
    )
    args = parser.parse_args()

    batch_extract(args.html_dir, rank=args.rank, all_ranks=args.all_ranks)