#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import argparse
import shutil
import sys
import os
from tqdm import tqdm

def get_optitype_samplesheet(paths, seq_type='dna', sample_naming='stem'):
    """Create a DataFrame from a list of paths."""
    data = []
    for p in paths:
        if sample_naming == 'stem':
            sample_name = p.stem
        elif sample_naming == 'full':
            sample_name = p.as_posix()
        else:
            raise ValueError("sample_naming must be 'stem' or 'full'")

        data.append([sample_name, p.as_posix(), seq_type])

    df = pd.DataFrame(data, columns=['sample', 'bam', 'seq_type'])
    return df

def find_bam_files(base_path: Path, seq_type='dna', tissue_code=None, wes_suffix='S000021.bam',rna_suffix='mRNA.bam') -> list[Path]:
    """Finds BAM files, optionally filtering by a tissue code in the filename."""
    suffix_map = {'dna': wes_suffix, 'rna':rna_suffix}
    suffix = suffix_map.get(seq_type)
    if not suffix:
        raise ValueError(f"seq_type must be one of {list(suffix_map.keys())}")

    glob_pattern = f"*/*/02_bam/*{suffix}"
    print(f"Searching for files with pattern: {glob_pattern}")
    all_bam_files = list(base_path.glob(glob_pattern))
    print(f"Found {len(all_bam_files)} total matching files.")

    if tissue_code:
        print(f"Filtering for tissue code '{tissue_code}' (filename starts with 'CH_{tissue_code}_')...")
        filter_prefix = f"CH_{tissue_code}_"
        filtered_files = [p for p in all_bam_files if p.name.startswith(filter_prefix)]
        print(f"Found {len(filtered_files)} files matching the tissue code.")
        return filtered_files

    return all_bam_files

def organize_files_and_get_path_map(source_paths: list[str], destination_dir: Path, method: str, dry_run: bool) -> dict:
    """Organizes files into a destination directory by copying or symlinking them."""
    action_verb = method.capitalize()

    if not dry_run:
        destination_dir.mkdir(parents=True, exist_ok=True)
        print(f"\n{action_verb}ing {len(source_paths)} unique files to '{destination_dir}'...")
    else:
        print(f"\n--- Dry Run: Would {method} {len(source_paths)} unique files to '{destination_dir}' ---")

    path_map = {}

    for src_str in tqdm(source_paths, desc=f"{action_verb}ing files", disable=dry_run):
        src_path = Path(src_str)
        dest_path = destination_dir / src_path.name
        path_map[src_str] = dest_path.as_posix()

        if not dry_run:
            if not src_path.exists():
                tqdm.write(f"Warning: Source file not found, skipping: {src_path}")
                continue

            if dest_path.exists() or dest_path.is_symlink():
                dest_path.unlink()

            if method == 'copy':
                shutil.copy2(src_path, dest_path)
            elif method == 'symlink':
                os.symlink(src_path.resolve(), dest_path)

    return path_map

def main(args):
    """Main execution function."""

    data_dir = Path(args.data_dir)
    results_dir = Path(args.results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    # --- 1. Find BAM files, with optional tissue filtering ---
    print("Gathering paths for all datasets...")
    bam_files = find_bam_files(data_dir, seq_type=args.seq_type, tissue_code=args.tissue_code)
    if not bam_files:
        print("No BAM files found for the specified criteria. Exiting.")
        sys.exit(0)

    df_with_duplicates = get_optitype_samplesheet(bam_files, seq_type=args.seq_type,sample_naming='stem')
    df_no_duplicates = get_optitype_samplesheet(bam_files, seq_type=args.seq_type,sample_naming='full')

    # --- 2. Define filenames based on whether a tissue was specified ---
    filename_suffix = f"_{args.tissue_code}" if args.tissue_code else ""

    out_local_dups = results_dir / f'samplesheet_local_with_duplicates{filename_suffix}.csv'
    out_local_unique = results_dir / f'samplesheet_local_unique_names{filename_suffix}.csv'

    df_with_duplicates.to_csv(out_local_dups, index=False)
    df_no_duplicates.to_csv(out_local_unique, index=False)
    print(f"\nLocal samplesheet with duplicate-prone names saved to:\n  {out_local_dups}")
    print(f"Local samplesheet with unique names saved to:\n  {out_local_unique}")

    # --- 3. Optionally organize files and create updated samplesheets ---
    if args.organize_in:
        destination_dir = Path(args.organize_in)
        unique_source_paths = df_no_duplicates['bam'].unique().tolist()

        path_map = organize_files_and_get_path_map(unique_source_paths, destination_dir, args.method, args.dry_run)

        df_organized_dups = df_with_duplicates.copy()
        df_organized_dups['bam'] = df_organized_dups['bam'].map(path_map)

        df_organized_unique = df_no_duplicates.copy()
        df_organized_unique['bam'] = df_organized_unique['bam'].map(path_map)

        out_organized_dups = results_dir / f'samplesheet_{args.method}_with_duplicates{filename_suffix}.csv'
        out_organized_unique = results_dir / f'samplesheet_{args.method}_unique_names{filename_suffix}.csv'

        if not args.dry_run:
            df_organized_dups.to_csv(out_organized_dups, index=False)
            df_organized_unique.to_csv(out_organized_unique, index=False)
            print(f"\n{args.method.capitalize()} path samplesheet with duplicate-prone names saved to:\n  {out_organized_dups}")
            print(f"{args.method.capitalize()} path samplesheet with unique names saved to:\n  {out_organized_unique}")
        else:
            print("\n--- Dry Run: Would have created the following samplesheets with updated paths ---")
            print(f"  - {out_organized_dups}")
            print(f"  - {out_organized_unique}")

    print("\n--- Process complete ---")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Find WES BAM files, create samplesheets, and optionally copy or symlink them to a new location.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '--data-dir',
        default='/data/PRTNR/CHUV/DIR/rgottar1/owkin_spatial/cegat',
        help="The base directory to search for BAM files."
    )
    parser.add_argument(
        '--results-dir',
        default='../../results/WES/hlatyping/',
        help="The directory where output samplesheets will be saved."
    )
    parser.add_argument(
        '--tissue-code',
        default='G',
        help="Optional: A single letter to filter files by tissue type.\n(e.g., 'G' for Glioblastoma files like 'CH_G_...')"
    )
    parser.add_argument(
        '--organize-in',
        default='../../results/WES/hlatyping/bam',
        help="Optional: A destination directory to organize all found BAM files into.\nIf specified, two additional samplesheets with updated paths will be created."
    )
    parser.add_argument(
        '--method',
        choices=['copy', 'symlink'],
        default='copy',
        help="Method to use for organizing files if --organize-in is specified.\n(default: 'copy')"
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help="Use with --organize-in to simulate the process without creating files or links."
    )
    parser.add_argument(
        '--seq-type',
        default='dna',
        help="dna or rna (to find mRNA or WES bam files)"
    )
    args = parser.parse_args()
    if not args.organize_in:
        print("Note: --organize-in not specified. Only local samplesheets will be created.")

    main(args)
