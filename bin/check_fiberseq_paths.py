import os
import pandas as pd
import glob
import sys
import re
from tqdm import tqdm
import argparse

tqdm.pandas()

def letter_index(letter):
    # Convert the letter to uppercase to ensure consistency
    letter = letter.upper()
    # Calculate the index by subtracting 64 from the ASCII value
    index = ord(letter) - 64
    return index

def sizeof_fmt(path, suffix="B"):
    num = os.path.getsize(path)
    for unit in ("", "Ki", "Mi", "Gi", "Ti", "Pi", "Ei", "Zi"):
        if abs(num) < 1024.0:
            return f"{num:3.1f}{unit}{suffix}"
        num /= 1024.0
    return f"{num:.1f}Yi{suffix}"

def check_wells(base_path, well_id, fname):
    d = letter_index(well_id[0])
    files = glob.glob(f"{base_path}/*_{well_id}/{fname}") 
    new_well_id = f"{d}_{well_id}"
    assert files == glob.glob(f"{base_path}/{new_well_id}/{fname}"), "Well ID does not match the expected pattern"
    assert len(files) <= 1, f"Expected 1 bam file, found {len(files)}: {files}. base_path: {base_path}, well_id: {well_id}, fname: {fname}"
    if len(files) != 1:
        raise ValueError(f"Expected 1 bam file, found {len(files)}: {files}. base_path: {base_path}, well_id: {well_id}, fname: {fname}")
    return new_well_id

    
def check_bam_files(row):
    platform = row['Instrument']
    well_id = row['Sample Well ID']
    base_path = row['base_path']
    barcode = None if pd.isna(row['Barcode']) else row['Barcode']
    assert platform in ('Sequel 2e', 'Revio'), f"Platform is not Sequel 2e or Revio: {platform}"
    
    reads_type = 'hifi'
    if platform == 'Sequel 2e':
        if barcode is not None:
            barcode = f"{barcode}--{barcode}"
            fname = f'{barcode}/*.hifi_reads.{barcode}.bam'
        else:
            fname = "*.hifi_reads.bam"
        
            
        if not re.match(r'^\d+_', well_id):
            try:
                well_id = check_wells(base_path, well_id, fname)
            except ValueError as e:
                if barcode is None:
                    print('hifi reads not found')
                    fname = "*.subreads.bam"
                    well_id = check_wells(base_path, well_id, fname)
                    reads_type = 'subreads'
                else:
                    raise e
    else:
        if barcode is not None:
            fname = f'hifi_reads/*.hifi_reads.{barcode}.bam'
        else:
            fname = 'hifi_reads/*.hifi_reads.bam'
    # Find all bam files matching the pattern
    files = glob.glob(f"{base_path}/{well_id}/{fname}")
    row['Sample Well ID fixed'] = well_id
    try:
        assert len(files) == 1, f"Expected 1 bam file, found {len(files)}: {files}. base_path: {base_path}, well_id: {well_id}, fname: {fname}"
        result = files[0]
        tmp = os.path.basename(result).split('.')
        expected_len = 3 if barcode is None else 4
        assert len(tmp) == expected_len, f"Expected {expected_len} parts in the filename, found {len(tmp)}: {tmp}"
        row['bam'] = result
        row['bam_size'] = sizeof_fmt(result)
        row['Flowcell ID'] = tmp[0]
        row['reads_type'] = reads_type
        
        return row
    except AssertionError as e:
        if len(files) == 0:
            print(f"No bam files found for {base_path}/{well_id}/{fname}")
            row['bam'] = None
            row['Flowcell ID'] = None
            row['bam_size'] = None
            row['reads_type'] = None
            return row
        raise e


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fill in missing bam files in the metadata file')
    parser.add_argument('input_file', type=str, help='Path to metadata file')
    parser.add_argument('output_file', type=str, help='Path to output file')
    parser.add_argument('--include_initial_columns', action='store_true', help='Include initial columns in the output file', default=False)
    args = parser.parse_args()

    df = pd.read_table(args.input_file, low_memory=False)
    base_path = '/net/seq/pacbio/runs/'
    df['base_path'] = base_path + df['Run ID (Data folder)']

    added_cols = ['Flowcell ID', 'Sample Well ID fixed', 'bam', 'bam_size']
    result_cols = list(df.columns) + added_cols if args.include_initial_columns else added_cols
    df.progress_apply(check_bam_files, axis=1)[result_cols].to_csv(args.output_file, sep='\t', index=False)
