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
    if not re.match(r'^\d+_', well_id):
        d = letter_index(well_id[0])
        files = glob.glob(f"{base_path}/*_{well_id}/{fname}") 
        new_well_id = f"{d}_{well_id}"
        new_well_id_files = glob.glob(f"{base_path}/{new_well_id}/{fname}")
        if len(files) == 1 and len(new_well_id_files) == 0:
            # hotfix for misc sheet
            print (f"Well_id doesn't match expected pattern {new_well_id} found {well_id} instead. base_path: {base_path}, fname: {fname}")
            new_well_id = files[0].split('/')[-2]
        else:
            assert files == new_well_id_files, f"Well ID does not match the expected pattern, multiple wells correspond to {well_id}. base_path: {base_path}, fname: {fname}"
    else:
        new_well_id = well_id
        files = new_well_id_files
    assert len(files) <= 1, f"Expected <1 bam file, found {len(files)}: {files}. base_path: {base_path}, well_id: {well_id}, fname: {fname}"
    if len(files) != 1:
        raise ValueError(f"Expected 1 bam file, found {len(files)}: {files}. base_path: {base_path}, well_id: {well_id}, fname: {fname}")
    return new_well_id

    
def check_bam_files(row):
    platform = row['Instrument']
    well_id = row['Sample Well ID']
    base_path = row['base_path']
    barcode = None if pd.isna(row['Barcode']) else row['Barcode']
    expected_len = 3 if barcode is None else 4
    assert platform in ('Sequel 2e', 'Revio'), f"Platform is not Sequel 2e or Revio: {platform}"
    
    reads_type = 'hifi'
    if platform == 'Sequel 2e':
        if barcode is not None:
            barcode = f"{barcode}--{barcode}"
            fname = f'{barcode}/*.hifi_reads.{barcode}.bam'
        else:
            fname = "*.hifi_reads.bam"
        
        try:
            well_id = check_wells(base_path, well_id, fname)
        except ValueError as e:
            if barcode is None:
                print(f'hifi reads not found for {base_path}. Trying subreads.')
                fname = "*.subreads.bam"
                try: 
                    well_id = check_wells(base_path, well_id, fname)
                    reads_type = 'subreads'
                except ValueError as e:
                    print(f'subreads not found for {base_path}. well_id: {well_id}. Setting reads to None.')
                    row['reads'] = None
                    row['reads_size'] = None
                    row['sample_id'] = None
                    row['reads_type'] = None
                    return row
            else:
                raise e
    else:
        if barcode is not None:
            fname = f'hifi_reads/*.hifi_reads.{barcode}.bam'
        else:
            fname = 'hifi_reads/*.hifi_reads.bam'
            # hotfix for misc sheet
            if len(glob.glob(f"{base_path}/{well_id}/{fname}")) == 0:
                fname = 'hifi_reads/*.hifi_reads.default.bam' 
                expected_len = 4

    # Find all bam files matching the pattern
    files = glob.glob(f"{base_path}/{well_id}/{fname}")
    row['Well ID fixed'] = well_id

    # FIXME repeated code. A lot of things don't work on misc sheet. All the format are weird there.
    if len(files) != 1:
        print(f"Expected 1 bam file, found {len(files)}: {files}. base_path: {base_path}, well_id: {well_id}, fname: {fname}")
        if len(files) == 0:
            row['reads'] = None
            row['reads_size'] = None
            row['sample_id'] = None
            row['reads_type'] = None
            return row
        raise AssertionError
    result = files[0]
    tmp = os.path.basename(result).split('.')
    assert len(tmp) == expected_len, f"Expected {expected_len} parts in the filename, found {len(tmp)}: {tmp}"
    row['reads'] = result
    row['reads_size'] = sizeof_fmt(result)
    row['sample_id'] = tmp[0]
    row['reads_type'] = reads_type
        
    return row


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fill in missing bam files in the metadata file')
    parser.add_argument('input_file', type=str, help='Path to metadata file')
    parser.add_argument('output_file', type=str, help='Path to output file')
    parser.add_argument('--include_initial_columns', action='store_true', help='Include initial columns in the output file', default=False)
    args = parser.parse_args()

    df = pd.read_table(args.input_file, low_memory=False)
    base_path = '/net/seq/pacbio/runs/'
    df['base_path'] = base_path + df['Run ID (Data folder)']

    added_cols = ['sample_id', 'Well ID fixed', 'reads', 'reads_type', 'reads_size']
    result_cols = list(df.columns) + added_cols if args.include_initial_columns else added_cols
    df.progress_apply(check_bam_files, axis=1)[result_cols].to_csv(args.output_file, sep='\t', index=False)
