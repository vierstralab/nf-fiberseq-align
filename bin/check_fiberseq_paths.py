import os
import pandas as pd
import glob
import sys
import re
from tqdm import tqdm

tqdm.pandas()

def letter_index(letter):
    # Convert the letter to uppercase to ensure consistency
    letter = letter.upper()
    # Calculate the index by subtracting 64 from the ASCII value
    index = ord(letter) - 64
    return index


def check_wells(base_path, well_id, fname):
    d = letter_index(well_id[0])
    files = glob.glob(f"{base_path}/*_{well_id}/{fname}") 
    new_well_id = f"{d}_{well_id}"
    assert files == glob.glob(f"{base_path}/{new_well_id}/{fname}"), "Well ID does not match the expected pattern"
    assert len(files) == 1, f"Expected 1 bam file, found {len(files)}: {files}. base_path: {base_path}, well_id: {well_id}, fname: {fname}"
    return new_well_id

    
def check_bam_files(row):
    platform = row['Instrument']
    well_id = row['Sample Well ID']
    base_path = row['base_path']
    barcode = None if pd.isna(row['Barcode']) else row['Barcode']
    assert platform in ('Sequel 2e', 'Revio'), f"Platform is not Sequel 2e or Revio: {platform}"

    if platform == 'Sequel 2e':
        if barcode is not None:
            barcode = f"{barcode}--{barcode}"
            fname = f'{barcode}/*.hifi_reads.{barcode}.bam'
        else:
            fname = "'*.subreads.bam'"
        if not re.match(r'^\d+_', well_id):
            well_id = check_wells(base_path, well_id, fname)

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
        row['Flowcell ID'] = tmp[0]
        
        return row[['Flowcell ID', 'Sample Well ID fixed', 'bam']]
    except AssertionError as e:
        if len(files) == 0:
            print(f"No bam files found for {base_path}/{well_id}/{fname}")
            row['bam'] = None
            row['Flowcell ID'] = None
            return row[['Flowcell ID', 'Sample Well ID fixed', 'bam']]
        raise e


if __name__ == '__main__':
    assert len(sys.argv) == 3, "Usage: python check_bam_files.py <input_file> <output_file>"
    df = pd.read_table(sys.argv[1])
    base_path = '/net/seq/pacbio/runs/'
    df['base_path'] = base_path + df['Run ID (Data folder)']
    df.progress_apply(check_bam_files, axis=1).to_csv(sys.argv[2], sep='\t', index=False)

    