import h5py
import numpy as np
import os
import gzip
import shutil

SRC_TYPES = ["high_coverage", "high_quality"]
SRC_BASENAMES = {"high_coverage": "HC", "high_quality": "HQ"}

SNPS_DIR = "snps"
INDELS_DIR = "indels"
OUTPUT_DIR = "merged"
CHROMOSOMES = [f"{n}H" for n in range(1, 8)]

def merge_hdf5_pair(snp_path, indel_path, output_path):
    with h5py.File(snp_path, 'r') as snp_file, h5py.File(indel_path, 'r') as indel_file:
        with h5py.File(output_path, 'w') as out_file:
            columns = list(snp_file.keys())
            arrays = {}
            for col in columns:
                snp_arr = snp_file[col][:]
                indel_arr = indel_file[col][:]
                arrays[col] = np.concatenate([snp_arr, indel_arr], axis=0)
            # Determine sorting keys
            pos = arrays['POS']
            # Decode POS values if necessary
            #if isinstance(pos[0], bytes):
            #    pos = np.array([int(p.decode('utf-8')) for p in pos])
            #else:
            pos = np.array([int(p) for p in pos])
            if 'CHROM' in arrays:
                chrom = arrays['CHROM']
                # decode chrom for proper string order
                #if isinstance(chrom[0], bytes):
                #    chrom = np.array([c.decode('utf-8') for c in chrom])
                chrom = np.array([c for c in chrom])
                sort_idx = np.lexsort((pos, chrom))
            else:
                sort_idx = np.argsort(pos)
            for col in columns:
                arr = arrays[col]
                arr_sorted = arr[sort_idx]
                out_file.create_dataset(col, data=arr_sorted, dtype=snp_file[col].dtype)

def gzip_file(infile, outfile):
    with open(infile, 'rb') as f_in, gzip.open(outfile, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(infile)

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for chrom in CHROMOSOMES:
        for typ in SRC_TYPES:
            snp = os.path.join(SNPS_DIR, f"chr{chrom}.{typ}.h5")
            indel = os.path.join(INDELS_DIR, f"chr{chrom}.{typ}.h5")
            tag = SRC_BASENAMES[typ]
            out_unzipped = os.path.join(OUTPUT_DIR, f"barley_chr{chrom}_{tag}.h5")
            out_gz = out_unzipped + ".gz"
            print(f"Merging {snp} and {indel} into {out_gz}")
            merge_hdf5_pair(snp, indel, out_unzipped)
            os.system(f"pigz {out_unzipped}")
            #gzip_file(out_unzipped, out_gz)

if __name__ == "__main__":
    main()
