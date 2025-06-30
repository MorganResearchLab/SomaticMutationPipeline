import logging
import pysam
import pandas as pd
import sys
import argparse

parser = argparse.ArgumentParser(
                    prog='split_by_celltype',
                    description='Splits BAM reads into cell type specific BAM files based on cell barcodes')
parser.add_argument("--bamfile",
                    help="Input BAM file")
parser.add_argument("--outdir",
                    help="Output directory for cell type split BAMs")
parser.add_argument("--chrom",
                    help="Chromosome identifier - check for chrNum or Num format in BAM header")
parser.add_argument("--metafile",
                    help="Meta data that maps cell barcode to cell type annotations")
parser.add_argument("--celltypes",
                    help="Comma-separated list of cell type annotations to split reads into")
parser.add_argument("--threads",
                    help="Threads to speed up cell type BAM splitting",
                    default=1, type=int)
parser.add_argument("--sra_meta", default=None,
                    help="Meta data that maps SRA ID to batch ID")
parser.add_argument("--log", default=sys.stdout,
                    help="logfile")
args = parser.parse_args()

if args.log != sys.stdout:
    logging.basicConfig(filename=args.log,
                        format='%(levelname)s:%(message)s', level=logging.DEBUG)
else:
    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

logging.info("Reading meta data file: {}".format(args.metafile))
metadata = pd.read_csv(args.metafile, sep="\t")

logging.info("Opening BAM file buffer: {}".format(args.bamfile))
inbam = pysam.AlignmentFile(args.bamfile, "rb", threads=args.threads)

# check the header reference naming convention
chr_contigs = any([rx.lower().startswith("chr") for rx in inbam.references])
is_upper = any([rx.lower().startswith("Chr") for rx in inbam.references])

chrome_st = args.chrom.startswith("chr")

if chr_contigs and is_upper:
    if not chrome_st:
        chrom_str = "Chr" + args.chrom
    else:
        chrom_str = "Chr" + args.chrom.lstrip("chr")
elif chr_contigs and not is_upper:
    if chrome_st:
        chrom_str = args.chrom
    else:
        chrom_str = "chr" + args.chrom
else:
    if chrome_st:
        chrom_str = args.chrom.lstrip("chr")
    else:
        chrom_str = args.chrom

if args.sra_meta is None:
    logging.info("no SRA meta-data provided - skipping")
    bam_name = args.bamfile.split('/')[-1].split('.')[0]
    batch_id = bam_name.split("_")[0]
else:
    logging.info("Reading sequencing sample meta data: {}".format(args.sra_meta))
    df = pd.read_csv(args.sra_meta)

    filename_column = df.columns[0]
    batch_id_column = df.columns[-1]

    bam_name = args.bamfile.split('/')[-1].split('.')[0]
    matching_row = df[df[filename_column] == bam_name]

    if not matching_row.empty:
        batch_id = matching_row[batch_id_column].values[0]
        logging.info(f"Batch ID for {bam_name}: {batch_id}")

    else:
        logging.info(f"No matching batch ID found for {bam_name}")
        batch_id = None
        

cell_types = args.celltypes.split(",")

logging.info("Extracting reads matching barcodes with {}".format(",".join(cell_types)))

outputbam = {}
for cell_type in cell_types:
   outputbam[cell_type] = {}
   if chr_contigs:
       outputbam[cell_type][chrom_str] = pysam.AlignmentFile("{}/{}_{}_{}.bam".format(args.outdir, bam_name, cell_type, chrom_str),
                                                             "wb", template=inbam, threads=args.threads)
   else:
       outputbam[cell_type][chrom_str] = pysam.AlignmentFile("{}/{}_{}_chr{}.bam".format(args.outdir, bam_name, cell_type, chrom_str),
                                                             "wb", template=inbam, threads=args.threads)

logging.info("Fetching reads on {}".format(chrom_str))
summ_dict = {"Reads_processed":0, "Valid": 0, "No_CB": 0, "CB_unmatched": 0, "CB_matched": 0}
count = 0

for x in (inbam.fetch(chrom_str)): 
    count += 1
    try:
     cellbarcode = x.get_tag("CB")
     summ_dict["Valid"] += 1
     # logging.info("CB tag: {}".format(cellbarcode)) # for debugging
     cb_id = batch_id + "_" + cellbarcode

     if not cellbarcode.endswith("-1"):
         cb_id = cb_id.rstrip("-1")
     
     x.set_tag("CB", cb_id)
     if metadata['Index'].isin([cb_id]).any():
         summ_dict["CB_matched"] += 1
         celltype = metadata.loc[metadata['Index'] == cb_id, 'Cell_type'].values[0]
         if count % 5000 == 0:
             logging.info("Read: {}. CB: {}. Celltype: {}.".format(x.query_name, cb_id, celltype))
             logging.info("Processed {} reads".format(count))
         outputbam[celltype][chrom_str].write(x)
     else:
         summ_dict["CB_unmatched"] += 1

    except KeyError :
        summ_dict["No_CB"] += 1

summ_dict["Reads_processed"] = count

for cell_type, chrom_bams in outputbam.items():
   for chrom, bam in chrom_bams.items():
   
       bam.close()
inbam.close()

logging.info("Writing summary results for {} processed reads".format(count))
sum_df = pd.DataFrame(summ_dict)
sum_df.to_csv("{}/{}_read_summary.tsv".format(args.outdir, bam_name))



