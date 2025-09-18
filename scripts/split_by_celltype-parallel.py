import concurrent.futures
import logging
import pysam
import pandas as pd
import sys
import argparse
import os
import time
from functools import partial

def setup_logging(log_file, log_level=logging.INFO):
    """Set up logging configuration"""
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    if log_file != "stdout":
        logging.basicConfig(filename=log_file, format=log_format, level=log_level)
    else:
        logging.basicConfig(format=log_format, level=log_level)

def process_region(region, bamfile, outdir, metadata, cell_types, batch_id, batch_position, bam_name, chr_contigs, chrom_str, log_file=None):
    """Process a specific genomic region"""
    # Set up logging for this process
    if log_file:
        setup_logging(log_file)
    
    region_id = f"{region[0]}:{region[1]}-{region[2]}"
    worker_id = os.getpid()  # Get the process ID for logging
    start_time = time.time()
    
    logging.info(f"Worker {worker_id} starting processing region: {region_id}")
    
    # Create temporary output files for this region
    region_str = f"{region[0]}_{region[1]}_{region[2]}"
    temp_outputbam = {}
    
    inbam = pysam.AlignmentFile(bamfile, "rb", threads=1)
    
    # Create temporary directory if it doesn't exist
    temp_dir = os.path.join(outdir, "temp")
    os.makedirs(temp_dir, exist_ok=True)
    
    for cell_type in cell_types:
        temp_outputbam[cell_type] = {}
        if chr_contigs:
            out_path = f"{temp_dir}/temp_{bam_name}_{cell_type}_{chrom_str}_{region_str}.bam"
        else:
            out_path = f"{temp_dir}/temp_{bam_name}_{cell_type}_chr{chrom_str}_{region_str}.bam"
        temp_outputbam[cell_type][chrom_str] = pysam.AlignmentFile(out_path, "wb", template=inbam, threads=1)
    
    # Process this region
    summ_dict = {"Reads_processed": 0, "Valid": 0, "No_CB": 0, "CB_unmatched": 0, "CB_matched": 0}
    reads_by_celltype = {}
    for cell_type in cell_types:
        reads_by_celltype[cell_type] = 0
    
    count = 0
    logging_interval = 10000  # Log every 10,000 reads
    
    for x in inbam.fetch(region[0], region[1], region[2]):
        count += 1
        try:
            cellbarcode = x.get_tag("CB")
            summ_dict["Valid"] += 1

            if batch_position == "prefix":
                cb_id = batch_id + "_" + cellbarcode
            elif batch_position == "suffix":
                cb_id = cellbarcode + "_" + batch_id
            elif batch_position == "None":
                cb_id = cellbarcode
            else:
                raise IOError(f"Batch position argument {batch_position} not recognised")

            if not cellbarcode.endswith("-1"):
                cb_id = cb_id.rstrip("-1")
                
            x.set_tag("CB", cb_id)
            if metadata['Index'].isin([cb_id]).any():
                summ_dict["CB_matched"] += 1
                celltype = metadata.loc[metadata['Index'] == cb_id, 'Cell_type'].values[0]

                if count % logging_interval == 0:
                    elapsed = time.time() - start_time
                    rate = count / elapsed if elapsed > 0 else 0
                    logging.info(f"Worker {worker_id} - Region {region_id} - "
                                f"Read: {x.query_name}. CB: {cb_id}. Celltype: {celltype}. "
                                f"Processed {count} reads in {elapsed:.2f}s ({rate:.1f} reads/s)")
                
                temp_outputbam[celltype][chrom_str].write(x)
                reads_by_celltype[celltype] = reads_by_celltype.get(celltype, 0) + 1
            else:
                summ_dict["CB_unmatched"] += 1
        except KeyError:
            summ_dict["No_CB"] += 1
    
    summ_dict["Reads_processed"] = count
    
    # Close all temporary files
    for cell_type, chrom_bams in temp_outputbam.items():
        for chrom, bam in chrom_bams.items():
            bam.close()
    
    inbam.close()
    
    # Calculate output file sizes
    output_file_sizes = {}
    for cell_type in cell_types:
        if chr_contigs:
            temp_path = f"{temp_dir}/temp_{bam_name}_{cell_type}_{chrom_str}_{region_str}.bam"
        else:
            temp_path = f"{temp_dir}/temp_{bam_name}_{cell_type}_chr{chrom_str}_{region_str}.bam"
        
        if os.path.exists(temp_path):
            output_file_sizes[cell_type] = os.path.getsize(temp_path) / (1024 * 1024)  # Size in MB
    
    # Calculate performance metrics
    end_time = time.time()
    elapsed_time = end_time - start_time
    processing_rate = count / elapsed_time if elapsed_time > 0 else 0
    
    # Log detailed completion information
    if count == 0:
        count += 1
        
    logging.info(f"Worker {worker_id} COMPLETED region {region_id}:")
    logging.info(f"  - Processing time: {elapsed_time:.2f} seconds")
    logging.info(f"  - Reads processed: {count} ({processing_rate:.1f} reads/second)")
    logging.info(f"  - Valid reads: {summ_dict['Valid']} ({summ_dict['Valid']*100/count:.1f}% of total)")
    logging.info(f"  - Reads with cell barcodes matched: {summ_dict['CB_matched']} ({summ_dict['CB_matched']*100/count:.1f}% of total)")
    
    # Log cell type distribution
    logging.info(f"  - Reads by cell type:")
    for cell_type, read_count in reads_by_celltype.items():
        if read_count > 0:
            logging.info(f"    * {cell_type}: {read_count} reads, output size: {output_file_sizes.get(cell_type, 0):.2f} MB")
    
    return region_id, summ_dict, elapsed_time, reads_by_celltype

def main():
    # Your existing argument parsing code
    parser = argparse.ArgumentParser(
                        prog='split_by_celltype',
                        description='Splits BAM reads into cell type specific BAM files based on cell barcodes')
    # Add existing arguments
    parser.add_argument("--bamfile", help="Input BAM file")
    parser.add_argument("--outdir", help="Output directory for cell type split BAMs")
    parser.add_argument("--chrom", help="Chromosome identifier - check for chrNum or Num format in BAM header")
    parser.add_argument("--metafile", help="Meta data that maps cell barcode to cell type annotations")
    parser.add_argument("--celltypes", help="Comma-separated list of cell type annotations to split reads into")
    parser.add_argument("--threads", help="Threads to speed up cell type BAM splitting", default=1, type=int)
    parser.add_argument("--sra_meta", default=None, help="Meta data that maps SRA ID to batch ID")
    parser.add_argument("--batch_column", default=None, help="Column of meta data containing batch information")
    parser.add_argument("--batch_pos", choices=["prefix", "suffix", "None"], nargs="?",
                        help="Where to place the batch ID on the CB", default="suffix", const="suffix")
    parser.add_argument("--log", default="stdout", help="logfile")
    # Add new arguments
    parser.add_argument("--chunk_size", default=10000000, type=int, 
                        help="Size of genomic chunks (in bp) for parallelization")
    parser.add_argument("--logging_interval", default=10000, type=int,
                        help="Number of reads between logging updates")
    args = parser.parse_args()
    
    # Set up logging - now using string instead of file object
    log_file = args.log
    setup_logging(log_file)
    
    start_time = time.time()
    logging.info(f"Starting parallel BAM splitting with {args.threads} threads")
    
    # Read metadata
    logging.info(f"Reading meta data file: {args.metafile}")
    metadata = pd.read_csv(args.metafile, sep="\t")
    logging.info(f"Loaded metadata with {len(metadata)} entries")
    
    # Open BAM to get header information
    logging.info(f"Opening BAM file buffer: {args.bamfile}")
    inbam = pysam.AlignmentFile(args.bamfile, "rb")
    
    # Check chromosome naming
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
    
    logging.info(f"Using chromosome identifier: {chrom_str}")
    
    # Get chromosome length from header
    chrom_length = 0
    for ref in inbam.references:
        if ref == chrom_str:
            chrom_length = inbam.get_reference_length(ref)
            break
    
    if chrom_length == 0:
        logging.error(f"Could not find chromosome {chrom_str} in BAM header")
        sys.exit(1)
    
    logging.info(f"Chromosome {chrom_str} length: {chrom_length:,} bp")
    
    # Determine batch ID
    if args.sra_meta is None:
        logging.info("No SRA meta-data provided - using filename for batch ID")
        bam_name = args.bamfile.split('/')[-1].split('.')[0]
        batch_id = bam_name.split("_")[0]
    else:
        logging.info(f"Reading sequencing sample meta data: {args.sra_meta}")
        df = pd.read_csv(args.sra_meta)
        
        filename_column = df.columns[0]
        if args.batch_column is not None:
            batch_id_column = args.batch_column
        else:           
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
    logging.info(f"Extracting reads matching barcodes with {len(cell_types)} cell types: {','.join(cell_types)}")
    
    # Create temp directory
    temp_dir = os.path.join(args.outdir, "temp")
    os.makedirs(temp_dir, exist_ok=True)
    
    # Divide the chromosome into chunks
    chunk_size = args.chunk_size
    regions = []
    for start in range(0, chrom_length, chunk_size):
        end = min(start + chunk_size, chrom_length)
        regions.append((chrom_str, start, end))
    
    num_regions = len(regions)
    logging.info(f"Divided {chrom_str} into {num_regions} regions for parallel processing")
    logging.info(f"Region size: {chunk_size:,} bp")
    
    # Process regions in parallel
    # Note: We're not using partial anymore to avoid serializing issues
    # Instead, we pass all parameters directly to the worker
    
    all_summaries = {}
    region_times = {}
    reads_by_region_celltype = {}
    
    # Create executor with the specified number of threads
    max_workers = min(args.threads, num_regions)  # Don't create more workers than regions
    logging.info(f"Using {max_workers} worker processes")
    
    completed_regions = 0
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks with explicit parameters (no partial function)
        future_to_region = {}
        for region in regions:
            future = executor.submit(
                process_region, 
                region=region,
                bamfile=args.bamfile,
                outdir=args.outdir,
                metadata=metadata,
                cell_types=cell_types,
                batch_id=batch_id,
                batch_position=args.batch_pos,
                bam_name=bam_name,
                chr_contigs=chr_contigs,
                chrom_str=chrom_str,
                log_file=log_file
            )
            future_to_region[future] = region
        
        for future in concurrent.futures.as_completed(future_to_region):
            region = future_to_region[future]
            region_id_str = f"{region[0]}:{region[1]}-{region[2]}"
            
            try:
                region_id, summary, elapsed_time, reads_by_celltype = future.result()
                all_summaries[region_id] = summary
                region_times[region_id] = elapsed_time
                reads_by_region_celltype[region_id] = reads_by_celltype
                
                completed_regions += 1
                logging.info(f"Main process received results from region {region_id} "
                            f"({completed_regions}/{num_regions} regions completed)")
                
                # Log some key statistics from the completed region
                logging.info(f"Region {region_id} processed {summary['Reads_processed']} reads "
                            f"in {elapsed_time:.2f} seconds "
                            f"({summary['Reads_processed']/elapsed_time:.1f} reads/second)")
                
            except Exception as exc:
                logging.error(f"Region {region_id_str} generated an exception: {exc}")
                raise
    
    logging.info(f"All {num_regions} regions have completed processing")
    
    # Merge the temporary BAM files for each cell type
    logging.info("Merging temporary BAM files")
    
    for cell_type in cell_types:
        merge_start_time = time.time()
        if chr_contigs:
            final_bam_path = f"{args.outdir}/{bam_name}_{cell_type}_{chrom_str}.bam"
        else:
            final_bam_path = f"{args.outdir}/{bam_name}_{cell_type}_chr{chrom_str}.bam"
        
        logging.info(f"Merging temporary files for cell type: {cell_type} to {final_bam_path}")
        
        # Get list of temp files for this cell type
        temp_files = []
        for region in regions:
            region_str = f"{region[0]}_{region[1]}_{region[2]}"
            if chr_contigs:
                temp_path = f"{temp_dir}/temp_{bam_name}_{cell_type}_{chrom_str}_{region_str}.bam"
            else:
                temp_path = f"{temp_dir}/temp_{bam_name}_{cell_type}_chr{chrom_str}_{region_str}.bam"
            if os.path.exists(temp_path):
                temp_files.append(temp_path)
        
        if temp_files:
            logging.info(f"Found {len(temp_files)} temporary files for cell type {cell_type}")
            # Use pysam to merge the files
            try:
                pysam.merge("-f", final_bam_path, *temp_files)
                merge_end_time = time.time()
                logging.info(f"Successfully merged {len(temp_files)} files for {cell_type} "
                           f"in {merge_end_time - merge_start_time:.2f} seconds")
                
                # Index the merged BAM file
                logging.info(f"Indexing {final_bam_path}")
                pysam.index(final_bam_path)
                
                # Get final file size
                final_size_mb = os.path.getsize(final_bam_path) / (1024 * 1024)
                logging.info(f"Final BAM file for {cell_type}: {final_size_mb:.2f} MB")
                
                # Remove temporary files
                logging.info(f"Removing temporary files for {cell_type}")
                for temp_file in temp_files:
                    os.remove(temp_file)
            except Exception as e:
                logging.error(f"Error merging files for {cell_type}: {e}")
        else:
            logging.warning(f"No temporary files found for cell type {cell_type}")
    
    # Combine summary dictionaries
    combined_summary = {
        "Reads_processed": 0,
        "Valid": 0,
        "No_CB": 0,
        "CB_unmatched": 0,
        "CB_matched": 0
    }
    
    for summary in all_summaries.values():
        for key in combined_summary:
            combined_summary[key] += summary[key]
    
    # Summarize cell type distribution across all regions
    all_cell_types_summary = {}
    for cell_type in cell_types:
        all_cell_types_summary[cell_type] = 0
        
    for region_celltype_data in reads_by_region_celltype.values():
        for cell_type, count in region_celltype_data.items():
            all_cell_types_summary[cell_type] += count
    
    # Write summary
    total_reads = combined_summary['Reads_processed']
    
    logging.info(f"Writing summary results for {total_reads:,} processed reads")
    
    # Get execution time and rate
    end_time = time.time()
    total_elapsed = end_time - start_time
    overall_rate = total_reads / total_elapsed if total_elapsed > 0 else 0
    
    logging.info("===== EXECUTION SUMMARY =====")
    logging.info(f"Total execution time: {total_elapsed:.2f} seconds")
    logging.info(f"Total reads processed: {total_reads:,}")
    logging.info(f"Overall processing rate: {overall_rate:.1f} reads/second")
    logging.info(f"Reads with valid cell barcodes: {combined_summary['Valid']:,} ({combined_summary['Valid']*100/total_reads:.1f}%)")
    logging.info(f"Reads matched to cell types: {combined_summary['CB_matched']:,} ({combined_summary['CB_matched']*100/total_reads:.1f}%)")
    
    # Log cell type distribution
    logging.info("Cell type distribution:")
    for cell_type, count in all_cell_types_summary.items():
        percentage = count * 100 / combined_summary['CB_matched'] if combined_summary['CB_matched'] > 0 else 0
        logging.info(f"  - {cell_type}: {count:,} reads ({percentage:.1f}% of matched reads)")
    
    # Write summary to file
    sum_df = pd.DataFrame([combined_summary])
    sum_df.to_csv(f"{args.outdir}/{bam_name}_{chrom_str}_read_summary.tsv", sep='\t', index=False)
    
    # Write detailed cell type distribution
    cell_type_df = pd.DataFrame([all_cell_types_summary])
    cell_type_df.to_csv(f"{args.outdir}/{bam_name}_{chrom_str}_celltype_distribution.tsv", sep='\t', index=False)
    
    # Clean up temp directory if empty
    try:
        os.rmdir(temp_dir)
        logging.info(f"Removed empty temporary directory: {temp_dir}")
    except:
        logging.info(f"Temporary directory {temp_dir} not empty or cannot be removed")
    
    logging.info("BAM splitting completed successfully")

if __name__ == "__main__":
    main()
