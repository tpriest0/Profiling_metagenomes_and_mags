import argparse
import pyrodigal_gv
from pysam.libcfaidx import FastxFile
from multiprocessing.pool import Pool
import gzip
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')

logging.info("Starting gene prediction")

parser = argparse.ArgumentParser()
parser.add_argument("-i","--in_fasta", type=str, required=True, help="Fasta file containing genome/element sequences")
parser.add_argument("-a","--out_genes_aa", type=str, required=True, help="Path and name of output that will contain gene predictions in amino acid format.")
parser.add_argument("-n","--out_genes_nc", type=str, required=True, help="Path and name of output that will contain gene predictions in nucleotide format.")
parser.add_argument("-g","--out_genes_gff", type=str, required=True, help="Path and name of output that that will contain gene predictions in GFF format.")
parser.add_argument("-s","--out_gene_scores", type=str, required=True, help="Path and name of output that will contain information on the predicted genes. The file will be a tab-delimited table.")
parser.add_argument("-t","--threads", type=int, required=False)
parser.add_argument("-c", "--chunk_size", type=int, default=10000, help="Number of sequences to process per chunk")

args = parser.parse_args()

orf_finder = pyrodigal_gv.ViralGeneFinder(meta=True)

# Function to process a chunk of sequences
def process_chunk(entries_chunk):
    results = []
    for entry in entries_chunk:
        result = orf_finder.find_genes(entry.sequence)
        results.append((entry, result))
    return results

chunk_count = 0

with open(args.out_genes_aa, "w") as genes_aa_file, open(args.out_genes_nc, "w") as genes_nc_file, open(args.out_gene_scores, "w") as stats_file, open(args.out_genes_gff, "w") as gff_file:
    stats_file.write("Contig\tGene\tStart\tStop\tGene_length\tGC_content\tRBS_motif\tScore\tStrand\tTranslation_table\n")

    # Process the FASTA file in chunks to minimize memory usage
    with FastxFile(args.in_fasta) as fasta:
        chunk = []
        for entry in fasta:
            chunk.append(entry)
            if len(chunk) == args.chunk_size:
                # Process the chunk
                with Pool(args.threads) as pool:
                    results = pool.map(process_chunk, [chunk])

                # Write the results to output files
                for chunk_results in results:
                    for entry, result in chunk_results:
                        for i, gene in enumerate(result):
                            gene_name = f"{entry.name}_{i+1}"

                            # Write nucleotide gene sequences
                            genes_nc_file.write(f">{gene_name}\n")
                            genes_nc_file.write(gene.sequence() + "\n")

                            # Write amino acid gene sequences
                            genes_aa_file.write(f">{gene_name}\n")
                            genes_aa_file.write(gene.translate() + "\n")

                            # Write GFF entry for the gene
                            gff_file.write(f"{entry.name}\tpyrodigal\tCDS\t{gene.begin}\t{gene.end}\t{gene.score:.2f}\t{gene.strand}\t0\tID={gene_name};Name={gene_name}\n")

                            # Write gene stats
                            gene_length = gene.end - gene.begin + 1
                            stats_file.write(f"{entry.name}\t{gene_name}\t{gene.begin}\t{gene.end}\t{gene_length}\t"
                                             f"{gene.gc_cont:.2f}\t{gene.rbs_motif}\t{gene.score:.2f}\t{gene.strand}\t{gene.translation_table}\n")
                
                chunk_count += 1
                print(f"Processed chunk {chunk_count} containing {args.chunk_size} sequences.")

                chunk = []

        # Process any remaining sequences
        if chunk:
            with Pool(args.threads) as pool:
                results = pool.map(process_chunk, [chunk])
            for chunk_results in results:
                for entry, result in chunk_results:
                    for i, gene in enumerate(result):
                        gene_name = f"{entry.name}_{i+1}"

                        # Write nucleotide gene sequences
                        genes_nc_file.write(f">{gene_name}\n")
                        genes_nc_file.write(gene.sequence() + "\n")

                        # Write amino acid gene sequences
                        genes_aa_file.write(f">{gene_name}\n")
                        genes_aa_file.write(gene.translate() + "\n")

                        # Write GFF entry for the gene
                        gff_file.write(f"{entry.name}\tpyrodigal\tCDS\t{gene.begin}\t{gene.end}\t{gene.score:.2f}\t{gene.strand}\t0\tID={gene_name};Name={gene_name}\n")

                        # Write gene stats
                        gene_length = gene.end - gene.begin + 1
                        stats_file.write(f"{entry.name}\t{gene_name}\t{gene.begin}\t{gene.end}\t{gene_length}\t"
                                         f"{gene.gc_cont:.2f}\t{gene.rbs_motif}\t{gene.score:.2f}\t{gene.strand}\t{gene.translation_table}\n")
