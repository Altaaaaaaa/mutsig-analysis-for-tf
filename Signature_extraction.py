import os
import argparse
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerMatrixGenerator import install as genInstall

parser = argparse.ArgumentParser(description='sigprofiler_matgen')

parser.add_argument('--ref_genome', required=True, help='select reference genome(e.g. GRCh37)')
parser.add_argument('--min', type = int, required=True, help='enter minimum_signatures')
parser.add_argument('--max', type = int, required=True, help='enter maximum_signatures')
parser.add_argument('--input_dir', required=True, help='directory name of input data')
parser.add_argumnet('--output_dir', required=True, help = 'desired directory name of output data')
parser.add_argument('--threads', type = int, default = -1, help = 'number of threads to use')

args = parser.parse_args()


genInstall.install(str(args.ref_genome), rsync=False, bash=True)

# Matrix generation
matGen.SigProfilerMatrixGeneratorFunc('DATA', str(args.ref_genome), str(args.input_dir), plot=False,
                                                     exome=False, bed_file=None, chrom_based=False, tsb_stat=False,
                                                     seqInfo=False, cushion=100)

# Signature extraction
path_to_example_table = sig.importdata("matrix")
data = path_to_example_table # you can put the path to your tab delimited file containing the mutational catalog matrix/table
sig.sigProfilerExtractor("matrix", args.output_dir,
                     f"{args.input_dir}/output/SBS/DATA.SBS96.all", reference_genome=str(args.ref_genome),
                         minimum_signatures=args.min, maximum_signatures=args.max, nmf_replicates=100, cpu=args.threads)
