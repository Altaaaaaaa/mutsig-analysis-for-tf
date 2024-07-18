import os
import argparse
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerMatrixGenerator import install as genInstall

# 인자값을 받을 수 있는 인스턴스 생성
parser = argparse.ArgumentParser(description='sigprofiler_matgen')

# 입력받을 인자값 등록
parser.add_argument('--ref_genome', required=True, help='select reference genome(e.g. GRCh37)')
parser.add_argument('--input_dir', required=True, help='directory name of input data')

# 입력받은 인자값을 args에 저장 (type: namespace)
args = parser.parse_args()




















genInstall.install(str(args.ref_genome), rsync=False, bash=True)

matGen.SigProfilerMatrixGeneratorFunc('DATA', str(args.ref_genome), str(args.input_dir), plot=False,
                                                     exome=False, bed_file=None, chrom_based=False, tsb_stat=False,
                                                     seqInfo=False, cushion=100)