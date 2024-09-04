#!/usr/bin/env python3

import pysam
import argparse
import os

def main(args):

        inputfile = args.input

        bam_file = pysam.AlignmentFile(inputfile, "rb")

        with open(args.output, mode = 'w') as f:
                for read in bam_file.fetch():
                #read_text = f"{first_read.query_name},{first_read.reference_name},{first_read.reference_start},{first_read.mapping_quality},{first_read.cigarstring},{first_read.flag},{first_read.query_sequence},{','.join(map(str, first_read.query_qualities))}"
                        read_text = read
                        f.write("{}\t{}\t{}\t{}\t{}\n".format(read.reference_name, read.reference_start,read.flag, read.mapping_quality,read.cigarstring))

        bam_file.close()

        folder= "/".join(args.output.split("/")[:-1])
        tar_command = "tar -C {} -czf {}.tar.gz {}".format(folder, args.output,args.output.split("/")[-1])
        os.system(tar_command)
        # Optional: Remove the original file if only the zipped version is needed
        os.remove(args.output)


def run():
        """
                Parse arguments and run
        """
        parser = argparse.ArgumentParser(description="program determine pangenome-allele types")
        parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
        parser.add_argument("-o", "--output", help="path to output file", dest="output",type=str, required=True)
        parser.set_defaults(func=main)
        args = parser.parse_args()
        args.func(args)


if __name__ == "__main__":
        run()
