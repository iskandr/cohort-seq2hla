
# Copyright (c) 2015. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
from glob import glob
import logging
from os.path import join, exists, isdir
from os import listdir
import subprocess

parser = argparse.ArgumentParser()

parser.add_argument(
    "--base-dir",
    default="/hdfs/datasets/martignetti_ovarian/",
    help="Directory which contains DNA & RNA samples for sample groups")

parser.add_argument(
    "--rna-subdir",
    default=join(
        "Illumina_RNA",
        "QC_S106.B202_PCT189_RNA.PE.RNASeqRibozero.RAPiD.Human"),
    help="Subdirectory which contains all RNA samples in a group")

parser.add_argument(
    "--rna-fastq-subdir",
    default="Raw/RNA.IlluminaHiSeq2500.RiboZero/",
    help="Within a particular RNA sample's dir, where are the FASTQ files")

parser.add_argument(
    "--left-reads-pattern",
    default="*_R1_*.fastq.gz",
    help="Wildcard pattern for FASTQ files containing left reads of mate-pairs")

parser.add_argument(
    "--right-reads-pattern",
    default="*_R2_*.fastq.gz",
    help="Wildcard pattern for FASTQ files containing left reads of mate-pairs")

parser.add_argument(
    "--seq2hla-path",
    default="../seq2hla/seq2HLA.py",
    help="Where is the seq2HLA repository located?")

parser.add_argument(
    "--temp-fastq-dir",
    default=".",
    help="Dir for concatenated FASTQ's (seq2HLA can't handle multiple lanes)")

args = parser.parse_args()

def paths_from_pattern(base_dir, pattern):
    full_pattern = join(base_dir, pattern)
    paths = [
        join(base_dir, filename)
        for filename in glob(full_pattern)
    ]
    if len(paths) == 0:
        raise ValueError("No FASTQ files found for %s" % full_pattern)
    return paths

def concat_compressed_fastq_files(
            group_name,
            sample_name,
            suffix,
            fastq_paths,
            target_dir):
        combined_fastq_name = "%s_%s_%s.fastq" % (
            group_name,
            sample_name,
            suffix)
        combined_fastq_path = join(target_dir, combined_fastq_name)
        logging.info(
            "Combining %d FASTQ files into %s",
            len(fastq_paths),
            combined_fastq_path)
        with open(combined_fastq_path, 'w') as output_file:
            subprocess.check_call(
                ["zcat"] + fastq_paths,
                stdout=output_file)
        assert exists(combined_fastq_path)
        return combined_fastq_path

if __name__ == "__main__":
    if not exists(args.base_dir):
        raise ValueError("Directory '%s' does not exist" % args.base_dir)
    # dictionary mapping sample group names to directory paths
    group_paths = {}
    for group_name in listdir(args.base_dir):
        group_path = join(args.base_dir, group_name)
        if isdir(group_path):
            logging.info("Sample group %s => %s", group_name, group_path)
            group_paths[group_name] = group_path

    # dictionary mapping (group_name, sample_name) pairs full paths
    sample_paths = {}
    for (group_name, group_path) in group_paths.items():
        rna_path = join(group_path, args.rna_subdir)
        if not exists(rna_path):
            raise ValueError(
                "Missing RNA subdirectory for sample group %s, expected %s" % (
                    group_name, rna_path))

        for sample_name in listdir(rna_path):
            sample_path = join(rna_path, sample_name)
            if isdir(sample_path):
                logging.info("Sample %s:%s => %s",
                    group_name,
                    sample_name,
                    sample_path)
                sample_paths[(group_name, sample_name)] = sample_path

    for ((group_name, sample_name), sample_path) in sample_paths.items():
        logging.info("Looking for FASTQ files for %s:%s" % (
            group_name, sample_name))
        fastq_path = join(sample_path, args.rna_fastq_subdir)
        if not exists(fastq_path):
            raise ValueError("Missing FASTQ subdirectory '%s' for %s:%s" % (
                fastq_path,
                group_name,
                sample_name))

        left_fastq_paths = paths_from_pattern(
            fastq_path,
            args.left_reads_pattern)

        left_combined_fastq_path = concat_compressed_fastq_files(
            group_name=group_name,
            sample_name=sample_name,
            suffix="R1",
            fastq_paths=left_fastq_paths,
            target_dir=args.temp_fastq_dir)

        right_fastq_paths = paths_from_pattern(
            fastq_path,
            args.right_reads_pattern)

        right_combined_fastq_path = concat_compressed_fastq_files(
            group_name=group_name,
            sample_name=sample_name,
            suffix="R2",
            fastq_paths=right_fastq_paths,
            target_dir=args.temp_fastq_dir)

        subprocess.check_call(
            [
                "python", args.seq2hla_path,
                "-1", left_combined_fastq_path,
                "-2", right_combined_fastq_path,
                "-r", "%s_%s" % (group_name, sample_name)
            ])