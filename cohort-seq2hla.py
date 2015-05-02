
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

import logging
import argparse
from os.path import join, exists, isdir
from os import listdir
from glob import glob

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

args = parser.parse_args()

def paths_from_pattern(base_dir, pattern):
    full_pattern = join(base_dir, pattern)
    paths = [
        join(base_dir, filename)
        for filename in glob(full_pattern)
    ]
    if len(paths) == 0:
        raise ValueError("No FASTQ files found for %s" % full_pattern)
    else:
        for path in paths:
            logging.info("Found FASTQ file %s" % path)
    return paths

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
        fastq_path = join(sample_path, args.rna_fastq_subdir)
        if not exists(fastq_path):
            raise ValueError("Missing FASTQ subdirectory '%s' for %s:%s" % (
                fastq_path,
                group_name,
                sample_name))

        left_fastq_paths = paths_from_pattern(
            fastq_path,
            args.left_reads_pattern)

        right_fastq_paths = paths_from_pattern(
            fastq_path,
            args.right_reads_pattern)
