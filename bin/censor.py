#!/usr/bin/env python

# Based on: https://github.com/Alexander-Sol/mzLib/blob/censorMzml/mzLib/Test/FileReadingTests/CensorMzmls.cs
#
# From the paper (https://doi.org/10.1021/acs.jproteome.5c00065):
#
# Data censoring was performed to estimate each PIP algorithm’s native
# peak matching error rate. This involved comparing peak traces that
# were unambiguously associated with MS2-detected peptides with
# corresponding peak traces identified by PIP using the same donor
# peptide but with the detecting MS2 spectra removed from the acceptor
# run. Data censoring was performed separately for each software tool
# that was evaluated as follows. Let MS2 Peptidesi be the list of all
# peptides that were detected in the pure-sample run i. MS2 Peptidesi
# was filtered according to the following criteria:
#
# 1. The peptide must have been detected at an FDR < 0.1%.
#
# 2. The peptide sequence must be derived from the H.  sapiens
#    proteome and must not occur in any other proteome used during the
#    search.
#
# 3. The peptide must have been detected in at least one other run in
#    the data set with a score greater than the score of the peptide
#    detection in run i.
#
# 4. The peptide must only be associated with one MS2 spectrum in run i.
#
# Then, Censored peptidesi was created by randomly selecting 500
# peptides from MS2 Peptidesi for censoring. Censoring was performed
# by editing the data file from run i. The MS2 spectrum associated
# with a peptide in Censored peptidesi was replaced with an MS2
# spectrum that contained only one peak at 150 m/z. This was performed
# for every data file associated with a pure-sample run, and then the
# resulting censored data files were searched and quantified alongside
# the original mixedsample data files. A different set of censored
# data files was generated for each software.
import argparse
import csv
import functools
import logging
import random
from pathlib import Path

import pyopenms as oms
from pipe import chain, groupby, select, take, where


################################################################################
class Censor:
    def __init__(self, peptide_file, psm_file, mzmls):
        self.peptide_file = peptide_file
        self.psm_file = psm_file
        self.output_dir = Path(".")

        self.mzml_files = {mzml.stem: mzml for mzml in mzmls}
        self.log = logging.getLogger(__name__)

        self.min_peptides = 100
        self.min_mzml_files = 2
        self.min_peptide_qval = 0.001
        self.min_psm_qval = 0.01
        self.max_psms_per_file = 500

        self.col_rt = "Scan Retention Time"
        self.col_seq = "Full Sequence"
        self.col_file = "File Name"

    def run(self):
        self._load_peptides()
        self._load_psms()
        self._restrict_psms()
        self._censor()

    def _load_peptides(self):
        self.log.info("reading peptides")
        self.parsed_peptides = self._read_psmtsv_file(
            self.peptide_file, self.min_peptide_qval
        )
        self.all_peptide_seqs = {row[self.col_seq] for row in self.parsed_peptides}
        self.peptides_by_file = dict(
            self.parsed_peptides | groupby(lambda p: p[self.col_file])
        )

        self.peptide_seq_by_file = {}

        for file, peptides in self.peptides_by_file.items():
            seq = {p[self.col_seq] for p in peptides}
            self.peptide_seq_by_file[file] = seq

        self.log.info(f"peptides that passed filtering: {len(self.parsed_peptides)}")

        assert len(self.parsed_peptides) >= self.min_peptides, (
            f"{len(self.parsed_peptides)} is not enough peptides"
        )

        self.log.info(
            f"number of mzML files in peptide list: {len(self.peptides_by_file)}"
        )

        assert len(self.peptides_by_file) >= self.min_mzml_files, (
            f"{len(self.peptides_by_file)} is not enough mzML files"
        )

    def _load_psms(self):
        self.log.info("reading PSMs")
        self.parsed_psms = self._read_psmtsv_file(self.psm_file, self.min_psm_qval)
        self.log.info(f"number of PSMs that passed filtering: {len(self.parsed_psms)}")

        self.psms_by_file = dict(
            self.parsed_psms
            | groupby(lambda psm: psm[self.col_file])
            | select(lambda group: (group[0], list(group[1])))
        )
        self.log.info(f"number of mzML files in PSMs list {len(self.psms_by_file)}")

        assert len(self.psms_by_file) >= self.min_mzml_files, (
            f"{len(self.psms_by_file)} is not enough mzML files"
        )

    def _restrict_psms(self):
        """
        Restrict PSMs that eluded more than once in a single run.
        """
        self.log.info("checking for multi-eluding PSMs")
        self.forbidden_sequences = set()

        def get_rt(psm):
            if self.col_rt in psm:
                return float(psm[self.col_rt])
            else:
                return 0.0

        for psms in self.psms_by_file.values():
            for sequence, group in psms | groupby(lambda psm: psm[self.col_seq]):
                group = list(group | select(get_rt))

                if len(group) > 1:
                    rt_min = min(group)
                    rt_max = max(group)

                    if (rt_max - rt_min) > 3:
                        self.forbidden_sequences.add(sequence)

        self.log.info(
            f"{len(self.forbidden_sequences)} PSMs forbidden from being censored"
        )

    def _read_psmtsv_file(self, file_name: str, qval: float):
        rows = []

        with open(file_name, "r") as file:
            reader = csv.DictReader(file, delimiter="\t")

            for row in reader:
                if (
                    float(row["PEP_QValue"]) <= qval
                    and "xidation" not in row[self.col_seq]
                ):
                    rows.append(row)

        return rows

    def _censor(self):
        self.log.info("censoring PSMs")
        censored_psms = []

        for file, psms in self.psms_by_file.items():
            self.log.info(f"censoring {len(psms)} PSMs in file {file}")
            assert len(psms) > 0, f"no PSMs in file {file}?"

            if "HeYe" in file:
                self.log.info(f"skipping HeYe file {file}")
                path = self.mzml_files[file]
                link = self._censored_mzml_name(path)
                link.symlink_to(path)
                continue

            def not_top_peptide(file, psm):
                # If this file contains the top peptide instance, we can't remove it
                return (
                    not (psm[self.col_seq] in self.peptide_seq_by_file[file])
                    and "Homo sapiens" in psm["Organism Name"]
                    and psm["Decoy/Contaminant/Target"] == "T"
                )

            def allowed_sequence(seq, psms):
                return (
                    seq in self.all_peptide_seqs
                    and len(psms) == 1
                    and (not (seq in self.forbidden_sequences))
                )

            to_remove = list(
                psms
                | where(functools.partial(not_top_peptide, file))
                | groupby(lambda psm: psm[self.col_seq])
                | select(lambda group: (group[0], list(group[1])))
                | where(lambda group: allowed_sequence(group[0], group[1]))
                | select(lambda group: group[1])
                | chain
                | groupby(lambda psm: round(float(psm["Precursor Mass"]), 2))
                | select(lambda group: next(iter(group[1])))
            )

            random.shuffle(to_remove)
            to_remove = list(to_remove | take(self.max_psms_per_file))
            assert len(to_remove) > 0, f"no censored PSMs from file {file}"

            censored_psms.extend(to_remove)
            self._remove_psms(self.mzml_files[file], to_remove)

        self.log.info(f"censored {len(censored_psms)} PSMs")
        assert len(censored_psms) > 0
        self._record_censored_psms(censored_psms)

    def _censored_mzml_name(self, path):
        return Path(Path(path).stem + "-censored.mzML")

    def _remove_psms(self, file, psms):
        base = self._censored_mzml_name(file)
        output_path = self.output_dir / base
        self.log.info(f"generating censored file file: {output_path}")

        exp = oms.MSExperiment()
        oms.MzMLFile().load(str(file), exp)
        spectra = exp.getSpectra()
        for psm in psms:
            scan_number = int(psm["Scan Number"])
            spectra[scan_number - 1].set_peaks([150.0], [10000.0])

        exp.setSpectra(spectra)
        oms.MzMLFile().store(str(output_path), exp)

    def _record_censored_psms(self, psms):
        file = self.output_dir / "censored.psmtsv"

        with open(file, "w") as file:
            writer = csv.DictWriter(file, fieldnames=psms[0].keys(), delimiter="\t")
            writer.writeheader()
            writer.writerows(psms)


################################################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="censor.py",
        description="Censor mzML files for FlashLFQ",
    )

    # https://docs.python.org/3/library/argparse.html
    parser.add_argument("--peptides")
    parser.add_argument("--psms")
    parser.add_argument("--mzml", action="extend", nargs="+", type=Path)

    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG)

    censor = Censor(args.peptides, args.psms, args.mzml)
    censor.min_peptide_qval = 0.01  # FIXME
    censor.run()
