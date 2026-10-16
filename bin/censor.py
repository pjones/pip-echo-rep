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
import re
from pathlib import Path

import pyopenms as oms
from pipe import chain, groupby, select, take, traverse, where


################################################################################
class Ident:
    """
    A peptide identification with some useful attributes.
    """

    @classmethod
    def from_peptide_identification(
        cls, pep: oms.PeptideIdentification, file_ref=None
    ) -> "Ident":
        """
        Create an Ident object from an OpenMS `PeptideIdentification`.
        """
        assert len(pep.getHits()) > 0
        assert pep.getScoreType() == "q-value"

        hit = pep.getHits()[0]
        aaseq = hit.getSequence()
        ident = Ident()

        ident.file_ref = file_ref
        ident.rt = pep.getRT()  # Seconds
        ident.sequence = aaseq.toUnmodifiedString()
        ident.sequence_with_mods = str(aaseq)
        ident.score = hit.getScore()
        ident.spectrum_ref = pep.getSpectrumReference()
        ident.scan = None  # Must use `spectrum_ref`
        ident.is_target = not hit.isDecoy()

        ident.is_human = any(
            hit.extractProteinAccessionsSet()
            | select(lambda s: re.search(r"human", s, re.IGNORECASE))
        )

        return ident

    def to_row(self):
        return {
            "File": self.file_ref,
            "Full Sequence": self.sequence_with_mods,
            "Base Sequence": self.sequence,
            "Score": self.score,
            "Is Target": self.is_target,
            "Is Human": self.is_human,
        }

    def sort_key(self) -> float:
        """
        Return the Q-value if this identification is not a decoy.
        """
        if self.is_target:
            return self.score
        else:
            return 1.0

    def is_valid(self, max_score: float) -> bool:
        """Can this identification be used for censoring?"""
        return (
            self.score <= max_score
            and self.is_human
            #
            # The paper doesn't mention this but the FlashLFQ censor
            # script includes this:
            and not "xidation" in self.sequence_with_mods
        )


################################################################################
class IdFile:
    peptides: list[Ident]

    def __init__(self, file: Path):
        self.id_file = file

    def name(self) -> str:
        return str(self.id_file.stem)

    def mzml_file_path(self) -> Path:
        raise RuntimeError("not implemented")

    def duplicates(self, max_seconds: float) -> list[str]:
        """
        Return a list of peptide sequences that eluded more than once.
        """

        def longer_than_max_seconds(group: list[Ident]):
            rts = list(group | select(lambda p: p.rt))

            min_rt = min(rts)
            max_rt = max(rts)

            return abs(min_rt - max_rt) > max_seconds

        return list(
            self.peptides
            | groupby(lambda pep: pep.sequence)
            | select(lambda group: [group[0], list(group[1])])
            | where(lambda group: len(group[1]) > 1)
            | where(lambda group: longer_than_max_seconds(group[1]))
            | select(lambda group: group[0])
        )


################################################################################
class Censor:
    def __init__(self):
        self.output_dir = Path(".")
        self.log = logging.getLogger(__name__)

        # When a peptide takes longer than this number of seconds to
        # elude it is considered to have eluded more than once in the
        # same run.
        self.max_rt_seconds = 180  # 3 minutes

        self.min_peptides = 100
        self.min_mzml_files = 2
        self.max_peptide_qval = 0.001
        self.max_psm_qval = 0.01
        self.max_psms_per_file = 500

    def censor(self, files: dict[str, IdFile], peptides: dict[str, list[Ident]]):
        peptides = self.remove_forbidden_sequences(files, peptides)
        self.log.info(f"reviewing {len(peptides)} sequences for censoring")

        to_remove_by_file = dict(
            peptides.values()
            | select(self.remove_best_ident)
            | where(lambda ids: len(ids) > 1)
            | chain
            | groupby(lambda ident: ident.file_ref)
            | select(lambda group: [group[0], list(group[1])])
            | select(lambda group: [group[0], (random.shuffle(group[1]), group[1])[1]])
            | select(
                lambda group: [group[0], list(group[1] | take(self.max_psms_per_file))]
            )
        )

        for name, file in files.items():
            mzml_file = file.mzml_file_path()

            if "HeYe" in str(mzml_file):
                self.log.info(f"skipping HeYe file {name}")
                link = self.file_to_censored_path(mzml_file)
                link.symlink_to(mzml_file)
                continue

            self.log.info(
                f"replacing {len(to_remove_by_file[name])} spectra from {name}"
            )
            self.replace_spectra(mzml_file, to_remove_by_file[name])

        self.dump_censored(
            list(to_remove_by_file.values() | chain | select(Ident.to_row))
        )

    def remove_best_ident(self, ids: list[Ident]) -> list[Ident]:
        """
        Given a list of identified peptides, return another list with
        the best scoring peptide removed.

        They are also grouped by their source file and only one
        peptide for each file is returned (the best scoring one).
        """
        by_score = sorted(
            (
                ids
                | groupby(lambda id: id.file_ref)
                | select(lambda group: sorted(group[1], key=Ident.sort_key))
                | select(lambda ids: ids[0])
            ),
            key=Ident.sort_key,
        )
        return by_score[1:]

    def remove_forbidden_sequences(
        self, files: dict[str, IdFile], peptides: dict[str, list[Ident]]
    ) -> dict[str, list[Ident]]:
        """
        Remove identifications from `peptides` that are listed as
        forbidden in the paper or the original FlashLFQ censoring
        script.
        """
        #
        # Peptides that eluded more than once in a single run:
        #
        multi_eluders = set(
            files.values()
            | select(lambda file: file.duplicates(self.max_rt_seconds))
            | chain
        )
        self.log.info(f"number of multi-eluding sequences: {len(multi_eluders)}")

        #
        # Peptides that only eluded in a single run:
        #
        only_one_file = set(
            peptides.values()
            | where(lambda peps: len(peps) <= 1)
            | select(lambda peps: peps[0].sequence)
            | chain
        )
        self.log.info(f"number of sequences from single runs: {len(only_one_file)}")

        forbidden_sequences = multi_eluders | only_one_file
        self.log.info(f"number of forbidden sequences: {len(forbidden_sequences)}")

        return dict(
            list(peptides.items())
            | where(lambda group: group[0] not in forbidden_sequences)
        )

    def file_to_censored_path(self, file: Path) -> Path:
        base = Path(file).stem + "-censored.mzML"
        return self.output_dir / base

    def replace_spectra(self, file: Path, idents: list[Ident]):
        """
        Censor the given mzML file by replacing some of the spectra
        with a dummy spectrum that only contains a single peak.
        """
        output_path = self.file_to_censored_path(file)
        self.log.info(f"generating censored file file: {output_path}")

        exp = oms.MSExperiment()
        oms.MzMLFile().load(str(file), exp)

        lookup = oms.SpectrumLookup()
        lookup.readSpectra(exp)

        spectra = exp.getSpectra()

        for ident in idents:
            index = None

            if hasattr(ident, "spectrum_ref"):
                index = lookup.findByNativeID(ident.spectrum_ref)
            elif hasattr(ident, "scan"):
                index = ident.scan - 1

            spectra[index].set_peaks([150.0], [10000.0])

        exp.setSpectra(spectra)
        oms.MzMLFile().store(str(output_path), exp)

    def dump_censored(self, ids: list[Ident]):
        with open(self.output_dir / "censored.tsv", "w") as file:
            writer = csv.DictWriter(file, fieldnames=ids[0].keys(), delimiter="\t")
            writer.writeheader()
            writer.writerows(ids)


################################################################################
class ProteomicsLFQ(Censor):
    class IdXmlFile(IdFile):
        def __init__(self, file: Path):
            IdFile.__init__(self, file)
            self.mzml_file = file.parent / (file.stem + ".mzML")

        def mzml_file_path(self) -> Path:
            return self.mzml_file

        def load(self, max_score: float):
            _, peptide_ids = oms.IdXMLFile().load(str(self.id_file))
            to_ident = lambda pep: Ident.from_peptide_identification(pep, self.name())

            self.peptides = list(
                peptide_ids
                | select(to_ident)
                | where(lambda ident: ident.is_valid(max_score))
            )

            return self

    def __init__(self, input_dir: Path):
        Censor.__init__(self)
        self.input_dir = input_dir

    def run(self):
        self.log.info("censoring spectra for ProteomicsLFQ")
        self._load_data()
        return self.censor(self.id_files, self.peptides)

    def _load_data(self):
        self.id_files = dict(
            self.input_dir.glob("*.idXML")
            | select(lambda path: ProteomicsLFQ.IdXmlFile(path))
            | select(lambda file: [file.name(), file.load(self.max_psm_qval)])
        )

        self.log.info(f"number of identity files: {len(self.id_files)}")
        assert len(self.id_files) >= self.min_mzml_files

        self.peptides = dict(
            self.id_files.values()
            | select(lambda file: file.peptides)
            | traverse
            | groupby(lambda ident: ident.sequence)
            | select(lambda group: [group[0], list(group[1])])
        )

        self.log.info(f"number of read peptides: {len(self.peptides)}")
        assert len(self.peptides) >= self.min_peptides


################################################################################
class FlashLFQ(Censor):
    def __init__(self, peptide_file, psm_file, mzmls):
        Censor.__init__(self)

        self.peptide_file = peptide_file
        self.psm_file = psm_file
        self.mzml_files = {mzml.stem: mzml for mzml in mzmls}

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
            self.peptide_file, self.max_peptide_qval
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
        self.parsed_psms = self._read_psmtsv_file(self.psm_file, self.max_psm_qval)
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
    logopts = {
        "format": "[%(asctime)s %(levelname)s]: %(message)s",
        "datefmt": "%Y-%m-%dT%H:%M:%S",
        "level": logging.DEBUG,
    }

    parser = argparse.ArgumentParser(
        prog="censor.py",
        description="Censor mzML files",
    )

    # https://docs.python.org/3/library/argparse.html
    parser.add_argument("--tool", required=True, choices=["flashlfq", "proteomicslfq"])
    parser.add_argument("files", action="extend", nargs="+", type=Path)

    args = parser.parse_args()
    logging.basicConfig(**logopts)

    censor = None

    if args.tool == "flashlfq":
        censor = FlashLFQ(args.files[0], args.files[1], args.files[2:])
        censor.max_peptide_qval = 0.01  # FIXME
    elif args.tool == "proteomicslfq":
        censor = ProteomicsLFQ(*args.files)

    censor.run()
