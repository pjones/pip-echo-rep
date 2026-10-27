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
# that was evaluated as follows. Let MS2 Peptides_i be the list of all
# peptides that were detected in the pure-sample run i. MS2 Peptides_i
# was filtered according to the following criteria:
#
# 1. The peptide must have been detected at an FDR < 0.1%.
#
#     (Checked with the peptide is loaded from an identification file via
#      `Ident.is_valid`.)
#
# 2. The peptide sequence must be derived from the H.  sapiens
#    proteome and must not occur in any other proteome used during the
#    search.
#
#     (Checked with the peptide is loaded from an identification file via
#      `Ident.is_valid`.)
#
# 3. The peptide must have been detected in at least one other run in
#    the data set with a score greater than the score of the peptide
#    detection in run i.
#
#    (Checked in the `Censor.remove_forbidden_sequences` function
#     which rejects peptides identified in only one file.  Then later
#     the `Censor.remove_best_ident` function ensures the best scoring
#     run is removed from the list of peptides that can be censored.)
#
# 4. The peptide must only be associated with one MS2 spectrum in run i.
#
#    (Also verified in `Censor.remove_forbidden_sequences`.)
#
# Then, Censored peptides_i was created by randomly selecting 500
# peptides from MS2 Peptides_i for censoring.  Censoring was performed
# by editing the data file from run i. The MS2 spectrum associated
# with a peptide in Censored peptides_i was replaced with an MS2
# spectrum that contained only one peak at 150 m/z.  This was
# performed for every data file associated with a [human] sample run,
# and then the resulting censored data files were searched and
# quantified alongside the original mixedsample data files.  A
# different set of censored data files was generated for each
# software.
import argparse
import csv
import logging
import random
import re
from operator import attrgetter
from pathlib import Path

import pyopenms as oms
from pipe import chain, groupby, select, take, traverse, where


################################################################################
class Ident:
    """
    A peptide identification with some useful attributes.
    """

    def to_row(self) -> dict[str, str]:
        """
        Return a dictionary for inserting into the manifest file.
        """
        raise RuntimeError("not implemented")

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
        self.peptides = []

    def name(self) -> str:
        return str(self.id_file.stem)

    def mzml_file_path(self) -> Path:
        raise RuntimeError("not implemented")

    def duplicates(self, max_minutes: float) -> list[str]:
        """
        Return a list of peptide sequences that eluded more than once.
        """

        def longer_than_max_minutes(group: list[Ident]):
            rts = list(group | select(lambda p: p.rt))

            min_rt = min(rts)
            max_rt = max(rts)

            return abs(min_rt - max_rt) > max_minutes

        return list(
            self.peptides
            | groupby(lambda pep: pep.sequence)
            | select(lambda group: [group[0], list(group[1])])
            | where(lambda group: len(group[1]) > 1)
            | where(lambda group: longer_than_max_minutes(group[1]))
            | select(lambda group: group[0])
        )


################################################################################
class Censor:
    def __init__(self):
        self.output_dir = Path(".")
        self.manifest_file = "censored.tsv"
        self.log = logging.getLogger(__name__)

        # When a peptide takes longer than this number of minutes to
        # elude it is considered to have eluded more than once in the
        # same run.
        self.max_rt_minutes = 3.0

        self.min_peptides = 100
        self.min_mzml_files = 2
        self.max_peptide_qval = 0.001
        self.max_psm_qval = 0.01
        self.max_psms_per_file = 500

    def censor(self, files: dict[str, IdFile]):
        self.log.info(f"number of identity files: {len(files)}")
        assert len(files) >= self.min_mzml_files

        peptides = dict(
            files.values()
            | select(lambda file: file.peptides)
            | traverse
            | groupby(lambda ident: ident.sequence_with_mods)
            | select(lambda group: [group[0], list(group[1])])
            | where(lambda group: len(group[1]) > 1)
        )

        self.log.info(f"number of unique peptide sequences: {len(peptides)}")
        assert len(peptides) >= self.min_peptides

        peptides = self.remove_forbidden_sequences(files, peptides)
        self.log.info(f"reviewing {len(peptides)} sequences for censoring")

        to_remove_by_file = self.filter_peptides(peptides)
        self.log.info(f"censoring {len(to_remove_by_file)} files")

        for name, file in files.items():
            mzml_file = file.mzml_file_path()

            if "HeYe" in str(mzml_file):
                self.log.info(f"skipping HeYe file {name}")
                link = self.file_to_censored_path(mzml_file)
                link.symlink_to(mzml_file)
                continue

            if name not in to_remove_by_file:
                self.log.warning(f"{name} has not peptides to censor!")
                continue

            self.log.info(
                f"replacing {len(to_remove_by_file[name])} peptides from {name}"
            )
            self.replace_spectra(mzml_file, to_remove_by_file[name])

        self.dump_censored(
            list(to_remove_by_file.values() | chain | select(lambda id: id.to_row()))
        )

    def filter_peptides(self, peptides: dict[str, list[Ident]]):
        def best(group: tuple[float, list[Ident]]) -> Ident:
            """
            Return the best peptide from `[mass, [peptides]]`.
            """
            by_key = sorted(group[1], key=lambda id: id.sort_key())
            by_rt_diff = sorted(by_key, key=attrgetter("rt_diff"))
            return by_rt_diff[0]

        def by_precursor(ids: [Ident]) -> [Ident]:
            """
            `ids` should be a list of peptides from the same file.
            Limits the results to one peptide per precursor.
            """
            return (
                ids  # [id]
                | groupby(lambda id: round(id.precursor_mass, 2))  # [[mass, [id]]]
                | select(best)  # [id]
            )

        def limit(ids: [Ident]) -> [Ident]:
            return list(ids | take(self.max_psms_per_file))

        return dict(
            peptides.values()
            | select(self.remove_best_ident)
            | where(lambda ids: len(ids) > 0)
            | traverse
            | groupby(lambda ident: ident.file_ref)
            | select(lambda group: [group[0], list(by_precursor(group[1]))])
            | select(lambda group: [group[0], (random.shuffle(group[1]), group[1])[1]])
            | select(lambda group: [group[0], limit(group[1])])
        )

    def remove_best_ident(self, ids: list[Ident]) -> list[Ident]:
        """
        Given a list of peptides that all have the same amino acid
        sequence, return list of peptides from the original list that
        can be censored.

        The peptides will be grouped by the mzML file from which they
        originated and the best scoring peptide for each group will be
        removed from the return list to ensure it can be used by the
        quantification process in MBR.

        Returned peptides will have a new attribute (`rt_diff`) that
        is the difference in retention time between it and the best
        scoring peptide.
        """

        def by_key(group: tuple[str, list[Ident]]) -> list[Ident]:
            return sorted(group[1], key=lambda id: id.sort_key())

        def rt_diff(best: Ident, others: list[Ident]) -> list[Ident]:
            for pep in others:
                pep.rt_diff = abs(best.rt - pep.rt)

            return others

        return list(
            ids  # [id]
            | groupby(lambda id: id.file_ref)  # [[file, [id]]]
            | select(by_key)  # [[id]]
            | select(lambda ids: rt_diff(ids[0], ids[1:]))  # [[id]] (remove the best)
            | traverse  # [id]
        )

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
            | select(lambda file: file.duplicates(self.max_rt_minutes))
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

            assert spectra[index].getMSLevel() == 2
            spectra[index].set_peaks([150.0], [10000.0])

        exp.setSpectra(spectra)
        oms.MzMLFile().store(str(output_path), exp)

    def dump_censored(self, ids: list[dict[str, str]]):
        if len(ids) == 0:
            raise RuntimeError("no peptides were censored!")

        with open(self.output_dir / self.manifest_file, "w") as file:
            writer = csv.DictWriter(file, fieldnames=ids[0].keys(), delimiter="\t")
            writer.writeheader()
            writer.writerows(ids)


################################################################################
class ProteomicsLFQ(Censor):
    class PepId(Ident):
        def __init__(self, pep: oms.PeptideIdentification, file_ref=None):
            """
            Create an Ident object from an OpenMS `PeptideIdentification`.
            """
            Ident.__init__(self)

            assert len(pep.getHits()) > 0
            assert pep.getScoreType() == "q-value"

            hit = pep.getHits()[0]
            aaseq = hit.getSequence()

            self.file_ref = file_ref
            self.rt = pep.getRT() / 60.0  # Seconds -> Minutes
            self.sequence = aaseq.toUnmodifiedString()
            self.sequence_with_mods = str(aaseq)
            self.score = hit.getScore()
            self.spectrum_ref = pep.getSpectrumReference()
            self.scan = None  # Must use `spectrum_ref`
            self.is_target = not hit.isDecoy()

            self.is_human = any(
                hit.extractProteinAccessionsSet()
                | select(lambda s: re.search(r"human", s, re.IGNORECASE))
            )

            self.precursor_mass = (
                pep.getMZ() * hit.getCharge() - hit.getCharge() * 1.007276466771
            )

        def to_row(self):
            return {
                "File": self.file_ref,
                "Full Sequence": self.sequence_with_mods,
                "Base Sequence": self.sequence,
                "Score": self.score,
                "Is Target": self.is_target,
                "Is Human": self.is_human,
            }

    class IdXmlFile(IdFile):
        def __init__(self, file: Path):
            IdFile.__init__(self, file)
            self.mzml_file = file.parent / (file.stem + ".mzML")

        def mzml_file_path(self) -> Path:
            return self.mzml_file

        def load(self, max_score: float):
            _, peptide_ids = oms.IdXMLFile().load(str(self.id_file))
            to_ident = lambda pep: ProteomicsLFQ.PepId(pep, self.name())

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
        return self.censor(self.id_files)

    def _load_data(self):
        self.id_files = dict(
            self.input_dir.glob("*.idXML")
            | select(lambda path: ProteomicsLFQ.IdXmlFile(path))
            | select(lambda file: [file.name(), file.load(self.max_psm_qval)])
        )

        for name, file in self.id_files.items():
            self.log.info(f"assigned {len(file.peptides)} to file {name}")


################################################################################
class FlashLFQ(Censor):
    class PepId(Ident):
        def __init__(self, row: dict[str, str]):
            Ident.__init__(self)
            row["File Name"] = row["File Name"].removesuffix("-calib")

            self.row = row
            self.file_ref = row["File Name"]
            self.rt = float(row["Scan Retention Time"])
            self.sequence = row["Base Sequence"]
            self.sequence_with_mods = row["Full Sequence"]
            self.score = float(row["PEP_QValue"])
            # self.spectrum_ref = None  # Must use `scan`
            self.scan = int(row["Scan Number"])
            self.is_target = row["Decoy/Contaminant/Target"] == "T"
            self.is_human = "Homo sapiens" in row["Organism Name"]
            self.precursor_mass = float(row["Precursor Mass"])

        def to_row(self):
            return self.row

    class MzMLFile(IdFile):
        def __init__(self, file: Path):
            file = file.parent / file.stem.removesuffix("-calib")
            IdFile.__init__(self, file)

        def mzml_file_path(self) -> Path:
            return self.id_file.parent / (self.id_file.stem + ".mzML")

    def __init__(self, psms_file, mzml_files):
        Censor.__init__(self)

        self.psms_file = psms_file
        self.mzml_files = mzml_files

    def run(self):
        self.log.info("censoring spectra for FlashLFQ")
        self._load_data()
        return self.censor(self.id_files)

    def _load_data(self):
        self.id_files = dict(
            self.mzml_files
            | select(FlashLFQ.MzMLFile)
            | select(lambda f: [f.id_file.stem, f])
        )

        valid_sequences = set()
        with open(self.psms_file.parent / "AllPeptides.psmtsv", "r") as file:
            for row in csv.DictReader(file, delimiter="\t"):
                pep = FlashLFQ.PepId(row)
                # NOTE: The paper says we should use
                # `max_peptide_qval` here but that leaves no peptides
                # to censor!
                if pep.is_valid(self.max_psm_qval):
                    valid_sequences.add(pep.sequence_with_mods)

        psms = []
        invalid_sequences = 0
        with open(self.psms_file, "r") as file:
            for row in csv.DictReader(file, delimiter="\t"):
                pep = FlashLFQ.PepId(row)

                if pep.is_valid(self.max_psm_qval):
                    if pep.sequence_with_mods in valid_sequences:
                        psms.append(pep)
                    else:
                        invalid_sequences += 1

        for group in psms | groupby(lambda p: p.file_ref):
            file = group[0]
            peps = list(group[1])

            self.log.info(f"assigned {len(peps)} peptides to file {group[0]}")
            assert file in self.id_files
            self.id_files[file].peptides.extend(peps)

        self.log.info(f"filtered out {invalid_sequences} invalid PSMS")


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
    parser.add_argument("--log", type=Path)
    parser.add_argument("--tool", required=True, choices=["flashlfq", "proteomicslfq"])
    parser.add_argument("files", action="extend", nargs="+", type=Path)

    args = parser.parse_args()

    if args.log is not None:
        logopts["filename"] = args.log

    logging.basicConfig(**logopts)

    censor = None

    if args.tool == "flashlfq":
        censor = FlashLFQ(args.files[0], args.files[1:])
    elif args.tool == "proteomicslfq":
        censor = ProteomicsLFQ(*args.files)

    censor.run()
