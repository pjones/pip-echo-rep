#!/usr/bin/env python


import argparse
import json
import logging
import re
from pathlib import Path

import pandas as pd
import pyopenms as oms
from pipe import select, where


################################################################################
class DataHelper:
    def __init__(
        self,
        data_dir: Path,
        dataset_name: str,
        censored_manifest: Path,
        orig_path: Path,
        censored_path: Path,
    ):
        self.data_dir = data_dir
        self.dataset_name = dataset_name
        self.entrapment_seqs_set = set()

        self.censored_manifest = censored_manifest
        self.orig_path = orig_path
        self.censored_path = censored_path

        # These need to be configured from the outside depending on
        # the data set:
        self.tool = None
        self.human_file_pattern = None
        self.foreign_species = None
        self.q_value = None

        # Max RT (in minutes) that two peptide identifications are
        # allowed to be.  Exceeding this threshold means the features
        # were linked by mistake.
        #
        # Each dataset has its own number:
        #
        # - 54s (Lim et al.)
        # - 36s (E. coli)
        # - 24s (Single-cell)
        self.rt_delta = None

    def path(self, file):
        return self.data_dir / "PXD057758" / file

    def entrapment_seqs(self):
        if len(self.entrapment_seqs_set) > 0:
            return self.entrapment_seqs_set
        else:
            file_name = self.path("EntrapmentProteinPeptideSequences_50percent.txt")
            with open(file_name) as file:
                for line in file:
                    self.entrapment_seqs_set.add(line.strip())

            return self.entrapment_seqs_set


################################################################################
# Adapted from: https://github.com/Alexander-Sol/PIP-ECHOanalysis.git
class Analysis:
    def __init__(self, data: DataHelper):
        self.data = data
        self.log = logging.getLogger(__name__)
        self.file_name_col = "File Name"
        self.sequence_col = "Full Sequence"

    def run(self, original, censored, manifest):
        mbr_original = len(original[original["Peak Detection Type"] == "MBR"].index)
        mbr_censored = len(censored[censored["Peak Detection Type"] == "MBR"].index)

        json = self.fdp(original, censored, manifest)

        json["mbr_original"] = mbr_original
        json["mbr_censored"] = mbr_censored

        return json

    def fdp(self, original, censored, manifest):
        key = [self.file_name_col, self.sequence_col]
        entrapment_seqs = self.data.entrapment_seqs()

        # Limit to human files.
        original = original.loc[
            original[key[0]].str.contains(self.data.human_file_pattern)
        ]

        msms_detected_foreign_peaks = original.loc[
            (original["Peak Detection Type"] == "MSMS")
            & (original["Organism"].str.contains(self.data.foreign_species, na=False))
        ]
        msms_detected_foreign_peaks = msms_detected_foreign_peaks.drop_duplicates(
            subset=[key[1]], keep="first"
        )
        msms_foreign_seqs = msms_detected_foreign_peaks[key[1]].tolist()

        msms_peaks = original.loc[
            (original["Peak Detection Type"] == "MSMS")
            & (original["Decoy Peptide"] == False)
            & (
                (original["Organism"].str.contains(self.data.foreign_species, na=False))
                | (original["Organism"].str.contains("Homo sapiens", na=False))
                | (original["Organism"].str.contains("HUMAN", na=False))
                | (original["Organism"].str.contains("Arabidopsis", na=False))
            )
        ]

        msms_count = msms_peaks.shape[0]

        mbr_peaks = original.loc[
            (original["Peak Detection Type"] == "MBR")
            & (original["PIP Q-Value"] <= self.data.q_value)
            & (original["Decoy Peptide"] == False)
            & (original["Random RT"] == False)
        ]

        human_peaks = mbr_peaks.loc[
            (mbr_peaks["Organism"].str.contains("Homo sapiens", na=False))
            | (mbr_peaks["Organism"].str.contains("HUMAN", na=False))
        ]

        entrapment_peaks = human_peaks.loc[
            human_peaks["Full Sequence"].isin(entrapment_seqs)
        ]

        human_peaks = human_peaks.loc[
            ~human_peaks["Base Sequence"].isin(entrapment_seqs)
        ]

        human_count = human_peaks.shape[0]
        arabida_count = entrapment_peaks.shape[0]

        foreign_peaks = mbr_peaks.loc[
            mbr_peaks["Organism"].str.contains(self.data.foreign_species, na=False)
        ]
        foreign_peaks = foreign_peaks.loc[
            ~(
                (foreign_peaks["Organism"].str.contains("Homo sapiens", na=True))
                | (foreign_peaks["Organism"].str.contains("HUMAN", na=True))
            )
        ]
        foreign_count = foreign_peaks.shape[0]

        msms_foreign_count = sum(foreign_peaks[key[1]].isin(msms_foreign_seqs))
        foreign_count = foreign_count - msms_foreign_count
        nper_results = self.native_error(original, censored, manifest)

        # Scaling factor for entrapment calculations, based on the
        # ratio on the number of entrapment peptide sequences and
        # human peptide sequences in the database
        #
        # scaling_factor = (arabida_count + human_count) / (1.0 + arabida_count)
        #
        # Using the value from the paper, which uses the count from
        # the FASTA file, not from the discovered peptides:
        scaling_factor = 1.0 + 3105275.0 / 3608159.0
        scaled_arabida = arabida_count * scaling_factor
        scaled_nper = human_count * nper_results["error_rate"]

        total = human_count + foreign_count + arabida_count
        total = max(1, total)

        fper = 100 * foreign_count / total
        nper = 100 * scaled_nper / total
        fder = 100 * scaled_arabida / total

        false_discovery_proportion = fper + nper + fder

        return {
            "tool": self.data.tool,
            "dataset": self.data.dataset_name,
            "q_value": self.data.q_value,
            "max_rt_delta_min": self.data.rt_delta,
            "max_rt_delta_sec": self.data.rt_delta * 60.0,
            "human": human_count,
            "foreign": foreign_count,
            "arabida": arabida_count,
            "total": total,
            "msms_count": msms_count,
            "scaling_factor": scaling_factor,
            "scaled_arabida": scaled_arabida,
            "scaled_nper": scaled_nper,
            "sensitivity": nper_results["sensitivity"],
            "error_rate": nper_results["error_rate"],
            "fper": fper,
            "nper": nper,
            "fder": fder,
            "false_discovery_proportion": false_discovery_proportion,
        }

    def native_error(self, original, censored, manifest):
        key = [self.file_name_col, self.sequence_col]

        # Only MSMS detected peptides should have been removed by the
        # censor script, but just in case:
        msms = original[original["Peak Detection Type"] == "MSMS"]

        # These are the MSMS peptides removed from the mzML files:
        removed = pd.merge(msms, manifest, on=key)
        removed = removed[[self.file_name_col, self.sequence_col, "RT"]]
        removed.to_csv("removed_df.csv")

        censored_mbr = censored.loc[
            (censored["Peak Detection Type"] == "MBR")
            & (censored["Decoy Peptide"] == False)
            & (censored["PIP Q-Value"] <= self.data.q_value)
        ]

        # These are the peptides that were removed from the original
        # run but then found via MBR:
        mbr = pd.merge(censored_mbr, removed, on=key)
        mbr["RT_d"] = mbr.apply(lambda row: abs(row["RT_x"] - row["RT_y"]), axis=1)
        mbr["Agreement"] = mbr.apply(lambda r: r["RT_d"] <= self.data.rt_delta, axis=1)
        mbr.sort_values(by="Agreement", ascending=False, inplace=True)
        mbr.to_csv("mbr_df_pre.csv")

        # FIXME: why isn't there more overlap between manifest and
        # removed?  This might indicate a bug in the censor script.
        self.log.info(f"    MSMS peptides: {msms.shape}")
        self.log.info(f"manifest peptides: {manifest.shape}")
        self.log.info(f" removed peptides: {removed.shape}")
        self.log.info(f"     MBR peptides: {mbr.shape}")

        mbr = mbr.drop_duplicates(subset=key, keep="first")
        manifest = manifest.drop_duplicates(subset=key, keep="first")

        mbr.sort_values(by="PIP Q-Value", ascending=True, inplace=True)
        mbr["good_transfers"] = (mbr["Agreement"] == True).cumsum()
        mbr["bad_transfers"] = (mbr["Agreement"] == False).cumsum()

        mbr["error_rate"] = mbr.apply(
            lambda row: (
                row["bad_transfers"] / (row["good_transfers"] + row["bad_transfers"])
            ),
            axis=1,
        )
        mbr.to_csv("mbr_df_post.csv")

        if mbr.shape[0] == 0 or manifest.shape[0] == 0:
            return {"sensitivity": 0.0, "error_rate": 0.0}
        else:
            return {
                "sensitivity": mbr.shape[0] / manifest.shape[0],
                "error_rate": mbr["error_rate"].iloc[-1],
            }


################################################################################
class ProteomicsLFQ(Analysis):
    class Peptide:
        def __init__(
            self,
            file_ref: str,
            detection_type: str,
            mbr_qval: float,
            pep: oms.PeptideIdentification,
            feat: oms.FeatureHandle | None = None,
        ):
            self.file_ref = file_ref
            self.detection_type = detection_type
            self.mbr_qvalue = mbr_qval

            hit = pep.getHits()[0]
            aaseq = hit.getSequence()

            self.sequence = aaseq.toUnmodifiedString()
            self.sequence_with_mods = str(aaseq)
            self.score = hit.getScore()
            self.is_decoy = hit.isDecoy()
            self.organism = str(hit.extractProteinAccessionsSet())

            # Note: Retention time in OpenMS is measured in seconds.
            # However, other tools use minutes so the remainder of
            # this script does too.
            if feat is not None:
                self.rt = feat.getRT() / 60.0
            else:
                self.rt = pep.getRT() / 60.0

        def to_row(self):
            return {
                "File": self.file_ref,
                "Full Sequence": self.sequence_with_mods,
                "Base Sequence": self.sequence,
                "Peak Detection Type": self.detection_type,
                "Decoy Peptide": self.is_decoy,
                "Organism": self.organism,
                "Score": self.score,
                "PIP Q-Value": self.mbr_qvalue,
                "PIP PEP": 1 - self.mbr_qvalue,
                "RT": self.rt,
                "Random RT": False,
            }

    class IdFile:
        def __init__(self, file):
            self.maps = {}
            self.peptides = []
            self.empty_features = 0
            self.mbr_events = 0

            self.log = logging.getLogger(__name__)
            self._load(file)

        def to_df(self):
            return pd.DataFrame(self.peptides | select(lambda p: p.to_row()))

        def _load(self, file):
            self.log.info(f"loading identifications from {file}")

            consensus_map = oms.ConsensusMap()
            oms.ConsensusXMLFile().load(str(file), consensus_map)

            for id, header in consensus_map.getColumnHeaders().items():
                base = str(Path(header.filename).stem)
                self.maps[id] = base.replace("-censored", "")

            for consensus_feature in consensus_map:
                peps = list(
                    consensus_feature.getPeptideIdentifications()
                    | where(lambda pep: len(pep.getHits()) > 0)
                )

                if len(peps) == 0:
                    self.empty_features += 1
                    continue

                assert len(peps) == 1

                self.peptides.extend(
                    self._peptides_from_feature(consensus_feature, peps[0])
                )

            self.log.info(
                f"loaded {len(self.peptides)} peptides with {self.mbr_events} MBR events"
            )

        def _peptides_from_feature(self, feature, pep):
            results = []

            seen_maps = set()
            mbr_maps = []
            mbr_qvals = []

            if "mbr_transfer_map_indices" in feature.getKeys():
                mbr_maps = feature.getMetaValue("mbr_transfer_map_indices")
                mbr_qvals = feature.getMetaValue("mbr_transfer_qvalues")
                self.mbr_events += len(mbr_maps)

                if mbr_qvals is None or mbr_qvals == []:
                    mbr_qvals = [1.0] * len(mbr_maps)

            for handle in feature.getFeatureList():
                detection = "MSMS"
                qval = 0.0

                idx = handle.getMapIndex()
                seen_maps.add(idx)

                assert idx in self.maps
                file_ref = self.maps[idx]

                if idx in mbr_maps:
                    place = mbr_maps.index(idx)
                    qval = mbr_qvals[place]
                    detection = "MBR"

                results.append(
                    ProteomicsLFQ.Peptide(file_ref, detection, qval, pep, handle)
                )

            # Ensure that we saw all of the maps recorded as accepting
            # an MBR transfer.
            for idx in mbr_maps:
                assert idx in seen_maps

            return results

    def __init__(self, data: DataHelper):
        Analysis.__init__(self, data)
        self.file_name_col = "File"

        # Match the names exported by OpenMS:
        if self.data.foreign_species == "Saccharomyces cerevisiae":
            self.data.foreign_species = "YEAST"
        else:
            raise RuntimeError(f"Does OpenMS know what {self.data.foreign_species} is?")

    def run(self):
        self.log.info("analyzing ProteomicsLFQ data")
        self._load()
        return Analysis.run(self, self.original_df, self.censored_df, self.manifest_df)

    def _load(self):
        file_name = "proteomicslfq.consensusXML"
        self.original = ProteomicsLFQ.IdFile(self.data.orig_path / file_name)
        self.censored = ProteomicsLFQ.IdFile(self.data.censored_path / file_name)

        self.original_df = self.original.to_df()
        self.censored_df = self.censored.to_df()
        self.manifest_df = pd.read_csv(self.data.censored_manifest, sep="\t")

        self.original_df.to_csv("original_df.csv")
        self.censored_df.to_csv("censored_df.csv")
        self.manifest_df.to_csv("manifest_df.csv")


################################################################################
class FlashLFQ(Analysis):
    def __init__(self, data: DataHelper):
        Analysis.__init__(self, data)

    def run(self):
        self.log.info("analyzing output from FlashLFQ")

        orig = self._load_peaks(self.data.orig_path)
        censored = self._load_peaks(self.data.censored_path)
        manifest = self._load_censored_psms()

        return Analysis.run(self, orig, censored, manifest)

    def _load_peaks(self, directory):
        file_name = directory / "QuantifiedPeaks.tsv"
        peaks = pd.read_csv(file_name, sep="\t")

        pattern = re.compile(r"(-censored)?-(calib)$")
        peaks["File Name"] = peaks["File Name"].str.replace(pattern, "", regex=True)

        return self._fixup_df(peaks)

    def _load_censored_psms(self):
        manifest = pd.read_csv(self.data.censored_manifest, sep="\t")
        return self._fixup_df(manifest)

    def _fixup_df(self, df):
        """
        Make the data look like what the `Analysis` class expects.
        """
        if "Peak RT Apex" in df.columns:
            df = df[df["Peak RT Apex"] != "-"]
            df["RT"] = df.apply(lambda r: float(r["Peak RT Apex"]), axis=1)

        return df


################################################################################
if __name__ == "__main__":
    logopts = {
        "format": "[%(asctime)s %(levelname)s]: %(message)s",
        "datefmt": "%Y-%m-%dT%H:%M:%S",
        "level": logging.DEBUG,
    }

    parser = argparse.ArgumentParser(
        prog="analysis.py",
        description="Analyze the performance for PIP FDR",
    )

    parser.add_argument("--id", required=True, type=str)
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--censored-manifest", required=True, type=Path)
    parser.add_argument("--output", default="analysis.json", type=Path)
    parser.add_argument("--log", type=Path)

    parser.add_argument("--tool", required=True, choices=["flashlfq", "proteomicslfq"])
    parser.add_argument("--qval", default=0.01, type=float)
    parser.add_argument("--human-file-pattern", default="_1x02nguL_")
    parser.add_argument("--foreign-species", default="Saccharomyces cerevisiae")
    parser.add_argument("--rt-delta", default=(24.0 / 60.0), type=float)

    parser.add_argument("file", action="extend", nargs="+", type=Path)

    args = parser.parse_args()

    if args.log is not None:
        logopts["filename"] = args.log

    logging.basicConfig(**logopts)

    data = DataHelper(args.data, args.id, args.censored_manifest, *args.file)
    data.tool = args.tool
    data.q_value = args.qval
    data.human_file_pattern = args.human_file_pattern
    data.foreign_species = args.foreign_species
    data.rt_delta = args.rt_delta

    results = None

    if args.tool == "flashlfq":
        flashlfq = FlashLFQ(data)
        results = flashlfq.run()
    elif args.tool == "proteomicslfq":
        proteomicslfq = ProteomicsLFQ(data)
        results = proteomicslfq.run()

    if results is not None:
        with open(args.output, "w") as file:
            file.write(json.dumps(results, indent=2))
