#!/usr/bin/env python


import argparse
import logging
import re
from pathlib import Path

import pandas as pd


################################################################################
class DataHelper:
    def __init__(self, data_dir: Path):
        self.data_dir = data_dir
        self.entrapment_seqs_set = set()

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
class FlashLFQ:
    def __init__(self, data_dir, orig_path, censored_path):
        self.data = DataHelper(data_dir)
        self.orig_path = orig_path
        self.censored_path = censored_path

        # These need to be configured from the outside depending on
        # the data set:
        self.human_file_pattern = "_1x02nguL_"
        self.foreign_species = "Yeast"
        self.q_value = 0.01
        self.rt_delta = 40 * 0.01

        # Scaling factor for entrapment calculations, based on the
        # ratio on the number of entrapment peptide sequences and
        # human peptide sequences in the database
        self.scaling_factor = 1 + 3105275 / 3608159

        self.log = logging.getLogger(__name__)

    def run(self):
        stats = self._fdp()
        print(stats["false_discovery_proportion"])

    def _fdp(self):
        peaks = self._load_peaks(self.orig_path)
        entrapment_seqs = self.data.entrapment_seqs()

        # Limit to the human mzML files:
        peaks = peaks.loc[peaks["File Name"].str.contains(self.human_file_pattern)]

        msms_detected_yeast_peaks = peaks.loc[
            (peaks["Peak Detection Type"] == "MSMS")
            & (peaks["Organism"].str.contains(self.foreign_species))
        ]

        msms_detected_yeast_peaks.drop_duplicates(
            subset=["Full Sequence"], keep="first", inplace=True
        )

        msms_yeast_seqs = msms_detected_yeast_peaks["Full Sequence"].tolist()

        msms_peaks = peaks.loc[
            (peaks["Peak Detection Type"] == "MSMS")
            & (peaks["Decoy Peptide"] == False)
            & (
                (peaks["Organism"].str.contains(self.foreign_species))
                | (peaks["Organism"].str.contains("Homo sapiens"))
                | (peaks["Organism"].str.contains("Arabidopsis"))
            )
        ]

        msms_count = msms_peaks.shape[0]

        peaks = peaks.loc[
            (peaks["Peak Detection Type"] == "MBR")
            & (peaks["PIP Q-Value"] < self.q_value)
            & (peaks["Decoy Peptide"] == False)
        ]

        peaks = peaks.loc[
            (peaks["Random RT"] == False) & (peaks["Decoy Peptide"] == False)
        ]

        human_peaks = peaks.loc[peaks["Organism"].str.contains("Homo sapiens")]

        # Our entrapment database has a bunch of scrambled sequences
        # appended to the end of the DB, which is why we're checking
        # for them
        entrapment_peaks = human_peaks.loc[
            human_peaks["Full Sequence"].isin(entrapment_seqs)
        ]

        human_peaks = human_peaks.loc[
            ~human_peaks["Base Sequence"].isin(entrapment_seqs)
        ]

        human_count = human_peaks.shape[0]
        arabida_count = entrapment_peaks.shape[0]

        yeast_peaks = peaks.loc[peaks["Organism"].str.contains(self.foreign_species)]
        yeast_peaks = yeast_peaks.loc[
            ~yeast_peaks["Organism"].str.contains("Homo sapiens")
        ]
        yeast_count = yeast_peaks.shape[0]

        msms_yeast_count = sum(yeast_peaks["Full Sequence"].isin(msms_yeast_seqs))
        yeast_count = yeast_count - msms_yeast_count

        score_threshold = peaks["PIP PEP"].max()

        # calculate the error rate for (presumably) native human transfers
        nper_results = self._get_native_peak_error_rate(score_threshold)

        scaled_arabida = arabida_count * self.scaling_factor
        scaled_nper = human_count * nper_results["error_rate"]

        total = human_count + yeast_count + arabida_count
        total = max(1, total)

        fper = 100 * yeast_count / total
        nper = 100 * scaled_nper / total
        fder = 100 * scaled_arabida / total

        false_discovery_proportion = fper + nper + fder

        return {
            "human": human_count,
            "foreign": yeast_count,
            "arabida": arabida_count,
            "total": total,
            "msms_count": msms_count,
            "scaled_arabida": scaled_arabida,
            "scaled_nper": scaled_nper,
            "sensitivity": nper_results["sensitivity"],
            "fper": fper,
            "nper": nper,
            "fder": fder,
            "false_discovery_proportion": false_discovery_proportion,
        }

    def _get_native_peak_error_rate(self, score_threshold):
        original_peaks = self._load_peaks(self.orig_path)
        original_peaks = original_peaks[original_peaks["Peak RT Apex"] != "-"]

        censored_psms = self._load_censored_psms()
        censored_peaks = pd.merge(
            censored_psms,
            original_peaks,
            on=["File Name", "Full Sequence"],
        )

        censored_peaks = censored_peaks[
            [
                "File Name",
                "Full Sequence",
                "Peak RT Start",
                "Peak RT Apex",
                "Peak RT End",
                "Peak Charge",
                "Peak intensity",
            ]
        ]

        new_peaks = self._load_peaks(self.censored_path)
        new_peaks = new_peaks.loc[
            (new_peaks["Peak Detection Type"] == "MBR")
            & (new_peaks["Random RT"] == False)
            & (new_peaks["Decoy Peptide"] == False)
        ]

        new_peaks = new_peaks.loc[(new_peaks["PIP Q-Value"] < self.q_value)]
        new_peaks = new_peaks[
            [
                "File Name",
                "Full Sequence",
                "Peak RT Start",
                "Peak RT Apex",
                "Peak RT End",
                "Peak Charge",
                "Peak intensity",
                "PIP PEP",
                "PIP Q-Value",
            ]
        ]

        peak_join = pd.merge(
            censored_peaks,
            new_peaks,
            on=["File Name", "Full Sequence"],
        )

        peak_join["Agreement"] = peak_join.apply(
            lambda x: self._check_overlap_time(x, self.rt_delta), axis=1
        )

        peak_join.sort_values(by="Agreement", ascending=False, inplace=True)

        # compare the old and new results (here, 1 means the match, 0 means they dont)
        peak_join.drop_duplicates(
            subset=["File Name", "Full Sequence"], keep="first", inplace=True
        )

        censored_peaks.drop_duplicates(
            subset=["File Name", "Full Sequence"], keep="first", inplace=True
        )
        sensitivity = peak_join.shape[0] / censored_peaks.shape[0]

        if peak_join.shape[0] == 0:
            return {"sensitivity": 0, "error_rate": 0}

        peak_join.sort_values(by="PIP Q-Value", ascending=True, inplace=True)
        peak_join["good_transfers"] = (peak_join["Agreement"] == 1).cumsum()
        peak_join["bad_transfers"] = (peak_join["Agreement"] == 0).cumsum()
        peak_join["error_rate"] = peak_join.apply(
            lambda x: self._calculate_error_rate(x), axis=1
        )

        return {
            "sensitivity": sensitivity,
            "error_rate": peak_join["error_rate"].iloc[-1],
        }

    def _check_overlap_time(self, table, rt_delta):
        pip_apex = float(table["Peak RT Apex_y"])
        ms_apex = float(table["Peak RT Apex_x"])

        if abs(pip_apex - ms_apex) < rt_delta:
            return 1
        else:
            return 0

    def _calculate_error_rate(self, table):
        return table["bad_transfers"] / (
            table["good_transfers"] + table["bad_transfers"]
        )

    def _load_peaks(self, directory):
        file_name = directory / "QuantifiedPeaks.tsv"
        peaks = pd.read_csv(file_name, sep="\t")

        pattern = re.compile(r"(-censored)?-(calib)$")
        peaks["File Name"] = peaks["File Name"].str.replace(pattern, "", regex=True)

        return peaks

    def _load_censored_psms(self):
        file_name = self.censored_path / "censored.psmtsv"
        return pd.read_csv(file_name, sep="\t")


################################################################################
if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)

    parser = argparse.ArgumentParser(
        prog="analysis.py",
        description="Analyze the performance for PIP FDR",
    )

    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--tool", required=True, choices=["flashlfq"])
    parser.add_argument("file", action="extend", nargs="+", type=Path)

    args = parser.parse_args()

    if args.tool == "flashlfq":
        flashlfq = FlashLFQ(args.data, *args.file)
        flashlfq.run()
