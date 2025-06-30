# AIMS ancestry inference pipeline integrated into wasptk
from __future__ import annotations

import os
import re
import subprocess
import tempfile
from itertools import islice
from pathlib import Path
from typing import Dict, List, Tuple

import pysam

_DATA_DIR = Path(__file__).resolve().parent / "ancestry"

# resource files
_AIM_POSITIONS = _DATA_DIR / "Capture96AIMPos.txt"
_REF_GT = _DATA_DIR / "1KGP_sup_AIM_GT.txt"
_POP_CODE = _DATA_DIR / "1KGP_sup_Pop_code.txt"
_EXTRA_PARAMS = _DATA_DIR / "extraparams"
_MAIN_PARAMS = _DATA_DIR / "mainparams"


def _read_aim_pos(path: Path) -> Tuple[List[str], List[str]]:
    chr_no: List[str] = []
    aim_pos: List[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            chr_no.append(line.split('_')[0])
            aim_pos.append(line.split('_')[1])
    all_aim_pos = [f"{c}_{p}" for c, p in zip(chr_no, aim_pos)]
    return all_aim_pos, aim_pos


def _read_vcf(vcf_file: str) -> Tuple[str, Dict[str, str]]:
    """Read a simple VCF (optionally gzipped) and return sample name and genotypes."""
    sample = ""
    gt_map: Dict[str, str] = {}
    opener = open
    if vcf_file.endswith(".gz"):
        import gzip

        opener = gzip.open
    with opener(vcf_file, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                fields = line.strip().split("\t")
                sample = fields[9]
                continue
            fields = line.strip().split("\t")
            chrom = fields[0]
            if not chrom.startswith("chr"):
                chrom = f"chr{chrom}"
            pos = fields[1]
            key = f"{chrom}_{pos}"
            gt = fields[9].split(":")[0]
            gt_map[key] = gt
    return sample, gt_map


def _add_sample_code(ref_sample_no: int, sample_name: str) -> Dict[str, int]:
    return {sample_name: ref_sample_no + 1}


def _convert_to_structure(
    aim_pos: List[str],
    gt_map: Dict[str, str],
    sample_code: Dict[str, int],
) -> Tuple[List[str], List[str], List[str]]:
    all_gt: List[str] = []
    high_cov: List[str] = []
    for pos in aim_pos:
        gt = gt_map.get(pos, "-9/-9")
        if gt not in {"-9/-9", "./."}:
            high_cov.append(pos)
        all_gt.append(gt)
    hap1: List[str] = []
    hap2: List[str] = []
    for gt in all_gt:
        a1, a2 = gt.split("/")
        a1 = a1.replace("1", "2").replace("0", "1") if a1 != "-9" else "-9"
        a2 = a2.replace("1", "2").replace("0", "1") if a2 != "-9" else "-9"
        hap1.append(a1)
        hap2.append(a2)
    code = list(sample_code.values())[0]
    hap1.insert(0, str(code))
    hap2.insert(0, str(code))
    hap1 = hap1[:1] + ["1", "0"] + hap1[1:]
    hap2 = hap2[:1] + ["1", "0"] + hap2[1:]
    return hap1, hap2, high_cov


def _write_input_file(hap1: List[str], hap2: List[str], out_path: Path) -> None:
    with open(out_path, "w") as out, open(_REF_GT) as ref:
        for line in ref:
            out.write(line)
        out.write(" ".join(hap1) + "\n")
        out.write(" ".join(hap2) + "\n")


def _run_structure(work_dir: Path) -> None:
    subprocess.run(
        [
            "structure",
            "-i",
            "1KGP_SuperPop_with_SampleGT.txt",
            "-m",
            str(_MAIN_PARAMS),
            "-e",
            str(_EXTRA_PARAMS),
        ],
        cwd=work_dir,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def _parse_structure_result(
    result_file: Path, pop_file: Path, sample_no: int, high_cov: List[str]
) -> Dict[str, object]:
    InputPCode: List[str] = []
    InputPop: List[str] = []
    with open(pop_file) as code:
        for line in code:
            line = line.strip("\n")
            InputPCode.append(line.split("\t")[1])
            InputPop.append(line.split("\t")[0])
    Clno = 5
    number = [str(i) for i in range(1, Clno + 1)]
    mymatrix: List[str] = []
    Popvalue: List[str] = []
    ResultLine: str | None = None
    with open(result_file) as old:
        for line in old:
            line = line.rstrip()
            if "Given    Inferred Clusters" in line:
                mymatrix.append(" ".join(islice(old, len(InputPop) + 2)))
            if line.startswith(str(sample_no)):
                ResultLine = line
    if ResultLine is None:
        raise RuntimeError("Failed to parse structure output")
    ResultLine = re.sub(r"\s+", " ", ResultLine)
    ResultLine = ResultLine.split(":")[1]
    result_values = [float(n) for n in ResultLine.split(" ")[1 : Clno + 1]]
    for n in mymatrix:
        char = n.split("\n")
    del char[0:2]
    del char[-1]
    for line in char:
        line = line.strip(" ")
        Popvalue.append(line.split(":")[1])
    AllPopValue: List[List[str]] = []
    for n in Popvalue:
        n = re.sub(r"\s+", " ", n)
        n = n.split(" ")
        AllPopValue.append(n[1 : Clno + 1])
    FinalPopValue = list(map(list, zip(*AllPopValue)))
    PopDict = {str(a): [InputPop, val] for a, val in zip(number, FinalPopValue)}
    for key, value in PopDict.items():
        PopDict[key] = dict(zip(*value))
    Mysampledict = {i + 1: v for i, v in enumerate(result_values)}
    Mysampledict = dict(
        sorted(Mysampledict.items(), key=lambda item: item[1], reverse=True)
    )
    ancestry_labels = {
        "1": "European",
        "2": "American",
        "3": "African",
        "4": "East Asian",
        "5": "South Asian",
    }
    max_ancestry = max(Mysampledict, key=Mysampledict.get)
    inferred_label = ancestry_labels.get(str(max_ancestry), str(max_ancestry))
    probabilities = {ancestry_labels.get(str(k), str(k)): v for k, v in Mysampledict.items()}
    return {
        "inferred_ancestry": inferred_label,
        "probabilities": probabilities,
        "aims_used": len(high_cov),
        "aims_missing": 96 - len(high_cov),
    }


def run_aims(vcf: str) -> Dict[str, object]:
    all_aim_pos, _ = _read_aim_pos(_AIM_POSITIONS)
    sample_name, gt_map = _read_vcf(vcf)
    sample_code = _add_sample_code(2504, sample_name)
    hap1, hap2, high_cov = _convert_to_structure(all_aim_pos, gt_map, sample_code)
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        input_file = tmpdir / "1KGP_SuperPop_with_SampleGT.txt"
        _write_input_file(hap1, hap2, input_file)
        try:
            _run_structure(tmpdir)
        except FileNotFoundError as exc:
            raise RuntimeError("structure command not found") from exc
        result_file = tmpdir / "OutfileK5withSuperPop_new_f"
        return _parse_structure_result(result_file, _POP_CODE, 2505, high_cov)
