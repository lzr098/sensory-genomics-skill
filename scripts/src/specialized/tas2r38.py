"""TAS2R38 Haplotype 分析器.

基于 rs713598、rs1726866、rs10246939 三个 SNP 判定 diplotype，
输出苦味感知能力评估。
"""

import json
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

from src.logger import get_logger
from src.models import TAS2R38Result, Variant

logger = get_logger(__name__)


# 内置 haplotype 定义（作为 fallback）
_DEFAULT_HAPLOTYPE_DEFS = {
    "PAV": {"rs713598": "C", "rs1726866": "T", "rs10246939": "G"},
    "AVI": {"rs713598": "G", "rs1726866": "C", "rs10246939": "A"},
    "AAI": {"rs713598": "G", "rs1726866": "C", "rs10246939": "G"},
    "PVI": {"rs713598": "C", "rs1726866": "T", "rs10246939": "A"},
}

_DEFAULT_DIPLOTYPE_PHENOTYPES = {
    "PAV/PAV": {"phenotype": "对苦味化合物高度敏感", "level": "苦味高敏感型"},
    "PAV/AVI": {"phenotype": "对苦味化合物中等敏感", "level": "苦味中等敏感型"},
    "PAV/AAI": {"phenotype": "对苦味化合物中等敏感", "level": "苦味中等敏感型"},
    "PAV/PVI": {"phenotype": "对苦味化合物高度敏感", "level": "苦味高敏感型"},
    "AVI/AVI": {"phenotype": "对苦味化合物不敏感", "level": "苦味不敏感型"},
    "AVI/AAI": {"phenotype": "对苦味化合物不敏感", "level": "苦味不敏感型"},
    "AVI/PVI": {"phenotype": "对苦味化合物中等敏感", "level": "苦味中等敏感型"},
    "AAI/AAI": {"phenotype": "对苦味化合物不敏感", "level": "苦味不敏感型"},
    "AAI/PVI": {"phenotype": "对苦味化合物中等敏感", "level": "苦味中等敏感型"},
    "PVI/PVI": {"phenotype": "对苦味化合物高度敏感", "level": "苦味高敏感型"},
}


class TAS2R38Analyzer:
    """TAS2R38 Haplotype 分析器."""

    def __init__(self, data_path: Optional[str] = None) -> None:
        """初始化分析器.

        Args:
            data_path: tas2r38_snps.json 路径，默认查找 data/ 目录。
        """
        if data_path is None:
            src_dir = Path(__file__).resolve().parent.parent
            data_path = src_dir.parent / "data" / "tas2r38_snps.json"
        else:
            data_path = Path(data_path)

        self.snp_defs: List[Dict[str, Any]] = []
        self.haplotype_defs: Dict[str, Any] = {}
        self.diplotype_phenotypes: Dict[str, Any] = {}
        self._load_data(data_path)

    def _load_data(self, path: Path) -> None:
        """加载 SNP 定义数据."""
        if not path.exists():
            logger.warning("tas2r38_snps.json not found, using built-in defaults")
            self.haplotype_defs = _DEFAULT_HAPLOTYPE_DEFS.copy()
            self.diplotype_phenotypes = _DEFAULT_DIPLOTYPE_PHENOTYPES.copy()
            self.snp_defs = [
                {"rsid": "rs713598", "chrom": "7", "pos": 141972755, "ref": "G", "alt": "C"},
                {"rsid": "rs1726866", "chrom": "7", "pos": 141973545, "ref": "T", "alt": "C"},
                {"rsid": "rs10246939", "chrom": "7", "pos": 141974933, "ref": "G", "alt": "A"},
            ]
            return

        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)

        self.snp_defs = data.get("snps", [])
        hap_def = data.get("haplotype_definitions", {})
        self.haplotype_defs = {
            k: {sk.split("_")[0]: sk.split("_")[1] for sk in v}
            for k, v in hap_def.items()
            if not k.endswith("note") and isinstance(v, list)
        }
        raw_diplotypes = data.get("diplotype_phenotypes", _DEFAULT_DIPLOTYPE_PHENOTYPES.copy())
        # 对 diplotype 键进行排序，确保 "PAV/AVI" 和 "AVI/PAV" 都映射到同一键
        self.diplotype_phenotypes = {}
        for key, value in raw_diplotypes.items():
            if key == "default":
                self.diplotype_phenotypes[key] = value
            else:
                sorted_key = "/".join(sorted(key.split("/")))
                self.diplotype_phenotypes[sorted_key] = value
        logger.info("Loaded TAS2R38 SNP definitions: %d SNPs", len(self.snp_defs))

    def analyze(self, tas2r38_variants: List[Variant]) -> TAS2R38Result:
        """分析 TAS2R38 变异，判定 diplotype.

        Args:
            tas2r38_variants: TAS2R38 基因上的变异列表。

        Returns:
            TAS2R38Result 分析结果。
        """
        # 建立位置 -> 变异映射
        pos_map: Dict[int, Variant] = {v.pos: v for v in tas2r38_variants}

        snp_gts: Dict[str, str] = {}
        for snp in self.snp_defs:
            rsid = snp["rsid"]
            pos = int(snp["pos"])
            variant = pos_map.get(pos)
            if variant:
                snp_gts[rsid] = variant.gt
            else:
                snp_gts[rsid] = "./."

        hap1, hap2 = self._call_haplotypes(snp_gts)
        diplotype = f"{hap1}/{hap2}"

        # 标准化 diplotype（排序）
        sorted_diplotype = "/".join(sorted(diplotype.split("/")))

        pheno = self.diplotype_phenotypes.get(
            sorted_diplotype,
            self.diplotype_phenotypes.get("default", {"phenotype": "未知", "level": "未知"}),
        )

        return TAS2R38Result(
            rs713598_gt=snp_gts.get("rs713598", "./."),
            rs1726866_gt=snp_gts.get("rs1726866", "./."),
            rs10246939_gt=snp_gts.get("rs10246939", "./."),
            diplotype=diplotype,
            phenotype_zh=pheno.get("phenotype", "未知"),
            phenotype_level=pheno.get("level", "未知"),
        )

    def _call_haplotypes(self, snp_gts: Dict[str, str]) -> tuple:
        """根据三个 SNP 的基因型调用两个单体型.

        Args:
            snp_gts: rsid -> gt 的字典。

        Returns:
            (haplotype1, haplotype2) 元组。
        """
        # 提取每个 SNP 的等位基因
        alleles_per_snp: Dict[str, List[str]] = {}
        for snp in self.snp_defs:
            rsid = snp["rsid"]
            gt = snp_gts.get(rsid, "./.")
            if "/" in gt:
                alleles = gt.split("/")
            elif "|" in gt:
                alleles = gt.split("|")
            else:
                alleles = [gt, gt]
            # 将 0/1 映射到 ref/alt
            ref = snp["ref"]
            alt = snp["alt"]
            mapped = []
            for a in alleles:
                if a == "0":
                    mapped.append(ref)
                elif a == "1":
                    mapped.append(alt)
                else:
                    mapped.append(".")
            alleles_per_snp[rsid] = mapped

        # 尝试推断单体型
        hap1_alleles = []
        hap2_alleles = []
        for snp in self.snp_defs:
            rsid = snp["rsid"]
            alleles = alleles_per_snp.get(rsid, [".", "."])
            hap1_alleles.append(alleles[0])
            hap2_alleles.append(alleles[1])

        hap1 = self._match_haplotype(hap1_alleles)
        hap2 = self._match_haplotype(hap2_alleles)

        # 如果有一个未匹配，尝试交换相位
        if hap1 == "?" or hap2 == "?":
            hap1_swapped = self._match_haplotype(hap2_alleles)
            hap2_swapped = self._match_haplotype(hap1_alleles)
            if hap1_swapped != "?" or hap2_swapped != "?":
                hap1, hap2 = hap1_swapped, hap2_swapped

        if hap1 == "?":
            hap1 = "未知"
        if hap2 == "?":
            hap2 = "未知"

        return hap1, hap2

    def _match_haplotype(self, alleles: List[str]) -> str:
        """根据等位基因列表匹配单体型.

        Args:
            alleles: 三个 SNP 的等位基因列表（顺序与 snp_defs 一致）。

        Returns:
            单体型名称，若无匹配返回 "?" .
        """
        rsids = [snp["rsid"] for snp in self.snp_defs]
        for hap_name, hap_alleles in self.haplotype_defs.items():
            match = True
            for i, rsid in enumerate(rsids):
                expected = hap_alleles.get(rsid, "")
                if expected and alleles[i] != expected:
                    match = False
                    break
            if match:
                return hap_name
        return "?"
