"""OR 基因分级展示逻辑.

Tier A: 已知配体 + 纯合 LoF/显著错义 → 正文详细展示
Tier B: 功能研究但配体未知 + 纯合 LoF → 表格列出
Tier C: 其余纯合 LoF → 附录折叠
"""

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from src.assessment.engine import ImpactEngine
from src.gene_sets.loader import GeneSetLoader
from src.logger import get_logger
from src.models import ImpactAssessment, ORTierResult, Sex, Variant

logger = get_logger(__name__)


class ORTierClassifier:
    """OR 基因分级分类器."""

    def __init__(self, data_path: Optional[str] = None) -> None:
        """初始化分类器.

        Args:
            data_path: or_ligands.json 路径。
        """
        if data_path is None:
            src_dir = Path(__file__).resolve().parent.parent
            data_path = src_dir.parent / "data" / "or_ligands.json"
        else:
            data_path = Path(data_path)

        self.ligand_map: Dict[str, Dict[str, str]] = {}
        self._load_data(data_path)

    def _load_data(self, path: Path) -> None:
        """加载 OR 配体映射表."""
        if not path.exists():
            logger.warning("or_ligands.json not found, using empty defaults")
            return

        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)

        for entry in data.get("ligands", []):
            gene = entry.get("gene_symbol", "")
            if gene:
                self.ligand_map[gene] = entry

        logger.info("Loaded OR ligand map: %d entries", len(self.ligand_map))

    def classify(
        self, or_variants: List[Variant], sex: Sex, gene_sets: GeneSetLoader
    ) -> List[ORTierResult]:
        """对 OR 基因变异进行分级分类（按基因维度聚合）.

        每个 OR 基因只保留最高 tier 的代表性变异，避免变异级别膨胀。

        Args:
            or_variants: OR 基因变异列表。
            sex: 样本性别。
            gene_sets: 感官基因集加载器。

        Returns:
            ORTierResult 列表（按基因去重）。
        """
        impact_engine = ImpactEngine()

        # 按基因分组聚合变异
        gene_variants: Dict[str, List[Variant]] = {}
        for variant in or_variants:
            gene = variant.gene_symbol
            if not gene or not gene.startswith("OR"):
                continue
            gene_variants.setdefault(gene, []).append(variant)

        results = []
        for gene, variants in gene_variants.items():
            ligand_info = self.ligand_map.get(gene, {})
            known_ligand = ligand_info.get("ligand_zh")

            # 评估该基因的所有变异，取最高 tier
            best_tier = "C"
            best_variant = variants[0]
            best_assessment = impact_engine.assess(best_variant, sex, gene_sets)

            for variant in variants:
                assessment = impact_engine.assess(variant, sex, gene_sets)
                tier = self._determine_tier(variant, assessment, known_ligand)
                tier_rank = {"A": 0, "B": 1, "C": 2}
                if tier_rank.get(tier, 99) < tier_rank.get(best_tier, 99):
                    best_tier = tier
                    best_variant = variant
                    best_assessment = assessment

            result = ORTierResult(
                tier=best_tier,
                gene_symbol=gene,
                known_ligand_zh=known_ligand,
                odor_description_zh=ligand_info.get("odor_description_zh"),
                variant=best_variant,
                assessment=best_assessment,
            )
            results.append(result)

        # 按 tier 排序：A > B > C
        tier_order = {"A": 0, "B": 1, "C": 2}
        results.sort(key=lambda x: tier_order.get(x.tier, 99))

        return results

    def _determine_tier(
        self, variant: Variant, assessment: ImpactAssessment, known_ligand: Optional[str]
    ) -> str:
        """确定 OR 基因的分级.

        Tier A: 已知配体 + 任何功能变异（纯合/杂合 LoF 或错义）
        Tier A+: 已知配体 + 杂合错义（单独标记展示）
        Tier B: 功能研究但配体未知 + 纯合 LoF
        Tier C: 其余不展示

        Args:
            variant: 变异记录。
            assessment: 功能影响评估。
            known_ligand: 已知配体（中文）。

        Returns:
            Tier 等级 "A" / "B" / "C"。
        """
        consequence = variant.consequence or ""
        is_lof = "frameshift" in consequence or "stop_gained" in consequence
        is_missense = "missense" in consequence

        # Tier A: 已知配体 + 纯合 LoF/错义
        if known_ligand and variant.is_homozygous:
            if is_lof or is_missense:
                return "A"

        # Tier A: 已知配体 + 杂合错义（降低展示优先级但仍展示）
        if known_ligand and not variant.is_homozygous and is_missense:
            return "A"

        # Tier B: 功能研究但配体未知 + 纯合 LoF
        if not known_ligand and variant.is_homozygous and is_lof:
            return "B"

        # Tier C: 其余纯合 LoF
        if variant.is_homozygous and is_lof:
            return "C"

        # 默认不展示
        return "C"
