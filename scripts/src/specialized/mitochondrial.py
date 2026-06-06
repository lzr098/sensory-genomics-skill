"""线粒体耳聋注释器.

识别 MT-RNR1/MT-TS1 已知致聋位点，输出药物风险警告和异质性水平。
"""

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from src.logger import get_logger
from src.models import MitochondrialResult, Variant

logger = get_logger(__name__)


class MitochondrialAnnotator:
    """线粒体耳聋注释器."""

    def __init__(self, data_path: Optional[str] = None) -> None:
        """初始化注释器.

        Args:
            data_path: mitochondrial_variants.json 路径。
        """
        if data_path is None:
            src_dir = Path(__file__).resolve().parent.parent
            data_path = src_dir.parent / "data" / "mitochondrial_variants.json"
        else:
            data_path = Path(data_path)

        self.known_variants: List[Dict[str, Any]] = []
        self._load_data(data_path)

    def _load_data(self, path: Path) -> None:
        """加载已知线粒体致聋变异库."""
        if not path.exists():
            logger.warning("mitochondrial_variants.json not found, using empty defaults")
            return

        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)

        self.known_variants = data.get("variants", [])
        logger.info("Loaded mitochondrial variant database: %d entries", len(self.known_variants))

    def annotate(self, mt_variants: List[Variant]) -> List[MitochondrialResult]:
        """注释线粒体变异.

        Args:
            mt_variants: MT 染色体上的变异列表。

        Returns:
            MitochondrialResult 列表。
        """
        results = []
        for variant in mt_variants:
            matched = self._match_variant(variant)
            if matched:
                heteroplasmy = self._extract_heteroplasmy(variant)
                result = MitochondrialResult(
                    variant_name=matched["variant_name"],
                    gene=matched["gene"],
                    heteroplasmy=heteroplasmy,
                    drug_warning_zh=matched.get("drug_warning", ""),
                    risk_level=matched.get("risk_level", "未知"),
                )
                results.append(result)

        return results

    def _match_variant(self, variant: Variant) -> Optional[Dict[str, Any]]:
        """在已知变异库中查找匹配."""
        for known in self.known_variants:
            if (
                int(known.get("position", 0)) == variant.pos
                and known.get("ref", "") == variant.ref
                and known.get("alt", "") == variant.alt
            ):
                return known
        return None

    @staticmethod
    def _extract_heteroplasmy(variant: Variant) -> float:
        """从变异中提取异质性水平.

        若 VCF 中无 AF/heteroplasmy 字段，返回 1.0（假设为纯质性）。

        Args:
            variant: 变异记录。

        Returns:
            异质性水平（0-1）。
        """
        raw = variant.raw_vep or {}
        # 尝试从 raw_vep 中提取
        if "allele_frequency" in raw:
            try:
                return float(raw["allele_frequency"])
            except (ValueError, TypeError):
                pass
        # 若 GT 为 1（无 / 分隔），可能为纯质性
        if variant.gt == "1":
            return 1.0
        return 1.0
