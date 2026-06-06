"""个人特征定性预测模块.

基于关键 SNP 基因型推断 + 基因层面变异数据，生成瞳孔颜色、毛发颜色、
皮肤颜色、酒精/咖啡因/乳糖代谢等定性预测。
"""

from typing import Dict, List, Optional

from src.logger import get_logger
from src.models import GeneCard, KeySNPResult, PersonalTraitPrediction

logger = get_logger(__name__)


class PersonalTraitPredictor:
    """个人特征定性预测器."""

    def __init__(
        self,
        key_snp_results: Optional[List[KeySNPResult]] = None,
        gene_cards: Optional[List[GeneCard]] = None,
    ):
        self.key_snps = {s.rsid: s for s in (key_snp_results or [])}
        self.gene_cards = {c.gene_symbol: c for c in (gene_cards or [])}

    def predict_all(self) -> List[PersonalTraitPrediction]:
        """生成所有个人特征的定性预测."""
        predictions = []
        predictions.append(self._predict_eye_color())
        predictions.append(self._predict_hair_color())
        predictions.append(self._predict_skin_color())
        predictions.append(self._predict_caffeine_metabolism())
        predictions.append(self._predict_alcohol_metabolism())
        predictions.append(self._predict_lactose_tolerance())
        return [p for p in predictions if p is not None]

    def _get_snp(self, rsid: str) -> Optional[KeySNPResult]:
        return self.key_snps.get(rsid)

    def _get_gene(self, gene: str) -> Optional[GeneCard]:
        return self.gene_cards.get(gene)

    def _predict_eye_color(self) -> Optional[PersonalTraitPrediction]:
        """预测瞳孔颜色.

        核心逻辑：
        - HERC2 rs12913832：G 等位基因（AG/GG）→ 蓝眼倾向；AA → 棕眼倾向
        - OCA2：功能缺失 → 蓝/绿眼倾向
        - SLC45A2 rs16891982：G（GG/AG）→ 浅色虹膜倾向
        - SLC24A4：功能缺失 → 蓝眼倾向
        - SLC24A5：功能缺失 → 浅色眼倾向
        - TYRP1：功能缺失 → 蓝/绿眼倾向
        - IRF4：功能缺失 → 浅色眼倾向
        """
        herc2 = self._get_snp("rs12913832")
        slc45 = self._get_snp("rs16891982")
        oca2 = self._get_gene("OCA2")
        slc24a4 = self._get_gene("SLC24A4")
        slc24a5 = self._get_gene("SLC24A5")
        tyrp1 = self._get_gene("TYRP1")
        irf4 = self._get_gene("IRF4")

        blue_score = 0
        brown_score = 0
        evidence = []
        key_genes = []
        key_snps = []

        # HERC2 是关键调控因子
        if herc2:
            key_snps.append(herc2.rsid)
            key_genes.append("HERC2")
            gt = herc2.inferred_genotype
            if "G" in gt:
                blue_score += 3
                evidence.append(f"HERC2 rs12913832 {gt} → 蓝眼风险等位基因")
            else:
                brown_score += 3
                evidence.append(f"HERC2 rs12913832 {gt} → 棕眼参考等位基因（东亚人群常见）")

        # OCA2 功能状态
        if oca2:
            key_genes.append("OCA2")
            if oca2.assessment.level in ("完全丧失", "显著影响"):
                blue_score += 2
                evidence.append("OCA2 检出功能缺失变异 → 蓝/绿眼倾向")
            elif oca2.assessment.level == "部分影响":
                blue_score += 1
                evidence.append("OCA2 检出部分功能影响变异")

        # SLC45A2
        if slc45:
            key_snps.append(slc45.rsid)
            key_genes.append("SLC45A2")
            gt = slc45.inferred_genotype
            if "G" in gt:
                blue_score += 1
                evidence.append(f"SLC45A2 rs16891982 {gt} → 浅色虹膜倾向")
            else:
                brown_score += 1
                evidence.append(f"SLC45A2 rs16891982 {gt} → 深色虹膜倾向（东亚参考型）")

        # SLC24A4
        if slc24a4 and slc24a4.assessment.level in ("完全丧失", "显著影响"):
            key_genes.append("SLC24A4")
            blue_score += 1
            evidence.append("SLC24A4 功能缺失 → 蓝眼倾向")

        # SLC24A5
        if slc24a5 and slc24a5.assessment.level in ("完全丧失", "显著影响"):
            key_genes.append("SLC24A5")
            blue_score += 1
            evidence.append("SLC24A5 功能缺失 → 浅色眼倾向")

        # TYRP1
        if tyrp1 and tyrp1.assessment.level in ("完全丧失", "显著影响"):
            key_genes.append("TYRP1")
            blue_score += 1
            evidence.append("TYRP1 功能缺失 → 蓝/绿眼倾向")

        # IRF4
        if irf4 and irf4.assessment.level in ("完全丧失", "显著影响"):
            key_genes.append("IRF4")
            blue_score += 1
            evidence.append("IRF4 功能缺失 → 浅色眼倾向")

        if blue_score == 0 and brown_score == 0:
            # 没有任何数据时，根据人群频率推断
            prediction = "深色瞳孔（棕色系）—— 东亚人群常见表型"
            confidence = "中"
            evidence.append("未检出蓝眼风险等位基因，参考东亚人群分布推断")
        elif blue_score > brown_score + 1:
            prediction = "浅色瞳孔倾向（蓝/绿/浅棕色）"
            confidence = "中"
        elif brown_score > blue_score + 1:
            prediction = "深色瞳孔（棕色系）—— 东亚典型表型"
            confidence = "中高"
        else:
            prediction = "中等深浅瞳孔（深棕/ Hazel）"
            confidence = "中"

        return PersonalTraitPrediction(
            trait="瞳孔颜色",
            subsystem="pigmentation",
            prediction=prediction,
            confidence=confidence,
            evidence="；".join(evidence) if evidence else "基于人群参考等位基因推断",
            key_genes=list(set(key_genes)),
            key_snps=list(set(key_snps)),
        )

    def _predict_hair_color(self) -> Optional[PersonalTraitPrediction]:
        """预测毛发颜色.

        核心逻辑：
        - SLC45A2 rs16891982：GG/AG → 浅色毛发；AA → 深色毛发
        - SLC24A5：功能缺失 → 浅色毛发
        - OCA2：功能缺失 → 浅色毛发
        - TYRP1：功能缺失 → 浅色毛发
        - IRF4：功能缺失 → 浅色毛发
        - EDAR：功能变异 → 东亚典型直黑发
        """
        slc45 = self._get_snp("rs16891982")
        slc24a5 = self._get_gene("SLC24A5")
        oca2 = self._get_gene("OCA2")
        tyrp1 = self._get_gene("TYRP1")
        irf4 = self._get_gene("IRF4")
        edar = self._get_gene("EDAR")

        light_score = 0
        dark_score = 0
        evidence = []
        key_genes = []
        key_snps = []

        if slc45:
            key_snps.append(slc45.rsid)
            key_genes.append("SLC45A2")
            gt = slc45.inferred_genotype
            if "G" in gt:
                light_score += 2
                evidence.append(f"SLC45A2 {gt} → 浅色毛发倾向")
            else:
                dark_score += 2
                evidence.append(f"SLC45A2 {gt} → 深色毛发倾向（东亚参考型）")

        if slc24a5 and slc24a5.assessment.level in ("完全丧失", "显著影响"):
            key_genes.append("SLC24A5")
            light_score += 1
            evidence.append("SLC24A5 功能缺失 → 浅色毛发")

        if oca2 and oca2.assessment.level in ("完全丧失", "显著影响"):
            key_genes.append("OCA2")
            light_score += 1
            evidence.append("OCA2 功能缺失 → 浅色毛发")

        if tyrp1 and tyrp1.assessment.level in ("完全丧失", "显著影响"):
            key_genes.append("TYRP1")
            light_score += 1
            evidence.append("TYRP1 功能缺失 → 浅色毛发")

        if irf4 and irf4.assessment.level in ("完全丧失", "显著影响"):
            key_genes.append("IRF4")
            light_score += 1
            evidence.append("IRF4 功能缺失 → 浅色毛发")

        if edar:
            key_genes.append("EDAR")
            # EDAR 的东亚等位基因与直黑发相关，但 daughter 数据以杂合/参考为主
            evidence.append("EDAR 基因存在，与毛发形态相关")

        if light_score == 0 and dark_score == 0:
            prediction = "深色毛发（黑色/深棕色）—— 东亚人群常见表型"
            confidence = "中"
            evidence.append("未检出浅色毛发风险等位基因")
        elif light_score > dark_score + 1:
            prediction = "浅色毛发倾向（金发/浅棕/红发）"
            confidence = "中"
        elif dark_score > light_score + 1:
            prediction = "深色毛发（黑色/深棕色）"
            confidence = "中高"
        else:
            prediction = "中等深浅毛发（深棕/棕色）"
            confidence = "中"

        return PersonalTraitPrediction(
            trait="毛发颜色",
            subsystem="pigmentation",
            prediction=prediction,
            confidence=confidence,
            evidence="；".join(evidence) if evidence else "基于人群参考等位基因推断",
            key_genes=list(set(key_genes)),
            key_snps=list(set(key_snps)),
        )

    def _predict_skin_color(self) -> Optional[PersonalTraitPrediction]:
        """预测皮肤颜色.

        核心逻辑：
        - SLC45A2 rs16891982：GG/AG → 浅色皮肤；AA → 深色皮肤
        - SLC24A5：功能缺失 → 浅色皮肤（最强效应）
        - OCA2：功能缺失 → 浅色皮肤
        - HERC2：通过调控 OCA2 影响色素沉着
        """
        slc45 = self._get_snp("rs16891982")
        slc24a5 = self._get_gene("SLC24A5")
        oca2 = self._get_gene("OCA2")
        herc2 = self._get_gene("HERC2")

        light_score = 0
        dark_score = 0
        evidence = []
        key_genes = []
        key_snps = []

        if slc45:
            key_snps.append(slc45.rsid)
            key_genes.append("SLC45A2")
            gt = slc45.inferred_genotype
            if "G" in gt:
                light_score += 2
                evidence.append(f"SLC45A2 {gt} → 浅色皮肤倾向")
            else:
                dark_score += 2
                evidence.append(f"SLC45A2 {gt} → 深色皮肤倾向（东亚参考型）")

        if slc24a5 and slc24a5.assessment.level in ("完全丧失", "显著影响"):
            key_genes.append("SLC24A5")
            light_score += 2
            evidence.append("SLC24A5 功能缺失 → 显著浅色皮肤倾向")
        elif slc24a5 and slc24a5.assessment.level == "部分影响":
            light_score += 1
            evidence.append("SLC24A5 部分功能影响")

        if oca2 and oca2.assessment.level in ("完全丧失", "显著影响"):
            key_genes.append("OCA2")
            light_score += 1
            evidence.append("OCA2 功能缺失 → 浅色皮肤")

        if herc2 and herc2.assessment.level in ("完全丧失", "显著影响"):
            key_genes.append("HERC2")
            light_score += 1
            evidence.append("HERC2 功能缺失 → 通过 OCA2 调控减轻色素沉着")

        if light_score == 0 and dark_score == 0:
            prediction = "中等偏深肤色（东亚典型肤色）"
            confidence = "中"
            evidence.append("未检出显著浅色皮肤风险等位基因")
        elif light_score > dark_score + 1:
            prediction = "浅色皮肤倾向"
            confidence = "中"
        elif dark_score > light_score + 1:
            prediction = "深色皮肤（东亚典型肤色）"
            confidence = "中高"
        else:
            prediction = "中等肤色"
            confidence = "中"

        return PersonalTraitPrediction(
            trait="皮肤颜色",
            subsystem="pigmentation",
            prediction=prediction,
            confidence=confidence,
            evidence="；".join(evidence) if evidence else "基于人群参考等位基因推断",
            key_genes=list(set(key_genes)),
            key_snps=list(set(key_snps)),
        )

    def _predict_caffeine_metabolism(self) -> Optional[PersonalTraitPrediction]:
        """预测咖啡因代谢能力.

        核心逻辑：
        - CYP1A2 rs762551：AA → 快代谢；AC/CC → 慢代谢
        - CYP1A2 基因层面功能缺失 → 慢代谢
        """
        cyp1a2_snp = self._get_snp("rs762551")
        cyp1a2_gene = self._get_gene("CYP1A2")

        evidence = []
        key_genes = []
        key_snps = []

        if cyp1a2_snp:
            key_snps.append(cyp1a2_snp.rsid)
            key_genes.append("CYP1A2")
            gt = cyp1a2_snp.inferred_genotype
            if gt == "AA":
                prediction = "咖啡因快代谢型 — 咖啡因耐受度高，可较快清除"
                confidence = "中高"
                evidence.append(f"CYP1A2 rs762551 {gt} → 快代谢等位基因纯合")
            elif "C" in gt:
                prediction = "咖啡因慢代谢型 — 咖啡因敏感，建议控制摄入量"
                confidence = "中高"
                evidence.append(f"CYP1A2 rs762551 {gt} → 慢代谢等位基因携带")
            else:
                prediction = "咖啡因代谢能力中等"
                confidence = "中"
                evidence.append(f"CYP1A2 rs762551 {gt}")
        elif cyp1a2_gene:
            key_genes.append("CYP1A2")
            if cyp1a2_gene.assessment.level in ("完全丧失", "显著影响"):
                prediction = "咖啡因慢代谢型 — CYP1A2 功能显著受损"
                confidence = "中"
                evidence.append("CYP1A2 检出功能显著受损变异")
            else:
                prediction = "咖啡因代谢能力正常范围"
                confidence = "中"
                evidence.append("CYP1A2 未检出关键功能变异")
        else:
            prediction = "咖啡因代谢能力正常范围（参考型）"
            confidence = "低"
            evidence.append("CYP1A2 数据不足")

        return PersonalTraitPrediction(
            trait="咖啡因代谢",
            subsystem="metabolism",
            prediction=prediction,
            confidence=confidence,
            evidence="；".join(evidence),
            key_genes=list(set(key_genes)),
            key_snps=list(set(key_snps)),
        )

    def _predict_alcohol_metabolism(self) -> Optional[PersonalTraitPrediction]:
        """预测酒精代谢能力.

        核心逻辑：
        - ALDH2 rs671：A（AA/AG）→ 乙醛脱氢酶活性降低 → 酒精不耐受/脸红
        - ADH1B rs1229984：T（TT/CT）→ 乙醇脱氢酶活性增强 → 代谢加快
        - 两者同时存在时：乙醛积累风险增加
        """
        aldh2_snp = self._get_snp("rs671")
        aldh2_gene = self._get_gene("ALDH2")
        adh1b = self._get_gene("ADH1B")

        evidence = []
        key_genes = []
        key_snps = []

        # 优先使用 SNP 数据
        if aldh2_snp:
            key_snps.append(aldh2_snp.rsid)
            key_genes.append("ALDH2")
            gt = aldh2_snp.inferred_genotype
            if "A" in gt:
                prediction = "酒精代谢能力降低 — 乙醛脱氢酶活性受损，饮酒易脸红"
                confidence = "中高"
                evidence.append(f"ALDH2 rs671 {gt} → 乙醛脱氢酶活性降低等位基因")
            else:
                prediction = "酒精代谢能力正常（ALDH2 酶活性正常）"
                confidence = "中高"
                evidence.append(f"ALDH2 rs671 {gt} → 正常活性等位基因纯合")
        elif aldh2_gene:
            key_genes.append("ALDH2")
            if aldh2_gene.assessment.level in ("完全丧失", "显著影响"):
                prediction = "酒精代谢能力降低 — 乙醛脱氢酶活性可能受损，饮酒易脸红"
                confidence = "中"
                evidence.append("ALDH2 检出功能显著受损变异 → 乙醛积累风险")
            elif aldh2_gene.assessment.level == "部分影响":
                prediction = "酒精代谢能力可能轻度降低"
                confidence = "中"
                evidence.append("ALDH2 检出部分功能影响变异")
            else:
                prediction = "酒精代谢能力正常（ALDH2 参考型）"
                confidence = "中"
                evidence.append("ALDH2 未检出显著功能变异")
        else:
            prediction = "酒精代谢能力正常范围（参考型推断）"
            confidence = "低"
            evidence.append("ALDH2 基因数据不足")

        if adh1b:
            key_genes.append("ADH1B")
            if adh1b.assessment.level in ("完全丧失", "显著影响"):
                evidence.append("ADH1B 检出功能显著受损变异 → 乙醇代谢可能受影响")
            else:
                evidence.append("ADH1B 未检出显著功能变异")

        # 综合判断
        if "降低" in prediction:
            pass  # 保持上面的判断
        elif adh1b and adh1b.assessment.level in ("完全丧失", "显著影响"):
            prediction = "酒精代谢可能异常 — 乙醇/乙醛代谢链存在变异"
            confidence = "中"

        return PersonalTraitPrediction(
            trait="酒精代谢",
            subsystem="metabolism",
            prediction=prediction,
            confidence=confidence,
            evidence="；".join(evidence),
            key_genes=list(set(key_genes)),
            key_snps=list(set(key_snps)),
        )

    def _predict_lactose_tolerance(self) -> Optional[PersonalTraitPrediction]:
        """预测乳糖耐受能力.

        核心逻辑：
        - LCT/MCM6 rs4988235：T（TT/CT）→ 乳糖耐受；CC → 乳糖不耐受
        - LCT 基因层面功能缺失 → 乳糖不耐受
        """
        lct_snp = self._get_snp("rs4988235")
        lct_gene = self._get_gene("LCT")

        evidence = []
        key_genes = []
        key_snps = []

        # 优先使用 SNP 数据
        if lct_snp:
            key_snps.append(lct_snp.rsid)
            key_genes.append("LCT")
            gt = lct_snp.inferred_genotype
            if "A" in gt:
                prediction = "乳糖耐受 — 乳糖酶持久型，可正常消化乳制品"
                confidence = "中高"
                evidence.append(f"LCT rs4988235 {gt} → 乳糖酶持久型等位基因携带")
            else:
                prediction = "乳糖不耐受倾向 — 乳糖酶非持久型，成年后乳糖酶活性可能下降"
                confidence = "中高"
                evidence.append(f"LCT rs4988235 {gt} → 乳糖酶非持久型参考等位基因（东亚常见）")
        elif lct_gene:
            key_genes.append("LCT")
            if lct_gene.assessment.level in ("完全丧失", "显著影响"):
                prediction = "乳糖不耐受倾向 — LCT 功能可能受损，建议控制乳制品摄入"
                confidence = "中"
                evidence.append("LCT 检出功能显著受损变异 → 乳糖酶活性可能降低")
            elif lct_gene.assessment.level == "部分影响":
                prediction = "乳糖耐受能力可能轻度降低"
                confidence = "中"
                evidence.append("LCT 检出部分功能影响变异")
            else:
                prediction = "乳糖耐受能力正常（LCT 参考型）"
                confidence = "中"
                evidence.append("LCT 未检出显著功能变异")
        else:
            prediction = "乳糖耐受能力正常范围（参考型推断）"
            confidence = "低"
            evidence.append("LCT 基因数据不足")

        return PersonalTraitPrediction(
            trait="乳糖耐受",
            subsystem="metabolism",
            prediction=prediction,
            confidence=confidence,
            evidence="；".join(evidence),
            key_genes=list(set(key_genes)),
            key_snps=list(set(key_snps)),
        )
