"""Markdown 报告生成器.

使用 Jinja2 模板引擎渲染人类可读的 Markdown 报告。
"""

import os
from pathlib import Path
from typing import Any, Dict, List, Optional

from jinja2 import Environment, FileSystemLoader, select_autoescape

from src.logger import get_logger
from src.models import GeneCard, SensoryReport, Variant

logger = get_logger(__name__)


class MarkdownReportGenerator:
    """Markdown 报告生成器."""

    def __init__(self, templates_dir: Optional[str] = None) -> None:
        """初始化生成器.

        Args:
            templates_dir: 模板目录路径，默认查找 ../templates/。
        """
        if templates_dir is None:
            # 模板与 src 代码同目录下的 templates/
            src_dir = Path(__file__).resolve().parent.parent
            templates_dir = src_dir / "templates"
        else:
            templates_dir = Path(templates_dir)

        self.templates_dir = Path(templates_dir)
        self.jinja_env = Environment(
            loader=FileSystemLoader(str(self.templates_dir)),
            autoescape=select_autoescape(["html", "xml"]),
            trim_blocks=True,
            lstrip_blocks=True,
        )

        # 注册自定义过滤器
        self.jinja_env.filters["level_badge"] = self._level_badge_filter
        self.jinja_env.filters["risk_badge"] = self._risk_badge_filter
        self.jinja_env.filters["sift_badge"] = self._sift_badge_filter
        self.jinja_env.filters["polyphen_badge"] = self._polyphen_badge_filter
        self.jinja_env.filters["freq_desc"] = self._freq_desc_filter
        self.jinja_env.filters["key_variants"] = self._key_variants_filter
        # Register 'in' test for selectattr/rejectattr
        self.jinja_env.tests["in"] = lambda x, y: x in y

    def generate(self, report: SensoryReport) -> str:
        """生成 Markdown 报告.

        Args:
            report: 报告数据模型。

        Returns:
            Markdown 格式报告字符串。
        """
        template = self.jinja_env.get_template("report.md.j2")
        context = self._build_context(report)
        return template.render(**context)

    def _build_context(self, report: SensoryReport) -> Dict[str, Any]:
        """构建模板上下文."""
        # 构建各子系统的基因卡片列表
        subsys_cards: Dict[str, List[GeneCard]] = {}
        for card in report.gene_cards:
            ss = card.subsystem or "unknown"
            subsys_cards.setdefault(ss, []).append(card)

        # 构建关键发现列表（用于执行摘要）
        key_findings = []
        for card in report.gene_cards:
            level = card.assessment.level
            if level in ("完全丧失", "显著影响", "部分影响"):
                # 收集关键变异信息
                top_vars = self._get_top_variants(card.variants, n=2)
                var_info = ", ".join([
                    f"{v.chrom}:{v.pos} {v.consequence[:20]}" if v.consequence else f"{v.chrom}:{v.pos}"
                    for v in top_vars
                ])
                key_findings.append({
                    "gene": card.gene_symbol,
                    "subsystem": card.subsystem,
                    "level": level,
                    "rationale": card.assessment.rationale_zh[:80] if card.assessment.rationale_zh else "",
                    "variants": var_info,
                })

        # 构建综合特征档案
        profile = self._build_comprehensive_profile(report)

        # 构建 gnomAD 频率参考表
        gnomad_refs = []
        for card in report.gene_cards:
            for v in card.variants:
                if v.rsid and (v.gnomad_af_exome is not None or v.gnomad_af_genome is not None):
                    gnomad_refs.append({
                        "gene": card.gene_symbol,
                        "variant": f"{v.chrom}:{v.pos} {v.ref}>{v.alt}",
                        "rsid": v.rsid,
                        "exome_af": v.gnomad_af_exome,
                        "genome_af": v.gnomad_af_genome,
                        "flags": v.lof_flags,
                    })
                    if len(gnomad_refs) >= 30:
                        break
            if len(gnomad_refs) >= 30:
                break

        # OR tier 分层
        or_tier_a = []
        or_tier_b = []
        or_tier_c = []
        if report.or_tiers:
            for ot in report.or_tiers:
                if ot.tier == "A":
                    or_tier_a.append(ot)
                elif ot.tier == "B":
                    or_tier_b.append(ot)
                else:
                    or_tier_c.append(ot)

        return {
            "report": report,
            "sample_id": report.sample_id,
            "sex": report.sex,
            "ref_genome": report.ref_genome,
            "analysis_date": report.analysis_date.strftime("%Y-%m-%d %H:%M:%S"),
            "subsystems": report.subsystems,
            "gene_cards": report.gene_cards,
            # 按影响级别分层：有影响 (>=部分影响) → 主报告，其他 → 附录
            "impactful_cards": [c for c in report.gene_cards 
                if c.assessment.level in ("完全丧失", "显著影响", "部分影响")],
            "mild_cards": [c for c in report.gene_cards 
                if c.assessment.level in ("可能轻微影响", "无影响")],
            "subsys_cards": subsys_cards,
            "tas2r38": report.tas2r38,
            "mitochondrial": report.mitochondrial,
            "or_tiers": report.or_tiers,
            "or_tier_a": or_tier_a,
            "or_tier_b": or_tier_b,
            "or_tier_c": or_tier_c,
            "executive_summary": report.executive_summary,
            "personal_traits": report.executive_summary.personal_traits,
            "disclaimer_zh": report.disclaimer_zh,
            "data_availability": report.data_availability,
            "key_findings": key_findings,
            "profile": profile,
            "gnomad_refs": gnomad_refs,
            "key_snps": report.key_snps if report.key_snps else [],
            # SNP hits: key SNPs actually found in VCF
            "snp_hits": [s for s in (report.key_snps or []) if s.found_in_vcf],
            # High-impact cards: significant or above
            "high_impact_cards": [c for c in report.gene_cards
                if c.assessment.level in ("完全丧失", "显著影响")],
            # LOF/GOF variants across all genes
            "lof_gof_variants": self._collect_lof_gof_variants(report.gene_cards),
            # All analyzed genes list for appendix
            "all_genes_list": sorted(report.gene_cards, key=lambda c: c.gene_symbol),
        }

    @staticmethod
    def _collect_lof_gof_variants(gene_cards: List[GeneCard]) -> List[Dict[str, Any]]:
        """收集所有 LOF/GOF 变异."""
        lof_gof = []
        for card in gene_cards:
            for v in card.variants:
                cons = v.consequence.lower() if v.consequence else ""
                is_lof = any(x in cons for x in [
                    "frameshift_variant", "stop_gained", "stop_lost", "start_lost",
                    "splice_acceptor_variant", "splice_donor_variant", "transcript_ablation"
                ])
                # GOF is harder to predict; flag missense with high CADD or known domain
                is_gof_candidate = (
                    "missense" in cons
                    and (v.cadd_score is not None and v.cadd_score > 20)
                )
                if is_lof or is_gof_candidate:
                    lof_gof.append({
                        "gene": card.gene_symbol,
                        "subsystem": card.subsystem,
                        "variant": v,
                        "type": "LOF" if is_lof else "GOF候选",
                    })
        # Sort: LOF first, then by consequence severity
        def sort_key(item: Dict[str, Any]) -> int:
            cons = item["variant"].consequence.lower() if item["variant"].consequence else ""
            if "frameshift" in cons:
                return 0
            if "stop_gained" in cons or "stop_lost" in cons or "start_lost" in cons:
                return 1
            if "splice_acceptor" in cons or "splice_donor" in cons:
                return 2
            if "missense" in cons:
                return 3
            return 4
        return sorted(lof_gof, key=sort_key)

    @staticmethod
    def _get_top_variants(variants: List[Variant], n: int = 3) -> List[Variant]:
        """从变异列表中选出最重要的 n 个变异.

        优先级：1) 后果严重度（frameshift > stop_gained > splice > missense > synonymous）
               2) 纯合 > 杂合
               3) 有 SIFT/PolyPhen 有害预测
               4) 质量值
        """
        if not variants:
            return []

        consequence_rank = {
            "frameshift": 6, "stop_gained": 5, "splice": 4,
            "missense": 3, "inframe": 2, "synonymous": 1, "": 0,
        }

        def score(v: Variant) -> int:
            s = 0
            cons = v.consequence.lower() if v.consequence else ""
            for key, rank in consequence_rank.items():
                if key in cons:
                    s = max(s, rank * 100)
            if v.is_homozygous:
                s += 50
            if v.sift and "deleterious" in v.sift.lower():
                s += 30
            if v.polyphen and ("probably_damaging" in v.polyphen.lower() or "possibly_damaging" in v.polyphen.lower()):
                s += 20
            s += min(int(v.qual / 10), 20)
            return s

        sorted_vars = sorted(variants, key=score, reverse=True)
        return sorted_vars[:n]

    @staticmethod
    def _build_comprehensive_profile(report: SensoryReport) -> List[Dict[str, Any]]:
        """构建综合特征档案."""
        profile = []
        subsys_names = {
            "vision": "视觉", "hearing": "听觉", "olfaction": "嗅觉",
            "taste": "味觉", "somatosensation": "体感/痛觉",
            "pigmentation": "色素沉着", "metabolism": "代谢特征",
            "muscle": "肌肉特征", "hair": "毛发特征",
        }

        # 按基因预设的特征映射
        feature_map = {
            "OR2B11": ("嗅觉", "花香类嗅觉受体", "⭐⭐"),
            "OR7D4": ("嗅觉", "雄烯酮敏感度", "⭐⭐⭐"),
            "OR6A2": ("嗅觉", "香菜/芫荽偏好", "⭐⭐⭐"),
            "TAS2R38": ("味觉", "苦味敏感度", "⭐⭐⭐"),
            "TAS1R2": ("味觉", "甜味感知", "⭐⭐"),
            "TAS1R3": ("味觉", "甜味感知", "⭐⭐"),
            "TAS1R1": ("味觉", "鲜味感知", "⭐⭐"),
            "SCNN1A": ("味觉", "咸味感知", "⭐⭐"),
            "OTOP1": ("味觉", "酸味感知", "⭐⭐"),
            "SCN9A": ("体感", "痛觉阈值", "⭐⭐⭐"),
            "SCN10A": ("体感", "痛觉产生", "⭐⭐"),
            "PIEZO2": ("体感", "触觉分辨与本体感觉", "⭐⭐"),
            "TRPV1": ("体感", "辣椒素/热痛耐受", "⭐⭐"),
            "TRPM8": ("体感", "冷感/薄荷醇敏感度", "⭐⭐"),
            "OCA2": ("视觉", "瞳孔颜色", "⭐⭐"),
            "HERC2": ("视觉", "瞳孔颜色", "❌"),
            "TYRP1": ("视觉", "瞳孔颜色", "⭐⭐"),
            "EYS": ("视觉", "夜间视力", "⭐⭐"),
            "OPN1LW": ("视觉", "红绿色觉", "⭐⭐⭐"),
            "OPN1MW": ("视觉", "红绿色觉", "⭐⭐⭐"),
            "OPN1SW": ("视觉", "蓝黄色觉", "⭐⭐⭐"),
            "MYO7A": ("听觉", "毛细胞静纤毛发育", "⭐⭐"),
            "CDH23": ("听觉", "毛细胞尖端链接", "⭐⭐"),
            "GJB2": ("听觉", "缝隙连接蛋白", "⭐⭐⭐"),
            "MT-RNR1": ("听觉", "抗生素耳毒性风险", "⭐⭐⭐⭐⭐"),
            "HERC2": ("色素", "瞳孔颜色", "⭐⭐⭐⭐"),
            "SLC45A2": ("色素", "皮肤/眼色色素深浅", "⭐⭐⭐"),
            "SLC24A4": ("色素", "眼色色素沉着", "⭐⭐"),
            "SLC24A5": ("色素", "皮肤色素沉着", "⭐⭐⭐"),
            "IRF4": ("色素", "眼色与毛发色素", "⭐⭐"),
            "CYP1A2": ("代谢", "咖啡因代谢速率", "⭐⭐⭐⭐"),
            "LCT": ("代谢", "乳糖耐受能力", "⭐⭐⭐⭐"),
            "ALDH2": ("代谢", "酒精代谢（乙醛脱氢酶）", "⭐⭐⭐⭐⭐"),
            "ADH1B": ("代谢", "酒精代谢（乙醇脱氢酶）", "⭐⭐⭐⭐"),
            "ACTN3": ("肌肉", "快肌纤维/爆发力", "⭐⭐⭐⭐"),
            "EDAR": ("毛发", "毛发密度与直卷", "⭐⭐⭐"),
        }

        seen_genes = set()
        for card in report.gene_cards:
            gene = card.gene_symbol
            if gene in seen_genes:
                continue
            seen_genes.add(gene)
            if gene in feature_map:
                subsystem, feature, confidence = feature_map[gene]
                # 根据评估等级调整推断
                level = card.assessment.level
                if level == "完全丧失":
                    inference = "功能完全丧失"
                elif level == "显著影响":
                    inference = "可能有显著个体差异"
                elif level == "部分影响":
                    inference = "可能有轻微个体差异"
                elif level == "可能轻微影响":
                    inference = "可能有细微差异"
                else:
                    inference = "正常范围"
                profile.append({
                    "feature": feature,
                    "subsystem": subsystem,
                    "gene": gene,
                    "inference": inference,
                    "confidence": confidence,
                    "level": level,
                })

        # TAS2R38 特殊处理
        if report.tas2r38:
            if report.tas2r38.diplotype:
                for p in profile:
                    if p["gene"] == "TAS2R38":
                        p["inference"] = f"苦味感知 {report.tas2r38.diplotype} → {report.tas2r38.phenotype_zh}"
            else:
                for p in profile:
                    if p["gene"] == "TAS2R38":
                        p["inference"] = "数据不足（关键SNP缺失）"
                        p["confidence"] = "❌"

        # 线粒体特殊处理
        if report.mitochondrial:
            for mt in report.mitochondrial:
                profile.append({
                    "feature": "抗生素耳毒性风险",
                    "subsystem": "听觉",
                    "gene": mt.gene,
                    "inference": mt.drug_warning_zh,
                    "confidence": "⭐⭐⭐⭐⭐",
                    "level": "显著影响",
                })
        else:
            profile.append({
                "feature": "抗生素耳毒性风险",
                "subsystem": "听觉",
                "gene": "MT-RNR1",
                "inference": "无风险（无线粒体致聋变异）",
                "confidence": "⭐⭐⭐⭐⭐",
                "level": "无影响",
            })

        return profile

    @staticmethod
    def _level_badge_filter(level: str) -> str:
        """将评估等级转换为 Markdown badge."""
        badges = {
            "完全丧失": "🔴 完全丧失",
            "显著影响": "🟠 显著影响",
            "部分影响": "🟡 部分影响",
            "可能轻微影响": "🔵 可能轻微影响",
            "无影响": "🟢 无影响",
        }
        return badges.get(level, level)

    @staticmethod
    def _risk_badge_filter(risk: str) -> str:
        """将风险等级转换为 Markdown badge."""
        badges = {
            "高风险": "🔴 高风险",
            "中风险": "🟠 中风险",
            "低风险": "🟢 低风险",
            "未知": "⚪ 未知",
        }
        return badges.get(risk, risk)

    @staticmethod
    def _sift_badge_filter(sift: Optional[str]) -> str:
        """将 SIFT 预测转换为可读文本."""
        if not sift:
            return "N/A"
        s = sift.lower()
        if "deleterious" in s:
            return "**有害**" if "low_confidence" not in s else "有害(低置信)"
        if "tolerated" in s:
            return "可耐受"
        return sift

    @staticmethod
    def _polyphen_badge_filter(polyphen: Optional[str]) -> str:
        """将 PolyPhen 预测转换为可读文本."""
        if not polyphen:
            return "N/A"
        s = polyphen.lower()
        if "probably_damaging" in s:
            return "**可能有害**"
        if "possibly_damaging" in s:
            return "可能有害"
        if "benign" in s:
            return "良性"
        if "unknown" in s:
            return "未知"
        return polyphen

    @staticmethod
    def _freq_desc_filter(af: Optional[float]) -> str:
        """将等位基因频率转换为人群频率描述."""
        if af is None:
            return "未知"
        if af > 0.4:
            return f"常见多态性 ({af*100:.1f}%)"
        if af > 0.05:
            return f"常见 ({af*100:.1f}%)"
        if af > 0.01:
            return f"低频 ({af*100:.1f}%)"
        if af > 0.001:
            return f"罕见 ({af*100:.2f}%)"
        return f"极罕见 ({af*100:.3f}%)"

    @staticmethod
    def _variant_impact_score(variant: Variant) -> int:
        """计算变异的蛋白影响分数（用于排序）."""
        score = 0
        cons = (variant.consequence or "").lower()
        consequence_scores = {
            "frameshift": 100, "stop_gained": 95, "stop_lost": 90,
            "start_lost": 85, "splice": 80, "missense": 50,
            "inframe": 40, "protein_altering": 35,
        }
        for key, val in consequence_scores.items():
            if key in cons:
                score = max(score, val)
        if variant.is_homozygous:
            score += 30
        if variant.sift and "deleterious" in variant.sift.lower():
            score += 20
        if variant.polyphen and "probably_damaging" in variant.polyphen.lower():
            score += 15
        if variant.polyphen and "possibly_damaging" in variant.polyphen.lower():
            score += 10
        if variant.protein_domain:
            score += 5
        score += min(int(variant.qual / 10), 10)
        return score

    @staticmethod
    def _key_variants_filter(variants: List[Variant], n: int = 3) -> List[Variant]:
        """筛选关键变异 — 优先使用已标记的关键变异."""
        # 如果 variants 中已有 is_key_variant 标记，直接使用
        key_vars = [v for v in variants if getattr(v, "is_key_variant", False)]
        if key_vars:
            # 按影响分数排序，取前 n 个
            return sorted(key_vars, key=MarkdownReportGenerator._variant_impact_score, reverse=True)[:n]
        # fallback：从全部变异中筛选
        return MarkdownReportGenerator._get_top_variants(variants, n)
