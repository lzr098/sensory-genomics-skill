"""Markdown 报告生成器.

使用 Jinja2 模板引擎渲染人类可读的 Markdown 报告。
"""

import os
from pathlib import Path
from typing import Any, Dict, List, Optional

from jinja2 import Environment, FileSystemLoader, select_autoescape

from src.enrichment.clinvar import query_clinvar_by_variant
from src.logger import get_logger
from src.models import GeneCard, SensoryReport, Variant

logger = get_logger(__name__)

# ClinVar significance rank (lower = more benign)
_CLINVAR_RANK = {
    "benign": 0,
    "likely_benign": 1,
    "uncertain_significance": 2,
    "vus": 2,
    "likely_pathogenic": 3,
    "pathogenic": 4,
}


def _parse_clinvar_significance(clnsig: str) -> str:
    """Normalize ClinVar CLNSIG to a simple category."""
    if not clnsig or clnsig == ".":
        return "unknown"
    s = clnsig.lower().replace(" ", "_").replace("-", "_")
    if "benign/likely_benign" in s or "benign" in s and "pathogenic" not in s:
        return "benign"
    if "likely_benign" in s:
        return "likely_benign"
    if "pathogenic/likely_pathogenic" in s or "pathogenic" in s and "benign" not in s:
        return "pathogenic"
    if "likely_pathogenic" in s:
        return "likely_pathogenic"
    if "uncertain_significance" in s or "vus" in s:
        return "vus"
    return "unknown"


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
        # Cross-check impactful genes against gnomAD/ClinVar
        cross_checks: Dict[str, Dict[str, Any]] = {}
        for card in report.gene_cards:
            if card.assessment.level in ("完全丧失", "显著影响", "部分影响"):
                cross_checks[card.gene_symbol] = self._cross_check_gene(card)

        # Downgraded gene symbols set
        downgraded_genes = {
            card.gene_symbol for card in report.gene_cards
            if card.assessment.level in ("完全丧失", "显著影响")
            and cross_checks.get(card.gene_symbol, {}).get("downgrade")
        }

        # 构建各子系统的基因卡片列表（排除已降级的显著影响/完全丧失基因）
        subsys_cards: Dict[str, List[GeneCard]] = {}
        for card in report.gene_cards:
            if card.gene_symbol in downgraded_genes:
                continue
            ss = card.subsystem or "unknown"
            subsys_cards.setdefault(ss, []).append(card)

        # 构建关键发现列表（用于执行摘要）
        key_findings = []
        downgraded_findings = []
        for card in report.gene_cards:
            level = card.assessment.level
            if level in ("完全丧失", "显著影响", "部分影响"):
                cc = cross_checks.get(card.gene_symbol, {})

                # For 显著影响/完全丧失: downgrade removes from key_findings
                if level in ("完全丧失", "显著影响") and cc.get("downgrade"):
                    downgraded_findings.append({
                        "gene": card.gene_symbol,
                        "subsystem": card.subsystem,
                        "original_level": level,
                        "downgrade_reasons": cc.get("downgrade_reasons", []),
                        "gnomad_af_max": cc.get("gnomad_af_max", 0.0),
                        "variant_checks": cc.get("variant_checks", []),
                    })
                    continue

                # 收集关键变异信息
                top_vars = self._get_top_variants(card.variants, n=2)
                var_info = ", ".join([
                    f"{v.chrom}:{v.pos} {v.consequence[:20]}" if v.consequence else f"{v.chrom}:{v.pos}"
                    for v in top_vars
                ])
                # Build per-variant quality summary for top variants
                var_quality = []
                for v in top_vars:
                    qparts = []
                    if v.dp:
                        qparts.append(f"DP={v.dp}")
                    if v.gq is not None and v.gq < 90:
                        qparts.append(f"GQ={v.gq}")
                    if v.ad and len(v.ad) >= 2:
                        qparts.append(f"AD={v.ad[0]}/{v.ad[1]}")
                    var_quality.append(", ".join(qparts) if qparts else "—")
                quality_summary = "; ".join(var_quality) if var_quality else "—"

                # Build cross-check evidence summary for display
                cc_notes = []
                if cc.get("gnomad_af_max"):
                    cc_notes.append(f"gnomAD AF={cc['gnomad_af_max']:.1%}")
                if cc.get("clinvar_sigs"):
                    cc_notes.append(f"ClinVar={', '.join(cc['clinvar_sigs'])}")
                cc_summary = "; ".join(cc_notes) if cc_notes else "—"

                key_findings.append({
                    "gene": card.gene_symbol,
                    "subsystem": card.subsystem,
                    "level": level,
                    "rationale": card.assessment.rationale_zh[:80] if card.assessment.rationale_zh else "",
                    "variants": var_info,
                    "variant_objects": top_vars,
                    "quality_summary": quality_summary,
                    "cross_check_summary": cc_summary,
                    "inheritance_pattern": card.assessment.inheritance_pattern or "未知",
                    "phenotypic_impact": self._infer_phenotypic_impact(card.gene_symbol, card.subsystem, level, top_vars),
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

        # 预计算 LOF/GOF 分组（用于 v5 模板分层展示）
        lof_gof_variants = self._collect_lof_gof_variants(report.gene_cards)
        # 2.1 高影响非 OR 基因：纯合 LOF 或评估 >= 显著影响
        non_or_lof = [item for item in lof_gof_variants
                      if item.get("subsystem") != "olfaction"
                      and item["variant"].is_homozygous]
        # 2.4 其他子系统功能丧失变异：杂合非 OR LOF
        other_lof = [item for item in lof_gof_variants
                     if item.get("subsystem") != "olfaction"
                     and item["variant"].is_heterozygous]
        or_heterozygous_lof = [item for item in lof_gof_variants
                               if item.get("subsystem") == "olfaction"
                               and item["variant"].is_heterozygous]
        or_homozygous_lof_unknown = [ot for ot in (or_tier_b + or_tier_c)]
        non_or_impactful = [c for c in report.gene_cards
                           if c.subsystem != "olfaction"
                           and c.assessment.level in ("完全丧失", "显著影响", "部分影响")]

        # Downgraded gene symbols set
        downgraded_genes = {d["gene"] for d in downgraded_findings}

        return {
            "report": report,
            "sample_id": report.sample_id,
            "sex": report.sex,
            "ref_genome": report.ref_genome,
            "analysis_date": report.analysis_date.strftime("%Y-%m-%d %H:%M:%S"),
            "subsystems": report.subsystems,
            "gene_cards": report.gene_cards,
            # 按影响级别分层：有影响 (>=部分影响) → 主报告，其他 → 附录
            # Downgraded 显著影响/完全丧失 genes are excluded from impactful_cards
            "impactful_cards": [c for c in report.gene_cards
                if c.assessment.level in ("完全丧失", "显著影响", "部分影响")
                and c.gene_symbol not in downgraded_genes],
            "mild_cards": [c for c in report.gene_cards
                if c.assessment.level in ("可能轻微影响", "无影响")
                or c.gene_symbol in downgraded_genes],
            "subsys_cards": subsys_cards,
            "tas2r38": report.tas2r38,
            "mitochondrial": report.mitochondrial,
            "or_tiers": report.or_tiers,
            "or_tier_a": or_tier_a,
            "or_tier_b": or_tier_b,
            "or_tier_c": or_tier_c,
            "non_or_lof": non_or_lof,
            "or_heterozygous_lof": or_heterozygous_lof,
            "or_homozygous_lof_unknown": or_homozygous_lof_unknown,
            "non_or_impactful": [c for c in non_or_impactful
                if c.gene_symbol not in downgraded_genes],
            "other_lof": other_lof,
            "executive_summary": report.executive_summary,
            "personal_traits": report.executive_summary.personal_traits,
            "disclaimer_zh": report.disclaimer_zh,
            "data_availability": report.data_availability,
            "key_findings": key_findings,
            "downgraded_findings": downgraded_findings,
            "cross_checks": cross_checks,
            "profile": profile,
            "gnomad_refs": gnomad_refs,
            "key_snps": report.key_snps if report.key_snps else [],
            # SNP hits: key SNPs actually found in VCF
            "snp_hits": [s for s in (report.key_snps or []) if s.found_in_vcf],
            # High-impact cards: significant or above (excluding downgraded)
            "high_impact_cards": [c for c in report.gene_cards
                if c.assessment.level in ("完全丧失", "显著影响")
                and c.gene_symbol not in downgraded_genes],
            # LOF/GOF variants across all genes
            "lof_gof_variants": self._collect_lof_gof_variants(report.gene_cards),
            # All analyzed genes list for appendix
            "all_genes_list": sorted(report.gene_cards, key=lambda c: c.gene_symbol),
            # Phenotypic impact mapping for gene cards
            "phenotypic_impacts": {
                card.gene_symbol: self._infer_phenotypic_impact(
                    card.gene_symbol, card.subsystem, card.assessment.level,
                    self._get_top_variants(card.variants, n=2)
                )
                for card in report.gene_cards
            },
        }

    # Protein-affecting consequence types (for downgrade check eligibility)
    _PROTEIN_AFFECTING_CONS = frozenset({
        "frameshift", "stop_gained", "stop_lost", "start_lost",
        "splice_acceptor", "splice_donor", "missense",
        "inframe_insertion", "inframe_deletion", "protein_altering",
        "transcript_ablation",
    })

    @classmethod
    def _is_protein_affecting(cls, variant: Variant) -> bool:
        """Check if a variant's consequence affects protein structure/function."""
        cons = (variant.consequence or "").lower()
        for key in cls._PROTEIN_AFFECTING_CONS:
            if key in cons:
                return True
        return False

    def _cross_check_gene(self, card: GeneCard) -> Dict[str, Any]:
        """Cross-check impactful gene variants against gnomAD AF and ClinVar.

        Rules (per user request — 2026-06-07 revision):
        - ONLY check HOMOZYGOUS variants that affect protein function
          (missense, frameshift, stop_gained, splice, inframe, start_lost, etc.)
        - OR logic: ClinVar benign/likely_benign → downgrade
                     gnomAD AF > 30% → downgrade
        - All other variants (heterozygous, synonymous, intronic, UTR, etc.) → NOT checked.

        Returns dict with keys:
            downgrade (bool), downgrade_reasons (List[str]),
            gnomad_af_max (float), clinvar_sigs (List[str]),
            variant_checks (List[dict])
        """
        result: Dict[str, Any] = {
            "downgrade": False,
            "downgrade_reasons": [],
            "gnomad_af_max": 0.0,
            "clinvar_sigs": [],
            "variant_checks": [],
        }

        variants = card.key_variants or card.variants
        if not variants:
            return result

        # Step 1: Filter to ONLY homozygous + protein-affecting variants
        eligible = [
            v for v in variants
            if v.is_homozygous and self._is_protein_affecting(v)
        ]

        if not eligible:
            # No eligible variants → nothing to downgrade
            return result

        any_downgrade = False
        for v in eligible:
            # gnomAD AF (variant-level from VEP)
            af = max(
                v.gnomad_af_exome or 0.0,
                v.gnomad_af_genome or 0.0,
                v.af_gnomad or 0.0,
            )
            result["gnomad_af_max"] = max(result["gnomad_af_max"], af)

            # ClinVar (variant-level via local VCF)
            cv = query_clinvar_by_variant(v.chrom, v.pos, v.ref, v.alt)
            sig = _parse_clinvar_significance(cv["clnsig"] if cv else "")
            if sig != "unknown":
                result["clinvar_sigs"].append(sig)

            cons_short = (v.consequence or "").split(",")[0]

            var_check = {
                "variant": f"{v.chrom}:{v.pos} {v.ref}>{v.alt}",
                "consequence": cons_short,
                "af": af,
                "clinvar": cv["clnsig"] if cv else "未注释",
                "downgrade": False,
                "downgrade_reason": "",
            }

            var_downgrade = False

            # OR Rule: ClinVar benign/likely_benign → downgrade
            if sig in ("benign", "likely_benign"):
                var_downgrade = True
                var_check["downgrade_reason"] = f"ClinVar={cv['clnsig'] if cv else 'benign'}"

            # OR Rule: gnomAD AF > 30% → downgrade
            elif af > 0.30:
                var_downgrade = True
                var_check["downgrade_reason"] = f"gnomAD AF={af:.1%}"

            var_check["downgrade"] = var_downgrade
            if var_downgrade:
                any_downgrade = True

            result["variant_checks"].append(var_check)

        # Gene-level: downgrade if ANY eligible homozygous protein-affecting
        # variant triggers ClinVar benign OR gnomAD AF > 30%
        if any_downgrade:
            result["downgrade"] = True
            result["downgrade_reasons"] = [
                vc["downgrade_reason"] for vc in result["variant_checks"]
                if vc["downgrade"]
            ]

        return result

    def _collect_lof_gof_variants(self, gene_cards: List[GeneCard]) -> List[Dict[str, Any]]:
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
                        "phenotypic_impact": self._infer_phenotypic_impact(
                            card.gene_symbol, card.subsystem,
                            "完全丧失" if is_lof else "部分影响", [v]
                        ),
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
    def _infer_phenotypic_impact(gene: str, subsystem: str, level: str, variants: List[Variant]) -> str:
        """根据基因和变异推断对表型的具体影响.

        返回一段中文描述，解释该基因功能改变在实际表型上可能意味着什么。
        """
        # 获取最严重后果类型
        worst_cons = ""
        for v in variants:
            c = v.consequence or ""
            if "frameshift" in c or "stop_gained" in c or "stop_lost" in c:
                worst_cons = "lof"
                break
            elif "missense" in c:
                worst_cons = "missense"
            elif worst_cons == "":
                worst_cons = c

        # 听觉系统
        if gene == "CDH23":
            return "耳蜗毛细胞机械传感器，双等位基因失活→先天性重度感音神经性聋（DFNB12/Usher 1D）"
        if gene == "MYO7A":
            return "毛细胞静纤毛结构蛋白，双等位基因严重变异→Usher 1B（先天性聋+视网膜色素变性）"
        if gene == "GJB2":
            return "耳蜗钾离子通道蛋白，功能丧失→DFNB1先天性重度感音神经性聋"
        if gene == "OTOF":
            return "内毛细胞突触囊泡释放蛋白，双等位失活→听觉神经病变（ANSD）"
        if gene == "SLC26A4":
            return "内耳离子平衡蛋白，功能受损→大前庭导水管综合征（波动性听力下降）"
        if gene in ("MT-RNR1", "MT-TS1"):
            return "线粒体变异→氨基糖苷类抗生素耳毒性敏感，低剂量可致不可逆聋"

        # 体感/痛觉
        if gene == "SCN9A":
            return "Nav1.7痛觉通道，变异罕见，高频率多态性通常无临床意义"
        if gene == "SCN10A":
            return "Nav1.8钠通道，功能改变→炎症性/慢性疼痛易感性变化"
        if gene == "OPRM1":
            return "μ阿片受体，纯合无义变异→内源性镇痛受损，阿片类镇痛药反应减弱"
        if gene == "TRPV1":
            return "辣椒素/热敏受体，变异→辛辣灼热感耐受度和高温感知阈值变化"
        if gene == "TRPM8":
            return "冷敏感/薄荷醇受体，变异→低温敏感度和薄荷清凉感变化"
        if gene == "PIEZO2":
            return "机械敏感通道，严重受损→触觉精细度下降、本体感觉减退"

        # 视觉
        if gene in ("OPN1LW", "OPN1MW"):
            return "红/绿视蛋白（X染色体），男性单倍体，功能缺失→红绿色觉异常"
        if gene == "OPN1SW":
            return "蓝视蛋白，功能丧失→蓝黄色觉异常，红绿色觉正常"
        if gene == "RHO":
            return "视紫红质，功能严重受损→视网膜色素变性/先天性静止性夜盲"
        if gene == "ABCA4":
            return "视网膜类视黄醇转运蛋白，双等位失活→Stargardt病/锥杆细胞营养不良"
        if gene == "EYS":
            return "感光细胞外节结构蛋白，双等位失活→视网膜色素变性（RP25）"
        if gene in ("GJA8", "CRYAA"):
            return "晶状体蛋白，功能严重受损→先天性白内障"

        # 味觉
        if gene == "TAS2R38":
            return "苦味受体，PAV/PAV=敏感型，AVI/AVI=迟钝型，杂合=中间表型"
        if gene in ("TAS1R2", "TAS1R3"):
            return "甜味受体，功能严重受损→糖类甜味感知下降"
        if gene == "TAS1R1":
            return "鲜味受体，功能严重受损→谷氨酸（味精）鲜味感知下降"
        if gene in ("SCNN1A", "SCNN1B", "SCNN1G"):
            return "咸味感知钠通道，功能受损→咸味敏感度下降"
        if gene == "OTOP1":
            return "酸味感知质子通道，功能受损→酸味感知下降"

        # 色素沉着
        if gene == "OCA2":
            return "黑色素体转运蛋白，功能严重受损→眼皮肤白化病II型（OCA2）"
        if gene == "HERC2":
            return "调控OCA2表达，功能受损→浅色虹膜（蓝/绿眼）和浅色皮肤"
        if gene in ("TYRP1", "TYR"):
            return "黑色素合成关键酶，严重受损→眼皮肤白化病I/III型"
        if gene in ("SLC45A2", "SLC24A5"):
            return "黑色素体pH/离子交换蛋白，受损→皮肤/毛发色素显著变浅"

        # 代谢
        if gene == "CYP1A2":
            return "咖啡因代谢酶，AA型=快代谢耐受好，AC/CC型=慢代谢易失眠"
        if gene == "ALDH2":
            return "酒精代谢关键酶，活性降低型（GA/AA）→饮酒面部潮红（东亚常见）"
        if gene == "ADH1B":
            return "乙醇脱氢酶，快速氧化型（AG/GG）→饮酒后乙醛快速生成"
        if gene == "LCT":
            return "乳糖酶基因，CT/TT型=成年后仍可消化乳制品，CC=乳糖不耐受"

        # 肌肉
        if gene == "ACTN3":
            return "快肌α肌动蛋白，R577X纯合→爆发力下降，耐力可能不变（耐力基因型）"
        # 毛发
        if gene == "EDAR":
            return "外胚层发育蛋白，东亚370A型→直发粗发，严重受损→毛发稀疏/卷曲"
        if gene == "ABCC11":
            return "大汗腺转运蛋白，AA=干耳垢/体味轻（东亚常见），GG=湿耳垢/体味明显"

        # 嗅觉 — OR 基因
        if gene.startswith("OR"):
            return "嗅觉受体，约400个冗余受体，单个失活通常不明显影响嗅觉"
        # 嗅觉信号转导
        if gene in ("CNGA2", "ADCY3"):
            return "嗅觉cAMP信号通路核心，严重受损→可能先天性嗅觉丧失"
        # 默认描述
        if level == "完全丧失":
            return f"{gene}功能完全丧失→{subsystem or '感官'}功能可能严重受损"
        elif level == "显著影响":
            return f"{gene}功能显著受损→{subsystem or '感官'}功能可能明显异常，建议关注"
        elif level == "部分影响":
            return f"{gene}功能部分受损→{subsystem or '感官'}功能可能轻微影响"
        return "表型影响待进一步分析。"


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
