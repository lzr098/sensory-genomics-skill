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
            "non_or_lof": non_or_lof,
            "or_heterozygous_lof": or_heterozygous_lof,
            "or_homozygous_lof_unknown": or_homozygous_lof_unknown,
            "non_or_impactful": non_or_impactful,
            "other_lof": other_lof,
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
            # Phenotypic impact mapping for gene cards
            "phenotypic_impacts": {
                card.gene_symbol: self._infer_phenotypic_impact(
                    card.gene_symbol, card.subsystem, card.assessment.level,
                    self._get_top_variants(card.variants, n=2)
                )
                for card in report.gene_cards
            },
        }

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
            return "可以把 CDH23 理解为耳蜗毛细胞上的'机械弹簧'——声音振动传来时，这个蛋白把毛细胞表面的纤毛连接在一起，形成感知机械力的'传感器'。双等位基因失活意味着这个弹簧断了，与先天性重度至极重度感音神经性听力损失（DFNB12 型）高度相关，患者通常在婴幼儿期即表现为双耳重度听力障碍。少数情况下若伴有视网膜色素变性，则符合 Usher 1D 型。"
        if gene == "MYO7A":
            return "MYO7A 参与毛细胞静纤毛的结构维护。该基因的双等位基因严重变异通常导致 Usher 1B 型（先天性聋 + 视网膜色素变性）。若检出的变异在人群中频率较高（如 >30%），则多为常见多态性，临床意义有限。"
        if gene == "GJB2":
            return "Connexin-26 是耳蜗钾离子循环的关键缝隙连接蛋白。功能丧失 → 耳蜗内电位无法维持 → 感音神经性聋（DFNB1），通常表现为语前重度-极重度听力损失。"
        if gene == "OTOF":
            return "Otoferlin 是内毛细胞突触囊泡释放的关键蛋白。双等位基因失活 → 声音信号无法从耳蜗传向听神经 → 听觉神经病变（ANSD），表现为听阈正常或轻度升高但言语识别极差。"
        if gene == "SLC26A4":
            return "Pendrin 负责内耳内淋巴液离子平衡。功能严重受损 → 内淋巴积水 → 大前庭导水管综合征（EVA），表现为波动性听力下降，头部外伤或气压变化可诱发/加重聋。"
        if gene in ("MT-RNR1", "MT-TS1"):
            return "线粒体基因变异 → 氨基糖苷类抗生素（庆大霉素、链霉素等）耳毒性敏感性显著增加 → 低剂量即可导致不可逆的感音神经性聋。"

        # 体感/痛觉
        if gene == "SCN9A":
            return "Nav1.7 是痛觉神经纤维上的'电闸'，控制疼痛信号是否向大脑传递。若检出的变异在东亚人群中频率 >80%，则属于极高频常见多态性，通常不具有临床致病意义。真正的 SCN9A 致病变异非常罕见。"
        if gene == "SCN10A":
            return "Nav1.8 钠通道主要表达于伤害性感受器。功能改变 → 炎症性疼痛和慢性疼痛的易感性可能发生变化。"
        if gene == "OPRM1":
            return "OPRM1 编码μ-阿片受体，是身体'天然止痛系统'的核心开关，也是吗啡、芬太尼等止痛药的作用靶点。纯合无义变异导致受体蛋白提前终止，产生截短的无功能蛋白。可能的影响：内源性镇痛系统受损，对阿片类镇痛药（术后吗啡、癌痛芬太尼）的反应可能显著减弱，常规剂量可能效果不佳。"
        if gene == "TRPV1":
            return "TRPV1 是辣椒素受体兼热敏感受器。功能变异 → 对辣椒等辛辣食物的灼热感耐受度、以及对高温伤害的感知阈值可能发生变化。"
        if gene == "TRPM8":
            return "TRPM8 是冷敏感受器和薄荷醇受体。功能变异 → 对低温环境的冷感敏感度、以及对薄荷清凉感的感知可能发生变化。"
        if gene == "PIEZO2":
            return "PIEZO2 是机械敏感离子通道，负责触觉分辨和本体感觉。功能严重受损 → 触觉精细度下降、本体感觉减退（闭眼时难以判断肢体位置）、共济失调步态。"

        # 视觉
        if gene in ("OPN1LW", "OPN1MW"):
            return "红/绿视蛋白基因位于 X 染色体。男性为单倍体，任何功能缺失变异都可能导致红绿色觉异常（红色弱/绿色弱或全色盲）。本样本为男性，需关注。"
        if gene == "OPN1SW":
            return "蓝视蛋白基因。功能丧失 → 蓝黄色觉异常（tritanopia），表现为蓝/黄颜色辨别困难，但红/绿色觉正常。"
        if gene == "RHO":
            return "视紫红质是杆状细胞光信号传导的核心蛋白。功能严重受损 → 视网膜色素变性（RP）或先天性静止性夜盲（CSNB），表现为夜盲、进行性视野缩小。"
        if gene == "ABCA4":
            return "ABCA4 负责视网膜色素上皮细胞中类视黄醇代谢产物的转运。双等位基因严重功能受损 → Stargardt 病（青少年黄斑变性）或锥杆细胞营养不良，表现为中心视力进行性下降。"
        if gene == "EYS":
            return "EYS 是视网膜感光细胞外节结构蛋白。双等位基因功能丧失 → 视网膜色素变性（RP25），表现为夜盲和进行性周边视野丧失。"
        if gene in ("GJA8", "CRYAA"):
            return "晶状体结构/代谢蛋白。功能严重受损 → 先天性白内障，表现为婴幼儿期即出现的晶状体混浊和视力障碍。"

        # 味觉
        if gene == "TAS2R38":
            return "苦味受体 TAS2R38 的基因型决定了对苯硫脲（PTC）等苦味物质的敏感度。PAV/PAV（味觉敏感型）→ 对苦味高度敏感；AVI/AVI（味觉迟钝型）→ 对苦味不敏感；杂合型 → 中间表型。"
        if gene in ("TAS1R2", "TAS1R3"):
            return "TAS1R2/TAS1R3 构成甜味受体。功能严重受损 → 对糖类甜味物质的感知能力下降，可能表现为'食不知甜'。"
        if gene == "TAS1R1":
            return "TAS1R1 与 TAS1R3 构成鲜味受体。功能严重受损 → 对谷氨酸（味精）等鲜味物质的感知能力下降。"
        if gene in ("SCNN1A", "SCNN1B", "SCNN1G"):
            return "上皮钠通道（ENaC）亚基，负责咸味感知。功能严重受损 → 对咸味的敏感度下降，可能偏好更咸的食物。"
        if gene == "OTOP1":
            return "OTOP1 是质子通道，负责酸味感知。功能严重受损 → 对酸性物质的酸味感知能力下降。"

        # 色素沉着
        if gene == "OCA2":
            return "OCA2 是黑色素体膜转运蛋白，影响黑色素合成量。功能严重受损 → 眼皮肤白化病 II 型（OCA2），表现为极浅色皮肤、浅色头发、蓝/灰/淡褐色虹膜，伴眼球震颤和视力低下。"
        if gene == "HERC2":
            return "HERC2 调控 OCA2 基因表达。功能严重受损 → 通过下调 OCA2 间接导致色素减少，表现为浅色虹膜（蓝/绿色眼）和浅色皮肤。"
        if gene in ("TYRP1", "TYR"):
            return "酪氨酸酶相关蛋白，黑色素合成通路的关键酶。功能严重受损 → 眼皮肤白化病 I/III 型，表现为极浅色皮肤和毛发，伴严重视力障碍。"
        if gene in ("SLC45A2", "SLC24A5"):
            return "黑色素体 pH 调节蛋白 / 钾离子交换蛋白，影响黑色素合成效率。功能严重受损 → 皮肤/毛发色素显著变浅，表现为白化病或极浅色表型。"

        # 代谢
        if gene == "CYP1A2":
            return "CYP1A2 是咖啡因代谢的主要酶。快代谢型（AA）→ 咖啡因清除快，耐受性好，不易因咖啡因导致心悸或失眠；慢代谢型（AC/CC）→ 咖啡因半衰期长，少量摄入即可能引起焦虑、失眠。"
        if gene == "ALDH2":
            return "乙醛脱氢酶2是酒精代谢的关键酶。正常活性型（GG）→ 饮酒后乙醛可迅速代谢为乙酸，不易出现面部潮红和不适；活性降低型（GA/AA）→ 乙醛蓄积 → 饮酒后快速面部潮红、心悸、恶心（东亚人群中常见）。"
        if gene == "ADH1B":
            return "乙醇脱氢酶2催化乙醇氧化为乙醛。快速氧化型（AG/GG）→ 饮酒后乙醛迅速生成，可能更快出现醉酒感；慢速氧化型（AA）→ 乙醇代谢较慢。"
        if gene == "LCT":
            return "LCT 编码乳糖酶。乳糖酶持久型（CT/TT）→ 成年后乳糖酶持续表达，可正常消化乳制品；非持久型（CC）→ 成年后乳糖酶活性下降 → 摄入乳制品后可能出现腹胀、腹泻、肠鸣（乳糖不耐受）。"

        # 肌肉
        if gene == "ACTN3":
            return "ACTN3（α-肌动蛋白-3）主要存在于快肌纤维 IIx 型。R577X 无义变异纯合（XX）→ 快肌纤维功能下降 → 爆发力/短跑/力量型运动能力降低，但耐力运动表现可能不受影响甚至略有优势（'耐力基因型'）。杂合（RX）→ 中间表型。"

        # 毛发
        if gene == "EDAR":
            return "EDAR 参与外胚层衍生物（毛发、牙齿、汗腺、乳腺）的发育。东亚典型变异（370A）→ 直发、粗发、汗腺发达、乳腺导管密度增加；功能严重受损 → 可能表现为毛发稀疏/卷曲、汗腺发育不全、牙齿发育异常。"
        if gene == "ABCC11":
            return "ABCC11 参与大汗腺分泌物的转运。功能丧失型（AA）→ 干耳垢（片状）、腋下体味较轻（大汗腺分泌减少），在东亚人群中极为常见；功能正常型（GG）→ 湿耳垢、体味较明显。"

        # 嗅觉 — OR 基因
        if gene.startswith("OR"):
            return "嗅觉受体（OR）基因负责识别特定气味分子。由于人类嗅觉系统存在大量冗余受体（约400个），单个 OR 基因失活通常不会导致完全无法感知某种气味，仅可能降低对该气味分子的敏感度或辨识度。"

        # 嗅觉信号转导
        if gene in ("CNGA2", "ADCY3"):
            return "CNGA2/ADCY3 是嗅觉信号转导通路的核心组件（cAMP 信号级联）。功能严重受损 → 可能导致先天性嗅觉丧失（anosmia）或嗅觉显著减退，因为所有 OR 信号都依赖此通路传递。"

        # 默认描述
        if level == "完全丧失":
            return f"{gene} 功能完全丧失 → 该基因所参与的{subsystem or '感官'}功能可能严重受损，具体表型取决于该基因在通路中的位置和冗余性。"
        elif level == "显著影响":
            return f"{gene} 功能显著受损 → 该基因所参与的{subsystem or '感官'}功能可能出现明显异常，建议在相关临床表型上加以关注。"
        elif level == "部分影响":
            return f"{gene} 功能部分受损 → 该基因所参与的{subsystem or '感官'}功能可能有轻微影响，但由于系统冗余或其他代偿机制，实际表型变化可能不明显。"
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
