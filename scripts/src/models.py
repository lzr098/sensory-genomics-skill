"""Pydantic 数据模型定义.

定义跨模块传递的结构化数据契约，包括变异、评估、报告等核心实体。
"""

from __future__ import annotations

from datetime import datetime
from typing import Any, Dict, List, Literal, Optional

from pydantic import BaseModel, Field

ImpactLevel = Literal[
    "完全丧失",
    "显著影响",
    "部分影响",
    "可能轻微影响",
    "无影响",
]

InheritancePattern = Literal[
    "隐性纯合",
    "隐性复合杂合",
    "X-连锁",
    "显性",
    "线粒体",
    "未知",
]

Sex = Literal["M", "F"]

Subsystem = Literal[
    "vision",
    "hearing",
    "olfaction",
    "taste",
    "somatosensation",
    "pigmentation",
    "metabolism",
    "muscle",
    "hair",
]


class Variant(BaseModel):
    """VCF 变异记录 + VEP 注释后的完整模型."""

    chrom: str = Field(..., description="染色体")
    pos: int = Field(..., description="位置（1-based）")
    ref: str = Field(..., description="参考等位基因")
    alt: str = Field(..., description="替代等位基因")
    gt: str = Field(..., description="基因型，如 0/1, 1/1, 0|1")
    dp: int = Field(0, description="测序深度")
    gq: Optional[int] = Field(None, description="基因型质量 (Phred-scaled)")
    ad: Optional[List[int]] = Field(None, description="等位基因深度 [ref, alt]")
    qual: float = Field(0.0, description="质量值")
    filter_status: str = Field("PASS", description="FILTER 字段")

    # VEP 注释字段
    gene_symbol: str = Field("", description="主基因符号（最优转录本）")
    gene_symbols: List[str] = Field(default_factory=list, description="该变异涉及的所有感官基因符号")
    consequence: str = Field("", description="VEP 后果类型，如 missense_variant")
    hgvsc: Optional[str] = Field(None, description="HGVSc")
    hgvsp: Optional[str] = Field(None, description="HGVSp")
    protein_domain: Optional[str] = Field(None, description="蛋白域")
    protein_topology: Optional[str] = Field(None, description="蛋白拓扑位置")
    af_gnomad: Optional[float] = Field(None, description="gnomAD 等位基因频率")
    raw_vep: Dict[str, Any] = Field(default_factory=dict, description="原始 VEP 响应")

    # 扩展注释字段
    rsid: Optional[str] = Field(None, description="dbSNP rsID")
    sift: Optional[str] = Field(None, description="SIFT 预测")
    polyphen: Optional[str] = Field(None, description="PolyPhen 预测")
    gnomad_af_exome: Optional[float] = Field(None, description="gnomAD exome AF")
    gnomad_af_genome: Optional[float] = Field(None, description="gnomAD genome AF")
    lof_flags: Optional[str] = Field(None, description="LoF flags (lc_lof, etc)")
    cadd_score: Optional[float] = Field(None, description="CADD 分数")
    amino_acid_change: Optional[str] = Field(None, description="氨基酸替换，如 D1919G")
    protein_position: Optional[int] = Field(None, description="蛋白序列位置")

    # 关键变异标记（由评估引擎填充）
    is_key_variant: bool = Field(False, description="是否为关键功能变异")
    key_variant_reason: str = Field("", description="关键变异的判定理由")

    def __hash__(self) -> int:
        return hash((self.chrom, self.pos, self.ref, self.alt, self.gt))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Variant):
            return NotImplemented
        return (
            self.chrom == other.chrom
            and self.pos == other.pos
            and self.ref == other.ref
            and self.alt == other.alt
            and self.gt == other.gt
        )

    @property
    def vcf_id(self) -> str:
        """返回变异的唯一标识字符串."""
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"

    @property
    def is_homozygous(self) -> bool:
        """判断是否为纯合变异."""
        return "1/1" in self.gt or "1|1" in self.gt

    @property
    def is_heterozygous(self) -> bool:
        """判断是否为杂合变异."""
        return "0/1" in self.gt or "0|1" in self.gt or "1/0" in self.gt or "1|0" in self.gt

    @property
    def is_hemizygous(self) -> bool:
        """判断是否为半合子（男性 X 染色体）."""
        return "1" in self.gt and "/" not in self.gt and "|" not in self.gt


class ImpactAssessment(BaseModel):
    """功能影响评估结果."""

    level: ImpactLevel = Field("无影响", description="影响程度五级")
    protein_impact: str = Field("", description="蛋白影响评估文本")
    gene_certainty: str = Field("", description="基因功能确定性评估")
    zygosity_match: bool = Field(False, description="基因型是否符合遗传模式")
    inheritance_pattern: InheritancePattern = Field("未知", description="遗传模式")
    rationale_zh: str = Field("", description="影响依据（中文）")
    limitation_note: Optional[str] = Field(None, description="局限说明")


class GeneCard(BaseModel):
    """基因卡片：每个感官基因的完整分析结果."""

    gene_symbol: str = Field(..., description="基因符号")
    subsystem: str = Field("", description="所属感官子系统")
    sensory_function_zh: str = Field("", description="该基因在感官系统中的功能（中文）")
    variants: List[Variant] = Field(default_factory=list, description="该基因上的全部变异列表")
    key_variants: List[Variant] = Field(default_factory=list, description="该基因上的关键功能变异列表（基因维度展示用）")
    assessment: ImpactAssessment = Field(default_factory=ImpactAssessment, description="功能影响评估")
    protein_impact_summary: str = Field("", description="综合蛋白影响摘要（基于关键变异）")
    special_logic_type: Optional[str] = Field(None, description="专用逻辑类型标记（如 tas2r38, mt, or_tier）")
    enrichment_data: Dict[str, Any] = Field(default_factory=dict, description="API 富集数据")


class TAS2R38Result(BaseModel):
    """TAS2R38 Haplotype 分析结果."""

    rs713598_gt: str = Field("", description="rs713598 基因型")
    rs1726866_gt: str = Field("", description="rs1726866 基因型")
    rs10246939_gt: str = Field("", description="rs10246939 基因型")
    diplotype: str = Field("", description="双体型，如 PAV/AVI")
    phenotype_zh: str = Field("", description="苦味感知表型（中文）")
    phenotype_level: str = Field("", description="苦味感知能力等级")


class MitochondrialResult(BaseModel):
    """线粒体耳聋注释结果."""

    variant_name: str = Field(..., description="变异名称，如 m.1555A>G")
    gene: str = Field(..., description="线粒体基因")
    heteroplasmy: float = Field(0.0, description="异质性水平（0-1）")
    drug_warning_zh: str = Field("", description="药物风险警告（中文）")
    risk_level: str = Field("", description="风险等级")


class ORTierResult(BaseModel):
    """OR 基因分级展示结果."""

    tier: Literal["A", "B", "C"] = Field(..., description="分级")
    gene_symbol: str = Field(..., description="OR 基因符号")
    known_ligand_zh: Optional[str] = Field(None, description="已知配体（中文）")
    odor_description_zh: Optional[str] = Field(None, description="气味描述（中文）")
    variant: Variant = Field(..., description="相关变异")
    assessment: ImpactAssessment = Field(default_factory=ImpactAssessment, description="功能影响评估")


class DataAvailability(BaseModel):
    """数据可用性状态."""

    gnomad_af: str = Field("N/A", description="gnomAD AF 状态")
    clinvar: str = Field("N/A", description="ClinVar 状态")
    spliceai: str = Field("N/A", description="SpliceAI 状态")
    cadd: str = Field("N/A", description="CADD 状态")
    topology: str = Field("N/A", description="蛋白拓扑状态")


class KeySNPResult(BaseModel):
    """关键性状 SNP 推断结果."""

    rsid: str = Field(..., description="dbSNP rsID")
    gene: str = Field(..., description="基因符号")
    chrom: str = Field(..., description="染色体")
    pos: int = Field(..., description="GRCh38 位置")
    ref: str = Field(..., description="参考等位基因")
    alt: List[str] = Field(default_factory=list, description="替代等位基因列表")
    inferred_genotype: str = Field(..., description="推断基因型（如 AA, AG, GG）")
    is_heterozygous: bool = Field(False, description="是否杂合")
    is_homozygous_alt: bool = Field(False, description="是否纯合变异")
    is_homozygous_ref: bool = Field(False, description="是否纯合参考")
    phenotype_label: str = Field(..., description="表型标签")
    phenotype_description: str = Field("", description="表型描述")
    notes: str = Field("", description="注释说明")
    found_in_vcf: bool = Field(False, description="是否在 VCF 中检出")
    dp: Optional[int] = Field(None, description="测序深度")
    gq: Optional[int] = Field(None, description="基因型质量 (Phred-scaled)")
    ad_ref: Optional[int] = Field(None, description="参考等位基因深度")
    ad_alt: Optional[int] = Field(None, description="替代等位基因深度")


class PersonalTraitPrediction(BaseModel):
    """个人特征定性预测结果."""

    trait: str = Field(..., description="特征名称，如 瞳孔颜色")
    subsystem: str = Field("", description="所属子系统")
    prediction: str = Field("", description="定性预测结论")
    confidence: str = Field("", description="置信度：高/中/低")
    evidence: str = Field("", description="推断依据简述")
    key_genes: List[str] = Field(default_factory=list, description="关键基因列表")
    key_snps: List[str] = Field(default_factory=list, description="关键 SNP 列表")


class ExecutiveSummary(BaseModel):
    """执行摘要统计."""

    subsystem_counts: Dict[str, Dict[str, int]] = Field(
        default_factory=dict,
        description="各子系统 × 影响程度计数",
    )
    key_findings: List[str] = Field(default_factory=list, description="值得关注的关键发现")
    personal_traits: List[PersonalTraitPrediction] = Field(
        default_factory=list, description="个人特征定性预测"
    )


class SensoryReport(BaseModel):
    """完整报告数据模型."""

    sample_id: str = Field("", description="样本 ID")
    sex: Sex = Field("M", description="样本性别")
    ref_genome: str = Field("GRCh38", description="参考基因组")
    analysis_date: datetime = Field(default_factory=datetime.now, description="分析日期")
    subsystems: List[Subsystem] = Field(default_factory=list, description="分析子系统列表")
    gene_cards: List[GeneCard] = Field(default_factory=list, description="全部基因卡片")
    tas2r38: Optional[TAS2R38Result] = Field(None, description="TAS2R38 分析结果")
    mitochondrial: Optional[List[MitochondrialResult]] = Field(None, description="线粒体分析结果")
    or_tiers: Optional[List[ORTierResult]] = Field(None, description="OR 分级结果")
    key_snps: Optional[List[KeySNPResult]] = Field(None, description="关键性状 SNP 推断结果")
    executive_summary: ExecutiveSummary = Field(default_factory=ExecutiveSummary, description="执行摘要")
    data_availability: Dict[str, DataAvailability] = Field(
        default_factory=dict, description="各基因数据可用性状态"
    )
    enrichment_summary: Dict[str, Any] = Field(
        default_factory=dict, description="Stage 6 API 富集策略与统计"
    )
    disclaimer_zh: str = Field("", description="免责声明（中文）")


class AnalysisConfig(BaseModel):
    """分析运行时配置."""

    vcf_path: str = Field(..., description="VCF 文件路径")
    sex: Sex = Field(..., description="样本性别 M/F")
    subsystems: List[Subsystem] = Field(
        default_factory=lambda: ["vision", "hearing", "olfaction", "taste", "somatosensation"],
        description="分析子系统子集",
    )
    known_phenotype: Optional[str] = Field(None, description="已知表型文本")
    show_reference_info: bool = Field(True, description="是否展示参考信息")
    vep_source: str = Field("rest_api", description="VEP 来源")
    output_dir: Optional[str] = Field(None, description="输出目录")
    strict_filter: bool = Field(False, description="使用精确外显子/CDS坐标预过滤（排除内含子变异）")
    output_format: str = Field("markdown", description="报告输出格式：markdown/html/pdf")
    precompute_db: Optional[str] = Field(None, description="预计算 VEP SQLite 数据库路径，运行时优先查库")
    original_vcf_path: Optional[str] = Field(None, description="原始 VCF 路径（BED 过滤前），用于 KeySNPInferrer 查询")
    auto_sex: bool = Field(False, description="从 VCF chrX/chrY 基因型模式自动推断性别")
    enrich_all_genes: bool = Field(False, description="对所有候选基因执行 API 富集（默认仅富集高影响基因和关键 OR 基因）")
