"""Markdown 报告生成器.

使用 Jinja2 模板引擎渲染人类可读的 Markdown 报告。
"""

import os
from pathlib import Path
from typing import Any, Dict, Optional

from jinja2 import Environment, FileSystemLoader, select_autoescape

from src.logger import get_logger
from src.models import SensoryReport

logger = get_logger(__name__)


class MarkdownReportGenerator:
    """Markdown 报告生成器."""

    def __init__(self, templates_dir: Optional[str] = None) -> None:
        """初始化生成器.

        Args:
            templates_dir: 模板目录路径，默认查找 ../templates/。
        """
        if templates_dir is None:
            src_dir = Path(__file__).resolve().parent.parent
            templates_dir = src_dir.parent / "templates"
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
        return {
            "report": report,
            "sample_id": report.sample_id,
            "sex": report.sex,
            "ref_genome": report.ref_genome,
            "analysis_date": report.analysis_date.strftime("%Y-%m-%d %H:%M:%S"),
            "subsystems": report.subsystems,
            "gene_cards": report.gene_cards,
            "tas2r38": report.tas2r38,
            "mitochondrial": report.mitochondrial,
            "or_tiers": report.or_tiers,
            "executive_summary": report.executive_summary,
            "disclaimer_zh": report.disclaimer_zh,
            "data_availability": report.data_availability,
        }

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
