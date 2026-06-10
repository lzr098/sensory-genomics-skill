"""HTML 报告生成器.

复用 Markdown 报告的内容与上下文，通过 markdown 库转换为 HTML，
并包装为带内联 CSS 的打印友好页面。
"""

import html
import os
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional

from src.logger import get_logger
from src.models import SensoryReport
from src.report.markdown_generator import MarkdownReportGenerator

logger = get_logger(__name__)


def _css() -> str:
    """返回报告用内联 CSS."""
    return """
    :root {
        --color-bg: #ffffff;
        --color-text: #1f2937;
        --color-muted: #6b7280;
        --color-primary: #2563eb;
        --color-danger: #dc2626;
        --color-warning: #d97706;
        --color-success: #059669;
        --color-info: #0891b2;
        --color-border: #e5e7eb;
        --font-sans: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
            "Helvetica Neue", Arial, "Noto Sans SC", "PingFang SC",
            "Microsoft YaHei", sans-serif;
        --font-mono: "SFMono-Regular", Consolas, "Liberation Mono", Menlo,
            "Courier New", monospace;
    }
    * { box-sizing: border-box; }
    body {
        margin: 0; padding: 0;
        font-family: var(--font-sans);
        font-size: 14px; line-height: 1.7;
        color: var(--color-text);
        background: var(--color-bg);
    }
    .report-container {
        max-width: 920px; margin: 0 auto; padding: 40px 32px;
    }
    .report-header {
        border-bottom: 3px solid var(--color-primary);
        padding-bottom: 24px; margin-bottom: 32px;
    }
    .report-header h1 {
        margin: 0 0 8px 0; font-size: 28px; font-weight: 700;
        color: var(--color-primary);
    }
    .report-meta {
        display: flex; flex-wrap: wrap; gap: 16px;
        color: var(--color-muted); font-size: 13px;
    }
    .disclaimer-box {
        background: #fffbeb; border-left: 4px solid var(--color-warning);
        padding: 14px 18px; margin: 24px 0; border-radius: 6px;
        font-size: 13px; color: #92400e;
    }
    h2 {
        margin-top: 36px; margin-bottom: 14px;
        font-size: 20px; font-weight: 600;
        border-bottom: 1px solid var(--color-border);
        padding-bottom: 8px;
    }
    h3 {
        margin-top: 26px; margin-bottom: 10px;
        font-size: 16px; font-weight: 600;
        color: var(--color-text);
    }
    h4 { margin-top: 18px; margin-bottom: 8px; font-size: 14px; font-weight: 600; }
    p { margin: 0 0 12px 0; }
    table {
        width: 100%; border-collapse: collapse; margin: 14px 0;
        font-size: 13px;
    }
    th, td {
        border: 1px solid var(--color-border);
        padding: 10px 12px; text-align: left; vertical-align: top;
    }
    th { background: #f9fafb; font-weight: 600; }
    tr:nth-child(even) td { background: #fafafa; }
    code {
        font-family: var(--font-mono); font-size: 12px;
        background: #f3f4f6; padding: 2px 5px; border-radius: 4px;
    }
    blockquote {
        margin: 14px 0; padding: 10px 16px;
        border-left: 4px solid var(--color-info);
        background: #f0f9ff; color: #0c4a6e;
    }
    ul, ol { margin: 10px 0; padding-left: 24px; }
    li { margin-bottom: 4px; }
    details {
        border: 1px solid var(--color-border); border-radius: 6px;
        padding: 12px 16px; margin: 12px 0; background: #fafafa;
    }
    summary {
        font-weight: 600; cursor: pointer; color: var(--color-primary);
    }
    .badge {
        display: inline-block; padding: 2px 8px; border-radius: 999px;
        font-size: 12px; font-weight: 600; white-space: nowrap;
    }
    .badge-danger { background: #fee2e2; color: #991b1b; }
    .badge-warning { background: #fef3c7; color: #92400e; }
    .badge-success { background: #d1fae5; color: #065f46; }
    .badge-info { background: #e0f2fe; color: #0c4a6e; }
    .status-bar { font-size: 15px; letter-spacing: 1px; }
    .gene-card {
        border: 1px solid var(--color-border); border-radius: 8px;
        padding: 18px; margin: 16px 0; background: #ffffff;
        box-shadow: 0 1px 2px rgba(0,0,0,0.03);
    }
    .gene-card-header {
        display: flex; justify-content: space-between; align-items: flex-start;
        flex-wrap: wrap; gap: 8px; margin-bottom: 10px;
    }
    .gene-card-header strong { font-size: 16px; color: var(--color-primary); }
    .gene-card-section {
        margin-top: 12px; padding-top: 12px;
        border-top: 1px dashed var(--color-border);
    }
    .gene-card-section-title {
        font-size: 12px; font-weight: 600; color: var(--color-muted);
        text-transform: uppercase; letter-spacing: 0.5px;
        margin-bottom: 6px;
    }
    .section-empty {
        color: var(--color-muted); font-style: italic; padding: 8px 0;
    }
    .footer {
        margin-top: 48px; padding-top: 16px;
        border-top: 1px solid var(--color-border);
        font-size: 12px; color: var(--color-muted); text-align: center;
    }
    @media print {
        body { font-size: 12px; }
        .report-container { padding: 24px; max-width: 100%; }
        h2 { page-break-after: avoid; }
        .gene-card { page-break-inside: avoid; }
        table { page-break-inside: avoid; }
        details { page-break-inside: avoid; }
    }
    """


def _html_wrapper(title: str, meta_html: str, body_html: str, generated_at: str) -> str:
    """将 body HTML 包装为完整 HTML 文档."""
    css = _css()
    return f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{html.escape(title)}</title>
    <style>{css}</style>
</head>
<body>
<div class="report-container">
    <header class="report-header">
        <h1>{html.escape(title)}</h1>
        <div class="report-meta">{meta_html}</div>
    </header>
    {body_html}
    <footer class="footer">
        报告生成时间：{generated_at} · Sensory Genomics Skill
    </footer>
</div>
</body>
</html>
"""


class HtmlReportGenerator:
    """HTML 报告生成器.

    复用 Markdown 报告生成器的内容，转换为带样式的 HTML 输出。
    """

    def __init__(self, markdown_generator: Optional[MarkdownReportGenerator] = None) -> None:
        """初始化生成器.

        Args:
            markdown_generator: 可选的 MarkdownReportGenerator 实例用于复用。
        """
        self.md_generator = markdown_generator or MarkdownReportGenerator()

    def generate(self, report: SensoryReport) -> str:
        """生成 HTML 报告字符串.

        Args:
            report: 报告数据模型。

        Returns:
            HTML 格式报告字符串。
        """
        try:
            import markdown as md_lib
        except ImportError as exc:
            logger.error("markdown library not installed: %s", exc)
            raise RuntimeError(
                "HTML 报告需要 'markdown' 库。请运行: "
                "pip install markdown"
            ) from exc

        # 1) 先生成 Markdown
        md_content = self.md_generator.generate(report)

        # 2) 转换为 HTML（支持表格、extra 扩展）
        extensions = [
            "tables",
            "fenced_code",
            "toc",
            "nl2br",
            "sane_lists",
        ]
        body_html = md_lib.markdown(md_content, extensions=extensions)

        # 3) 构建页头信息
        sample_id = report.sample_id or "Unknown"
        title = f"感官基因组学分析报告 — {html.escape(sample_id)}"
        meta_parts = [
            f"<span>样本 ID：<strong>{html.escape(sample_id)}</strong></span>",
            f"<span>性别：<strong>{'男性' if report.sex == 'M' else '女性' if report.sex == 'F' else report.sex}</strong></span>",
            f"<span>参考基因组：<strong>{html.escape(report.ref_genome or 'GRCh38')}</strong></span>",
        ]
        if report.executive_summary and report.executive_summary.personal_traits:
            traits_count = len(report.executive_summary.personal_traits)
            meta_parts.append(f"<span>推断特征：<strong>{traits_count} 项</strong></span>")

        meta_html = "\n".join(meta_parts)
        generated_at = datetime.now().strftime("%Y-%m-%d %H:%M")

        return _html_wrapper(title, meta_html, body_html, generated_at)

    def generate_to_file(self, report: SensoryReport, output_path: str) -> str:
        """生成并保存 HTML 报告到文件.

        Args:
            report: 报告数据模型。
            output_path: 输出 HTML 文件路径。

        Returns:
            输出文件路径。
        """
        html_content = self.generate(report)
        Path(output_path).write_text(html_content, encoding="utf-8")
        logger.info("HTML report saved to %s", output_path)
        return output_path
