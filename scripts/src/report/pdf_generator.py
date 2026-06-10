"""PDF 报告生成器.

优先通过 pandoc 将 Markdown 转换为 PDF；
当 pandoc 不可用时 fallback 到 weasyprint；
两者都不可用时生成 HTML 并提示用户手动转换.
"""

import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

from src.logger import get_logger
from src.models import SensoryReport
from src.report.html_generator import HtmlReportGenerator

logger = get_logger(__name__)


class PdfReportGenerator:
    """PDF 报告生成器.

    使用以下优先级生成 PDF:
    1. pandoc (最佳, 通过系统命令调用)
    2. weasyprint (纯 Python, 依赖复杂)
    3. 降级为 HTML 输出并提示用户
    """

    def __init__(self, html_generator: Optional[HtmlReportGenerator] = None) -> None:
        """初始化生成器.

        Args:
            html_generator: 可选的 HtmlReportGenerator 实例。
        """
        self.html_generator = html_generator or HtmlReportGenerator()

    def generate(self, report: SensoryReport) -> str:
        """生成 PDF 二进制内容（字节串）.

        Args:
            report: 报告数据模型。

        Returns:
            PDF 文件字节串；若降级为 HTML，则返回空串并记录日志。
        """
        raise NotImplementedError(
            "请使用 generate_to_file(report, output_path) 方法，"
            "因为 PDF 生成器依赖外部工具写入文件。"
        )

    def generate_to_file(self, report: SensoryReport, output_path: str) -> str:
        """生成 PDF 文件.

        按优先级尝试：pandoc → weasyprint → HTML 降级。

        Args:
            report: 报告数据模型。
            output_path: 期望的 PDF 输出路径。

        Returns:
            实际生成的文件路径。若降级为 HTML，则返回 HTML 文件路径并在日志中提示。
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # 1) 尝试 pandoc
        pandoc_path = shutil.which("pandoc")
        if pandoc_path:
            try:
                return self._generate_with_pandoc(report, output_path)
            except Exception as exc:
                logger.warning("pandoc PDF generation failed: %s", exc)

        # 2) 尝试 weasyprint
        try:
            import weasyprint  # noqa: F401
            return self._generate_with_weasyprint(report, output_path)
        except ImportError:
            logger.info("weasyprint not installed, skipping")
        except Exception as exc:
            logger.warning("weasyprint PDF generation failed: %s", exc)

        # 3) 降级：输出 HTML 并提示
        html_path = output_path.with_suffix(".html")
        self.html_generator.generate_to_file(report, str(html_path))
        logger.warning(
            "PDF 生成工具不可用（未安装 pandoc 或 weasyprint），"
            "已降级输出 HTML 报告：%s。如需 PDF，请安装 pandoc: brew install pandoc",
            html_path,
        )
        return str(html_path)

    def _generate_with_pandoc(self, report: SensoryReport, output_path: Path) -> str:
        """使用 pandoc 从 Markdown 生成 PDF."""
        from src.report.markdown_generator import MarkdownReportGenerator

        md_gen = MarkdownReportGenerator()
        md_content = md_gen.generate(report)

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".md", encoding="utf-8", delete=False
        ) as md_file:
            md_file.write(md_content)
            md_path = md_file.name

        try:
            cmd = [
                "pandoc",
                md_path,
                "-o",
                str(output_path),
                "--pdf-engine=xelatex",
                "-V",
                "CJKmainfont=PingFang SC",
                "-V",
                "geometry:margin=2.5cm",
                "--toc",
            ]
            logger.info("Running pandoc: %s", " ".join(cmd))
            subprocess.run(
                cmd,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            logger.info("PDF report saved to %s", output_path)
            return str(output_path)
        finally:
            try:
                os.unlink(md_path)
            except OSError:
                pass

    def _generate_with_weasyprint(self, report: SensoryReport, output_path: Path) -> str:
        """使用 weasyprint 从 HTML 生成 PDF."""
        import weasyprint

        html_content = self.html_generator.generate(report)
        weasyprint.HTML(string=html_content).write_pdf(str(output_path))
        logger.info("PDF report saved to %s", output_path)
        return str(output_path)
