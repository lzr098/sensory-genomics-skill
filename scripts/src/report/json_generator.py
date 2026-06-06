"""JSON 结构化输出生成器.

将 SensoryReport 序列化为机器可读的 JSON 格式。
"""

import json
from typing import Any, Dict

from src.logger import get_logger
from src.models import SensoryReport

logger = get_logger(__name__)


class JsonReportGenerator:
    """JSON 报告生成器."""

    def generate(self, report: SensoryReport) -> str:
        """生成 JSON 报告.

        Args:
            report: 报告数据模型。

        Returns:
            JSON 格式字符串（缩进 2 空格）。
            默认排除 raw_vep 以控制文件体积。
        """
        return report.model_dump_json(
            indent=2,
            ensure_ascii=False,
            exclude={"gene_cards": {"__all__": {"variants": {"__all__": {"raw_vep"}}}}},
        )

    def generate_dict(self, report: SensoryReport, *, include_raw_vep: bool = False) -> Dict[str, Any]:
        """生成 JSON 可序列化的字典.

        Args:
            report: 报告数据模型。
            include_raw_vep: 是否包含原始 VEP 数据（默认关闭以控制体积）。

        Returns:
            字典格式报告数据。
        """
        if include_raw_vep:
            return report.model_dump(mode="json")
        return report.model_dump(
            mode="json",
            exclude={"gene_cards": {"__all__": {"variants": {"__all__": {"raw_vep"}}}}},
        )
