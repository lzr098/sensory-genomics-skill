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
        """
        return report.model_dump_json(indent=2, ensure_ascii=False)

    def generate_dict(self, report: SensoryReport) -> Dict[str, Any]:
        """生成 JSON 可序列化的字典.

        Args:
            report: 报告数据模型。

        Returns:
            字典格式报告数据。
        """
        return report.model_dump(mode="json")
