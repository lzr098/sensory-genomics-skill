"""统一日志配置.

提供统一格式的日志输出到文件和控制台，按运行时间命名日志文件。
"""

import logging
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional


class SkillLogFormatter(logging.Formatter):
    """自定义日志格式器，包含毫秒时间戳."""

    def __init__(self, fmt: str = "%(asctime)s [%(levelname)s] %(name)s: %(message)s") -> None:
        super().__init__(fmt, datefmt="%Y-%m-%d %H:%M:%S")

    def formatTime(self, record: logging.LogRecord, datefmt: Optional[str] = None) -> str:
        ct = datetime.fromtimestamp(record.created)
        if datefmt:
            s = ct.strftime(datefmt)
        else:
            s = ct.strftime("%Y-%m-%d %H:%M:%S")
        return f"{s}.{int(record.msecs):03d}"


def setup_logger(
    name: str = "sensory_genomics",
    level: int = logging.INFO,
    log_dir: Optional[str] = None,
    file_prefix: str = "sensory_genomics",
    console: bool = True,
) -> logging.Logger:
    """配置并返回统一格式的 Logger.

    若该 name 的 Logger 已配置过 Handler，直接返回已有实例，
    避免重复添加 Handler 导致日志文件分裂和重复输出。

    Args:
        name: Logger 名称.
        level: 日志级别.
        log_dir: 日志文件保存目录，默认 ~/.workbuddy/skills/sensory-genomics/logs/.
        file_prefix: 日志文件名前缀.
        console: 是否同时输出到控制台.

    Returns:
        配置好的 Logger 实例.
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)

    # 避免重复添加 Handler（同一进程内多次初始化时）
    if logger.handlers:
        return logger

    formatter = SkillLogFormatter()

    # 文件 Handler
    if log_dir is None:
        log_dir = os.path.expanduser("~/.workbuddy/skills/sensory-genomics/logs/")
    os.makedirs(log_dir, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = os.path.join(log_dir, f"{file_prefix}_{timestamp}.log")
    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setLevel(level)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # 控制台 Handler
    if console:
        stream_handler = logging.StreamHandler(sys.stdout)
        stream_handler.setLevel(level)
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)

    logger.info("Logger initialized. Log file: %s", log_path)
    return logger


def get_logger(name: str = "sensory_genomics") -> logging.Logger:
    """获取已配置的 Logger，若未配置则返回默认 Logger."""
    logger = logging.getLogger(name)
    if not logger.handlers:
        return setup_logger(name)
    return logger
