"""YAML 配置加载与校验.

将 config.yaml 加载为 Pydantic 模型，支持路径展开（~ → HOME）。
"""

import os
from pathlib import Path
from typing import Any, Dict, Optional

import yaml
from pydantic import BaseModel, Field, field_validator

from src.exceptions import ConfigError
from src.logger import get_logger

logger = get_logger(__name__)


class VepConfig(BaseModel):
    """VEP 配置子模型."""

    source: str = Field("auto", description="VEP 来源: auto | rest_api | local_docker")
    base_url: str = Field("https://rest.ensembl.org")
    batch_size: int = Field(200)
    rate_limit: int = Field(15)
    max_retries: int = Field(3)
    timeout: int = Field(60)
    cache_dir: str = Field("~/.workbuddy/tools/vep/cache", description="本地 VEP cache 目录路径")


class FilterConfig(BaseModel):
    """过滤配置子模型."""

    min_qual: int = Field(30)
    min_dp: int = Field(10)
    pass_only: bool = Field(True)


class CacheConfig(BaseModel):
    """缓存配置子模型."""

    db_path: str = Field("~/.workbuddy/skills/sensory-genomics/cache.sqlite")
    default_ttl_days: int = Field(30)


class RateLimitConfig(BaseModel):
    """限速配置子模型."""

    ensembl: int = Field(15)
    ncbi: int = Field(3)
    uniprot: int = Field(10)


class OutputConfig(BaseModel):
    """输出配置子模型."""

    default_dir: str = Field("~/.workbuddy/skills/sensory-genomics/output/")
    formats: list = Field(default_factory=lambda: ["markdown", "json"])


class ReferenceInfoConfig(BaseModel):
    """参考信息展示配置."""

    show_gnomad_af: bool = Field(True)
    show_clinvar: bool = Field(True)
    show_spliceai: bool = Field(True)
    show_cadd: bool = Field(True)
    show_topology: bool = Field(True)


class LoggingConfig(BaseModel):
    """日志配置子模型."""

    level: str = Field("INFO")
    file_prefix: str = Field("sensory_genomics")
    format: str = Field("%(asctime)s [%(levelname)s] %(name)s: %(message)s")


class SkillConfig(BaseModel):
    """Skill 顶层配置模型."""

    skill: Dict[str, Any] = Field(default_factory=dict)
    vep: VepConfig = Field(default_factory=VepConfig)
    filter: FilterConfig = Field(default_factory=FilterConfig)
    cache: CacheConfig = Field(default_factory=CacheConfig)
    rate_limits: RateLimitConfig = Field(default_factory=RateLimitConfig)
    output: OutputConfig = Field(default_factory=OutputConfig)
    subsystems: list = Field(
        default_factory=lambda: ["vision", "hearing", "olfaction", "taste", "somatosensation"]
    )
    reference_info: ReferenceInfoConfig = Field(default_factory=ReferenceInfoConfig)
    logging: LoggingConfig = Field(default_factory=LoggingConfig)
    precompute_db: Optional[str] = Field(None, description="预计算 VEP SQLite 数据库路径")

    @field_validator("subsystems", mode="before")
    @classmethod
    def validate_subsystems(cls, v: Any) -> Any:
        valid = {"vision", "hearing", "olfaction", "taste", "somatosensation", "pigmentation", "metabolism", "muscle", "hair"}
        if isinstance(v, list):
            invalid = set(v) - valid
            if invalid:
                raise ConfigError(f"Invalid subsystems: {invalid}")
        return v


def load_config(config_path: Optional[str] = None) -> SkillConfig:
    """加载 YAML 配置文件并校验.

    Args:
        config_path: 配置文件路径，默认查找当前目录 config.yaml。

    Returns:
        校验后的 SkillConfig 实例。
    """
    if config_path is None:
        config_path = "config.yaml"

    path = Path(config_path)
    if not path.exists():
        # 尝试从工作目录查找
        work_dir = Path(__file__).resolve().parent.parent
        path = work_dir / "config.yaml"

    if not path.exists():
        logger.warning("Config file not found at %s, using defaults.", path)
        return SkillConfig()

    try:
        with open(path, "r", encoding="utf-8") as f:
            raw = yaml.safe_load(f)
    except yaml.YAMLError as exc:
        raise ConfigError(f"Failed to parse YAML config: {exc}") from exc

    if raw is None:
        raw = {}

    try:
        config = SkillConfig(**raw)
    except Exception as exc:
        raise ConfigError(f"Config validation failed: {exc}") from exc

    # 展开路径中的 ~
    config.cache.db_path = os.path.expanduser(config.cache.db_path)
    config.output.default_dir = os.path.expanduser(config.output.default_dir)

    logger.info("Config loaded from %s", path)
    return config
