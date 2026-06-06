"""测试 YAML 配置加载与校验.

覆盖 config.yaml 解析、默认值、路径展开、子系统校验、错误处理。
"""

from pathlib import Path

import pytest
import yaml

from src.config_loader import SkillConfig, load_config
from src.exceptions import ConfigError


class TestSkillConfig:
    """SkillConfig Pydantic 模型测试."""

    def test_defaults(self) -> None:
        cfg = SkillConfig()
        assert cfg.vep.source == "rest_api"
        assert cfg.vep.base_url == "https://rest.ensembl.org"
        assert cfg.vep.batch_size == 200
        assert cfg.filter.min_qual == 30
        assert cfg.filter.min_dp == 10
        assert cfg.filter.pass_only is True
        assert cfg.cache.default_ttl_days == 30
        assert cfg.rate_limits.ensembl == 15
        assert cfg.rate_limits.ncbi == 3
        assert cfg.rate_limits.uniprot == 10
        assert cfg.output.formats == ["markdown", "json"]
        assert cfg.reference_info.show_gnomad_af is True
        assert cfg.logging.level == "INFO"

    def test_validate_subsystems_valid(self) -> None:
        cfg = SkillConfig(subsystems=["vision", "hearing"])
        assert cfg.subsystems == ["vision", "hearing"]

    def test_validate_subsystems_invalid(self) -> None:
        with pytest.raises(ConfigError):
            SkillConfig(subsystems=["vision", "invalid"])

    def test_path_expansion_not_applied_in_model(self) -> None:
        """路径中的 ~ 应在 load_config 中展开，不在模型层面."""
        cfg = SkillConfig()
        assert "~" in cfg.cache.db_path
        assert "~" in cfg.output.default_dir


class TestLoadConfig:
    """load_config 函数测试."""

    def test_load_default_config(self) -> None:
        """从项目根目录加载 config.yaml."""
        cfg = load_config()
        assert cfg.vep.source == "rest_api"
        assert cfg.subsystems == ["vision", "hearing", "olfaction", "taste", "somatosensation"]
        # 路径已展开
        assert "~" not in cfg.cache.db_path
        assert "~" not in cfg.output.default_dir

    def test_load_specific_path(self) -> None:
        cfg = load_config("config.yaml")
        assert cfg.vep.batch_size == 200

    def test_load_nonexistent_returns_defaults(self) -> None:
        cfg = load_config("/nonexistent/path/config.yaml")
        assert cfg.vep.source == "rest_api"

    def test_load_invalid_yaml(self, tmp_path: Path) -> None:
        bad_yaml = tmp_path / "bad.yaml"
        bad_yaml.write_text("{invalid yaml:::")
        with pytest.raises(ConfigError):
            load_config(str(bad_yaml))

    def test_load_empty_yaml(self, tmp_path: Path) -> None:
        empty_yaml = tmp_path / "empty.yaml"
        empty_yaml.write_text("")
        cfg = load_config(str(empty_yaml))
        assert cfg.vep.source == "rest_api"  # 应回退到默认值

    def test_custom_values(self, tmp_path: Path) -> None:
        custom = tmp_path / "custom.yaml"
        data = {
            "vep": {"source": "local", "batch_size": 50},
            "filter": {"min_qual": 50, "pass_only": False},
        }
        custom.write_text(yaml.safe_dump(data))
        cfg = load_config(str(custom))
        assert cfg.vep.source == "local"
        assert cfg.vep.batch_size == 50
        assert cfg.filter.min_qual == 50
        assert cfg.filter.pass_only is False
        # 未指定的字段保持默认值
        assert cfg.filter.min_dp == 10
