"""测试自定义异常体系.

覆盖所有异常类型的构造、属性继承及消息传递。
"""

import pytest

from src.exceptions import (
    ApiError,
    AssessmentError,
    ConfigError,
    SensoryGenomicsError,
    VcfError,
    VepError,
)


class TestExceptions:
    """异常体系测试."""

    def test_base_exception(self) -> None:
        exc = SensoryGenomicsError("base error")
        assert str(exc) == "base error"
        assert exc.message == "base error"

    def test_vcf_error(self) -> None:
        exc = VcfError("parse failed", record="chr1:100:A>T")
        assert "parse failed" in str(exc)
        assert exc.record == "chr1:100:A>T"
        assert isinstance(exc, SensoryGenomicsError)

    def test_vep_error(self) -> None:
        exc = VepError("timeout", status_code=504)
        assert exc.status_code == 504
        assert isinstance(exc, SensoryGenomicsError)

    def test_config_error(self) -> None:
        exc = ConfigError("invalid field", field="subsystems")
        assert exc.field == "subsystems"
        assert isinstance(exc, SensoryGenomicsError)

    def test_api_error(self) -> None:
        exc = ApiError("rate limited", api_name="gnomAD", status_code=429)
        assert exc.api_name == "gnomAD"
        assert exc.status_code == 429
        assert isinstance(exc, SensoryGenomicsError)

    def test_assessment_error(self) -> None:
        exc = AssessmentError("unknown gene", gene_symbol="UNKNOWN")
        assert exc.gene_symbol == "UNKNOWN"
        assert isinstance(exc, SensoryGenomicsError)

    def test_exception_catch_base(self) -> None:
        """所有具体异常都应能被基类捕获."""
        exceptions = [
            VcfError("test"),
            VepError("test"),
            ConfigError("test"),
            ApiError("test"),
            AssessmentError("test"),
        ]
        for exc in exceptions:
            with pytest.raises(SensoryGenomicsError):
                raise exc
