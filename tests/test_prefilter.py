"""测试 VCF 预过滤规则.

覆盖 QUAL、DP、FILTER 字段的阈值过滤及边界条件。
"""

import pytest

from src.models import Variant
from src.vcf.prefilter import Prefilter


class TestPrefilter:
    """Prefilter 单元测试."""

    def test_pass_all(self) -> None:
        p = Prefilter(min_qual=30, min_dp=10, pass_only=True)
        v = Variant(chrom="1", pos=100, ref="A", alt="T", gt="0/1", qual=50, dp=20, filter_status="PASS")
        assert p.apply(v) is True

    def test_fail_qual(self) -> None:
        p = Prefilter(min_qual=30, min_dp=10, pass_only=True)
        v = Variant(chrom="1", pos=100, ref="A", alt="T", gt="0/1", qual=20, dp=20, filter_status="PASS")
        assert p.apply(v) is False

    def test_fail_dp(self) -> None:
        p = Prefilter(min_qual=30, min_dp=10, pass_only=True)
        v = Variant(chrom="1", pos=100, ref="A", alt="T", gt="0/1", qual=50, dp=5, filter_status="PASS")
        assert p.apply(v) is False

    def test_fail_filter(self) -> None:
        p = Prefilter(min_qual=30, min_dp=10, pass_only=True)
        v = Variant(chrom="1", pos=100, ref="A", alt="T", gt="0/1", qual=50, dp=20, filter_status="LowQual")
        assert p.apply(v) is False

    def test_pass_not_only(self) -> None:
        p = Prefilter(min_qual=30, min_dp=10, pass_only=False)
        v = Variant(chrom="1", pos=100, ref="A", alt="T", gt="0/1", qual=50, dp=20, filter_status="LowQual")
        assert p.apply(v) is True

    def test_boundary_qual(self) -> None:
        p = Prefilter(min_qual=30, min_dp=10, pass_only=True)
        v = Variant(chrom="1", pos=100, ref="A", alt="T", gt="0/1", qual=30, dp=10, filter_status="PASS")
        assert p.apply(v) is True

    def test_boundary_dp(self) -> None:
        p = Prefilter(min_qual=30, min_dp=10, pass_only=True)
        v = Variant(chrom="1", pos=100, ref="A", alt="T", gt="0/1", qual=30, dp=10, filter_status="PASS")
        assert p.apply(v) is True

    def test_apply_with_logging(self) -> None:
        p = Prefilter(min_qual=30, min_dp=10, pass_only=True)
        v = Variant(chrom="1", pos=100, ref="A", alt="T", gt="0/1", qual=20, dp=20, filter_status="PASS")
        assert p.apply_with_logging(v) is False
