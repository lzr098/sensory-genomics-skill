"""End-to-end test for sensory-genomics pipeline.

Runs the full pipeline on a small pre-filtered VCF and verifies:
- Pipeline completes without exceptions
- Markdown and JSON outputs are generated and non-empty
- Markdown contains key expected sections (disclaimer, CDH23, etc.)
- No NameError or traceback remnants in output
"""

import asyncio
import os
import sys
import tempfile
from pathlib import Path

# Add scripts to path so 'src' is importable
scripts_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(scripts_dir))

from src.main import run_analysis
from src.models import AnalysisConfig
from src.report.markdown_generator import MarkdownReportGenerator
from src.report.json_generator import JsonReportGenerator


VCF_PATH = "/Users/zhaorongli/WorkBuddy/2026-06-06-21-14-49/LSP_WGS_20260417.filtered2.sensory.vcf.gz"


def test_pipeline_completes():
    """Test that the full pipeline runs without exception."""
    config = AnalysisConfig(
        vcf_path=VCF_PATH,
        sex="M",
        subsystems=["vision", "hearing", "olfaction", "taste", "somatosensation"],
        output_dir=tempfile.mkdtemp(),
        show_reference_info=False,  # skip API enrichment for speed
        strict_filter=False,
    )

    report = asyncio.run(run_analysis(config))

    assert report is not None
    assert report.sample_id == "E-S22506309071"
    assert report.sex == "M"
    assert len(report.gene_cards) > 0
    assert report.tas2r38 is not None
    print(f"✅ Pipeline completed: {len(report.gene_cards)} gene cards")
    return report


def test_markdown_generation(report):
    """Test that Markdown report generates without NameError."""
    md_gen = MarkdownReportGenerator()
    md_content = md_gen.generate(report)

    assert md_content
    assert len(md_content) > 1000
    assert "NameError" not in md_content
    assert "Traceback" not in md_content
    print(f"✅ Markdown generated: {len(md_content)} chars")
    return md_content


def test_json_generation(report):
    """Test that JSON report generates successfully."""
    json_gen = JsonReportGenerator()
    json_content = json_gen.generate(report)

    assert json_content
    assert len(json_content) > 1000
    assert "error" not in json_content.lower() or "\"error\"" not in json_content.lower()
    print(f"✅ JSON generated: {len(json_content)} chars")
    return json_content


def test_markdown_content(md_content: str):
    """Verify key sections are present in the Markdown report."""
    # Disclaimer
    assert "重要声明" in md_content, "Missing disclaimer"
    # CDH23 with specific follow-up advice
    assert "CDH23" in md_content, "Missing CDH23"
    assert "OAE" in md_content or "耳声发射" in md_content, "Missing OAE follow-up advice"
    assert "ABR" in md_content or "听性脑干反应" in md_content, "Missing ABR follow-up advice"
    # OPRM1
    assert "OPRM1" in md_content, "Missing OPRM1"
    # OR section 2.3 should show count, not a table
    assert "2.3" in md_content, "Missing 2.3 section"
    # Appendix G should contain the OR homozygous LoF table
    assert "G. 配体未知的纯合 OR 功能丧失汇总" in md_content, "Missing Appendix G"
    # Subsystem sections
    assert "视觉系统" in md_content, "Missing vision section"
    assert "听觉系统" in md_content, "Missing hearing section"
    print("✅ Markdown content validation passed")


def test_api_on_demand(report):
    """Verify that API enrichment was skipped (show_reference_info=False)."""
    # When show_reference_info=False, enrichment_data should be empty
    for card in report.gene_cards:
        assert not card.enrichment_data, f"Unexpected enrichment data for {card.gene_symbol}"
    print("✅ API on-demand behavior verified (skipped when disabled)")


def main():
    print("=" * 60)
    print("Sensory Genomics End-to-End Test")
    print("=" * 60)

    if not os.path.exists(VCF_PATH):
        print(f"❌ Test VCF not found: {VCF_PATH}")
        sys.exit(1)

    # 1. Pipeline
    report = test_pipeline_completes()

    # 2. Markdown generation (critical: must not hit NameError)
    md_content = test_markdown_generation(report)

    # 3. JSON generation
    json_content = test_json_generation(report)

    # 4. Content validation
    test_markdown_content(md_content)

    # 5. API on-demand behavior
    test_api_on_demand(report)

    print("=" * 60)
    print("ALL TESTS PASSED ✅")
    print("=" * 60)


if __name__ == "__main__":
    main()
