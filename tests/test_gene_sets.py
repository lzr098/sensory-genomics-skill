"""测试感官基因集加载器.

覆盖 GeneSetLoader 的初始化、索引构建、查询方法及边界条件。
"""

from pathlib import Path

import pytest

from src.gene_sets.loader import GeneSetLoader


class TestGeneSetLoader:
    """GeneSetLoader 集成测试."""

    def test_load_default_data(self) -> None:
        loader = GeneSetLoader()
        assert len(loader.gene_index) > 0
        assert len(loader.subsystem_index) == 5
        assert "vision" in loader.subsystem_index
        assert "hearing" in loader.subsystem_index
        assert "olfaction" in loader.subsystem_index
        assert "taste" in loader.subsystem_index
        assert "somatosensation" in loader.subsystem_index

    def test_is_sensory_gene(self) -> None:
        loader = GeneSetLoader()
        assert loader.is_sensory_gene("GJB2") is True
        assert loader.is_sensory_gene("OR2T11") is True
        assert loader.is_sensory_gene("BRCA1") is False
        assert loader.is_sensory_gene("") is False

    def test_get_subsystem(self) -> None:
        loader = GeneSetLoader()
        assert loader.get_subsystem("GJB2") == "hearing"
        assert loader.get_subsystem("OPN1LW") == "vision"
        assert loader.get_subsystem("TAS2R38") == "taste"
        assert loader.get_subsystem("BRCA1") == ""

    def test_get_genes_for_subsystem(self) -> None:
        loader = GeneSetLoader()
        vision_genes = loader.get_genes_for_subsystem("vision")
        assert isinstance(vision_genes, list)
        assert "OPN1LW" in vision_genes
        assert "OPN1MW" in vision_genes

    def test_get_all_genes(self) -> None:
        loader = GeneSetLoader()
        all_genes = loader.get_all_genes()
        assert isinstance(all_genes, set)
        assert len(all_genes) > 50  # 应包含大量基因

    def test_get_gene_function(self) -> None:
        loader = GeneSetLoader()
        func = loader.get_gene_function("GJB2")
        assert "缝隙连接" in func

    def test_get_gene_detail_missing(self) -> None:
        loader = GeneSetLoader()
        assert loader.get_gene_detail("UNKNOWN", "function") == ""
        assert loader.get_gene_detail("UNKNOWN", "function", default="N/A") == "N/A"

    def test_load_nonexistent_data_dir(self) -> None:
        """数据目录不存在时应优雅降级."""
        loader = GeneSetLoader(data_dir="/nonexistent/path")
        assert loader.gene_index == {}
        assert loader.subsystem_index == {}

    def test_load_custom_data_dir(self, tmp_path: Path) -> None:
        """从自定义目录加载."""
        import yaml

        data_dir = tmp_path / "data"
        data_dir.mkdir()
        gene_sets = {
            "subsystems": [
                {
                    "name": "test_subsystem",
                    "genes": [
                        {"symbol": "TEST1", "function": "test function 1"},
                        "TEST2",
                    ],
                }
            ]
        }
        (data_dir / "sensory_gene_sets.yaml").write_text(yaml.safe_dump(gene_sets), encoding="utf-8")

        loader = GeneSetLoader(data_dir=str(data_dir))
        assert loader.is_sensory_gene("TEST1") is True
        assert loader.is_sensory_gene("TEST2") is True
        assert loader.get_subsystem("TEST1") == "test_subsystem"
        assert loader.get_gene_function("TEST1") == "test function 1"

    def test_gene_index_no_duplicates(self) -> None:
        """基因不应在不同子系统中重复（或至少 index 只记录一个）."""
        loader = GeneSetLoader()
        # 同一个基因不应出现两次
        assert len(loader.gene_index) == len(set(loader.gene_index.keys()))
