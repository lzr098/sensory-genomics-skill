"""关键性状 SNP 基因型推断模块.

从 VCF 中查询已知性状关联 SNP，对于 VCF 中不存在的位点推断为 REF/REF，
并根据预定义的基因型-表型映射返回推断结果。
"""

import subprocess
import yaml
from pathlib import Path
from typing import Dict, List, Optional
from pydantic import BaseModel, Field

from src.logger import get_logger
from src.models import KeySNPResult

logger = get_logger(__name__)


class KeySNPInferrer:
    """查询 VCF 中的关键性状 SNP 并推断基因型."""

    def __init__(self, vcf_path: str):
        self.vcf_path = vcf_path
        self.snps = self._load_snps()

    def _load_snps(self) -> Dict:
        """加载关键 SNP 数据库."""
        # 尝试多个路径
        candidates = [
            Path(__file__).resolve().parent.parent.parent / "assets" / "data" / "key_trait_snps.yaml",
            Path("assets/data/key_trait_snps.yaml"),
            Path("/Users/zhaorongli/.workbuddy/skills/sensory-genomics/assets/data/key_trait_snps.yaml"),
        ]
        for path in candidates:
            if path.exists():
                with open(path, "r") as f:
                    return yaml.safe_load(f)
        logger.warning("key_trait_snps.yaml not found, using empty SNP set")
        return {"snps": {}}

    def _query_vcf(self, chrom: str, pos: int) -> Optional[Dict]:
        """使用 bcftools 查询 VCF 中特定位置的变异.

        自动处理 chr 前缀差异（尝试 chrom 和 chr+chrom 两种形式）。

        Returns:
            dict with keys: chrom, pos, ref, alt, gt 或 None if not found.
        """
        # 尝试两种染色体命名格式
        chrom_candidates = [chrom]
        if not chrom.startswith("chr"):
            chrom_candidates.append(f"chr{chrom}")
        else:
            # 如果已有 chr 前缀，也尝试不带前缀的
            chrom_candidates.append(chrom[3:])

        for chrom_query in chrom_candidates:
            try:
                cmd = [
                    "bcftools", "view", "-H",
                    "-r", f"{chrom_query}:{pos}",
                    self.vcf_path,
                ]
                result = subprocess.run(
                    cmd, capture_output=True, text=True, timeout=30
                )
                if result.returncode != 0 or not result.stdout.strip():
                    continue

                line = result.stdout.strip().split("\n")[0]
                parts = line.split("\t")
                if len(parts) < 10:
                    continue

                ref = parts[3]
                alt = parts[4].split(",")
                format_fields = parts[8].split(":")
                sample_fields = parts[9].split(":")
                gt_field = sample_fields[0]  # e.g. "0/1" or "0|0"

                # Build dict of FORMAT -> sample value
                fmt = dict(zip(format_fields, sample_fields))

                # Normalize GT to alleles
                gt = self._gt_to_alleles(gt_field, ref, alt)

                # Extract quality metrics
                dp = self._parse_int(fmt.get("DP"))
                gq = self._parse_int(fmt.get("GQ"))
                ad_parts = fmt.get("AD", "").split(",") if fmt.get("AD") else []
                ad_ref = self._parse_int(ad_parts[0]) if len(ad_parts) > 0 else None
                ad_alt = self._parse_int(ad_parts[1]) if len(ad_parts) > 1 else None

                return {
                    "chrom": parts[0],
                    "pos": int(parts[1]),
                    "ref": ref,
                    "alt": alt,
                    "gt": gt,
                    "gt_raw": gt_field,
                    "dp": dp,
                    "gq": gq,
                    "ad_ref": ad_ref,
                    "ad_alt": ad_alt,
                }
            except Exception as e:
                logger.warning("Failed to query VCF for %s:%d: %s", chrom_query, pos, e)
                continue

        return None

    @staticmethod
    def _parse_int(value: Optional[str]) -> Optional[int]:
        """安全解析整数字段，返回 None 表示缺失或无效."""
        if value is None or value == "." or value == "":
            return None
        try:
            return int(value)
        except ValueError:
            return None

    @staticmethod
    def _gt_to_alleles(gt_field: str, ref: str, alt_list: List[str]) -> str:
        """将 VCF GT 字段转换为等位基因字符串（如 'AA', 'AG'）."""
        # Handle phased or unphased
        sep = "|" if "|" in gt_field else "/"
        indices = gt_field.split(sep)

        alleles = []
        for idx in indices:
            if idx == ".":
                alleles.append("?")
            elif idx == "0":
                alleles.append(ref)
            else:
                alt_idx = int(idx) - 1
                if alt_idx < len(alt_list):
                    alleles.append(alt_list[alt_idx])
                else:
                    alleles.append("?")

        return "".join(sorted(alleles))

    def infer_all(self) -> List[KeySNPResult]:
        """查询所有关键 SNP 并返回推断结果."""
        results = []
        for rsid, info in self.snps.get("snps", {}).items():
            chrom = str(info["chrom"])
            pos = info["pos"]
            ref = info["ref"]
            alt = info.get("alt", [])

            vcf_record = self._query_vcf(chrom, pos)

            if vcf_record:
                gt = vcf_record["gt"]
                found_in_vcf = True
                dp = vcf_record.get("dp")
                gq = vcf_record.get("gq")
                ad_ref = vcf_record.get("ad_ref")
                ad_alt = vcf_record.get("ad_alt")
            else:
                # Not in VCF -> homozygous reference
                gt = ref + ref
                found_in_vcf = False
                dp = gq = ad_ref = ad_alt = None

            # Look up phenotype
            pmap = info.get("phenotype_map", {})
            # Build a normalized lookup map with sorted genotype keys
            sorted_pmap = { "".join(sorted(k)): v for k, v in pmap.items() }
            pheno = sorted_pmap.get("".join(sorted(gt)), {})
            if not pheno:
                # Fallback: try original key (for backward compatibility)
                pheno = pmap.get(gt, {})
            if not pheno:
                logger.debug("No exact phenotype for %s genotype %s", rsid, gt)
                pheno = {
                    "label": "未知表型",
                    "description": f"基因型 {gt} 的表型数据未定义",
                }

            # Determine zygosity
            unique_alleles = set(gt)
            is_hom_ref = len(unique_alleles) == 1 and ref in unique_alleles
            is_hom_alt = len(unique_alleles) == 1 and ref not in unique_alleles
            is_het = len(unique_alleles) == 2

            result = KeySNPResult(
                rsid=rsid,
                gene=info["gene"],
                chrom=chrom,
                pos=pos,
                ref=ref,
                alt=alt,
                inferred_genotype=gt,
                is_heterozygous=is_het,
                is_homozygous_alt=is_hom_alt,
                is_homozygous_ref=is_hom_ref,
                phenotype_label=pheno.get("label", "未知"),
                phenotype_description=pheno.get("description", ""),
                notes=info.get("notes", ""),
                found_in_vcf=found_in_vcf,
                dp=dp,
                gq=gq,
                ad_ref=ad_ref,
                ad_alt=ad_alt,
            )
            results.append(result)

        return results

    def get_by_gene(self, gene_symbol: str) -> List[KeySNPResult]:
        """获取特定基因的所有关键 SNP 结果."""
        all_results = self.infer_all()
        return [r for r in all_results if r.gene == gene_symbol]
