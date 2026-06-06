"""VEP 注释模块.

提供 Ensembl VEP REST API 异步客户端和本地 VEP 预留接口。
"""

from src.vep.batcher import VariantBatcher
from src.vep.client import VepClient

__all__ = ["VepClient", "VariantBatcher"]
