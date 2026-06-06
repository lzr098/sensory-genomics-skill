"""API 富集与缓存模块.

提供外部生物信息学 API 的统一查询和本地缓存功能。
"""

from src.enrichment.cache import CacheManager

__all__ = ["CacheManager"]
