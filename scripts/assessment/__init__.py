"""功能影响评估模块.

提供基于规则的感官基因功能影响评估引擎。
"""

from src.assessment.engine import ImpactEngine
from src.assessment.inheritance import InheritanceMatcher
from src.assessment.rules import GeneCertaintyRule, ImpactRule, ProteinImpactRule

__all__ = [
    "ImpactEngine",
    "InheritanceMatcher",
    "ImpactRule",
    "ProteinImpactRule",
    "GeneCertaintyRule",
]
