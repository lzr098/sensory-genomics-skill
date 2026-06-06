"""专用逻辑模块.

提供 TAS2R38、线粒体耳聋、OR 基因分级等特殊基因的分析逻辑。
"""

from src.specialized.mitochondrial import MitochondrialAnnotator
from src.specialized.or_tiers import ORTierClassifier
from src.specialized.tas2r38 import TAS2R38Analyzer

__all__ = ["TAS2R38Analyzer", "MitochondrialAnnotator", "ORTierClassifier"]
