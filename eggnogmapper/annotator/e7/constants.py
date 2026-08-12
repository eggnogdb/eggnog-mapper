"""Shared constants for the eggNOG v7 annotation engine.

Extracted here so lightweight consumers (e.g. eggnogmapper.annotation.output)
can import these values without pulling in the full AnnotationEngine.
"""

# Canonical mapping from cascade tier integer to human-readable confidence
# label. Single source of truth for tier names. AnnotationEngine.TIER_CONFIDENCE
# must be assigned from this dict, never independently redefined.
TIER_CONFIDENCE: dict[int, str] = {0: "high", 1: "medium", 2: "low"}
