"""RNA-only pocket prediction components extracted from RiboPoseDiff."""

from __future__ import annotations

class SampleValidationError(ValueError):
    """An explicit scientific rejection of an individual input sample."""

    def __init__(self, message: str, *, stage: str = "sample_validation", details=None):
        super().__init__(message)
        self.stage = stage
        self.details = details or {}
