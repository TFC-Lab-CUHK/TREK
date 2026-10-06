"""TREK: transcription end site and alternative polyadenylation identification from long reads."""

__version__ = "1.0.0"
__author__ = "Jizhou Zhang"

from .gtf_processor import GTFProcessor
from .alignment_processor import AlignmentProcessor
from .apa_finder import TESFinder

__all__ = ['GTFProcessor', 'AlignmentProcessor', 'TESFinder']
