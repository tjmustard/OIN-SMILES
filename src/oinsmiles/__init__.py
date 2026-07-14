import logging as _logging

from .core.chirality import OINStereoWarning
from .core.translator import SMILESToXYZ, XYZToSMILES

# Library logging best practice: attach a NullHandler to the package logger so
# importing oinsmiles never emits "No handler found" warnings and never writes to
# the consumer's stderr on its own. Applications opt in to output by configuring
# logging themselves (e.g. logging.basicConfig(level=logging.DEBUG)).
_logging.getLogger(__name__).addHandler(_logging.NullHandler())

__all__ = ["XYZToSMILES", "SMILESToXYZ", "OINStereoWarning"]
