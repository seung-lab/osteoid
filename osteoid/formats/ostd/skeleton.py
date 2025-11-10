import numpy as np
import numpy.typing as npt

from .header import (
  OstdHeader, 
  OstdTransform,
  OstdTransformSection,
  OstdAttribute,
  OstdAttributeSection,
)

# represents one skeleton section
# in a possibly multipart file
class OstdSkeletonPart:
  pass

# represents a full skeleton including
# multiple parts
class OstdSkeleton:
  def __init__(self):
    self.parts = []

  def append(self, skel:"OstdSkeleton"):
    self.parts.append(skel)

  def save(self, filename:str):
    pass

  @classmethod
  def load(kls, filename:str) -> "OstdSkeleton":
    pass

  def to_bytes(self) -> bytes:
    pass

  @classmethod
  def from_bytes(kls, binary:bytes) -> "OstdSkeleton":
    pass