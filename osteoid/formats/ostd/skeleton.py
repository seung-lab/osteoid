import numpy as np

from .header import (
	OstdHeader, 
	OstdTransform,
	OstdTransformSection,
	OstdAttribute,
	OstdAttributeSection,
)


class OstdSkeleton:
	def __init__(self):
		self.header = None
		self.vertices = None
		self.edges = None
		self.attributes = {}

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