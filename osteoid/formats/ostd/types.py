from enum import IntEnum

class PhysicalUnit(IntEnum):
  VOXEL = 0

  ANGSTROM = 1
  FEMTOMETER = 2
  PICOMETER = 3
  NANOMETER = 4
  MICROMETER = 5
  MILLIMETER = 6
  CENTIMETER = 7
  METER = 8
  KILOMETER = 9
  MEGAMETER = 10
  LIGHTYEAR = 11
  PARSEC = 12

  MIL = 13
  INCH = 14
  FOOT = 15
  YARD = 16
  STATUTE_MILE = 17
  NAUTICAL_MILE = 18

TO_UNIT = {
  "vx": PhysicalUnit.VOXEL,
  "voxel": PhysicalUnit.VOXEL,

  "A": PhysicalUnit.ANGSTROM,
  "angstrom": PhysicalUnit.ANGSTROM,
  
  "fm": PhysicalUnit.FEMTOMETER,
  "femtometer": PhysicalUnit.FEMTOMETER,
  
  "pm": PhysicalUnit.PICOMETER,
  "picometer": PhysicalUnit.PICOMETER,
  
  "nm": PhysicalUnit.NANOMETER,
  "nanometer": PhysicalUnit.NANOMETER,
  
  "um": PhysicalUnit.MICROMETER,
  "micrometer": PhysicalUnit.MICROMETER,
  "micron": PhysicalUnit.MICROMETER,
  
  "mm": PhysicalUnit.MILLIMETER,
  "millimeter": PhysicalUnit.MILLIMETER,

  "cm": PhysicalUnit.CENTIMETER,
  "centimeter": PhysicalUnit.CENTIMETER,
  
  "m": PhysicalUnit.METER,
  "meter": PhysicalUnit.METER,
  
  "km": PhysicalUnit.KILOMETER,
  "kilometer": PhysicalUnit.KILOMETER,

  "Mm": PhysicalUnit.MEGAMETER,
  
  "ly": PhysicalUnit.LIGHTYEAR,
  "lightyear": PhysicalUnit.LIGHTYEAR,

  "pc": PhysicalUnit.PARSEC,
  "parsec": PhysicalUnit.PARSEC,

  "mil": PhysicalUnit.MIL,
  
  "in": PhysicalUnit.INCH,
  "inch": PhysicalUnit.INCH,
  "inches": PhysicalUnit.INCH,
  
  "ft": PhysicalUnit.FOOT,
  "foot": PhysicalUnit.FOOT,
  "feet": PhysicalUnit.FOOT,

  "yd": PhysicalUnit.YARD,
  "yard": PhysicalUnit.YARD,
  
  "mi": PhysicalUnit.STATUTE_MILE,
  "mile": PhysicalUnit.STATUTE_MILE,
  
  "nmi": PhysicalUnit.NAUTICAL_MILE,
}

class CompressionType(IntEnum):
  NONE = 0
  GZIP = 1
  BZIP2 = 2
  ZSTD = 3

class DataType(IntEnum):
  F8 = 0
  F16 = 1
  F32 = 2
  F64 = 3
  U8 = 4
  U16 = 5
  U32 = 6
  U64 = 7
  I8 = 8
  I16 = 9
  I32 = 10
  I64 = 11
  BOOLEAN = 12
  PACKED_BOOLEAN = 13

TO_DATATYPE = {
  np.float16: DataType.F16,
  np.float32: DataType.F32,
  np.float64: DataType.F64,
  np.uint8: DataType.U8,
  np.uint16: DataType.U16,
  np.uint32: DataType.U32,
  np.uint64: DataType.U64,
  np.int8: DataType.I8,
  np.int16: DataType.I16,
  np.int32: DataType.I32,
  np.int64: DataType.I64,
}

FROM_DATATYPE = { v,k for k,v in TO_DATATYPE.items() }

class EdgeRepresentation(IntEnum):
  PAIR = 0
  PARENT = 1

class GraphType(IntEnum):
  GRAPH = 0
  TREE = 1

class SpaceType(IntEnum):
  VOXEL = 0
  PHYSICAL = 1
