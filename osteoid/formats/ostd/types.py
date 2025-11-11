from dataclasses import dataclass
from enum import IntEnum

import numpy as np

SI_PREFIXES = {}
LENGTH_SYMBOLS = {}

class SIPrefixType(IntEnum):
  NONE = 0 
  QUECTO = 1 
  RONTO = 2 
  YOCTO = 3 
  ZEPTO = 4 
  ATTO = 5 
  FEMTO = 6 
  PICO = 7 
  NANO = 8 
  MICRO = 9 
  MILLI = 10
  CENTI = 11
  DECI = 12
  DEKA = 13
  HECTO = 14
  KILO = 15
  MEGA = 16
  GIGA = 17
  TERA = 18
  PETA = 19
  EXA = 20
  ZETTA = 21
  YOTTA = 22
  RONNA = 23
  QUETTA = 24

  def __str__(self):
    return SI_PREFIXES[self.value]

SI_PREFIXES = {
  SIPrefixType.NONE: "",
  SIPrefixType.QUECTO: "q",
  SIPrefixType.RONTO: "r",
  SIPrefixType.YOCTO: "y",
  SIPrefixType.ZEPTO: "z",
  SIPrefixType.ATTO: "a",
  SIPrefixType.FEMTO: "f",
  SIPrefixType.PICO: "p",
  SIPrefixType.NANO: "n",
  SIPrefixType.MICRO: "u",
  SIPrefixType.MILLI: "m",
  SIPrefixType.CENTI: "c",
  SIPrefixType.DECI: "d",
  SIPrefixType.DEKA: "D",
  SIPrefixType.HECTO: "h",
  SIPrefixType.KILO: "k",
  SIPrefixType.MEGA: "M",
  SIPrefixType.GIGA: "G",
  SIPrefixType.TERA: "T",
  SIPrefixType.PETA: "P",
  SIPrefixType.EXA: "E",
  SIPrefixType.ZETTA: "Z",
  SIPrefixType.YOTTA: "Y",
  SIPrefixType.RONNA: "R",
  SIPrefixType.QUETTA: "Q",
}

class LengthType(IntEnum):
  VOXEL = 0
  METER = 1
  ANGSTROM = 2
  ASTRONOMICAL_UNIT = 3
  LIGHTYEAR = 4
  PARSEC = 5
  MIL = 6
  INCH = 7
  FOOT = 8
  YARD = 9
  STATUTE_MILE = 10
  NAUTICAL_MILE = 11

  def __str__(self):
    return LENGTH_SYMBOLS[self.value]

class AreaType(IntEnum):
  VOXEL = 0
  METER = 1
  ANGSTROM = 2
  ASTRONOMICAL_UNIT = 3
  LIGHTYEAR = 4
  PARSEC = 5
  MIL = 6
  INCH = 7
  FOOT = 8
  YARD = 9
  STATUTE_MILE = 10
  NAUTICAL_MILE = 11

  def __str__(self):
    return f"{LENGTH_SYMBOLS[self.value]}^2"

class VolumeType(IntEnum):
  VOXEL = 0
  METER = 1
  ANGSTROM = 2
  ASTRONOMICAL_UNIT = 3
  LIGHTYEAR = 4
  PARSEC = 5
  MIL = 6
  INCH = 7
  FOOT = 8
  YARD = 9
  STATUTE_MILE = 10
  NAUTICAL_MILE = 11
  LITER = 12

  def __str__(self):
    if self.value == self.LITER:
      return "L"
    return f"{LENGTH_SYMBOLS[self.value]}^3"

LENGTH_SYMBOLS = {
  LengthType.VOXEL: "vx",
  LengthType.METER: "m",
  LengthType.ANGSTROM: "Å",
  LengthType.ASTRONOMICAL_UNIT: "amu",
  LengthType.LIGHTYEAR: "ly",
  LengthType.PARSEC: "pc",
  LengthType.MIL: "mil",
  LengthType.INCH: "in",
  LengthType.FOOT: "ft",
  LengthType.YARD: "yd",
  LengthType.STATUTE_MILE: "mi",
  LengthType.NAUTICAL_MILE: "nmi",
}

TEMPERATURE_SYMBOLS = {}

class TemperatureType(IntEnum):
  UNKNOWN = 0
  CELSIUS = 1
  FAHRENHEIT = 2
  RANKINE = 3
  KELVIN = 4

  def __str__(self):
    return TEMPERATURE_SYMBOLS[self.value]

TEMPERATURE_SYMBOLS = {
  TemperatureType.UNKNOWN: "",
  TemperatureType.CELSIUS: "°C",
  TemperatureType.FAHRENHEIT: "°F",
  TemperatureType.RANKINE: "°R",
  TemperatureType.KELVIN: "K",
}

class TimeType(IntEnum):
  UNKNOWN = 0
  SECOND = 1
  MINUTE = 2
  HOUR = 3
  DAY = 4
  MONTH = 5
  YEAR = 6
  HERTZ = 7

class LuminosityType(IntEnum):
  UNKNOWN = 0
  CANDELA = 1
  LUMEN = 2
  LUX = 3
  PHOTON = 4
  PHOTONS_PER_SECOND = 5

class ElectricalType(IntEnum):
  UNKNOWN = 0
  VOLT = 1
  AMPERE = 2
  OHM = 3
  SIEMEN = 4
  FARAD = 5
  HENRY = 6
  COULOMB = 7

class MassType(IntEnum):
  UNKNOWN = 0
  GRAM = 1
  DALTON = 2

class SubstanceAmount(IntEnum):
  UNKNOWN = 0
  MOLE = 1

class EnergyType(IntEnum):
  UNKNOWN = 0
  JOULE = 1
  WATT = 2

class CompressionType(IntEnum):
  NONE = 0
  GZIP = 1
  BZIP2 = 2
  ZSTD = 3
  DRACO = 4

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
  BOOL = 12
  PACKED_BOOL = 13

class EdgeRepresentationType(IntEnum):
  PAIR = 0
  PARENT = 1
  LINKED_PATHS = 2

class GraphType(IntEnum):
  GRAPH = 0
  TREE = 1

class SpaceType(IntEnum):
  GENERIC = 0
  PHYSICAL = 1
  SCANNER = 2
  ATLAS = 3
  ALIGNED = 4
  WORLD = 5
  SOMA = 6
  BASE = 7
  JOINT = 8
  TOOL = 9
  MODEL = 10
  CAMERA = 11

@dataclass
class CoordinateFrame:
  sign_x:bool
  sign_y:bool
  sign_z:bool
  permutation:AxisPermutationType

class AxisPermutationType(IntEnum):
  XYZ = 0
  XZY = 1
  YXZ = 2
  YZX = 3
  ZXY = 4
  ZYX = 5
  XY = 6
  YX = 7

class AttributeType(IntEnum):
  VERTEX = 0
  EDGE = 1

TO_AXIS_PERMUTATION = {
  'XYZ': AxisPermutationType.XYZ,
  'XZY': AxisPermutationType.XZY,
  'YXZ': AxisPermutationType.YXZ,
  'YXZ': AxisPermutationType.YZX,
  'ZXY': AxisPermutationType.ZXY,
  'ZYX': AxisPermutationType.ZYX,
  'XY': AxisPermutationType.XY,
  'YX': AxisPermutationType.YX,
}


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
  np.bool_: DataType.BOOL,

}
FROM_DATATYPE = { v:k for k,v in TO_DATATYPE.items() }

QUANTITY_TYPE = {
  0: LengthType,
  1: AreaType,
  2: VolumeType,
  3: TemperatureType,
  4: TimeType,
  5: LuminosityType,
  6: MassType,
  7: ElectricalType,
  8: SubstanceAmount,
  9: EnergyType,
}

def human_readable(prefix:SIPrefixType, unit:IntEnum) -> str:
  return f"{prefix}{unit}"

TO_LENGTH_UNIT = {
  "vx": (SIPrefixType.NONE, LengthType.VOXEL),
  "voxel": (SIPrefixType.NONE, LengthType.VOXEL),

  "A": (SIPrefixType.NONE, LengthType.ANGSTROM),
  "angstrom": (SIPrefixType.NONE, LengthType.ANGSTROM),
  
  "fm": (SIPrefixType.FEMTO, LengthType.METER),
  "femtometer": (SIPrefixType.FEMTO, LengthType.METER),
  
  "pm": (SIPrefixType.PICO, LengthType.METER),
  "picometer": (SIPrefixType.PICO, LengthType.METER),
  
  "nm": (SIPrefixType.NANO, LengthType.METER),
  "nanometer": (SIPrefixType.NANO, LengthType.METER),
  
  "um": (SIPrefixType.MICRO, LengthType.METER),
  "micrometer": (SIPrefixType.MICRO, LengthType.METER),
  "micron": (SIPrefixType.MICRO, LengthType.METER),
  
  "mm": (SIPrefixType.MILLI, LengthType.METER),
  "millimeter": (SIPrefixType.MILLI, LengthType.METER),

  "cm": (SIPrefixType.CENTI, LengthType.METER),
  "centimeter": (SIPrefixType.CENTI, LengthType.METER),
  
  "m": (SIPrefixType.NONE, LengthType.METER),
  "meter": (SIPrefixType.NONE, LengthType.METER),
  
  "km": (SIPrefixType.KILO, LengthType.METER),
  "kilometer": (SIPrefixType.KILO, LengthType.METER),

  "Mm": (SIPrefixType.MEGA, LengthType.METER),
  
  "ly": (SIPrefixType.NONE, LengthType.LIGHTYEAR),
  "lightyear": (SIPrefixType.NONE, LengthType.LIGHTYEAR),

  "pc": (SIPrefixType.NONE, LengthType.PARSEC),
  "parsec": (SIPrefixType.NONE, LengthType.PARSEC),

  "mil": (SIPrefixType.NONE, LengthType.MIL),
  
  "in": (SIPrefixType.NONE, LengthType.INCH),
  "inch": (SIPrefixType.NONE, LengthType.INCH),
  "inches": (SIPrefixType.NONE, LengthType.INCH),
  
  "ft": (SIPrefixType.NONE, LengthType.FOOT),
  "foot": (SIPrefixType.NONE, LengthType.FOOT),
  "feet": (SIPrefixType.NONE, LengthType.FOOT),

  "yd": (SIPrefixType.NONE, LengthType.YARD),
  "yard": (SIPrefixType.NONE, LengthType.YARD),
  
  "mi": (SIPrefixType.NONE, LengthType.STATUTE_MILE),
  "mile": (SIPrefixType.NONE, LengthType.STATUTE_MILE),
  
  "nmi": (SIPrefixType.NONE, LengthType.NAUTICAL_MILE),
}
FROM_LENGTH_UNIT = {
  v:k for k,v in TO_LENGTH_UNIT.items()
}





