from dataclasses import dataclass
from enum import IntEnum

import numpy as np

SI_PREFIXES = {}
LENGTH_SYMBOLS = {}

class SIPrefixType(IntEnum):
  ZEPTO = 0
  ATTO = 1 
  FEMTO = 2 
  PICO = 3
  NANO = 4 
  MICRO = 5
  MILLI = 6
  CENTI = 7
  NONE = 8
  KILO = 9
  MEGA = 10
  GIGA = 11
  TERA = 12
  PETA = 13
  EXA = 14
  ZETTA = 15

  def __str__(self):
    return SI_PREFIXES[self.value]

SI_PREFIXES = {
  SIPrefixType.NONE: "",
  SIPrefixType.ZEPTO: "z",
  SIPrefixType.ATTO: "a",
  SIPrefixType.FEMTO: "f",
  SIPrefixType.PICO: "p",
  SIPrefixType.NANO: "n",
  SIPrefixType.MICRO: "u",
  SIPrefixType.MILLI: "m",
  SIPrefixType.CENTI: "c",
  SIPrefixType.KILO: "k",
  SIPrefixType.MEGA: "M",
  SIPrefixType.GIGA: "G",
  SIPrefixType.TERA: "T",
  SIPrefixType.PETA: "P",
  SIPrefixType.EXA: "E",
  SIPrefixType.ZETTA: "Z",
}

SI_PREFIX_VALUE = {
  SIPrefixType.ZEPTO: 1e-21,
  SIPrefixType.ATTO: 1e-18,
  SIPrefixType.FEMTO: 1e-15,
  SIPrefixType.PICO: 1e-12,
  SIPrefixType.NANO: 1e-9,
  SIPrefixType.MICRO: 1e-6,
  SIPrefixType.MILLI: 0.001,
  SIPrefixType.CENTI: 0.01,
  SIPrefixType.NONE: 1.0,
  SIPrefixType.KILO: 1000.0,
  SIPrefixType.MEGA: 1e6,
  SIPrefixType.GIGA: 1e9,
  SIPrefixType.TERA: 1e12,
  SIPrefixType.PETA: 1e15,
  SIPrefixType.EXA: 1e18,
  SIPrefixType.ZETTA: 1e21,
}

class LengthType(IntEnum):
  UNKNOWN = 0
  VOXEL = 1
  ANGSTROM = 2
  MIL = 3
  INCH = 4
  FOOT = 5
  YARD = 6
  METER = 7
  STATUTE_MILE = 8
  NAUTICAL_MILE = 9
  ASTRONOMICAL_UNIT = 10
  LIGHTYEAR = 11
  PARSEC = 12

LENGTH_SYMBOLS = {
  LengthType.UNKNOWN: "",
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

LengthType.__str__ = lambda self: LENGTH_SYMBOLS[self]

LENGTH_CONVERSION_FACTORS = {
  (LengthType.UNKNOWN, LengthType.METER): float('NaN'),
  (LengthType.VOXEL, LengthType.METER): float('NaN'),
  (LengthType.METER, LengthType.METER): 1.0,
  (LengthType.ANGSTROM, LengthType.METER): 1e-10,
  (LengthType.ASTRONOMICAL_UNIT, LengthType.METER): 149597870691,
  (LengthType.LIGHTYEAR, LengthType.METER): 9.460528405e15,
  (LengthType.PARSEC, LengthType.METER): 3.0856778570831e16,
  (LengthType.MIL, LengthType.METER): 0.0000254,
  (LengthType.INCH, LengthType.METER): 0.0254,
  (LengthType.FOOT, LengthType.METER): 0.3048,
  (LengthType.YARD, LengthType.METER): 0.9144,
  (LengthType.STATUTE_MILE, LengthType.METER): 1609.344,
  (LengthType.NAUTICAL_MILE, LengthType.METER): 1852,
}

def length_conversion_factor(unit1:LengthType, unit2:LengthType) -> float:
  if unit1 == unit2:
    return 1.0

  u1_meters = LENGTH_CONVERSION_FACTORS[(unit1, LengthType.METER)]
  u2_meters = LENGTH_CONVERSION_FACTORS[(unit2, LengthType.METER)]

  return u2_meters / u1_meters

class AreaType(IntEnum):
  UNKNOWN = 0
  VOXEL = 1
  ANGSTROM = 2
  MIL = 3
  INCH = 4
  FOOT = 5
  YARD = 6
  METER = 7
  STATUTE_MILE = 8
  NAUTICAL_MILE = 9
  ASTRONOMICAL_UNIT = 10
  LIGHTYEAR = 11
  PARSEC = 12

AreaType.__str__ = lambda self: f"{LENGTH_SYMBOLS[self]}^2"

class VolumeType(IntEnum):
  UNKNOWN = 0
  VOXEL = 1
  ANGSTROM = 2
  MIL = 3
  INCH = 4
  FOOT = 5
  YARD = 6
  METER = 7
  STATUTE_MILE = 8
  NAUTICAL_MILE = 9
  ASTRONOMICAL_UNIT = 10
  LIGHTYEAR = 11
  PARSEC = 12
  LITER = 13

VOLUME_SYMBOLS = {**LENGTH_SYMBOLS, VolumeType.LITER: "L"}
VolumeType.__str__ = lambda self: "L" if self == VolumeType.LITER else f"{LENGTH_SYMBOLS[self]}^3"

class TemperatureType(IntEnum):
  UNKNOWN = 0
  CELSIUS = 1
  FAHRENHEIT = 2
  RANKINE = 3
  KELVIN = 4

TEMPERATURE_SYMBOLS = {
  TemperatureType.UNKNOWN: "",
  TemperatureType.CELSIUS: "°C",
  TemperatureType.FAHRENHEIT: "°F",
  TemperatureType.RANKINE: "°R",
  TemperatureType.KELVIN: "K",
}

TemperatureType.__str__ = lambda self: TEMPERATURE_SYMBOLS[self]

class TimeType(IntEnum):
  UNKNOWN = 0
  SECOND = 1
  MINUTE = 2
  HOUR = 3
  DAY = 4
  MONTH = 5
  YEAR = 6
  HERTZ = 7

TIME_SYMBOLS = {
  TimeType.UNKNOWN: "",
  TimeType.SECOND: "s",
  TimeType.MINUTE: "min",
  TimeType.HOUR: "h",
  TimeType.DAY: "d",
  TimeType.MONTH: "mo",
  TimeType.YEAR: "y",
  TimeType.HERTZ: "Hz",
}

TimeType.__str__ = lambda self: TIME_SYMBOLS[self]

class LuminosityType(IntEnum):
  UNKNOWN = 0
  CANDELA = 1
  LUMEN = 2
  LUX = 3
  PHOTON = 4
  PHOTONS_PER_SECOND = 5

LUMINOSITY_SYMBOLS = {
  LuminosityType.UNKNOWN: "",
  LuminosityType.CANDELA: "cd",
  LuminosityType.LUMEN: "lm",
  LuminosityType.LUX: "lx",
  LuminosityType.PHOTON: "photons",
  LuminosityType.PHOTONS_PER_SECOND: "pps",
}

LuminosityType.__str__ = lambda self: LUMINOSITY_SYMBOLS[self]


class ElectricalType(IntEnum):
  UNKNOWN = 0
  VOLT = 1
  AMPERE = 2
  OHM = 3
  SIEMEN = 4
  FARAD = 5
  HENRY = 6
  COULOMB = 7

ELECTRICAL_SYMBOLS = {
  ElectricalType.UNKNOWN: "",
  ElectricalType.VOLT: "V",
  ElectricalType.AMPERE: "A",
  ElectricalType.OHM: "Ω",
  ElectricalType.SIEMEN: "S",
  ElectricalType.FARAD: "F",
  ElectricalType.HENRY: "H",
  ElectricalType.COULOMB: "C",
}

ElectricalType.__str__ = lambda self: ELECTRICAL_SYMBOLS[self]

class MassType(IntEnum):
  UNKNOWN = 0
  GRAM = 1
  DALTON = 2

MASS_SYMBOLS = {
  MassType.UNKNOWN: "",
  MassType.GRAM: "g",
  MassType.DALTON: "da",
}

MassType.__str__ = lambda self: MASS_SYMBOLS[self]

class SubstanceAmount(IntEnum):
  UNKNOWN = 0
  MOLE = 1

SUBSTANCE_AMOUNT_SYMBOLS = {
  SubstanceAmount.UNKNOWN: "",
  SubstanceAmount.MOLE: "mol",
}

SubstanceAmount.__str__ = lambda self: SUBSTANCE_AMOUNT_SYMBOLS[self]

class EnergyType(IntEnum):
  UNKNOWN = 0
  JOULE = 1
  WATT = 2

ENERGY_SYMBOLS = {
  EnergyType.UNKNOWN: "",
  EnergyType.JOULE: "J",
  EnergyType.WATT: "W",
}

EnergyType.__str__ = lambda self: ENERGY_SYMBOLS[self]

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
  CYCLIC = 2

class SpaceType(IntEnum):
  GENERIC = 0
  VOXEL = 1
  PHYSICAL = 2
  SCANNER = 3
  ATLAS = 4
  ALIGNED = 5
  WORLD = 6
  SOMA = 7
  BASE = 8
  JOINT = 9
  TOOL = 10
  MODEL = 11
  CAMERA = 12

class AxisPermutationType(IntEnum):
  XYZ = 0
  XZY = 1
  YXZ = 2
  YZX = 3
  ZXY = 4
  ZYX = 5
  XY = 6
  YX = 7

@dataclass
class CoordinateFrame:
  sign_x:bool
  sign_y:bool
  sign_z:bool
  permutation:AxisPermutationType

  def __str__(self) -> str:
    def sgn(x):
      return "+" if not x else "-"

    signs = {
      'X': sgn(self.sign_x),
      'Y': sgn(self.sign_y),
      'Z': sgn(self.sign_z),
    }

    axes = FROM_AXIS_PERMUTATION[self.permutation]
    out = ""
    for i in range(len(axes)):
      out += f"{signs[axes[i]]}{axes[i]}"
    return out

  @classmethod
  def parse(kls, orientation:str) -> "CoordinateFrame":
    if len(orientation) > 6:
      raise ValueError(f"Unable to parse orientation: {orientation[:100]}")

    orientation = orientation.upper()
    normalized = orientation.replace('+', '').replace('-', '')

    if not (2 <= len(normalized) <= 3):
      raise ValueError(f"Unable to parse orientation: {normalized}")

    POSITIVE = 0
    NEGATIVE = 1

    signs = [ POSITIVE, POSITIVE, POSITIVE ]
    mapping = { "X": 0, "Y": 1, "Z": 2 }

    for i in range(len(orientation) - 1):
      if orientation[i] == "-":
        signs[mapping[orientation[i+1]]] = NEGATIVE

    permutation = TO_AXIS_PERMUTATION[normalized]

    return CoordinateFrame(*signs, permutation)

  def __eq__(self, other) -> bool:
    if isinstance(other, str):
      return str(self) == other.upper()

    return (
      self.sign_x == other.sign_x and
      self.sign_y == other.sign_y and
      self.sign_z == other.sign_z and
      self.permutation == other.permutation
    )


class AttributeType(IntEnum):
  VERTEX = 0
  EDGE = 1

TO_AXIS_PERMUTATION = {
  'XYZ': AxisPermutationType.XYZ,
  'XZY': AxisPermutationType.XZY,
  'YXZ': AxisPermutationType.YXZ,
  'YZX': AxisPermutationType.YZX,
  'ZXY': AxisPermutationType.ZXY,
  'ZYX': AxisPermutationType.ZYX,
  'XY': AxisPermutationType.XY,
  'YX': AxisPermutationType.YX,
}
FROM_AXIS_PERMUTATION = {
  v:k for k,v in TO_AXIS_PERMUTATION.items()
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

class DimensionlessType(IntEnum):
  UNKNOWN = 0

TO_QUANTITY_TYPE = {
  0: AreaType,
  1: DimensionlessType,
  2: ElectricalType,
  3: EnergyType,
  4: LengthType,
  5: LuminosityType,
  6: MassType,
  7: SubstanceAmount,
  8: TemperatureType,
  9: TimeType,
  10: VolumeType,
}
FROM_QUANTITY_TYPE = { v:k for k,v in TO_QUANTITY_TYPE.items() }

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

@dataclass
class PhysicalUnit:
  prefix:SIPrefixType
  base:IntEnum

  def __eq__(self, other) -> bool:
    return self.prefix == other.prefix and self.base == other.base

  def __str__(self) -> str:
    return f"{self.prefix}{self.base}"


