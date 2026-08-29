from dataclasses import dataclass
from enum import IntEnum
import math
import struct

import numpy as np

class AttributeType(IntEnum):
  VERTEX = 0
  EDGE = 1

class CompressionType(IntEnum):
  NONE = 0
  GZIP = 1
  BZIP2 = 2
  ZSTD = 3
  DRACO = 4

class DataType(IntEnum):
  F16 = 0
  F32 = 1
  F64 = 2
  U8 = 3
  U16 = 4
  U32 = 5
  U64 = 6
  I8 = 7
  I16 = 8
  I32 = 9
  I64 = 10
  BOOL = 11
  PACKED_BOOL = 12

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

class GraphType(IntEnum):
  GRAPH = 0
  TREE = 1
  CYCLIC = 2

class LogEnum(IntEnum):
  LINEAR = 0
  LOG10 = 1
  LOG2 = 2
  LN = 3

SI_PREFIX = {
  -30: 'q', # quecto
  -27: 'r', # ronto
  -24: 'y', # yocto
  -21: 'z', # zepto
  -18: 'a', # atto
  -15: 'f', # femto
  -12: 'p', # pico
  -9: 'n', # nano
  -6: 'u', # micro
  -3: 'm', # milli
  -2: 'c', # centi
  0: '',
  3: 'k', # kilo
  6: 'M', # mega
  9: 'G', # giga
  12: 'T', # tera
  15: 'P', # peta
  18: 'E', # exa
  21: 'Z', # zetta
  24: 'Y', # yotta
  27: 'R', # ronna
  30: 'Q', # quetta
}

@dataclass
class SIUnit:
  si_prefix:int = 0
  amperes:int = 0
  kelvin:int = 0
  kilograms:int = 0
  meters:int = 0
  moles:int = 0
  seconds:int = 0
  scale:LogEnum = LogEnum.LINEAR

  def __str__(self) -> str:
    dimensions = [
      (self.amperes, f"A"),
      (self.kelvin, f"K"),
      (self.kilograms, f"kg"),
      (self.meters, f"m"),
      (self.moles, f"mol"),
      (self.seconds, f"s"),
    ]

    rendered = ""
    for exp, symbol in dimensions:
      if exp == 0:
        continue
      elif exp == 1:
        rendered += f" * {symbol}"
      elif exp == -1:
        rendered += f" / {symbol}"
      elif exp < 0:
        rendered += f" / {symbol}^{abs(exp)}"
      else:
        rendered += f" * {symbol}^{abs(exp)}"

    if self.scale == LogEnum.LOG10:
      rendered += " (log10)"
    elif self.scale == LogEnum.LOG2:
      rendered += " (log2)"
    elif self.scale == LogEnum.LN:
      rendered += " (ln)"

    if len(rendered) > 0:
      rendered = rendered[3:] # strip first * or /

    prefix = SI_PREFIX.get(self.si_prefix, '?')

    if rendered == '':
      rendered = 'dimensionless'

    return f"{prefix}{rendered}"

  @classmethod
  def from_bytes(kls, binary:bytes) -> "SIUnit":
    return kls(
      si_prefix=int.from_bytes(binary[0:1], 'little', signed=True),
      amperes=int.from_bytes(binary[1:2], 'little', signed=True),
      kelvin=int.from_bytes(binary[2:3], 'little', signed=True),
      kilograms=int.from_bytes(binary[3:4], 'little', signed=True),
      meters=int.from_bytes(binary[4:5], 'little', signed=True),
      moles=int.from_bytes(binary[5:6], 'little', signed=True),
      seconds=int.from_bytes(binary[6:7], 'little', signed=True),
      scale=LogEnum(int.from_bytes(binary[7:8], 'little', signed=False)),
    )

  def to_bytes(self) -> bytes:
    return b''.join([
      int(self.si_prefix).to_bytes(1, 'little', signed=True),
      int(self.amperes).to_bytes(1, 'little', signed=True),
      int(self.kelvin).to_bytes(1, 'little', signed=True),
      int(self.kilograms).to_bytes(1, 'little', signed=True),
      int(self.meters).to_bytes(1, 'little', signed=True),
      int(self.moles).to_bytes(1, 'little', signed=True),
      int(self.seconds).to_bytes(1, 'little', signed=True),
      int(self.scale).to_bytes(1, 'little'),
    ])

class CompressionType(IntEnum):
  NONE = 0
  GZIP = 1
  BZIP2 = 2
  ZSTD = 3
  DRACO = 4

class GraphType(IntEnum):
  GRAPH = 0
  TREE = 1
  CYCLIC = 2

class SpaceType(IntEnum):
  UNKNOWN = 0
  VOXEL = 1
  PHYSICAL = 2
  SCANNER = 3
  ATLAS = 4
  WORLD = 5
  ANCHOR = 6

SPACE_SYMBOLS = {
  SpaceType.UNKNOWN: "UNKNOWN",
  SpaceType.VOXEL: "VOXEL",
  SpaceType.PHYSICAL: "PHYSICAL",
  SpaceType.SCANNER: "SCANNER",
  SpaceType.ATLAS: "ATLAS",
  SpaceType.WORLD: "WORLD",
  SpaceType.WORLD: "ANCHOR",
}

SpaceType.__str__ = lambda self: SPACE_SYMBOLS[self]

@dataclass
class CoordinateFrame:
  signs:list[bool]
  lehmer_code:int
  num_space_like:int
  voxel_centered:bool

  @classmethod
  def from_bytes(kls, binary:bytes) -> "CoordinateFrame":
    if len(binary) != 4:
        raise ValueError(f"Expected 4 bytes, got {len(binary)}")

    # B = unsigned char
    header, sign_bits, lehmer_lo, lehmer_hi = struct.unpack('BBBB', binary)
    lehmer_code = lehmer_lo | (lehmer_hi << 8)

    n_axes = (header & 0b111) + 1
    num_space_like = (header >> 3) & 0b111
    voxel_centered = bool((header >> 6) & 1)
    # reserved bit (header >> 7) is ignored

    signs = []
    for i in range(n_axes):
        signs.append(bool((sign_bits >> i) & 1))

    return kls(
      signs=signs,
      lehmer_code=lehmer_code,
      num_space_like=num_space_like,
      voxel_centered=voxel_centered,
    )

  def to_bytes(self) -> bytes:
    n_axes = len(self.signs)
    if n_axes > 8:
        raise ValueError("Supports up to 8 axes")

    header = ((n_axes - 1) & 0b111)          # 3 bits for axes
    header |= ((self.num_space_like & 0b111) << 3)  # next 3 bits for space-like
    header |= ((int(self.voxel_centered) & 1) << 6) # bit 6: voxel_centered
    # bit 7 reserved as 0

    sign_bits = 0
    for i, s in enumerate(self.signs):
        if s:
            sign_bits |= (1 << i)

    # pack Lehmer code as 2 bytes (assuming <= 8 axes, max 8! = 40320 fits in 2 bytes)
    lehmer_bytes = struct.pack('<H', self.lehmer_code) # H = unsigned short

    # B = unsigned char
    return struct.pack('BB', header, sign_bits) + lehmer_bytes

  @classmethod
  def rank_permutation(kls, permutation:list[int]) -> int:
    """Convert an axis permutation to an index."""
    axes = list(range(len(permutation)))
    n = len(axes)
    rank = 0

    for i in range(n):
        idx = axes.index(permutation[i])
        rank += idx * math.factorial(n - i - 1)
        axes.pop(idx)

    return rank

  @classmethod
  def unrank_permutation(kls, num_axes:int, k:int) -> list[int]:
    """Convert an index into an axis permutation for a given number of axes."""
    axes = list(range(num_axes))
    perm = []

    for i in range(num_axes, 0, -1):
        f = math.factorial(i - 1)
        idx = k // f
        k %= f
        perm.append(axes.pop(idx))

    return perm

  def __str__(self) -> str:
    def sgn(x):
      return "+" if not x else "-"

    signs = [ sgn(x) for x in self.signs ]
    axes = self.unrank_permutation(len(signs), self.lehmer_code)

    if self.num_space_like == 1:
      convention = [ 'X', 'T', '2', '3', '4', '5', '6', '7' ]
    elif self.num_space_like == 2:
      convention = [ 'X', 'Y', 'T', '3', '4', '5', '6', '7' ]  
    else:
      convention = [ 'X', 'Y', 'Z', 'T', '4', '5', '6', '7' ]

    out = ""
    for i in range(len(axes)):
      out += f"{signs[axes[i]]}{convention[axes[i]]}"
    return out

  @classmethod
  def parse(kls, orientation:str, voxel_centered:bool) -> "CoordinateFrame":
    if len(orientation) > 6:
      raise ValueError(f"Unable to parse orientation: {orientation[:100]}")

    orientation = orientation.upper()
    normalized = orientation.replace('+', '').replace('-', '')

    if not (1 <= len(normalized) <= 8):
      raise ValueError(f"Unable to parse orientation: {normalized}")

    POSITIVE = 0
    NEGATIVE = 1

    signs = [ POSITIVE, POSITIVE, POSITIVE, POSITIVE ]
    mapping = { 
      "X": 0, "Y": 1, "Z": 2, "T": 3, 
      "0": 0, "1": 1, "2": 2, "3": 3, "4": 4, "5": 5, "6": 6, "7": 7 
    }

    space_like = {'X', 'Y', 'Z'}
    num_space_like = 0

    for i in range(len(orientation) - 1):
      if orientation[i] == "-":
        signs[mapping[orientation[i+1]]] = NEGATIVE

    code = kls.rank_permutation([ mapping[axis] for axis in normalized ])

    for axis in normalized:
      num_space_like += int(axis in space_like)

    return kls(
      signs=signs, 
      lehmer_code=code, 
      num_space_like=num_space_like,
      voxel_centered=voxel_centered,
    )

  def __eq__(self, other) -> bool:
    if isinstance(other, str):
      return str(self) == other.upper()

    return (
      self.signs == other.signs and
      self.lehmer_code == other.lehmer_code
    )

TRANSLATE_UNIT = {
  # LENGTH
  "vx": SIUnit(),
  "voxel": SIUnit(),
  
  "fm": SIUnit(si_prefix=-15, meters=1),
  "femtometer": SIUnit(si_prefix=-5, meters=1),
  
  "pm": SIUnit(si_prefix=-12, meters=1),
  "picometer": SIUnit(si_prefix=-4, meters=1),
  
  "nm": SIUnit(si_prefix=-9, meters=1),
  "nanometer": SIUnit(si_prefix=-9, meters=1),
  
  "um": SIUnit(si_prefix=-6, meters=1),
  "micrometer": SIUnit(si_prefix=-6, meters=1),
  "micron": SIUnit(si_prefix=-6, meters=1),
  
  "mm": SIUnit(si_prefix=-3, meters=1),
  "millimeter": SIUnit(si_prefix=-3, meters=1),

  "cm": SIUnit(si_prefix=-2, meters=1),
  "centimeter": SIUnit(si_prefix=-2, meters=1),

  "m": SIUnit(si_prefix=0, meters=1),
  "meter": SIUnit(si_prefix=0, meters=1),
  
  "km": SIUnit(si_prefix=3, meters=1),
  "kilometer": SIUnit(si_prefix=3, meters=1),

  "Mm": SIUnit(si_prefix=6, meters=1),

  # Area
  "fm^2": SIUnit(si_prefix=-15, meters=2),
  "pm^2": SIUnit(si_prefix=-12, meters=2),
  "nm^2": SIUnit(si_prefix=-9, meters=2),  
  "um^2": SIUnit(si_prefix=-6, meters=2),  
  "mm^2": SIUnit(si_prefix=-3, meters=2),
  "cm^2": SIUnit(si_prefix=-2, meters=2),
  "m^2": SIUnit(si_prefix=0, meters=2),  
  "km^2": SIUnit(si_prefix=3, meters=2),
  "Mm^2": SIUnit(si_prefix=6, meters=2),

  # VOLUME
  "cc": SIUnit(si_prefix=-6, meters=3), # cubic centimeter
  "ccm": SIUnit(si_prefix=-6, meters=3), # cubic centimeter
  "mL": SIUnit(si_prefix=-6, meters=3),
  "milliliters": SIUnit(si_prefix=-6, meters=3),
  
  "L": SIUnit(si_prefix=-3, meters=3),
  "liters": SIUnit(si_prefix=-3, meters=3),

  # MASS
  "ug": SIUnit(si_prefix=-9, kilograms=1),
  "micrograms": SIUnit(si_prefix=-9, kilograms=1),

  "mg": SIUnit(si_prefix=-6, kilograms=1),
  "milligrams": SIUnit(si_prefix=-6, kilograms=1),
  
  "g": SIUnit(si_prefix=-3, kilograms=1),
  "grams": SIUnit(si_prefix=-3, kilograms=1),
  
  "kg": SIUnit(si_prefix=0, kilograms=1),
  "kilograms": SIUnit(si_prefix=0, kilograms=1),

  # FORCE
  "N": SIUnit(kilograms=1, meters=1, seconds=-2), # newtons

  # ENERGY
  "J": SIUnit(kilograms=1, meters=2, seconds=-2), # joules

  "mW": SIUnit(si_prefix=-3, kilograms=1, meters=2, seconds=-3), 
  "W": SIUnit(kilograms=1, meters=2, seconds=-3), # watts
}

