# Originally from CloudVolume's lib.py

from typing import Union

import re

import numpy as np

def map2(fn, a, b):
  assert len(a) == len(b), "Vector lengths do not match: {} (len {}), {} (len {})".format(a[:3], len(a), b[:3], len(b))

  result = np.empty(len(a))

  for i in range(len(result)):
    result[i] = fn(a[i], b[i])

  if isinstance(a, Vec) or isinstance(b, Vec):
    return Vec(*result)

  return result

def max2(a, b):
  return map2(max, a, b).astype(a.dtype)

def min2(a, b):
  return map2(min, a, b).astype(a.dtype)

def clamp(val, low, high):
  return min(max(val, low), high)

def check_bounds(val, low, high, index):
  if val > high or val < low:
    raise OutOfBoundsError(
      f'Value {val} (index={index}) cannot be outside of inclusive range {low} to {high}'
    )
  return val

class Vec(np.ndarray):
    def __new__(cls, *args, **kwargs):
      dtype = kwargs.pop('dtype', None)
      if dtype is None:
        if floating(args):
          dtype = float
        else:
          dtype = int
      
      return super().__new__(
        cls, shape=(len(args),), 
        buffer=np.array(args, dtype=dtype), 
        dtype=dtype
      )

    @classmethod
    def clamp(cls, val, minvec, maxvec):
      x = np.minimum.reduce([
        np.maximum.reduce([val,minvec]), 
        maxvec
      ])
      return Vec(*x)

    def clone(self):
      return Vec(*self[:], dtype=self.dtype)

    def null(self):
        return self.length() <= 10 * np.finfo(np.float32).eps

    def dot(self, vec):
      return sum(self * vec)

    def length2(self):
        return self.dot(self)

    def length(self):
        return math.sqrt(self.dot(self))

    def rectVolume(self):
        return reduce(operator.mul, self)

    def __hash__(self):
      return int(''.join(map(str, self)))

    def __repr__(self):
      values = ",".join([ str(x) for x in self ])
      return f"Vec({values}, dtype={self.dtype})"

def __assign(self, val, index):
  self[index] = val

Vec.x = property(lambda self: self[0], lambda self,val: __assign(self,val,0))
Vec.y = property(lambda self: self[1], lambda self,val: __assign(self,val,1))
Vec.z = property(lambda self: self[2], lambda self,val: __assign(self,val,2))
Vec.w = property(lambda self: self[3], lambda self,val: __assign(self,val,3))

Vec.r = Vec.x
Vec.g = Vec.y
Vec.b = Vec.z
Vec.a = Vec.w

def floating(lst):
  for x in lst:
    if isinstance(x, float):
      return True
  return False

FLT_RE = r'(-?\d+(?:\.\d+)?)' # floating point regexp
FILENAME_RE = re.compile(fr'{FLT_RE}-{FLT_RE}_{FLT_RE}-{FLT_RE}_{FLT_RE}-{FLT_RE}(?:\.gz|\.br|\.zstd)?$')

UNIT_SCALES = {
  "pm": 1e-12,
  "nm": 1e-9,
  "um": 1e-6,
  "mm": 1e-3,
  "m": 1.0,
  "km": 1e3,
  "Mm": 1e6,
  "Gm": 1e9,
  "Tm": 1e12,
}

class Bbox(object):
  __slots__ = [ 'minpt', 'maxpt', '_dtype', 'unit' ]

  def __init__(self, a, b, dtype=None, unit='vx'):
    """
    Represents a three dimensional cuboid in space. 
    Ex: `bbox = Bbox((xmin, ymin, zmin), (xmax, ymax, zmax))`

    Example Units:
      vx: voxels
      nm: nanometers
      um: micrometers
      mm: millimeters
      m: meters
    """
    if dtype is None:
      if floating(a) or floating(b):
        dtype = np.float32
      else:
        dtype = np.int32

    self.minpt = Vec(*[ min(ai,bi) for ai,bi in zip(a,b) ], dtype=dtype)
    self.maxpt = Vec(*[ max(ai,bi) for ai,bi in zip(a,b) ], dtype=dtype)

    self._dtype = np.dtype(dtype)
    self.unit = unit

  def convert_units(
    self, unit, 
    resolution=[1,1,1], resolution_unit="nm"
  ):
    """
    Convert the units of this bounding box either 
    from voxels to physical units (e.g. nanometers) or
    vice-versa, or convert between differing physical 
    scales such as um to nm.

    To convert either way between voxels ("vx") and 
    a physical dimension, you must provide a resolution.
    For voxel to voxel or physical units to physical units,
    the resolution is ignored.

    Allowed Units:
      dimensionless: vx (means voxels)
      physical: pm, nm, um, mm, m, km, Mm, Gm, Tm
    """
    if self.unit == "vx" and unit == "vx":
      return self.clone()
    elif self.unit == "vx" and unit != "vx":
      bbx = self.clone()
      bbx *= np.array(resolution)
      bbx.unit = unit
      return bbx
    elif self.unit != "vx" and unit == "vx":
      bbx = self.clone()
      bbx *= np.round(UNIT_SCALES[resolution_unit] / UNIT_SCALES[self.unit])
      bbx //= np.array(resolution) 
      bbx.unit = unit
      return bbx
    else:
      bbx = self.astype(np.float32)
      scale_factor = UNIT_SCALES[unit] / UNIT_SCALES[self.unit]
      bbx /= scale_factor
      bbx.unit = unit

      if scale_factor > 0:
        return bbx
      else:
        return bbx.astype(self.dtype)

  @property
  def is_physical(self):
    return self.unit != 'vx'

  @classmethod
  def deserialize(cls, bbx_data):
    bbx_data = json.loads(bbx_data)
    return Bbox.from_dict(bbx_data)

  def serialize(self):
    return json.dumps(self.to_dict())

  @property
  def ndim(self):
    return len(self.minpt)

  @property 
  def dtype(self):
    return self._dtype

  @property
  def dx(self):
    return self.maxpt.x - self.minpt.x

  @property
  def dy(self):
    return self.maxpt.y - self.minpt.y

  @property
  def dz(self):
    return self.maxpt.z - self.minpt.z  

  @classmethod
  def intersection(cls, bbx1, bbx2):
    result = Bbox( [ 0 ] * bbx1.ndim, [ 0 ] * bbx2.ndim )

    if not Bbox.intersects(bbx1, bbx2):
      return result
    
    for i in range(result.ndim):
      result.minpt[i] = max(bbx1.minpt[i], bbx2.minpt[i])
      result.maxpt[i] = min(bbx1.maxpt[i], bbx2.maxpt[i])

    return result

  @classmethod
  def intersects(cls, bbx1, bbx2):
    a = bbx1.minpt < bbx2.maxpt
    b = bbx1.maxpt > bbx2.minpt
    c = a & b
    if len(c) == 3:
      return c[0] & c[1] & c[2]
    return np.all(c)

  @classmethod
  def near_edge(cls, bbx1, bbx2, distance=0):
    return (
         np.any( np.abs(bbx1.minpt - bbx2.minpt) <= distance )
      or np.any( np.abs(bbx1.maxpt - bbx2.maxpt) <= distance )
    )

  @classmethod
  def create(cls, obj, context=None, bounded=False, autocrop=False):
    typ = type(obj)
    if typ is Bbox:
      obj = obj
    elif typ in (list, tuple, slice):
      obj = Bbox.from_slices(obj, context, bounded, autocrop)
    elif typ is Vec:
      obj = Bbox.from_vec(obj)
    elif typ in string_types:
      obj = Bbox.from_filename(obj)
    elif typ is dict:
      obj = Bbox.from_dict(obj)
    else:
      raise NotImplementedError("{} is not a Bbox convertible type.".format(typ))

    if context and autocrop:
      obj = Bbox.intersection(obj, context)

    if context and bounded:
      if not context.contains_bbox(obj):
        raise OutOfBoundsError(
          "{} did not fully contain the specified bounding box {}.".format(
            context, obj
        ))

    return obj

  @classmethod
  def from_delta(cls, minpt, plus):
    return Bbox( minpt, Vec(*minpt) + plus )

  @classmethod
  def from_dict(cls, data):
    dtype = data['dtype'] if 'dtype' in data else np.float32
    return Bbox( data['minpt'], data['maxpt'], dtype=dtype)

  @classmethod
  def from_vec(cls, vec, dtype=int):
    return Bbox( (0,0,0), vec, dtype=dtype)

  @classmethod
  def from_filename(cls, filename, dtype=int):
    fname = os.path.basename(filename)
    match = FILENAME_RE.search(fname)

    if match is None:
      raise ValueError(f"Unable to decode bounding box from: {filename}")

    root, ext = os.path.splitext(fname)
    parse_type = float if '.' in root else int

    (xmin, xmax,
     ymin, ymax,
     zmin, zmax) = map(parse_type, match.groups())

    return Bbox( (xmin, ymin, zmin), (xmax, ymax, zmax), dtype=dtype)

  @classmethod
  def from_slices(cls, slices, context=None, bounded=False, autocrop=False):
    if context:
      slices = context.reify_slices(
        slices, bounded=bounded, autocrop=autocrop
      )

    for slc in slices:
      if slc.step not in (None, 1):
        raise ValueError("Non-unitary steps are unsupported. Got: " + str(slc.step))

    return Bbox(
      [ slc.start for slc in slices ],
      [ (slc.start if slc.stop < slc.start else slc.stop) for slc in slices ]
    )

  @classmethod
  def from_list(cls, lst):
    """
    from_list(cls, lst)
    
    the first half of the values are the minpt, 
    the last half are the maxpt
    """
    half = len(lst) // 2 
    return Bbox( lst[:half], lst[half:] )

  @classmethod
  def from_points(cls, arr):
    """Create a Bbox from a point cloud arranged as:
      [
        [x,y,z],
        [x,y,z],
        ...
      ]
    """
    arr = np.array(arr, dtype=np.float32)
    
    mins = np.min(arr, axis=0)
    maxes = np.max(arr, axis=0) + 1

    return Bbox( mins, maxes, dtype=np.int64)

  def to_filename(self, precision=None):
    """
    Renders the Bbox as a string. For example:
    
    >>> Bbox([0,2,4],[1,3,5]).to_filename()
    > '0-1_2-3_4-5'

    If the data is floating point, adding a precision
    allows will round the numbers to that decimal place.
    """
    def render(x):
      if precision:
        return f"{round(x, precision):.{precision}f}"
      return str(x)

    return '_'.join(
      ( render(self.minpt[i]) + '-' + render(self.maxpt[i]) for i in range(self.ndim) )
    )

  def to_slices(self):
    return tuple([
      slice(int(self.minpt[i]), int(self.maxpt[i])) for i in range(self.ndim)
    ])

  def to_list(self):
    return list(self.minpt) + list(self.maxpt)

  def to_dict(self):
    return {
      'minpt': self.minpt.tolist(),
      'maxpt': self.maxpt.tolist(),
      'dtype': np.dtype(self.dtype).name,
    }

  def reify_slices(self, slices, bounded=True, autocrop=False):
    """
    Convert free attributes of a slice object 
    (e.g. None (arr[:]) or Ellipsis (arr[..., 0]))
    into bound variables in the context of this
    bounding box.

    That is, for a ':' slice, slice.start will be set
    to the value of the respective minpt index of 
    this bounding box while slice.stop will be set 
    to the value of the respective maxpt index.

    Example:
      bbx = Bbox( (-1,-2,-3), (1,2,3) )
      bbx.reify_slices( (np._s[:],) )
      
      >>> [ slice(-1,1,1), slice(-2,2,1), slice(-3,3,1) ]

    Returns: [ slice, ... ]
    """
    if isinstance(slices, integer_types) or isinstance(slices, floating_types):
      slices = [ slice(int(slices), int(slices)+1, 1) ]
    elif type(slices) == slice:
      slices = [ slices ]
    elif type(slices) == Bbox:
      slices = slices.to_slices()
    elif slices == Ellipsis:
      slices = []

    slices = list(slices)

    for index, slc in enumerate(slices):
      if slc == Ellipsis:
        fill = self.ndim - len(slices) + 1
        slices = slices[:index] +  (fill * [ slice(None, None, None) ]) + slices[index+1:]
        break

    while len(slices) < self.ndim:
      slices.append( slice(None, None, None) )

    # First three slices are x,y,z, last is channel. 
    # Handle only x,y,z here, channel seperately
    for index, slc in enumerate(slices):
      if isinstance(slc, integer_types) or isinstance(slc, floating_types):
        slices[index] = slice(int(slc), int(slc)+1, 1)
      elif slc == Ellipsis:
        raise ValueError("More than one Ellipsis operator used at once.")
      else:
        start = self.minpt[index] if slc.start is None else slc.start
        end = self.maxpt[index] if slc.stop is None else slc.stop 
        step = 1 if slc.step is None else slc.step

        if step < 0:
          raise ValueError('Negative step sizes are not supported. Got: {}'.format(step))

        if autocrop:
          start = clamp(start, self.minpt[index], self.maxpt[index])
          end = clamp(end, self.minpt[index], self.maxpt[index])
        # note: when unbounded, negative indicies do not refer to
        # the end of the volume as they can describe, e.g. a 1px
        # border on the edge of the beginning of the dataset as in
        # marching cubes.
        elif bounded:
          # if start < 0: # this is support for negative indicies
            # start = self.maxpt[index] + start         
          check_bounds(start, self.minpt[index], self.maxpt[index], index)
          # if end < 0: # this is support for negative indicies
          #   end = self.maxpt[index] + end
          check_bounds(end, self.minpt[index], self.maxpt[index], index)

        slices[index] = slice(start, end, step)

    return slices

  @classmethod
  def expand(cls, *args):
    result = args[0].clone()
    for bbx in args:
      result.minpt = min2(result.minpt, bbx.minpt)
      result.maxpt = max2(result.maxpt, bbx.maxpt)
    return result

  @classmethod
  def clamp(cls, bbx0, bbx1):
    result = bbx0.clone()
    result.minpt = Vec.clamp(bbx0.minpt, bbx1.minpt, bbx1.maxpt)
    result.maxpt = Vec.clamp(bbx0.maxpt, bbx1.minpt, bbx1.maxpt)
    return result

  def size(self):
    return Vec(*(self.maxpt - self.minpt), dtype=self.dtype)

  def size3(self):
    return Vec(*(self.maxpt[:3] - self.minpt[:3]), dtype=self.dtype)

  def subvoxel(self):
    """
    Previously, we used bbox.volume() < 1 for testing
    if a bounding box was larger than one voxel. However, 
    if two out of three size dimensions are negative, the 
    product will be positive. Therefore, we first test that 
    the maxpt is to the right of the minpt before computing 
    whether conjunctioned with volume() < 1.

    Returns: boolean
    """
    return (not self.valid()) or self.volume() < 1

  def empty(self):
    """
    Previously, we used bbox.volume() <= 0 for testing
    if a bounding box was empty. However, if two out of 
    three size dimensions are negative, the product will 
    be positive. Therefore, we first test that the maxpt 
    is to the right of the minpt before computing whether 
    the bbox is empty and account for 20x machine epsilon 
    of floating point error.

    Returns: boolean
    """
    return (not self.valid()) or (self.volume() < (20 * MACHINE_EPSILON))

  def valid(self):
    return np.all(self.minpt <= self.maxpt)

  def volume(self):
    if np.issubdtype(self.dtype, np.integer):
      return self.size3().astype(np.int64).rectVolume()
    else:
      return self.size3().astype(np.float64).rectVolume()

  def center(self):
    return (self.minpt + self.maxpt) / 2.0

  def adjust(self, amt: Union[int, tuple, list]):
    if isinstance(amt, tuple) or isinstance(amt, list):
      assert len(amt) == 3
      amt = np.asarray(amt)

    self.minpt -= amt
    self.maxpt += amt
    
    if not self.valid():
      raise ValueError("Cannot shrink bbox below zero volume.")

    return self

  def shrink(self, amt: Union[int, tuple, list]):
    if isinstance(amt, int):
      assert amt > 0
    elif isinstance(amt, tuple) or isinstance(amt, list):
      amt = np.array(amt)
      assert np.all(amt > 0) # type: ignore

    # make it negative for shrink
    return self.adjust(-amt) # type: ignore

  def grow(self, amt: Union[int, tuple, list]):
    return self.adjust(amt)

  def expand_to_chunk_size(self, chunk_size, offset=Vec(0,0,0, dtype=int)):
    """
    Align a potentially non-axis aligned bbox to the grid by growing it
    to the nearest grid lines.

    Required:
      chunk_size: arraylike (x,y,z), the size of chunks in the 
                    dataset e.g. (64,64,64)
      offset: arraylike (x,y,z) the origin of the coordinate system
        so that this offset can be accounted for in the grid line 
        calculation.
    Optional:
      offset: arraylike (x,y,z), the starting coordinate of the dataset
    """
    chunk_size = np.array(chunk_size, dtype=np.float32)
    result = self.clone()
    result -= offset
    result.minpt = np.floor(result.minpt / chunk_size) * chunk_size
    result.maxpt = np.ceil(result.maxpt / chunk_size) * chunk_size 
    return (result + offset).astype(self.dtype)

  def shrink_to_chunk_size(self, chunk_size, offset=Vec(0,0,0, dtype=int)):
    """
    Align a potentially non-axis aligned bbox to the grid by shrinking it
    to the nearest grid lines.

    Required:
      chunk_size: arraylike (x,y,z), the size of chunks in the 
                    dataset e.g. (64,64,64)
      offset: arraylike (x,y,z) the origin of the coordinate system
        so that this offset can be accounted for in the grid line 
        calculation.
    Optional:
      offset: arraylike (x,y,z), the starting coordinate of the dataset
    """
    chunk_size = np.array(chunk_size, dtype=np.float32)
    result = self.clone()
    result = result - offset
    result.minpt = np.ceil(result.minpt / chunk_size) * chunk_size
    result.maxpt = np.floor(result.maxpt / chunk_size) * chunk_size 

    # If we are inside a single chunk, the ends
    # can invert, which tells us we should collapse
    # to a single point.
    if np.any(result.minpt > result.maxpt):
      result.maxpt = result.minpt.clone()

    return (result + offset).astype(self.dtype)

  def round_to_chunk_size(self, chunk_size, offset=Vec(0,0,0, dtype=int)):
    """
    Align a potentially non-axis aligned bbox to the grid by rounding it
    to the nearest grid lines.

    Required:
      chunk_size: arraylike (x,y,z), the size of chunks in the 
                    dataset e.g. (64,64,64)
      offset: arraylike (x,y,z) the origin of the coordinate system
        so that this offset can be accounted for in the grid line 
        calculation.
    Optional:
      offset: arraylike (x,y,z), the starting coordinate of the dataset
    """
    chunk_size = np.array(chunk_size, dtype=np.float32)
    result = self.clone()
    result = result - offset
    result.minpt = np.round(result.minpt / chunk_size) * chunk_size
    result.maxpt = np.round(result.maxpt / chunk_size) * chunk_size
    return (result + offset).astype(self.dtype)

  def num_chunks(self, chunk_size):
    """Computes the number of chunks inside this bbox for a given chunk size."""
    Nfn = lambda i: math.ceil((self.maxpt[i] - self.minpt[i]) / chunk_size[i])
    return reduce(operator.mul, map(Nfn, range(len(self.minpt))))

  def contains(self, point):
    """
    Tests if a point on or within a bounding box.

    Returns: boolean
    """
    for x in (point < self.minpt):
      if x:
        return False
    for x in (point > self.maxpt):
      if x:
        return False
    return True

  def contains_bbox(self, bbox):
    return (
      self.contains(bbox.maxpt)
      and self.contains(bbox.minpt) 
    )

  def overlaps_bbox(self, bbox):
    return not (
      np.any(self.maxpt < bbox.minpt)
      or np.any(self.minpt > bbox.maxpt)
    )

  def clone(self):
    return Bbox(self.minpt, self.maxpt, dtype=self.dtype, unit=self.unit)

  def astype(self, typ):
    tmp = self.clone()
    tmp.minpt = tmp.minpt.astype(typ)
    tmp.maxpt = tmp.maxpt.astype(typ)
    tmp._dtype = tmp.minpt.dtype 
    return tmp

  def transpose(self):
    return Bbox(self.minpt[::-1], self.maxpt[::-1])

  # note that operand can be a vector 
  # or a scalar thanks to numpy
  def __isub__(self, operand): 
    if isinstance(operand, Bbox):
      self.minpt = np.subtract(self.minpt, operand.minpt, casting="safe")
      self.maxpt = np.subtract(self.maxpt, operand.maxpt, casting="safe")
    else:
      self.minpt = np.subtract(self.minpt, operand, casting="safe")
      self.maxpt = np.subtract(self.maxpt, operand, casting="safe")

    return self.astype(self.minpt.dtype)

  def __sub__(self, operand):
    tmp = self.clone()
    return tmp.__isub__(operand)

  def __iadd__(self, operand):
    if isinstance(operand, Bbox):
      self.minpt = np.add(self.minpt, operand.minpt, casting="safe")
      self.maxpt = np.add(self.maxpt, operand.maxpt, casting="safe")
    else:
      self.minpt = np.add(self.minpt, operand, casting="safe")
      self.maxpt = np.add(self.maxpt, operand, casting="safe")

    return self

  def __add__(self, operand):
    tmp = self.clone()
    return tmp.__iadd__(operand)

  def __imul__(self, operand):
    self.minpt = np.multiply(self.minpt, operand, casting="safe")
    self.maxpt = np.multiply(self.maxpt, operand, casting="safe")
    self._dtype = self.minpt.dtype 
    return self

  def __mul__(self, operand):
    tmp = self.clone()
    tmp.minpt = np.multiply(tmp.minpt, operand, casting="safe")
    tmp.maxpt = np.multiply(tmp.maxpt, operand, casting="safe")
    return tmp.astype(tmp.minpt.dtype)

  def __idiv__(self, operand):
    if (
      isinstance(operand, float) \
      or self.dtype in (float, np.float32, np.float64) \
      or (hasattr(operand, 'dtype') and operand.dtype in (float, np.float32, np.float64))
    ):
      return self.__itruediv__(operand)
    else:
      return self.__ifloordiv__(operand)

  def __div__(self, operand):
    if (
      isinstance(operand, float) \
      or self.dtype in (float, np.float32, np.float64) \
      or (hasattr(operand, 'dtype') and operand.dtype in (float, np.float32, np.float64))
    ):

      return self.__truediv__(operand)
    else:
      return self.__floordiv__(operand)

  def __ifloordiv__(self, operand):
    self.minpt //= operand
    self.maxpt //= operand
    return self

  def __floordiv__(self, operand):
    tmp = self.astype(float)
    tmp.minpt //= operand
    tmp.maxpt //= operand
    return tmp.astype(int)

  def __itruediv__(self, operand):
    res = self.__truediv__(operand)
    self.minpt[:] = res.minpt[:]
    self.maxpt[:] = res.maxpt[:]
    return self

  def __truediv__(self, operand):
    tmp = self.clone()

    if isinstance(operand, int):
      operand = float(operand)

    tmp.minpt = Vec(*( tmp.minpt.astype(float) / operand ), dtype=float)
    tmp.maxpt = Vec(*( tmp.maxpt.astype(float) / operand ), dtype=float)
    return tmp.astype(tmp.minpt.dtype)

  def __ne__(self, other):
    return not (self == other)

  def __eq__(self, other):
    return np.array_equal(self.minpt, other.minpt) and np.array_equal(self.maxpt, other.maxpt) and self.unit == other.unit

  def __hash__(self):
    return int(''.join(map(str, map(int, self.to_list()))))

  def __repr__(self):
    normfn = int
    if np.issubdtype(self.dtype, np.floating):
      normfn = float
    return f"Bbox({[ normfn(x) for x in self.minpt ]},{[ normfn(x) for x in self.maxpt ]}, dtype=np.{self.dtype}, unit='{self.unit}')"

# From SO: https://stackoverflow.com/questions/14313510/how-to-calculate-rolling-moving-average-using-python-numpy-scipy
def moving_average(a:np.ndarray, n:int, mode:str = "symmetric") -> np.ndarray:
  if n <= 0:
    raise ValueError(f"Window size ({n}), must be >= 1.")
  elif n == 1:
    return a

  if len(a) == 0:
    return a

  if a.ndim == 2:
    a = np.pad(a, [[n, n],[0,0]], mode=mode)
  else:
    a = np.pad(a, [n, n], mode=mode)

  ret = np.cumsum(a, dtype=float, axis=0)
  ret = (ret[n:] - ret[:-n])[:-n]
  ret /= float(n)
  return ret