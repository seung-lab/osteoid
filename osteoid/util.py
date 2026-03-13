from typing import Union, IO

import mmap
import os
import gzip
import lzma

from .skeleton import Skeleton

def _load(filelike, size:int = -1, allow_mmap:bool = False) -> IO[bytes]:
  if hasattr(filelike, 'read'):
    binary = filelike.read(size)
  elif (
    isinstance(filelike, str) 
    and os.path.splitext(filelike)[1] == '.gz'
  ):
    with gzip.open(filelike, 'rb') as f:
      binary = f.read(size)
  elif (
    isinstance(filelike, str) 
    and os.path.splitext(filelike)[1] in ('.lzma', '.xz')
  ):
    with lzma.open(filelike, 'rb') as f:
      binary = f.read(size)
  else:
    with open(filelike, 'rb') as f:
      if allow_mmap:
        binary = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
      else:
        binary = f.read(size)
  
  return binary

def load(filename:str, allow_mmap:bool = False) -> Skeleton:
  binary = _load(filename, allow_mmap=allow_mmap)

  if filename.endswith("swc"):
    return Skeleton.from_swc(bytes(binary).decode("utf8"))
  else:
    return Skeleton.from_precomputed(binary)

def save(
  filename:str,
  skeleton:Skeleton,
  **kwargs
):
  """Save labels into the file-like object or file path."""
  if filename.endswith("swc"):
    binary = skeleton.to_swc().encode("utf8")
  else:
    binary = skeleton.to_precomputed()

  if (
    isinstance(filename, str) 
    and os.path.splitext(filename)[1] == '.gz'
  ):
    with gzip.open(filename, 'wb') as f:
      f.write(binary)
  elif (
    isinstance(filename, str) 
    and os.path.splitext(filename)[1] in ('.lzma', '.xz')
  ):
    with lzma.open(filename, 'wb') as f:
      f.write(binary)
  else:
    with open(filename, 'wb') as f:
      f.write(binary)
