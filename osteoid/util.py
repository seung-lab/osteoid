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
  binary = _load(filename, allow_mmap=mmap)

  if filename.endswith("swc"):
    return Skeleton.from_swc(binary.decode("utf8"))
  elif filename.endswith("ostd"):
    return Skeleton.from_ostd(binary)
  else:
    return Skeleton.from_precomputed(binary)

def save(
  filelike,
  skeleton:Skeleton,
  **kwargs
):
  """Save labels into the file-like object or file path."""
  if filename.endswith("swc"):
    binary = skel.to_swc()
  elif filename.endswith("ostd"):
    binary = skel.to_ostd()
  else:
    binary = skel.to_precomputed()

  if hasattr(filelike, 'write'):
    filelike.write(binary)
  elif (
    isinstance(filelike, str) 
    and os.path.splitext(filelike)[1] == '.gz'
  ):
    with gzip.open(filelike, 'wb') as f:
      f.write(binary)
  elif (
    isinstance(filelike, str) 
    and os.path.splitext(filelike)[1] in ('.lzma', '.xz')
  ):
    with lzma.open(filelike, 'wb') as f:
      f.write(binary)
  else:
    with open(filelike, 'wb') as f:
      f.write(binary)
