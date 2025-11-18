import os
from pathlib import Path
import time

import click
import numpy as np

@click.group()
def main():
  """
  Tools for viewing, converting, and reading
  skeleton formats such as swc, ostd.
  """
  pass


@main.command()
@click.argument("src")
@click.option('-u', '--uniq', is_flag=True, default=False, help="Print unique ids in sorted order.", show_default=True)
def ids(src, uniq):
  """Prints all unique ids in the potentially multi-part ostd file."""
  import osteoid.formats
  from osteoid.formats.ostd import OstdSkeleton, OstdHeader

  oskel = OstdSkeleton.load(src, allow_mmap=True)

  ids = oskel.ids()
  if uniq:
    ids = np.unique(ids)

  ids = [ str(x) for x in ids ]
  click.echo(",".join(ids))

@main.command()
@click.argument("src")
def info(src):
  import osteoid.formats
  from osteoid.formats.ostd import OstdSkeleton, OstdHeader

  if src.endswith(".swc"):
    click.echo(osteoid.formats.swc.read_header(src))
  elif src.endswith(".ostd"):
    with open(src, "rb") as f:
      binary = f.read(OstdHeader.HEADER_BYTES)
    header = OstdHeader.from_bytes(binary, skip_total_length_check=True)
    click.echo(header.details())

    oskel = OstdSkeleton.load(src, allow_mmap=True)
    attrs = [ f"{attr} ({unit})" for attr, unit in oskel.attributes ]
    click.echo("\nattributes: " + ", ".join(attrs))

@main.command()
@click.argument("src")
@click.argument("dest")
def convert(src:str, dest:str):
  """Convert a skeleton from one format to another."""
  from .util import load, save
  skel = load(src, allow_mmap=True)
  save(dest, skel)

@main.command()
@click.option('-v', '--verbose', is_flag=True, default=False, help="Print information about conversion progress.", show_default=True)
@click.option('-d', '--dir', default=".", help="Where to write the new skeletons.", show_default=True)
@click.option('-f', '--force', is_flag=True, default=False, help="Overwrite files without asking.", show_default=True)
@click.option('-t', '--type', required=True, type=click.Choice(['swc', 'ostd']), help="Convert all file specified to this type. If the type is the same leave the file alone.", show_default=True)
@click.argument("src", nargs=-1)
def convertall(src, type, verbose, force, dir):
  """Convert many skeletons from one format to another."""
  from .util import load, save

  srcs = []
  for path in src:
    if path[-1] == "*":
      path = path[:-1]
      path = [ os.path.join(path, x) for x in os.listdir(path) ]
      srcs.extend(path)
    else:
      srcs.append(path)

  new_dir = Path(dir)

  for path_str in srcs:
      path = Path(path_str)
      if not path.exists() and verbose:
        click.echo(f"File {path} does not exist, skipping.")
        continue
      
      dest_path = new_dir / path.with_suffix(f".{type}").name

      if dest_path.exists() and not force:
        click.echo(f"{path.name} exists, aborting.")
        return

      if path.suffix.lower() == f".{type}" and verbose:
          click.echo(f"{path.name} is already {type}, skipping.")
          continue
      
      s = time.perf_counter()
      skel = load(str(path), allow_mmap=True)
      save(str(dest_path), skel)      
      t = time.perf_counter() - s
      if verbose:
        click.echo(f"Converted {path.name} to {dest_path.name} in {t:.2f} sec")

@main.command()
@click.argument("filename")
@click.option('-c', '--color-by', type=click.Choice(['r', 'c', 'x']), default='r', help="For skeleton visualization. r = radius, c = components, x = cross sectional area (if available).", show_default=True)
def view(filename, color_by):
  """Visualize a .swc or .npy file."""
  import microviewer
  from .util import load

  skel = load(filename, allow_mmap=True)

  error_text = {
    'r': "radius",
    'x': "cross sectional area",
  }

  try:
    microviewer.objects([ skel ], skeleton_color_by=color_by)
  except AttributeError:
    click.echo(f"ostd: skeleton does not have a {error_text[color_by]} attribute.")
   
@main.command()
def license():
  """Prints the license for this library and cli tool."""
  path = os.path.join(os.path.dirname(__file__), 'LICENSE')
  with open(path, 'rt') as f:
    print(f.read())

