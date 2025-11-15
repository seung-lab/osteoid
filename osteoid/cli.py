import os

import click

@click.group()
def main():
  """
  Tools for viewing, converting, and reading
  skeleton formats such as swc, ostd.
  """
  pass

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

