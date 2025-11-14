import os

import click

@click.group()
def main():
  """
  View and convert skeleton formats.
  """
  pass

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

