import click
from time import time
from .ftcs import FTCS
from .upwind15 import UpWind15
from .crankn import CrankN
from . import __version__


@click.version_option(__version__)
@click.command()
@click.option('--runame', required=False, default="heat_results", show_default=True,
              type=click.STRING,
              help="name to give run and results dir.")
@click.option('--prec', required=False, default="double", show_default=True,
              type=click.Choice(["half", "float", "double", "quad"]),
              help="precision of the solution convergence.")
@click.option("--alpha", required=False, default=0.2, show_default=True,
              type=click.FLOAT,
              help="material thermal diffusivity (sq-meters/second).")
@click.option("--lenx", required=False, default=1.0, show_default=True,
              type=click.FLOAT,
              help="material length (meters).")
@click.option("--dx", required=False, default=0.1, show_default=True,
              type=click.FLOAT,
              help="x-incriment. Best if lenx/dx==int. (meters).")
@click.option("--dt", required=False, default=0.004, show_default=True,
              type=click.FLOAT,
              help="t-incriment (seconds).")
@click.option("--maxt", required=False, default=2.0, show_default=True,
              type=click.FLOAT,
              help=">0:max sim time (seconds) | <0:min l2 change in soln.")
@click.option("--bc0", required=False, default=0, show_default=True,
              type=click.FLOAT,
              help="boundary condition @ x=0: u(0,t) (Kelvin)")
@click.option("--bc1", required=False, default=1, show_default=True,
              type=click.FLOAT,
              help="boundary condition @ x=lenx: u(lenx,t) (Kelvin)")
@click.option('--ic', required=False, default="const(1)", show_default=True,
              type=click.STRING,
              help="initial condition @ t=0: u(x,0) (Kelvin)")
@click.option('--alg', required=False, default="ftcs", show_default=True,
              type=click.Choice(["ftcs", "upwind15", "crankn"]),
              help="algorithm")
@click.option("--savi", required=False, default=0, show_default=True,
              type=click.INT,
              help="save every i-th solution step")
@click.option("--save", required=False, default=0, show_default=True,
              type=click.INT,
              help="save error in every saved solution")
@click.option("--outi", required=False, default=100, show_default=True,
              type=click.INT,
              help="output progress every i-th solution step")
@click.option("--noout", required=False, default=0, show_default=True,
              type=click.INT,
              help="disable all file outputs")
def main(runame: str, prec: str, alpha: float, lenx: float,
         dx: float, dt: float, maxt: float, bc0: float,
         bc1: float, ic: str, alg: str, savi: int,
         save: int, outi: int, noout: int) -> None:
    """Main entry point for heateq_design."""
    click.echo('Invoking heat equation solver...')
    t0 = time()
    if alg == 'ftcs':
        heat_solver = FTCS(lenx, maxt, alpha, dx, dt, bc0, bc1, ic, outi)
    elif alg == 'upwind15':
        heat_solver = UpWind15(lenx, maxt, alpha, dx, dt, bc0, bc1, ic, outi)
    else:
        heat_solver = CrankN(lenx, maxt, alpha, dx, dt, bc0, bc1, ic, outi)
    heat_solver.solve()
    t1 = time() - t0
    click.echo('Solver complete. Results generated here:' + runame)
    click.echo("Time elapsed: " + str(t1))
