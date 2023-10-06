from distutils.command.config import config
import pysam
import numpy as np
import matplotlib.pyplot as plt
from BenchPlot import *

plot_params(snakemake.input['a'], snakemake.output, snakemake.config)