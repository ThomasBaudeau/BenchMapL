from distutils.command.config import config
import pysam
import numpy as np
import matplotlib.pyplot as plt
from BenchPlot import *



plot_simple(snakemake.input, snakemake.output, snakemake.config)
