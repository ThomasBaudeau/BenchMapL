from distutils.command.config import config
import pysam
import numpy as np
import matplotlib.pyplot as plt
from BenchPlot import *
import pandas as pd

files_stats(snakemake.input['a'],snakemake.input['b'], snakemake.output, snakemake.config)