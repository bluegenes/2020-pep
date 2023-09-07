import os
import numpy as np
import seaborn as sns
from sklearn.manifold import TSNE

from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

compare_matrix = snakemake.input
outplot = snakemake.output



#TSNE (
#n_components: int=2,
#perplexity: float=30,
#early_exaggeration: float=12,
#learning_rate: float=200,
#n_iter: int=1000,
#n_iter_without_progress: int=300,
#min_grad_norm: float=0.0000001,
#metric: str=builtins.str,
#init: =random,
#verbose: int=0,
#random_state: NoneType=None,
#method: str=builtins.str,
#angle: float=0.5
#)
