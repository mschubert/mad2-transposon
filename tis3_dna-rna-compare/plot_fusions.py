# requirements:
# https://github.com/jrderuiter/imfusion-analyses
# https://github.com/jrderuiter/geneviz.git
# https://github.com/jrderuiter/pybiomart.git
# https://github.com/jrderuiter/ngs-tk.git (+ link ngs_tk to ngstk) //ERR: ngstk.tabix, ngstk.util.genomic
import sys
sys.path.append('imfusion-analyses/src')
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from pathlib import Path
from imfusion.model import Insertion
from nbsupport.insertions import plot_insertion
sns.set_style('white')

#def plot_insertion(key):
#    insertion = insertion_lookup[key]


de_ctgs = (pd.read_csv("../data/rnaseq_imfusion/merged_ctgs.txt", sep='\t')
           .query('de_pvalue < 0.05'))
de_ctgs.head()

insertions = Insertion.read_csv("../data/rnaseq_imfusion/insertions.txt", sep='\t')
insertion_lookup = {ins.id: ins for ins in Insertion.from_frame(insertions)}

base_dir = Path("../data/rnaseq_imfusion")
junction_paths = {fp.parent.parent.name: fp 
                  for fp in base_dir.glob('**/SJ.out.tab')
                  if '_STARpass1' not in str(fp)}
bam_paths = {fp.parent.name: fp 
             for fp in base_dir.glob('**/alignment.bam')}

