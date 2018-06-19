## plotting utils for network analyses
## Author: Zichen Wang
## 3/27/2014

import os
import operator
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.lines import Line2D
from math import sqrt
from itertools import combinations

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 ## Output Type 3 (Type3) or Type 42 (TrueType)
rcParams['font.sans-serif'] = 'Arial'
rcParams['lines.markersize'] = 4
from matplotlib_venn import venn2, venn3

def plot_venn(t1=None, t2=None, t3=None, ax=None, set_colors=('r', 'b', 'k')):
	"""input: 2 or 3 tuples: (list/set, name_to_display) """
	assert len(t1) == len(t2) == 2
	if t3:
		venn3( [set(t[0]) for t in [t1,t2,t3]], tuple( ['%s\n(%s)'%(t[1], len(set(t[0])) ) for t in [t1,t2,t3]]) , set_colors=set_colors, alpha=0.5,ax=ax)
	else:
		venn2( [set(t[0]) for t in [t1,t2]], tuple( ['%s\n(%s)'%(t[1], len(set(t[0])) ) for t in [t1,t2]]), set_colors=set_colors[0:2],alpha=0.5, ax=ax)


COLORS8 = [
'#1F77B4',
'#26A9E0',
'#75DBA7',
'#2CA02C',
'#9467BD',
'#FF0000',
'#FF7F0E',
'#E377C2',
]

COLORS = ['#1F78B4','#E31A1C','#FF7F00','#6A3D9A','#B15928','#666666', 'k']
COLORS10 = [
'#1f77b4',
'#ff7f0e',
'#2ca02c',
'#d62728',
'#9467bd',
'#8c564b',
'#e377c2',
'#7f7f7f',
'#bcbd22',
'#17becf',
]
COLORS20 = [
'#1f77b4',
'#aec7e8',
'#ff7f0e',
'#ffbb78',
'#2ca02c',
'#98df8a',
'#d62728',
'#ff9896',
'#9467bd',
'#c5b0d5',
'#8c564b',
'#c49c94',
'#e377c2',
'#f7b6d2',
'#7f7f7f',
'#c7c7c7',
'#bcbd22',
'#dbdb8d',
'#17becf',
'#9edae5',
]

COLORS20b = [
'#393b79',
'#5254a3',
'#6b6ecf',
'#9c9ede',
'#637939',
'#8ca252',
'#b5cf6b',
'#cedb9c',
'#8c6d31',
'#bd9e39',
'#e7ba52',
'#e7cb94',
'#843c39',
'#ad494a',
'#d6616b',
'#e7969c',
'#7b4173',
'#a55194',
'#ce6dbd',
'#de9ed6',
]

def enlarge_tick_fontsize(ax,fontsize):
	for tick in ax.xaxis.get_major_ticks():
		tick.label.set_fontsize(fontsize)
	for tick in ax.yaxis.get_major_ticks():
		tick.label.set_fontsize(fontsize)

def scatter_plot(coords, samples, sep='_', legend_size=14,
	marker='o', alpha=0.8):

	fig = plt.figure(figsize=(10,8))
	ax = fig.add_subplot(111)
	conditions = [s.split(sep)[0] for s in samples]
	conditions_uniq = list(set(conditions))
	for idx, condition in enumerate(conditions_uniq):
		mask = np.in1d(conditions, [condition])
		ax.scatter(coords[mask, 0], coords[mask, 1], label=condition, color=COLORS10[idx], 
			marker=marker, alpha=alpha)

	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	enlarge_tick_fontsize(ax, 14)
	fig.tight_layout()
	return fig

