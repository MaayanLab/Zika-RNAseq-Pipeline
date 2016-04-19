## plotting utils for network analyses
## Author: Zichen Wang
## 3/27/2014

import os
import operator
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.lines import Line2D
from matplotlib_venn import venn2, venn3
from math import sqrt
from itertools import combinations

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 ## Output Type 3 (Type3) or Type 42 (TrueType)
rcParams['font.sans-serif'] = 'Arial'

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

def plot_venn(t1=None, t2=None, t3=None, ax=None, set_colors=('r', 'b', 'k')):
	"""input: 2 or 3 tuples: (list/set, name_to_display) """
	assert len(t1) == len(t2) == 2
	if t3:
		venn3( [set(t[0]) for t in [t1,t2,t3]], tuple( ['%s\n(%s)'%(t[1], len(set(t[0])) ) for t in [t1,t2,t3]]) , set_colors=set_colors, alpha=0.5,ax=ax)
	else:
		venn2( [set(t[0]) for t in [t1,t2]], tuple( ['%s\n(%s)'%(t[1], len(set(t[0])) ) for t in [t1,t2]]), set_colors=set_colors[0:2],alpha=0.5, ax=ax)

def plot_venn3(t1,t2,t3, fn):
	assert len(t1) == len(t2) == len(t3) == 2
	pwd = os.getcwd()
	pwd = '/'.join( pwd.split('\\') )

	sets = [ set(t[0]) for t in [t1, t2, t3] ]
	names = [ t[1] for t in [t1, t2, t3] ]

	areas = [str(len(s)) for s in sets]
	names_areas = [n +'('+ a +')' for n, a in zip(names, areas)]
	n2s = [len(sets[i1] & sets[i2]) for i1, i2 in combinations(range(3), 2)]
	n123 = len(sets[0] & sets[1] & sets[2])
	all_numbers = areas + n2s + [n123]
	all_values = [pwd] + all_numbers + names_areas + [fn]
	with open ('temp_venn.R','w') as out:
		out.write('''
setwd("%s")
library("grid")
library("VennDiagram")
venn.plot <- draw.triple.venn(
	area1 = %s,
	area2 = %s,
	area3 = %s,
	n12 = %s,
	n13 = %s,
	n23 = %s,
	n123 = %s,
	category = c("%s", "%s", "%s"),
	fill = c("red", "magenta", "blue"),
	lty = "dashed",
	cex = 2,
	cat.cex = 2,
	cat.col = c("red", "magenta", "blue")
	);
pdf("%s", width=6, height=6);
grid.draw(venn.plot);
dev.off();
		'''%tuple(all_values))
	parameters = ['Rscript', 'temp_venn.R']
	print parameters
	p = subprocess.Popen(parameters)
	p.wait()
	return


def plot_venn4(t1,t2,t3,t4, fn):
	assert len(t1) == len(t2) == len(t3) == len(t4) == 2
	pwd = os.getcwd()
	pwd = '/'.join( pwd.split('\\') )

	sets = [ set(t[0]) for t in [t1, t2, t3, t4] ]
	names = [ t[1] for t in [t1, t2, t3, t4] ]

	areas = [str(len(s)) for s in sets]
	names_areas = [n +'\n('+ a +')' for n, a in zip(names, areas)]

	n2s = [len(sets[i1] & sets[i2]) for i1, i2 in combinations(range(4), 2)]
	n3s = [len(sets[i1] & sets[i2] & sets[i3]) for i1, i2, i3 in combinations(range(4), 3)]
	n1234 = len(sets[0] & sets[1] & sets[2] & sets[3])

	all_numbers = areas + n2s + n3s + [n1234]
	all_values = [pwd] + all_numbers + names_areas + [fn]

	with open ('temp_venn.R','w') as out:
		out.write('''
setwd("%s")
library("grid")
library("VennDiagram")
venn.plot <- draw.quad.venn(
	area1 = %s,
	area2 = %s,
	area3 = %s,
	area4 = %s,
	n12 = %s,
	n13 = %s,
	n14 = %s,
	n23 = %s,
	n24 = %s,
	n34 = %s,
	n123 = %s,
	n124 = %s,
	n134 = %s,
	n234 = %s,
	n1234 = %s,
	category = c("%s", "%s", "%s", "%s"),
	fill = c("orange", "red", "magenta", "blue"),
	lty = "dashed",
	cex = 2,
	cat.cex = 2,
	cat.col = c("orange", "red", "magenta", "blue")
	);
pdf("%s", width=14, height=8);
grid.draw(venn.plot);
dev.off();
		'''%tuple(all_values))
	parameters = ['Rscript', 'temp_venn.R']
	print parameters
	p = subprocess.Popen(parameters)
	p.wait()
	return

