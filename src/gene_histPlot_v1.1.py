#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Zhihao Xie"
__version__ = "v1.1.0"

import os,sys
import re

if len(sys.argv) < 2:
    sys.stderr.write("Usage: python3 %s <gene_seq[fasta]> <out_prefix>\n" % sys.argv[0])
    sys.stderr.write("\t<out_prefix> default is \"gene\"\n")
    sys.stderr.flush()
    sys.exit()

import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO

if len(sys.argv) <= 2:
    out_tag = "gene"
else:
    out_tag = sys.argv[2]

if re.search(r"^\.|^/", out_tag):
    out_dir = os.path.dirname(out_tag)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    else:
        sys.stderr.write("Warning! The out directory had exist.\n")
        sys.stderr.flush()

lengthFile = out_tag+".length.txt"
lengthFile_out = open(lengthFile, "w")
lengthFile_out.write("#SeqID\tLength\n")
for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
    seqID = seq_record.id
    seqLength = len(seq_record)
    lengthFile_out.write("{0}\t{1}\n".format(seqID, seqLength))
    lengthFile_out.flush()
lengthFile_out.close()

file = open(lengthFile,'r')
x = []
for line in file:
    linarr = line.strip().split("\t",1)
    if re.match(r'^#|^#?SeqID', linarr[0]):
        continue
    if float(linarr[1]) > 3000:
        linarr[1] = 3000
    x.append(int(float(linarr[1])))

file.close()

# hist plot
fig, ax1 = plt.subplots()

num_bins = 30
# the histogram of the data
n, bins, patches = ax1.hist(x, num_bins, normed=False, histtype='bar', facecolor='steelblue', alpha=0.75, edgecolor='black', linewidth=0.5, label='Frequence(#)')

ax1.set_xlabel("Gene length(bp)", fontsize = 10)
ax1.set_ylabel("Frequence(#)", fontsize = 10)
ax1.set_title("Gene Length Distribution")

xticks_label = list(ax1.get_xticks())
xticks_label = [int(x) for x in xticks_label]
xticks_label[-2] = ">3000"
ax1.set_xticklabels(xticks_label)
#ax.legend(["Frequence"], loc="upper right", fontsize='x-small', bbox_to_anchor=(0.982,0.98))

#print(len(n))
#print(len(bins))
#print(patches)

percent = []
for item in n:
    percent.append(item/sum(n)*100)
percent.insert(0, 0)
max_item = max(percent)

center = (bins[:-1] + bins[1:]) / 2
center = np.insert(center, 0, 0)

ax2 = ax1.twinx()
ax2.plot(center, percent, 'r-', linewidth=0.7, label='Percentage(%)')
ax2.set_ylim(bottom=0, top=max_item*1.2)
ax2.set_ylabel("Percentage(%)",fontsize = 10)
#ax2.legend(["Percentage"], loc="upper left", fontsize='x-small', bbox_to_anchor=(0.8,0.9))

# set legend
legend1=ax1.legend(loc=(0.8,0.94),fontsize=5,shadow=None)
legend2=ax2.legend(loc=(0.8,0.9),fontsize=5,shadow=None)
legend1.get_frame().set_facecolor('#FFFFFF')
legend2.get_frame().set_facecolor('#FFFFFF')

#fig.set_size_inches(6,4)

#plt.grid(True)
#fig.tight_layout()
#plt.show()

fig.savefig("%s.length_distribution.pdf" % (out_tag), dpi=300)

cmd = r"convert -density 300 -resize 70% " + out_tag + ".length_distribution.pdf " + out_tag + ".length_distribution.png"
os.system(cmd)

