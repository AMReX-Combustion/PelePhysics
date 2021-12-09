import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def prettyLabels(xlabel,ylabel,fontsize,title=None):
    plt.xlabel(xlabel, fontsize=fontsize, fontweight='bold', fontname="Times New Roman")
    plt.ylabel(ylabel, fontsize=fontsize, fontweight='bold', fontname="Times New Roman")
    if not title==None:
        plt.title(title, fontsize=fontsize, fontweight='bold', fontname="Times New Roman")
    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontname("Times New Roman")
        tick.label1.set_fontweight('bold')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontname("Times New Roman")
        tick.label1.set_fontweight('bold')
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)
        ax.spines[axis].set_color("black")
    plt.grid(color='k', linestyle='-', linewidth=0.5)
    plt.tight_layout()

def plotLegend():
    fontsize = 16
    plt.legend()
    leg=plt.legend(prop={'family':'Times New Roman','size': fontsize-3,'weight':'bold' })
    leg.get_frame().set_linewidth(2.0)
    leg.get_frame().set_edgecolor('k')    

