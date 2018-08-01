import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_table("MasterTable.txt", sep="\t", index_col=0)

binBoundaries = range(0, 2200000, 50000)

mean = df['length'].mean()
median = df['length'].median()
print("Plasmid mean length:{}bp\nPlasmid median length:{}bp".format(mean, median))

df.replace(to_replace='st.tic', value='non-mosaic', regex=True, inplace=True)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3, 3))
plt.boxplot([df[df['character'] == 'mosaic']['length'], df[df['character'] == 'non-mosaic']['length']],
            labels=['mosaic', 'non-mosaic'], boxprops=dict(color='c'),
            flierprops=dict(marker='.', markerfacecolor='white', markeredgecolor='black'),
            whiskerprops=dict(linestyle='-', color='c'), medianprops=dict(color='black'))

ax.set_ylabel("Plasmid length (bp)")
plt.gca().set_yscale("log")
plt.tight_layout()
plt.savefig("LengthBoxplot.png")
