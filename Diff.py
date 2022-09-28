#!/usr/bin/env python3
import numpy as np 
import warnings
import matplotlib.pyplot as plt
from numpy import savetxt, loadtxt
import freud


b=np.loadtxt("LJ.txt")
N=len(b)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    box_data = np.genfromtxt("Box.txt", skip_header=0, max_rows=3)
    data = np.genfromtxt("LJ.txt", skip_header=0, invalid_raise=False)

data = data[~np.isnan(data).all(axis=1)].reshape(-1, N, 3)

box = freud.box.Box.from_box(box_data[:, 1] - box_data[:, 0])

data[..., :3] -= box.L / 2


for frame in data:    
    dp = freud.diffraction.DiffractionPattern(grid_size=1024, output_size=1024)
    fig, ax = plt.subplots(figsize=(3, 3), dpi=150)
    dp.compute((box, frame[:, :3]), view_orientation=[1, 0, 0,0])
    dp.plot(ax,cmap='bone')
    plt.savefig("DP.png")
    


