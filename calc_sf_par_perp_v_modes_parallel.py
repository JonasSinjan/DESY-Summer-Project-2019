"""
TODO:
  look throuhgh all lines and understand how it works
  find out how data files work/ do I need some examples to run?
FIXME:
"""

"""
Code File that interprets synthetic data by generating the structure function
Structure Function is analagous to Energy spectrum, but in real space, not in fourier space
"""

import struct  # performs conversions between Python values and C structs represented as Python strings
import numpy as np
# import matplotlib.pyplot as plt
from multiprocessing import Pool  # Pool class represents a pool of worker processes
# its methods allows tasks to be offloaded to the worker processes in a few ways
import math
import scipy.interpolate as spint
from numpy.random import seed
from numpy.random import randint
from numpy.random import rand

# from pylab import *

xpt = 512
ypt = 512
zpt = 512
lent = 512
Lx = 1.0
Ly = 1.0
Lz = 1.0
t_start = 5
t_stop = 15
step = 1
dir_data = "../512_runs_4/C4/data/decomped_modes/"  # data files
dir_output = "../512_runs_4/C4/data/decomped_modes/"  # data files
seed(1)
n_avg_bfield_pts = 5
nrandpts = 50000
mode = 'F'
nprocs = 16

# initliasing 1D arrays
magk = np.zeros(lent / 2)
mag_power_spec = np.zeros(lent / 2)
kin_power_spec = np.zeros(lent / 2)
mag_pspec_tavg = np.zeros(lent / 2)
kin_pspec_tavg = np.zeros(lent / 2)
# sf_snapshot = np.zeros((lent/4,3))

ntstp = 0
sf_par = np.zeros(lent / 2)
sf_perp = np.zeros(lent / 2)
npts = np.zeros(lent / 2)


def struc_funk(ff):
    ll = ff * 1.0
    print(ll)

    numpt = 0.0
    sf_pare = 0.0
    sf_perpe = 0.0
    for kup in range(0, nrandpts):
        # print(kup)

        # choose a random point
        ri = randint(0, lent, size=3)

        # calculate the average b field direction in a sphere of radius ll around random point (xi,yi,zi)
        # do this by looping over npts_avg_field
        lr = rand(n_avg_bfield_pts) * ll / 2.0
        theta = rand(
            n_avg_bfield_pts) * np.pi  # not sure if this is random, think you must take
        # random value of cos(theta) to get uniform distribution over a sphere
        phi = rand(n_avg_bfield_pts) * 2.0 * np.pi  # rand(5) - random nummber in certain shape array
        lx = lr * np.sin(theta) * np.cos(phi)
        ly = lr * np.sin(theta) * np.sin(phi)
        lz = lr * np.cos(theta)
        xis = np.int_(np.floor(ri[0] + lx))
        yis = np.int_(np.floor(ri[1] + ly))
        zis = np.int_(np.floor(ri[2] + lz))

        # renormalize indexes if they are going out of box
        xis = xis % lent
        yis = yis % lent
        zis = zis % lent

        bhat = np.array([np.mean(bx[xis, yis, zis]), np.mean(by[xis, yis, zis]), np.mean(bz[xis, yis, zis])])
        bhat = bhat / (np.sqrt(np.sum(bhat * bhat)))

        # bhat contains the unit vector along the local magnetic field
        # now take 2 points separated by distance ll along this direction with center at ri
        r1 = np.int_(ri + (ll / 2.0) * bhat)
        r2 = np.int_(ri - (ll / 2.0) * bhat)

        # renormalize indices
        r1 = r1 % lent
        r2 = r2 % lent

        b1 = np.array([vx[r1[0], r1[1], r1[2]], vy[r1[0], r1[1], r1[2]], vz[r1[0], r1[1], r1[2]]])
        b2 = np.array([vx[r2[0], r2[1], r2[2]], vy[r2[0], r2[1], r2[2]], vz[r2[0], r2[1], r2[2]]])

        sf_pare = sf_pare + np.sum((b1 - b2) * (b1 - b2))

        # next calculating the perpendicular structure function
        # make a unit vector perpendicular to bhat
        # generate another random vector
        vrand = rand(3)
        # taking cross product with bhat to generate perpendicular vector
        perpb_x = bhat[1] * vrand[2] - bhat[2] * vrand[1]
        perpb_y = bhat[2] * vrand[0] - bhat[0] * vrand[2]
        perpb_z = bhat[0] * vrand[1] - bhat[1] * vrand[0]
        perpb = np.array([perpb_x, perpb_y, perpb_z])
        perpb = perpb / (np.sqrt(np.sum(perpb * perpb)))

        # now take 2 points separated by distance ll along the perpendicular direction with center at ri
        r1 = np.int_(ri + (ll / 2.0) * perpb)
        r2 = np.int_(ri - (ll / 2.0) * perpb)

        # renormalize indices
        r1 = r1 % lent
        r2 = r2 % lent

        b1 = np.array([vx[r1[0], r1[1], r1[2]], vy[r1[0], r1[1], r1[2]], vz[r1[0], r1[1], r1[2]]])
        b2 = np.array([vx[r2[0], r2[1], r2[2]], vy[r2[0], r2[1], r2[2]], vz[r2[0], r2[1], r2[2]]])

        sf_perpe = sf_perpe + np.sum((b1 - b2) * (b1 - b2))

        numpt = numpt + 1.0

    # print(ll,numpt,sf_pare,sf_perpe)
    return [numpt, sf_pare, sf_perpe]


for t in range(t_start, t_stop + 1, step):  # the time loop

    filename = dir_data + 'B' + mode + str(t) + '.BIN'
    print(filename)
    fd = open(filename, 'rb')
    fd.read(4)
    nx = np.fromfile(file=fd, dtype=np.int32, count=1)[0]
    ny = np.fromfile(file=fd, dtype=np.int32, count=1)[0]
    nz = np.fromfile(file=fd, dtype=np.int32, count=1)[0]
    print(nx, ny, nz)
    fd.read(4)
    fd.read(4)
    abx = np.fromfile(file=fd, dtype=np.float64, count=nx * ny * nz)
    fd.read(4)
    fd.read(4)
    aby = np.fromfile(file=fd, dtype=np.float64, count=nx * ny * nz)
    fd.read(4)
    fd.read(4)
    abz = np.fromfile(file=fd, dtype=np.float64, count=nx * ny * nz)
    fd.read(4)

    temp = np.reshape(abx, (lent, lent, lent))
    bx = temp.transpose()
    temp = np.reshape(aby, (lent, lent, lent))
    by = temp.transpose()
    temp = np.reshape(abz, (lent, lent, lent))
    bz = temp.transpose()

    filename = dir_data + 'V' + mode + str(t) + '.BIN'
    print(filename)
    fd = open(filename, 'rb')
    fd.read(4)
    nx = np.fromfile(file=fd, dtype=np.int32, count=1)[0]
    ny = np.fromfile(file=fd, dtype=np.int32, count=1)[0]
    nz = np.fromfile(file=fd, dtype=np.int32, count=1)[0]
    fd.read(4)
    fd.read(4)
    abx = np.fromfile(file=fd, dtype=np.float64, count=nx * ny * nz)
    fd.read(4)
    fd.read(4)
    aby = np.fromfile(file=fd, dtype=np.float64, count=nx * ny * nz)
    fd.read(4)
    fd.read(4)
    abz = np.fromfile(file=fd, dtype=np.float64, count=nx * ny * nz)
    fd.read(4)

    temp = np.reshape(abx, (lent, lent, lent))
    vx = temp.transpose()
    temp = np.reshape(aby, (lent, lent, lent))
    vy = temp.transpose()
    temp = np.reshape(abz, (lent, lent, lent))
    vz = temp.transpose()

    if __name__ == '__main__':
        pool = Pool(processes=nprocs)
        sf_snapshot = pool.map(struc_funk, range(lent / 4))
        sff = np.asarray(sf_snapshot)
        pool.terminate()

        # pool = Pool(processes=nprocs)              # start 4 worker processes

    # sf_snapshot = pool.map(struc_funk, range(lent/4))

    # sff = np.asarray(sf_snapshot)

    npts[0:lent / 4] = npts[0:lent / 4] + sff[:, 0]
    sf_par[0:lent / 4] = sf_par[0:lent / 4] + sff[:, 1]
    sf_perp[0:lent / 4] = sf_perp[0:lent / 4] + sff[:, 2]

    # print(np.shape(sf_snapshot))
    # for q in range (0,lent/4) :
    #  print("sf_snapshot= ",q, sf_snapshot[q])

sf_par = sf_par / npts
sf_perp = sf_perp / npts

# writing the spectra to a file
f = open(dir_output + 'sf_par_perp_v_' + mode + '.txt', 'w')
for i in range(0, lent / 2):
    value = str(i * 1.0) + " " + str(sf_par[i]) + " " + str(sf_perp[i])
    f.write(value + "\n")
f.close()
