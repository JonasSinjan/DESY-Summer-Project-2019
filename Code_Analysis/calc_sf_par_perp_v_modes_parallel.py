"""
TODO:
    check alternative method for 2D
    check that data from displacement + squares method generates correct spectrum, try 2D, 128 and 256 only
    create different dir for 128 and 256
FIXME:
    maybe set up json to use for setup instead of ugly setup at beginning
    setup correct data input and output paths
    not sure if second argument passed through correctly to the struc function
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
from functools import partial
import math
import scipy.interpolate as spint
from numpy.random import seed
from numpy.random import randint
from numpy.random import rand

# from pylab import *

# testing git properly setup

##############################################################################################
#####                                      SETUP                                       #######
##############################################################################################

# NUMBER OF POINTS: OPTIONS 128, 256, 512 ETC
size = 64

sq_bool = False

if sq_bool:
    nx = size
    ny = size
    nz = size

else:
    nx = size + 1
    ny = size + 1
    nz = size + 1

# DATA INPUT AND OUTPUT PATH
dir_data = "/lustre/fs23/group/that/jonas/Github_repo/DESY/3d_displacement/64run3D/"  # data files
dir_output = "/lustre/fs23/group/that/jonas/Github_repo/DESY/3d_displacement/64run3D/"  # data files
#dir_data = "c:/Users/jonas/DESY/2d_displacement/256run2D_73_frac/"  # data files
#dir_output = "c:/Users/jonas/DESY/2d_displacement/256run2D_73_frac/"  # data files

# IF DISPLACEMENT MUST PULL FROM PHI.BIN FOR SQUARES RHO.BIN -  CHECK!!!

# NUMBER OF DIMENSIONS
twoD_bool = False  # if set to true, will assume data in 2D, otherwise when false defaults to 3D

if twoD_bool is True:
    shape = (lent + 1, lent + 1)  # for 2D
else:
    shape = (lent + 1, lent + 1, lent + 1)

###############################################################################################

xpt = size
ypt = size
zpt = size

Lx = 1.0
Ly = 1.0
Lz = 1.0
t_start = 5
t_stop = 6  # only want one loop
step = 1

seed(1)
n_avg_bfield_pts = 5
nrandpts = 10000
mode = 'F'
nprocs = 16

# initliasing 1D arrays
magk = np.zeros(int(lent / 2))
mag_power_spec = np.zeros(int(lent / 2))
kin_power_spec = np.zeros(int(lent / 2))
mag_pspec_tavg = np.zeros(int(lent / 2))
kin_pspec_tavg = np.zeros(int(lent / 2))
# sf_snapshot = np.zeros((lent/4,3))

ntstp = 0
sf_par = np.zeros(int(lent / 2))
sf_perp = np.zeros(int(lent / 2))
npts = np.zeros(int(lent / 2))


def read_files(dir_data):
    filename = dir_data + 'PHI0' + '.BIN'
    # print(filename)
    fd = open(filename, 'rb')

    abx = np.fromfile(file=fd, dtype=np.float64, count=nx * ny)

    temp = np.reshape(abx, (nx, ny))
    phi = temp.transpose()  # missed the empty brackets here
    # print(phi[22,:]) - working correctly

    filename = dir_data + 'BX' + '.BIN'  # 'B' + mode + str(t) + '.BIN' not sure why this was used:
    # print(filename)
    fd = open(filename, 'rb')

    abx = np.fromfile(file=fd, dtype=np.float64, count=nx * ny)

    temp = np.reshape(abx, (nx, ny))
    bx = temp.transpose()

    bx.fill(1)

    filename = dir_data + 'BY' + '.BIN'
    # print(filename)
    fd = open(filename, 'rb')

    aby = np.fromfile(file=fd, dtype=np.float64, count=nx * ny)

    temp = np.reshape(aby, (nx, ny))
    by = temp.transpose()
    print(bx[:, 1])
    print(np.mean(bx), np.mean(by))
    return phi, bx, by

def read_files3D(dir_data):
    filename = dir_data + 'PHI0' + '.BIN'
    # print(filename)
    fd = open(filename, 'rb')

    abx = np.fromfile(file=fd, dtype=np.float64, count=nx * ny * nz)

    temp = np.reshape(abx, (nx, ny, nz))
    phi = temp.transpose()  # missed the empty brackets here
    # print(phi[22,:]) - working correctly

    filename = dir_data + 'BX' + '.BIN'  # 'B' + mode + str(t) + '.BIN' not sure why this was used:
    # print(filename)
    fd = open(filename, 'rb')

    abx = np.fromfile(file=fd, dtype=np.float64, count=nx * ny * nz)

    temp = np.reshape(abx, (nx, ny, nz))
    bx = temp.transpose()

    bx.fill(1)

    filename = dir_data + 'BY' + '.BIN'
    # print(filename)
    fd = open(filename, 'rb')

    aby = np.fromfile(file=fd, dtype=np.float64, count=nx * ny * nz)

    temp = np.reshape(aby, (nx, ny, nz))
    by = temp.transpose()

    filename = dir_data + 'BZ' + '.BIN'
    # print(filename)
    fd = open(filename, 'rb')

    aby = np.fromfile(file=fd, dtype=np.float64, count=nx * ny * nz)

    temp = np.reshape(aby, (nx, ny, nz))
    bz = temp.transpose()

    print(bx[:, :, 1])
    print(np.mean(bx), np.mean(by), np.mean(by), np.mean(bz+by))
    return phi, bx, by, bz

def struc_funk3D(ff, phi, bx, by, bz):
    ll = ff * 1.0

    numpt = 0.0
    sf_pare = 0.0
    sf_perpe = 0.0

    for kup in range(0, nrandpts):
        # 3D method
        # choose a random point
        ri = randint(0, lent, size=3)  # 1x3 array

        # calculate the average b field direction in a sphere of radius ll around random point (xi,yi,zi)
        # do this by looping over npts_avg_field
        lr = rand(n_avg_bfield_pts) * ll / 2.0
        theta = rand(n_avg_bfield_pts) * np.pi  # want random theta not random costheta
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
        # for global bhat becomes x, y, z = 1, 0, 0 bhat is unit vector pointing in the x direction
        bhat = bhat / (np.sqrt(np.sum(bhat * bhat)))

        # bhat contains the unit vector along the local magnetic field
        # now take 2 points separated by distance ll along this direction with center at ri
        r1 = np.int_(ri + (ll / 2.0) * bhat)
        r2 = np.int_(ri - (ll / 2.0) * bhat)

        # renormalize indices
        r1 = r1 % lent
        r2 = r2 % lent

        b1 = np.array([phi[r1[0], r1[1], r1[2]]])
        b2 = np.array([phi[r2[0], r2[1], r2[2]]])

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

        # b1 = np.array([vx[r1[0], r1[1], r1[2]], vy[r1[0], r1[1], r1[2]], vz[r1[0], r1[1], r1[2]]])
        # b2 = np.array([vx[r2[0], r2[1], r2[2]], vy[r2[0], r2[1], r2[2]], vz[r2[0], r2[1], r2[2]]])

        # b1 will be phi at point r1 and b2 will be phi at point r2
        b1 = np.array([phi[r1[0], r1[1], r1[2]]])
        b2 = np.array([phi[r2[0], r2[1], r2[2]]])

        sf_perpe = sf_perpe + np.sum((b1 - b2) * (b1 - b2))

        numpt = numpt + 1.0

    print(ll, numpt, sf_pare, sf_perpe, )
    return [numpt, sf_pare, sf_perpe]


    def struc_funk2D(ff, phi, bx, by):
        ll = ff * 1.0

        numpt = 0.0
        sf_pare = 0.0
        sf_perpe = 0.0

        for kup in range(0, nrandpts):
            # 2D method
            # choose a random point
            ri = randint(0, lent, size=2)  # 1x2 array

            # 2D polars
            lr = rand(n_avg_bfield_pts) * ll / 2.0
            theta = rand(n_avg_bfield_pts) * np.pi  # want random theta not random costheta
            lx = lr * np.cos(theta)
            ly = lr * np.sin(theta)
            xis = np.int_(np.floor(ri[0] + lx))
            yis = np.int_(np.floor(ri[1] + ly))

            # renormalize indexes if they are going out of box
            xis = xis % lent
            yis = yis % lent

            bhat = np.array([np.mean(bx[xis, yis]), np.mean(by[xis, yis])])
            bhat = bhat / (np.sqrt(np.sum(bhat * bhat)))

            # bhat contains the unit vector along the local magnetic field
            # now take 2 points separated by distance ll along this direction with center at ri
            r1 = np.int_(ri + (ll / 2.0) * bhat)
            r2 = np.int_(ri - (ll / 2.0) * bhat)

            # renormalize indices
            r1 = r1 % lent
            r2 = r2 % lent

            # dont need v, only need phi now
            b1 = np.array([phi[r1[0], r1[1]]])
            b2 = np.array([phi[r2[0], r2[1]]])

            sf_pare = sf_pare + np.sum((b1 - b2) * (b1 - b2))

            # find one of the 2D perpendicular vectors
            perpb = np.array([bhat[1], -bhat[0]])

            # print(ri, type(ri), ll, type(ll), perpb, type(perpb))
            # now take 2 points separated by distance ll along the perpendicular direction with center at ri
            r1 = np.int_(ri + (ll / 2.0) * perpb)
            r2 = np.int_(ri - (ll / 2.0) * perpb)

            # renormalize indices
            r1 = r1 % lent
            r2 = r2 % lent

            b1 = np.array([phi[r1[0], r1[1]]])
            b2 = np.array([phi[r2[0], r2[1]]])

            sf_perpe = sf_perpe + np.sum((b1 - b2) * (b1 - b2))

            numpt = numpt + 1.0

    print(ll, numpt, sf_pare, sf_perpe, )
    return [numpt, sf_pare, sf_perpe]


for t in range(0, 1, 1):  # the time loop

    if __name__ == '__main__':
        #phi, bx, by = read_files(dir_data)  # will these be recognised by the struc funk function?
        phi0, bx,by,bz = read_files3D(dir_data)

        # pool = Pool(processes=nprocs)
        # sf_snapshot = pool.map(struc_funk, range(lent / 4)) #ff/ll is the distance taken
        sf_snapshot = []
        sff = np.zeros((3, int(lent / 4)))
        for i in range(int(lent / 4)):
            numpt_tmp, par_tmp, perp_tmp = struc_funk3D(i, phi0, bx, by, bz)
            sff[0, i] = numpt_tmp
            sff[1, i] = par_tmp
            sff[2, i] = perp_tmp

        # pool.terminate()
        print(np.shape(sff))
        print("The Process has Completed")

        # pool = Pool(processes=nprocs)              # start 4 worker processes

    # sf_snapshot = pool.map(struc_funk, range(lent/4))

    # sff = np.asarray(sf_snapshot)

    npts[0:int(lent / 4)] = npts[0:int(lent / 4)] + sff[0, :]
    sf_par[0:int(lent / 4)] = sf_par[0:int(lent / 4)] + sff[1, :]
    sf_perp[0:int(lent / 4)] = sf_perp[0:int(lent / 4)] + sff[2, :]

    # print(np.shape(sf_snapshot))
    # for q in range (0,lent/4) :
    #  print("sf_snapshot= ",q, sf_snapshot[q])

sf_par = sf_par / npts
sf_perp = sf_perp / npts

# writing the spectra to a file
f = open(dir_output + 'sf_par_perp_v_phi0' + mode + '.txt', 'w')
for i in range(0, int(lent / 2)):
    value = str(i * 1.0) + " " + str(sf_par[i]) + " " + str(sf_perp[i])
    f.write(value + "\n")
f.close()
