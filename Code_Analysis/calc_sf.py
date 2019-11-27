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
import time

def read_files(dir_data):
    filename = dir_data + 'PHI' + '.BIN'
    fd = open(filename, 'rb')

    abx = np.fromfile(file=fd, dtype=np.float64, count=nx * ny)

    temp = np.reshape(abx, (nx, ny))
    phi = temp.transpose()

    filename = dir_data + 'BX' + '.BIN'
    fd = open(filename, 'rb')

    abx = np.fromfile(file=fd, dtype=np.float64, count=nx * ny)

    temp = np.reshape(abx, (nx, ny))
    bx = temp.transpose()

    #bx.fill(1)

    filename = dir_data + 'BY' + '.BIN'
    fd = open(filename, 'rb')

    aby = np.fromfile(file=fd, dtype=np.float64, count=nx * ny)

    temp = np.reshape(aby, (nx, ny))
    by = temp.transpose()

    print(bx[:, 1])
    print(np.mean(bx), np.mean(by))
    return phi, bx, by

def read_files_phi0(dir_data):
    filename = dir_data + 'PHI0' + '.BIN'
    fd = open(filename, 'rb')

    abx = np.fromfile(file=fd, dtype=np.float64, count=nx * ny)

    temp = np.reshape(abx, (nx, ny))
    phi0 = temp.transpose() 

    filename = dir_data + 'BX' + '.BIN'
    fd = open(filename, 'rb')

    abx = np.fromfile(file=fd, dtype=np.float64, count=nx * ny)

    temp = np.reshape(abx, (nx, ny))
    bx = temp.transpose()

    bx.fill(1)

    filename = dir_data + 'BY' + '.BIN'
    fd = open(filename, 'rb')

    aby = np.fromfile(file=fd, dtype=np.float64, count=nx * ny)

    temp = np.reshape(aby, (nx, ny))
    by = temp.transpose()

    by.fill(0)

    print(bx[:, 1])
    print(np.mean(bx), np.mean(by))
    return phi0, bx, by

def read_files_sq(dir_data):
    filename = dir_data + 'RHO' + '.BIN'
    fd = open(filename, 'rb')

    abx = np.fromfile(file=fd, dtype=np.float64, count=nx * ny)

    temp = np.reshape(abx, (nx, ny))
    phi = temp.transpose()
  
    filename = dir_data + 'BX' + '.BIN' 
    fd = open(filename, 'rb')

    abx = np.fromfile(file=fd, dtype=np.float64, count=nx * ny)

    temp = np.reshape(abx, (nx, ny))
    bx = temp.transpose()

    filename = dir_data + 'BY' + '.BIN'
    fd = open(filename, 'rb')

    aby = np.fromfile(file=fd, dtype=np.float64, count=nx * ny)

    temp = np.reshape(aby, (nx, ny))
    by = temp.transpose()

    print(bx[:, 1])
    print(np.mean(bx), np.mean(by))
    return phi, bx, by

def read_files3D_phi0(dir_phi0, dir_B, local = False):
    filename = dir_phi0 + 'PHI0' + '.BIN'
    fd = open(filename, 'rb')
    abx = np.fromfile(file=fd, dtype=np.float64, count=nx * ny * nz)
    temp = np.reshape(abx, (nx, ny, nz))
    phi0 = temp.transpose()

    if local:
    
        filename = dir_B + 'BX' + '.BIN' 
        fd = open(filename, 'rb')
        abx = np.fromfile(file=fd, dtype=np.float64, count=nx * ny * nz)
        temp = np.reshape(abx, (nx, ny, nz))
        bx = temp.transpose()

        filename = dir_B + 'BY' + '.BIN'
        fd = open(filename, 'rb')
        aby = np.fromfile(file=fd, dtype=np.float64, count=nx * ny * nz)
        temp = np.reshape(aby, (nx, ny, nz))
        by = temp.transpose()

        filename = dir_B + 'BZ' + '.BIN'
        fd = open(filename, 'rb')
        aby = np.fromfile(file=fd, dtype=np.float64, count=nx * ny * nz)
        temp = np.reshape(aby, (nx, ny, nz))
        bz = temp.transpose()

        dbx = bx - 1
        dby = by
        dbz = bz

        tmp = dbx**2 + dby**2 + dbz**2 #calculating the mach alfven number
        tmp = np.mean(tmp)
        mach_alfven = np.sqrt(tmp)

    else:
        bx = np.ones((nx, ny, nz))
        by = np.zeros((nx, ny, nz))
        bz = np.zeros((nx, ny, nz))
    print(bx[:, :, 1])
    print(np.mean(bx), np.mean(by), np.mean(by), np.mean(bz+by))
    return phi0, bx, by, bz, mach_alfven

def read_files3D_phi(dir_data):
    filename = dir_data + 'PHI' + '.BIN'
    fd = open(filename, 'rb')

    abx = np.fromfile(file=fd, dtype=np.float64, count=nx * ny * nz)

    temp = np.reshape(abx, (nx, ny, nz))
    phi = temp.transpose()

    filename = dir_data + 'BX' + '.BIN'  
    fd = open(filename, 'rb')

    abx = np.fromfile(file=fd, dtype=np.float64, count=nx * ny * nz)

    temp = np.reshape(abx, (nx, ny, nz))
    bx = temp.transpose()

    dbx = bx - 1
    
    #bx.fill(1)

    filename = dir_data + 'BY' + '.BIN'
    fd = open(filename, 'rb')

    aby = np.fromfile(file=fd, dtype=np.float64, count=nx * ny * nz)

    temp = np.reshape(aby, (nx, ny, nz))
    by = temp.transpose()

    dby = by
    
    #by.fill(0)

    filename = dir_data + 'BZ' + '.BIN'
    fd = open(filename, 'rb')

    aby = np.fromfile(file=fd, dtype=np.float64, count=nx * ny * nz)

    temp = np.reshape(aby, (nx, ny, nz))
    bz = temp.transpose()

    dbz = bz
    
    #bz.fill(0)

    tmp = dbx**2 + dby**2 + dbz**2 #calculating the mach alfven number
    tmp = np.mean(tmp)
    mach_alfven = np.sqrt(tmp)

    print(bx[:, :, 1])
    print(np.mean(bx), np.mean(by), np.mean(by), np.mean(bz+by))
    return phi, bx, by, bz #, mach_alfven

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
        phi_ang = rand(n_avg_bfield_pts) * 2.0 * np.pi  # rand(5) - random nummber in certain shape array
        lx = lr * np.sin(theta) * np.cos(phi_ang)
        ly = lr * np.sin(theta) * np.sin(phi_ang)
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



if __name__ == '__main__':

    #---------------------------------------------------------------------------------------------
    #  SETUP                                       
    #----------------------------------------------------------------------------------------------

    # data input and output path

    #desy cluster path
    # dir_data = '/lustre/fs23/group/that/jonas/Github_repo/DESY/phi0init/Runs/512_15_kpara/'  # data files
    # dir_output = '/lustre/fs23/group/that/jonas/Github_repo/DESY/phi0init/Runs/512_15_kpara/'  # data files

    dir_phi0 = '/lustre/fs23/group/that/jonas/Github_repo/DESY/phi0init/Runs/512_test/'
    dir_B = '/lustre/fs23/group/that/jonas/Github_repo/DESY/localB/Runs/512_B_amp1/'  # data files
    dir_output = '/lustre/fs23/group/that/jonas/Github_repo/DESY/localB/Runs/512_B_amp1/'  # data files
    
    #windows laptop
    # dir_data = "c:/Users/jonas/DESY/2d_displacement/256run2D_73_frac/"  # data files
    # dir_output = "c:/Users/jonas/DESY/2d_displacement/256run2D_73_frac/"  # data files

    #linux home pc
    # dir_data = '/home/jonas/Documents/VSCode/DESY/3d_disp_mem/Runs/256_test_5-2/'
    # dir_output = '/home/jonas/Documents/VSCode/DESY/3d_disp_mem/Runs/256_test_5-2/'
    # dir_data = '/home/jonas/Documents/VSCode/DESY/phi0init/Runs/512_test/'
    # dir_output = '/home/jonas/Documents/VSCode/DESY/phi0init/Runs/512_test/'

    # resolution size must be specified
    size = 512
    lent = size

    # dimensions
    twoD_bool = False # if set to true, will assume data in 2D, otherwise when false defaults to 3D

    #set to true if decoding data files from squares method - if false defaults to displacment method results
    sq_bool = False

    if sq_bool:
        nx = size
        ny = size
        nz = size

    else:
        nx = size + 1
        ny = size + 1
        nz = size + 1

    xpt, ypt, zpt = size, size, size

    Lx, Ly, Lz = 1.0, 1.0, 1.0

    seed(1)
    n_avg_bfield_pts = 5
    nrandpts = 10000
    mode = 'F'

    if twoD_bool:
        shape = (lent + 1, lent + 1)  # for 2D

        sf_par = np.zeros(int(lent / 2))
        sf_perp = np.zeros(int(lent / 2))
        npts = np.zeros(int(lent / 2))

        #2D PHI0
        #phi0, bx, by = read_files_phi0(dir_data)

        sf_snapshot = []
        sff = np.zeros((3, int(lent / 4)))
        for i in range(int(lent / 4)):
            numpt_tmp, par_tmp, perp_tmp = struc_funk2D(i, phi0, bx, by)
            sff[0, i] = numpt_tmp
            sff[1, i] = par_tmp
            sff[2, i] = perp_tmp
            
        npts[0:int(lent / 4)] = npts[0:int(lent / 4)] + sff[0, :]
        sf_par[0:int(lent / 4)] = sf_par[0:int(lent / 4)] + sff[1, :]
        sf_perp[0:int(lent / 4)] = sf_perp[0:int(lent / 4)] + sff[2, :]

        sf_par = sf_par / npts
        sf_perp = sf_perp / npts

        # writing the spectra to a file
        f = open(dir_output + 'sf_par_perp_v_phi0_wrt_global' + mode + '.txt', 'w')
        for i in range(0, int(lent / 2)):
            value = str(i * 1.0) + " " + str(sf_par[i]) + " " + str(sf_perp[i])
            f.write(value + "\n")
        f.close()

    else:
        shape = (lent + 1, lent + 1, lent + 1)

        sf_par_2 = np.zeros(int(lent / 2))
        sf_perp_2 = np.zeros(int(lent / 2))
        npts_2 = np.zeros(int(lent / 2))

        #3D PHI0
        localbool = True
        phi0,bx,by,bz,mach_alfven= read_files3D_phi0(dir_phi0, dir_B, localbool)
        print('Mach Alfven = ', mach_alfven)
        #time.sleep(30)
	sf_snapshot = []
        sff_2 = np.zeros((3, int(lent / 4)))
        for i in range(int(lent / 4)):
            numpt_tmp, par_tmp, perp_tmp = struc_funk3D(i, phi0, bx, by, bz)
            sff_2[0, i] = numpt_tmp
            sff_2[1, i] = par_tmp
            sff_2[2, i] = perp_tmp

        npts_2[0:int(lent / 4)] = npts_2[0:int(lent / 4)] + sff_2[0, :]
        sf_par_2[0:int(lent / 4)] = sf_par_2[0:int(lent / 4)] + sff_2[1, :]
        sf_perp_2[0:int(lent / 4)] = sf_perp_2[0:int(lent / 4)] + sff_2[2, :]

        sf_par_2 = sf_par_2 / npts_2
        sf_perp_2 = sf_perp_2 / npts_2

        # writing the spectra to a file - must change name of output file depending on phi0 or phi & if wrt global or local frame
        f = open(dir_output + 'sf_par_perp_v_phi0_wrt_local_amp1' + mode + '.txt', 'w')
        for i in range(0, int(lent / 2)):
            value = str(i * 1.0) + " " + str(sf_par_2[i]) + " " + str(sf_perp_2[i]) #+ " " + str(mach_2)
            f.write(value + "\n")
        f.close()

print("The Process has Completed")
