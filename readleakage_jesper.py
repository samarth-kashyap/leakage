# {{{ Library Importsj
from heliosPy import iofuncs as cio
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import time
# }}} imports

# {{{ Global variables
dl_mat = 6
dm_mat = 15
lmax = 249
# }}} global vars


# {{{ def sign(a, b):
def sign(a, b):
    """Returns the the value of a with sign of b

    Inputs:
    -------
    a - float
        variable whose value is returned
    b - float
        variable whose sign is returned

    Returns:
    --------
    abs(a) if b >= 0
    -abs(a) if b < 0

    """
    if b < 0:
        return -abs(a)
    else:
        return abs(a)
# }}} def sign(a, b):


if __name__ == "__main__":
    t1 = time.time()
    data_dir = "/scratch/seismogroup/code/leakvw0/default/"
    write_dir = "/scratch/g.samarth/HMIDATA/leakmat/"
    new_dir = "/scratch/g.samarth/leakage/HMI/"
    new_dir2 = "/scratch/g.samarth/leakage/MDI/jsoc.stanford.edu/SUM9/D367602322/S00000/leakvw.ref"

    jesper_data = np.loadtxt("/home/g.samarth/leakage/shravan_200812.txt")

    rmat1 = fits.open(data_dir + "vradsum/leakrlist.vradsum.fits")[0].data
    rmat2 = fits.open(data_dir + "vradsum/leakilist.vradsum.fits")[0].data

    # rmat1new = fits.open(new_dir2 + "leakvw.ap90w05/leakrr1.fits")[0].data
    # rmat2new = fits.open(new_dir2 + "leakvw.ap90w05/leakri1.fits")[0].data
    # rmat1new = fits.open(new_dir + "leakfd.ref/leakrr1.fits")[0].data
    # rmat2new = fits.open(new_dir + "leakfd.ref/leakri1.fits")[0].data
    rmat1new = fits.open(new_dir2 + "/leakrr1.fits")[0].data
    rmat2new = fits.open(new_dir2 + "/leakri1.fits")[0].data

    # rmat1new = fits.open(new_dir + "leakfd.ref/leakrr1.fits")[0].data
    # rmat2new = fits.open(new_dir + "leakfd.ref/leakri1.fits")[0].data

    lp = 100
    dm = -14
    ell = 106
    dl = ell - lp

    count = 0

    cc = np.zeros(2*lp+1)
    cj = np.zeros(2*lp+1)
    cc1 = np.zeros(2*lp+1)

    ind0 = int(lp*(lp+1)/2)
    for mp in range(-lp, lp+1):
        ind = ind0 + abs(mp)
        m = mp + dm
        if abs(m) <= ell:
            cj[count] = jesper_data[count, 4]

            # old format
            helpr = rmat1[ind, dl+dl_mat, abs(m)-abs(mp)+dm_mat].copy()
            helpi = rmat2[ind, dl+dl_mat, abs(m)-abs(mp)+dm_mat].copy()
            helpi = -helpi if m*mp < 0 else helpi
            cc[count] = (helpr + helpi)/2

            # new format
            helpr1 = rmat1new[lp, abs(mp), dl+dl_mat, abs(m)-abs(mp)+dm_mat].copy()
            helpi1 = rmat2new[lp, abs(mp), dl+dl_mat, abs(m)-abs(mp)+dm_mat].copy()
            helpi1 = -helpi1 if m*mp < 0 else helpi1
            cc1[count] = (helpr1 + helpi1)/2
            # print(f"{ell:4d} {m:5d} {lp:4d} {mp:5d} {cj[count]:15.5e} {cc[count]:15.5e}" +
                  # f" {cc1[count]:15.5e}")
            count += 1

    fig = plt.figure(figsize=(8, 3))
    plt.subplot(121)
    plt.plot(cc, '-', alpha=0.7, color='black', label='old format')
    plt.plot(cj, '--', alpha=0.5, color='red', label='jesper data')
    plt.legend()

    plt.subplot(122)
    plt.plot(cc-cj, color='black', label='difference')
    plt.legend()
    plt.tight_layout()
    fig.show()

    fig = plt.figure(figsize=(8, 3))
    plt.subplot(121)
    plt.plot(cc1, '-', alpha=0.5, color='black', label='new format')
    plt.plot(cj, '--', alpha=0.5, color='red', label='jesper data')
    plt.legend()

    plt.subplot(122)
    plt.plot(cc1-cj, color='black', label='difference')
    plt.legend()
    plt.tight_layout()
    fig.show()
