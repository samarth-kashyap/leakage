# {{{ Library Importsj
from heliosPy import iofuncs as cio
from astropy.io import fits
import numpy as np
import time
# }}} imports

# {{{ Global variables
dl_mat = 6
dm_mat = 15
lmax = 249
# }}} global vars


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


if __name__ == "__main__":
    t1 = time.time()
    data_dir = "/scratch/seismogroup/code/leakvw0/default/"
    write_dir = "/scratch/g.samarth/HMIDATA/leakmat/"
    rmat1 = fits.open(data_dir + "vradsum/leakrlist.vradsum.fits")[0].data
    rmat2 = fits.open(data_dir + "vradsum/leakilist.vradsum.fits")[0].data
    hormat1 = fits.open(data_dir + "vhorsum/leakrlist.vhorsum.fits")[0].data
    hormat2 = fits.open(data_dir + "vhorsum/leakilist.vhorsum.fits")[0].data

    rleaks = np.zeros((2*dm_mat+1, 2*dl_mat+1,
                       lmax+1, 2*lmax + 1), dtype=complex)
    horleaks = np.zeros((2*dm_mat+1, 2*dl_mat+1,
                         lmax+1, 2*lmax + 1), dtype=complex)

    for lp in range(lmax+1):
        ind0 = int(lp*(lp+1)/2)
        for mp in range(lp+1):
            ind = ind0 + mp
            for j in range(2*dl_mat + 1):
                ell = lp + j - dl_mat
                for i in range(2*dm_mat + 1):
                    m = mp + i - dm_mat
                    if abs(m) <= ell:
                        L_ii = 0.5 * sign(1, m) *\
                            rmat2[ind, ell-lp+dl_mat, abs(m)-abs(mp)+dm_mat]
                        L_rr = 0.5 * rmat1[ind, ell-lp+dl_mat,
                                           abs(m)-abs(mp)+dm_mat]
                        rleaks[i, j, lp, mp+lmax] = L_rr + L_ii
                        rleaks[i, j, lp, -mp+lmax] = L_rr + sign(1, mp)*L_ii

                        L_ii = 0.5 * sign(1, m) *\
                            hormat2[ind, ell-lp+dl_mat, abs(m)-abs(mp)+dm_mat]
                        L_rr = 0.5 * hormat1[ind, ell-lp+dl_mat,
                                             abs(m)-abs(mp)+dm_mat]
                        horleaks[i, j, lp, mp+lmax] = L_rr + L_ii
                        horleaks[i, j, lp, -mp+lmax] = L_rr + sign(1, mp)*L_ii
    t2 = time.time()
    print(f" Time taken for computation = {(t2 - t1)/60:.3f} minutes")
    cio.writefitsfile(rleaks.real, write_dir + 'rleaks1.fits')
    cio.writefitsfile(rleaks.imag, write_dir + 'rleaks2.fits')
    cio.writefitsfile(horleaks.real, write_dir + 'horleaks1.fits')
    cio.writefitsfile(horleaks.imag, write_dir + 'horleaks2.fits')
    print("writing complete")
