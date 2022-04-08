from astropy.io import fits
import numpy as np
import time

dl_mat = 6
dm_mat = 15
lmax = 249


def writefitsfile(a, fname):
    """Writes a fits file with a given filename

    Parameters:
    -----------
    a - np.ndarray(dtype=float)
        array that needs to be stored in the FITS file
    fname - string
        filename containing full path

    Returns:
    --------
    None
    """
    print("Writing "+fname)
    hdu = fits.PrimaryHDU()
    hdu.data = a
    hdu.writeto(fname, overwrite=True)
    return None


def sign(a, b):
    """ Returns value of a with sign of b """
    return abs(a) if b >= 0 else -abs(a)


def minus1pow(m):
    """ Returns (-1)^m """
    return 1 if m % 2 == 0 else -1


# This is given as a reference as we've already established the correctness
# of this using some of the element values you sent previously
def make_leaks_VW():
    rmat1 = fits.open(vw_dir + "vradsum/leakrlist.vradsum.fits")[0].data
    rmat2 = fits.open(vw_dir + "vradsum/leakilist.vradsum.fits")[0].data
    hormat1 = fits.open(vw_dir + "vhorsum/leakrlist.vhorsum.fits")[0].data
    hormat2 = fits.open(vw_dir + "vhorsum/leakilist.vhorsum.fits")[0].data

    rleaks = np.zeros((2*dm_mat+1, 2*dl_mat+1,
                       lmax+1, 2*lmax + 1), dtype=complex)
    horleaks = np.zeros((2*dm_mat+1, 2*dl_mat+1,
                         lmax+1, 2*lmax + 1), dtype=complex)

    for lp in range(lmax+1):
        ind0 = int(lp*(lp+1)/2)
        for mp in range(-lp, lp+1):
            ind = ind0 + abs(mp)
            for j in range(2*dl_mat + 1):
                ell = lp + j - dl_mat
                for i in range(2*dm_mat + 1):
                    m = mp + i - dm_mat
                    if abs(m) <= ell:
                        csphase = 1
                        if m > 0:
                            csphase *= minus1pow(m)
                        if mp > 0:
                            csphase *= minus1pow(mp)
                        pm = sign(1, m*mp)
                        Lrr = rmat1[ind, j, abs(m)-abs(mp)+dm_mat]*csphase
                        Lii = rmat2[ind, j, abs(m)-abs(mp)+dm_mat]*csphase
                        rleaks[i, j, lp, mp+lmax] = (Lrr + pm*Lii)/2

                        Lrr = hormat1[ind, j, abs(m)-abs(mp)+dm_mat]*csphase
                        Lii = hormat2[ind, j, abs(m)-abs(mp)+dm_mat]*csphase
                        horleaks[i, j, lp, mp+lmax] = (Lrr + pm*Lii)/2
    t2 = time.time()
    print(f" Time taken for computation = {(t2 - t1)/60:.3f} minutes")
    writefitsfile(rleaks.real, write_dir + 'rleaks1.fits')
    writefitsfile(horleaks.real, write_dir + 'horleaks1.fits')
    print("writing complete")
    return 0


def make_leaks_FD():
    Lrr = fits.open(fd_dir + "leakfd.ref/leakrr1.fits")[0].data
    Lri = fits.open(fd_dir + "leakfd.ref/leakri1.fits")[0].data
    Lhr = fits.open(fd_dir + "leakfd.ref/leakhr1.fits")[0].data
    Lhi = fits.open(fd_dir + "leakfd.ref/leakhi1.fits")[0].data

    rleaks = np.zeros((2*dm_mat+1, 2*dl_mat+1,
                       lmax+1, 2*lmax + 1), dtype=np.float)
    horleaks = np.zeros((2*dm_mat+1, 2*dl_mat+1,
                         lmax+1, 2*lmax + 1), dtype=np.float)

    for lp in range(lmax+1):
        for mp in range(-lp, lp+1):
            for j in range(2*dl_mat + 1):
                ell = lp + j - dl_mat
                for i in range(2*dm_mat + 1):
                    m = mp + i - dm_mat
                    if abs(m) <= ell:
                        csphase = 1
                        if m > 0:
                            csphase *= minus1pow(m)
                        if mp > 0:
                            csphase *= minus1pow(mp)
                        pm = sign(1, m*mp)
                        lr_rr = Lrr[lp, abs(mp), j, i]
                        lr_ii = Lri[lp, abs(mp), j, i]
                        lh_rr = Lhr[lp, abs(mp), j, i]
                        lh_ii = Lhi[lp, abs(mp), j, i]
                        rleaks[i, j, lp, mp+lmax] = (lr_rr + pm*lr_ii) * csphase / 2
                        horleaks[i, j, lp, mp+lmax] = (lh_rr + pm*lh_ii) * csphase / 2
    t2 = time.time()
    print(f" Time taken for computation = {(t2 - t1)/60:.3f} minutes")
    writefitsfile(rleaks, write_dir + 'rleaks1fd.fits')
    writefitsfile(horleaks, write_dir + 'horleaks1fd.fits')
    print("writing complete")
    return 0


if __name__ == "__main__":
    t1 = time.time()
    vw_dir = "/scratch/seismogroup/code/leakvw0/default/"
    fd_dir = "/scratch/g.samarth/leakage/HMI/"
    write_dir = "/scratch/g.samarth/HMIDATA/leakmat/"
    make_leaks_FD()

