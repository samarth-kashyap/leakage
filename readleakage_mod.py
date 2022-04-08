# {{{ Library Imports
from astropy.io import fits
import numpy as np
import time
# }}} imports

# {{{ Global variables
dl_mat = 6
dm_mat = 15
lmax = 300
# }}} global vars


# {{{ def writefitsfile(a, fname, overwrite=True):
def writefitsfile(a, fname, overwrite=True):
    """Writes a fits file with a given filename
    WARNING: If a file exists with the same filename, it will be deleted.

    Parameters:
    -----------
    a - np.ndarray(dtype=float)
        array that needs to be stored in the FITS file
    fname - string
        filename containing full path
    overwrite - bool (default = True)
        overwrites existing fits file if true

    Returns:
    --------
    None
    """
    print("Writing "+fname)
    hdu = fits.PrimaryHDU()
    hdu.data = a
    hdu.writeto(fname, overwrite=overwrite)
    return None
# }}} writefitsfile(a, fname, overwrite=True)


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
# }}} sign(a, b)


# {{{ def minus1pow(m):
def minus1pow(m):
    if (m % 2) == 0:
        return 1
    else:
        return -1
# }}} minus1pow(m)


# {{{ def make_leaks_old():
def make_leaks_old():
    """ Read leakage matrices in the old format """
    rmat1 = fits.open(old_dir + "vradsum/leakrlist.vradsum.fits")[0].data
    rmat2 = fits.open(old_dir + "vradsum/leakilist.vradsum.fits")[0].data
    hormat1 = fits.open(old_dir + "vhorsum/leakrlist.vhorsum.fits")[0].data
    hormat2 = fits.open(old_dir + "vhorsum/leakilist.vhorsum.fits")[0].data

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
    writefitsfile(rleaks.imag, write_dir + 'rleaks2.fits')
    writefitsfile(horleaks.real, write_dir + 'horleaks1.fits')
    writefitsfile(horleaks.imag, write_dir + 'horleaks2.fits')
    print("writing complete")
    return 0
# }}} make_leaks_old()


# {{{ def make_leaks_new():
def make_leaks_new():
    """ Read leakage matrices in the new format """
    if VW:
        # VW leaks
        Lrr = fits.open(new_dir + "leakvw.ap90w05/leakrr1.fits")[0].data
        Lri = fits.open(new_dir + "leakvw.ap90w05/leakri1.fits")[0].data
        Lhr = fits.open(new_dir + "leakvw.ap90w05/leakhr1.fits")[0].data
        Lhi = fits.open(new_dir + "leakvw.ap90w05/leakhi1.fits")[0].data
    else:
        # FD leaks
        Lrr = fits.open(new_dir + "leakfd.ref/leakrr1.fits")[0].data
        Lri = fits.open(new_dir + "leakfd.ref/leakri1.fits")[0].data
        Lhr = fits.open(new_dir + "leakfd.ref/leakhr1.fits")[0].data
        Lhi = fits.open(new_dir + "leakfd.ref/leakhi1.fits")[0].data

    rleaks = np.zeros((2*dm_mat+1, 2*dl_mat+1,
                       lmax+1, 2*lmax + 1), dtype=complex)
    horleaks = np.zeros((2*dm_mat+1, 2*dl_mat+1,
                         lmax+1, 2*lmax + 1), dtype=complex)

    for lp in range(lmax+1):
        print(lp)
    # for lp in range(142, 143):
        for mp in range(-lp, lp+1):
        # for mp in range(-92, -91):#lp+1):
            for j in range(2*dl_mat + 1):
                ell = lp + j - dl_mat
                for i in range(2*dm_mat + 1):
                    m = mp + i - dm_mat
                    if abs(m) <= ell:
                        csphase = 1
                        if m > 0:  # M > 0
                            csphase *= minus1pow(m)
                        if mp > 0:
                            csphase *= minus1pow(mp)
                        pm = sign(1, m*mp)
                        i1 = abs(m) - abs(mp) + dm_mat
                        lr_rr, lr_ii = Lrr[lp, abs(mp), j, i1], Lri[lp, abs(mp), j, i1]
                        lh_rr, lh_ii = Lhr[lp, abs(mp), j, i1], Lhi[lp, abs(mp), j, i1]
                        rleaks[i, j, lp, mp+lmax] = (lr_rr + pm*lr_ii)/2*csphase
                        horleaks[i, j, lp, mp+lmax] = (lh_rr + pm*lh_ii)/2*csphase
                        """
                        print(f"lp = {lp:03d}, mp = {mp:03d}, l = {ell:04d}, m = {m:04d}")
                        print(f"lp = {lp:03d}, abs(mp) = {mp:03d}, l-lp+dl_mat = {j:02d}, abs(m)-abs(mp)+dm_mat = {i1:02d}")
                        print(f"leakrr = {lr_rr}, leakri = {lr_ii}")
                        print(f"leakhr = {lh_rr}, leakhi = {lh_ii}")
                        print(f"csphase = {csphase}, pm={pm}")
                        print(f"rleaks = (leakrr+pm*leakri)/2*csphase = {(lr_rr + pm*lr_ii)/2*csphase}")
                        print(f"horleaks = (leakhr+pm*leakhi)/2*csphase = {(lh_rr + pm*lh_ii)/2*csphase}")
                        print("====================================================================================================")
                        """
    t2 = time.time()
    print(f" Time taken for computation = {(t2 - t1)/60:.3f} minutes")
    if VW:
        writefitsfile(rleaks.real, write_dir + 'rleaks1vw.fits')
        writefitsfile(rleaks.imag, write_dir + 'rleaks2vw.fits')
        writefitsfile(horleaks.real, write_dir + 'horleaks1vw.fits')
        writefitsfile(horleaks.imag, write_dir + 'horleaks2vw.fits')
    else:
        writefitsfile(rleaks.real, write_dir + 'rleaks1fd.fits')
        writefitsfile(rleaks.imag, write_dir + 'rleaks2fd.fits')
        writefitsfile(horleaks.real, write_dir + 'horleaks1fd.fits')
        writefitsfile(horleaks.imag, write_dir + 'horleaks2fd.fits')
    print("writing complete")
    return 0
# }}} make_leaks_new()


# {{{ def make_leaks_new2():
def make_leaks_new2(writenpy=False):
    """ Read leakage matrices in the new format - vectorized version """
    if VW:
        # VW leaks
        Lrr = fits.open(new_dir + "leakvw.ap90w05/leakrr1.fits")[0].data
        Lri = fits.open(new_dir + "leakvw.ap90w05/leakri1.fits")[0].data
        Lhr = fits.open(new_dir + "leakvw.ap90w05/leakhr1.fits")[0].data
        Lhi = fits.open(new_dir + "leakvw.ap90w05/leakhi1.fits")[0].data
    else:
        # FD leaks
        Lrr = fits.open(new_dir + "leakfd.ref/leakrr1.fits")[0].data
        Lri = fits.open(new_dir + "leakfd.ref/leakri1.fits")[0].data
        Lhr = fits.open(new_dir + "leakfd.ref/leakhr1.fits")[0].data
        Lhi = fits.open(new_dir + "leakfd.ref/leakhi1.fits")[0].data

    rleaks = np.zeros((2*dm_mat+1, 2*dl_mat+1,
                       lmax+1, 2*lmax + 1), dtype=complex)
    horleaks = np.zeros((2*dm_mat+1, 2*dl_mat+1,
                         lmax+1, 2*lmax + 1), dtype=complex)

    t1 = time.time()
    for lp in range(lmax+1):
        ell = np.arange(lp-dl_mat,lp+dl_mat+1)
        for mp in range(-lp, lp+1):
            for i in range(2*dm_mat + 1):
                m = mp + i - dm_mat
                irrelevant_ind = np.ones(ell.shape)
                irrelevant_ind[abs(m)>ell] = 0
                csphase = 1
                if m > 0:
                    csphase *= minus1pow(m)
                if mp > 0:
                    csphase *= minus1pow(mp)
                pm = sign(1, m*mp)

                i1 = abs(m) - abs(mp) + dm_mat
                lr_rr, lr_ii = Lrr[lp, abs(mp), :, i1], Lri[lp, abs(mp), :, i1]
                lh_rr, lh_ii = Lhr[lp, abs(mp), :, i1], Lhi[lp, abs(mp), :, i1]
                rleaks[i, :, lp, mp+lmax] = (lr_rr + pm*lr_ii)/2*csphase
                horleaks[i, :, lp, mp+lmax] = (lh_rr + pm*lh_ii)/2*csphase
                rleaks[i, :, lp, mp+lmax] *= irrelevant_ind
                horleaks[i, :, lp, mp+lmax] *= irrelevant_ind
    t2 = time.time()
    print(f" Time taken for computation = {(t2 - t1)/60:.3f} minutes")

    if writenpy:
        if VW:
            np.save(write_dir + 'rleaks1vw.npy', rleaks.real)
            np.save(write_dir + 'horleaks1vw.npy', horleaks.real)
        else:
            np.save(write_dir + 'rleaks1fd.npy', rleaks.real)
            np.save(write_dir + 'horleaks1fd.npy', horleaks.real)
        print("writing complete")
        return 0
    else:
        if VW:
            writefitsfile(rleaks.real, write_dir + 'rleaks1vw.fits')
            writefitsfile(rleaks.imag, write_dir + 'rleaks2vw.fits')
            writefitsfile(horleaks.real, write_dir + 'horleaks1vw.fits')
            writefitsfile(horleaks.imag, write_dir + 'horleaks2vw.fits')
        else:
            writefitsfile(rleaks.real, write_dir + 'rleaks1fd.fits')
            writefitsfile(rleaks.imag, write_dir + 'rleaks2fd.fits')
            writefitsfile(horleaks.real, write_dir + 'horleaks1fd.fits')
            writefitsfile(horleaks.imag, write_dir + 'horleaks2fd.fits')
        print("writing complete")
        return 0
# }}} make_leaks_new2()


if __name__ == "__main__":
    t1 = time.time()
    old_dir = "/scratch/seismogroup/code/leakvw0/default/"
    new_dir = "/scratch/seismogroup/leakage/HMI/"
    # write2dir = "/scratch/g.samarth/HMIDATA/leakmat/"
    write_dir = "/scratch/g.samarth/HMIDATA/leakmat_new/"

    VW = False

    make_leaks_new2(writenpy=True)
