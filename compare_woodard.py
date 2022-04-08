# {{{ Library imports
import numpy as np
from math import pi
from astropy.io import fits
from heliosPy import datafuncs as cdata
import matplotlib.pyplot as plt
# }}} imports

# {{{ global vars
lmax = 249
dl, dm = 4, 4
ell, emm = 200, 120
dl_mat, dm_mat = 6, 15
rsun = 6.9598e10
twopiemin6 = 2*pi*1e-6
leak_dir = "/scratch/g.samarth/HMIDATA/leakmat/"

plot_cmp = True
VWleaks = False
leak_dict = {"VW": False,
             "FD": True,
             "WD": False,
             "JS": False}
# }}} global vars


def get_leaks():
    """Get both leakage matrices.

    Inputs:
    -------
    None

    Outputs:
    --------
    """
    wl = np.loadtxt("leak_woodard_py.dat")
    jlmix = np.zeros_like(wl, dtype=np.float)

    if leak_dict['VW']:
        rleaks = fits.open(f"{leak_dir}rleaks1.fits")[0].data
        horleaks = fits.open(f"{leak_dir}horleaks1.fits")[0].data
        radlk = rleaks[:, :, ell, emm+lmax]
        horlk = horleaks[:, :, ell, emm+lmax]
        del rleaks, horleaks

    if leak_dict['FD']:
        rleaks = fits.open(f"{leak_dir}rleaks1fd.fits")[0].data
        horleaks = fits.open(f"{leak_dir}horleaks1fd.fits")[0].data
        radlk = rleaks[:, :, ell, emm+lmax]
        horlk = horleaks[:, :, ell, emm+lmax]
        del rleaks, horleaks

    count = 0
    for dell in range(dl, -dl-1, -1):
        l1 = ell - dell
        for dem in range(dm, -dm-1, -1):
            m1 = emm - dem
            if abs(m1) <= l1:
                omeganl1, fwhmnl1, amp1 = cdata.findfreq(mode_data,
                                                         l1, 0, m1)
                mix = (274.8*1e2*(l1+0.5) /
                    ((twopiemin6)**2*rsun))/omeganl1**2
                compensate_CS = 1
                if (m1 > 0) and (m1 % 2 == 1):
                    compensate_CS *= -1
                if leak_dict["FD"]:
                    compensate_CS = 1
                r1 = radlk[dem+dm_mat, dell+dl_mat]
                h1 = horlk[dem+dm_mat, dell+dl_mat]
                jlmix[count, 0] = l1
                jlmix[count, 1] = m1
                jlmix[count, 2] = compensate_CS*(r1 + mix*h1)
                count += 1
    return jlmix, wl


if __name__ == "__main__":
    mode_data = np.loadtxt(f"{leak_dir}hmi.6328.36")

    # if leak_dict["VW"]:
    #     rleaks = fits.open(f"{leak_dir}rleaks1.fits")[0].data
    #     horleaks = fits.open(f"{leak_dir}horleaks1.fits")[0].data
    #     radlk = rleaks[:, :, ell, emm+lmax]
    #     horlk = horleaks[:, :, ell, emm+lmax]
    #     del rleaks, horleaks

    # if leak_dict["FD"]:
    #     rleaks = fits.open(f"{leak_dir}rleaks1fd.fits")[0].data
    #     horleaks = fits.open(f"{leak_dir}horleaks1fd.fits")[0].data
    #     radlk = rleaks[:, :, ell, emm+lmax]
    #     horlk = horleaks[:, :, ell, emm+lmax]
    #     del rleaks, horleaks

    jlmix, wl = get_leaks()

    print(f" sum(abs(diff_ll)) = {abs(jlmix[:, 0] - wl[:, 0]).sum()}")
    print(f" sum(abs(diff_mm)) = {abs(jlmix[:, 1] - wl[:, 1]).sum()}")
    print(f" sum(abs(leak_r)) = {abs(jlmix[:, 2] - wl[:, 2]).sum()}")
    print(f" sum(abs(leak_mix)) = {abs(jlmix[:, 2] - wl[:, 2]).sum()}")

    norm = abs(jlmix[:, 2]).max()/abs(wl[:, 2]).max()
    print(f"norm = {norm}")

    if plot_cmp:
        plt.figure(figsize=(10, 10))
    #    plt.plot(jl[:, 2], 'b', alpha=0.6, label="Jesper (r)")
        plt.plot(jlmix[:, 2]*0.3, '--', color='black',
                alpha=0.8, label="Jesper (r + mix*h)")
        plt.plot(wl[:, 2], 'r', alpha=0.8, label="Woodard")
    #    plt.plot(jlmix[:, 2] + wl[:, 2], '-.b', label="Difference")
        plt.legend()
        plt.show()
