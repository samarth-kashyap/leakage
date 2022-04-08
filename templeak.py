from astropy.io import fits

dl_mat = 6
dm_mat = 15
def print_leakmat(l, lp, m, mp):
    lrr, lri = leakmatr[lp, abs(mp), l-lp+dl_mat, abs(m)-abs(mp)+dm_mat], leakmati[lp, abs(mp), l-lp+dl_mat, abs(m)-abs(mp)+dm_mat]
    print(f"l = {l:03d}; lp = {lp:03d}; m = {m:03d}; mp = {mp:03d}; leak_r = {lrr}, leak_i = {lri}")
    return None

leakmatr = fits.open("/scratch/seismogroup/leakage/HMI/leakvw.ap90w05/leakrr1.fits")[0].data
leakmati = fits.open("/scratch/seismogroup/leakage/HMI/leakvw.ap90w05/leakri1.fits")[0].data
print_leakmat(140, 146, -140, -136)
