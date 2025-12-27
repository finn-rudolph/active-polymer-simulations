import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

# -------- user params (edit) ----------
vx_file = "vx_profile.dat"     # output from ave/chunk (time-averaged)
# dump with columns: id x y z vx vy vz c_peratom[4]
dump_file = "dump.stress"
Lx = 50.0                      # box lengths
Ly = 50.0
Lz = 50.0
nbins = 100
gdot_imposed = 0.01            # if you used erate
exclude_edge_frac = 0.1        # fraction of top/bottom bins to exclude from "bulk" fit

timestep = 66000
# --------------------------------------

# --- helper: read vx_profile.dat (adapt if file header differs) ---
# Expect lines like: <bin> <count> <vx> ...
# If ave/chunk file has headings, skip them appropriately. Try to read all floats and find columns.
with open(vx_file) as f:
    for i, line in enumerate(f):
        if line.startswith(str(timestep)):
            data = np.loadtxt(vx_file, skiprows=i+1, max_rows=nbins)
            break

print(data)
# Typical ave/chunk file format: bin_center count vx ...
# Try to find vx as the last column:
bin_centers = data[:, 0]
vx = data[:, -1]

# Exclude edges for bulk fit
n = len(bin_centers)
i0 = int(n*exclude_edge_frac)
i1 = int(n*(1-exclude_edge_frac))
x_bulk = bin_centers[i0:i1]
vx_bulk = vx[i0:i1]

# Linear fit: vx = slope * y + intercept
slope, intercept, r_value, p_value, stderr = stats.linregress(x_bulk, vx_bulk)
r2 = r_value**2
measured_gdot = slope  # because vx = gdot * y (affine)
print(
    f"Measured shear rate (slope) = {measured_gdot:.6e}, R^2 = {r2:.6f}, stderr = {stderr:.3e}")

# Residuals and normalized residuals
vx_fit = slope * bin_centers + intercept
residuals = vx - vx_fit
rel_resid = residuals / (slope * Ly)  # normalized by typical velocity range
print(f"Max abs residual (absolute) = {np.max(np.abs(residuals)):.3e}")
print(
    f"Max rel residual (to total shear-range) = {np.max(np.abs(rel_resid)):.3e}")

# Plot velocity profile and residuals
plt.figure()
plt.subplot(2, 1, 1)
plt.plot(bin_centers, vx, 'o', label='avg vx (bins)')
plt.plot(bin_centers, vx_fit, '-', label=f'linear fit slope={slope:.3e}')
plt.xlabel('y')
plt.ylabel('vx')
plt.legend()
plt.subplot(2, 1, 2)
plt.plot(bin_centers, rel_resid, 'o-')
plt.axhline(0, color='k', linewidth=0.8)
plt.xlabel('y')
plt.ylabel('rel residual (to global shear range)')
plt.tight_layout()
plt.savefig("vx_profile_check.png", dpi=200)
print("Saved vx_profile_check.png")

# --- Process dump.stress to compute sigma_xy(y) per snapshot and time-average ---
# Expect dump with header lines; try to parse generically
# Quick parsing assuming LAMMPS custom dump with columns: id x y z vx vy vz c_peratom[4]
cols = ["id", "x", "y", "z", "vx", "vy", "vz", "sxy"]
# read as whitespace delim, skip comment lines starting with ITEM or loop
# Fast but simple parser (may need adaptation):
frames = []
with open(dump_file, 'r') as f:
    block = []
    for line in f:
        line = line.strip()
        if line == "":
            continue
        # detect start of a snapshot by non-numeric token "ITEM:" or "timestep"
        if line.startswith("ITEM: TIMESTEP") or line.startswith("ITEM: TIMESTEPS") or line.startswith("ITEM:"):
            # ignore header tokens
            block = []
            continue
        parts = line.split()
        # lines with 8 floats -> data
        if len(parts) == len(cols):
            block.append([float(x) for x in parts])
        else:
            # skip unexpected lines
            continue
        # rough heuristic: if block reaches expected #atoms, treat as a frame
        # If you know natoms, set natoms variable above and check length
    # block parsing may fail for different dump formats; if so adapt parser

# If very simple single-frame dump, convert block to array:
if len(block) > 0:
    arr = np.array(block)
    df = pd.DataFrame(arr, columns=cols)
    # Bin by y
    ys = df['y'].values
    sxy = df['sxy'].values
    bins = np.linspace(0, Ly, nbins+1)
    bin_idx = np.digitize(ys, bins) - 1
    sigma_bins = np.zeros(nbins)
    counts = np.zeros(nbins)
    for i in range(nbins):
        mask = (bin_idx == i)
        counts[i] = np.sum(mask)
        if counts[i] > 0:
            # LAMMPS per-atom stress units: be careful; we compute bin-averaged stress
            # Sigma_xy(bin) = - (sum over atoms of sxy_atom) / V_bin
            Vbin = Lx * (bins[i+1]-bins[i]) * Lz
            sigma_bins[i] = - np.sum(sxy[mask]) / Vbin
        else:
            sigma_bins[i] = np.nan

    # plot sigma_xy(y)
    bin_centers_s = 0.5*(bins[:-1]+bins[1:])
    plt.figure()
    plt.plot(bin_centers_s, sigma_bins, 'o-')
    plt.xlabel('y')
    plt.ylabel(r'$\sigma_{xy}(y)$')
    plt.savefig("sigma_xy_profile.png", dpi=200)
    print("Saved sigma_xy_profile.png")

    # Compare sigma uniformity
    finite = np.isfinite(sigma_bins)
    sig_mean = np.nanmean(sigma_bins[finite])
    sig_std = np.nanstd(sigma_bins[finite])
    print(
        f"Mean sigma_xy = {sig_mean:.3e}, std = {sig_std:.3e}, rel std = {sig_std/abs(sig_mean):.3f}")

    # Optionally compute local eta(y) = sigma_xy(y) / shear_rate_local
    # If you assume shear_rate_local ~= measured_gdot, compute eta(y):
    eta_local = sigma_bins / measured_gdot
    plt.figure()
    plt.plot(bin_centers_s, eta_local, 'o-')
    plt.xlabel('y')
    plt.ylabel(r'$\eta(y)$')
    plt.savefig("eta_profile.png")
    print("Saved eta_profile.png")
else:
    print("No dump frames parsed. If your dump format differs, adapt parser. "
          "Alternatively dump in a very simple whitespace format (no ITEM headers) or give natoms to stop per-frame parsing.")
