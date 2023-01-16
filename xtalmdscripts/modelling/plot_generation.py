import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# This script will generate plots of results for different OpenFF result comparision
off_box = pd.read_csv('minimization_results.csv')

# convert RMSD values to numpy arrays
off100_box_rmsd = off_box[off_box['COD ID'].str.endswith("1.0.0")]['RMSD20 Mean'].to_numpy(dtype='float64')
off110_box_rmsd = off_box[off_box['COD ID'].str.endswith("1.1.0")]['RMSD20 Mean'].to_numpy(dtype='float64')
off120_box_rmsd = off_box[off_box['COD ID'].str.endswith("1.2.0")]['RMSD20 Mean'].to_numpy(dtype='float64')
off130_box_rmsd = off_box[off_box['COD ID'].str.endswith("1.3.0")]['RMSD20 Mean'].to_numpy(dtype='float64')
off200_box_rmsd = off_box[off_box['COD ID'].str.endswith("2.0.0")]['RMSD20 Mean'].to_numpy(dtype='float64')

# Axis setting
fig, axs = plt.subplots(2, 1)
axs[0].hist((off100_box_rmsd, off110_box_rmsd, off120_box_rmsd, off130_box_rmsd, off200_box_rmsd))
axs[0].set_xlabel('RMSD20')
axs[0].set_ylabel('Frequency')
axs[0].legend(['OpenFF 1.0.0', 'OpenFF 1.1.0', 'OpenFF 1.2.0', 'OpenFF 1.3.0', 'OpenFF 2.0.0'])
axs[0].set_title('Full Box Minimization')
axs[1].hist((off100_box_rmsd, off200_box_rmsd))
axs[1].set_xlabel('RMSD20')
axs[1].set_ylabel('Frequency')
axs[1].legend(['OpenFF 1.0.0', 'OpenFF 2.0.0'])
axs[1].set_title('Full Box Minimization')
fig.tight_layout(h_pad=3)
plt.savefig('RMSD_comparison.png')

# Mean value calculation
print('RMSD20 (box min) OFF 1.0.0 Mean: %s nm' % np.mean(off100_box_rmsd))
print('RMSD20 (box min) OFF 1.1.0 Mean: %s nm' % np.mean(off110_box_rmsd))
print('RMSD20 (box min) OFF 1.2.0 Mean: %s nm' % np.mean(off120_box_rmsd))
print('RMSD20 (box min) OFF 1.3.0 Mean: %s nm' % np.mean(off130_box_rmsd))
print('RMSD20 (box min) OFF 2.0.0 Mean: %s nm' % np.mean(off200_box_rmsd))
plt.show()
