import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import chisquare
from scipy.optimize import curve_fit

# Assuming these values based on context
C0 = 3e8  # Speed of light in meters/second
RND = 0
dis = True
datetime = 20231206021552

# Read  data
trigger_time_data = pd.read_csv(f'D:\\STUDY\\GRAND\\result\\reconstruction\\reconstruction_exact\\trigger_{datetime}_exact.txt', delimiter='\t')
position_data = pd.read_csv(f'D:\\STUDY\\GRAND\\result\\reconstruction\\reconstruction_exact\\position_{datetime}_exact.txt', delimiter='\t')

# Merge dataframes on 'event_id' and 'DU_number'
merged_data = pd.merge(trigger_time_data, position_data, on=['event_id', 'DU_number'])

# Sort by trigger time
sorted_data = merged_data.sort_values(by='trigger_time')

# Extract relevant columns
detId = sorted_data['DU_number']
detPos = sorted_data[['position_x', 'position_y', 'position_z']].values
expdelays = sorted_data['trigger_time'].apply(lambda x: x.split('+')[1].split('*')[0]).astype(float)

detId = detId.tolist()
detPos = detPos.tolist()
expdelays = expdelays.tolist()

# reconstruction results
x = np.array([1.57390e+05])  
y = np.array([-4.74582e+05])  
z = np.array([-3.14097e+04]) 

DetectorType = np.zeros(len(detId))
thetap = np.array([8.61200e+01])  
phip = np.array([1.08335e+02])

FSAMPLING = np.array([5e8])  
Lant = np.array([4])
ind = 0  
ErrorTrig = 1

# Plane Reconstruction
k = np.array([np.cos(np.radians(phip[ind])) * np.sin(np.radians(thetap[ind])),
              np.sin(np.radians(phip[ind])) * np.sin(np.radians(thetap[ind])),
              np.cos(np.radians(thetap[ind]))]).T
plandelays = np.column_stack([detId, np.dot(detPos, k)])
plandelays[:, 1] -= np.min(plandelays[:, 1])
plandelays[:, 1] /= C0 / 1e9  # in ns

# Spherical Reconstruction
Xs = np.column_stack([x[ind], y[ind], z[ind]])
sphdelays = np.column_stack([detId, np.sqrt(np.sum((detPos - np.ones((len(detId), 1)) * Xs) ** 2, axis=1))])
sphdelays[:, 1] -= np.min(sphdelays[:, 1])
sphdelays[:, 1] /= C0 / 1e9  # in ns

# set the minimum time to 0
expdelays -= np.min(expdelays)

# Select relevant data based on DetectorType
isAnt = np.where(DetectorType == 0)[0]
isSci = np.where(DetectorType == 1)[0]
plandelays_ant = plandelays[isAnt, :]
plandelays_sci = plandelays[isSci, :]
sphdelays_ant = sphdelays[isAnt, :]
sphdelays_sci = sphdelays[isSci, :]
expdelays_ant = expdelays[isAnt]
expdelays_sci = expdelays[isSci]

# linear_fit
def linear_fit(x, slope, intercept):
    return slope * x + intercept

errort = 1 * ErrorTrig / FSAMPLING * 1e9  # ns

# Calculate chi2 for plane reconstruction # not used
chi2_plan, p_plan = chisquare((expdelays_ant - plandelays_ant[:, 1]) ** 2 / errort ** 2)
chi2_sph, p_sph = chisquare((expdelays_ant - sphdelays_ant[:,1]) **2 / errort ** 2)

# Fit a line to the scatter plot
fit_params_plan, _ = curve_fit(linear_fit, plandelays_ant[:, 1], expdelays_ant, sigma=np.ones(len(isAnt)) * errort)
fit_params_sph, _ = curve_fit(linear_fit, sphdelays_ant[:, 1], expdelays_ant, sigma=np.ones(len(isAnt)) * errort)

# chi2 analysis
# Plane reconstruction
LinParametersPlan, _ = curve_fit(lambda x, a: a * x, plandelays_ant[:, 1], expdelays_ant)
slope_plan = LinParametersPlan[0]
khi2_plan = np.sum((expdelays_ant - slope_plan * plandelays_ant[:, 1])**2 / errort**2)
n_plan = len(expdelays_ant) - 2
khi2n_plan = khi2_plan / n_plan

# Spherical reconstruction
LinParametersSph, _ = curve_fit(lambda x, a: a * x, sphdelays_ant[:, 1], expdelays_ant)
slope_sph = LinParametersSph[0]
khi2_sph = np.sum((expdelays_ant - slope_sph * sphdelays_ant[:, 1])**2 / errort**2)
n_sph = len(expdelays_ant) - 2
khi2n_sph = khi2_sph / n_sph

# Plotting
fig, ax = plt.subplots(1, 2, figsize=(12, 6), dpi = 100)


print(plandelays_ant[:, 1])
print(expdelays_ant)
print(np.ones(len(isAnt)) * errort,end='\n\n')

print(sphdelays_ant[:,1])
print(expdelays_ant)
print(np.ones(len(isAnt)) * errort)



ax[0].scatter(plandelays_ant[:, 1], expdelays_ant, c='blue', alpha=0.7)
ax[0].plot(plandelays_ant[:, 1], linear_fit(plandelays_ant[:, 1], *fit_params_plan), color='red')
ax[0].set_ylabel('Measured Delays (ns)')
ax[0].set_xlabel('Expected Delays (ns)')
ax[0].set_title('Plane Reconstruction')
ax[0].annotate(f'Theta: {np.mean(thetap):.2f}\nPhi: {np.mean(phip):.2f}\nChi2/ndf: {khi2n_plan:.2f}', xy=(0.15, 0.9), xycoords='axes fraction', ha='center')
for det_id, txt in zip(plandelays_ant[:, 0], isAnt):
    ax[0].annotate(int(det_id), (plandelays_ant[txt, 1], expdelays_ant[txt]),fontsize=12)
  

ax[1].scatter(sphdelays_ant[:, 1], expdelays_ant)
ax[1].plot(sphdelays_ant[:, 1], linear_fit(sphdelays_ant[:, 1], *fit_params_sph), color='red', label='Fit Line')
ax[1].set_ylabel('Measured Delays (ns)')
ax[1].set_xlabel('Expected Delays (ns)')
ax[1].set_title('Spherical Reconstruction')
ax[1].annotate(f'Theta: {np.mean(thetap):.2f}\nPhi: {np.mean(phip):.2f}\nChi2/ndf: {khi2n_sph:.2f}', xy=(0.15, 0.9), xycoords='axes fraction', ha='center')
for det_id, txt in zip(sphdelays_ant[:, 0], isAnt):
    ax[1].annotate(int(det_id), (sphdelays_ant[txt, 1], expdelays_ant[txt]),fontsize=12)

plt.suptitle(f'Expected Delays vs Measured Delays  {datetime}')
#plt.show()
plt.savefig(f'D:\\STUDY\\GRAND\\plot\\reconstruction\\{datetime}.png')

