anapath = '/home/olivier/GRAND/GP300/ana/GP13/'
galsim = ["VoutRMS2_NSgalaxy.npy","VoutRMS2_EWgalaxy.npy","VoutRMS2_Zgalaxy.npy"]
## FFT
plt.figure(12)
u = "(V$^2$/MHz)" # 
for j in range(3):
    plt.subplot(311+j) # 3*1
    a = np.load(anapath+galsim[j])
    sig_gal18 = a[:,17]/gainlin/gainlin  #Take 18h LST - correct for 20dB gain
    sig_gal6 = a[:,5]/gainlin/gainlin  #Take 6h LST - correct for 20dB gain
    freq_gal = range(10,250) # MHz
    plt.semilogy(freq_gal,sig_gal18[0:240],"--",label="Galaxy - 18hLST (Sandra)")