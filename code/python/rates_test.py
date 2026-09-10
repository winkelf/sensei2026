import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from scipy.stats import chi2
plt.style.use('dark_background')

def upperLimit(obs_number_input):
    CL = 0.9
    p_on_the_left = 1-CL
    degrees_of_freedom = 2* (obs_number_input+1)

    # Calculate the x = double of mu
    x = chi2.ppf(1 - p_on_the_left, degrees_of_freedom)
    mu_output = x / 2
    #return mu_output
    return obs_number_input 

# Load all data
filenames = [
    ("../cxx/new_test_results/txts/more_survElec.txt",    "../cxx/new_test_results/txts/LEC_survElec.txt", 'Number of unmasked n-electron events vs q cut', "Number of unmasked n-electron events"),
    ("../cxx/new_test_results/txts/more_survPix.txt",     "../cxx/new_test_results/txts/LEC_survPix.txt",  "Number of unmasked pixels vs q cut",            "Number of unmasked pixels"),
    ("../cxx/new_test_results/txts/more_rates.txt",       "../cxx/new_test_results/txts/LEC_rates.txt",    "Rates vs q cut",                                 "n-Electron rate")
]

data = []
for f, fLEC, _, _ in filenames: data.append((np.loadtxt(f), np.loadtxt(fLEC)))

# Setup the layout
fig = plt.figure(figsize=(14, 8))
gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])
ax0 = fig.add_subplot(gs[0, 0])  # Top left
ax1 = fig.add_subplot(gs[1, 0])  # Bottom left
ax2 = fig.add_subplot(gs[:, 1])  # Right (spans both rows)

# Plot 1: survElec
x, y1, y2, y3       = data[0][0][:, 0], data[0][0][:, 1], data[0][0][:, 2], data[0][0][:, 3]
y1LEC, y2LEC, y3LEC = data[0][1][:, 1], data[0][1][:, 2], data[0][1][:, 3]

y1_UL    = np.array([upperLimit(val) for val in y1])
y2_UL    = np.array([upperLimit(val) for val in y2])
y3_UL    = np.array([upperLimit(val) for val in y3])
y1LEC_UL = np.array([upperLimit(val) for val in y1LEC])
y2LEC_UL = np.array([upperLimit(val) for val in y2LEC])
y3LEC_UL = np.array([upperLimit(val) for val in y3LEC])

#for val in y2: print(val/upperLimit(val))

ax0.plot(x, y1_UL,   "-o", label="1 electron count", color="blue")
ax0.plot(x, y2_UL,   "-o", label="2 electrons count", color="green")
ax0.plot(x, y3_UL,   "-o", label="3 electrons count", color="orange")
ax0.plot(x, y1LEC_UL, "-", label="1 electron count (LEC)", color="yellow")
ax0.plot(x, y2LEC_UL, "-", label="2 electrons count (LEC)", color="red")
ax0.plot(x, y3LEC_UL, "-", label="3 electrons count (LEC)", color="magenta")
#ax0.set_title("Upper limit on number of surviving n-electron events vs cut")
ax0.set_title("Number of surviving n-electron events vs cut")
ax0.set_xlabel('Cut in $k/N$')
#ax0.set_ylabel("Upper limit on number of surviving n-electron events")
ax0.set_ylabel("Number of surviving n-electron events")
ax0.grid(True, linestyle=':')
ax0.legend()

# Plot 2: survPix
x_pix, y1_pix, y2_pix, y3_pix = data[1][0][:, 0], data[1][0][:, 1], data[1][0][:, 2], data[1][0][:, 3]
y1LEC_pix, y2LEC_pix, y3LEC_pix = data[1][1][:, 1], data[1][1][:, 2], data[1][1][:, 3]
ax1.plot(x_pix, y1_pix, "-o", label="1 electron pixels", color="blue")
ax1.plot(x_pix, y2_pix, "-o", label="2 electrons pixels", color="green")
ax1.plot(x_pix, y3_pix, "-o", label="3 electrons pixels", color="orange")
ax1.plot(x_pix, y1LEC_pix, "-", label="1 electron pixels (LEC)", color="yellow")
ax1.plot(x_pix, y2LEC_pix, "-", label="2 electrons pixels (LEC)", color="red")
ax1.plot(x_pix, y3LEC_pix, "-", label="3 electrons pixels (LEC)", color="magenta")
ax1.set_title(filenames[1][2])
#ax1.set_xlabel('Cut in $q=-2\\ln(\\lambda)$')
ax1.set_xlabel('Cut in $k/N$')
ax1.set_ylabel(filenames[1][3])
ax1.grid(True, linestyle=':')
ax1.legend()

# Plot 3: rates (computed as survElec / survPix)
# Avoid division by zero with np.where or masking
rates1 = np.divide(y1, y1_pix, out=np.zeros_like(y1), where=y1_pix != 0)
rates2 = np.divide(y2, y2_pix, out=np.zeros_like(y2), where=y2_pix != 0)
rates3 = np.divide(y3, y3_pix, out=np.zeros_like(y3), where=y3_pix != 0)
rates1LEC = np.divide(y1LEC, y1LEC_pix, out=np.zeros_like(y1LEC), where=y1LEC_pix != 0)
rates2LEC = np.divide(y2LEC, y2LEC_pix, out=np.zeros_like(y2LEC), where=y2LEC_pix != 0)
rates3LEC = np.divide(y3LEC, y3LEC_pix, out=np.zeros_like(y3LEC), where=y3LEC_pix != 0)

## Desired x value
#x_val = 15
#
## Find index of the given x
#idx = np.where(x == x_val)[0]
#
#if len(idx) > 0:
#    print(f"x = {x_val}")
#    print(f"rates1 = {rates1[idx[0]]}")
#    print(f"rates1LEC = {rates1LEC[idx[0]]}")
#    print(f"Rates ratio = {rates1[idx[0]]/rates1LEC[idx[0]]}")
#else:
#    print("x value not found in the array")

ax2.plot(x, rates1, "-o", label="1 electron rate", color="blue")
ax2.plot(x, rates2, "-o", label="2 electrons rate", color="green")
ax2.plot(x, rates3, "-o", label="3 electrons rate", color="orange")
ax2.plot(x, rates1LEC, "-", label="1 electron rate (LEC)", color="yellow")
ax2.plot(x, rates2LEC, "-", label="2 electrons rate (LEC)", color="red")
ax2.plot(x, rates3LEC, "-", label="3 electrons rate (LEC)", color="magenta")
ax2.set_title(filenames[2][2])
#ax2.set_xlabel('Cut in $q=-2\\ln(\\lambda)$')
ax2.set_xlabel('Cut in $k/N$')
ax2.set_ylabel("Ratio")
#ax2.set_xlim([0,30]) 
#ax2.set_xlim([10,60]) 
#ax2.set_ylim([0,0.01]) 
#ax2.set_ylim([0,0.015])
ax2.set_yscale('symlog', linthresh=1e-7)
ax2.set_ylim(bottom=0)
ax2.grid(True, linestyle=':')
ax2.legend()

pdfName = "../cxx/new_test_results/rates_test.pdf"

plt.tight_layout()
#plt.savefig("summary_plots.pdf")
#print("Generated summary_plots.pdf")
plt.savefig(pdfName)
print(f"Generated {pdfName}")
