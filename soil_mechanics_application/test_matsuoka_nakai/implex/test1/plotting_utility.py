#import matplotlib
#import matplotlib.pyplot as plt
#from matplotlib import rc
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
#matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

def plot_deviatoric_strain_vs_deviatoric_stress(e_dev,q_dev):
	plt.plot(e_dev,q_dev)
	plt.grid(color='black', linestyle='--', linewidth=0.05)
	plt.xticks( fontsize=36)
	plt.yticks(fontsize=36)
	plt.xlabel('$\\varepsilon_{q} [-]$',size=36)
	plt.ylabel('$q [Pa] $',size=36)
	plt.show()


def plot_strain_zz_vs_deviatoric_stress(e_zz,q_dev):
	plt.plot(e_zz,q_dev)
	plt.grid(color='black', linestyle='--', linewidth=0.05)
	plt.xticks( fontsize=36)
	plt.yticks(fontsize=36)
	plt.xlabel('$\\varepsilon_{zz} [-]$',size=36)
	plt.ylabel('$q [Pa] $',size=36)
	plt.show()

def plot_strain_axial_vs_deviatoric_stress(e_axial,q_dev):
	plt.plot(e_axial,q_dev)
	plt.grid(color='black', linestyle='--', linewidth=0.05)
	plt.xticks( fontsize=36)
	plt.yticks(fontsize=36)
	plt.xlabel('$\\epsilon_{axial} [-]$',size=36)
	plt.ylabel('$q [Pa] $',size=36)
	plt.show()

def plot_volumetric_strain_vs_deviatoric_stain(e_dev,e_vol):
	plt.plot(e_dev,e_vol)
	plt.grid(color='black', linestyle='--', linewidth=0.05)
	plt.xticks( fontsize=36)
	plt.yticks(fontsize=36)
	plt.xlabel('$\\epsilon_{zz} [-]$',size=36)
	plt.ylabel('$\epsilon_{vol} [-]$',size=36)
	plt.show()

def plot_volumetric_strain_vs_strain_zz(e_zz,e_vol):
	plt.plot(e_zz,e_vol)
	plt.grid(color='black', linestyle='--', linewidth=0.05)
	plt.xticks( fontsize=36)
	plt.yticks(fontsize=36)
	plt.xlabel('$\\epsilon_{zz} [-]$',size=36)
	plt.ylabel('$\epsilon_{vol} [-]$',size=36)
	plt.show()

def plot_volumetric_strain_vs_strain_axial(e_axial,e_vol):
	plt.plot(e_axial,e_vol)
	plt.grid(color='black', linestyle='--', linewidth=0.05)
	plt.xticks( fontsize=36)
	plt.yticks(fontsize=36)
	plt.xlabel('$\\epsilon_{axial} [-]$',size=36)
	plt.ylabel('$\epsilon_{vol} [-]$',size=36)
	plt.show()