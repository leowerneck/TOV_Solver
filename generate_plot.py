# .-----------------------------------------.
# |  Copyright (c) 2019, Leonardo Werneck   |
# | Licensed under the BSD 2-Clause License |
# .-----------------------------------------.

from numpy import loadtxt, linspace
import matplotlib.pyplot as plt

# Names of the EOSs used
EOSnames = ["WFF1","WFF2","APR4","SLy","ENG","APR3","MPA1","ALF2","H4","MS1b","MS1"]

# Colors to be used in the plot
EOScolors = ["Blue","Red","Green","Black","Orange","Magenta"]

# Linestyles
EOSlines = ["-","-","-","-","-","-","-.","-.","-.","-.","-."]

fig = plt.figure(dpi=300)
ax = fig.add_subplot(111)
ax.set_color_cycle(EOScolors)
ax.set_xlim(9,17)
ax.set_ylim(0,3)
ax.set_xlabel(r"$R$ (km)")
ax.set_ylabel(r"$M$ (solar masses)")
ax.minorticks_on()
ax.grid(which='major',ls='-',lw=0.5)

for EOS in range(len(EOSnames)):
    filename = EOSnames[EOS]+"_mass_vs_radius.dat"
    dummy, dummy, M, R = loadtxt(filename).T
    ax.plot(R,M,lw='2',label=EOSnames[EOS],ls=EOSlines[EOS])

ax.legend()

plt.savefig("Mass_vs_Radius.png",dpi=300)
plt.close(fig)
