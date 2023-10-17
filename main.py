'''
    Arya Koureshi
    MSc student at Sharif University of Technology
    Researcher at Sharif Brain Center
    Researcher at GhazizadehLab
    
    arya.koureshi@gmail.com
'''

# %% Imports
import numpy as np
import matplotlib.pyplot as plt

# %% Model I
# rate constants k1 and ke
K1 = [0.9776, 0.7448, 0.3293, 0.2213]
Ke = [0.05]

# initial drug concentration
c0 = 500

t = np.linspace(0, 6, 1000)

plt.figure(figsize=(10, 10))
for k1 in K1:
    for ke in Ke:
        # Calculate the drug concentration in the GI tract and blood stream compartments
        c1 = c0 * np.exp(-k1 * t)
        c2 = (c0 * k1 / (k1 - ke)) * (np.exp(-ke * t) - np.exp(-k1 * t))
        plt.plot(t, c1, label=str(k1)+" / hr")
plt.legend()
plt.ylim([0, c0])
plt.xlim([0, t[-1]])
plt.xlabel('Time (hours)')
plt.ylabel('Concentration (units)')
plt.title('Drug concentration pattern in single compartment model')
plt.show()
plt.savefig("Model1_1.png", dpi=312)

plt.figure(figsize=(10, 10))
ke = 0.05
k1 = 0.9776
C0 = [600, 250, 200, 100]
for c0 in C0:
        # Calculate the drug concentration in the GI tract and blood stream compartments
        c1 = c0 * np.exp(-k1 * t)
        c2 = (c0 * k1 / (k1 - ke)) * (np.exp(-ke * t) - np.exp(-k1 * t))
        plt.plot(t, c1, label=r'$c_0$ = ' + str(c0))
plt.legend()
plt.ylim([0, C0[0]])
plt.xlim([0, t[-1]])
plt.xlabel('Time (hours)')
plt.ylabel('Concentration (units)')
plt.title(r'Variation of drug concentration in single compartment with $k_1 = 9776 / hr$')
plt.show()
plt.savefig("Model1_2.png", dpi=312)

plt.figure(figsize=(10, 10))
Ke = [0.2213, 0.9776, 0.3293, 0.5228, 0.7448]
K1 = [0.9776]
c0 = 500
for k1 in K1:
    for ke in Ke:
        # Calculate the drug concentration in the GI tract and blood stream compartments
        c1 = c0 * np.exp(-k1 * t)
        if ke == 0.9776:
            plt.plot(t, c1, label=str(k1)+" / hr")
        else:
            if ke == 0.2213:
                plt.text(1.2, 350, "Drug concentration in bloodstream")
            c2 = (c0 * k1 / (k1 - ke)) * (np.exp(-ke * t) - np.exp(-k1 * t))
            plt.plot(t, c2, label=str(ke)+" / hr")
plt.legend()
plt.ylim([0, c0])
plt.xlim([0, t[-1]])
plt.xlabel('Time (hours)')
plt.ylabel('Concentration (units)')
plt.title('Drug distribution pattern in oral administration - Model-I')
plt.show()
plt.savefig("Model1_3.png", dpi=312)

# %% Model II
# case 1
# rate constants kb, ke, and kt
kb = 0.9776
ke = 0.2213
kt = 0.3293

# initial drug concentration
c0 = 500

t = np.linspace(0, 2.5, 1000)

# Calculate the roots of the equation
a = 1
b = kb + ke + kt
c = ke * kt
D = np.sqrt(b**2 - 4*a*c)
u1 = -1 * (-b + D)/(2*a)
u2 = -1 * (-b - D)/(2*a)

# Calculate the drug concentration in the blood and tissue compartments
cb = (c0 / (u2 - u1)) * ((-u1 + kt) * np.exp(-u1 * t) - (-u2 + kt) * np.exp(-u2 * t))
ct = (c0 * kb / (u1 - u2)) * (np.exp(-u1 * t) - np.exp(-u2 * t))

plt.figure(figsize=(10, 10))
plt.plot(t, cb, label='Blood')
plt.plot(t, ct, label='Tissue')
plt.xlim([0, t[-1]])
plt.xlabel('Time (hours)')
plt.ylabel('Concentration (units)')
plt.suptitle('Drug concentration pattern in Model-II - case 1')
plt.title(r'$k_b$ = {} , $k_e$ = {} , $k_t$ = {}'.format(kb, ke, kt))
plt.legend()
plt.show()
plt.savefig("Model2_case1.png", dpi=312)

# case 2
# rate constants kb, ke, and kt
kb = 0.5
ke = 0.05
kt = 0.25

# initial drug concentration
c0 = 500

t = np.linspace(0, 6, 1000)

# Calculate the roots of the equation
a = 1
b = kb + ke + kt
c = ke * kt
D = np.sqrt(b**2 - 4*a*c)
u1 = -1* (-b + D)/(2*a)
u2 = -1* (-b - D)/(2*a)

# Calculate the drug concentration in the blood and tissue compartments
cb = (c0 / (u2 - u1)) * ((-u1 + kt) * np.exp(-u1 * t) - (-u2 + kt) * np.exp(-u2 * t))
ct = (c0 * kb / (u1 - u2)) * (np.exp(-u1 * t) - np.exp(-u2 * t))

plt.figure(figsize=(10, 10))
plt.plot(t, cb, label='Blood')
plt.plot(t, ct, label='Tissue')
plt.xlim([0, t[-1]])
plt.xlabel('Time (hours)')
plt.ylabel('Concentration (units)')
plt.suptitle('Drug concentration pattern in Model-II - case 2')
plt.title(r'$k_b$ = {} , $k_e$ = {} , $k_t$ = {}'.format(kb, ke, kt))
plt.legend()
plt.show()
plt.savefig("Model2_case2.png", dpi=312)

# %% Model III
# case 1
# rate constants kb, kt, and ke
kb = 0.9776
ke = 0.2213
kt = 0.3293

# initial drug concentration
c0 = 500

# Define the time range
t = np.linspace(0, 6, 1000)

# Calculate the drug concentration in the arterial blood, tissue, and venous blood compartments
cab = c0 * np.exp(-kb * t)
ct = (c0 * kb / (kb - kt)) * (np.exp(-kt * t) - np.exp(-kb * t))
cvb = c0 * kb * kt * ((np.exp(-kt * t)) / ((kb - kt) * (ke - kt)) - (np.exp(-kb * t)) / ((kb - kt) * (ke - kb)) + (np.exp(-ke * t)) / ((ke - kt) * (ke - kb)))

plt.figure(figsize=(10, 10))
plt.plot(t, cab, 'r', label='Arterial Blood')
plt.plot(t, ct, 'g', label='Tissue')
plt.plot(t, cvb, 'b', label='Venous Blood')
plt.xlabel('Time (hours)')
plt.ylabel('Concentration (units)')
plt.ylim([0, c0])
plt.xlim([0, t[-1]])
plt.suptitle('Drug distribution pattern in intravenous administration with 500 units of initial drug dosage for Model-III - case 1')
plt.title(r'$k_b$ = {} , $k_e$ = {} , $k_t$ = {}'.format(kb, ke, kt))
plt.legend()
plt.show()
plt.savefig("Model3_case1.png", dpi=312)

# case 2
# rate constants kb, kt, and ke
kb = 0.5
ke = 0.05
kt = 0.25

# initial drug concentration
c0 = 500

t = np.linspace(0, 6, 1000)

# Calculate the drug concentration in the arterial blood, tissue, and venous blood compartments
cab = c0 * np.exp(-kb * t)
ct = (c0 * kb / (kb - kt)) * (np.exp(-kt * t) - np.exp(-kb * t))
cvb = c0 * kb * kt * ((np.exp(-kt * t)) / ((kb - kt) * (ke - kt)) - (np.exp(-kb * t)) / ((kb - kt) * (ke - kb)) + (np.exp(-ke * t)) / ((ke - kt) * (ke - kb)))

plt.figure(figsize=(10, 10))
plt.plot(t, cab, 'r', label='Arterial Blood')
plt.plot(t, ct, 'g', label='Tissue')
plt.plot(t, cvb, 'b', label='Venous Blood')
plt.xlabel('Time (hours)')
plt.ylabel('Concentration (units)')
plt.ylim([0, c0])
plt.xlim([0, t[-1]])
plt.suptitle('Drug distribution pattern in intravenous administration with 500 units of initial drug dosage for Model-III - case 2')
plt.title(r'$k_b$ = {} , $k_e$ = {} , $k_t$ = {}'.format(kb, ke, kt))
plt.legend()
plt.show()
plt.savefig("Model3_case2.png", dpi=312)