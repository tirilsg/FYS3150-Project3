import numpy as np
import matplotlib.pyplot as plt

#we load and store the data for the simulations of 100 particles in the trap with names that tells us what they contain
nointertest = np.loadtxt("nointertest.txt") 
wv_nointerest = nointertest[:, 0]
q_nointerest = nointertest[:, 1:]

data_intertest = np.loadtxt("intertest.txt")
wv_interest = data_intertest[:, 0]
q_interest = data_intertest[:, 1:]

nointertestzoom = np.loadtxt("nointertestzoom.txt")
wv_nointerestzoom = nointertestzoom[:, 0]
q_nointerestzoom = nointertestzoom[:, 1:]

intertestzoom = np.loadtxt("intertestzoom.txt")
wv_interestzoom = intertestzoom[:, 0]
q_interestzoom = intertestzoom[:, 1:]

amplitudes = [0.1, 0.4, 0.7] #amplitudes used in our simulations, with and without consideration for particle interaction
#plots the fraction of total particles within the trap as function of w_v
plt.figure(figsize=(12, 6))
plt.suptitle("Fraction of trapped particles as a funciton of resonance")
plt.subplot(1, 2, 1)
for idx, amp in enumerate(amplitudes):
    plt.plot(wv_nointerest, q_nointerest[:, idx] / 100, label=f"Amplitude f {amp}")
plt.ylabel("% trapped particles")
plt.xlabel(r"$w_v$ (MHz)")
plt.title("Without interactions between particles")
plt.legend()
plt.grid()
plt.subplot(1, 2, 2)
for idx, amp in enumerate(amplitudes):
    plt.plot(wv_interest, q_interest[:, idx] / 100, label=f"Amplitude f {amp}") 
plt.xlabel(r"$w_v$ (MHz)")
plt.ylabel("% trapped particles")
plt.title("Interactions between particles")
plt.legend()
plt.grid()
plt.savefig(r"Trappedfrac.pdf")
plt.show()

#uses the zoomed-in data to create similar plots
plt.figure(figsize=(12, 6))
plt.suptitle(r"Plot zoomed into resonance of interest, interval $w_v \in$ [1.0,1.8]")
plt.subplot(1, 2, 1)
for idx, amp in enumerate(amplitudes):
    plt.plot(wv_nointerestzoom, q_nointerestzoom[:, idx] / 100, label=f"Amplitude f {amp}")
plt.xlabel(r"$w_v$ (MHz)")
plt.ylabel("% trapped particles")
plt.title("No interactions between particles")
plt.legend()
plt.grid()
plt.subplot(1, 2, 2)
for idx, amp in enumerate(amplitudes):
    plt.plot(wv_interestzoom, q_interestzoom[:, idx] / 100, label=f"Amplitude f {amp}")
plt.xlabel(r"$w_v$ (MHz)")
plt.ylabel("% trapped particles")
plt.title("Interactions between particles")
plt.legend()
plt.grid()
plt.savefig(r"Trappedfraczoom.pdf")
plt.show()