import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

file_path = "speedup.data"
data = pd.read_csv(file_path, delim_whitespace=True)

ideal_speedup = np.array(data["n_proc"])

plt.figure(figsize=(8, 5))

plt.plot(data['n_proc'], data['speedup'], marker='o', linestyle='-', color='b', label='Speedup Actuel')
plt.plot(data['n_proc'], ideal_speedup, linestyle='--', color='r', label='Speedup Idéal (S(p) = p)')

plt.xlabel("Nombre de processus")
plt.ylabel("Speedup")
plt.title("Speedup vs. Nombre de processus")
plt.legend()
plt.grid(True)

plt.savefig("speedup.png", dpi=300, bbox_inches='tight')
