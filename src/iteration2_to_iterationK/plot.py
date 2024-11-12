from collections import Counter
import matplotlib.pyplot as plt

x_iterations = [1, 2, 3]
l2_clade_counts_iterations = [
    Counter({'Viridiplantae': 11, 'Opisthokonta': 6, 'NA': 6}),
    Counter({'Viridiplantae': 14, 'NA': 12, 'Opisthokonta': 7, 'Sar': 5}),
    Counter({'Viridiplantae': 14, 'NA': 13, 'Sar': 13, 'Opisthokonta': 7})
]
l3_clade_counts_iterations = [
    Counter({'Opisthokonta': 19, 'Viridiplantae': 13, 'NA': 9, 'Sar': 5, 'Amoebozoa': 1}),
    Counter({'Opisthokonta': 23, 'Viridiplantae': 15, 'Sar': 15, 'NA': 8, 'Amoebozoa': 1}),
    Counter({'Opisthokonta': 19, 'Viridiplantae': 15, 'Sar': 14, 'NA': 7})
]

# Extract unique clades from both L2 and L3
clades = sorted(set(key for counts in l2_clade_counts_iterations + l3_clade_counts_iterations for key in counts))

# Define unique markers for each clade
markers = ['o', 's', '^', 'D', 'x', 'P', '*', 'h', '+']
clade_marker_map = {clade: markers[i % len(markers)] for i, clade in enumerate(clades)}

# Plotting L2 and L3 clade frequencies over iterations
plt.figure(figsize=(15, 8))

# Set a larger marker size
marker_size = 10

# Plot for L2 Clades (Blue Color Family)
for clade in clades:
    l2_counts = [iteration_counts.get(clade, 0) for iteration_counts in l2_clade_counts_iterations]
    plt.plot(x_iterations, l2_counts, label=f'L2 {clade}', linestyle='-', marker=clade_marker_map[clade], color='blue', markersize=marker_size)

# Plot for L3 Clades (Red Color Family)
for clade in clades:
    l3_counts = [iteration_counts.get(clade, 0) for iteration_counts in l3_clade_counts_iterations]
    plt.plot(x_iterations, l3_counts, label=f'L3 {clade}', linestyle='-', marker=clade_marker_map[clade], color='red', markersize=marker_size)

# Ensuring x-axis only has integer ticks
plt.xticks(x_iterations)

plt.xlabel('Iterations', fontsize = 20)
plt.ylabel('Count', fontsize = 20)
plt.title('Clade Frequency Over Iterations for L2 and L3', fontsize = 20)
plt.legend()
# plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0., fontsize = 20)
plt.grid(True)
plt.savefig("../../result/iteration2_to_iterationK/line_plot.png")
