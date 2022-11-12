import matplotlib.pyplot as plt


def energy_plot(plot_y, prefix):
    plt.ylabel('Energy', fontsize=14)
    plt.xlabel('Iteration', fontsize=14)
    plt.title('Energy vs Iteration', fontsize=15)
    plt.plot([i for i in range(len(plot_y))], plot_y)
    plt.tight_layout()
    plt.savefig("." + prefix + '.jpg')
    plt.clf()
