# The effect of migration on population models with bottlenecks
A Python package to simulate the Site Frequency Spectrum (SFS) of structured population models with bottlenecks.

We consider coalescent processes describing the limiting genealogy of a population divided in several islands and undergoing recurrent demographic bottlenecks. Such coalescent processes present three kinds of transition:
- binary mergers: happening at rate 1 for each pair of individual inside the same island;
- migration: happening at constant rate between two islands;
- simultaneous multiple  mergers and migrations:  happening at the occurrence of a bottlenecks, their effect depends on the type of bottleneck.

We consider two types of bottlenecks:

### Drastic Bottlenecks 
In this case the size of the affected island is reduced to a finite number of individual for a random, but finite, amount of generations.

### Soft Bottlenecks
In this case we run, on the affected island, a Kingman coalescent for a given amount of time and then we instantly update the genealogy based on the final state of the Kingman coalescent.

One can find in the code how to produce the figures in the paper, but can also produce further examples (see Usage).

## Usage

This code can be used to simulate the SFS of structured population models undergoing migration and recurring bottlenecks, based on the following parameters:

``sample`` array of the number of individuals sampled from each island (size is number of islands)

``migr`` rate at which migration between islands occur

``bott`` rate at which bottlenecks occur

``bott_size`` size of the bottleneck; the size of the bottlenecks is the length in the soft type while it is the actual number of surviving individuals in the drastic case.

``bott_type`` type of bottleneck: 'd' for drastic; 's' for soft

``n`` number of iterations

The function ``soft(sample,migr,bott,bott_size)``, resp. ``drastic(sample,migr,bott,bott_size)`` generates one realization of the SFS of the model with soft bottlenecks with the chosen parameters. The function ``SFS(bott_type,sample,migr,bott,bott_size,n)`` generates the average SFS for the model with ``bott_type`` bottlenecks, over ``n`` iterations. The code can also be used to produce the plot of the logit transform of the average SFS.

### Example
```
# EXAMPLE USAGE
bott_type = 'd' # choses the type of bottleneck: 'd' for drastic; 's' for soft
sample = [25,25,25,25]  # array of initial sample; in this case 4 islands
migr = 1    # migration rate
bott = 0.1  # bottleneck rate
bott_size = 4    # bottleneck size
n = 3000    # number of iterations

print SFS(bott_type,sample,migr,bott,bott_size,n)
# print the average SFS for the selected parameters

# plot the logit transform of the average SFS
a = [bott_type,sample,migr,bott,bott_size,n]
fig = plt.figure(figsize=(10, 10))
plt.plot(
    logit([i / sum(a[1]) for i in range(1, sum(a[1]))]),
    logit(SFS(a[0], a[1], a[2], a[3], a[4], a[5])),
    "black",
    linestyle="-",
    linewidth=2,
    label="Example plot",
)

plt.legend(loc="upper right", fontsize=16)
plt.xlabel("logit mutant frequency", fontsize=17)
plt.ylabel("logit SFS", fontsize=17)
plt.show()
```

## License
This project is licensed under the GNU Affero General Public License v3.0.
