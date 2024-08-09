import random
import math
import numpy
import scipy.special
from matplotlib import pyplot as plt
from scipy.special import logit, expit


## this function divides "set_size" individuals randomly in "blocks_number" sets each one going randomly in other "islands_number" boxes
def random_partition(set_size, blocks_number, islands_number):
    box_size = [
        k
        for k in numpy.random.multinomial(
            set_size, [1 / blocks_number] * blocks_number, size=None
        )
    ]
    box_migr = numpy.random.multinomial(
        blocks_number, [1 / islands_number] * islands_number, size=None
    )
    boxes = [[] for k in range(islands_number)]

    for k in range(islands_number):
        box_index = random.sample(range(len(box_size)), box_migr[k])
        boxes[k] = [box_size[r] for r in box_index]
        for r in range(len(box_index)):
            box_size.pop(max(box_index))
            box_index.remove(max(box_index))

    return boxes


def merging_soft(
    lineages,
    size,
    parent,
    next_parent,
    flags,
    node_time,
    i,
    t,
):
    if len(lineages) >= 2:
        next_children = random.sample(range(len(lineages)), size)
        for r in range(size):
            parent[lineages[next_children[r]]] = next_parent
        # children[next_parent] = [lineages[next_children[r]] for r in range(size)]
        for r in range(size):
            flags[next_parent] = flags[next_parent] + flags[lineages[next_children[r]]]
        for r in range(size):
            lineages.pop(max(next_children))
            next_children.remove(max(next_children))
        lineages.append(next_parent)
        node_time[next_parent] = t
        next_parent = next_parent + 1
        return (lineages, parent, next_parent, flags, node_time)


def merging_drastic(
    lineages,
    size,
    parent,
    next_parent,
    flags,
    gen_parent,
    node_time,
    i,
    t,
):
    if len(lineages) >= 2:
        next_children = random.sample(range(len(lineages)), size)
        for r in range(size):
            parent[lineages[next_children[r]]] = next_parent
        for r in range(size):
            flags[next_parent] = flags[next_parent] + flags[lineages[next_children[r]]]
        for r in range(size):
            lineages.pop(max(next_children))
            next_children.remove(max(next_children))
        gen_parent.append(next_parent)
        node_time[next_parent] = t
        next_parent = next_parent + 1
        return (gen_parent, parent, next_parent, flags, node_time)


def soft(sample, migr, bott, bott_size):
    ## sample=[size isalnd 1,size island 2,...,size island n]    len(sample)=number of islands
    active_lineages = [[] for i in range(len(sample))]
    active_lineages[0] = [x for x in range(sample[0])]
    size = 0
    for i in range(0, len(sample)):
        active_lineages[i] = [x for x in range(size, size + sample[i])]
        size = size + sample[i]
    parent = [-1 for i in range(0, 2 * size - 1)]
    node_time = [0 for i in range(0, 2 * size - 1)]
    flags = [1 for i in range(0, size)] + [0 for i in range(size - 1)]
    next_parent = size
    t = 0

    while sum(len(active_lineages[i]) for i in range(len(active_lineages))) > 1:
        product = 1
        total_coal = sum(
            scipy.special.binom(len(active_lineages[i]), 2)
            for i in range(len(active_lineages))
        )
        total_migr = sum(len(active_lineages[i]) for i in range(len(active_lineages)))

        if max(len(active_lineages[i]) for i in range(len(active_lineages))) > 1:
            coal = random.expovariate(total_coal)
        else:
            coal = math.inf
        if migr > 0 and len(sample) > 1:
            migr = random.expovariate(migr * total_migr)
        else:
            migr = math.inf
        if bott > 0:
            bott = random.expovariate(bott)
        else:
            bott = math.inf
        t = t + min(coal, migr, bott)

        if min(coal, migr, bott) == coal:
            total_coal = sum(
                scipy.special.binom(len(active_lineages[r]), 2)
                for r in range(len(active_lineages))
            )
            i = numpy.random.choice(
                len(active_lineages),
                1,
                p=[
                    scipy.special.binom(len(active_lineages[r]), 2) / total_coal
                    for r in range(len(active_lineages))
                ],
            )[
                0
            ]  # random choice of the island
            [active_lineages[i], parent, next_parent, flags, node_time] = merging_soft(
                active_lineages[i], 2, parent, next_parent, flags, node_time, i, t
            )

        elif min(coal, migr, bott) == migr:
            total_migr = sum(
                len(active_lineages[r]) for r in range(len(active_lineages))
            )
            i = 1
            j = 1

            while i == j:
                i = numpy.random.choice(
                    len(active_lineages),
                    1,
                    p=[
                        len(active_lineages[r]) / total_migr
                        for r in range(len(active_lineages))
                    ],
                )[0]
                j = random.choice(range(len(active_lineages)))

            next_migrant = random.choice(range(len(active_lineages[i])))
            active_lineages[j].append(active_lineages[i][next_migrant])
            active_lineages[i].pop(next_migrant)

        elif min(coal, migr, bott) == bott:
            i = random.choice(range(len(active_lineages)))  # random choice of an island
            s = 0
            if len(active_lineages[i]) >= 2:
                while s < bott_size and len(active_lineages[i]) >= 1:
                    if len(active_lineages[i]) > 1:
                        coal_bott = random.expovariate(
                            scipy.special.binom(len(active_lineages[i]), 2)
                        )
                    else:
                        coal_bott = math.inf
                    if migr > 0 and len(sample) > 1:
                        migr_bott = random.expovariate(
                            migr * len(active_lineages[i])
                        )
                    else:
                        migr_bott = math.inf

                    s = s + min(coal_bott, migr_bott)
                    if (
                        min(coal_bott, migr_bott) == coal_bott
                        and len(active_lineages[i]) > 1
                    ):
                        [
                            active_lineages[i],
                            parent,
                            next_parent,
                            flags,
                            node_time,
                        ] = merging_soft(
                            active_lineages[i],
                            2,
                            parent,
                            next_parent,
                            flags,
                            node_time,
                            i,
                            t,
                        )

                    elif (
                        min(coal_bott, migr_bott) == migr_bott and len(sample) > 1
                    ):
                        j = random.choice(range(len(active_lineages)))
                        while i == j:
                            j = random.choice(range(len(active_lineages)))

                        next_migrant = random.choice(range(len(active_lineages[i])))
                        active_lineages[j].append(active_lineages[i][next_migrant])
                        active_lineages[i].pop(next_migrant)

    frequency = [0 for i in range(size - 1)]
    for k in range(1, size):
        for i in range(len(flags)):
            if flags[i] == k:
                frequency[k - 1] = (
                    frequency[k - 1] + node_time[parent[i]] - node_time[i]
                )
    tot_frequency = 0
    for k in range(len(frequency)):
        tot_frequency = tot_frequency + frequency[k]
    for k in range(len(frequency)):
        frequency[k] = frequency[k] / tot_frequency
    return frequency


def drastic(sample, migr, bott, bott_size):
    # sample=[size island 1,size island 2,...,size island n]    len(sample)=number of islands
    # size = sum(sample[i] for i in range(len(sample)))
    active_lineages = [[] for i in range(len(sample))]
    active_lineages[0] = [x for x in range(sample[0])]
    size = 0
    for i in range(0, len(sample)):
        active_lineages[i] = [x for x in range(size, size + sample[i])]
        size = size + sample[i]
    parent = [-1 for i in range(0, 2 * size - 1)]
    node_time = [0 for i in range(0, 2 * size - 1)]
    flags = [1 for i in range(0, size)] + [0 for i in range(size - 1)]
    next_parent = size
    t = 0

    while sum(len(active_lineages[i]) for i in range(len(active_lineages))) > 1:
        product = 1
        total_coal = sum(
            scipy.special.binom(len(active_lineages[i]), 2)
            for i in range(len(active_lineages))
        )
        total_migr = sum(len(active_lineages[i]) for i in range(len(active_lineages)))

        if max(len(active_lineages[i]) for i in range(len(active_lineages))) > 1:
            coal = random.expovariate(total_coal)
        else:
            coal = math.inf
        if migr > 0 and len(sample) > 1:
            migr = random.expovariate(migr * total_migr)
        else:
            migr = math.inf
        if bott > 0:
            bott = random.expovariate(bott)
        else:
            bott = math.inf
        t = t + min(coal, migr, bott)

        if min(coal, migr, bott) == coal:  # COALESCENCE
            total_coal = sum(
                scipy.special.binom(len(active_lineages[r]), 2)
                for r in range(len(active_lineages))
            )
            i = numpy.random.choice(
                len(active_lineages),
                1,
                p=[
                    scipy.special.binom(len(active_lineages[r]), 2) / total_coal
                    for r in range(len(active_lineages))
                ],
            )[
                0
            ]  # random choice of the island
            [active_lineages[i], parent, next_parent, flags, node_time] = merging_soft(
                active_lineages[i], 2, parent, next_parent, flags, node_time, i, t
            )

        elif min(coal, migr, bott) == migr:
            total_migr = sum(
                len(active_lineages[r]) for r in range(len(active_lineages))
            )
            i = 1
            j = 1

            while i == j:
                i = numpy.random.choice(
                    len(active_lineages),
                    1,
                    p=[
                        len(active_lineages[r]) / total_migr
                        for r in range(len(active_lineages))
                    ],
                )[0]
                j = random.choice(range(len(active_lineages)))

            next_migrant = random.choice(range(len(active_lineages[i])))
            active_lineages[j].append(active_lineages[i][next_migrant])
            active_lineages[i].pop(next_migrant)

        elif min(coal, migr, bott) == bott:
            i = random.choice(range(len(active_lineages)))  # random choice of an island
            lenght_bott = (
                numpy.random.poisson(bott_size) + 1
            )  
            for g in range(lenght_bott):
                gen_parent = [
                    [] for k in range(len(active_lineages))
                ]  # matrix where to put the new parents created during the bottleneck
                boxes = random_partition(
                    len(active_lineages[i]), bott_size, len(active_lineages)
                )  # create a box (boxes[j]) for every island wich is made of boxes (boxes[j][k]) which represent the mergings
                for j in range(len(active_lineages)):
                    for k in range(len(boxes[j])):
                        if (
                            boxes[j][k] > 1
                        ):  # the merging is between individual in island i but the parent goes to island j
                            [
                                gen_parent[j],
                                parent,
                                next_parent,
                                flags,
                                node_time,
                            ] = merging_drastic(
                                active_lineages[i],
                                boxes[j][k],
                                parent,
                                next_parent,
                                flags,
                                gen_parent[j],
                                node_time,
                                j,
                                t,
                            )

                        elif boxes[j][k] == 1:
                            next_migrant = random.choice(range(len(active_lineages[i])))
                            gen_parent[j].append(active_lineages[i][next_migrant])
                            active_lineages[i].pop(next_migrant)

                for j in range(len(active_lineages)):
                    active_lineages[j] = active_lineages[j] + gen_parent[j]

    frequency = [0 for i in range(size - 1)]
    for k in range(1, size):
        for i in range(len(flags)):
            if flags[i] == k:
                frequency[k - 1] = (
                    frequency[k - 1] + node_time[parent[i]] - node_time[i]
                )
    tot_frequency = 0
    for k in range(len(frequency)):
        tot_frequency = tot_frequency + frequency[k]
    for k in range(len(frequency)):
        frequency[k] = frequency[k] / tot_frequency
    return frequency


def SFS(bott_type, sample, migr, bott, bott_size, n):
    p = [0 for i in range(sum(sample) - 1)]
    if bott_type == "s":
        for i in range(n):
            print(i)
            s = soft(
                sample, migr, bott, bott_size
            )  # s is the frequency obtained at the current step i for i=1,...,n
            for i in range(sum(sample) - 1):
                p[i] = p[i] + s[i]
        for i in range(sum(sample) - 1):
            p[i] = p[i] / n
    elif bott_type == "d":
        for i in range(n):
            print(i)
            s = drastic(
                sample, migr, bott, bott_size
            )  # s is the frequency obtained at the current step i for i=1,...,n
            for i in range(sum(sample) - 1):
                p[i] = p[i] + s[i]
        for i in range(sum(sample) - 1):
            p[i] = p[i] / n
    else:
        p = 0
        print("error in bottlenck type")
    return p


# EXAMPLE USAGE
bott_type = 'd' # choses the type of bottleneck: 'd' for drastic; 's' for soft
sample = [25,25,25,25]  # array of initial sample; in this case 4 islands
migr = 1    # migration rate
bott = 0.1  # bottleneck rate
bott_size = 4    # bottleneck size
##the size of the bottlenecks is the length in the soft type while it is the actual number of surviving individuals in the drastic case.
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

# FIGURE6 left-hand side
a=['d',[100,0,0,0],0.1,10,4,30000]
b=['d',[100,0,0,0],10,10,4,30000]
c=['s',[100,0,0,0],0.1,100,1,30000]
d=['s',[100,0,0,0],10,100,1,30000]


# FIGURE6 right-hand side
##a=['d',[25,25,25,25],0.1,10,1,30000]
##b=['d',[25,25,25,25],10,10,1,30000]
##c=['s',[25,25,25,25],0.1,100,1,30000]
##d=['s',[25,25,25,25],10,100,1,30000]


# FIGURE7
##a=['d',[100],0,100,5,30000]
##b=['d',[100],0,100,20,30000]



fig = plt.figure(figsize=(10, 10))
plt.plot(
    logit([i / sum(a[1]) for i in range(1, sum(a[1]))]),
    logit(SFS(a[0], a[1], a[2], a[3], a[4], a[5])),
    "black",
    linestyle="-",
    linewidth=2,
    label="Drastic (low migration)",
)
plt.plot(
    logit([i / sum(b[1]) for i in range(1, sum(b[1]))]),
    logit(SFS(b[0], b[1], b[2], b[3], b[4], b[5])),
    "grey",
    linestyle="-",
    linewidth=2,
    label="Drastic (high migration)",
)
plt.plot(
    logit([i / sum(c[1]) for i in range(1, sum(c[1]))]),
    logit(SFS(c[0], c[1], c[2], c[3], c[4], c[5])),
    "black",
    linestyle="dashdot",
    linewidth=2,
    label="Soft (low migration)",
)
plt.plot(
    logit([i / sum(d[1]) for i in range(1, sum(d[1]))]),
    logit(SFS(d[0], d[1], d[2], d[3], d[4], d[5])),
    "grey",
    linestyle="dashdot",
    linewidth=2,
    label="Soft (high migration)",
)
plt.legend(loc="upper right", fontsize=16)
plt.xlabel("logit mutant frequency", fontsize=17)
plt.ylabel("logit SFS", fontsize=17)


plt.show()
