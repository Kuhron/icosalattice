import matplotlib.pyplot as plt

from icosalattice.Adjacency import get_adjacency_from_point_code


def plot_adjacency_of_point_code(pc):
    rt3 = 3**0.5
    adj = get_adjacency_from_point_code(pc)
    xs = [0]
    ys = [0]
    labels = [pc]
    direction_to_coords = {
        1: (-2, 0),
        2: (-1, -rt3),
        3: (1, -rt3),
        -1: (2, 0),
        -2: (1, rt3),
        -3: (-1, rt3),
    }
    ax = plt.gca()
    for direction, pc2 in adj.items():
        x,y = direction_to_coords[direction]
        if pc2 is not None:
            xs.append(x)
            ys.append(y)
            labels.append(pc2)
            plt.arrow(x*0.1, y*0.1, x*0.75, y*0.75, head_width=0.05, head_length=0.05, color="k")
            ax.annotate(f"{direction:+}", xy=(x/2+0.1, y/2+0.1/2))
    plt.scatter(xs, ys)
    ax.set_aspect("equal")
    ax.set_xlim(-2.5, 2.5)
    ax.set_ylim(-2.5, 2.5)
    for x,y,label in zip(xs, ys, labels):
        ax.annotate(label, xy=(x,y), xytext=(x+0.1,y+0.1/2))
    plt.show()
