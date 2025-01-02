import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection



R3 = 3**0.5
Y_BOTTOM = -1/3 * R3/2
Y_TOP = 2/3 * R3/2
Y_SUB_BOTTOM = Y_TOP + 2*(Y_BOTTOM - Y_TOP)
X_LEFT = -0.5
X_MID = 0
X_RIGHT = 0.5
pA = (X_MID, Y_TOP)
pC = (X_RIGHT, Y_BOTTOM)
pK = (X_LEFT, Y_BOTTOM)
pL = (X_MID, Y_SUB_BOTTOM)


def ld_to_xy_centered_upward_triangle(l, d):
    x0, y0 = pC
    xl, yl = pA
    xd, yd = pL
    dxl = xl - x0
    dyl = yl - y0
    dxd = xd - x0
    dyd = yd - y0
    x = x0 + l * dxl + d * dxd
    y = y0 + l * dyl + d * dyd
    return x, y


def toggle_visible_on_plot(event, graph_objects_by_label, fig):
    text_item = event.artist
    assert type(text_item) is matplotlib.text.Text
    label_with_checkbox = text_item.get_text()
    is_visible = label_with_checkbox.startswith("[X] ")
    if not is_visible:
        assert label_with_checkbox.startswith("[O] "), f"label without checkbox: {label_with_checkbox}"
    
    label = label_with_checkbox[len("[X] "):]
    graph_objects_by_label[label].set_visible(not is_visible)
    # legend_item.set_visible(not is_visible)
    if is_visible:
        # making it invisible
        text_item.set_text("[O] " + label)
    else:
        # making it visible
        text_item.set_text("[X] " + label)
    fig.canvas.draw()


def xy_distance(xy1, xy2):
    x1, y1 = xy1
    x2, y2 = xy2
    return ((x1-x2)**2 + (y1-y2)**2)**0.5


def transform_xy_1(x, y):
    # doesn't preserve edges
    r_raw = (x**2 + y**2)**0.5
    max_r = Y_TOP
    r = r_raw / max_r
    a = 0.8
    c = 0.5
    r1 = 1/(1+np.exp(a * np.log(1/r - 1) + np.log(1/c - 1)))
    new_r_raw = r1 * max_r
    x *= new_r_raw / r_raw
    y *= new_r_raw / r_raw
    return x, y


def transform_xy_2(x, y):
    a = xy_distance((x,y), pA)
    c = xy_distance((x,y), pC)
    k = xy_distance((x,y), pK)
    print(a, c, k)
    


# get the coordinates for lattice points in (l,d) space
iterations = 5
# transform_xy = lambda x, y: (x, y)
transform_xy = transform_xy_1

step = 2 ** -iterations
nums = np.arange(0, 1 + step, step)
ls = nums
ds = nums

ijs_of_points = [(i,j) for i in range(len(nums)) for j in range(len(nums)) if j >= i]

ijs_on_l_lines = [[(i,j) for j in range(i, len(nums))] for i in range(1, len(nums)-1)]
ijs_on_dl_lines = [[(j-i,j) for j in range(i, len(nums))] for i in range(1, len(nums)-1)]
ijs_on_d_lines = [[(i,j) for i in range(j+1)] for j in range(1, len(nums)-1)]

ij_to_xy = {}
xs = []
ys = []
for i,j in ijs_of_points:
    x, y = ld_to_xy_centered_upward_triangle(ls[j], ds[i])
    x, y = transform_xy(x, y)
    xs.append(x)
    ys.append(y)
    ij_to_xy[(i,j)] = (x,y)

# plot stuff

fig, ax = plt.subplots()

graph_objects_by_label = {}
for ijs_on_lines_collection, color, label in zip([ijs_on_l_lines, ijs_on_dl_lines, ijs_on_d_lines], ["indianred", "goldenrod", "mediumseagreen"], ["L", "DL", "D"]):
    lines = []
    for ijs in ijs_on_lines_collection:
        xys_this_line = []
        for i,j in ijs:
            x,y = ij_to_xy[(i,j)]
            xys_this_line.append((x,y))
        lines.append(xys_this_line)
    collection = LineCollection(lines, colors=[color for _ in lines], alpha=0.7, label=label, zorder=0)
    # plt.plot(xs_this_line, ys_this_line, c=color, label=label, zorder=0)
    ax.add_collection(collection)
    graph_objects_by_label[label] = collection

plt.plot([X_RIGHT, X_MID, X_LEFT, X_RIGHT], [Y_BOTTOM, Y_TOP, Y_BOTTOM, Y_BOTTOM], c="k", zorder=1)
lattice_points = plt.scatter(xs, ys, c="gray", label="lattice", alpha=0.5, zorder=2)
graph_objects_by_label["lattice"] = lattice_points

plt.gca().set_aspect("equal")

# make legend where things can be toggled using their labels: https://www.youtube.com/watch?v=QcRE8dinls4 (Jie Jenn)
legend = plt.legend(loc="upper right")
for text_item in legend.get_texts():
    text_item.set_picker(True)
    text_item.set_text("[X] " + text_item.get_text())
plt.connect("pick_event", lambda event: toggle_visible_on_plot(event, graph_objects_by_label=graph_objects_by_label, fig=fig))

plt.show()
