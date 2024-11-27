from icosalattice.PlotDataOnMap import plot_variable_interpolated_from_dict


fp = "/home/kuhron/elegen/scripts/test_data.txt"

with open(fp) as f:
    lines = f.readlines()
lines = [x.strip().split(",") for x in lines]
pcs, vals = zip(*lines)
vals = [float(x) for x in vals]
d = dict(zip(pcs, vals))

plot_variable_interpolated_from_dict(d, dots_per_degree=10)
