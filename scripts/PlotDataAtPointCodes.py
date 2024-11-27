import icosalattice.PlottingUtil as pu


fp = "/home/kuhron/elegen/scripts/test_data.txt"

with open(fp) as f:
    lines = f.readlines()
lines = [x.strip().split(",") for x in lines]
pcs, vals = zip(*lines)
vals = [float(x) for x in vals]

print(pcs, vals)