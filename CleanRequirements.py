import re

req_fp = "requirements.txt"
ignore_fp = "dev_requirements.txt"

with open(ignore_fp) as f:
    ignore_patterns = [x.strip() for x in f.readlines()]
ignore_patterns = [x for x in ignore_patterns if x != ""]

with open(req_fp) as f:
    lines = f.readlines()

lines = [l for l in lines if not any(re.match(pattern, l) for pattern in ignore_patterns)]

with open(req_fp, "w") as f:
    for l in sorted(lines):
        f.write(l)


