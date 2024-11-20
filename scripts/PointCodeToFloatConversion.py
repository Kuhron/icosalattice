import random

import icosalattice.StartingPoints as sp


# STARTING_POINT_CODE_TO_BINARY = {
#     "A": "0", "B": "1", "C": "10", "D": "11",
#     "E": "100", "F": "101", "G": "110", "H": "111",
#     "I": "1000", "J": "1001", "K": "1010", "L": "1011",
# }
# STARTING_POINT_BINARY_TO_CODE = {v:k for k,v in STARTING_POINT_CODE_TO_BINARY.items()}
# DIRECTION_CODE_TO_BINARY = {"0": "00", "1": "10", "2": "11", "3": "01"}
# DIRECTION_BINARY_TO_CODE = {v:k for k,v in DIRECTION_CODE_TO_BINARY.items()}


# def point_code_to_binary(pc):
#     return STARTING_POINT_CODE_TO_BINARY[pc[0]] + "." + "".join(DIRECTION_CODE_TO_BINARY[x] for x in pc[1:])


# def binary_to_point_code(pc_bin):
#     spc_bin, xs = pc_bin.split(".")
#     assert len(xs) % 2 == 0
#     return STARTING_POINT_BINARY_TO_CODE[spc_bin] + "".join(DIRECTION_BINARY_TO_CODE[xs[2*i: 2*i+2]] for i in range(len(xs) // 2))


STARTING_POINT_CODE_TO_FLOAT = {pc: float(i) for i, pc in enumerate(sp.STARTING_POINT_CODES)}
STARTING_POINT_FLOAT_TO_CODE = {x: pc for pc, x in STARTING_POINT_CODE_TO_FLOAT.items()}
DIRECTION_CODE_TO_BASE_FOUR = {
    "0": 0b00, "1": 0b10,
    "3": 0b01, "2": 0b11,
}
DIRECTION_BASE_FOUR_TO_CODE = {x:y for y,x in DIRECTION_CODE_TO_BASE_FOUR.items()}

def point_code_to_float(pc):
    # a representation of point code as a float object rather than string
    # hopefully without data corruption by float rounding because denominators are always powers of 2
    res = STARTING_POINT_CODE_TO_FLOAT[pc[0]]
    for i, x in enumerate(pc[1:]):
        denom = 4**(i+1)
        res += DIRECTION_CODE_TO_BASE_FOUR[x]/denom
    return res


def point_float_to_code(x):
    x_orig = x
    n, x = divmod(x, 1)
    s = STARTING_POINT_FLOAT_TO_CODE[n]
    iterations = 0
    max_iterations = 128
    while x > 0:
        x *= 4
        y,x = divmod(x, 1)
        s += DIRECTION_BASE_FOUR_TO_CODE[y]
        iterations += 1
        if iterations > max_iterations:
            raise ValueError(f"invalid float representation of point code; remainder found: {x_orig}")
    return s


if __name__ == "__main__":
    i = 0
    while True:
        pc = random.choice(sp.STARTING_POINT_CODES[2:]) + "".join(random.choice("0123") for i in range(random.randrange(9)))
        while pc[-1] == "0":
            pc = pc[:-1]
        
        fpc = point_code_to_float(pc)
        print(pc)
        print(fpc)
        pc2 = point_float_to_code(fpc)
        assert pc == pc2, f"{pc} != {pc2}"

        i += 1
        # if i % 100000 == 0:
        #     print(i)

        print()
