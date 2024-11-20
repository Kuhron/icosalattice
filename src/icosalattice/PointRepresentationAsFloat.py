import icosalattice.StartingPoints as sp


DIRECTION_CODE_TO_BASE_FOUR = {
    "0": 0b00, "1": 0b10,
    "3": 0b01, "2": 0b11,
}
DIRECTION_BASE_FOUR_TO_CODE = {x:y for y,x in DIRECTION_CODE_TO_BASE_FOUR.items()}



def point_code_to_float(pc):
    # a representation of point code as a float object rather than string
    # hopefully without data corruption by float rounding because denominators are always powers of 2
    res = sp.STARTING_POINT_CODE_TO_FLOAT[pc[0]]
    for i, x in enumerate(pc[1:]):
        denom = 4**(i+1)
        res += DIRECTION_CODE_TO_BASE_FOUR[x]/denom
    return res


def point_float_to_code(x):
    x_orig = x
    n, x = divmod(x, 1)
    s = sp.STARTING_POINT_FLOAT_TO_CODE[n]
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
    import random
    
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
