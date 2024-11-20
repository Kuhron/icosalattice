# another number system for lattice points

# for converting between {point code, number array, lookup number}
LETTER_TO_NUMBER_DICT = {c:i for i,c in enumerate("CDEFGHIJKL")}
LETTER_TO_NUMBER_DICT["A"] = -2
LETTER_TO_NUMBER_DICT["B"] = -3
NUMBER_TO_LETTER_DICT = {i:c for c,i in LETTER_TO_NUMBER_DICT.items()}
INITIAL_POINT_LOOKUP_NUMBERS = [-2, -3, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]


