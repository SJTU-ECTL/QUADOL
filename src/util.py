import math
import heapq
import itertools


def bit_majority(tt1: list, tt2: list, tt3: list) -> list:
    # Check that all truth tables have the same length
    if len(tt1) != len(tt2) or len(tt1) != len(tt3):
        raise ValueError('All truth tables must have the same length')
    # Initialize the majority truth table
    majority_tt = [0] * len(tt1)
    # Iterate over the truth tables
    for i in range(len(tt1)):
        # If the number of 1s is greater than or equal to 2, the majority bit is 1
        if tt1[i] + tt2[i] + tt3[i] >= 2:
            majority_tt[i] = 1
    return majority_tt


def list_combinations(lst: list, k: int = 2) -> list:
    """
    Find all possible combinations of `k` elements from `lst`
    Sample Input:
    lst = [1, 2, 3, 4]
    k = 2
    Output: [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
    """
    return list(itertools.combinations(lst, k))


def has_common_element(g_structure: set) -> bool:
    all_combine = list_combinations([0, 1, 2, 3])
    if len(g_structure) != 2:
        raise ValueError('Graph structure must have exactly two in Graph Type-22')
    return len(set(all_combine[list(g_structure)[0]]) & set(all_combine[list(g_structure)[1]])) > 0


def type2_combinations(A: list, B: list, C: list, D: list) -> list:
    return [[(A, B), (C, D)], [(A, D), (B, C)], [(A, C), (B, D)]]


def pair_elements(lst: list, k: int) -> dict:
    """
    Pair all elements, except the `k`-th element.
    Sample Input:
    lst = ['a', 'b', 'c', 'd', 'e', 'f']
    k = 2
    Output: [('a', 4), ('b', 3), ('c', 5), ('d', 2), ('e', 1), ('f', 0)]
    """
    n = len(lst)
    pairs = {}
    count = n - 1
    for i in range(n):
        if i == k:
            pairs[n - 1] = lst[i]
        else:
            pairs[count - 1] = lst[i]
            count -= 1
    return pairs


def least_elements(nums: list, size: int = 1) -> list[tuple[int, int]]:
    """
    Find the least two elements in a list
    Sample Input:
    nums = [5, 2, 9, 1, 7, 4]
    size = 2
    Output: [(3, 1), (1, 2)]
    Output Format: [(idx, val), ...]
    """
    return heapq.nsmallest(size, enumerate(nums), key=lambda x: x[1])


def reorder_tt(f_tt: list, f_order: list, g_order: list) -> list:
    """
    Reorder the truth table of f_tt to g_tt
    Sample Input:
    f_tt = ['w', 'x', 'y', 'z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l']
    f_order = ['a', 'b', 'c', 'd']
    g_order = ['b', 'd', 'c', 'a']
    Output: ['w', 'y', 'x', 'z', 'a', 'c', 'b', 'd', 'e', 'g', 'f', 'h', 'i', 'k', 'j', 'l']
    """
    # Determine the mapping from f_order to g_order
    mapping = [f_order.index(b) for b in g_order]

    g_tt = [0] * len(f_tt)
    for i in range(len(f_tt)):
        # Convert index to binary and pad with zeros to the left to get n bits
        binary_index = list(format(i, '0' + str(len(f_order)) + 'b'))
        # Reorder the binary index according to the mapping
        reordered_index = [binary_index[m] for m in mapping]
        # Convert the reordered binary index back to decimal
        new_index = int(''.join(reordered_index), 2)
        # Assign the value from the old truth table to the new position in the reordered truth table
        g_tt[new_index] = f_tt[i]
    return g_tt


def list_complement(tt: list) -> list:
    return [1 - i for i in tt]


def HD_2tt(tt1: list, tt2: list) -> float:
    """
    Function to calculate the Hamming distance between two truth tables.
    """
    if len(tt1) != len(tt2):
        raise ValueError("Truth tables must have the same length")
    return sum(bit1 != bit2 for bit1, bit2 in zip(tt1, tt2))


def HD_3tt(tt1: list, tt2: list, tt3: list) -> float:
    """
    Function to calculate the Hamming distance between two truth tables.
    """
    if len(tt1) != len(tt2) or len(tt1) != len(tt3):
        raise ValueError('All truth tables must have the same length')
    bit_major = bit_majority(tt1, tt2, tt3)
    return sum([HD_2tt(tt1, bit_major),
                HD_2tt(tt2, bit_major),
                HD_2tt(tt3, bit_major)])


def decompose_tt(tt: list, input_idx: int, parameter: bool = True) -> list:
    """
    Function to decompose the truth table.
    If the parameter is true, decompose `idx = 1`.
    If the parameter is false, decompose `idx = 0`.
    Sample Input:
    tt = [1, 0, 1, 0, 1, 0, 1, 1]
    input_idx = 1, parameter = True
    Output: [1, 0, 1, 1]
    """
    new_tt = []
    for i in range(len(tt)):
        # Convert i to binary and pad with zeros to get num_inputs bits
        bin_str = format(i, '0' + str(int(math.log2(len(tt)))) + 'b')
        # Check if the bit at input_index is 1 or 0
        if bin_str[input_idx] == ('1' if parameter else '0'):
            # If it is, add the corresponding value from the truth table to the new truth table
            new_tt.append(tt[i])
    return new_tt


def rebuild_tt(tt_0: list, tt_1: list, input_idx: int) -> list:
    """
    Function to rebuild the truth table according to the index
    Sample Input:
    tt_0 = [1, 0, 0, 1], tt_1 = [0, 1, 0, 0]
    input_idx = 0
    Output: [1, 0, 0, 1, 0, 1, 0, 0]
    """
    if len(tt_0) != len(tt_1):
        raise ValueError("Truth tables must have the same length")
    if 2**int(math.log2(len(tt_0))) != len(tt_0):
        raise ValueError("Truth tables must be 2 to the power of n")
    input_size = int(math.log2(len(tt_0))) + 1
    f_order = [i for i in range(input_size) if i != input_idx]
    f_order.insert(0, input_idx)
    return reorder_tt(tt_0 + tt_1, f_order, list(range(input_size)))


def bin_to_verilog_hex(tt: list) -> str:
    """
    Convert a binary array to a Verilog Hex string
    Sample Input:
    tt = [0, 0, 0, 0, 1, 1, 1, 1]
    Output: "8'h0F"
    """
    tt_len = len(tt)
    if 2**int(math.log2(tt_len)) != tt_len:
        raise ValueError("Truth tables must be 2 to the power of n")
    # Convert list of numbers into a binary string
    binary_str = ''.join(str(x) for x in tt)
    # Convert binary string into an integer
    binary_int = int(binary_str, 2)
    # Convert integer into a hex string and remove the '0x' prefix
    hex_str = hex(binary_int)[2:].upper()
    return str(tt_len) + "'h" + hex_str.zfill(int(tt_len / 4))


def bool_list_to_hex(bool_list: list) -> str:
    # Convert boolean list to binary string
    binary_str = ''.join('1' if b else '0' for b in bool_list)
    integer = int(binary_str, 2)
    hex_str = hex(integer)[2:].upper()
    return hex_str


def verilog_hex_to_bin(verilog_str: str) -> list:
    """
    Convert a Verilog Hex string to binary array
    Sample Input:
    verilog_str = "8'hAB"
    Output: [1, 0, 1, 0, 1, 0, 1, 1]
    """
    # Extract the hexadecimal number from the string
    hex_num = verilog_str.split("'h")[-1]
    list_len = int(verilog_str.split("'h")[0])
    # Convert the hexadecimal number to binary
    bin_num = bin(int(hex_num, 16))[2:]
    # Pad the binary number with zeros to the left to get 8 bits
    bin_num = bin_num.zfill(list_len)
    # Convert the binary number to an array of bits
    tt = [int(bit) for bit in bin_num]
    return tt


def same_inputs(f1_inputs: list, f2_inputs: list) -> list:
    """
    Function to find common inputs between two functions
    Sample Input:
    f1_inputs = ['a', 'b', 'c', 'd']
    f2_inputs = ['x', 'd', 'c', 'a']
    Output: ['a', 'c', 'd']
    """
    # Convert the input lists to sets
    f1_inputs_set = set(f1_inputs)
    f2_inputs_set = set(f2_inputs)
    # Get the intersection of the two sets
    same_inputs_set = f1_inputs_set & f2_inputs_set
    # Convert the set back to a list
    same_inputs_list = list(same_inputs_set)
    return same_inputs_list


def diff_inputs(f1_inputs: list, f2_inputs: list) -> tuple[list, list]:
    """
    Function to find different inputs between two functions
    Sample Input:
    f1_inputs = ['a', 'k', 'b', 'c', 'd', 'm']
    f2_inputs = ['x', 'd', 'c', 'a']
    Output: ['k', 'b', 'm'], ['x']
    """
    # Convert the input lists to sets
    f1_inputs_set = set(f1_inputs)
    f2_inputs_set = set(f2_inputs)
    # Get the difference of the two sets
    func1_diff = list(f1_inputs_set - f2_inputs_set)
    func2_diff = list(f2_inputs_set - f1_inputs_set)
    # Sort the lists in the same order as the original input lists
    func1_diff.sort(key=f1_inputs.index)
    func2_diff.sort(key=f2_inputs.index)
    return func1_diff, func2_diff


def sop2tt(sop: list) -> list:
    """
    Sample Input:
    sop = [
        '--001 1',
        '-1-01 1',
        '-10-1 1',
        '1---1 1',
        '1-00- 1',
        '11-0- 1',
        '110-- 1'
    ]
    Output: [0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0,
             1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1]
    """
    if len(sop) == 0:
        raise ValueError("Empty list")
    takeComplement = True if int(sop[0][-1]) == 0 else False
    minterms, tt = [], []
    for row in sop:
        inputs, output = row[:-1].strip(), row[-1]
        indices = [i for i, x in enumerate(inputs) if x == '-']
        for values in itertools.product('01', repeat=len(indices)):
            new_inputs = list(inputs)
            for index, value in zip(indices, values):
                new_inputs[index] = value
            minterms.append(''.join(new_inputs))
    num_inputs = len(minterms[0])
    if num_inputs == 0:
        tt.append(1 if sop == ['1'] else 0)
        return tt
    else:
        for i in range(2 ** num_inputs):
            pattern = format(i, '0' + str(num_inputs) + 'b')
            tt.append(1 if pattern in minterms else 0)
    if takeComplement:
        tt = list_complement(tt)
    return tt


def tt2sop(tt: list) -> list:
    n = int(math.log2(len(tt)))  # number of input variables
    if n == 0:
        inputs = ['']
    else:
        inputs = [bin(i)[2:].zfill(n) for i in range(2 ** n)]  # input patterns
    tt_str = [str(i) for i in tt]
    empty_str = [" "] * len(tt)
    rst_all = [''.join(t) for t in zip(inputs, empty_str, tt_str)]
    rst_1 = [s for s in rst_all if not s.endswith('0')]
    rst_0 = [s for s in rst_all if not s.endswith('1')]
    return rst_1 if len(rst_1) > 0 else rst_0
