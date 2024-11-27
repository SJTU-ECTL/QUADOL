import definition
import util
import writer


FLOAT_INF = float('inf')


def HD_4tt(tt: list, forbidden_combine: list = None) -> definition.LUT62:
    # Assume: tt = [A, B, C, D]
    if len(tt) != 4:
        raise ValueError('Truth table must have 4 sub-truth tables')

    # Calculate Type22 Min HD
    # Assume: tt = [A, B, C, D]
    # [(A, B), (A, C), (A, D), (B, C), (B, D), (C, D)]
    tt_list_22 = util.list_combinations(tt)
    # Calculate the HD between cofactors
    HD_list_22 = [util.HD_2tt(i[0], i[1]) for i in tt_list_22]
    # Format of `rst_22`: [(idx, HD), (idx, HD)]
    rst_22 = util.least_elements(HD_list_22, 2)
    HD_22 = sum(x[1] for x in rst_22)
    # Skip Forbidden Structures
    if {x[0] for x in rst_22} in forbidden_combine:
        rst_22 = util.least_elements(HD_list_22, 3)
        del rst_22[1]
        HD_22 = sum(x[1] for x in rst_22)

    # Define a dictionary to map the set to the corresponding assignments
    rst_22_set = frozenset(x[0] for x in rst_22)
    mapping_type22 = {
        frozenset({0, 5}): ([tt[0], tt[1]], [tt[2], tt[3]], definition.MUX_TYPE.NO_MUX),
        frozenset({1, 4}): ([tt[0], tt[2]], [tt[1], tt[3]], definition.MUX_TYPE.MUX_0),
    }

    # Calculate Type 31 Min HD
    # [(A, B, C), (A, B, D), (A, C, D), (B, C, D)]
    tt_list_31 = util.list_combinations(tt, 3)
    # Calculate the HD between cofactors
    bit_maj_list = [util.bit_majority(i[0], i[1], i[2]) for i in tt_list_31]
    HD_list_31 = [util.HD_3tt(i[0], i[1], i[2]) for i in tt_list_31]
    # Format of `rst_31`: [(idx, HD)]
    rst_31 = util.least_elements(HD_list_31)
    if len(rst_31) != 1:
        raise ValueError('Function `least_elements` must return a single element')
    HD_31 = rst_31[0][1]

    # Define a dictionary to map isolated truth table to the corresponding assignments
    isolated_tt_ID = 3 - rst_31[0][0]
    mapping_type31 = {
        0: (tt[0], bit_maj_list[3], definition.MUX_TYPE.MUX_A, False),
        1: (tt[1], bit_maj_list[2], definition.MUX_TYPE.MUX_A, True),
        2: (bit_maj_list[1], tt[2], definition.MUX_TYPE.MUX_B, False),
        3: (bit_maj_list[0], tt[3], definition.MUX_TYPE.MUX_B, True),
    }

    rst = definition.LUT62()
    if HD_31 <= HD_22:
        # Graph Type #1: 3-1
        # For example: A-B-C, D
        rst.HD = HD_31
        A_tt, B_tt, rst.mux_type, rst.mux_inv = mapping_type31[isolated_tt_ID]
        rst.A_tt.append(A_tt)
        rst.B_tt.append(B_tt)
    else:
        # Graph Type #2: 2-2
        # For example: A-B, C-D
        rst.HD = HD_22
        rst.A_tt.extend(mapping_type22[rst_22_set][0])
        rst.B_tt.extend(mapping_type22[rst_22_set][1])
        rst.mux_type = mapping_type22[rst_22_set][2]
    return rst


def MinHD_helper(tt1: list, tt2: list, idx: int, forbidden_combine: list = None) -> definition.LUT62:
    if len(tt1) != len(tt2):
        raise ValueError("LUTs must have same size")
    # Decompose the truth table according to `idx`
    # Assume: tt1 = 'G', tt2 = F, tt = [('G', True), ('G', False), ('F', True), ('F', False)]
    tt_p = [util.decompose_tt(cofactor, idx, flag)
            for cofactor in [tt1, tt2] for flag in [True, False]]
    rst_p = HD_4tt(tt_p, forbidden_combine)
    # Take the complement
    tt1_inv, tt2_inv = util.list_complement(tt1), util.list_complement(tt2)
    tt1_comp = [util.decompose_tt(cofactor, idx, flag)
                for cofactor in [tt1_inv, tt2] for flag in [True, False]]
    tt2_comp = [util.decompose_tt(cofactor, idx, flag)
                for cofactor in [tt1, tt2_inv] for flag in [True, False]]
    rst_1_comp = HD_4tt(tt1_comp, forbidden_combine)
    rst_2_comp = HD_4tt(tt2_comp, forbidden_combine)

    min_HD = min(rst_p.HD, rst_1_comp.HD, rst_2_comp.HD)
    if min_HD == rst_p.HD:
        rst = rst_p
    elif min_HD == rst_1_comp.HD:
        rst = rst_1_comp
        rst.O1_inv = True
    elif min_HD == rst_2_comp.HD:
        rst = rst_2_comp
        rst.O2_inv = True
    else:
        raise ValueError("There must exist a min HD")
    return rst


def MinHD_LUT6_LUT6_S6(tt1: list, tt2: list, input_list: list) -> definition.LUT62:
    if len(tt1) != 2**6 or len(tt2) != 2**6:
        raise ValueError("LUTs must have pre-defined size")
    HD, RST = [], []
    for idx in range(6):
        # Sum the least two HD
        rst = MinHD_helper(tt1, tt2, idx, [{2, 3}])
        HD.append(rst.HD)
        RST.append(rst)
    HD_least_element = util.least_elements(HD)
    rst = RST[HD_least_element[0][0]]
    if rst.mux_type == definition.MUX_TYPE.MUX_E:
        raise ValueError("MUX Type Error")
    rst.mux_input = input_list.pop(HD_least_element[0][0])
    rst.port_list = input_list
    return rst


def MinHD_LUT6_LUT6_S5(tt1: list, tt2: list, input_list: list) -> definition.LUT62:
    if len(tt1) != 2**6 or len(tt2) != 2**6:
        raise ValueError("LUTs must have pre-defined size")
    rst = MinHD_helper(tt1, tt2, 5, [{2, 3}, {1, 4}])
    if rst.mux_type == definition.MUX_TYPE.MUX_A:
        rst.mux_input = input_list[5]
    elif rst.mux_type == definition.MUX_TYPE.MUX_B:
        rst.mux_input = input_list[6]
    elif rst.mux_type == definition.MUX_TYPE.NO_MUX:
        rst.mux_input = input_list[5]
    else:
        raise ValueError("MUX Type Error")
    rst.port_list = input_list[:5]
    return rst


def MinHD_LUT5_LUT6_S5(tt1: list, tt2: list, input_list: list) -> definition.LUT62:
    if len(tt1) != 2**5 or len(tt2) != 2**6:
        raise ValueError("LUTs must have pre-defined size")
    tt = [tt1, util.decompose_tt(tt2, 5, True),
          util.decompose_tt(tt2, 5, False)]
    hd_list = [util.HD_2tt(tt[0], tt[1]), util.HD_2tt(tt[0], tt[2]), util.HD_2tt(tt[1], tt[2])]
    hd_min = util.least_elements(hd_list)
    # Define a dictionary to map the value of hd_min[0][0] to the corresponding assignments
    mapping = {
        0: ([tt[0], tt[1]], [tt[2]], definition.MUX_TYPE.MUX_B, False),
        1: ([tt[0], tt[2]], [tt[1]], definition.MUX_TYPE.MUX_B, True),
        2: ([tt[0]], [tt[1], tt[2]], definition.MUX_TYPE.NO_MUX, False),
    }

    rst = definition.LUT62()
    rst.HD = hd_min[0][1]
    rst.mux_input = input_list[5]
    rst.port_list = input_list[:5]
    rst.A_tt.extend(mapping[hd_min[0][0]][0])
    rst.B_tt.extend(mapping[hd_min[0][0]][1])
    rst.mux_type = mapping[hd_min[0][0]][2]
    rst.mux_inv = mapping[hd_min[0][0]][3]
    return rst


def MinHD(lut1: definition.LUT, lut2: definition.LUT) -> definition.LUT62:
    # Obtain original LUT input ports & truth table
    f1_input, f2_input = lut1.port_list, lut2.port_list
    f1_tt, f2_tt = lut1.tt, lut2.tt
    # Analyze same & different ports between lut1 & lut2
    same_input_list = util.same_inputs(f1_input, f2_input)

    def accelerator_1():
        f1_diff, f2_diff = util.diff_inputs(f1_input, f2_input)
        # New port list for lut1 & lut2
        f1_new_input_ = same_input_list + f1_diff
        f2_new_input_ = same_input_list + f2_diff
        new_input_ = same_input_list + f1_diff + f2_diff
        return f1_new_input_, f2_new_input_, new_input_

    def accelerator_2(f1_new_input_, f2_new_input_):
        # Reorder truth table due to new port order
        f1_new_tt_1 = util.reorder_tt(f1_tt, f1_input, f1_new_input_)
        f2_new_tt_1 = util.reorder_tt(f2_tt, f2_input, f2_new_input_)
        return f1_new_tt_1, f2_new_tt_1

    # Calculate Hamming Distance (HD) between two LUTs

    # Type #1: LUT6 + LUT6 => LUT6_2 (Similarity 6)
    if lut1.size == 6 and lut2.size == 6 and len(same_input_list) == 6:
        f1_new_input, f2_new_input, _ = accelerator_1()
        if f1_new_input != f2_new_input:
            raise ValueError("Inputs do not match")
        # Calculate Min Hamming Distance
        f1_new_tt, f2_new_tt = accelerator_2(f1_new_input, f2_new_input)
        rst = MinHD_LUT6_LUT6_S6(f1_new_tt, f2_new_tt, f1_new_input)
        rst.lut1, rst.lut2 = lut1, lut2
        return rst

    # Type #2: LUT6 + LUT6 => LUT6_2 (Similarity 5)
    if lut1.size == 6 and lut2.size == 6 and len(same_input_list) == 5:
        f1_new_input, f2_new_input, new_input = accelerator_1()
        if f1_new_input[:5] != f2_new_input[:5]:
            raise ValueError("Inputs do not match")
        # Calculate Min Hamming Distance
        f1_new_tt, f2_new_tt = accelerator_2(f1_new_input, f2_new_input)
        rst = MinHD_LUT6_LUT6_S5(f1_new_tt, f2_new_tt, new_input)
        rst.lut1, rst.lut2 = lut1, lut2
        return rst

    # Type #3: LUT5 + LUT6 => LUT6_2 (Similarity 5)
    if lut1.size == 5 and lut2.size == 6 and len(same_input_list) == 5:
        f1_new_input, f2_new_input, _ = accelerator_1()
        f1_new_tt = util.reorder_tt(f1_tt, f1_input, f1_new_input)
        f2_new_tt = util.reorder_tt(f2_tt, f2_input, f2_new_input)
        rst = MinHD_LUT5_LUT6_S5(f1_new_tt, f2_new_tt, f2_new_input)
        rst.lut1, rst.lut2 = lut1, lut2
        return rst
    if lut1.size == 6 and lut2.size == 5 and len(same_input_list) == 5:
        f1_new_input, f2_new_input, _ = accelerator_1()
        f1_new_tt = util.reorder_tt(f1_tt, f1_input, f1_new_input)
        f2_new_tt = util.reorder_tt(f2_tt, f2_input, f2_new_input)
        rst = MinHD_LUT5_LUT6_S5(f2_new_tt, f1_new_tt, f1_new_input)
        rst.lut1, rst.lut2 = lut2, lut1
        return rst

    return definition.LUT62()


class SimilarPair:
    def __init__(self, lut1: definition.LUT, lut2: definition.LUT):
        # Calculate the min Hamming Distance
        self.MinHD = MinHD(lut1, lut2)
        self.value = self.MinHD.HD
        if self.value != float('inf'):
            self.MinHD.implement = writer.LUT62_BLIF(self.MinHD)

    def __str__(self):
        return str(self.value)
