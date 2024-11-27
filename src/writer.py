import itertools
import definition
import util
import os
import shutil


def append_inverse(templet: list, lut: str, lut_name: str, inv: bool):
    inverse = definition.LUT(lut_name, [lut])
    inverse.tt = [1, 0] if inv else [0, 1]
    templet.append(inverse)


def define_mux(lut62: definition.LUT62, templet: list,
               lut_a: str, lut_b: str, lut_M: str):
    if lut62.mux_type == definition.MUX_TYPE.NO_MUX:
        append_inverse(templet, lut_a, lut62.lut1.name, lut62.O1_inv)
        append_inverse(templet, lut_b, lut62.lut2.name, lut62.O2_inv)
    elif lut62.mux_type in {definition.MUX_TYPE.MUX_A, definition.MUX_TYPE.MUX_B, definition.MUX_TYPE.MUX_0}:
        mux = definition.LUT(lut_M, [lut_a, lut_b, lut62.mux_input])
        if lut62.mux_inv:
            mux.tt = [0, 0, 0, 1, 1, 0, 1, 1]
        else:
            mux.tt = [0, 0, 1, 0, 0, 1, 1, 1]
        templet.append(mux)
        if lut62.mux_type == definition.MUX_TYPE.MUX_A:
            append_inverse(templet, lut_M, lut62.lut1.name, lut62.O1_inv)
            append_inverse(templet, lut_b, lut62.lut2.name, lut62.O2_inv)
        elif lut62.mux_type == definition.MUX_TYPE.MUX_B:
            append_inverse(templet, lut_M, lut62.lut2.name, lut62.O2_inv)
            append_inverse(templet, lut_a, lut62.lut1.name, lut62.O1_inv)
        elif lut62.mux_type == definition.MUX_TYPE.MUX_0:
            append_inverse(templet, lut_M, lut62.lut1.name, lut62.O1_inv)
            append_inverse(templet, lut_M, lut62.lut2.name, lut62.O2_inv)


def LUT62_BLIF(lut62: definition.LUT62) -> list:
    rst, templet = [], []
    # Define wire names
    lut_name = lut62.lut1.name + '_' + lut62.lut2.name
    lut_a_str = lut_name + '_A'
    lut_b_str = lut_name + '_B'
    # MUX output
    lut_M_str = lut_name + '_M'

    # Define MUX Connection
    define_mux(lut62, templet, lut_a_str, lut_b_str, lut_M_str)

    tt_combine = list(itertools.product(lut62.A_tt, lut62.B_tt))
    for tt in tt_combine:
        lut_a = definition.LUT(lut_a_str, lut62.port_list)
        lut_b = definition.LUT(lut_b_str, lut62.port_list)
        lut_a.tt, lut_b.tt = tt[0], tt[1]
        tmp_implement = [lut_a, lut_b]
        tmp_implement.extend(templet)
        rst.append(tmp_implement)
    return rst


def approx_lut_list(design: definition.BLIF_design, approx: definition.approx_config, impl) -> list:
    rst = design.luts
    removed_luts, used_luts = set(), []
    for idx, impl_ in zip(approx.config_idx, impl):
        removed_luts.add(design.lut62[idx].lut1)
        removed_luts.add(design.lut62[idx].lut2)
        used_luts.extend(design.lut62[idx].implement[impl_])
    rst_lst = [x for x in rst if x not in removed_luts]
    rst_lst.extend(used_luts)
    return rst_lst


def write_BLIF(filename: str, design: definition.BLIF_design, approx: definition.approx_config, impl: list):
    with open(filename, 'w') as file:
        file.write(f'.model {design.name}\n')
        file.write('.inputs ')
        for port in design.inputs:
            file.write(f'{port} ')
        file.write('\n.outputs ')
        for outputs in design.outputs:
            file.write(f'{outputs} ')
        for luts in approx_lut_list(design, approx, impl):
            file.write(f'\n.names ')
            for port in luts.port_list:
                file.write(f'{port} ')
            file.write(f'{luts.name}')
            for sop_str in util.tt2sop(luts.tt):
                file.write(f'\n{sop_str}')
        file.write('\n.end\n')


def initial_test_folder(design: definition.BLIF_design):
    if os.path.exists(design.test_folder) and os.path.isdir(design.test_folder):
        shutil.rmtree(design.test_folder)
    os.makedirs(design.test_folder)
