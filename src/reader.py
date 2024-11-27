import util
import definition
from enum import Enum


class READ_MODE(Enum):
    INPUT = 'INPUT'
    OUTPUT = 'OUTPUT'
    LUTS_TT = 'LUTS_TT'
    LUTS_INFO = 'LUTS_INFO'
    NO_MODE = 'NONE'


def read_BLIF(filename: str) -> definition.BLIF_design:
    with open(filename, 'r') as file:
        lines = file.readlines()

    design_name = ''
    luts, luts_input, luts_output = [], [], []
    read_mode, tmp_line, tmp_lines = READ_MODE.NO_MODE, '', []

    for line in lines:
        # Remove trailing newline
        line = line.strip()

        # Model Start
        if line.startswith('.model'):
            design_name = line.split()[1]
            continue
        # Model End
        if line.startswith('.end'):
            if len(tmp_lines) != 0:
                tt = util.sop2tt(tmp_lines)
                luts[-1].tt = tt
            break

        # Model IOs
        if line.startswith(('.inputs', '.outputs')):
            read_mode = READ_MODE.INPUT if line.startswith('.inputs') else READ_MODE.OUTPUT
        if read_mode in (READ_MODE.INPUT, READ_MODE.OUTPUT):
            tmp_line += line
            if tmp_line.endswith('\\'):
                tmp_line = tmp_line.split('\\', 1)[0]
                continue
            if read_mode == READ_MODE.INPUT:
                luts_input = tmp_line.split()[1:]
            elif read_mode == READ_MODE.OUTPUT:
                luts_output = tmp_line.split()[1:]
            read_mode, tmp_line = READ_MODE.NO_MODE, ''
            continue

        # Model LUTs
        if line.startswith('.names'):
            read_mode = READ_MODE.LUTS_INFO
        if read_mode == READ_MODE.LUTS_INFO:
            tmp_line += line
            if tmp_line.endswith('\\'):
                tmp_line = tmp_line.split('\\', 1)[0]
                continue
            if len(tmp_lines) != 0:
                tt = util.sop2tt(tmp_lines)
                luts[-1].tt = tt
                tmp_lines = []
            lut_ports = tmp_line.split()[1:]
            lut_name = lut_ports[-1]
            lut = definition.LUT(lut_name, lut_ports[:-1])
            luts.append(lut)
            read_mode, tmp_line = READ_MODE.LUTS_TT, ''
            continue
        if read_mode == READ_MODE.LUTS_TT:
            tmp_lines.append(line)

    return definition.BLIF_design(design_name, luts_input, luts_output, luts)
