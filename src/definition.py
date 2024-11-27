import networkx as nx
from enum import Enum
from itertools import product


class MUX_TYPE(Enum):
    NO_MUX = 'NO_MUX'
    MUX_A = 'MUX_A'
    MUX_B = 'MUX_B'
    MUX_0 = 'MUX_DUPLICATE'
    MUX_E = 'ERROR'


class CONST(Enum):
    CONST = 'CONST'


class LUT:
    def __init__(self, name: str, port_list: list):
        self.name = name
        self.port_list = port_list
        self.tt = None
        self.size = len(self.port_list)


class LUT62:
    def __init__(self):
        self.lut1, self.lut2 = None, None
        self.HD = float('inf')
        self.port_list, self.mux_input = None, None
        self.mux_type, self.mux_inv = MUX_TYPE.MUX_E, False
        self.A_tt, self.B_tt = [], []
        self.O1_inv, self.O2_inv = False, False
        self.implement = []


class BLIF_design:
    def __init__(self, name: str, inputs: list, outputs: list, luts: list):
        self.name = name
        self.inputs = inputs
        self.outputs = outputs
        self.luts = luts
        self.lut_dic = dict()
        for lut in self.luts:
            self.lut_dic[lut.name] = lut
        # Ignore Constant Signals
        self.size = sum(1 for lut in self.luts if lut.size != 0)
        # Ignore Direct Copied Signals
        wire_num = sum(1 for lut in self.luts if lut.size == 1 and lut.tt == [0, 1])
        self.size -= wire_num
        self.implement_num = None
        self.best_implement = None
        self.lut62 = None
        # Define the working space
        self.origin_blif = None
        self.approx_blif = None
        self.test_folder = None
        self.test_blif = f'{self.name}.blif'
        # Define Error Simulation Constraints
        self.simulator = None
        self.e_th = None
        self.e_mode = None

    def validation(self, config: list):
        test_set = set()
        for d, l in zip(config, self.lut62):
            if d:
                if l.lut1 in test_set or l.lut2 in test_set:
                    return False
                test_set.add(l.lut1)
                test_set.add(l.lut2)
        return True

    def updateLUT62(self, lut62: list):
        self.lut62 = lut62
        self.implement_num = [len(x.implement) for x in self.lut62]
        self.best_implement = [None] * len(self.lut62)


class approx_config:
    def __init__(self, design: BLIF_design, config: list):
        self.config_idx = [idx for idx in range(len(config)) if config[idx]]
        self.impl_id = []
        for idx in self.config_idx:
            if design.best_implement[idx] is None:
                IDs = range(design.implement_num[idx])
            else:
                IDs = [design.best_implement[idx]]
            self.impl_id.append(IDs)
        self.implement_idx = [list(x) for x in list(product(*self.impl_id))]
