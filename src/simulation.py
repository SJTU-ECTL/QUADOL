import os
import subprocess
import definition
import writer
import uuid


def synthesis(design: definition.BLIF_design, config):
    # Synthesize the design and return its cost.
    return design.size - config.count(True)


def accuracy_sim(design: definition.BLIF_design, config, ID: int):
    approx = definition.approx_config(design, config)
    best_error = float('inf')
    best_implement, test_folder = None, None
    # Define the working dir
    path_str = str(uuid.uuid4())
    test_folder = os.path.join(design.test_folder, path_str)
    if os.path.exists(test_folder):
        raise ValueError(f'{test_folder} already exists')
    os.makedirs(test_folder)
    approx_file = os.path.join(test_folder, design.test_blif)
    for impl in approx.implement_idx:
        # Write BLIF files
        writer.write_BLIF(approx_file, design, approx, impl)
        # Run error simulation
        test_file = os.path.join(test_folder, 'error.txt')
        simulator = os.path.join(os.getcwd(), design.simulator)
        command = [simulator, '-a', design.origin_blif, '-b', approx_file, '-f', test_file, '-m', design.e_mode]
        process = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if process.returncode != 0:
            print(test_folder)
        # Find best implement
        with open(test_file, 'r') as file:
            error = float(file.read())
            if error < best_error:
                best_error = error
                best_implement = impl
    # Record best implement to a txt file
    if best_implement:
        implement_file = os.path.join(test_folder, 'implement.txt')
        with open(implement_file, 'w') as file:
            file.write(str(best_implement))
    # Return Format:
    # 0. Config of the best design
    # 1. Simulation results: Measured Error
    # 2. Which LUT pair is Changed in this round
    # 3. The best implementation for this config
    return config, best_error, ID, best_implement


def find_best_design(rst_lst: list, skip_sets: set, e_max: float) -> list:
    # Sort the designs by accuracy, then reverse so that it's in ascending order.
    sorted_results = sorted(rst_lst, key=lambda x: x[1], reverse=False)
    if len(sorted_results) == 0:
        raise ValueError("Items not found in Design List")
    # Add config IDs with large error to skipped set.
    for i in sorted_results:
        if i[1] > e_max:
            skip_sets.add(i[2])
    return list(sorted_results[0])


def find_best_implement(rst_lst: list, design: definition.BLIF_design):
    for rst in rst_lst:
        if len(rst[3]) != 1:
            raise ValueError("Invalid implementation")
        design.best_implement[rst[2]] = rst[3][0]
