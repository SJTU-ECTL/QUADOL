import sys
import os
import re
import shutil

script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)

overall_dir = os.path.join(script_dir, '..')
results_dir = os.path.join(overall_dir, 'results')
ALSRAC_dir = os.path.join(results_dir, 'ALSRAC')


# Read ALSRAC Log file & Return the best results in the file
def single_ALSRAC_log(filename: str):
    # Define patterns
    original_size_pattern = r"Original circuit: size = (\d+)"
    round_pattern = r"--------------- round (\d+) ---------------"
    error_pattern = r"current error = ([-+]?\d*\.?\d+(?:e[-+]?\d+)?)"
    size_pattern = r"size = (\d+)"

    with open(filename, 'r') as file:
        content = file.read()
        # Find all matches in the file
        original_size = re.findall(original_size_pattern, content)
        if len(original_size) == 0:
            return float('inf'), float('inf')
        original_size = int(original_size[0])
        rounds = re.findall(round_pattern, content)
        errors = re.findall(error_pattern, content)
        errors_float = [float(num) for num in errors]
        sizes = re.findall(size_pattern, content)
        sizes = sizes[1:]
        # Zip together the results from each pattern
        results = list(zip(rounds, errors_float, sizes))
    filtered_tuples = results[:-1]
    # if len(filtered_tuples) == 0:
    #     return float('inf'), float('inf')
    min_size = min(int(t[2]) for t in filtered_tuples)
    return original_size, min_size


# First Step: Find the best results from ALSRAC
def restruct_ALSRAC_results():
    example_list = dict()
    for filename in os.listdir(ALSRAC_dir):
        full_path = os.path.join(ALSRAC_dir, filename)
        if os.path.isdir(full_path):
            continue
        term = os.path.splitext(filename)[0]
        term_lst = term.split('_')
        if term_lst[0] in example_list:
            example_list[term_lst[0]].append(full_path)
        else:
            example_list[term_lst[0]] = [full_path]
    for key, values in example_list.items():
        rst = None
        for value in values:
            _, min_size = single_ALSRAC_log(value)
            if rst is None or min_size < rst[1]:
                rst = (_, min_size, value)

        filename = os.path.basename(rst[2])
        term = os.path.splitext(filename)[0]

        old_folder_name = os.path.join(ALSRAC_dir, term)
        if os.path.isdir(old_folder_name):
            old_file_name = rst[2]
            new_folder_name = os.path.join(ALSRAC_dir, key)
            new_file_name = os.path.join(ALSRAC_dir, f'{key}.log')

            os.rename(old_folder_name, new_folder_name)
            os.rename(old_file_name, new_file_name)

    for filename in os.listdir(ALSRAC_dir):
        term = os.path.splitext(filename)[0]
        term_lst = term.split('_')
        if len(term_lst) > 1:
            full_path = os.path.join(ALSRAC_dir, filename)
            if os.path.isdir(full_path):
                shutil.rmtree(full_path)
            if os.path.isfile(full_path):
                os.remove(full_path)


# Second Step: Remove half of the results in ALSRAC
# Reason: These results cannot defeate ALSRAC after LUT Merging
def remove_ALSRAC_redundancies():
    error_th = float(sys.argv[1])
    for filename in os.listdir(ALSRAC_dir):
        full_path = os.path.join(ALSRAC_dir, filename)
        if os.path.isfile(full_path):
            continue
        log_file = os.path.join(ALSRAC_dir, f'{filename}.log')
        o_size, m_size = single_ALSRAC_log(log_file)
        if abs(o_size - m_size) < 40:
            continue
        m_size = min(o_size, m_size)
        mid_size_1 = (o_size + m_size) // 2
        mid_size_2 = m_size * 2
        mid_size = min(mid_size_1, mid_size_2)
        for blif in os.listdir(full_path):
            term_lst = blif.split('_')
            if int(term_lst[-2]) > mid_size:
                os.remove(os.path.join(full_path, blif))
            if float(term_lst[-3]) > error_th:
                os.remove(os.path.join(full_path, blif))


def main():
    restruct_ALSRAC_results()
    remove_ALSRAC_redundancies()


if __name__ == '__main__':
    main()
