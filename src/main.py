import argparse

import binarySearch
import reader
import distanceCal
import writer


def parse_arg() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Approximate LUT Merging')
    parser.add_argument('-i', type=str, help='Exact BLIF file', required=False)
    parser.add_argument('-a', type=str, help='Approximate BLIF file', required=True)
    parser.add_argument('-o', type=str, help='Target dict', required=True)
    parser.add_argument('-s', type=str, help='Simulator', required=True)
    parser.add_argument('-e', type=float, help='Error Threshold', required=True)
    parser.add_argument('-m', type=str, help='Error Metric', required=True)

    args = parser.parse_args()
    return args


def main():
    args = parse_arg()
    design = reader.read_BLIF(args.a)
    # Define working space
    if args.i:
        design.origin_blif = args.i
    else:
        design.origin_blif = args.a
    design.approx_blif = args.a
    design.test_folder = args.o
    # Define error constraints
    design.simulator = args.s
    design.e_th = args.e
    design.e_mode = args.m

    # Build working space
    writer.initial_test_folder(design)

    lut62 = []
    for i in range(len(design.luts)):
        for j in range(i + 1, len(design.luts)):
            # Calculate the min HD between node[i] & node[j]
            pair = distanceCal.SimilarPair(design.luts[i], design.luts[j])
            if pair.value == float('inf'):
                continue
            lut62.append(pair.MinHD)
    design.updateLUT62(lut62)
    # Binary Search Algorithm:
    rst = binarySearch.binary_search(design)
    # Format: Original Size, Error, Approximate Size
    print(design.size, rst[0], rst[1], sep=',')


if __name__ == "__main__":
    main()
