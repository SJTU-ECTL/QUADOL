# QUADOL: A Quality-Driven Approximate Logic Synthesis Method Exploiting Dual-Output LUTs for Modern FPGAs

## Requirements
### ALSRAC Dependency

QUADOL uses [ALSRAC](https://github.com/SJTU-ECTL/ALSRAC) for comparative experiments.
Several end-to-end comparison scripts also utilize ALSRAC by default.
Please ensure that ALSRAC is properly compiled and deployed in your environment.

The location of the ALSRAC installation can be specified in the `scripts/env.sh` file.
This is done by modifying the path on line 13. By default, ALSRAC is expected to be located at `QUADOL/../ALSRAC`.

### Python Environment

The project requires Python 3.7 or higher. The following additional libraries are needed:

* matplotlib
* pandas
* networkx

* To install the libraries, use the following command:
   
```[bash]
pip install matplotlib pandas networkx
```

### Building the Error Simulator

Navigate to the `scripts` directory and run the `build.sh` script:

```bash
bash build.sh
```

Upon successful compilation, the error simulator will be located at `utils/build/simulator.out`.

### Using QUADOL

To use QUADOL, navigate to the `src` directory and run the `main.py` script with the appropriate arguments:

```bash
python main.py -i [Exact BLIF file] -o [Target dict] -s [Simulator] -e [Error Threshold] -m [Error Metric]
```

The `-m` (Error Metric) argument supports the following options:

* `er` - Error Rate
* `med` - Mean Error Distance
* `mred` - Mean Relative Error Distance
* `nmed` - Normalized Mean Error Distance

Please replace `[Exact BLIF file]`, `[Target dict]`, `[Simulator]`, `[Error Threshold]`, and `[Error Metric]` with your specific inputs.

For example, in the root folder of the project:
```bash
python src/main.py -i Inputs/Arith/log2_size_2018.blif -o results -s utils/build/simulator.out -e 0.001 -m mred
```

### Using QUADOL+

We also provide a framework for QUADOL+.
The sample enhanced version of QUADOL strengthens ALSRAC under the constraint of \(MRED \leq 0.001\) on the EPFL Arithmetic Benchmark.
To use QUADOL+, run the following command:

```bash
bash QUADOL-plus.sh
```
