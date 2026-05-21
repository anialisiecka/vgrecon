## About the Project

This repository provides an algorithm for generating a **reconciled graph** and computing the **homology relation cardinality** (size) for GFA graphs. 

## Installation

To run `vgrecon`, you need a C++ compiler supporting C++11 or later (e.g., `g++`) and the `gfalace` tool installed in your system PATH.

1. **Clone the repository:**
   ```bash
   git clone https://github.com
   cd vgrecon
   ```

2. **Compile the core algorithm:**
   ```bash
   g++ -O3 vgrecon.cpp -o vgrecon
   ```

3. **Compile the algorithm for homology relation cardinality calculation**
   ```bash
   g++ -O3 homRel.cpp -o homrel
   ```
## Usage
### Step 1: Run the core reconciliation algorithm by passing two GFA files, a FASTA file, and the output path:
```bash
./vgrecon input1.gfa input2.gfa input.fasta output_rec.gfa
```

### Step 2: Use `gfalace` to insert trivial nodes into the `output_rec.gfa` graph obtained in Step 1, using the FASTA sequences to fill the gaps:
```bash
gfalace -g output_rec.gfa -o output.gfa --fill-gaps 2 --fasta-files input.fasta
```

### Step 3 (optional): Compute the cardinality of the homology relation induced by the reconciled graph `output.gfa` obtained in Step 2:
```bash
./homrel output.gfa
```
