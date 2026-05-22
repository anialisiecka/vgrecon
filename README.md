## About the Project

This repository provides an algorithm for generating a **reconciled graph** and computing the **homology relation cardinality** (size) for GFA graphs. 

## Installation

To run `vgrecon`, you need a C++ compiler supporting C++11 or later (e.g., `g++`) and the [gfalace](https://github.com/pangenome/gfalace) tool installed in your system PATH.

1. **Clone the repository:**
   ```bash
   git clone https://github.com
   cd vgrecon
   ```

2. **Compile the core algorithm and homology calculation algorithm:**
   ```bash
   g++ -O3 rec.cpp -o rec
   g++ -O3 homRel.cpp -o homrel
   ```

3. **Make the main pipeline script executable:**
   ```bash
   chmod +x vgrecon
   ```

## Usage

Ensure that the compiled `rec` binary remains in the same directory as the `vgrecon` script, then run:

```bash
./vgrecon --gfa1 graph1.gfa --gfa2 graph2.gfa --fasta sequences.fa --output reconciled.gfa
```

Additionally, if you want to measure the cardinality (size) of the homology relation induced by any variation graph (e.g., a reconciled graph), you can run the following:

```bash
./homrel graph.gfa
```
