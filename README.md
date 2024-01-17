### Install and run MA and PA:

```
cd MA (or cd PA)
make -f Makefile.gcc
./ssa <sequence_file> <suffix_list> <output_filename>
```

<sequence_file> is a textfile containing the input text and
<suffix_list> is a textfile containing the starting positions of the suffixes to be sorted (separated with a new line).

Two files are output: <output_filename>.ssa which contains the sparse suffix array; and <output_filename>.lcp which contains the lcp values between consecutive suffixes. Every entry is on a new line.

The program can be run using the following example, where the files out.ssa and out.lcp will be output:

```
./ssa ./data/text/text1 ./data/suffixes/suffixes1 out
```
When publishing work that is based on the results from SSA please cite:
```
Lorraine A. K. Ayad, Grigorios Loukides, Solon P. Pissis, Hilde Verbeek:
Sparse Suffix and LCP Array: Simple, Direct, Small, and Fast. 
LATIN 2024
```
________________________________

### Install and run SSA-LCE:

```
cd SSA-LCE
cmake .
make
./sa-rk <sequence_file> <output_filename> <suffix_list> 
```

<sequence_file> is a textfile containing the input text and
<suffix_list> is a textfile containing the starting positions of the suffixes to be sorted (separated with a new line).

Two files are output: <output_filename>.ssa which contains the sparse suffix array; and <output_filename>.lcp which contains the lcp values between consecutive suffixes. Every entry is on a new line.

The program can be run using the following example, where the files out.ssa and out.lcp will be output:

```
./sa-rk ./data/text/text1 out ./data/suffixes/suffixes1
```
