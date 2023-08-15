### Install and run MA and PA using:

```
make -f Makefile.gcc
./ssa <sequence_file> <suffix_list> <output_filename>
```

Two output files are output, namely; <output_filename>.ssa which contains the sparse suffix array and <output_filename>.lcp which contains the lcp values between consecutive suffixes. The tool can be run using the following example, where the files out.ssa and out.lcp will be output:

```
./ssa ./data/text/text1 ./data/suffixes/suffixes1 out
```
________________________________

### Install and run SSA-LCE:

```
cmake .
make
./sa-rk <sequence_file> <output_filename> <suffix_list> 
```

Two output files are output, namely; <output_filename>.ssa which contains the sparse suffix array and <output_filename>.lcp which contains the lcp values between consecutive suffixes. The tool can be run using the following example, where the files out.ssa and out.lcp will be output:

```
./sa-rk ./data/text/text1 ./data/suffixes/suffixes1 out
```
