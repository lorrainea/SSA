### Install and run PA:

```
make -f Makefile.gcc
./ssa <sequence_file> <suffix_list> <output_filename>
```

<sequence_file> is a textfile containing the input text and
<suffix_list> is a textfile containing the starting positions of the suffixes to be sorted (separated by a new line).

Two files are output: <output_filename>.ssa which contains the sparse suffix array; and <output_filename>.lcp which contains the lcp values between consecutive suffixes. Every entry is on a new line.

The program can be run using the following example, where the files out.ssa and out.lcp will be output:

```
./ssa ./data/text/text1 ./data/suffixes/suffixes1 out
```

The time complexity is O(n + b log b + (b'n/b) log b) and the extra space used is O(b), 
where n is the number of bytes in <sequence_file>, 
b is the number of positions in <suffix_list>, 
and b' is the number of suffixes with an LCP value at least 2n/b.

________________________________

