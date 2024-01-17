### Install and run MA:

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

The time complexity is  <b> O(n log b) </b> and the extra space used is <b> O(b) </b>, 
where <b> n </b> is the number of bytes in <sequence_file> and <b> b </b> is the number of positions in <suffix_list>.

MA is a space-efficient simulation of the Monte Carlo algorithm, for constructing the sparse suffix tree, presented in:

```
Tomohiro I, Juha Kärkkäinen, Dominik Kempa: Faster Sparse Suffix Sorting. STACS 2014: 386-396
```
MA rather constructs the sparse suffix and LCP arrays <b>directly</b>.
________________________________

