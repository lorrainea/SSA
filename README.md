Run using:

./ssa <sequence_file> <suffix_list> <output_filename>

Two output files are output, namely; <output_filename>.ssa which contains the sparse suffix array and <output_filename>.lcp which contains the lcp values between consecutive suffixes. The tool can be run using the following example, where the files out.ssa and out.lcp will be output:

./ssa text suffixes out
