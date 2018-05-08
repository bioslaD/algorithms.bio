# Guide to solve the NW and SW challenges

## Working solutions in C
To get an idea of how the real world applications in C implemented these algorithms, one can visit this page to run the tools online:
https://www.ebi.ac.uk/Tools/sss/
The relevant tools are GGSEARCH (optimal global-global alignment searches using the Needleman-Wunsch algorithm), and SSEARCH (local alignment search tool using the Smith-Waterman algorithm).

One can also clone this github repo to look at the C source code: https://github.com/wrpearson/fasta36
The relevant files, but not limited to, are `src/global_sse2.*` for GGSEARCH and `src/smith_waterman_sse2.*`. It will be OK to extract and translate these file to D. However, please remeber to make the translated code more idiomatic D than C-ish D.
The C code uses SSE2, it is not a requirement that D code also utilizes SSE2, but it will be nice if you can do it :)

To get a sense of how GGSEARCH and SSEARCH work, use the above sequence and run the tools on EBI page (Tools/SSS). If you can compile and run `fasta36` command locally on your computer, use the following additional files as databases for respective commands below. See information about input files
in the "Expected input and output" section.


<a name="ggsearch"></a>GGSEARCH command:
```
ggsearch36 -T 8 -n -r +5/-4 -f -14 -g -4 -E "5.0 -1.0" -F 0.0 -b 10 -d 10 -m "F9 globalout.m9" -z 1  input.txt miniglobal.fasta
# Runs with 8 threads (-T 8) for DNA sequence (-n)
# Score (-r) for match is +5 and mismatch is -4 
# First gap gets penalty (-f) is -14, extending the gap gets -4. 
# -E 
# In essence:
# -r, -f, and -g are parameters for calculating and selecting alignment with highest score when considering two sequences.
# -E, -F, -b, and -d form constraints for how many hits to give back from database for a query sequence.
```
SSEARCH command:

```
ssearch36 -T 8 -n -r +5/-4 -f -14 -g -4 -E "5.0 -1.0" -F 0.0 -b 10 -d 10 -m "F9 localout.m9" -z 1  input.txt minilocal.fasta
```
Below are explainations for the commanline parameters:
```
  % ssearch36 -help
USAGE
 ssearch36 [-options] query_file library_file
 "@" query_file uses stdin; query_file:begin-end sets subset range
 library file formats: 0:FASTA; 1:GenBankFF; 3:EMBL_FF; 7:FASTQ; 10:subset; 12:NCBI blastdbcmd;
   alternate library formats: "library_file 7" for 7:FASTQ

DESCRIPTION
 SSEARCH performs a Smith-Waterman search
 version: 36.3.8g Feb, 2018

OPTIONS (options must preceed query_file library_file)
 -3   compare forward strand only
 -a   show complete Query/Sbjct sequences in alignment
 -b:  high scores reported (limited by -E by default);
      =<int> forces <int> results;
 -C:  length of the query/sbjct name in alignments
 -d:  number of alignments shown (limited by -E by default)
 -D   enable debugging output
 -e:  expand_script to extend hits
 -E:  E()-value,E()-repeat threshold
 -f:  gap-open penalty
 -F:  min E()-value displayed
 -g:  gap-extension penalty
 -h   help - show options, arguments
 -H   show histogram
 -i   search with reverse-complement
 -I   interactive mode
 -k:  number of shuffles
 -l:  FASTLIBS abbreviation file
 -L   long library descriptions
 -m:  Output/alignment format;
      0 - standard ":. " alignment; 1 - " xX"; 2 - ".MS.."; 3 - separate >fasta entries;
      4 - "---" alignment map; 5 - 0+4; 6 - <html>;
      8 - BLAST tabular; 8C commented BLAST tabular; 8CC BLAST tab CIGAR, 8CD BLAST tab CIGAR ext; 8CB BLAST tab BTOP
      B - BLAST Query/Sbjct alignments; BB - complete BLAST output;
      9 - FASTA tabular; 9c - FASTA tabular encoded; 9C FASTA tabular CIGAR encoded; 9B FASTA tabular BTOP encoded
     10 - parseable key:value; 11 - lav for LALIGN;
      A - aligned residue score
      F - 'F0,6,9c out_file' - alternate output formats to files;
 -M:  filter on library sequence length
 -n   DNA/RNA query
 -N:  max library length before overlapping
 -o:  offset coordinates of query/subject
 -O:  write results to file
 -p   protein query
 -P:  PSSM file
 -q   quiet [default] -- do not prompt
 -Q   quiet [default] -- do not prompt
 -r:  [+0/0]  +match/-mismatch for DNA/RNA
 -R:  raw score file
 -s:  Scoring matrix: (protein)
      BL50, BP62 (sets -f -11 -g -1); P250, OPT5, VT200,
      VT160, P120, VT120, BL80, VT80, MD40, VT40, MD20, VT20, MD10, VT10;
      scoring matrix file name; -s ?BL50 adjusts matrix for short queries;
 -S   filter lowercase (seg) residues
 -T:  max threads/workers
 -U   RNA query
 -v:  shuffle window size
 -V:  annotation characters in query/library for aligments
 -w:  width of alignment display
 -W:  alignment context length (surrounding unaligned sequence)
 -X:  Extended options
 -z:  Statistics estimation method:
      1 - regression; -1 - no stats.; 0 - no scaling; 2 - Maximum Likelihood Est.;
      3 - Altschul/Gish; 4 - iter. regress.; 5 - regress w/variance;
      6 - MLE with comp. adj.;
     11 - 16 - estimates from shuffled library sequences;
     21 - 26 - E2()-stats from shuffled high-scoring sequences;
 -Z:  [library entries] database size for E()-value

```
## Expected input 
Use this input as query sequence for testing.

To be closer to real-world usage, and still simple, let's keep the input format as FASTA as shown in input.txt.

```
# This file is used for testing GGSEARCH, SSEARCH and your own programs
cat input.txt
>EB01
GTTAGACTCATTGGTAGAGCGCGTGCTTAGCATGTACGAGGtCCGGGTTCAATCCGGGCA
CCTCCA
>EB02
GTTAGTACGATTGGTAGAGCGCGTAATGCGAGTGTACGAGGtCCGGGTTCAATCCATGAT
CCTCCA

```
And search databases below:
```
# Data for use with GGSEARCH or your own equivalent
cat miniglobal.fasta

>KH527 KH780527.1
GGGGGTGTAGCTCAGTGGTAGAGCGCGTGCTTAGCATGTACGAGGTCCCGGGTTCAATCC
CCGGCACCTCCA
>KH528 KH780528.1
GGGGGTGTAGCTCAGTGGTAGAGCGCGTGCTTAGCATGCACGAGGCCCCGGGTTCAATCC
CCGGCACCTCCA
>KH523 KH780523.1
GGGGGTGTAGCTCAGTGGTAGAGCGCGTGCTTAGCATGCACGAGGCCCCGGGTTCAATCC
CTGGCACCTCCA
>HG549 HG983549.1
GGGGGTGTAGCTCAGTGGTAGAGCGCGTGCTTAGCATGCACGAGGCCCCGGGTTCAATCC
CCGGCACCTCCA
>HG587 HG983587.1
ATTGCCGAATGGACAGTCAGTCCAGAATGCCAGTCAGTCGGTCCCAAACCCTTGTCACTG
CCGGCATCTCCA
>HG193 HG984193.1
GGGGATGTAGCTCAGTGGTAGAGCGCGTATGCTTCGCATGTATGAGGCCCCGGGTTCGATCC
CCGGCATCTCCA
>KH526 KH780526.1
GGGGATGTAGCTCAGTGGTAGCGCATGCTTAGCATGCATGAGGTCCCGGGTTCGATCC
CCAGCATCTCCA
```

And this is input database for local search
```
# this file is for local search test by SSEARCH or your own program
cat minilocal.fasta

>CR570 CR163570.1
    CTTCGAAAGCTGTTGGTTCGGACCCGGGAGCGCGCAGATGGTCAGATGCCATGGATTTGG
    ACTTGATTTTCAAAAAAATCGGTGAAAATAATTATGGCCTTCGGCTTCTGTTCATACTTC
    CAAAACCTAGTAAAAGTCAGATTTCTGTCTTGAATGTATTCCCAAGAGCTTAGCTGTGTC
    GGGGGTGTAGCTCAGTGGTAGAGCGCGTGCTTAGCATGCACGAGGCCCCGGGTTCAATCC
    CCGGCACCTCCAGTGTTTTCAAGCCTCCATCTGCTCGCTTTTCTCCAACAAATATACTTA
    CACTTTTCCTTGCTCGAAGAGAAAGTCAAACCCTTCCATACAATATTATGCTGCGGGACA
    AACAGCAAACCTTCACTTCACTTCACTGCTGCCGAGCCGATACCATCTCATAGCATCAGA
    AGAAGATGGAAAATCAAAATACAGCTAGGAACAGCGGTTCTCTCAGTCCAAAAACGGCAA
    GGTTCAGGTTCTGAAGGTAACGTGACCATCTACCTTCATTTCCAGGAAAGAAAAGGATCA
    ATGGTAACAGGAGAATGGCACAGGTTCTCTAAATTTTCCTGTTGAAATAGCTGTTTTTAC
    CGCATCGCCCTACAGTTCGGCACGTGGGTTCCCGGTTTGTTTTCAGATGATCCCTGGAAT
    CCACTCACCTGAGACAGCAGCGCCGTCTGTCGGAGAGACGGTGTCGAGACCCTAGAGAGA
    AATCGCAATCTGGGTTAGCAATGTACGAAAATTAGAACCGGATTTTAAATGGTAATCCAT
    GCCGATCCTTCG
>AG587 AG736587.1
    GTATGCTCAGAGTCGTCAATTTGGCGACGCCATTGAAGTTTTGGAATCGGGAGCCGTCTT
    TGTCTTCCAGCTCCATCTTTTCCACCTTTTGCTTGGGCAGCCCCCGTGAGTCGTGTCAAG
    GCTGACGAGTAGAAATGGAACAGCACTAATATTAATGGCGAAACCTTTGTAAAATACGTT
    TGGTTTCTGTTTAACAAGGAAAAGAAAGTAAAGCAATGGAAAAGAATTTAAACACAAAAA
    CAACTGGAGGTGCCGGGGATTGAACCCGGGGCCTCGTGCATGCTAAGCACGCGCTCTACC
    ACTGAGCTACACCCCCTTATCGAAACTACTCTCTCGAGAATATATTCAAGACAGAATCTG
    ACCATTTTGCTACGTTTCAGAACAATTAGTTGTAACAAGCCAAAGTCTACTTTATTTAGC
    TATTTCTTATATCTTAAATTTAGGTTTTGTGTCCCTGTTTGCTGATAGCTGAGCAAACCG
    CATTCTACACCGAAGGCCTCTATTGATGGTCCTGGGATTTTTCTGCTCGTCAGTCCGGAG
    TCACTTCCCGGCCACCACTAGAAGAACCAGGGATGAAATTTTCTCCTGAGTTTTGACTGT
    CTCCTTTCTTTCTCCGCTGCTTTAAAGGGCTGCGAGAAAGTCACAAAGTACAAAGCGAAG
    CATTTAGAGACCATATTAGATGCAGATGGCGAGGGAAGACAGGTGGACAANAGCAGACAG
    GTTCGTGTCGGTGCAGCCACTGCTTTGGACCCGAGCGTCCGTCC
>HE540 HE318540.1
    TCTGGTGTCCTATTCTTCTTTCACCCAGGAGGTTGTAATCCTTGGCGGTTTACTCTTCAA
    GTCAATGCTCAGTGTTTGCAAGAGTGGACTCGAAGACACTGTGCTCTTTGGGGCACGCCA
    TCCTCTTTTGAGCACTCGAAAGACTGTCTAGGGCATCTTGTACTTTTCCTGCCCTAGTGC
    TTCTAGTGACCAGTTTTGCCCGGAGCTTTGGTTCCATCCGGGGGGATGTGTATTAAGCTG
    CTGAGATGGGGGACCAGCTCTCCTCCTTGCCCCTTGGGTGTCTTTTTGAGTCTAGGCTTT
    TTTTAACAGGCAGAACTAGGAAAGGGCTCTTCTTTGCCTTCCTCCATTCAAGTTGTCTGT
    AGGCCTGACACTGGGGTAGGCTGACAAGCGGGACTAGATCCTAGGAAGGGCGCCTAAATT
    TTAGGCCGGCATTCGGTTGCCTTTCTCCTCATGTCCTCGGCGAAAAGGGAATAGGAGTTT
    CGCAAGTAACCGGGTTGTACACAACGGAGCTTCGAACCAACACGTGAGGAGTTGCCATAC
    CATTTTTTGTTTGTGTGTTTGGTTGTTTGTTTGCTTTCCCTTTTCCTGCCTCTGAGGTTT
    CGGCAGCCGCGGGGGTGTAGCTCAGTGGTAGAGCGCGTGCTTAGCATGCACGAGGCCCCG
    GGTTCAATCCCCGGCACCTCCAACCGCCACGTCCACTTGCTCCGCCCTTCCTTGCTCTC
>HE391 HE161391.1
    AGAATCTCCTGTGCACGAACACCGCCAGGCTGAACTCTCACCACCGTGTCCTCAACGGAC
    AGGAGAGCAACTGCCGCAAACATGCGAAAATTTGAAAAGGCGTCACTCGTCACAAACGAC
    GTTCCTCGCGCGAAACATCCTTTAAACGTTCTCCGAATCCAGGAACAACGCGGTGGCGGT
    GGCTGCCGACGCCGGTTTGGTGAAAGATGAGACGAGGGCGTGGTGAGAGTGATTTCCTGC
    TCTCAGAAGTCCTGGTAGTTGATTTTTGTGGCGTATCGATCGGCCTTGGCCAGGAGAGGA
    GTGATTGCGCAGAGTGACTAGAATGGATGGGTTTAAAGAAAAAGAAAAGAAAAGAAAAGA
    AAAGAAAAGAAGAAAAGAAAAAAAAAGAAGAAAAGAAAGAAAGACAAAGTTGGACGCATA
    CTACTTTCGCGGCAGTCAGCAGGGACTGTGAGCTCCCGCCTCGCACTGTCTCTCCAGAGT
    TGAGAGCAAGGAAGGGCGGAGCAAGTGGACGTGGCGGTTGGAGGTGCCGGGGATTGAACC
    CGGGGCCTCGTGCATGCTAAGCACGCGCTCTACCACTGAGCTACACCCCCGCGGCTGCCG
    AAACCTCAGAGGCAGGAAAAGGGAAAGCAAACAACCAACCCAACACACAAACAAAAAATG
    GTAT
```

## Alignment score and gap panelty
Alignment score is to be calculated as following:
`AlignmentScore = 5*NumberOfMatches -4*NumberOfMismatches - (15 + 4*(GapLenght - 1))`
We use affine gap penalty method as described sucintly here: http://rosalind.info/problems/gaff/ You can read more about [Gap penalty](https://en.wikipedia.org/wiki/Gap_penalty), and more about [affine](http://homepage.usask.ca/~ctl271/857/affine_gap_penalties.shtml) gap panelty.

## Expected output

We gonna mimic GGSEARCH for global search and SSEARCH for local search. It may be quite complicate to reproduce output exactly as the output of the [FASTA commands](#ggsearch). So make the output and output constraints as close as possible to these tools.
Check this again about their scoring and output constraints:

```
# In essence:
# -r, -f, and -g are parameters for calculating and selecting alignment with highest score when considering two sequences.
# -E, -F, -b, and -d form constraints for how many hits to give back from database for a query sequence.

```
It maybe complicated to impletment constraints based on -E, -F, -b and -d parameters. This is from the [author](https://github.com/wrpearson/fasta36/issues/8) of `fasta36` repo:

>For ggsearch36, the alignment score is calculated by the do_work() function in the dropnnw2.c file, for ssearch36, the score is calculated by do_work() in the dropgsw2.c file. 
>
>The E()-value calculation for local scores is done by functions in scaleswn.c. ssearch and other local comparison methods use the extreme value distribution. ggsearch uses a normal distribution to calculate the significance of global alignment scores.

You these information as much as you can, or device a method yourself to limit the number of hits for a query sequence. If you do so, please write documentation to explain it.
