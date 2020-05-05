# MultiGeneBlastParser

### A script to simplify the output of MultiGeneBlast


This script will convert the output of MultiGeneBlast (Medema, M. *et al* 2013) in to distinct entries for each hit, for easier downstream processing,

It depends on the `pandas` module as the only non-standard module.
This can be installed as follows:

     $ python -m pip install pandas

or

     $ conda install pandas

## Full Script Options

```
usage: MGBparser.py [-h] [-v] [-r] [-q] [-s] [-b] [-o OUTFILE] [-m MAX_RESULT]
                    clusterfile

Parse and extract data from the output of MultiGeneBlast.

positional arguments:
  clusterfile           The text file of hits output by MGB. By default this
                        is called 'clusterblast_output.txt'.

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Verbose behaviour, printing parameters of the script.
  -r, --references      Display relevant references. This option cannot be
                        used with any others.
  -q, --query           Output the query sequence information (and write a
                        file). [Boolean]
  -s, --sighits         Output the table of significant hits (and write a
                        file). [Boolean]
  -b, --blastfile       Output the table of BLAST hits for each match (and
                        write a file). [Boolean]
  -o OUTFILE, --outfile OUTFILE
                        Name stem for all output files (extensions are set
                        internally).
  -m MAX_RESULT, --max_result MAX_RESULT
                        Return results for just the top n details sections
                        [Def: 50].
```

### N.B. Python2 _only_ at the present time. Pull requests to port to Py3 are welcome.

###### A simple invocation may look like:

    python MGBparser.py --outfile myclusters clusterblast_output.txt

This would simply return as many hits as are in the input file, but will silently write files with the BLAST infomation and operon locations with the filename prefix `myoperon*`.
Note, the input file is a _positional_ argument, and must come last. Note also that the output file is a _name stem_, and not the full filename. Output files will be created where the input file originated.


A more comprehensive invocation might look like:

    python MGBparser.py  --outfile myclusers -m 10 -qsbv clusterblast_output.txt

Which would run the script as before, but this time only return the top 10 hits (`-m 10`), whilst
also printing and writing files with the details of the input sequence, BLAST results, operon locations,
significant hit descriptions and would print many more status and information messages.


## Output

Here is some example output (as generated from the `/test_data/single_clusterblast_output.txt` file):

(Note there is no need for the `-m` flag as there is only a single entry.)

```$ python MGBparser.py --clusterfile test_data/single_clusterblast_output.txt --out single -vqbs
Opening test_data/single_clusterblast_output.txt for reading...
Parsing Query details section...
Parsing Significant hit section...
Parsing Hit details section...
Splitting Hit details section...
Establishing Hit classes...
Writing query details to file: test_data/single_queryinfo.tsv
Query sequence information:
===========================
Input file:PVCcif_ATCC43949.gbk
        Locus  Start   Stop Strand                               Annotation       Comment
0   PAU_01961      8    457      +     T4-like_virus_tail_tube_protein_gp19  no_locus_tag
1   PAU_01962    521   1618      +                major_tail_sheath_protein  no_locus_tag
2   PAU_01963   1799   3280      +                      tail_sheath_protein  no_locus_tag
3   PAU_01964   3334   4533      +                      tail_sheath_protein  no_locus_tag
4   PAU_01965   4547   5005      +     T4-like_virus_tail_tube_protein_gp19  no_locus_tag
5   PAU_01966   5002   5181      +                     hypothetical_protein  no_locus_tag
6   PAU_01967   5168   5851      +                     hypothetical_protein  no_locus_tag
7   PAU_01968   5848   7449      +                  Rhs_element_Vgr_protein  no_locus_tag
8   PAU_01969   7462   7905      +                  baseplate_wedge_subunit  no_locus_tag
9   PAU_01970   7902   8318      +                     hypothetical_protein  no_locus_tag
10  PAU_01971   8527  11253      +                     hypothetical_protein  no_locus_tag
11  PAU_01972  11246  14197      +                     hypothetical_protein  no_locus_tag
12  PAU_01973  14333  15184      +                     hypothetical_protein  no_locus_tag
13  PAU_01974  15247  17145      +                     hypothetical_protein  no_locus_tag
14  PAU_01975  17155  19227      +  ATP-dependent_zinc_metalloprotease_FtsH  no_locus_tag
15  PAU_01976  19252  20166      +                     hypothetical_protein  no_locus_tag
16  PAU_01977  20327  21223      +                     hypothetical_protein  no_locus_tag
17  PAU_01978  21308  22291      +                     hypothetical_protein  no_locus_tag
18  PAU_01979  22788  23684      +                     hypothetical_protein  no_locus_tag
19  PAU_01980  23656  24114      +                     hypothetical_protein  no_locus_tag
Writing Significant hit details to file: test_data/single_sighitinfo.tsv
Significant Hit information:
============================
   Rank      ID                                Description
0     1   PAU_1  Photorhabdus asymbiotica strain ATCC43949
Writing Hit details for: 1. PAU_1 to test_data/single_PAU_1_locationinfo.tsv (1 of 1)
Hit location information:
============================
        Locus    Start     Stop Strand                            Annotation       Comment
0   PAU_01961  2233799  2234248      +  T4-like_virus_tail_tube_protein_gp19  no_locus_tag
1   PAU_01962  2234312  2235409      +             major_tail_sheath_protein  no_locus_tag
2   PAU_01963  2235590  2237071      +                   tail_sheath_protein  no_locus_tag
3   PAU_01964  2237125  2238324      +                   tail_sheath_protein  no_locus_tag
4   PAU_01965  2238338  2238796      +  T4-like_virus_tail_tube_protein_gp19  no_locus_tag
5   PAU_01966  2238793  2238972      +                  hypothetical_protein  no_locus_tag
6   PAU_01967  2238959  2239642      +                  hypothetical_protein  no_locus_tag
7   PAU_01968  2239639  2241240      +               Rhs_element_Vgr_protein  no_locus_tag
8   PAU_01969  2241253  2241696      +               baseplate_wedge_subunit  no_locus_tag
9   PAU_01970  2241693  2242109      +                  hypothetical_protein  no_locus_tag
10  PAU_01971  2242318  2245044      +                  hypothetical_protein  no_locus_tag
11  PAU_01972  2245037  2247988      +                  hypothetical_protein  no_locus_tag
12  PAU_01973  2248124  2248975      +                  hypothetical_protein  no_locus_tag
13  PAU_01974  2249038  2250936      +                  hypothetical_protein  no_locus_tag
14  PAU_01976  2253043  2253957      +                  hypothetical_protein  no_locus_tag
15  PAU_01977  2254118  2255014      +                  hypothetical_protein  no_locus_tag
16  PAU_01978  2255099  2256082      +                  hypothetical_protein  no_locus_tag
17  PAU_01979  2256579  2257475      +                  hypothetical_protein  no_locus_tag
18  PAU_01980  2257447  2257905      +                  hypothetical_protein  no_locus_tag
Writing Hit BLAST details for: 1. PAU_1 to test_data/single_PAU_1_blastinfo.tsv (1 of 50)
Hit location information:
============================
        Locus    Start     Stop Strand                            Annotation       Comment
0   PAU_01961  2233799  2234248      +  T4-like_virus_tail_tube_protein_gp19  no_locus_tag
1   PAU_01962  2234312  2235409      +             major_tail_sheath_protein  no_locus_tag
2   PAU_01963  2235590  2237071      +                   tail_sheath_protein  no_locus_tag
3   PAU_01964  2237125  2238324      +                   tail_sheath_protein  no_locus_tag
4   PAU_01965  2238338  2238796      +  T4-like_virus_tail_tube_protein_gp19  no_locus_tag
5   PAU_01966  2238793  2238972      +                  hypothetical_protein  no_locus_tag
6   PAU_01967  2238959  2239642      +                  hypothetical_protein  no_locus_tag
7   PAU_01968  2239639  2241240      +               Rhs_element_Vgr_protein  no_locus_tag
8   PAU_01969  2241253  2241696      +               baseplate_wedge_subunit  no_locus_tag
9   PAU_01970  2241693  2242109      +                  hypothetical_protein  no_locus_tag
10  PAU_01971  2242318  2245044      +                  hypothetical_protein  no_locus_tag
11  PAU_01972  2245037  2247988      +                  hypothetical_protein  no_locus_tag
12  PAU_01973  2248124  2248975      +                  hypothetical_protein  no_locus_tag
13  PAU_01974  2249038  2250936      +                  hypothetical_protein  no_locus_tag
14  PAU_01976  2253043  2253957      +                  hypothetical_protein  no_locus_tag
15  PAU_01977  2254118  2255014      +                  hypothetical_protein  no_locus_tag
16  PAU_01978  2255099  2256082      +                  hypothetical_protein  no_locus_tag
17  PAU_01979  2256579  2257475      +                  hypothetical_protein  no_locus_tag
18  PAU_01980  2257447  2257905      +                  hypothetical_protein  no_locus_tag
Writing Hit coordinates for: 1. PAU_1 to test_data/single_PAU_1_coords.tsv (1 of 1)
Hit coordinate information:
===========================
Hit No  ID      Start Locus     End Locus       Start Index     End Index       Main Strand     Source
1       PAU_1   PAU_01961       PAU_01980       2233799 2257905 +       Photorhabdus asymbiotica strain ATCC43949
```



MultiGeneBlast has been previously published here:

> Medema MH, Takano E, Breitling R.
> "Detecting Sequence Homology at the Gene Cluster Level with MultiGeneBlast."
> Molecular Biology and Evolution. 2013;30(5):1218-1223. doi:10.1093/molbev/mst025.
