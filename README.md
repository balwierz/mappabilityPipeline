
Output: a bigWig file with signal value equal to *the shortest read length necessary to unique map a read to a given locus*.
[currently a pseudo-bedGraph file, you would need to sort, conver to bedGraph, and then to bigWig manually]

Uses bowtie2. Iteratively maps a chopped genome to itself to test if the mapping is unique.
Filters for mapping quality q>=10 (can be easily changed).
Keeps a hash table of genomic coordinates in the memory. For the mouse genome at least 90GB or RAM is necessary.

Note: Bowtie2 sometimes reports a very good alignment (q=42) of a sequence to a wrong locus.
Probably somehow does not discover the second (intended) locus and reports only one perfect match somewhere else.
This is really wrong, but it is bowtie2 where this should be fixed.
The vast majority of positions and mappings are OK.

Piotr Balwierz
