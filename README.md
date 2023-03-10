# motif-mark
This is a Python script based on object-oriented programming that generates a .png image indicating binding sites for up to 5 motifs (each motif being 10 bases or less) on up to 10 sequences (each sequence being 1000 bases or less). The resulting image portrays introns and exons and also shows motif, sequence, intron, and exon sizes to scale. The prefix of the created .png file matches the prefix used for the FASTA file.

The script requires a FASTA file that has a maximum of 10 sequences, where each sequence should not exceed 1000 bases. Introns are to be represented using lowercase nucleotides, while uppercase nucleotides are to be used for exons. Additionally, a txt file is needed containing one motif per line. Motifs can include any IUPAC degenerate base symbol, and they can be denoted in uppercase or lowercase.

Usage: ./motif-mark-oop.py -f [Fasta File] -m [Motif File]


Example Output:
![alt text](https://github.com/pranavmuthu/motif-mark/blob/main/Motifs.png?raw=true)
