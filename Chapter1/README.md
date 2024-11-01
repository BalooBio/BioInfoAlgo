# Finding-Ori
This repository contains my solutions to the code challenges found in Chapter 1 of Bioinformatics Algorithms: An Active Learnign Approach, 3rd Edition, by Phillip Compeau and Paval Pevzner, as well as my solution to the final challenge of finding the origin of replication (OriC) and DnaA Boxes in the Salmonella enterica genome. The code challenges (1A-1N) can be viewed as individual python files, or in a single file called chapter1.py, in the code_challenges directory. The final code for finding the OriC and DnaA Boxes in the Salmonella enterica genome is located in the DnaA_box_finder directory and consists of two python files; finding_ori.py and supporting_functions.py.

finding_ori.py will prompt the user for a file name containing the genome to anlyze. Bacterial genomes for S. enterica and Escherichia Coli are stored as fasta files in the genomes directory. The program will return the locations and sequences of the OriC and DnaA boxes for the provided genome. NOTE: This program will only work for bacterial genomes closely related to E. coli having DnaA boxes similar to the consensus binding sequence 5' TTATC[CA]A[CA]A 3'.

A breif desctiption of the finding ori solution is provided below. All information herein has been adapted from, or inspired by, Chapter 1 of Bioinformatics Algorithms: An Active Learnign Approach, which I recommend reading for a more complete understanding of the problem and this solution set.

## STEP 1: Identifying the Approximate Position of OriC

The approximate position of OriC is identified using the function minimum_skew(seq), which accepts a DNA sequence, seq, as an arument (in this case, a bacterial genome) and returns the index position(s) for which GC-skew reaches it mimimum. What is GC-skew? Due to the nature of DNA replication, DNA on the 3' side of OriC is subject to a higher rate of deanimation mutions, which changes Cytosine to Thymine (C -> T), compared to DNA on the 5' side of the OriC. As a result, DNA on the 5' side of OriC will have more C's than G's and DNA on the 3' side of OriC will have more G's than C's. The funciton, minimum_skew(seq) traverses the genome 1 nucleotide at a time, updating a tracker, GC-skew, by +1 when it encounters a G and -1 when it encounters a C. The genome position(s) with the lowest (most negative) GC-skew value represent the loaction(s) where the DNA transitions from having more C's than G's to having more G's than C's, thus providing an approximate estimation of OriC's locaiton. Typically, more than one index position is returned by the minimum_skew() function, but for the bacterial geneomes included in the genomes directory, these positions are within a few base-pairs(bp) of one another and thus represent a single OriC, which average 500bp in length. Two minima are returned for the Salmonella enterica genome: [4084376, 4084378]. finding_ori.py will use the first minima to pinpoint the location of OriC.

## STEP 2: Pinpointing the Location of OriC Using the DNA Box Consensus Sequence of E. coli
If OriC were centered on the GC minima, we could simply take the 500 bp sequence centered on the GC-minima (genome[GC_minima - 250: GC_minima + 249]) and search for DnaA boxes. Unfortunately, GC-minima is only the approximate position of OriC, so we should expect that OriC could be centered up to 500 bp on either side of the GC minima. Instead of searching one 500bp sequence for DnaA Boxes, finding_ori.py itterates through 1000, 500bp sequeneces, starting with the sequence centered on GC_minima - 500 and ending with the sequence centered on GC_minima + 500.

An earlier version of finding_ori.py passed each 500bp to a function called final_frequent_words(text, k, d), which returns the set of k-mers with the most frequent occurance in text (a DNA sequence) allowing for up to d number of mismatches and considering reverse complements. If any of the most common 9-mers in a given sequence matched the DnaA box consensus sequence  5' TTATC[CA]A[CA]A 3, the index position at the center of that sequence was added to a list of potential oriC index positons. This approach takes some time to compute, and returns the following 77 consecutive index positions for the S. enterica genome:

[4084137, 4084138, 4084161, 4084162, 4084163, 4084164, 4084165, 4084166, 4084167, 4084168, 4084169, 4084170, 4084171, 4084172, 4084173, 4084174, 4084175, 4084176, 4084177, 4084178, 4084179, 4084180, 4084181, 4084182, 4084183, 4084184, 4084185, 4084186, 4084187, 4084188, 4084189, 4084190, 4084191, 4084192, 4084193, 4084194, 4084195, 4084196, 4084197, 4084198, 4084199, 4084200, 4084201, 4084202, 4084203, 4084204, 4084205, 4084206, 4084207, 4084208, 4084209, 4084210, 4084211, 4084212, 4084213, 4084214, 4084215, 4084216, 4084217, 4084218, 4084219, 4084220, 4084221, 4084222, 4084223, 4084224, 4084225, 4084226, 4084227, 4084228, 4084229, 4084230, 4084231, 4084232, 4084233, 4084234, 4084235]

From here, the most reasonable approach would be to select the index at the center of this list as the proposed OriC/ To save time, I've modified finding_ori.py to pass every 100th 500bp sequence to final_frequent_words(text, k, d), meaning that this is done with 10, 500 bp sequences rather than 1000. This returns a single index, 4084176.

## STEP 3: Identify Specific DnaA Boxes in the Proposed OriC and Print Results

finding_ori_py prints the starting and ending position of Oric, identified in step 2, and prints the full 500 bp OriC seuqence. It then iterates through each 9-mer in the 500bp sequence to identify DnaA Boxes that are represented by the consensus sequence  5' TTATC[CA]A[CA]A 3'. The starting position and the 9bp sequence of each DnaA box is printed.

Running fining_ori.py with the S. enterica genome prints the following:

    Proposed Origin:

    Position = 4083926:4084425

    Sequence = GATCTGTTCTATTGTGATCTCTTATTAGGATCGCGCCAGGCTGTGGATAACCCGGATCCTGTAATAAAGATCAATGCGTTGGAAAGGATCACTAGCTGTGAATGATCGGTGATCGTGGTCCGTATAAGCTGGGATCAAAACGGGTACTTATACACAACTCAAAAAGTGAACAACGGTTATTCTTTGGATAACTACCGGTTGATCCAAGCTTTCCACCAGATTTATCCACAATGGATCGCACGATCTTTACACTTATTTGAGTAAATTAATCCAGGATCCGAGCCAAATCTCCGCTGGATCTTCCGGAATCTCATGTTCAAGGATGTTGATCTTCAGTGTTTCCCCAACCTGTTTTGCGCCAGCGCCTTTCAGTTCCGCTTCTATTTTCTCAATCGCGCCGCAAAACGTGTCGTATTCTCGACTGCCAATGCCAATTGCGCCGAAACGTACCGCGGAAAGATCGGGTTTCTGCGTCTGAAGGTCTTCATAGAAAGGGGTC

    Proposed DnaA Boxes:

    1: Position = 4083967, Sequence = TGTGGATAA

    2: Position = 4084073, Sequence = TTATACACA

    3: Position = 4084108, Sequence = TTTGGATAA

    4: Position = 4084147, Sequence = TTATCCACA


Running fining_ori.py with the E. coli genome prints the following:

    Proposed Origin:

    Position = 3923569:3924068

    Sequence = ACCCATACGCGCCGCGGCCATCGCGGCCTCGGTGCCTGCATGACCCCCGCCAATGATGATGACGTCAAAAGGATCCGGATAAAACATGGTGATTGCCTCGCATAACGCGGTATGAAAATGGATTGAAGCCCGGGCCGTGGATTCTACTCAACTTTGTCGGCTTGAGAAAGACCTGGGATCCTGGGTATTAAAAAGAAGATCTATTTATTTAGAGATCTGTTCTATTGTGATCTCTTATTAGGATCGCACTGCCCTGTGGATAACAAGGATCCGGCTTTTAAGATCAACAACCTGGAAAGGATCATTAACTGTGAATGATCGGTGATCCTGGACCGTATAAGCTGGGATCAGAATGAGGGGTTATACACAACTCAAAAACTGAACAACAGTTGTTCTTTGGATAACTACCGGTTGATCCAAGCTTCCTGACAGAGTTATCCACAGTAGATCGCACGATCTGTATACTTATTTGAGTAAATTAACCCACGATCCCAGCCAT

    Proposed DnaA Boxes:

    1: Position = 3923823, Sequence = TGTGGATAA

    2: Position = 3923929, Sequence = TTATACACA

    3: Position = 3923964, Sequence = TTTGGATAA

    4: Position = 3924003, Sequence = TTATCCACA
