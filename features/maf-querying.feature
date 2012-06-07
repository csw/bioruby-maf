@milestone_3
Feature: Filter results from MAF files
  In order to work with only relevant data from a MAF file
  Such as only species recognized by PhyloCSF
  I want to filter the results of MAF queries

  Scenario: Return only specified species
    Given MAF data:
    """
    ##maf version=1
    a score=10542.0
    s mm8.chr7                 80082334 34 + 145134094 GGGCTGAGGGC--AGGGATGG---AGGGCGGTCC--------------CAGCA- 
    s rn4.chr1                136011785 34 + 267910886 GGGCTGAGGGC--AGGGACGG---AGGGCGGTCC--------------CAGCA- 
    s oryCun1.scaffold_199771     14021 43 -     75077 -----ATGGGC--AAGCGTGG---AGGGGAACCTCTCCTCCCCTCCGACAAAG- 
    s hg18.chr15               88557580 27 + 100338915 --------GGC--AAGTGTGGA--AGGGAAGCCC--------------CAGAA- 
    s panTro2.chr15            87959837 27 + 100063422 --------GGC--AAGTGTGGA--AGGGAAGCCC--------------CAGAA- 
    s rheMac2.chr7             69864714 28 + 169801366 -------GGGC--AAGTATGGA--AGGGAAGCCC--------------CAGAA- 
    s canFam2.chr3             56030570 39 +  94715083 AGGTTTAGGGCAGAGGGATGAAGGAGGAGAATCC--------------CTATG- 
    s dasNov1.scaffold_106893      7435 34 +      9831 GGAACGAGGGC--ATGTGTGG---AGGGGGCTGC--------------CCACA- 
    s loxAfr1.scaffold_8298       30264 38 +     78952 ATGATGAGGGG--AAGCGTGGAGGAGGGGAACCC--------------CTAGGA 
    s echTel1.scaffold_304651       594 37 -     10007 -TGCTATGGCT--TTGTGTCTAGGAGGGGAATCC--------------CCAGGA 
    """
    When I open it with a MAF reader
    And filter for only the species
    | hg18    |
    | mm8     |
    | rheMac2 |
    Then an alignment block can be obtained
    And the alignment block has 3 sequences

  Scenario: Return only blocks having all specified species
    Given a MAF source file "mm8_chr7_tiny.maf"
    When I open it with a MAF reader
    And build an index on the reference sequence
    And filter for blocks with the species
    | panTro2 |
    | loxAfr1 |
    And search for blocks between positions 80082471 and 80082730 of mm8.chr7
    Then 1 block is obtained

  Scenario: Return only blocks having a certain number of sequences
    Given a MAF source file "mm8_chr7_tiny.maf"
    When I open it with a MAF reader
    And build an index on the reference sequence
    And filter for blocks with at least 6 sequences
    And search for blocks between positions 80082767 and 80083008 of mm8.chr7
    Then 1 block is obtained

