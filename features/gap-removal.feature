Feature: Remove gaps from MAF files
  In order to work with only the alignment data involving sequences
  Which can be used by downstream software
  We may want to filter out certain species
  Which can leave gap regions where sequence data was only present
  For removed species
  So it is useful to be able to remove those gaps

  Background:
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
    | mm8     |
    | rn4     |
    | hg18    |
    | canFam2 |
    | loxAfr1 |

  Scenario: Detect filtered blocks
    When an alignment block can be obtained
    Then the alignment block is marked as filtered
    And the alignment block has 5 sequences

  Scenario: Detect gaps
    When an alignment block can be obtained
    Then 1 gap is found with length [14]

  Scenario: Remove gaps
    When an alignment block can be obtained
    And gaps are removed
    Then the text size of the block is 40

  Scenario: Remove gaps in the parser
    When I enable the :remove_gaps parser option
    And an alignment block can be obtained
    Then the text size of the block is 40
