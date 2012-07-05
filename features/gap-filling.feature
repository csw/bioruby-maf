Feature: Join alignment blocks with reference data
  In order to produce FASTA output with one sequence per species
  For use in downstream tools
  We need to join adjacent MAF blocks together
  And fill gaps in the reference sequence from reference data

  Scenario: Non-overlapping MAF blocks in region of interest
    Given MAF data:
    """
    ##maf version=1
    a score=20.0
    s sp1.chr1        10 13 +      50 GGGCTGAGGGC--AG
    s sp2.chr5     53010 13 +   65536 GGGCTGACGGC--AG
    s sp3.chr2     33010 15 +   65536 AGGTTTAGGGCAGAG

    a score=21.0
    s sp1.chr1        30 10 +      50 AGGGCGGTCC
    s sp2.chr5     53030 10 +   65536 AGGGCGGTGC
    """
    And chromosome reference sequence:
    """
    > sp1.chr1
    CCAGGATGCT
    GGGCTGAGGG
    CAGTTGTGTC
    AGGGCGGTCC
    GGTGCAGGCA
    """
    When I open it with a MAF reader
    And build an index on the reference sequence
    And tile sp1.chr1:0-50 with the chromosome reference
    And tile with species [sp1, sp2, sp3]
    And write the tiled data as FASTA
    Then the FASTA data obtained should be:
    """
    > sp1
    CCAGGATGCTGGGCTGAGGGC--AGTTGTGTCAGGGCGGTCCGGTGCAGGCA
    > sp2
    **********GGGCTGACGGC--AG*******AGGGCGGTGC**********
    > sp3
    **********AGGTTTAGGGCAGAG***************************
    """

  Scenario: Non-overlapping MAF blocks with species map
    Given MAF data:
    """
    ##maf version=1
    a score=20.0
    s sp1.chr1        10 13 +      50 GGGCTGAGGGC--AG
    s sp2.chr5     53010 13 +   65536 GGGCTGACGGC--AG
    s sp3.chr2     33010 15 +   65536 AGGTTTAGGGCAGAG

    a score=21.0
    s sp1.chr1        30 10 +      50 AGGGCGGTCC
    s sp2.chr5     53030 10 +   65536 AGGGCGGTGC
    """
    And chromosome reference sequence:
    """
    > sp1.chr1
    CCAGGATGCT
    GGGCTGAGGG
    CAGTTGTGTC
    AGGGCGGTCC
    GGTGCAGGCA
    """
    When I open it with a MAF reader
    And build an index on the reference sequence
    And tile sp1.chr1:0-50 with the chromosome reference
    And tile with species [sp1, sp2, sp3]
    And map species sp1 as mouse
    And map species sp2 as hippo
    And map species sp3 as squid
    And write the tiled data as FASTA
    Then the FASTA data obtained should be:
    """
    > mouse
    CCAGGATGCTGGGCTGAGGGC--AGTTGTGTCAGGGCGGTCCGGTGCAGGCA
    > hippo
    **********GGGCTGACGGC--AG*******AGGGCGGTGC**********
    > squid
    **********AGGTTTAGGGCAGAG***************************
    """

  Scenario: Subset of non-overlapping MAF blocks in region
    Given MAF data:
    """
    ##maf version=1
    a score=20.0
    s sp1.chr1        10 13 +      50 GGGCTGAGGGC--AG
    s sp2.chr5     53010 13 +   65536 GGGCTGACGGC--AG
    s sp3.chr2     33010 15 +   65536 AGGTTTAGGGCAGAG

    a score=21.0
    s sp1.chr1        30 10 +      50 AGGGCGGTCC
    s sp2.chr5     53030 10 +   65536 AGGGCGGTGC
    """
    And chromosome reference sequence:
    """
    > sp1.chr1
    CCAGGATGCT
    GGGCTGAGGG
    CAGTTGTGTC
    AGGGCGGTCC
    GGTGCAGGCA
    """
    When I open it with a MAF reader
    And build an index on the reference sequence
    And tile sp1.chr1:12-36 with the chromosome reference
    And tile with species [sp1, sp2, sp3]
    And write the tiled data as FASTA
    Then the FASTA data obtained should be:
    """
    > sp1
    GCTGAGGGC--AGTTGTGTCAGGGCG
    > sp2
    GCTGACGGC--AG*******AGGGCG
    > sp3
    GTTTAGGGCAGAG*************
    """
  Scenario: Overlapping MAF blocks in region of interest
    Given MAF data:
    """
    ##maf version=1
    a score=20.0
    s sp1.chr1        10 13 +      50 GGGCTGAGGGC--AG
    s sp2.chr5     53010 13 +   65536 GGGCTGACGGC--AG
    s sp3.chr2     33010 15 +   65536 AGGTTTAGGGCAGAG

    a score=21.0
    s sp1.chr1        20 10 +      50 AGGGCGGTCC
    s sp2.chr5     53020 10 +   65536 AGGGCGGTGC
    """
    And chromosome reference sequence:
    """
    > sp1.chr1
    CCAGGATGCT
    GGGCTGAGGG
    CAGTTGTGTC
    AGGGCGGTCC
    GGTGCAGGCA
    """
    When I open it with a MAF reader
    And build an index on the reference sequence
    And tile sp1.chr1:0-50 with the chromosome reference
    And tile with species [sp1, sp2, sp3]
    And write the tiled data as FASTA
    Then the FASTA data obtained should be:
    """
    > sp1
    CCAGGATGCTGGGCTGAGGGAGGGCGGTCCAGGGCGGTCCGGTGCAGGCA
    > sp2
    **********GGGCTGACGGAGGGCGGTGC********************
    > sp3
    **********AGGTTTAGGG******************************
    """



