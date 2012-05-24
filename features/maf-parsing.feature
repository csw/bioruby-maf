Feature: Parse MAF files
  In order to extract information from a MAF file
  I want to read it and pull out information

  Scenario: Read MAF header
    Given MAF data:
    """
    ##maf version=1 scoring=humor.v4
    # humor.v4 R=30 M=10 /cluster/data/hg15/bed/blastz.mm3/axtNet25/chr22.maf /cluster/data/hg15/bed/blastz.rn3/axtNet25/chr22.maf

    a score=0.128
    s human_hoxa 100  8 + 100257 ACA-TTACT
    s horse_hoxa 120  9 -  98892 ACAATTGCT
    s fugu_hoxa   88  7  + 90788 ACA--TGCT
    """
    And in a temp file
    When I open it with a MAF reader
    Then the MAF version should be "1"
    And the scoring scheme should be "humor.v4"
    # third line a continuation
    And the alignment parameters should be "humor.v4 R=30 M=10 /cluster/data/hg15/bed/blastz.mm3/axtNet25/chr22.maf /cluster/data/hg15/bed/blastz.rn3/axtNet25/chr22.maf" 

  Scenario: Read alignment block
    Given MAF data:
    """
    ##maf version=1 scoring=humor.v4
    # humor.v4 R=30 M=10 /cluster/data/hg15/bed/blastz.mm3/axtNet300/chr1.maf
    # /cluster/data/hg15/bed/blastz.rn3/axtNet300/chr1.maf

    a score=0.128
    s human_hoxa 100  8 + 100257 ACA-TTACT
    s horse_hoxa 120  9 -  98892 ACAATTGCT
    s fugu_hoxa   88  7  + 90788 ACA--TGCT
    """
    And in a temp file
    When I open it with a MAF reader
    Then an alignment block can be obtained
    And the alignment block has 3 sequences
    And sequence 0 has source "human_hoxa"
    And sequence 0 has start 100
    And sequence 0 has size 8
    And sequence 0 has strand "+"
    And sequence 0 has source size 100257
    And sequence 0 has text "ACA-TTACT"
    And sequence 1 has strand "-"

