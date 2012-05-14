Feature: Convert MAF file to FASTA
  In order to use multiple alignment data with other tools
  I want to read a Multiple Alignment Format (MAF) file and write out its data as FASTA

  Scenario: Convert simple MAF file
    Given a MAF source file "t1.maf"
    When I select FASTA output
    And process the file
    Then the output should match "t1.fasta"

  Scenario: Convert simple MAF data
    Given MAF data:
    """
    ##maf version=1 scoring=humor.v4
    # humor.v4 R=30 M=10 /cluster/data/hg15/bed/blastz.mm3/axtNet300/chr1.maf
    # /cluster/data/hg15/bed/blastz.rn3/axtNet300/chr1.maf

    a score=0.128
    s human_hoxa 100  8 + 100257 ACA-TTACT
    s horse_hoxa 120  9 -  98892 ACAATTGCT
    s fugu_hoxa   88  7  + 90788 ACA--TGCT


    a score=0.071
    s human_unc 9077 8 + 10998 ACAGTATT
    # Comment
    s horse_unc 4555 6 -  5099 ACA--ATT
    s fugu_unc  4000 4 +  4038 AC----TT
    """
    When I select FASTA output
    And process the file
    Then the output should be:
    """
    >human_hoxa:100-108
    ACA-TTACT
    >horse_hoxa:120-129
    ACAATTGCT
    >fugu_hoxa:88-95
    ACA--TGCT
    >human_unc:9077-9085
    ACAGTATT
    >horse_unc:4555-4561
    ACA--ATT
    >fugu_unc:4000-4004
    AC----TT
    """

