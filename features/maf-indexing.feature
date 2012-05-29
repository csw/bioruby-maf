@milestone_2
Feature: Indexed access to MAF files
  In order to extract alignment blocks from MAF files
  By chromosomal ranges matching a source sequence
  I want to have a way to build indexes on MAF files
  And use indexes to efficiently find alignment blocks
  Because linear searches of a 200 GB file are impractical

  Scenario: Index a MAF file
    Given a MAF source file "mm8_chr7_tiny.maf"
    When I open it with a MAF reader
    And build an index on the reference sequence in "mm8_chr7_tiny.maf.index"
    Then the index has 8 entries
    
  Scenario: Extract alignment blocks by chromosomal range
    Given a MAF source file "mm8_chr7_tiny.maf"
    When I open it with a MAF reader
    And search for blocks between positions 80082592 and 80082766 of mm8.chr7
    Then 2 sequences are obtained
    And sequence 0 has start 80082592
    And sequence 1 has start 80082713
