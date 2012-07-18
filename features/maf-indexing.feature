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
    And build an index on the reference sequence
    Then the index has at least 8 entries
    
  Scenario: Extract alignment blocks by chromosomal range
    Given a MAF source file "mm8_chr7_tiny.maf"
    When I open it with a MAF reader
    And build an index on the reference sequence
    And search for blocks between positions 80082592 and 80082766 of mm8.chr7
    Then 2 blocks are obtained
    And sequence mm8.chr7 of block 0 has start 80082592
    And sequence mm8.chr7 of block 1 has start 80082713
    
  Scenario: Extract alignment blocks by chromosomal range from index file
    Given a MAF source file "mm8_chr7_tiny.maf"
    And a Kyoto Cabinet index file "mm8_chr7_tiny.kct"
    When I open it with a MAF reader
    And search for blocks between positions 80082592 and 80082766 of mm8.chr7
    Then 2 blocks are obtained
    And sequence mm8.chr7 of block 0 has start 80082592
    And sequence mm8.chr7 of block 1 has start 80082713

  @no_jruby
  Scenario: Build MAF index with CLI tool
    Given test files:
    | mm8_chr7_tiny.maf |
    When I run `maf_index mm8_chr7_tiny.maf mm8_chr7_tiny.kct`
    Then it should pass with:
    """
    """
    And a file named "mm8_chr7_tiny.kct" should exist

  @no_jruby
  Scenario: Dump MAF index with CLI tool
    Given test files:
    | mm8_chr7_tiny.maf |
    | mm8_chr7_tiny.kct |
    When I run `maf_index -d mm8_chr7_tiny.kct`
    Then it should pass with regex:
    """
    0 \[bin 1195\] 80082334:80082368
    """

  