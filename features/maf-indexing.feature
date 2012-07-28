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

  Scenario: Extract alignment blocks by chromosomal range on non-ref sequence
    Given a MAF source file "mm8_chr7_tiny.maf"
    When I open it with a MAF reader
    And build an index on all sequences
    And search for blocks between positions 136011819 and 136012026 of rn4.chr1
    Then 2 blocks are obtained
    And sequence mm8.chr7 of block 0 has start 80082368
    And sequence mm8.chr7 of block 1 has start 80082471

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
  Scenario: Build MAF index on BGZF file with CLI tool
    Given test files:
    | mm8.chrM.maf.gz |
    When I run `maf_index mm8.chrM.maf.gz mm8.chrM.kct`
    Then it should pass with:
    """
    """
    And a file named "mm8.chrM.kct" should exist

  @no_jruby
  Scenario: Build MAF index on all sequences with CLI tool
    Given test files:
    | mm8_chr7_tiny.maf |
    When I run `maf_index --all mm8_chr7_tiny.maf mm8_chr7_tiny.kct`
    And I run `maf_index -d mm8_chr7_tiny.kct`
    Then it should pass with regex:
    """
    9 \[bin 585\] 594:631
    """

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


