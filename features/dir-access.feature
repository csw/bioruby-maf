Feature: Provide access to multiple MAF files in a directory
  In order to efficiently work with many MAF files
  We need to provide a convenient interface to them

  Scenario: Query for several chromosomes at once
    Given indexed MAF files in "test/data"
    When I query for the genomic intervals
    | chrom    | start    | end      |
    | mm8.chr7 | 80082580 | 80082612 |
    | mm8.chrM | 1400     | 1590     |
    Then 5 blocks are obtained

  Scenario: Apply block filters
    Given indexed MAF files in "test/data"
    When I filter for blocks with text size at most 200
    And I query for the genomic intervals
    | chrom    | start    | end      |
    | mm8.chr7 | 80082580 | 80082612 |
    | mm8.chrM | 1400     | 1590     |
    Then 3 blocks are obtained

  Scenario: Apply sequence filters
    Given indexed MAF files in "test/data"
    When I filter for only the species
    | mm8  |
    | rn4  |
    | hg18 |
    And I query for the genomic intervals
    | chrom    | start    | end      |
    | mm8.chr7 | 80082580 | 80082612 |
    | mm8.chrM | 1400     | 1590     |
    Then 5 blocks are obtained
    And block 0 has 3 sequences

  Scenario: Set parse options
    Given indexed MAF files in "test/data"
    When I enable the :remove_gaps parser option
    And I filter for only the species
    | mm8  |
    | rn4  |
    | hg18 |
    And I query for the genomic intervals
    | chrom    | start    | end      |
    | mm8.chr7 | 80082580 | 80082612 |
    Then 2 blocks are obtained
    And the text size of block 1 is 121
