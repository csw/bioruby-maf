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

