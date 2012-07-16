Feature: Join adjacent alignment blocks
  After filtering out sequences
  The sequence that caused two blocks to be separate may be removed
  So it can be desirable to join such blocks together

  Scenario: Two blocks natively in indexed access
    Given indexed MAF files in "test/data"
    When I query for the genomic intervals
    | chrom    | start    | end      |
    | mm8.chr7 | 80082334 | 80082471 |
    Then 2 blocks are obtained
    And the text size of block 0 is 54
    And the text size of block 1 is 156

  Scenario: Two blocks joined in indexed access
    Given indexed MAF files in "test/data"
    When I enable the :join_blocks parser option
    And I filter for only the species
    | mm8     |
    | rn4     |
    | oryCun1 |
    | hg18    |
    | panTro2 |
    | rheMac2 |
    | canFam2 |
    | loxAfr1 |
    | echTel1 |
    And I query for the genomic intervals
    | chrom    | start    | end      |
    | mm8.chr7 | 80082334 | 80082471 |
    Then 1 block is obtained
    And the text size of block 0 is 210
