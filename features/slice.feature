Feature: MAF slicing
  In order to obtain just the alignment data covering a given region
  I want to be able to take slices of alignment blocks over
  A given interval

  Scenario: Interval covering two blocks
    Given a MAF source file "mm8_chr7_tiny.maf"
    And a Kyoto Cabinet index file "mm8_chr7_tiny.kct"
    When I open it with a MAF reader
    And I enable the :remove_gaps parser option
    And open a new MAF writer
    And write the header from the original MAF file
    And filter for only the species
      | mm8 |
      | rn4 |
    And search for blocks between positions 80082350 and 80082380 of mm8.chr7
    And slice the resulting blocks according to the given interval
    And write all the matched blocks
    Then the output should match, except whitespace, "mm8_chr7_tiny_slice1.maf"

  Scenario: Interval covering two blocks, using directory access, counting
    Given indexed MAF files in "test/data"
    When I enable the :remove_gaps parser option
    And filter for only the species
      | mm8 |
      | rn4 |
    And I extract a slice over the genomic interval
      | chrom    |    start |      end |
      | mm8.chr7 | 80082350 | 80082380 |
    Then 2 blocks are obtained

  Scenario: Interval covering two blocks, using directory access
    Given indexed MAF files in "test/data"
    When I enable the :remove_gaps parser option
    And open a new MAF writer
    And write a default header
    And filter for only the species
      | mm8 |
      | rn4 |
    And I extract a slice over the genomic interval
      | chrom    |    start |      end |
      | mm8.chr7 | 80082350 | 80082380 |
    And write all the matched blocks
    Then the output should match, except whitespace, "mm8_chr7_tiny_slice1.maf"
    
  Scenario: Interval in block subset
    Given indexed MAF files in "test/data"
    When I open a new MAF writer
    And write a default header
    And I extract a slice over the genomic interval
      | chrom    |    start |      end |
      | mm8.chr7 | 80082718 | 80082728 |
    And write all the matched blocks
    Then the output should match, except whitespace, "mm8_chr7_tiny_slice2.maf"
    
  Scenario: Interval to end of block
    Given indexed MAF files in "test/data"
    When I open a new MAF writer
    And write a default header
    And I extract a slice over the genomic interval
      | chrom    |    start |      end |
      | mm8.chr7 | 80082757 | 80082767 |
    And write all the matched blocks
    Then the output should match, except whitespace, "mm8_chr7_tiny_slice3.maf"
    
