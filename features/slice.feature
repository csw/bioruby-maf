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
