Feature: MAF output
  In order to output modified MAF files or subsets of them
  I want to be able to write out parsed MAF data

  Scenario: Reproduce simple test data
    Given a MAF source file "mm8_single.maf"
    When I open it with a MAF reader
    And open a new MAF writer
    And write the header from the original MAF file
    And write all the parsed blocks
    Then the output should match, except whitespace, "mm8_single.maf"

  Scenario: Reproduce test data with i, e, q lines
    Given a MAF source file "chr22_ieq.maf"
    When I enable the :parse_extended parser option
    And I enable the :parse_empty parser option
    And I open it with a MAF reader
    And open a new MAF writer
    And write the header from the original MAF file
    And write all the parsed blocks
    Then the output should match, except whitespace, "chr22_ieq.maf"
