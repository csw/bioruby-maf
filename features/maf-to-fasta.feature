Feature: Convert MAF file to FASTA
  In order to use multiple alignment data with other tools
  I want to read a Multiple Alignment Format (MAF) file and write out its data as FASTA

  Scenario: Convert simple MAF file
    Given a MAF source file "t1.maf"
    When I select FASTA output
    And process the file
    Then the output should match "t1.fasta"
