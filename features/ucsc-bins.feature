Feature: Computation of UCSC bins
  In order to efficiently use indexes
  We will use the UCSC bin indexing system
  Per http://genomewiki.ucsc.edu/index.php/Bin_indexing_system

  Scenario Outline: Compute smallest containing bin
    Given I have a region with start <Start> and end <End>
    When I compute the smallest containing bin
    Then the bin should be <Bin>

    Examples:
    |    Start |      End |  Bin |
    | 25079603 | 25079787 |  776 |
    | 25128173 | 25128248 |  776 |
    | 50312474 | 50312703 |  968 |
    | 41905591 | 41906101 |  904 |
    | 16670899 | 16673060 |  712 |
    | 75495356 | 75495494 | 1160 |
    | 92259501 | 92261053 | 1288 |
    | 83834063 | 83838132 | 1224 |
    |  7309597 |  7310411 |  640 |
    |  6190410 |  6190999 |  632 |
    # from https://github.com/polyatail/biopython/blob/af34c033d78c4c72dffbb500e513e568a2ba5e29/Tests/test_MafIO_index.py#L48
    
