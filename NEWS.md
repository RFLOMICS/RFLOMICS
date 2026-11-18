CHANGES IN VERSION 1.3.3
-------------------------
      
NEW FEATURES

    o Added option to filter missing values in abundance data (metabolomics and proteomics).
    o Added default edgeR-based filtering for count data.

SIGNIFICANT USER-VISIBLE CHANGES

    o Updated export of results, including archive and report outputs.

BUG FIXES

    o Updated compatibility with recent versions of clusterProfiler.
    
CHANGES IN VERSION 1.5.0
-------------------------

BUG FIXES

    o Replace expect_identical with expect_equal in tests for BioC 3.24 devel compatibility
    o Revise RAM warning and add note on resources
    