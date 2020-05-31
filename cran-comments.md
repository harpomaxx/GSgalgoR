## Test environments

* local Debian GNU/Linux bullseye/sid R 4.0.0
* local ubuntu  18.04 R 3.6.1
* local windows 8 (9200) R 3.5.2
* win-builder (devel, release and old release)
* rhub
    - macOS 10.13.6 High Sierra, R-release, CRAN's setup
    - Windows Server 2008 R2 SP1, R-release, 32/64 bit 
    - Debian Linux, R-devel, clang, ISO-8859-15 locale

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 2 NOTES:

* Possibly mis-spelled words in DESCRIPTION:  galgo (9:50), omics (10:54), subtype (8:67), subtypes (12:46)

A: All of the identified words are spelled correctly.

* Found the following (possibly) invalid URLs: URL: https://arxiv.org/abs/1810.05691 From: man/cluster_algorithm.Rd Status: Error Message: libcurl error code 60: server certificate verification failed. CAfile: none CRLfile: none (Status without verification: OK) URL: https://doi.org/10.1371/journal.pone.0036530 From: man/use_rna_luad.Rd Status: Error Message: libcurl error code 60: server certificate verification failed. CAfile: none CRLfile: none (Status without verification: OK)

A: The flagged URLs are correct.







* checking package dependencies ... NOTE
Package suggested but not available for checking: ‘gpuR’

A: The package gpuR for the moment has some issues for Windows, however is not necesary for running galgoR

* checking examples ... NOTE
Examples with CPU (user + system) or elapsed time > 10s
            user system elapsed
plot_pareto 1.95   0.31   18.16

A:the galgo functions runs a genetic algorithm that demands considerable time. We have tried to reduce examples' time to a mininum.


