## Test environments

* local Debian GNU/Linux bullseye/sid R 4.0.0
* local ubuntu  18.04 R 3.6.1
* local windows 8 (9200) R 3.5.2
* Fedora Linux, R-devel, clang, gfortran via rhub
* Ubuntu Linux 16.04 LTS, R-release, GCC via rhub
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 4  NOTES found:

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Carlos Catania <harpomaxx@gmail.com>’

New submission

Found the following (possibly) invalid URLs:
  URL: (https://arxiv.org/abs/1810.05691)
    From: man/cluster_algorithm.Rd
    Message: Invalid URI scheme

We have double checked and the above URL seems to be OK.

* checking package dependencies ... NOTE
Package suggested but not available for checking: ‘gpuR’

A: The package gpuR for the moment has some issues for Windows, however is not necesary for running galgoR

* checking Rd cross-references ... NOTE
Package unavailable to check Rd xrefs: ‘gpuR’

* checking examples ... NOTE
Examples with CPU or elapsed time > 5s
                       user system elapsed
classify_multiple     4.968  0.184  35.988
create_centroids      3.960  0.176  38.356
non_dominated_summary 2.556  0.088  36.214
toList                1.064  0.096  34.006
toDataFrame           1.060  0.092  32.206
galgo                 1.056  0.088  33.225

A:the galgo functions runs a genetic algorithm that demans considerable time. We have tried to reduce examples' time to a mininum.


