# mod_seq
Scripts for analyzing high-throughput sequencing data generated from DMS modification of in vivo E. coli ribosomes.

If many sequences are being aligned, there is a potential issue with the number of files open. This can be fixed at the command line.

$ ulimit -a

open files                      (-n) 256

Use $ ulimit -Sn 10000  (or some other large value) to increase the number of open files.
