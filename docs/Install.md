##Installation instructions

```
git clone https://github.com/Illumina/akt.git
cd akt/
make
```
If you get a warning about '<omp.h> not found' (especially on Macs)
you can try
```
make no_omp
```
Everything will be run on a single thread, so the `-n` option does nothing
in `akt kin` and `akt ibd`.
