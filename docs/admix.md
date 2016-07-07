##admix

The `admix` function attempts to use the output of `akt pca` to assign admixture fractions to data based on known populations.

**-c** *VALUE* Which columns in input file to use.
**-C** *VALUE* file with vectors whose transformation to admixture fractions is known.

The input file should be in the form
```
e1 e2 e3 ... ed POP a1 a2 a3 ... ad
...
```
where `e1...` are coordinates in the space of principal components and `a1...` are coordinates in the space
of admixture fractions. The mapping from PCA space to admixture space is
```
Ta + v = e
```
with inverse
```
a = T^{-1}(e-v)
```
`admix` calculates this inverse mapping.

```
./akt admix pcadata -c 2-6 -C pca_to_admix.txt > admixtures
```

The file `admixtures` contains

```
SAMPLE_ID %POP1 %POP2 ... %POPD
...
```

A vector must be chosen to fix the zero point of the transformation (the vector `v` above). 
For K populations the `-C` input file can contain K or K+1 lines.
If K+1 vectors are given the one with the lowest norm in the space of admixtures equals `v`.
If K vectors and their images in admixture space are specified then the zero point is fixed by 
choosing an arbitrary vector orthogonal to the sides of the tetrahedron made by the K inputs
and roughly equidistant from them. There are two possible choices of sign: the positive one is always chosen.



