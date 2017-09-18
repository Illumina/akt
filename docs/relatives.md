## relatives

Takes the output from `akt kin` and detects/reconstructs pedigrees from the information. Can also flag duplicated samples and create lists of unrelated samples.

**-k** *VALUE* Only keep links with kinship above this threshold (searches in this set for duplicate, parent-child and sibling links).  
**-i** *VALUE* Iteration parameter for unrelated set output.  
**-p** *PREFIX* Prefix for output files.  
**-g** If present output graphviz files. These can be visualised using e.g. `neato -Tpng -O out.allgraph` or for family pedigrees `dot -Tpng -O out.Fam0.graph`.  

```
./akt relatives allibd -g > allrelatives
```

The output contains duplicates, families and relationship types.

```
grep ^Dup allrelatives
Dup0 Sample0
Dup0 Sample1
...
grep ^Fam allrelatives
Fam0 Sample2
Fam0 Sample3
...
...
grep ^Type allrelatives
Type Fam0 Sample2 Sample3 Parent/Child
...
grep ^Unrel allrelatives
Sample0
Sample2
...
```

The file `out.allgraph` can be viewed with gviz e.g. `fdp out.allgraph -Tpng -O` and the families can be viewed using 
e.g. `dot out.Fam0.graph -Tpng -O`. The parent child relations are also recorded in PLINK fam format in `out.fam`. If
e.g. a sibling pair, is found the samples will appear in `out.fam` without parents. If the direction of the relationship
can't be determined e.g. for parent/child duos a random sample is assigned to be the parent in `out.fam`. The final column
in the `.fam` file specifies how many potential parents the sample had.

Note that `relatives` is quite a aggressive in its pedigree search, and can make errors when founders 
are missing (for example a mother and two children). We can remove false pedigrees via a simple Mendel consistency check:

```
akt kin --force -M 1 test.bcf > kinship.txt
akt relatives kinship.txt
akt mendel -p out.fam test.bcf > mendel.txt
python ~/workspace/akt/scripts/check_pedigree.py -fam out.fam -m mendel.txt > corrected.fam
``` 
