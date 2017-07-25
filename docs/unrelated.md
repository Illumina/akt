## unrelated

This takes the output from "akt kin" and creates a list of nominally unrelated individuals. 

**-k** *value* individuals with kinship coefficient > *value* are considered related (default 0.025)

**-i** *value* setting *value*>0 enables stochastic approach (default 0)

The algorithm has two options:

### simple greedy algorithm

1. Select individual with smallest number of relatives (defined as kinship coefficient > k) and remove all their relatives.
2. Repeat 1. until remaining individuals are unrelated.

### stochastic approach

For each sub-graph:

1. Randomly select individuals within each sub-graph and remove their relatives
2. Repeat 1. until all individuals are unrelated

Repeat this *i* times, storing the largest unconnected set found. If the stochastic approach yields a larger unconnected set than the greedy approach then that is returned, else the greedy result is returned.

Note this [maximal independent set problem](https://en.wikipedia.org/wiki/Maximal_independent_set) is NP-hard.
