# (G)enetic (A)gortihm for Vehicle (R)outing (P)roblems

This project hopes to refresh my knowledge of heuristic's and meta heuristics for solving NP-complete problems. It also seeks to experiment with rust's advantages in the space of heuristic problems.

_A meta heuristic approach using genetic algorithm for the [vehicle routing problem](https://en.wikipedia.org/wiki/Vehicle_routing_problem)_

## This specific implementation uses:
- [ordered cross-over](https://en.wikipedia.org/wiki/Crossover_(genetic_algorithm)#Order_crossover_(OX1))
- [haversine distance](https://en.wikipedia.org/wiki/Haversine_formula) for distances
- [tournamenta](https://en.wikipedia.org/wiki/Tournament_selection) & [elitism](https://en.wikipedia.org/wiki/Selection_(genetic_algorithm)#Elitist_Selection) for selection
- [swaps](https://www.geeksforgeeks.org/mutation-algorithms-for-string-manipulation-ga/) for mutation

## Roadmap

1. Use Polars and investigate the vrp crate 
2. Investigate benchmarks
3. Set up a way to use _real world_ examples
4. Investigate if we can generate realistic examples

## Inspirations:

- GE - https://www.mdpi.com/2079-9292/10/24/3147
- state of the art - https://www.mdpi.com/2076-3417/11/21/10295
