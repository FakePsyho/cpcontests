
## Code
[Final Solution](https://github.com/FakePsyho/cpcontests/blob/master/atcoder/ahc011/main.cpp) - If you ignore debug/stat gathering, it should be very readable as I have refactored the code near the end of the contest.

## Approach
There are three main parts:
* (1) Generate valid tilings (placement of tiles that create single connected component)
	- Simple unpolished HC with kicks (forced changes when stuck)
	- Transition: swap two random pieces
	- Eval: ∑ **component_size** * log(**component_size**+1)
	- ~200 tilings for N=6 and ~30 for N=10
	- Tilings are completely random.
* (2) Create matching between initial tiles and their final position
	- Weighted bipartite matching via MCMF (Min-Cost Max Flow)
	- Weight: Euclidean squared (it was a significant improvement over Manhattan/Euclidean)
* (3) Based on the matching, find the sequence of movements
	- Single beam search (BS) with fixed width. Unfotunately, it fails pretty often for N=9 & N=10
	- Eval: Euclidean squared
	- BS is fairly optimized in terms of memory usage / execution time

Those parts are run in sequential order. (2) is repeated until a valid matching is found (50% of the matchings are incorrect due to parity). (3) is run mutiple times (each time for a different matching since those are very fast to calculate) since (3) will often fail for bigger tests. IIRC the success rate for N=10 WIDTH=1000 was around 10-15%.


## Comments
* There's a leftover backup code for tile-by-tile solver. It's used as a backup solver in case I won't find a good solution for N=9..10. It solves the edges in order to reduce the grid to 8x8 and then normal BS solves the rest. 
* While the approach is simple conceptually, there's a lack of synergy between the different parts. Selecting random tilings is a particular poor use of resources.
* My guess is that improving eval function in BS would drastically improve the quality of results. BS was not able to solve a single case for N=10 with manhattan distance. 95% of bad runs end up in almost fully solved grid except for 1-3 groups with a 2-5 cycle of permuted tiles. I'm almost sure there's a clever way of fixing that.
* There was a significant correlation between the value of the matching and the number of moves that BS generated. Easy way of getting another 0.5M for my solution would be to generate more tilings, filter them and run BS only for the best ones.
