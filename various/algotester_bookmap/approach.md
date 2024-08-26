### Links
[Problem Statement](https://algotester.com/en/ContestProblem/DisplayWithEditor/135403)

[Final Solution](https://github.com/FakePsyho/cpcontests/blob/master/various/algotester_bookmap/main.cpp) - one of the worst codes I have produced: very messy, lots of dead code, many functions and variables are incorrectly named, etc.

### Approach

The core observation was that, due to the test generation constraints, the alpha and beta parts of the price calculation formula can be ignored without sacrificing too much accuracy. Not only it drastically speed up score calculation, but it also makes the final price calculation for each buyer completely independent of other buyers.

My final approach is glued together from multiple parts, each solving particular subproblem.
* **A**) Single buyer no alpha and no beta: Optimize a single buyer using Hill Climbing (HC) by performing random moves. Eval ignores alpha and beta, so it works in O(T).
* **B**) All buyers no alpha and no beta: Same as **A** but optimize all buyers at once allowing moves between different buyers.
* **C**) All buyers full eval: Global HC that doesn't ignore alpha and beta. Since the full eval is very slow O(NT(N+T)), a small optimization here is to go over the buyers in increasing order and consider only partial scores up to the current buyer. This reduces the sim calculation time down to O(DT(N+T)) where D is the maximum difference between buyers indices. In other words, if I'm currently at buyer B, I'm only considering moves between buyers in [B-D+1, B] range.
* **D**) Global allocation of boxes between buyers: If you can estimate the best possible score for a particular number of boxes given to a buyer, you can use very naive Simulated Annealing (SA) to solve the allocation problem (how many boxes are given to each buyer) nearly optimally. If I have several data points (number of boxes -> score) for each buyer, I can estimate the score for every possible number of boxes by performing linear interpolation between two data points.

Almost all of the difficulty in the problem comes from not having enough boxes to properly influence the prices. At the same time, it's usually possible to dump excessive boxes without much trouble. Because of this, we can estimate the difficulty of test by taking sum(L) / C ratio. The higher the ratio, the harder the test is and it's less likely to get a "perfect" 9999999.

For easy tests with sum(L) / C < 8.0:
* For the first 50% of the time, generate lots of random solutions using **B**, each run uses 2% of total time limit or if 5 * T * N moves happened without any improvement. Evaluate all of the solutions on full eval (with alpha & beta) and take the best one.
* Run **C** for the remaining time.

For hard tests with sum(L) / C > 8.0:
* Generate data points for **D**. Each data point is a single run of **A** with specified buyer and number of boxes. For each buyer, perform binary search (using 10 steps) to find the approximate number of boxes needed to "solve" it (get reasonably low score). Let's call this number X. Add three more data points: X/3, 2*X/3, X+5. This gives a total of 13 data points for each buyer. This exploits the fact that usually (always?) we can easily dump excessive number of boxes. As soon as we know the minimum number of boxes that gives us a very low score, we can assume that more boxes won't increase the score.
* Run **D** using generated data points. It's very fast as score evals take O(1). 
* Using generated box allocations, run **A** for 5% of total time, then **B** for 10% of total time and then finally **C** till the end.


### Possible improvements

* Some of the easier tests would benefit from **D** as well, but it would either require drastically speeding up gathering data points or some quick test to find out that the particular test is harder. A possible approach would be to perform a few runs of **B** (as normal), and then, if none them produced very low scores, switch to the other solution. Unfortunately, I didn't have enough time to try this.

* As mentioned above, generating data points is a major bottleneck. For the biggest tests (N=100, T=100), it's barely able to work in time. I'm pretty sure that finding a more effective method than **A** should be possible. Either a more analytical (or greedy?) solution or maybe something more clever than performing random transitions, considering the "optimal" solution usually have mostly zeros in most time slots. 

