## Links
[Problem Statement](https://www.topcoder.com/challenges/30285365)

[Final Solution](https://github.com/FakePsyho/cpcontests/blob/master/topcoder/mm139/PipeConnector.cpp) - My final submission minus dead code.

## Solution
My solution is a very simple application of Simulated Annealing. I'd imagine that most (if not all) solutions are SA-based.

The SA is very straightforward:
- State: All of the connections (in the form of paths)
- Scoring function: Same as the problem statement
- Transition: First choose a color, then try to remove random 0-2 paths (where paths of the same color have their probability drastically increased) and then try to add 1-3 random paths of the same color
- Temp schedule is exponential from 30 to 4. It seemed that N/C/P didn't matter too much for choosing correct temperatures.

Finding paths:
- When searching for paths I randomize a starting point (where the higher value points have higher chance of being chosen) and I ran a BFS to all of the remaining points of the color same color. 
- Penalty for collision is randomized for each BFS (range = [1, 1+N*P])
- End point is chosen via evaluation function: `point_value^3 / sqrt(N+path_len) * sqrt(random())` where `random()` has uniform [0..1] distribution
- The order of the BFS (i.e. up, down, left & right) is randomized, but it's fixed for the entire BFS.

I experimented a bit with adding penalty to my scoring function for each visited cell (where this penalty would be gradually reduced and hit 0 at the end). Unfortunately, I wasn't able to get any improvement out of that.

Giving 5x time is only a 0.15 improvement, so clearly I was missing quite a lot.

## Example Results
```
{"id": 1, "score": 309.0, "time": 9823, "step": 10713856, "acc": 5193766}
{"id": 2, "score": 1592.0, "time": 9803, "step": 4552448, "acc": 1153718}
{"id": 3, "score": 2676.0, "time": 9826, "step": 4681984, "acc": 1174857}
{"id": 4, "score": 3403.0, "time": 9803, "step": 5960960, "acc": 1506137}
{"id": 5, "score": 3227.0, "time": 9826, "step": 5434368, "acc": 1902134}
{"id": 6, "score": 589.0, "time": 9826, "step": 5355776, "acc": 1618754}
{"id": 7, "score": 1393.0, "time": 9825, "step": 6110208, "acc": 2274100}
{"id": 8, "score": 1153.0, "time": 9824, "step": 9171200, "acc": 1698628}
{"id": 9, "score": 2249.0, "time": 9824, "step": 5746176, "acc": 2262110}
{"id": 10, "score": 2022.0, "time": 9826, "step": 5201920, "acc": 1865512}
```

`step` is the number of transitions, `acc` number of accepted changes

