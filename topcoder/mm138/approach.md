## Code

[Final Solution](https://github.com/FakePsyho/cpcontests/blob/master/topcoder/mm138/DiceRoller.cpp) - My final submission minus dead code.

## Solution
Since I only had half a day to do the contest, my solution is very simple & unoptimized Beam Search.

Scoring function: `(score + sum_in_connected_component * 0.5) * (B if there's a path to start, otherwise 1.0) - distance_to_edge * (2V+1)`

subtracting `distance_to_edge * (2V+1)` forced the path to go on a spiral.

The sides are calculated after the BS is done, so score is just the highest sum on each side. `sum_in_connected_component` is recomputed every step, so everything is quite slow. For comparison I get 4.6M state evals for seed=1 and 1.8M for seed=2.

Instead of doing 1 big BS, I do run 20 smaller ones. This reduces the issue of being stuck with bad values on sides, due to early unlucky few moves. Each run uses a different subset of starting positions.

## Example Results
```
Test Case #1:  Score = 55.46853322492484
Test Case #2:  Score = 1707.0038998728953
Test Case #3:  Score = 1301.7160949886
Test Case #4:  Score = 270.34505701228136
Test Case #5:  Score = 82.31861128802927
Test Case #6:  Score = 320.02144439905953
Test Case #7:  Score = 394.42886363867115
Test Case #8:  Score = 1406.2838808378442
Test Case #9:  Score = 1143.6227799868311
Test Case #10: Score = 1446.8247099797002
```





