
## Links
[Problem Statement](https://www.topcoder.com/challenges/30278915)

[Final Solution](https://github.com/FakePsyho/cpcontests/blob/master/topcoder/mm138/DiceRoller.cpp) - My final submission minus dead code.

[Post Contest Solution](https://github.com/FakePsyho/cpcontests/blob/master/topcoder/mm138/DiceRoller_post.cpp) - Final Solution + ideas from RafbillFr solution (more info below)

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

## Post Contest Solution
Out of curiosity, I tried adding ideas [mentioned by RafbillFr](https://discussions.topcoder.com/discussion/20051/post-your-approach) (forum is blocked if you haven't registered for the contest). This was very easy, since both of our solutions had the same core concept (Beam Search).

Turns out that adding "G" component (multiplying each cell value by `G[number_of_neighbors]` where G[2]=0.13, G[3]=0.45, G[4]=0.74 was much smarter than my blanket 0.5. This change alone gave me around 8.5 point improvement. I have also added final run with fixed sides (+0.5 point) & bonus for cells that are unused, but are not present on the dice (another +0.5).

I don't do any proper execution time handling in my code, so it's hard to say how it would be perform on TC servers. That being said, according to my local tests this version would be very close (~0.5 point) below RafbillFR. Considering mine solution is ~10x slower, I'd count it as a success. 

Scores on the first 10 seeds (rounded to integers, as I never downloaded fixed local tester):
```
{"id": 1, "score": 55.0}
{"id": 2, "score": 2080.0}
{"id": 3, "score": 1495.0}
{"id": 4, "score": 297.0}
{"id": 5, "score": 82.0}
{"id": 6, "score": 339.0}
{"id": 7, "score": 432.0}
{"id": 8, "score": 1536.0}
{"id": 9, "score": 1275.0}
{"id": 10, "score": 1560.0}
```

