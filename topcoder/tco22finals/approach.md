## Links
[Problem Statement](https://www.topcoder.com/challenges/30314404)

[Final Solution](https://github.com/FakePsyho/cpcontests/blob/master/topcoder/tco22finals/BlockGame.cpp) - Final submission minus dead code.

## Approach
My solution is a very basic greedy (with lookahead for T > 1) with a completely unmanagable eval function. 

Every turn, I simulate all possible placements of up to X tiles (where X <= T, i.e. I don't simulate future draws) and choose the one that maximizes the scoring function. For X > 1, I only consider moves where the Xth move makes a 2+ line clear, while the previous tiles didn't make any clears. Before the simulation, I estimate the number of possible states (upper bound on the possible number of different placements) and run the lookahead if the number of states is under a certain threshold. This is so I won't waste execution time on turns where there are too many possibilities.

If no moves are possible, I discard a tile with the most fruits on it.


## Evaluation Function
This is where it gets a little bit weird. The eval is slightly simplified (avoiding few micro-optimizations that had minor impact) in order to make it more readable.

* `Eval(state)` = `(row_bonus - hole_penalty + clear_bonus) * lines_mod + edge_bonus + tile_bonus - move_penalty * (X - 1)`, where X is the lookahead
* `row_bonus` = `dominant_fruit_cells ^ (3.8 - P * 3.5 + C * 0.05)`, summed over all rows/columns
* `hole_penalty` = `patterns_matched * (4.0 - C * 0.5)`, where `patterns_matched` is the count for all F.F patterns (empty cell surrounded by two fruits or fruit + edge of the grid); patterns counted on the board had slightly higher weight
* `clear_bonus` = `actual_score_for_the_clear * N * N * 0.2`
* `lines_mod` = `1.0 for no clear, 0.1 for single clear, LINES for clearing LINES simultenously`
* `edge_bonus` = `max(abs(row+1 - N*.5 - .5), abs(column+1 - N*.5 - .5)) * 1.4`, where row & column are the top-left corner of the placed tile
* `tile_bonus` = `N * number_of_fruits ^ 2`
* `move_penalty` = `N ^ 2.25`

There are two big downsides to this approach:
* There's no planning ahead. If I can't make a clear within X moves I'm choosing the next move without any future simulation. 
* Penalizing holes is too simple. The ultimate goal is to create a setup where with a single tile you can clear multiple lines, i.e. holes should be aligned in such a way that consecutive rows/columns have holes in the same 3x3 region.

My guess is that because of the above, I have rather low scores on big N / low P tests. T=1 is probably a bit weak as well, since the main advantage is the ability of searching through significant number of possible moves.


## Random Thoughts
This was a very unusual TCO final for me. For the past few months I'm struggling with long covid. While I'm close to being relatively healthy now, I figured that I probably have only like 6-8 effective hours instead of more expected 14-18. Thus, I approached the final as it was a 8-12h contest.

I believe I had my main solution (including lookahead up to 2 and most of the components of the scoring function) after around 2 hours. For the rest of the contest I just abused the fact that I'm getting 1K tests in 90-130 seconds, so I heavily tweaked the parameters without any plan. That's the main reason why the eval function is such a mess. I was just changing random parameters and looking for tiny improvements. I believe the parameters are near optimal as during the last 5 hours I have only climbed around 0.5-1%.

In the last several hours I was contemplating on writing a different type of solution that would specialize in low P / high T / high N cases: focus on making a cross pattern (while minimizing the holes) and leave one tile to finish the cross for a guaranteed 4-6 line clear. But, I wasn't able to figure out any easy way of doing it (and at the same time I felt too tired to implement something more ambitious). In hindsight, I should have focused on rewriting my "hole penalty" since this feels like the lowest hanging fruit. 

Overall, I think that was probably the most interesting problem we had in the final that used 24h format. The test cases were vastly different and allowed for wide range of approaches. At the same time, the problem felt rather light on the required implementation so there was a big room for experimentation. Also, a nice mix of randomness and provided knowledge. 

## Ablation Study
Results based on the first 1K seeds. That's the output from my [mmtester](https://github.com/FakePsyho/mmtester), that probably is the only reason why I'm in the top half of the leaderboard.
Each column represents a subset of the tests.

```
Tests             1000      215      232      207      153      193      197      227      226      176      174      213      209      211      182      185       167       177       177       168       165       147
Run            Overall    N=6-8   N=9-12  N=13-15  N=16-17  N=18-20      C=1      C=2      C=3      C=4      C=5      T=1      T=2      T=3      T=4      T=5  P=.2-.25  P=.25-.3  P=.3-.35  P=.35-.4  P=.4-.45  P=.45-.5
-------------  -------  -------  -------  -------  -------  -------  -------  -------  -------  -------  -------  -------  -------  -------  -------  -------  --------  --------  --------  --------  --------  --------
final          93.3245  93.7745  94.3903  93.4742  93.4524  91.2802  94.2712  93.1681  93.7148  92.9125  92.3667  94.5185  92.8429  93.1914  92.9582  93.0063   94.2820   94.3981   93.9381   92.6380   92.5316   91.9255
tl10x          94.5209  94.2671  95.0152  94.6437  94.8562  93.8118  94.6948  93.8725  94.7229  94.9897  94.4333  94.5185  92.8429  93.7191  95.1841  96.6813   95.4926   95.2169   95.4931   94.2024   93.1718   93.3214
tl.1x          88.7581  92.4481  91.2174  89.2441  86.6280  82.8586  90.0336  89.4523  89.6530  87.5711  86.4465  94.5185  92.8429  90.2456  84.2035  80.2952   90.1867   89.6177   88.0150   87.3826   88.8901   88.4635
noedge         93.1085  93.3991  93.7363  93.7386  92.9635  91.4694  94.2539  93.1875  93.3816  92.2386  92.2338  93.5722  92.9719  92.7324  93.5188  92.7544   94.0077   93.9713   94.1365   93.2970   91.8409   90.9970
notile         91.2269  93.0685  91.8089  90.8855  90.5079  89.4117  93.0587  91.2276  91.3162  90.4431  89.8287  94.5222  90.9863  88.6309  90.8432  91.0429   94.5373   93.6958   91.6028   90.4752   88.4841   88.0271
nomove         91.4339  93.2432  92.3019  91.3442  90.1922  89.4555  94.0793  92.4254  91.1489  89.8607  89.1067  94.5185  92.2849  91.2599  88.6039  89.9034   93.3599   93.8370   92.0288   90.2857   88.6305   90.1165
nohole         84.2173  85.0514  85.1931  83.8302  83.8631  82.8109  86.3884  84.3036  85.1703  82.6117  82.0327  86.9659  84.0387  83.0686  83.1836  83.5813   93.0657   90.1058   85.9980   83.4716   78.9415   71.7515
noholeweight   92.6119  92.7160  92.8070  92.9090  92.4147  92.0991  92.9359  93.6518  92.9528  92.2198  90.8422  93.9573  92.7317  92.1282  92.4379  91.6504   94.0732   93.8495   93.5211   92.0642   90.7738   91.1008
lookahead2     86.0055  83.3221  85.2247  86.9646  87.4881  87.7292  87.9276  86.1487  86.6448  85.5888  83.2335  94.5185  92.8429  85.2170  78.3904  76.8703   88.6088   88.1499   86.7553   84.5802   84.9894   82.1712
lookahead3     92.4673  91.8531  92.9734  92.6322  93.3166  91.6932  93.3075  91.9007  92.9631  92.3152  91.7653  94.5185  92.8429  93.1835  90.8961  90.4103   93.4752   94.0592   93.6335   91.5756   91.1228   90.4771
nolookahead    72.1608  65.0913  70.9594  74.6379  75.2951  76.3389  71.2865  73.6473  72.4050  73.3129  69.7290  94.5185  80.4182  67.9327  58.6280  55.2264   75.6772   74.1686   71.6416   68.6785   71.3789   70.9828
randomdiscard  92.5485  92.9395  93.2862  92.5307  92.2259  91.5011  93.5936  93.1582  92.2101  92.0982  91.4651  94.5185  92.1452  92.2142  91.8262  91.8279   94.0143   93.9938   93.8389   91.6722   91.0937   90.2495
rowpow2.725    89.5335  89.8004  89.3126  89.8110  89.4736  89.2514  92.9267  89.4766  89.0756  88.3319  87.5759  92.3233  91.2263  89.2247  87.7125  86.5527   85.5944   90.6506   92.0944   92.2690   90.3388   85.5951
```

* `final`: final submit
* `tl10x`: 10x time limit & 10x max state 
* `tl.1x`: .1x time limit & .1x max state
* `noedge`: `edge_bonus = 0`
* `notile`: `tile_bonus = 0`
* `nomove`: `move_bonus = 0`
* `nohole`: `hole_penalty = 0`
* `noholeweight`: no bonus weight for holes near the edges
* `nolookahead`: no lookahead
* `lookahead2`: lookahead (X) set to 2
* `lookahead3`: lookahead (X) set to 3
* `randomdiscard`: discard random tile instead of one with the most fruits 
* `rowpow2.725`: `row_bonus = dominant_fruit_cells ^ 2.725` (this removes fine tuning vs different C & P values)


## Example Results
```
{"id": 1, "score": 28135.0, "time": 1038, "dt": 5, "dc": 0.00912964, "m1": 12, "m2": 54, "m3": 65, "m4": 57, "m5": 13, "m6": 2}
{"id": 2, "score": 10792.0, "time": 140, "dt": 332, "dc": 0.39643, "m1": 116, "m2": 3, "m3": 1, "m4": 0, "m5": 0, "m6": 0}
{"id": 3, "score": 33876.0, "time": 189, "dt": 136, "dc": 0.198903, "m1": 29, "m2": 53, "m3": 1, "m4": 0, "m5": 0, "m6": 0}
{"id": 4, "score": 15009.0, "time": 246, "dt": 3, "dc": 0.0057732, "m1": 7, "m2": 94, "m3": 49, "m4": 20, "m5": 5, "m6": 0}
{"id": 5, "score": 68118.0, "time": 3856, "dt": 0, "dc": 0, "m1": 0, "m2": 16, "m3": 24, "m4": 5, "m5": 1, "m6": 0}
{"id": 6, "score": 28336.0, "time": 128, "dt": 32, "dc": 0.0566111, "m1": 18, "m2": 77, "m3": 10, "m4": 0, "m5": 0, "m6": 0}
{"id": 7, "score": 61631.0, "time": 2115, "dt": 0, "dc": 0, "m1": 1, "m2": 53, "m3": 25, "m4": 3, "m5": 0, "m6": 0}
{"id": 8, "score": 54013.0, "time": 6915, "dt": 11, "dc": 0.0173833, "m1": 0, "m2": 71, "m3": 16, "m4": 1, "m5": 1, "m6": 0}
{"id": 9, "score": 34584.0, "time": 82, "dt": 54, "dc": 0.096084, "m1": 31, "m2": 64, "m3": 4, "m4": 0, "m5": 0, "m6": 0}
{"id": 10, "score": 75520.0, "time": 113, "dt": 8, "dc": 0.0150041, "m1": 9, "m2": 67, "m3": 2, "m4": 0, "m5": 0, "m6": 0}
```

`dt` = dropped tiles; `dc` = dropped cells (total % of fruits in dropped tiles); `mX` = the number of clears with X lines
