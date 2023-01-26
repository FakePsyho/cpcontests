## Links
[Problem Statement](https://www.topcoder.com/challenges/30328975)

[Final Solution](https://github.com/FakePsyho/cpcontests/blob/master/topcoder/mm143/TreeMaker.cpp) - Final submission minus dead code. Beam search code is very clean, but `sim()` and `sim_score()` are hard to read due to optimizations and dynamic nature.

## Approach

My solution uses exactly N\*N moves and pushes all of the tiles from the same side. Every column (or row if pushes are sideways) will receive exactly N tiles. 

![](https://github.com/FakePsyho/cpcontests/blob/master/topcoder/mm143/seed1665.gif)

Whole solution is constructed via a single run of beam search that takes all of the available time. For the biggest test cases (N=30), my beam width is around 2500-3000. At every state I evaluate all of the possible N moves (unless column is completely full). Score evals are nearly constant time (with a very large constant) and I believe most of the execution time goes to state copying. I'm using union-find for calculating components.

Complexity is hidden within the beam search evaluation function. I'm using scoring function from the problem statement only at the very end, when selecting the final solution. 

## Eval Function

My eval function is a linear combination of all parts mentioned in the first list.

Component = set of connected tiles of a single color

Things that worked:

* Number of cells that are impossible to satisfy * 35; empty cells that have either two incoming edges from the same component or two incoming edges from different colors
* Broken connections * (20 + P); In order to minimize the number of components you have to minimize the number of broken connections, but broken connections have more immediate feedback
* Number of components with zero edges open * 17; Those are "dead" components that we can't extend anymore
* Number of components with a single edge open * 3; 
* Number of components * [0, 4, 8, 13] (based on C); Component weight is mostly to discourage placing tiles far apart
* Relative height difference: 1.2 * MSE (Mean Squared Difference); "Spiky" shape limits our options as cells in "concave" spots are going to be highly restricted; Note that in my final submission it's a little bit more complex, but 1.2 * MSE is almost as good as the final one
* Number of times when a tile is adjacent to an empty tile or tile of another color * 3; My intuition is that this helps avoiding the problem of creating "spiky" shapes, while also trying to keep colors close to each other


Things that didn't work:

* Ten different ways of trying to add penalty for spreading out tiles of different colors. Things like: pre-assigning columns for specific colors and giving penalties for placing tiles outside of those, counting inversions in each row, distance to the closest tile of the same color, forcing tiles into columns for the first X rows, etc. All of those changes drastically reduced my scores. 
* Bonus/Penalty for clearing out complete columns. My intuition was that clearing out complete column is bad because it limits our options. At the same clearing out complete column is good because it's hard to do that properly without incurring heavy penalty somewhere. Alternative that could've worked: giving a bonus for clearing out "outer" columns.
* Adding a temporary bonus based on the number of possible placements for the next tile.
* Customizing weights based on C/P/N. I found it surprising that my "locally optimal" set of weights was essentially independant of C/P/N parameters (with some very rare exceptions).


Promising things that I didn't try:

* Scoring how well the current board can accomodate all of the potential tiles and penalize situations when a lot of tiles cannot be put safely without immediately creating a bad situation.
* Scoring "effectiveness" of each component. While my eval tries to avoid creating situations where the cycle is forced, it has no way of penalizing cycles. Each cycle indirectly increases the number of components by 1. 

## Insert Facepalm Emoji

For nearly the whole duration of the contest, I have assumed that there's a hard limit on the number of moves and it's N\*N instead of 10K. I was very confused about the results, because I thought that with such limit, approach that I was using is probably the only viable one. And yet, for some reason it was clear (after making a submit for only C=1) that my solution has most of the bests of C=1, while it struggles with everything else. I just assumed that I'm missing some simple formula for C>1 tests that everyone else had.

I had "discovered" my mistake few hours before the submission deadline and in that short timespan I came up with 3 ideas that would utilize additional moves:

1. Perform some additional setup before the main beam search starts. I'm not sure what's the best strategy here, but the only one that worked for me was when I tried minimizing the number of outside edges. It gave a small ~.5% improvement. Unfortunately I wasn't convinced about adding 200-300 lines of code 20 minutes before the deadline, for a mere half a point gain.

2. Perform post-optimization after the main beam search is done. The last row in my solution is usually very bad / nearly random. If you perform UP c1 UP c2 DOWN c1 DOWN c2, you're essentially rotating the set of 3 tiles: [tile in hand, N-1:c1, N-1:c2]. This looked promising but it had 2 problems: Calculating components from scratch is very slow and while my solution could be improved by a lot, I'd probably need a proper beam search to find those improvements. Simple greedy was probably an .5% improvement again (but again, wasn't included in the final submit). With full beam search (or very clever heuristics) my guess is that it could be a further 2-3% gain.

3. During beam search, add an option to swap the current tile by performing perpendicular push that doesn't affect any formed rows. While this doesn't help with the poor "end game" of my solution, this should solve nearly all of the problems while the last row is still available. Sadly, I haven't tried this one, as it felt the most complex to implement. That being said, I wouldn't be surprised if that's a 5-15% improvement if properly implemented. 

## Ablation Study
Results based on the seeds 2001-3000. Output from [psytester](https://github.com/FakePsyho/psytester). For this ablation study, I have used beam search with fixed width in order to get comparable results among different runs.

The impact on C=1 is generally very small, because my score was already very close to N*N on many of the tests. 

```
Tests                 1000      275      237      243      245      195      171      223      194      217      497      503
Run                Overall      C=1      C=2      C=3      C=4   N=6-10  N=11-15  N=16-20  N=21-25  N=26-30    P=1-5   P=6-10  Bests  Uniques    Gain
-----------------  -------  -------  -------  -------  -------  -------  -------  -------  -------  -------  -------  -------  -----  -------  ------
final              93.6761  97.0913  92.8212  92.3060  92.0288  97.6062  94.0410  93.3293  92.0858  91.6351  94.0328  93.3237    256       75  0.2339
tlx0.1             82.4624  92.5623  80.0655  77.9337  77.9363  91.8811  85.1869  81.4593  78.0274  76.8476  83.2144  81.7195     45        1  0.0016
tlx0.5             90.9802  96.0663  89.5963  89.1053  88.4695  96.0267  91.8297  90.7453  88.9259  87.8536  91.5111  90.4556    147       29  0.0936
tlx1.2             94.7232  97.3884  93.9998  93.3462  93.7973  97.8423  95.3251  94.0402  93.9212  92.8651  95.1149  94.3363    286      106  0.3647
tlx2.0             96.5419  97.6742  96.1257  95.7783  96.4311  98.6600  96.5323  96.3953  96.2189  95.0857  96.7582  96.3282    448      231  1.0717
no_rotations       91.7919  95.2874  91.6842  89.3951  90.3497  91.7062  91.9655  92.8617  90.4755  91.8095  92.3220  91.2680    220      147  0.5937
no_hashing         88.6922  93.7506  86.1037  86.3542  87.8373  95.7272  90.0515  87.9828  85.6664  84.7334  89.1001  88.2891    123       20  0.0480
eval_og            47.4467  47.8335  39.2915  48.9033  53.4568  88.0282  64.0578  42.6530  28.7192  19.5584  45.6781  49.1942     28        0  0.0000
eval_no_relh       74.2523  95.4830  68.7305  64.9614  64.9785  94.5778  82.7033  72.7578  65.1190  59.0289  74.6082  73.9006    153       30  0.0764
eval_no_bad        78.1544  91.3060  71.9235  72.3227  75.2039  96.0966  88.1871  80.0144  69.0815  60.3250  77.7735  78.5307     92       14  0.0488
eval_no_01edge     86.5862  87.2276  83.6503  85.8217  89.4647  94.9255  89.6712  86.1723  82.9962  80.2962  87.2835  85.8973     90       24  0.0620
eval_no_adjacency  90.5303  94.7830  88.2132  89.3401  89.1787  97.5549  93.6518  89.4356  87.7074  85.4066  91.4758  89.5960    176       50  0.1677
```

* `final`: final submit
* `tl.1x`: 10% of final's time limit / beam width
* `tl.5x`: 50% of final's time limit / beam width
* `tl1.2x`: 120% of final's time limit / beam width; this is a maximum realistic improvement from code optimization alone
* `tl2x`: 200% of final's time limit / beam width
* `no_rotations`: beam search considers only a single direction of pushes instead of all 4
* `no_hashing`: no hashing for removal duplicate states
* `eval_og`: eval: default scoring function from problem statement
* `eval_no_relh`: eval: no relative height
* `eval_no_bad`: eval: no penalty for empty cells where it's impossible to place tiles without breaking connections
* `eval_no_01edge`: eval: no penalty based on compenents' open edges
* `eval_no_adjacency`: eval: removed penalty based on adjacency


## Example Results
```
{"id": 1, "score": 48.0, "time": 7286, "N": 6, "C": 1, "P": 1}
{"id": 2, "score": 2052.0, "time": 9737, "N": 30, "C": 4, "P": 2}
{"id": 3, "score": 727.0, "time": 9770, "N": 24, "C": 2, "P": 3}
{"id": 4, "score": 252.0, "time": 9684, "N": 11, "C": 2, "P": 10}
{"id": 5, "score": 57.0, "time": 9339, "N": 7, "C": 1, "P": 2}
{"id": 6, "score": 148.0, "time": 9725, "N": 12, "C": 1, "P": 3}
{"id": 7, "score": 498.0, "time": 9735, "N": 15, "C": 4, "P": 2}
{"id": 8, "score": 1228.0, "time": 9755, "N": 26, "C": 3, "P": 8}
{"id": 9, "score": 794.0, "time": 9740, "N": 28, "C": 1, "P": 2}
{"id": 10, "score": 682.0, "time": 9756, "N": 26, "C": 1, "P": 1}
```

