## Links
[Problem Statement](https://www.topcoder.com/challenges/9e14555a-a853-4137-8902-cd118f7ea5dc)

[Final Solution](https://github.com/FakePsyho/cpcontests/blob/master/topcoder/mm146/HappyGrid.cpp) - Final submission. Rather bloated due to excessive redundancy.


## Overall Approach

The approach is split into two parts: Building target board and finding moves that transform original board into target board.

Creating target board is handled by Simulated Annealing (SA) while finding the moves is done by a greedy method that moves balls one by one in a predefined order of cells.

The solution itself consists of 3 separate phases.

Phase 1: SA on the target board that tries to maximize the number of connected components (CCs) equal to K.

Phase 2: SA on the target board that tries to maximize the score 

Phase 3: HC on the predefined order of cells in order to minimize the number of moves (turns)


## Building target board

**Phase 1:**

State: board (with no additional data)

Transition: 1) 5%: Swap two random cells 2) 95%: Swap a random cell and a cell that is connected to that CC

Eval: for each CC: 1.0 if size == K, -0.1 if size < K, K-size if size > K

**Phase 2:**

State: board (with no additional data)

Transition: 1) 4%: Swap two random cells 2) 76%: Swap a random cell and a cell that is next to the CC of the original cell 3) 20%: Random cell A, Two cells B&C that are next to the CC of A. Perform A->B->C->A rotation

Eval: score from the problem statement (i.e. `CCs_of_size_K * P - moves`)

The sizes of CCs are calculated dynamically after each transition. This is fairly quick since the problem forces small CCs.

SA is split into two phases because (a) eval function used in the final scoring does not encourage enough exploration and (b) dynamic CC calculation is fast (few short BFS runs) while calculating `moves` requires me to simulate the whole process which is `O(board_size)`. For Phase 2, only the transitions that do not reduce the number of CCs are taken into consideration.


## Finding moves

Finding moves is performed via a rather simple greedy:

* Mark a single cell (closest to the board center) as the origin and then order the remaining cells by the distance from the marked cell.

* Starting with the cell that is the furthest away from the origin: find the closest ball that has the correct color and move it to that target cell (via a chain of swaps)

* Repeat until the board is solved. Sorting cells by the distance from origin ensures (I think) that the "unsolved" part of the board is never disconnected which means that greedy is able to find correct set of moves for every target board.

This method is deterministic and is used in Phase 2 for calculating the number of turns to perform all of the moves.

Phase 3 uses a slightly more generalized (and slower) greedy where I'm allowed to modify the order of the cells that I'm solving.


## Details & Comments

* My solution is still very far from the optimum. As you can see from the attached ablation study, increasing the time limit twice, increases the score by 0.8. Despite that, I couldn't tweak the constants to gain even 0.1 more.

* On the very last day, I tried to add beam search for finding moves. Unfortunately, it turned out that achieving decent results with beam search is far from trivial and I was unable to get results even remotely comparable to the greedy method that I'm using.

* At some point I ran a flawed experiment that convinced me that maximizing the number of CCs is often not desired. My (probably incorrect) intuition was that this is due to the limited possibilities for useful transitions, so it's easier to get stuck in a local optimum. Because of that, my solution doesn't always select the board with the highest CC count. 


## Ablation Study

Results based on the seeds 10001-12000. Output from [psytester](https://github.com/FakePsyho/psytester). 

You can see the impact of various features of my solution. r is the P/K ratio.

Note that the variance on the 2K tests is significant and it was very hard to test tiny improvements. Changing seed alone could move my score by +/- 0.05.


```
Tests                   2000      709      706      585      396      405      415      397      387      290      289      298      296      280      284      263      681         694      625
Run                  Overall   N=8-15  N=16-23  N=24-30      C=2      C=3      C=4      C=5      C=6      K=2      K=3      K=4      K=5      K=6      K=7      K=8   r=-6.7  r=6.7-9.51  r=9.51-  Bests  Uniques    Gain
-------------------  -------  -------  -------  -------  -------  -------  -------  -------  -------  -------  -------  -------  -------  -------  -------  -------  -------  ----------  -------  -----  -------  ------
final                97.2766  97.6808  97.2694  96.7955  95.9647  98.1452  97.7736  97.5479  96.8989  97.4596  97.9259  97.7160  97.2993  97.3430  96.6133  96.4835  96.5169     97.3952  97.9727     95       47  0.0112
tl.2x                94.3521  96.1770  93.9893  92.5782  94.9487  95.9672  94.7160  93.9403  92.0838  96.3971  96.1883  95.5001  94.3093  93.7054  92.6903  91.3100  93.2060     94.4889  95.4490     26        4  0.0030
tl.5x                96.3438  97.2502  96.2579  95.3489  95.8242  97.2813  96.8430  96.2912  95.4130  97.0957  97.3229  97.0035  96.1295  96.1806  95.5733  94.9381  95.3991     96.5183  97.1792     63       28  0.0106
tl2x                 98.0700  98.0180  98.1604  98.0239  96.3706  98.6716  98.6339  98.5504  98.0819  97.8115  98.4290  98.1821  97.9718  98.3439  98.0251  97.7009  97.4911     98.1509  98.6109    258      188  0.0561
tl5x                 98.9411  98.5191  99.0289  99.3467  96.7220  99.3701  99.5152  99.5047  99.5692  98.1364  98.7764  98.9991  98.8794  99.2996  99.3936  99.1431  98.6217     98.9920  99.2327    910      808  0.5470
corner               96.7694  97.3977  96.6709  96.1266  95.7874  97.8004  97.3754  96.9119  95.8990  97.2443  97.5170  97.0412  96.7755  96.6360  96.1593  95.9100  95.8101     97.0312  97.5238    136      108  0.0479
no_p1                94.8828  96.5817  94.4863  93.3024  96.0306  97.6082  95.9021  93.7017  90.9747  97.4806  97.9304  97.9757  96.5872  94.1347  91.3521  87.8559  94.6979     94.7396  95.2434    132       82  0.0320
p1_eval_1            95.6234  96.9322  95.3511  94.3659  96.1011  97.0280  96.2124  95.1625  93.5061  97.4412  97.8793  96.5946  96.0639  95.1448  94.0658  91.7356  94.9653     95.7054  96.2496    106       60  0.0248
no_p3                96.3482  96.7821  96.1787  96.0270  94.9174  97.1762  97.0028  96.6594  95.9247  96.4878  96.8601  96.5267  96.2333  96.6157  95.9301  95.7257  95.1110     96.6427  97.3694     15        4  0.0011
p2_hc                87.1470  87.1543  86.9434  87.3839  79.5464  92.6083  90.1184  87.8502  85.3013  91.9771  85.4246  89.7045  87.2748  86.7626  85.0221  83.3758  84.5437     87.8017  89.2565      2        2  0.0016
maxcomp              97.3088  97.6905  97.3537  96.7920  96.0893  98.0313  97.7614  97.6266  96.9891  97.4383  97.8792  97.5799  97.1718  97.4877  96.8120  96.7321  96.4948     97.4668  98.0202    111       64  0.0172
p2_no_triple         97.1983  97.5906  97.1818  96.7427  96.0936  97.6573  97.7014  97.4559  97.0445  97.0712  97.4808  97.6682  97.1605  97.4578  96.7966  96.6954  96.4837     97.3131  97.8494     97       46  0.0104
p2_allow_reduce_cc   95.2411  97.0978  95.0677  93.2004  97.0574  96.8257  95.3573  94.1718  92.6967  97.9359  96.6557  95.9411  94.7790  94.3162  93.6275  93.1696  93.7523     95.5713  96.4967    171      138  0.1706
p1p2_random          93.2309  96.3224  92.7663  90.0446  94.4702  94.7848  93.5794  92.8238  90.3804  95.8599  96.0065  92.5574  92.0139  92.4430  92.2521  91.3102  91.0180     93.5112  95.3307     46       30  0.0193
no_dynamic_cc        94.7133  96.7556  94.4966  92.4998  94.9616  95.9587  95.2749  94.3655  92.9107  95.9866  95.9076  95.7375  94.8264  94.2677  93.5314  92.4602  93.6867     94.8572  95.6722     52       22  0.0073
```

* `final`: final submit
* `tl.2x`: 20% of final's time limit
* `tl.5x`: 50% of final's time limit
* `tl2x`: 200% of final's time limit
* `tl5x`: 500% of final's time limit
* `corner`: Order of the cells is created by running BFS from the corner instead of middle of the board
* `no_p1`: Phase 1 skipped
* `p1_eval_1`: Simplified Phase 1 eval: CC has value 1 if it's equal to K, otherwise it's 0
* `no_p3`: Phase 3 skipped
* `p2_hc`: Phase 2 uses HC instead of SA
* `maxcomp`: After Phase 1 is done, take the solution with the maximum number of CCs
* `p2_no_triple`: Phase 2 has A->B->C->A transition removed
* `p2_allow_reduce_cc`: Allow reducing the number of CCs in Phase 2. This efffectively removes an effective pruning method and reduces the number of iterations.
* `p1p2_random`: Phase 1 & 2 have only uniform random swaps as transitions. Drastically reduces the proportion of useful transitions.
* `no_dynamic_cc`: Remove dynamic version of recalculating CCs. Drastically reduces the number of iterations for larger boards.



## Example Results

```
{"id": 1, "score": 372.0, "time": 9861, "N": 8, "C": 3, "K": 7, "P": 57, "r": 8.14286, "p1_steps": 4096, "p2_steps": 8736768, "p3_steps": 1107968, "comps": 7, "turns": 27, "max_comps": 7}
{"id": 2, "score": 7440.0, "time": 9872, "N": 30, "C": 6, "K": 2, "P": 22, "r": 11, "p1_steps": 1, "p2_steps": 4145152, "p3_steps": 63488, "comps": 374, "turns": 788, "max_comps": 384}
{"id": 3, "score": 6221.0, "time": 9868, "N": 27, "C": 3, "K": 5, "P": 57, "r": 11.4, "p1_steps": 1361920, "p2_steps": 4747264, "p3_steps": 86016, "comps": 118, "turns": 505, "max_comps": 123}
{"id": 4, "score": 3617.0, "time": 9843, "N": 29, "C": 4, "K": 8, "P": 54, "r": 6.75, "p1_steps": 946176, "p2_steps": 2872320, "p3_steps": 78848, "comps": 80, "turns": 703, "max_comps": 87}
{"id": 5, "score": 5275.0, "time": 9845, "N": 28, "C": 4, "K": 7, "P": 75, "r": 10.7143, "p1_steps": 1114112, "p2_steps": 3735552, "p3_steps": 82944, "comps": 82, "turns": 875, "max_comps": 87}
{"id": 6, "score": 1687.0, "time": 9839, "N": 19, "C": 3, "K": 8, "P": 47, "r": 5.875, "p1_steps": 900096, "p2_steps": 4931584, "p3_steps": 194560, "comps": 40, "turns": 193, "max_comps": 40}
{"id": 7, "score": 4257.0, "time": 9829, "N": 22, "C": 4, "K": 3, "P": 33, "r": 11, "p1_steps": 3, "p2_steps": 10829824, "p3_steps": 156672, "comps": 139, "turns": 330, "max_comps": 140}
{"id": 8, "score": 1513.0, "time": 9830, "N": 19, "C": 6, "K": 5, "P": 29, "r": 5.8, "p1_steps": 1425408, "p2_steps": 6008832, "p3_steps": 182272, "comps": 64, "turns": 343, "max_comps": 66}
{"id": 9, "score": 2676.0, "time": 9826, "N": 25, "C": 2, "K": 5, "P": 36, "r": 7.2, "p1_steps": 9, "p2_steps": 8031232, "p3_steps": 159744, "comps": 80, "turns": 204, "max_comps": 92}
{"id": 10, "score": 997.0, "time": 9827, "N": 21, "C": 2, "K": 6, "P": 24, "r": 4, "p1_steps": 9, "p2_steps": 11063296, "p3_steps": 306176, "comps": 48, "turns": 155, "max_comps": 54}
```

Those are scores from the local run, but they should be fairly comparable. `pX_steps` means the number of steps performed for phase X.


