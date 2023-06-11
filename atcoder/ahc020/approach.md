## Links
[Problem Statement](https://atcoder.jp/contests/ahc020)

[Final Solution](https://github.com/FakePsyho/cpcontests/blob/master/atcoder/ahc020/main.cpp) - Final submission minus dead code.

## Approach
There are two main components that solve each of the subproblems:

1) In order to find output power for each node, I'm using a greedy algorithm that in each step adds a person that's not yet covered and minimizes the additional increase in score. There's additional randomized penalty for using each node. This is called `kgreedy` in my code and it also allows to start with predefined power output for each node.

2) In order to choose edges, I'm just running a simple MST (Minimum Spanning Tree) algorithm. The implementation is straightforward and gets the list of nodes that we have to connect. Unsurprisingly, this is called `mst` in my code.

Global solution is a HC (Hill Climbing) that combines both of those components:
- randomly reduce power output from one of the nodes and recalculate required power output using (1)
- in order to find edges, run (2) for the complete set of nodes and then greedily remove/add nodes that don't need power; my guess is that this produces nearly-optimal result, but unfortunately you have to run (2) hundreds of times per each step and suddenly it's the main bottleneck
- compute new score and keep the new solution if the score is improved

The whole process is very slow and overall I get only around 200-800 steps which is enough for HC to get stuck in a local minimum, but I felt it's not enough to switch HC to SA.

## Random Comments
* As it's my tradition, I started the contest late (~40 mins) and very sleepy. Would be nice if those 4h AHCs didn't always start at the same time.
* I toyed around with my solution after the contest for a bit. Switching HC to SA and adding 10x bigger time limit (which should be doable considering how inefficient my code is) gives around 502M score locally. 
