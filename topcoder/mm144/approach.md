## Links

[Problem Statement](https://www.topcoder.com/challenges/30337643)

[Final Solution](https://github.com/FakePsyho/cpcontests/blob/master/topcoder/mm144/Flood.cpp) - Final submission. Complete mess. No educational value.


## Intro

I remember the web puzzle game that this problem is inspired by (or at least I'm pretty it is). Despite having some prior knowledge, I wasn't confident that there's any simple approach that would work decently. This, coupled with the TCO news definitely killed my motivation. I'd imagine I'm not the only one considering this might have been the record low participation rate, despite having an interesting and well-thought problem. In the end, I decided to fight for the 3rd place and write some quick & ugly solution.


## Approach

Core part is a full simulation where each builder gets an ordered list of walls to build and tries to build them one by one. I'm using A* (with Manhattan distance) instead of the bfs but otherwise it's a straightforward without any tricks.

There two phases: 

Phase One (P1) is a global SA where in each step I generate a random list of walls for either one or all of the builders (each builder gets a separate and independent list). There are following "paths" possible:
- Build 4 walls around a single cell (same as example solution)
- Build 8 walls around a single cell (essentially forming a + sign)
- Choose a starting point and a random direction (including diagonals). Build walls until the edge of the map is reached.
- Choose a starting point and a random direction (including diagonals). Build 1-10 walls, choose a perpendicular direction and the continue until the edge of the map is reached.
- Same as above but change direction twice.

Phase Two (P2) is a local HC where in each step we perform one of the following transitions:
- Choose a builder and remove a random wall
- Choose a builder and move a random wall in a random direction
- Choose a builder and add a random wall anywhere on the map
- Choose a builder and swap the order of the random walls

P2 is 7 seconds. P2 is 2.5 seconds and uses the solution produced by P1.

There are some minor tweaks, but nothing worth mentioning. I tried other transitions for P2, but I ran into the problem where I make random changes that would produce the same score, but would require more turns. 

I'm honestly surprised that this approach is able to produce reasonable results for a lot of testcases. Especially considering it's very unoptimized. Giving 2x time gives ~2% improvement, which may mean that this is a viable approach somehow (at least for smaller tests) with a little bit more polish.

## Example Results
```
{"id": 1, "score": 47.0, "time": 9518, "H": 8, "W": 8, "HW": 64, "T": 2, "B": 2, "S": 5}
{"id": 2, "score": 336.0, "time": 9524, "H": 30, "W": 30, "HW": 900, "T": 5, "B": 5, "S": 14}
{"id": 3, "score": 283.0, "time": 9519, "H": 27, "W": 29, "HW": 783, "T": 3, "B": 4, "S": 7}
{"id": 4, "score": 79.0, "time": 9523, "H": 29, "W": 22, "HW": 638, "T": 5, "B": 3, "S": 7}
{"id": 5, "score": 569.0, "time": 9523, "H": 28, "W": 23, "HW": 644, "T": 2, "B": 3, "S": 22}
{"id": 6, "score": 49.0, "time": 9521, "H": 19, "W": 13, "HW": 247, "T": 3, "B": 1, "S": 15}
{"id": 7, "score": 152.0, "time": 9526, "H": 22, "W": 22, "HW": 484, "T": 2, "B": 2, "S": 17}
{"id": 8, "score": 83.0, "time": 9529, "H": 19, "W": 13, "HW": 247, "T": 3, "B": 4, "S": 7}
{"id": 9, "score": 313.0, "time": 9521, "H": 25, "W": 14, "HW": 350, "T": 2, "B": 3, "S": 10}
{"id": 10, "score": 325.0, "time": 9518, "H": 21, "W": 16, "HW": 336, "T": 1, "B": 2, "S": 9}
```

