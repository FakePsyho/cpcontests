## Links
[Problem Statement](https://www.topcoder.com/challenges/30322331)

[Final Solution](https://github.com/FakePsyho/cpcontests/blob/master/topcoder/mm142/HelpTheElves.cpp) - Final submission (bad code with a lot of redundancy but easy to follow)

## Intro
The concept was really cool (especially from the design perspective), but as others have already mentioned the problem turned out to be way too complex. It also gave massive advantage to java coders, considering you had to translate AI solutions to your own language. AI #3 used rng which made it even worse. It would have been nice to get at least the AI's id as a parameter.

Even though I really like problems like this, the scope here is so large that I had trouble finding motivation to write a decent solution. In the end, I did most of the work on the very last day. I have implemented AI #2-#5 in my code, but I didn't have time to do any heavy exploitation, so I'm only using this to detect the opponent. I have a single solution that handles all of the opponents, but (as usual?) I abused testing (5K tests per 5 mins) for some parameter tuning.

My code is quite buggy and my elves are performing a ton of useless operations. I didn't have time to implement explicit exploitation of #1 (precompute boxes' placement to put elves on the path) and #4 (which IMHO was very easy to get nearly-optimal number of presents, but much harder to exploit the maximum number of wasted money).

## Solution
Definitions:
- reachable: means that the cell is reachable from the border
- alive: means that the cell is reachable from any of the presents
- guarding: special property of the elf with a box, which generally means that the elf will stand still and guard it's position
- sidestep: special move for guarding elves to make space so that other elves can pass; after performing the move, it tries to go back to the original guarding position
- basic elf: elf that doesn't carry anything

Order of each turn:
1. Update opponent AI prediction by running the code on the previous state
2. For all guarding elves: remove "guarding" property if the cell is not alive or all neighbors are reachable
3. For all guarding elves: if there's an elf with a present on a non-reachable cell, make a sidestep to make space for that elf
4. For all elves with presents: go along the fastest path to the border
5. For all guarding elves: if there's an elf in a non-alive cell, make a sidestep to make space for that elf
6. For all non-guarding elves with boxes: go along the fastest path to the border
7. For all basic elves that can reach any present: go along the fastest path to the present
8. For all basic elves that can reach any guarding elf that has at least 2 empty cells around: go to the closest one
9. For all basic elves that can reach an alive box: choose the one that has the most empty cells around it; in case of a tie, take the closest one
10. For all basic elves that can reach any box: go to the closest box
11. Perform all of the planned sidesteps or go to the border if unsuccessful after X turns

In each step, I'm only considering elves that have not moved yet. Except for the very few places, the logic is exactly the same for all opponents.

Apart from the above, there's also a special "waiting" logic that:
- for AI #2 & #3: in tests with smaller `elfP * C`, disallow all moves for the first 75-80% turns where we would pick up a box
- for AI #4: in tests with smaller `elfP * C`, don't perform any moves for the 70% turns.

Waiting is the main logic used for exploiting AIs #2-4, since they are allowed to place only a single box per turn. It was counterproductive for AI #1 as the guarding mechanic was good enough. 

I also have a cheap trick of temporary marking some close-to-the-border presents as trees so that the elves prefer the presents in the center at the start of the simulation.

## Example Results
```
{"id": 1, "AI": 1, "score": 412.0, "time": 20, "N": 10, "C": 1, "elfP": 0.1, "M": 36}
{"id": 2, "AI": 5, "score": 18750.0, "time": 354, "N": 30, "C": 10, "elfP": 0.2, "M": 9}
{"id": 3, "AI": 4, "score": 8620.0, "time": 256, "N": 29, "C": 7, "elfP": 0.111968, "M": 14}
{"id": 4, "AI": 4, "score": 221.0, "time": 14, "N": 10, "C": 3, "elfP": 0.117386, "M": 12}
{"id": 5, "AI": 1, "score": 5336.0, "time": 283, "N": 27, "C": 8, "elfP": 0.10486, "M": 8}
{"id": 6, "AI": 3, "score": 376.0, "time": 36, "N": 13, "C": 2, "elfP": 0.161918, "M": 14}
{"id": 7, "AI": 1, "score": 8053.0, "time": 302, "N": 28, "C": 3, "elfP": 0.130583, "M": 3}
{"id": 8, "AI": 4, "score": 845.0, "time": 31, "N": 14, "C": 5, "elfP": 0.150601, "M": 5}
{"id": 9, "AI": 3, "score": 5138.0, "time": 131, "N": 20, "C": 6, "elfP": 0.1832, "M": 42}
{"id": 10, "AI": 5, "score": 1911.0, "time": 22, "N": 14, "C": 1, "elfP": 0.10214, "M": 10}
```

My guess is that I do better on #3 than others, since I don't believe there is any clear exploitative strategy. I'd imagine the final standings can wildly differ from the provisional if the test distribution was highly skewed.

Out of curiosity, I made a run where I replaced all of the AIs with #5 to figure out how many points I'm missing vs each AI:
```
Tests                  5000     1021     1018     1008      981      972
Run                 Overall     AI=1     AI=2     AI=3     AI=4     AI=5
------------------  -------  -------  -------  -------  -------  -------
all_ai5             98.5909  98.4543  99.0673  98.8439  98.6716  97.8915
final_submit        85.3864  88.9982  77.8269  73.3607  89.5334  98.0376
```

