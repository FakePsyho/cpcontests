## Code

[Final Solution](https://github.com/FakePsyho/cpcontests/blob/master/topcoder/mm137/StoneSeries.cpp) - My final submission minus dead code plus improved readability. Surprisingly clean except for the `greedy()`/`notgreedy()`functions.


## The Problem

It's generally a good idea to understand how optimal/nearly-optimal solutions should look like. I spent some time toying around with the problem and convinced myself that snake-like solution should be very close to optimal. Your achievable score directly depends on the highest stone value, so you want to maximize the length of the snake. The easiest way of achieving that is having a snake (path) going through cells that only have a single 1 as a neighbor. 
```
..1oooooooo.
....1..1..1o
.oooooooooo.
o1..........
.oooooooooo.
....1..1..1o
.oooooooooo.
```
Zigzag pattern pictured above works quite well, spirals should be just as good if not easier to execute.

The snake can be done in a single segment, but it's also fine to do it in multiple segments (probably slightly less efficient, but might be the best option in tests with high % of blocked cells). You also want to move horizontally/vertically since you want to minimize the number of neighbors you have on the way.

Having a good idea of optimal solution helps with evaluating approaches, because if your approach is not able to converge to a nearly optimal solution, it's usually not a viable approach. 


## Base Solution

After a few initial failed approaches, I tried greedy (at each moment try to place a highest-valued stone) and it seemed it was producing quite good results. The two natural ways of extending it would be (a) beam search (instead of generating whole solution in one go, try to generate small chunks and keep a family of multiple solutions) and (b) simulated annealing (iteratively remove and add stones in various places). The problem with SA is that it's not great because moving stone X destroys all of the stones that are dependent on it. In particular moving 1s is very problematic. This leaves beam search as an obvious candidate. Unfortunately for me, SA is my comfort zone so I went with SA.

The vanilla SA is as following:
* State: actual stone placement
* Evaluation: actual score
* Transition: randomly select few stones and remove all of the stones that are dependent on them; after that, run the already mentioned greedy method to place new stones

With carefully implemented greedy method (you need a fast way of finding the next stone to build) and some very crude parameter tuning, this should get something around 65-70 points on the leaderboard.

This may look like a decent start, but SA+greedy won't be able to achieve long snakes: You can't rearrange 1s without destroying half of your solution. At the same time, placement of 1s heavily impacts the maximum achievable score on a particular test case.


## Improvements

The provided impact is calculated by comparing my final solution and final solution without that particular feature. This means that impact of those features don't necessarily add up, but hopefully such information is still useful. The rerun happened only on 100 tests (compared to 2K tests I have used near the end of the contest), so those values can be quite off, but probably within 0.5.

* Various tricks for dealing with "immovable 1s" problem:
	* Instead of starting with a random solution, generate multiple greedy solutions and take the best one. The intuition here is that there's a positive correlation between the solution quality from  greedy and the potential (i.e. maximum achievable score) for a particular set of 1s. After initializing with greedy, most of the non-1s solution is going to be destroyed anyway. Impact: +2.5.
	* The is the trick that made this whole approach viable: Based on my analysis of how optimal solutions should look like, it's clear that the easiest way of achieving high scores, is to have 1s placed in such a way, that you maximize the number of empty cells that neighbor with a single 1. This makes building long snake segments easy. Instead of running greedy on an empty grid, I generate a solution that consists of only 1 stones that maximize this property (see `build_ones()`). I do this via HC (transition: place a new 1 and remove "colliding" 1s with the newly places stone. Collision means that there is an empty cell with 2+ neighbors). Impact: +10 (on top of greedy). Since this process is costly, precomputing those configurations and reusing them, gives another +1.
	* Multiple hard restarts for smaller tests. Considering greedy creates a new solution in a single pass, it's simply cheaper to generate a new solution from scratch than make a very aggressive temperature schedule. Impact: +1.
	* Temperature schedule that never really stabilize. This is a very common thing configuration when you combine SA with a very complex transition using some greedy approach. At some point, my temperature schedule for the largest tests was actually increasing, which means I had SA without the "annealing" part. Converting SA to HC reduces the score by 5 points.
	* My SA+greedy has a tendency to rebuild 1s in bad spots. This is especially true in the beginning when the whole solution is very unoptimized. As a quick fix, I disallow removing 1s during the first 10% of the run. After few initial moves, the cost of removing 1 and placing it in a worse spot is much higher, so SA is less likely to accept such moves. Impact: +0.5.
* Greedy is very efficient, but in many cases, it's nearly deterministic so it's very easy to end up being stuck in a local minima. In order to remedy this, 25% of the time, I use a variation of greedy method (see `notgreedy()`) with a different priorities for stones. Instead of (highest value -> lowest value) it uses (highest value, one, randomized order for everything else). It's still important to prioritize highest value as it helps to build chains that increase the maximum value. Impact: +3.
* My normal routine when removing stones is to select few of them and then remove all of the stones that are dependent on them. The stones are selected at random disregarding their placement. Sometimes you do want remove multiple stones that are close to each other but are independent. To fix this, 30% of the time I remove a big patch (2x2 - 5x5) of stones instead. Impact: +1.


## Things that didn't work

* Tweaking building initial 1s placement. More conservative, less conservative, allowing cells with 2+ neighbors, etc. Considering how impactful this addition was, even a small improvement here should be significant. I also imagined that different tests will require slightly different strategies. Unfortunately, no luck here. 
* Multiple restarts in order to find the best candidate solution and then finish one of them with the remaining time. I thought it's going to be a certain improvement, but it turned out to be false. Probable reason is that for larger tests, I'm never getting enough moves in order to arrive at a good state, while for smaller tests I'm getting stuck very fast so it's easier to just add more full restarts.
* Increasing probability of extending the path sideways instead of diagonally. Increasing this probability only slightly almost doesn't affect the score, so there's probably more subtle way of doing this. This should be particularly useful for big tests where we don't get enough moves to arrive at the local maximum.
* Optimizing for maximum stone first and then later switching to the actual scoring. I used `max_stone * max_stone * 0.5 * (1 - moment) + score * moment` where moment was gradually changing from 0 to 1. I tried different schedules for the moment variable, but all of those produced significantly worse results.


## Example Results
```
Test Case #1:  Score = 296.0     Maximum = 21
Test Case #2:  Score = 254224.0  Maximum = 582
Test Case #3:  Score = 159700.0  Maximum = 475
Test Case #4:  Score = 167950.0  Maximum = 464
Test Case #5:  Score = 141206.0  Maximum = 446
Test Case #6:  Score = 10912.0   Maximum = 130
Test Case #7:  Score = 1431.0    Maximum = 43
Test Case #8:  Score = 5929.0    Maximum = 97
Test Case #9:  Score = 12754.0   Maximum = 142
Test Case #10: Score = 92651.0   Maximum = 372
```

## Random Remarks

Giving my solution 20x bigger time limit (with no other adjustments) increases my score by around 8 points. Mostly for larger tests. I was surprised that my solution was still able to scale quite well with more time, despite having a feeling that it easily gets stuck in a local minima. Immediate increase of time limit by 30% gave around 1 point, so I tried my best to optimize the code. At the same time I thought that all of my constants were already highly optimized and I ran out of simple tweaks to implement. My guess is that my core solution is inferior to other solutions at the top, but I went to extremes with squeezing every last bit of performance out of it.

I believe my solution works fine for tests with big number of blocked cells as the optimal solutions for those tests tend to be very chaotic and the placement of 1s shouldn't affect the scores too much. I'm also pretty sure that I'm losing a ton of points on big tests with few blocked cells (long snakes should triumph here) and very small tests (beam search based solutions should be much better at fully exploring the state space). 

My optimistic plan was to develop an secondary approach that would be very good at building long snakes and then just run my SA+greedy on top of that snake (possibly fixing it and/or adding perturbations). The nice part about building snake is that you don't really care about the actual values of stones as long you're only using cells that have a single 1 as a neighbor. This allows for building the snake via iterative methods (HC/SA). The snake tends to split the whole diagram into many independent diagrams (if you ignore the maximum stone value), which might help with quick calculation of maximum score for each particular snake. Combined with small snake perturbations this should be able to achieve extremely high scores on tests like #2.

I have to say, I really liked the problem. It's always cool to see problems where different archetypes of solutions (SA / beam search / greedy) can coexist with each other and none of them is a clear winner. I do feel that the problem quality skyrocketed throughout the last year or two, so it's a shame that the marathon participation is at the all time low :(

I definitely invested much more time in this marathon than in the last one, although the main reason was probably the fact that it was the very last qualifier for the TCO. Both @Stonefeang and @blackmath had almost a sure spot for stage 4. This meant that a few people (me, @colun, @tanzaku & @JacoCronje ) were fighting for a single spot and each one of us had to win in order to take the last available place at the finals. It felt weird that essentially I need to get 200 points in order to qualify for the TCO. Also, if I understand correctly, there's a small chance that @blackmath might drop to 7th place which takes away his TCO ticket and gives it to whoever placed 2nd. Whatever happens, the system tests are going to be stressful, so I just hope there won't be any close calls.

Oh, and speaking of close calls. I did the shameful thing where I decided to be a bit deceptive and artificially lowered my score by ~4% on the very last day. I believe I had around a 2 point lead going into the last 24 hours. That being said, both @colun and @tanzaku improved greatly, while I couldn't find any improvements in the last few hours. I decided that I'm too exhausted to try the "long snake" approach and instead focused on a bunch of tiny optimizations (mostly code optimizations, hence the large number of submits during the last few hours from me). Unfortunately, almost none of them worked and maybe I was able to get a mere 0.5-1 improvement out of those. Around 50 minutes till the end of the contest I found a silly bug in my stone removal DFS. It was the classic bug of adding multiple copies of the same vertex onto the stack and here the result was that quite often I was incorrectly removing almost the whole grid. Fixing it gave me +2.5 points, which (I believe) moved me from #3 position to the top spot with the comfortable lead. The very last submit just removes the 4% score reduction. Sorry!


