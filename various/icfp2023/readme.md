## Team: is_topcoder_dead

eulerscheZahl

Illedan

nika

Psyho

sullyper

## Intro

Our team had a rather laid-back approach to the competition. We have set up a discord before the contest, but our cooperation was mostly limited to exchanging ideas and (later) sharing output files. No shared codebase etc. Hence, this post mortem is going to focus on my personal experience. 

This was the first time I competed in the ICFP challenge. I looked at it, many years before, and both of those were closer to puzzle hunts. This year, the challenge was a rather straightforward optimization problem. This was very close to your standard Marathon/AHC except the tests were given to you (so you had to produce results on your own) & there were few changes introduced during the contest: more tests & additional rules.

## Approach

If you have any experience with the optimization problems, your inner voice was screaming "local search" during reading [the specification](https://github.com/FakePsyho/cpcontests/tree/master/various/icfp2023/specification.pdf). So obviously I went with SA (Simulated Annealing). Even if SA turns out to not be the best approach, it should be able to generate near-optimal solutions to small tests and you can use this knowledge to find a more structured approach later.

**Initialization:** Random state

**Transitions:**

* (A) Swap two musicians

* (B) Move a random musician to a different spot. Moves are accepted only if there are no collission with other musicians. All options have equal probability. 
    * (B0) Random place
    * (B1) Random place within 50x50 square centered on the musician
    * (B2) By a random vector of length 25<sup>[0..1]</sup> (uniform probability)
    * (B3) Random place touching the edge of the stage
    * (B4) Random place 10 units away from the edge of the stage
    * (B5) Choose another random musician and then choose a random place that touches this musician
    * (B6) Choose a random direction and move towards it until we hit another musician or edge of the stage

**Eval:** Same as in specification

**Complexity:**

Naive implementation is **O(M<sup>2</sup>A)** per eval. You can make it dynamic to reduce the complexity to **O(M+A)** for swap transitions (A) & **O(MA)** for move transitions (B). 

&nbsp;

After implementing initial SA (basically everything that was mentioned above except B6), my next goal was to increase the number of iterations in SA. I had two main ideas: reducing complexity to **O((M+A)log(M+A))** by doing some sort of radial sweep or to add multithreading (MT). My end goal was to do both, but I thought that adding MT first is going to reduce the amount of work, since I'll have to do a major refactor for this, so the sooner the better.

I already had generic SA code that does MT (originally developed for Reply/Hashcode contests), but it's a huge mess that's optimized for speed. It is very much untested and has essentially no debug messages if anything goes wrong. I have wasted a lot of time tracking esoteric issues originating from flawed design. For example, swap moves are much faster to execute than transition moves. Because of that, transition moves were never executed when acceptance rate was high. The simple fix was to sync transition types in all workers.

In the end, I spent so much time fighting those weird bugs and micromanaging my SA runs, that I haven't improved the time complexity of evals.

## Random Thoughts:

* I'm not sure which transitions are helpful. B3, B4 & B5 speed up initial convergence, but it's quite possible that they increase the chance of getting stuck in a local optimum and they make it much harder to achieve global one. It's nearly impossible to test those things out when your runs take several hours to complete in a contest that lasts 3 days.

* Since swaps were much faster to evaluate, swaps are executed more frequently. It definitely speeds up initial convergence but I'm not sure if that was a good idea in the long run.

* Originally, I thought about allowing overlapping musicians in my SA and this decision made my solution a bit more complicated than it should. Penalty for collision would start at 0 and gradually grow to infinity. This generally helps with achieving dense packing and smoothes out gradients when doing transitions. Unfortunately, this doesn't interact well with the sound blocking mechanic and potentially increases the number of required iterations in order to achieve a semi-decent solution. I thought about revisiting this idea when my solution is very fast and I still have some spare time.

* Towards the end, I thought about alternative method of speeding up the evals: Partition whole stage into small rectangular areas and then precompute potential blockings for each pair of partitions. You can do the same thing with (partition, pillar) pairs. While it's not the fastest pruning method, it's very easy to implement (despite that, I haven't done it).

* Due to how quickly the team "All your sound are belong to us" had a very good score, I had wrongly assumed that there's a non-SA solution that gives very good scores quickly. In many contests, you can extract a lot of info from the leaderboard, especially if you know the strenghts and the weaknesses of other people. I had no idea who hides behind that team name, other than it's probably a non-Japanese team. All being said, this assumption completely backfired for us, because we were all heavily focused on exploring non-SA approaches. IIRC, there were only 2 tests out of 90 that weren't solved entirely by the SA in our team.

* I have used 16 core machine during majority of the contest and then 96 core machine for the last several hours.

* I recommend reading [RafBill's approach](https://gitlab.com/rafaelbocquet-cpcontests/icfpc23) (the aforementioned "All your sound are belong to us" team) that most probably is going to win the contest with a huge lead. The SA explained there is very similar. I believe the main difference comes the number of iterations.

* Our final result was 132,227,609,016. [scores.txt](https://github.com/FakePsyho/cpcontests/tree/master/various/icfp2023/scores.txt) file contains scores for each test case.







