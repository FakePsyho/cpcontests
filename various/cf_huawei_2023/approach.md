### Links
[Problem Statement](https://codeforces.com/contest/1885/problem/A)

[Final Solution](https://github.com/FakePsyho/cpcontests/blob/master/various/cf_huawei_2023/main.cpp)

[Provisional Results](https://github.com/FakePsyho/cpcontests/blob/master/various/cf_huawei_2023/results.txt) (it's around 15 points lower than my best)

[Stats from provisional Tests](https://github.com/FakePsyho/cpcontests/blob/master/various/cf_huawei_2023/provisional.csv) (data scraped via memory usage)

### Intro

This was probably the most interesting Huawei contest I've done so far. The problem statement was also much better and easier to understand than the previous ones. 

Sadly, most of the complexity was hidden within the tests again and even though we got few tests to download, it would still be highly beneficial if we could get at least 20-30 more next time. Another thing that would be nice to have is implementation of scoring function. It's silly that the optimization contest with the massive prizes has worse standard than the average contest without any prizes.


### Approach

Conceptually my approach is very simple:
- Implemented all of the formulas (easier said than done)
- Defined state as power usage (exactly as the output)
- Implemented Simulated Annealing (SA) on top of it

Eval function for each frame is 1.0 if it's finished and `0.1 * finished_percentage ^ 3.5` otherwise. Complete eval function is just a sum of all frames.

I experimented with different transition functions, but in the end I was left with:
- Change power of a random `(t,k,r,n)`
- Move power from `(t,k,r,n)` to `(t,k,r,n')`
- Move power from `(t,k,r,n)` to `(t,k,r',n)`
- Move power from `(t,k,r,n)` to `(t,k',r,n)`
- Swap entire `(t,k,r,...)` with `(t,k',r,...)`
- Swap entire `(t,...,r,...)` with `(t,...,r',...)`

The eval code is very optimized (it's mostly `State::set_power` and `State::calc_user_data`), but sadly the solution doesn't get enough steps to obtain a good score in bigger / more complex tests. Reducing time by 50%, reduces my score by roughty ~150 points. Overall, this was enough for ~17k.

For tests with high d (2-30), I start with zero power. For tests with high d (31-51), I try to evenly split the power among all used users for that particular `t`, in such a way that for every `r` only a single `k` is used. This gave around 200 points, mostly from two big test cases (33 & 38).

There's a lot of tiny tweaks and tricks in my code in order to squeeze maximum value out of each SA step. If you're going to read the code, it's probably best to ignore all of the weird continues in my transitions.


### Additional Comments

- Most probably there will be a significant shakeup between provisional and final standings. The provisional test cases can be split in 4 different groups: one split based on average d (high vs low) and one split based on len=1 and len>1. If the final tests have a different distribution, the final standings might be wildly different from the current one (without taking extreme overfitting and failed tests into account).

- My solution needs way more time (steps) to be competitive. My plan was to develop decent SA solution, then run it for very long time and see what the near-optimal solution looks like in order to develop further optimizations/simplifications of the problem. Sadly, I never got to that point.

- Alternative approach, would be to add a greedy initialization (which I believe everyone else did at the top), and use SA only in order to squeeze the last few points. The reason why SA is so useful here, is that the formula feels very "random" and it's hard to find the optimum without performing any random perturbations

- Worth mentioning is that I don't have any transitions between different `t`. They were way more costly than others + I thought I'm mostly losing points on len=1 tests, so I didn't explore those.

- I really hope that for the next Huawei contest, we could get a significant number of offline tests. **It's really hard to develop a good solution, where most of the complexity is hidden within the tests itself.** Submitting hundreds of solutions is a huge waste of time and resources. To be clear, it's great that we at least got those 3, but it would be much better to get more.