Final Ranking: #33 in Legend


### Links

[Problem Statement](https://www.codingame.com/ide/puzzle/seabed-security)

[Final Results](https://www.codingame.com/contests/fall-challenge-2023/leaderboard/global)

[Other Approaches](https://www.codingame.com/forum/t/fall-challenge-2023-feedbacks-strategies/202891)

Note: Code is not included as CodinGame prohibits its posting


### Intro

My plan was to start with something very simple and add some basic logic to get a feel for the game. After that, I aimed to implement full tracking and focus on the strategy.

In reality, I only implemented tracking for fish on the very last day.  I also somehow missed that you receive scans for enemy drones (I guess I simply forgot about it after the holiday break), so my strategy almost completely ignores the enemy drones


### Approach

* Basic tracking using rectangles: expand them every turn & use data from the radar to "shrink" them
* Movement:
    - fixed moves for the first few turns
    - if we're within the first 20 turns and we can scare the fish out of bounds (+ some sketchy logic to figure out if the enemy drone will scan it next turn)
    - if we have 9 fish and we're not further from the surface than our opponent -> go to surface with both drones
    - if we have scanned all the fish and a drone have some unique scans -> go to surface with that drone
    - choose the closest unscanned fish (based on the center of rectangles) and go to it in a straight line
    - if all fish are scanned and our drone has no unique scans, try scaring a bottom fish out of bounds
    - if the move would collide with the monster (including depth 2 search), make depth 2 brute force and find the move that doesn't collide + minimizes distance to our target fish (or surface)
* Light:
    - if battery >= 5 and time was not used in the last turn, use it
    - if there are no potential fish nearby (based on distance to the closest rectangle), force not using light
    - if there is a rectangle (fish) that is very close (I'm using distance to the furthest corner) force using light

Overall, a really lazy solution that shouldn't work too well.


### Some thoughts

* I found the game quite interesting: low complexity but has a lot of non-trivial strategy. It has way more depth than it initially seems, which is always a plus in a longer contest. Unfortunately, significant randomness means that it requires more games for accurate ranking calculation than CodinGame offers. At one point, my submits could end up anywhere between 15th and 100th place.

* While I didn't spend much time on the contest itself, I spent significant time on developing [psyleague](https://github.com/FakePsyho/psyleague). I'm baffled that most people don't perform any local tests: they are almost always more reliable than the leaderboard and allow for better time management (you don't have to spam submits and stare at the noisy results). Even if you plan to spend only a few hours on the contest, it's going to save you some time.

* I think I tested around 150-200 different versions locally. Most of the local runs are simple constant changes to see which version performs slightly better. Without that, I doubt I would be anywhere close to top 100.