Links: [contest](https://atcoder.jp/contests/ahc012) & [code](https://github.com/FakePsyho/cpcontests/blob/master/atcoder/ahc012/main.cpp)
## Approach
* Allow only horizontal & vertical lines; lines are slightly skewed to allow making a line between consecutive points; the lines don't through the points
* Use all 100 lines: around 8-10 vertical - remaining are horizontal; the split between horizontal & vertical is fixed for each test; I do allow lines in the same spot as it allows to effectively "delete" a line if that's required
* Basic solution is SA:
	* State: set of horizontal & vertical lines
	* Transition: Take a random line and move it: 35% to a random spot, 65% to within [-13,+13] relative position
	* Score: `∑ min(a_d, b_d) * pow(d / 10.0, time_elapsed * 2.75)` where `time_elapsed` starts at 0 and goes to 1; this increases the weight of groups of 10 at the start of SA and then slowly converges to the scoring from the problem statement; was crucial when going from 98M to 99M.
	* Temp schedule: standard exponential with heavily-tweaked constants: `0.8 * pow(0.05 / 0.8, time_elapsed)`
* Around 500k evals within time limit.

Overall that's a very straightforward solution except for the trick with the scoring function. Without that trick, I usually end up with way too many of regions of size 1, sizes 2-9 mostly correctly and very few of size 10.

I though about adding an additional phase after the main SA where I would allow non-horizontal lines in order to make very small adjustments. This should drastically improve the score, but unfortunately it's not easy to implement.

## Random Comments
* As usual, I started the contest late (this time only 30 mins) & extremely sleepy; wouldn't mind if AHC started slightly later ;)
* I saved the screencast (3,5h long boring video without any commentary) from the contest. If I ever upload it, I will link it here.
* I tried to setup a VM, but I failed with installing rust on linux, so I gave up. With external VM this would be a completely different experience for me considering I was mostly tweaking constants during the last 2 hours.
* Overall, really cool contest and a good fit for 4 hours. 
