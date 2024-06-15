## Links

[Contest](https://codeweekend.dev/#/)

[Problem Statement (pdf)](https://drive.google.com/file/d/1ogJtACkj3MIGnIE9Ar3eeEW_EFPil3yX/view)

[Code](https://github.com/FakePsyho/cpcontests/tree/master/various/codeweekend_1/main.cpp) - My final code; note that the way I run my code is a bit more complicated (it's explained in the last two sections)


## Approach Summary

I developed a single solution for all 50 given tests using Simulated Annealing (SA), where the state represents the order of the monsters. Evaluation is done via a complete simulation (killing all monsters in the given order) with the score calculated as per the contest's exact scoring system (i.e., gold after `num_turns`). I used 3 transitions with equal probability: (a) move monster (uniform) to another position (uniform) (b) take a random monster (uniform) take another monster (uniform from the 10 closest monsters) and move it next to the first monster (c) take two monsters (uniform) and swap their positions.

During the simulation, the hero greedily moves towards the next monster. I calculate the expected position if the hero were to move in a straight line, then check all points in a 3x3 grid centered around that position, selecting the closest one within the speed range. If this move puts the hero in attacking range, I further optimize by pushing the destination towards the next monster. Generally, if the hero's destination is within range of the next X monsters, I adjust the destination to be closer to the X+1 monster.

For Day 2, the only big difference was adding the new rules to the simulation.


## Timeline

I expected that I won't have too much time available throughout the weekend and didn’t prepare any tools beforehand (which I realized was a huge mistake). The timeline is a very rough estimate, it's only to give you an idea of how I approached the contest. For the record, contest started at 11pm of my local time.

- **Day1 Sprint**. I read the problem statement, briefly checked the tests to understand the style and constraints, and started implementing SA. While it was clear that SA would perform well, it probably wasn’t the best method for any of the tests. I wasted some time manually submitting before realizing the provided API notebook implemented all necessary functions. Overall, I felt that I was doing well (in terms of ideas & implementation), but somehow my scores were significantly far from the best. If I remember correctly, it was nearly 10%, which is a lot for this particular type of a problem. Four hours in, I noticed I wasn't updating the hero's range after leveling up. Fixing this and rerunning the tests put me in the lead... which I lost 2 minutes later. I was already sleepy, so I decided to run some tests overnight on my laptop. I was too exhausted to do it properly, so I just wrote a single batch file with multiple runs and went to sleep. When I was going to sleep, I was in the 2nd place.

- **Day1 Morning**. I was surprised to wake up with being in the lead again. Based on [qwerty787788's video](https://www.youtube.com/watch?v=Vs7qYHiajSo) my & nika's solution were exchanging places throughout the whole night. I had a short coding session where I implemented a simple tool that would allow me to schedule multiple tasks (runs) and split them onto different cores. I thought that I'm still doing well despite running everything overnight on a single core so I might use all 6 cores of my laptop.

- **Day2 T-1 hour**. Not sure why I did it so late, but at that time I finally made decision to add a small 16-core VM machine. My code worked 3x as fast so I had an 8x improvement over my local setup. It was barely enough to snipe the day1 lead with a tiny margin. If I implemented an intermediate output in my code, I'd probably have 20-30 point lead, as a lot of runs finished not long after the deadline.

- **Day2 Sprint**. Seeing the new changes, I thought my solution wouldn’t perform too well. After some deliberation, I decided to just add the new rules to the simulation and found that SA solved most of the new test cases really well. Within an hour, I had a working solution except for two hard tests (36 & 37), for which I hardcoded a starting move (go straight down) to gain some points. I scheduled overnight tests (with gradually increasing `SA_STEPS`) on the 16-core machine, and went to sleep around 1 AM.

- **Day2 Morning**. When I woke up I was still in the first spot, but my lead was melting with every minute. In every test where I scored far below 1K, hero was not able to evade all of the monsters at the start of its path. I thought about potential solutions and figured out that probably the best way is to add any reasonable beam search (BS) that would solve the first 10-15% of the test and then solve the remaining part via already-working SA, considering it works surprisingly well. The only problem was that it felt like a lot of work and I didn't have enough time or energy to do it. At the time, wata got a huge lead, so it was clear that I won't be able to regain top spot with a tiny change. I decided to move on with my day with hope that I'll figure out a simpler solution.

- **T-5 hours**. Several hours later, I still haven't found an easy way to fix the bad start for harder tests. I could probably go with randomized greedy instead of BS, but that's about it. First place was far away (more or less equal to how much I was losing on 36, 37 & 39), but second place looked possible without solving the issue. I decided to add a 96-core machine for the last few hours. I thought that it will easily pay for itself if I rerun my code on some day1 tests where I lost the best. I added two additional features `AVOID` (when adjusting destination avoid running into more monsters) & `XPATHFIND` (if the destination is under monster's attack, rerun pathfinding step analyzing all potential movement destinations; extremely slow O(speed^2) and honestly no idea if it actually works) with very minimal impact. Other than that, I haven't done anything useful. Around 90 mins before the deadline, I've noticed that during simulation, the hero always goes as far as possible towards the monster disregarding monster attack range. After a quick solid facepalm, I fixed it by adding the `XSTOP` feature. It did increase some of my scores, but the impact was much smaller than I anticipated. My guess is that with enough iterations, SA was always able to overcome this. So effectively not having `XSTOP`, made SA "slower" by up to (I guess) 5x. The remaining of the contest was spent refreshing the results page.


## Tools

* **Runner**: A small script allowing me to compile and run my code with adjusted parameters via a single command. For example, to run test 30 with `AVOID=1` and `SA_STEPS=1e7`, I use `r 30 @AVOID=1 @SA_STEPS=1e7`. For every argument starting with @, it scans the source code and replaces the original declaration with the given one. This script enables easy testing of different solution versions without editing the source code or creating complex command-line argument parsing. I also use this in other heuristic contests and if it ever matures, I might open-source it so it joins the [psytester](https://github.com/FakePsyho/psytester) & [psyleague](https://github.com/FakePsyho/psyleague) family.

* **Task Scheduler**: During the contest, I created a basic tool to run predefined tasks on multiple cores. Post-contest, I learned this could be easily simulated via the [parallel command](https://www.gnu.org/software/parallel/). I wish I had a more complex tool to edit tasks on the fly, view outputs, save run history, etc. Managing runs during contests with downloadable tests is a huge timesink for me, and a GUI would be highly useful. If Codeweekend occurs more frequently, it might be worthwhile to create such a tool.


## Random Thoughts

* The contest was very well organized. It was really nice to have an API, everything was on time, and the website was very responsive. I don't think there were any hiccups during the 48 hours. The problem was nice, though a little bit basic. My guess is that they wanted to make something more approachable to people with less exposure to heuristic contests. 

* I didn’t do many experiments to test the effectiveness of my temperature schedule, initial state, or state transitions. With a limited number of tests, it's challenging to draw conclusive results. Most inefficiencies can be mitigated by running the solution longer, so having a "safe" method that eventually produces good results is crucial. I'm mentioning this, so it's clear that the choices made were more of an educated guess.

* As mentioned in the Tools section, I run my code through the runner. Which means, all runs had some parameters changed. Unfortunately, I didn't save the run history. SA_STEPS directly impacts runtime and result quality. `MAX_DEPTH` was usually omitted (set to maximum). For tests 26-50, I experimented with `MAX_DEPTH` and `AVOID`, but I don't know what produced the best results. `XPATHFIND` was tried only on 36 & 37 with mixed results. Some of my longest runs lasted up to 5 hours, but most of my bests came from runs between 15 minutes and 2 hours.