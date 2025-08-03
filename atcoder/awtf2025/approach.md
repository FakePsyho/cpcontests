## Links

[Humans vs AI](https://github.com/FakePsyho/cpcontests/blob/master/atcoder/awtf2025/humansvsai.md) - separate write-up about AI

[Problem Statement](https://atcoder.jp/contests/awtf2025heuristic/tasks/awtf2025heuristic_a)

[Final Submission](https://github.com/FakePsyho/cpcontests/blob/master/atcoder/awtf2025/main.cpp) - final submission with dead code removed



## Intro

This write-up is going to be different from what I usually write. In most cases, I only talk about the final solution and provide data on ablation runs. But since this post is likely to reach a lot of people who don't have much exposure to competitive programming, it's going to focus more on the psychological aspect. 

~~I decided to not include my comments about the AI model to keep this shorter. I might do a separate post about it if there's demand for that (so please let me know).~~


## What is competitive programming

Competitive programming has two main "branches". The most popular one is often just called competitive programming, but usually people refer to it as algorithmic contests. Those contests last for a few hours where you'll get multiple problems to solve and each one of them is something like a puzzle where you have to write a computer program that correctly solves it. Each problem can be either solved or not. Unless the problem is split into subproblems, you won't see any partial scores. This is the more mainstream category, as people often train this solely for the interview process. If you've met anyone doing competitive programming, they probably competed on one of the many platforms (Codeforces, Leetcode, Topcoder🪦). IOI falls into this category as well.

The other category are optimization problems, where you get a single very complex problem and the goal is to produce as good solution as possible. The classic example of an optimization problem is [TSP](https://en.wikipedia.org/wiki/Travelling_salesman_problem). Those contests are closer to a small-scale research and often last for one week or even longer. 

This contest was an optimization contest and was organized by [AtCoder](https://atcoder.jp/), which is currently the largest platform of this kind. They frequently hold online contests and this was the first time they organized an onsite finals: they have invited the 12 best-performing competitors from the 2024.

For the record, AtCoder organizes algorithmic contests as well, but that's not why you're here.


## Exact format of the contest

- individual contest, no communication allowed
- 10 hour long
- access to the internet is fine (except for social media)
- screen sharing is used for detecing collusion 
- we use our own laptops
- we write and submit programs that are remotely executed on the judge's machine (everyone's submission is executed under the same condition: same machine, time & memory limits)
- the provisional ranking (one used during the contest) is a fixed set of only 50 tests; final ranking published during the award ceremony used 2000 tests in order to minimize the variance
- all of the tests are generated according to the specification; we also get scripts for scoring our solutions and generating more tests, which is very useful for local testing
- there's no limit on external compute used for local testing; I believe that most people have rented some external VM for the duration of the contest
- use of genAI is allowed, except automatically generating a large number of solutions: [link](https://atcoder.jp/posts/1495); that being said, I believe no one used genAI except for simple autocomplete


## About me

I'm Psyho and I enjoy problem solving in various forms. I've been doing competitive programming as a hobby for around 20 years now. My main specialization is optimization contests, where I've [won quite a few](https://cphof.org/profile/topcoder:Psyho) major contests. For a few years, I was making a living out of winning machine learning contests. I also enjoy a bunch of random stuff loosely revolving around mental challenges: [competitive logic puzzles](https://en.wikipedia.org/wiki/World_Puzzle_Championship), escape rooms, puzzle hunts, [video game blind racing](https://mysteryfun.house/) and a few others. 

My other passion is game design. I've participated in many game jams, I've built a few escape rooms (in Warsaw). I even had my own indie game studio for a bit, but it turned out that finishing and releasing games requires actual work 😅


## Before the contest

You don't really have to do anything in order to prepare for such contest, except for some [tools for local testing](https://x.com/FakePsyho/status/1945096985792667840) (which is optional). Based on my experience, the only thing that truly matters is to [get enough sleep](https://cdn-ak.f.st-hatena.com/images/fotolife/a/atcoder/20250715/20250715100640.png) (those are my answers to the questionnaire that we filled out a few weeks before the event). I had plenty of contests in the past that I lost solely due to running out of energy at some point. If it's not obvious: contests like this are incredibly demanding in terms of focus. They drain your energy like crazy and when I'm done, I often feel like I just had an intense day full of hiking. Of course it didn't work out. I had a rather busy schedule before Japan and I was already slightly sleep deprived. I couldn't fall asleep on the plane (16 hour flight total) and I slept only [one hour](https://x.com/FakePsyho/status/1945230603688816912) on the night before the contest.

If I had been well rested, I would estimate my chances of winning at 15-20%. I'm definitely way worse at implementation than I was 5+ years ago. It's very rare for me to do any coding outside of the contests now, and I've barely done any contests this year. The thing that works in my favor is that during onsites (unlike online) I really give 100% of what I have and that all of my competitors have far fewer onsites under their belts than I do.

Another thing that helps me is that (at least historically) at onsites most people go for the high-risk high-reward approach: complex solution that if it works it's going to perform very well. The issues are that (a) everyone underperforms at onsites, so you overestimate how quickly you'll implement something and (b) it's hard to focus at onsites and complex solutions often require deep analytical thinking. On the other hand, my goal is to optimize expected position (I just want to be in top 3) and I'm the laziest competitive programmer you can find - I refuse to implement anything until I'm sure that there's no easier way of doing it. Surprisingly it works. 

All that said, I gave myself a 5% due to the sleep deprivation. And I thought that 5% was very generous.


## 10 long hours

From this moment on, I'll assume that you know the [problem statement](https://atcoder.jp/contests/awtf2025heuristic/tasks/awtf2025heuristic_a). Instead of looking at the constraints/test generation, it's enough to view the [web visualizer](https://img.atcoder.jp/awtf2025heuristic/sJKH3KO4.html?lang=en). You can change the seed number to quickly see different tests.

Feel free to look at the [official livestream](https://www.youtube.com/live/TG3ChQH61vE) of the contest. While the livestream is in Japanese and the English<>Japanese automatic translation is, to put it mildly, not good. It's useful for looking at how the solutions actually work. There's also an AI-summary available after each new submission and from what I gathered those are quite good.

The summary of the problem is as follows: There's a 30x30 grid with 10-100 robots, each with start & destination. You can divide robots into disjoint groups (i.e. each robot belongs to a single group). In each time step, you can issue a command in order to move the robots (either whole group or a single one) for a single cell in one of the 4 cardinal directions (up, down, left, right). The robots cannot share the same cells, so this means that they can block each other. The last thing is that you can build walls and some of the tests will have some pre-built walls. The goal is to minimize the number of commands issued to the robots, where each command counts as 1 no matter if you issue it to a single robot or the whole group. There's no cost for placing additional walls. If the robots fail to reach their target destinations, then you get a huge penalty and most probably you'll get nearly no score for the test.

My first impression was that there's a lot of things you have control over: groups, moves & walls. It might be impossible to effectively handle all three things within only 10 hours. In other words, I thought I had to simplify the problem. Usually, the best way to do it, is to look at the tests and figure out how a good solution is going to look like. The main observation I had is as follows: if you have a robot that starts at let's say R2C5 (row 2, column 5) and its target is at R20C24, then you have to issue at least 18 (= 20-2) commands of moving down and 19 commands of moving right (= 24-5). We can generalize this and see that for each direction we have to take the maximum distance between the target and source cells out of all robots. This is our lower bound and in order to achieve it, we have to move all robots as a single group. I've looked at the the few simple tests (low number robots & no walls) and it looked like I could easily find optimal (or near optimal) solutions by hand. For the record, [here](https://github.com/FakePsyho/cpcontests/blob/master/atcoder/awtf2025/notes.jpg) is the only thing that I have scribbled throughout the whole contest. Most of my "notes" were my hand subconsciously drawing something while I thought about the problem. 

So to sum things up: I have figured out that if I move all of the robots as a single group, it should be relatively easy to find very good solutions for small (or maybe even medium) tests. Now I have to find the simplest approach that is able to achieve this. During onsites, I always go for simple solutions, because it allows me to iterate more quickly. If something turns out to not work well, I'll minimize my wasted time. Plus, with my sleep deprived mental state, I'm not able to code anything more complex anyway. I believe this was about 20-30 minutes into the contest. 

The next 15-20 minutes were spent on figuring out, what's the simplest and most comfortable way of using my new observation and I have arrived at something like this:
- always use a single group
- let's use a fixed path for the whole group (easier to implement, maybe good enough)
- the goal is to find the placement of the walls (which is a very hard problem), so that when all robots move along the same path, they finish in positions that are somewhat close to their final position
- as a final step, use [BFS](https://en.wikipedia.org/wiki/Breadth-first_search) for each robot to move them individually to its final destination

This simplified the whole complex problem to just placing walls. I still had no idea how to do it efficiently, but I decided to do (yet again) the simplest possible thing: perform the [hill climbing](https://en.wikipedia.org/wiki/Hill_climbing) (HC) and in each step simulate whole movement and use the sum of manhattan distances as an eval function. Hill climbing here just means that I randomly add/remove walls and I keep the changes that improved the score. I remember I had my first working version around 1:40 into the contest. Results were barely an improvement over the basic BFS.

Over the next few hours I was making constant progress with some incremental changes:
- speed up simulation code (40K -> 120K sims on average)
- switch hill climbing to [simulated annealing](https://en.wikipedia.org/wiki/Simulated_annealing) (SA)
- different fixed paths (for some reason double spiral outperformed single spiral, even though single spiral requires fewer walls)
- experiment with adding/removing individual moves at the start via SA (didn't really work despite my intuition telling me it should be huge)
- in SA, reduce probability of adding a wall if it makes the score worse (in hindsight, that was stupid since better approach is to reduce probability of choosing wall addition as my transition in SA, which is my default approach)
- remove unused walls at the end, right before the BFS so that BFS has an easier job
- improve BFS at the end by performing HC on the BFS order throughout the last ~10% of the run; it drastically helped with harder tests since my robots couldn't get to theirs' destinations quite often

At this point, it was around 6 (out of 10) hours into the contest. While I was hovering around the top 3 human competitors most of the time, I knew I was progressing very slowly. I was constantly fighting urge to sleep (at this point I felt as if I had been awake for 30+ hours), it was hard for me to keep focus and I was introducing bugs I normally wouldn't. Wasted over two hours on a really dumb stuff. And the cherry on top was that OpenAI's model had a significant 20% lead on the top human, which felt like a lot.

At some point (probably right before toilet break), I decided to run my solution with a 10x larger time limit. That is 20s instead of 2s. Performing local run with an increased time limit is a common technique that gives you information about how well your solution scales with an increased time budget. To my surprise, the 10x time limit, increased the score by 25-30%. IIRC, that was slightly more than I needed to match OpenAI's score. In my code, the main bottleneck was the simulation code, everything else was negligible. 

Using this newly obtained information, I decided to focus the rest of my energy on squeezing maximum performance out of my sim code. Up to that point, my sim code had one major drawback: it simulated a single move at a time, while my fixed path usually contained 10-15 moves in a row in the same direction. If we want to move a robot in a single direction and we have information where is the closest wall in that direction, we can move that robot in a constant time, i.e. O(1). Maintaining that information is simple, considering I only change a single wall at each SA step, so I only have to update a single row/column after each step. 10 minutes later (based on my logs) I was done with implementation and my sims improved from 120K to ~700K, this gave me a 15% lead over the next human competitor, but I was still a little bit behind OpenAI. I remember, I was almost sure I'd be able to catch it up as long as I won't fall asleep. 

I was in a weird spot mentally. There was a mix of joy (I made a huge improvement), excitement (I should be able to overtake AI's top spot), relief (I was sure that my final rank would be decent, although I fully expected that someone might pass me right at the end) and exhaustion (I seriously considered taking a 20-60 minute nap, but I was afraid that it would completely mess up with my flow). Since I was short of only a few percent behind the model's score (and I have assumed it wouldn't improve much from that point), I decided to focus only on tiny changes. 

The remaining 4 hours were just a very slow grind on my side:
- At around 7:20, I finally managed to pass OpenAI's score. Something that I haven't mentioned before is that we had an opportunity to go to a "Confession Room", in order to give a monologue that will be shown on the livestream. I just needed to vent, so I decided to go as soon as I got the first place. I had no plan what to say, but I do think it turned out to be quite hilarious. You can see the whole thing [here](https://www.youtube.com/live/TG3ChQH61vE?t=26408s). For the record, I was so tired that I didn't even remember what I said until I looked at the recording.
- I had a lot of ambitious ideas, but due to tiredness I decided to skip anything that I couldn't implement in less than 5 minutes. Some of those ideas are mentioned in the "Possible improvements" section. In those contests, it's always the case that you want to do way more than your time allows. The hardest part is (a) to do it fast and (b) have a good intuition about what would work well. I'm usually good with both, but here I did terribly on (a) and exceptionally well on (b). 
- Most of my improvements came from very tiny changes: further speed improvement (I got ~1.5M sims in the end), small tweaks in the BFS phase, mixing linear & exponential temperature schedule in SA and some hyperparameter optimization. I have only gained around 15% from the 6 hour mark. 
- For the whole time, I was anxiously waiting for my human opponents to submit something that would suddenly overtake me. In the past, I've been in this situation many times and I always struggle with my mindset here: I always lose motivation as soon as I'm in the lead and then I'm just stressing out and nervously refresh the leaderboard looking at my shrinking lead. Looking at my logs, I did very little work throughout the last 4 hours, and I'm not sure if sleep deprivation is the only thing to blame. 

And that's it. When the timer hit the 10 hour mark, I was just happy that it's over and that soon I'll be able to recover some of my lost brain cells.


## Closing thoughts

In no particular order:
- I'd like to thank organizers and all of the remaining finalists. While I'm unbelievably happy that I managed to win, in the end it's important to remember that it's just one of the many possible outcomes. All of the finalists were incredibly strong and each one of them could easily win. I was lucky to get a problem that was well suited to my skillset. I was lucky that (AFAIK) all of my competitors fell into the trap of using multiple groups and thus struggled with overcomplicated approaches.
- While I'd like to think that my sleep deprivation effectively costed me 3-4 hours of wasted time (being slow & bugs) and without that my score would be 15-40% higher, it's also quite possible that if I was well-rested, I would try to implement a much more complex solution that would fail miserably. We will never know.
- If you've got this far and this whole contest sparked your imagination, consider participating in an online contest. It's fun if you like mental challenges.


## Possible improvements

A few of my ideas that could significantly improve my solution:
- My solution doesn't deal with pre-built walls at all. My fixed path doesn't guarantee that it will be possible to get to the targets. One way to fix this is that we can evaluate a fixed path as follows: for each robot, find how close we can get to its target if we only move along that fixed path; our eval is the sum of all of those distances. The goal is to take a path that has this sum equal to zero, so that (in theory) it's possible to find a placement of walls that every robot will end up in its final position without even doing additional BFS at the end.  
- When performing BFS, do A* instead to speed it up. This should allow for slightly more iterations of the final path-finding phase.
- Experiment with eval for the SA. Sum of manhattan distances is not the best, because robots in the corners get stuck more often. Making a more complex function that would take into account the absolute positions instead of just relative sounds promising. Another idea is to precalculate distances using pre-built walls instead of pure manhattan distance.
- The most ambitious idea I had was to add an additional greedy phase between the main SA and the final BFS thing. The greedy works as follows: Each iteration tries adding a single individual move at every possible moment (i.e. in the middle of the fixed path) for every robot in every direction. Eval is the sum of distances to targets, but since our walls are fixed now we can calculate exact distances. This should drastically improve hard tests, assuming it works fast enough (so I wouldn't have to sacriface time for the main SA). 

And if you want to see polished solution, check out the post-contest rankings [here](https://atcoder.jp/contests/awtf2025heuristic/standings/extended). Rafbill (one of the finalists) spend some time after the contest ended to create a program that is over 2x better than my final submission. A round score of 50M means, that currently (at the time of writing this) he holds the best solution for every test in the provisional scoreboard.














