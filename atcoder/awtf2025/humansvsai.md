## Intro

You saw the headlines that AI is now the 2nd best competitive programmer, but you don't know what to make of it? This is the place for you. While it's hard to draw any definite conclusions out of a single contest, there's quite a bit we can infer from this.

I highly recommend reading the [write-up](https://github.com/FakePsyho/cpcontests/blob/master/atcoder/awtf2025/approach.md) for the contest first, otherwise you won't get the full picture. As an absolute minimum, you should know the [exact format of the contest](https://github.com/FakePsyho/cpcontests/blob/master/atcoder/awtf2025/approach.md#exact-format-of-the-contest).


## Detailed look at heuristics contests

Heuristic contests resemble research. You get an open-ended multi-layered optimization problem. Unlike research, the problems are well defined: you get the full specification of the problem, scoring formula and complete process of creating the test data. You also get the code for scoring & test generation so you can fully recreate the judging system on your local machine.

Another way of looking at it is: we get the full specification of a single player (puzzle) game and our goal is to write a computer program that will play that game. Each test is just a randomly generated level. In fact, many contests are directly inspired by computer games. 

Properties of a well-designed heuristic problem:
- Finding a good approach should be hard even for the most experienced competitors
- Wide variety of possible approaches, i.e. there's no single solution archetype that completely dominates the top rankings
- Practically impossible to obtain optimal answers (except maybe for a tiny subset of the tests)
- The problem consists of multiple layers and has a lot of depth

For example, we could use the [classic NES tetris](https://www.youtube.com/watch?v=RlnlDKznIaw) as a starting point for the problem. Obviously it's way too simple, but we can add complexity via more sophisticated test generation. Instead of using classic tetromino set, each test would use randomly generated polyomino set. We can randomize the width instead of having it fixed at 10. Randomize the number of next visible polyominos. As a final touch we may add some powerups. Is it good enough? Probably not 😅 Experienced competitors would immediately know that either [MCTS](https://en.wikipedia.org/wiki/Monte_Carlo_tree_search) or [Beam Search](https://en.wikipedia.org/wiki/Beam_search) would be the best here (depending on the exact specifics of the problem). But at least, you have a rough idea of what you can expect.

Being good at heuristic contests means that you can perform well in a very wide variety of problems. It's almost pure problem solving. 

Solving heuristic problems is an iterative process and varies from problem to problem and from person to person, but usually it incorporates some of the following:
- write a few simple approaches, analyze their performance in order to build intuition for the problem
- do some deep analysis using pen & paper
- isolate a subproblem and find existing research about it
- stare at the visualizer to find problematic behavior
- optimize your code: design better algorithm, better efficiency, pruning, precomputing common calculations, reduce memory usage, etc.
- refactor / simplify your code; doesn't help directly, but makes the future work easier
- perform some tweaks, either by implementing tiny changes or through [hyperparameter optimization](https://en.wikipedia.org/wiki/Hyperparameter_optimization)
- benchmark your current solution

I hope this gives you a clearer idea of the complexity of these problems. And so, if these problems are so complex how did OpenAI manage to get such a good result? 


## Expectations

The current SOTA models are very good at some of those tasks (tiny tweaks, basic solutions), okay with others (applying existing research to a subproblem, code optimization) and meh with the rest (they can't use a visualizer, deep analysis for open-problems is hard when you don't know what you're looking for). The place where the agent shines is its ability to perform all of those things at the same time. It can simultaneously explore many different changes and it's only limited by the compute budget.

While I don't know exactly how OpenAI's agent works, I can safely assume it's quite similar to [AlphaEvolve](https://deepmind.google/discover/blog/alphaevolve-a-gemini-powered-coding-agent-for-designing-advanced-algorithms/) from DeepMind and [ALE-Agent](https://sakana.ai/ale-bench/) from Sakana AI (their ALE-Agent would [place 4th](https://x.com/SakanaAILabs/status/1945845946115661873) in that contest, if you don't count OpenAI's submission). 

The way it (most likely) works:
- maintain a set of diverse solutions
- take: problem statement + one of the solutions + its local testing results (scores, distribution of the execution times) + prompt motivating it to improve the solution + (optionally) other solutions (a summary is probably enough); you can have a wide variety of hardcoded prompts, but it's also possible to have another agent that will act as the master
- use the above to generate a new solution & run local testing for that solution
- every once in a while, submit your best solution
- repeat until the end of the contest

It's easy to see that the fact that local testing gives you a precise score, allows you to generate thousands of different versions, test them and reliably select the best ones. Based on some intuition and Sakana's results, I was expecting for OpenAI to quickly get some decent scores and then essentially plateau. The bigger the compute budget, the faster it's going to get to its maximum potential.

Given all of that, [this](https://x.com/FakePsyho/status/1945109137379221713) was my take before the contest started.


## vs Reality

Let's go over each aspect:

> - problem itself

Since the contest was 10 hours long instead of 10 days + it had to be easy to understand for the livestream spectators, it was expected to be on the simple side. There was no need for deep analysis and you could safely ignore the visualizer. At the same time, the problem was not straightforward: after all, I won because I managed to successfully simplify the problem (thanks to my intuition). I'd say that this problem slightly favored AI compared to the pool of possible problems we might have gotten at the onsite finals. If you compare it to the pool of possible problems that we usually encounter during contests that are 10 days long, then it heavily favored AI.

> - compute budget

Most probably, we will never know 😢 If I'm not mistaken, Sakana AI's system had a compute budget of a few thousand dollars (I recall seeing this somewhere, correct me if I hallucinated this information).

> - testing budget

Testing budget is our budget for running local tests. It's separate because it's a limiting factor for humans as well. This limits how quickly the agent is able to evaluate solutions.

> - variance in human performance (non-JP will be sleep deprived)

I do feel that human competitors drastically underperformed in this contest. Historically, in order to win I needed both a good approach and a good implementation of that approach. Even during the contest I was surprised by how big of a lead I had. The idea that ultimately secured my victory (in order to achieve near-optimal results for small tests, you have to use a single group + fixed path is okay-ish + SA works well for the wall placement) was relatively quick to figure out and is a rather standard concept. If I had competed in isolation (without access to the live leaderboard), my guess would be that I have finished somewhere between 2nd and 4th positions.

> - amount of work put towards finetuning agent for AHCs

As far as I know, only a few people were involved with this project directly. I'd imagine that after building the basic system for maintaining diverse set of solutions, additional work gives you diminishing results. Good and varied prompting might help, but if the underlying reasoning model is very strong, perhaps it isn’t necessary 🤷


## OpenAIAHC submits

I have downloaded all of the OpenAI's submissions, auto-formatted them and run them locally on 2000 tests (same number that atcoder used for their final results). Below is the output from [psytester](https://github.com/FakePsyho/psytester/), which is my tool for performing local tests.

| Run          | Overall |    W=0 |    W=1 |    W=2 | K=10-31 | K=32-56 | K=57-77 | K=78-100 | Bests | Uniques | Fails | Time    | Base |
| ------------ | ------- | ------ | ------ | ------ | ------- | ------- | ------- | -------- | ----- | ------- | ----- | ----    | ---- |
| ai00         |  13.307 | 11.848 | 13.165 | 14.962 |  24.945 |  11.679 |   9.656 |    6.833 |     0 |       0 |     0 | 0:15:01 | - |
| ai01         |  33.457 | 31.881 | 33.262 | 35.286 |  33.457 |  28.767 |  34.947 |   36.824 |     0 |       0 |     0 | 0:24:42 | - |
| ai02         |  51.678 | 53.813 | 51.02  | 50.073 |  46.598 |  43.814 |  56.206 |   60.455 |     0 |       0 |     0 | 0:49:10 | - |
| ai03         |  55.448 | 55.775 | 55.054 | 55.479 |  50.248 |  47.141 |  60.128 |   64.651 |     0 |       0 |     0 | 1:08:17 | - |
| ai04         |  58.356 | 58.396 | 58.109 | 58.547 |  49.054 |  50.055 |  64.846 |   69.915 |     0 |       0 |     2 | 1:41:13 | - |
| ai05         |  62.077 | 62     | 61.694 | 62.52  |  50.736 |  52.547 |  69.616 |   75.93  |     0 |       0 |     5 | 2:02:11 | ai04 (10%) |
| ai06         |  62.652 | 61.888 | 62.361 | 63.725 |  52.063 |  54.007 |  70.362 |   74.675 |     0 |       0 |     3 | 2:35:41 | ai05 (90%) |
| ai07         |  65.282 | 65.193 | 65.301 | 65.357 |  55.204 |  54.37  |  72.317 |   79.771 |     0 |       0 |     2 | 3:07:35 | ai05 (10%) |
| ai08         |  66.417 | 65.496 | 66.406 | 67.39  |  56.435 |  55.505 |  73.61  |   80.658 |     0 |       0 |     1 | 3:37:43 | ai06 (50%) |
| ai09         |  66.973 | 66.509 | 66.837 | 67.587 |  59.612 |  56.435 |  73.095 |   79.239 |     0 |       0 |     0 | 4:06:33 | ai03 (90%) |
| ai10         |  66.203 | 65.659 | 66.099 | 66.869 |  56.236 |  54.899 |  73.548 |   80.682 |     0 |       0 |     1 | 4:16:31 | ai08 (10%) |
| ai11         |  67.232 | 66.02  | 67.159 | 68.565 |  57.141 |  55.822 |  74.216 |   82.295 |     1 |       1 |     3 | 4:39:07 | ai08/ai10 (10%) |
| ai12         |  68.332 | 68.269 | 68.14  | 68.58  |  58.883 |  57.902 |  75.321 |   81.742 |     0 |       0 |     0 | 5:07:36 | ai09 (50%) |
| ai13         |  71.777 | 69.262 | 72.147 | 74.05  |  58.609 |  59.731 |  79.965 |   89.414 |     9 |       8 |     2 | 5:35:20 | ai07 (95%) |
| ai14         |  72.739 | 70.236 | 73.416 | 74.71  |  57.969 |  60.22  |  81.755 |   91.668 |    23 |       7 |     2 | 6:04:44 | ai13 (5 lines) |
| ai15         |  72.728 | 70.127 | 73.334 | 74.869 |  57.967 |  60.232 |  81.845 |   91.526 |    22 |       6 |     2 | 6:23:03 | same as ai14 |
| ai16         |  73.222 | 70.791 | 74.01  | 75.012 |  58.648 |  60.955 |  82.143 |   91.785 |    17 |      15 |     3 | 6:31:12 | ai15 (10%) |
| ai17         |  73.322 | 71.53  | 73.658 | 74.874 |  58.352 |  60.76  |  82.365 |   92.469 |    42 |      35 |     1 | 7:17:06 | ai16 (3 lines) |
| ai18         |  73.441 | 70.842 | 74.172 | 75.461 |  61.994 |  61.068 |  81.09  |   90.209 |    11 |      10 |     5 | 7:25:46 | ai07/ai13? (95%) |
| ai19         |  73.404 | 71.663 | 73.61  | 75.026 |  61.685 |  61.146 |  81.246 |   90.138 |    19 |       1 |    10 | 7:51:03 | ai18 (1 line) |
| ai20         |  73.783 | 71.965 | 73.834 | 75.631 |  62.024 |  61.205 |  81.437 |   91.067 |    25 |       6 |    12 | 7:56:14 | ai19 (10%) |
| ai21         |  74.396 | 72.326 | 74.562 | 76.399 |  62.511 |  61.259 |  82.541 |   91.907 |    23 |      17 |     0 | 8:30:58 | ai17 (30%) |
| ai22         |  74.686 | 72.057 | 75.433 | 76.723 |  62.198 |  61.032 |  83.252 |   92.923 |    61 |      51 |     2 | 9:11:35 | ai21 (10%) |
| ai23         |  75.652 | 73.183 | 76.137 | 77.77  |  62.93  |  61.776 |  83.806 |   94.753 |   123 |       9 |     0 | 9:20:58 | ai22 (50%) |
| ai24         |  74.988 | 72.324 | 75.046 | 77.712 |  62.64  |  62.528 |  82.643 |   92.743 |    20 |      14 |    13 | 9:27:42 | ai20 (10%) |
| ai25         |  74.964 | 72.478 | 75.144 | 77.389 |  61.735 |  61.128 |  83.779 |   93.892 |   103 |      92 |     1 | 9:32:50 | ai22 (10%) |
| ai26         |  74.569 | 71.553 | 74.658 | 77.631 |  62.639 |  62.792 |  83.084 |   90.37  |    82 |      70 |    31 | 9:37:55 | ai24 (10%) |
| ai27         |  75.621 | 73.175 | 76.091 | 77.727 |  62.931 |  61.754 |  83.8   |   94.654 |   121 |      10 |     0 | 9:42:56 | same as ai23 |
| psyho        |  83.131 | 88.348 | 82.341 | 78.435 |  89.921 |  85.45  |  85.239 |   71.89  |   302 |     273 |     0 |         | - |
| ai_psyho     |  90.653 | 89.82  | 91.076 | 91.123 |  96.884 |  92.056 |  91.125 |   82.501 |   592 |     514 |     0 |         | psyho |
| ai_psyho_fix |  91.806 | 91.594 | 91.62  | 92.204 |  97.054 |  93.29  |  93.199 |   83.665 |   689 |     606 |     0 |         | ai_psyho |

Columns: Overall is total relative score (scaled to 100.0 max; I will never understand why atcoder likes to use huge numbers), W=... & K=... represent the total relative score for the subset of tests with variables in the corresponding range, Fails is the number of failed tests (crashes or invalid output), Base is the previous version that this version is based upon and the number in the parenthesis is the rough estimate of the changed code (so 95% essentially means a near complete rewrite)

Rows: aiXX are submits made by the OpenAI during the contest, psyho is my final submit, ai_psyho is the OpenAI's submit that improved my solution that will be described in details in the next section

Something that my local tester doesn't support is that I don't detect when the program goes over the time limit. During the contest, many of OpenAI's later submits went over the stated 2s time-limit (which would normally result in 0 score for those tests). As a result, some OpenAI submissions appear slightly better than they actually were.

Few important observations about the submits:
- The introduced changes vary wildly in terms of quantity. Many submits are near-complete rewrites. Even when the changes are reasonable (10-30%), they usually touch multiple unrelated sections of code. Additionally, there are a few submits that are equivalent to changing a constant. 
- After a while, there are only two main solutions for the agent: ai26 & ai27 are both based on beam search, but they have no shared code; this supports maintaining a wider set of solutions that are developed independently
- There are two things that do not make sense: submits ai24-26 are strictly worse than ai23 & ai15 is exactly the same as ai14. Similarly, the very last submit ai27 is just a resubmit of ai23 (which is indeed the best submit). The simplest approach to submitting solutions would be to periodically (e.g. every 15 minutes) take the best scoring solution and submit it, but this should never result in submitting the same code twice and it should never result in submitting something strictly worse (1% difference is quite a lot in this contest and easily noticeable with 1K+ local tests). I can find only two explanations for such behavior: either sloppy programming by the OpenAI team or the agent itself made the decision about which solution should be submitted. Another possible explanation could be manual submissions by the OpenAI team (which would explain ai15 as human error), but my impression was that the system was fully automated.
- The code gets progressively more complex and less maintainable. Most of the late improvements are rather small tweaks. The agent nearly plateaued during second half of the contest. There's only a 4% improvement during the last 4 hours.
- Later submits had very little whitespace. My wild guess is that it's just a weird temporary side effect of rewriting the same code many times in the row: at some point the agent randomly "eats" some of the whitespace and then it maintains consistent coding style.
- All of the magic numbers in the code tend to be very round. Rather than finding the best-fitting constants, agent had the tendency of addings more complex code.

Based on my logs it looks like the version I have submitted right before going to the "Confession Room" (so around 7:15 into the contest) would've narrowly won against ai27. It wasn't as back-and-forth as the provisional standings indicated during the contest.


## Agent 🤝 my final submit

Around 24 hours after the contest ended, OpenAI submitted an improved version of my code. I carefully analyzed the changes and... I was quite impressed (something that doesn't happen too often 😅). [At the end of my write-up](https://github.com/FakePsyho/cpcontests/blob/master/atcoder/awtf2025/approach.md#possible-improvements) I mentioned several possible improvements and the agent implemented two of them 😲

Summary of changes (from the most impressive to the least impressive):
- In my main SA, the agent changed eval to use precomputed distances using pre-built walls instead of pure manhattan distance. This is why ai_psyho gets much better score for W=1 & W=2.
- Changed my BFS to A*, thus (on average) it increases the number of BFS hill climbing iterations at the end.
- Changed `TIME_SCALE` to 0.9. This constant controls how much time my solution uses. My original solution stops at around 1.94s (max is 1.944s on the atcoder server). Changing this constant to 0.9 means using about 10% less runtime. ai_psyho_fix is ai_psyho with `TIME_SCALE` set back to 1.0. I suspect that OpenAI hardcoded a guideline for the agent to keep the execution time below some threshold (most probably somewhere between 90-95% range). While my solution never gets close to 2s, it was probably still above that threshold. 
- Unrolled some of my REP macros for no reason.
- Added a lot of dead code. There's a giant block of code wrapped around something that's essentially an `if (false)` block.
- In some places I carefully handled edge cases to avoid pointless ifs. Agent added "guards" for those edge cases, making the code both slower and reducing the "entropy" of the code (fancy way of saying that it reduced the implied understanding of the problem).
- Several of the introduced changes add little to no benefit, while substantially increasing the complexity of the code.

All of those changes are rather isolated (unlike submits during the contest), which hints that (a) very little compute was used for this (b) the agent that was used during the contest, probably used way more compute than we might think 😅

It might be worth noting that the agent's code tends to be slow & unoptimized. My guess is that's what the training set looks like and thus the over-reliance on STL is the default coding style of the agent.  


## Conclusions

- As of July 2025, the AI agent can quickly create a decent solution, but then quickly reaches a plateau. I'm confident nearly all of the remaining finalists would easily beat the agent's score given more time. The problem used at the onsite finals heavily favored the AI when compared to the average problem used in an online contest that spans 10 days. To be clear, this comparison involves only the very best competitors. The gap between the top and the median is very wide. Your average competitor will lose on most of the problems, most of the time, even when given a lot of time to develop their solution. 
- A big issue with the agent is that it heavily favors short-term goal (score) while disregarding long-term planning (simplicity, maintainability). This is an interesting case of exploration vs exploitation dilemma. Based on my experience, beginner competitors always struggle with this, and it takes a lot of experience to even become decent at balancing the two. One could argue that this is one of the core skills of an experienced programmer. I wouldn't be surprised if "solving" this issue would result in a dramatic increase in the agent's performance. If the exploration doesn't work correctly, adding more compute won’t help significantly, since the solution space will eventually collapse.
- Tied to the point above, heuristic contests might be an interesting benchmark for AI agents. This is because it's one of the rare cases where you have a complex open-ended problem, but in the end you receive precise and immediate feedback on your performance.
- It seems the optimal combination for the not-so-distant future would be an AI agent doing all of the grunt work, guided by the human expert. The main obstacle is that most of the AI-generated code tends to be messy and bloated. If you look at the scoreboard, it's clear that the number of failed tests increased throughout the contest. That's not a good sign.
- Assuming we have more heuristic contests in the future in which frontier AI labs showcase their agents, I hope we will at least establish clear rules regarding the testing budget. I also hope this won't be the last time humans win against AI in heuristic contests.

And that's it. Hopefully this was a useful datapoint for our AI timeline.
