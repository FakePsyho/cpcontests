## Links
[Problem Statement](https://codeforces.com/contest/1751/problem/A)

[Final Solution](https://github.com/FakePsyho/cpcontests/blob/master/various/huawei2/main.cpp) - I have removed a lot of nearly useless (~0.4 gain) code to make it a little bit more readable.
## Problem
This problem is a fairly standard task scheduling problem. The biggest issue here is that the best solving method highly depends on the input. Large tests, small tests, unbalanced allowed machines, very strict disk capacities, narrow critical path are some of the potential properties that can drastically affect the preferred solving method. You could either randomly try a bunch of methods or do some heavy work with the test case extraction to make more educated guesses.

## Approach
Out of 51 tests, only two were maximal (n=10000) and every other test was n<=800. Because of that, I wanted to find a decent approach based on the local search. I tried a weird mix of SA (Simulated Annealing) and greedy. 
My state is: `task_order` + `task_machine` + `task_disk_preference`, where each variable is a vector of length t (number of tasks).
In order to evaluate the state I use a greedy method: 
- Use `task_disk_preference` to pre-assign a disk to each task, earlier in the list gets a faster disk.
- In each step, take the first (based on `task_order`) available (= all dependencies fulfilled) task and place it on the designated machine (based on `task_machine`) using the pre-assigned disk.
- It's allowed to set machine to "unknown" for a task. This version tries all machines and chooses one that has the fastest finish time.

This state representation is fairly stable, but the state evaluation is incredibly inefficient (since it restarts the whole process from scratch). For n=800 I got around 50K state evals (which is rather low, but not terrible), but for n=10000 I got only around 2K which is pretty much nothing and thus SA has no chance of improving it.

The main improvement came from initializing the `task_order` based on the "critical path". For each task, I calculate the longest possible chain of tasks that have to follow it. The tasks with the longest chains are given the highest priority. For some reason using "longest_chain_starting_from_this_task - longest_chain_ending_with_task" gave even better results so I kept it. Initializing `task_disk_preference` was simpler and I just used data*number_of_times_this_data_is_used. This + some minor tweaks gave around 2827 points. The remaining 10 points came from running the same solution many times and choosing the best set of parameters for each test case.

My greedy initialization is far from good and my SA does nothing for the 2 largest test cases, so I'd imagine that's where the remaining 10 points mostly come from. 

## Additional Comments
If someone started early, it was possible to extract some of smaller tests (all n<=300 and maybe even some/all of the n=800). You get a total of 8064 (24 per hour * 24 hours * 14 days) submits throughout the 2 weeks of the competition. With each submit you can extract 24-32 bits of information. I'd guess that should be achievable with a decent compression (removing redundant input + LZ compression). This would leave only the last 2 cases for online solving.




