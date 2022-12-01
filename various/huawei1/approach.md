## Links
[Problem Statement](https://codeforces.com/contest/1752/problem/A)

[Final Solution](https://github.com/FakePsyho/cpcontests/blob/master/various/huawei1/main.cpp) - the code is a complete mess since there was no reason for any refactoring
## Problem
The problem is fairly bloated, but it's turns out it's an interesting and rather unusual puzzle. Each test case consists of only 8 numbers, so you can immediately extract them with just a few submits (and pretty much everyone did exactly that). This means you can save CF some server costs and just solve the problem locally. It was also rather clear that it should be reasonable easy to get perfect/near-perfect solution.
Test cases can be separated into two groups: n=4, high k and high n, high k. Apart from that, there are a few additional "outliers".
## Approach
- As a starter, simple brute force to solve three small tests & verify that I can correctly calculate all of the stats.
- For high n (n>100), high k, it looked like any reasonable greedy approach should work:
	-  Let's build our graph sequentially. We start with an "empty" graph and we greedily add components. I have two kinds components:
		- one-vertex component: X->X,X+1; this creates 0...01/1...10 paths, which should solve "right" & "change" stats.
		- two-vertex component: X->X+1,X+1 & X+1->X,X+2; this allows any possible path for the first vertex & and 0...01/1...10 for the second one; the goal was to lower correlation stat, so the path in the first vertex is completely random (random sequence should have very low correlation)
	- At each step I randomly generate multiple random components and I take the one that minimizes evaluation function. The last component uses all of the remaining length of the route (k input).
	- This solves all of the big n tests on the first try; The smaller ones (n=100-300) require multiple runs and some luck (this could be solved by using beam search instead of greedy if there was a need for this).
- Tests with small n (n=4) are conceptually much simpler. I randomly generate graphs and run beam search on them, where each step is an addition of a single segment. This worked well, but unfortunately was quite slow. It turns out that you don't really need beam search as greedy (again) is enough as long as you had a decent start & end (so you can finish your path in the correct vertex). The fix was to reduce the beam width to 1 for everything except for the first and last ~100 steps. 
- Evaluation function in both greedy and beam search is MSE (Mean Squared Error). MSE works well because you need to give more priority for stats that are the furthest away from the goal.
- Outlier tests:
	- weirdly enough, #64 (n=152, k=198) was solved with brute force (unless I remember it wrong)
	- #65 (n=407, k=535) was solved via greedy with some minor tweaking
	- #100 (n=8, k=79), #101 (n=2, k=75) & #103 (n=17, k=252) were solved with beam search with a larger beam width. My guess is that #103 was probably the hardest cases for most people
	- #98 & #105 were too large to solve via beam search/greedy methods within 15s; By eyeballing the stats, I figured out that a full cycle graph (i.e. X->X+1,X+1) with uniform random probably solves those tests. It turned out to be correct.
## Additional Comments
If I had to solve this problem once again I would completely change my approach. #98 & #105 were the very last tests that I have solved. After solving them, it struck me that this "cycle" graph structure should trivialize every other test case. Since we always move 1 vertex "forward" with segment, we have a fixed path through vertices that doesn't depend on our route. This means we can apply simple hill climbing with O(P) evaluation time for every change. This should be very easy to implement and would solve all of the tests except maybe 1 or 2.

Overall, the problem turned out to be way more interesting than I anticipated (especially based on a rather bland problem statement). I only wish the testcases would be much harder (or maybe required higher precision?). I started 4 days late. Having someone finished half a day after starting was a bit demotivating. 



