## Links
[Problem Statement](https://atcoder.jp/contests/ahc017)

[Final Solution](https://github.com/FakePsyho/cpcontests/blob/master/topcoder/ahc017/main.cpp) - Final submission. Relatively clean.

## Intro

My initial thought was that this problem looks very straightforward: State representation is obvious and transitions should be very stable. This means that we need to speed up evals and perform some sort of local search. The problem is, that the constraints are unusually high and it's very hard to even fit a simple greedy within the time limit.

Additionally, pathfinding in general graphs is a somewhat cursed problem. There's a ton of potential algorithms, implementations are non-existant, every algorithm has it's pros and cons. It gets even more complicated if you consider that the graphs can be dynamic (recalculating distances if you add/remove edges). You can easily spent days on reading papers without really learning anything useful.

My initial goal was to limit myself to several hours with the contest and try to find an easy and efficient solution. This meant that I gave up on the pandora's box with pathfinding algorithms and focused on a simple greedy solution. Luckily, I got all of the important pieces very early on.

## Final Approach

1. Sort edges based on the added frustration.
2. Iterate over those edges and greedily add them to the day that minimizes the added frustration. 
3. The actual eval function that is used is `added_frustration * (1.0 - 0.3 * (edges_left/M)^3.5 * (day_edges_left/K)^12 * empty_day_bonus`, where `edges_left` = number of edges that we still need to use, `day_edges_left` = number of edges that we can still fit into this particular day, `empty_day_bonus` = 1.1 if the day has no edges so far, otherwise 1.0; the main reason for this eval function is that added frustration is generally smaller if particular day is already filled out, the goal of this function was to spread the edges somewhat evenly at the start of the process.
4. If there's still some time after the greedy finishes, perform a very simple hill-climbing: at each step I try to move a random edge to a different random day.


## Optimizations

* Each time frustration is recalculated, we only need to add/remove a single edge. For those cases, I'm using dynamic Dijkstra. This is by far the most important optimization that makes this solution possible. See `calc_all_distances_add_edge` & `calc_all_distances_rem_edge` for details.
* Divide all distances by `DIST_DIV` (1000 in my final solution) and use O(E + max_dist) Dijkstra insead of the "standard" O(V lg E). This sacrifices accuracy for improving speed. I'm not 100% sure if this indeed improves my score.
* If the test is big, use only a subset of starting vertices when calculating frustration.
* When calculating `added_frustation` in (2), progressively increase the accuracy (number of starting vertices). We start with a small number of starting vertices to limit the potential day candidates down to 5. Then we increase the number of vertices and reduce the candidates down to 4. Repeat until we're left with the best day according to an eval function in (3). This is called `beam` in my code, since this behavior resembles beam search.
* If the subset of starting vertices is small, vertices are selected in a way that they are evenly spread out throughout the whole graph. When this subset is large I'm just using uniform selection.
* In order to efficiently use all of the available time, I'm allowed to shrink/expand `beam` in the middle of the greedy. I haven't found any good way of doing this, so I'd say it's closer to "panic mode" when I suddenly shrink the size of the beam in order to avoid going over the time limit.


## Ablation Study

Results based on the seeds 2001-3000. Output from [psytester](https://github.com/FakePsyho/psytester). All of the tests' subgroups heavily correlate with overall score, since all tests were very similar.


```
Tests             1000        200        200        200        201         199         200          200          201          199          200     216      191      205      195      193
Run            Overall  N=500-623  N=624-723  N=724-814  N=816-909  N=910-1000  M=848-1258  M=1259-1491  M=1492-1674  M=1676-1959  M=1960-2914  D=5-10  D=11-15  D=16-20  D=21-25  D=26-30  Bests  Uniques   Gain
-------------  -------  ---------  ---------  ---------  ---------  ----------  ----------  -----------  -----------  -----------  -----------  ------  -------  -------  -------  -------  -----  -------  -----
final           99.370     99.525     99.448     99.435     99.316      99.127      99.503       99.474       99.352       99.332       99.192  99.299   99.321   99.375   99.443   99.421     18       18  0.001
tlx0.1          80.904     86.914     84.683     82.233     78.111      72.548      90.532       86.262       82.514       76.602       68.579  79.064   80.829   80.983   82.224   81.618      0        0  0.000
tlx0.5          98.678     99.305     99.101     98.781     98.362      97.841      99.236       99.090       98.837       98.542       97.686  98.745   98.658   98.612   98.739   98.633     14       14  0.003
tlx2.0          99.623     99.683     99.634     99.619     99.610      99.570      99.674       99.644       99.592       99.606       99.601  99.586   99.594   99.603   99.649   99.689     37       35  0.004
tlx4.0          99.805     99.867     99.815     99.764     99.785      99.793      99.853       99.799       99.780       99.798       99.793  99.815   99.780   99.786   99.807   99.834     77       60  0.005
tlx8.0          99.958     99.966     99.974     99.946     99.960      99.943      99.962       99.955       99.955       99.961       99.956  99.954   99.955   99.960   99.958   99.963    797      780  0.124
simple_eval     93.424     92.452     93.402     93.058     93.767      94.443      94.468       93.777       93.505       92.882       92.484  98.254   96.462   93.696   90.420   87.757      9        9  0.003
no_sorting      85.977     88.723     87.303     86.481     84.358      83.010      87.552       86.977       86.302       85.620       83.428  87.527   86.481   85.788   85.459   84.466      0        0  0.000
no_bs           98.379     99.099     98.814     98.291     98.018      97.672      99.060       98.690       98.362       98.060       97.722  99.035   98.465   98.221   98.132   97.977     14       14  0.005
uniformsubset   99.239     99.502     99.431     99.395     99.163      98.703      99.480       99.463       99.301       99.203       98.750  99.283   99.211   99.191   99.271   99.239     19       19  0.003
fixed_beam      99.219     99.526     99.442     99.333     99.028      98.765      99.504       99.403       99.271       99.151       98.763  99.257   99.180   99.212   99.229   99.212     18       18  0.004
no_hc           99.340     99.462     99.412     99.417     99.298      99.112      99.438       99.445       99.331       99.315       99.172  99.218   99.302   99.359   99.430   99.405      4        4  0.000
distdiv5k       97.066     97.920     97.617     97.228     96.695      95.866      98.769       98.040       97.651       96.760       94.107  97.090   96.894   97.040   97.151   97.153     12       12  0.004
```

* `final`: final submit
* `tlx?`: final submit with time limit multiplied by ?
* `simple_eval`: using `added_frustration` as eval instead of the one described in (3)
* `no_sorting`: removes (1)
* `no_bs`: no beam-search like progressive accuracy (mentioned in the Optimizations section)
* `uniformsubset`: removes "starting vertices are selected in a way that they are evenly spread out throughout the whole graph"
* `fixed_beam`: `beam` is fixed throughout the run. This includes changing few constants in order to avoid any TLEs
* `no_hc`: no hill-climbing phase at the end; the impact of HC is nearly negligible as I tried to use all of my allowed time on greedy
* `distdiv5k`: `DIST_DIV`=5000; shows the potential loss due to inaccuracy of distances


## Comments

* As you can see from the ablation study, it was critical to have: very fast pathfinding, correctly sorted edges and a good eval function. Missing even one meant that you were not able to achieve a very high score.

* I have not explored dynamic edge order in my greedy. I think it's possible to do that within the time limit. At least for the smaller tests.

* I don't know how much I'm losing because of the linear Dijkstra. I wouldn't be surprised if in the end I'm actually losing points due to the lost accuracy.

* I think I have a small of timing out, although it's hard to estimate probability of that. Also, every once in a while (1 test out of 20K?) my code produces a disconnected graph which essentially gives 0 for that particular test.

* It was relatively easy to get a complete set of bests in the provisional 50 tests, because (a) all tests are similar, there are no corner cases (b) (I think) greedy (+ optional hill-climbing) is the only viable approach here. My big jump from ~45M to 50M was when I switched from normal Dijkstra to dynamic version. At that point I had every other component of my solution that is mentioned above. Remaining improvements came mostly from small tweaks in eval and better time management.


