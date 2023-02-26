## Links

[Problem Statement](https://atcoder.jp/contests/ahc018)

[Final Solution](https://github.com/FakePsyho/cpcontests/blob/master/atcoder/ahc018/main.cpp) - Final submission. 

[Clean Solution](https://github.com/FakePsyho/cpcontests/blob/master/atcoder/ahc018/main_19.cpp) - 19th submit. Score is ~3% lower but the code is much simpler (around 1/3 shorter). I'd recommend looking here first and then looking at the diff between this and the final one.


## Intro

That was a really cool problem with a simple concept and yet a lot of depth. It was heavy biased towards problem solving in the sense that algorithmic knowledge and code optimization skills were not very useful here. I won't lie, I was expecting very tightly packed top so I'm highly surprised about the results. I believe the variance is actually somewhat small here (my leaderboard results were usually within 1% of what I expected), so the final ranking should be relatively close to the provisional one (I hope!).

I'm pretty sure that my approach is almost the same (or simpler) when compared to everyone else in the top 20. My guess is that the difference comes from optimizing parameters in order to minimize the stamina usage. Overall, you have to balance between spending stamina on probing/scanning the terrain and building the actual path. Those two are tightly connected and unfortunately there's no other way than brutally optimizing your hyperparameters.

I wanted to have a simple and fast approach because it translates into quick testing. On a 24-core VM, 10K tests (which was enough to find tiny improvements) took me 30-40 seconds. In order to optimize a single parameter I usually made a batch with 10-20 versions and then used optimal values for each C. I believe I made around 300-400 total runs of 10K tests. 

If you want to understand tiny details, I highly recommend looking at the code since it's short & simple. It would be nightmare to explain every little magic number here.


## Final Approach

![](https://github.com/FakePsyho/cpcontests/blob/master/topcoder/mm143/seed100047.mp4)

1. Create a graph consisting of water sources, houses and predetermined grid points. Cost of an edge is based on the manhattan distance between those two points and dig stats so far on those two specific cells.
2. Find the shortest path between any source (either water source or a tile where water already flows) and any house. On that path, choose the vertex that's not yet fully excavated. Tie-breakers: (a) smaller number of digs so far (b) closer distance to the target.
3. Repeat (2) until we get a path that only contains fully excavated vertices. After that, we build a full path from water to the house and we excavate tiles one by one sequentially. 
4. Repeat (1-3) until we're done.

And that's it. Most of the heavy lifting is done via parameter optimization. That is, instead of writing complex logic I opted for very optimized parameters. Functions for calculating edge costs in the graph and finding best power when digging the actual paths look very random, but that's what I got by trying to extract maximum value out of those.


## Potential improvements

- My parameters only depend on C. There are no adjustments towards the way Perlin noise is generated, which is probably the biggest potential improvement I could do. The problem here is that it's hard to extract exact hardness of each cell, especially for high C. At the same time, this improvement would mostly benefit high C tests. If I had to continue improving my solution, that would be my next step
- My grid points are fixed and there's no "smoothing-out" phase after the path is computed. This means that paths around water sources and houses are highly unoptimal. If I were to fix this, I'd probably "distort" grid around each water source/house so that each water source/house would translate into an actual vertex of the grid. TBH, I'm not sure if this would yield any substantial improvement.
- It's unclear to me how much my solution would improve by switching from greedy to searching for rectilinear Steiner trees, but that would require a significant rewrite of my solution. My intuition was that this wasn't worth it considering this would significanltly slow down testing.
- Ideally, I'd have a function that is able to estimate the distribution of the hardness for every existing cell instead of using my `edge_cost` function. This either requires applying Bayes (I'm bad at math) or using an ML approach (a lot of work), so I decided to leave it as it is. And again, this would slow down the testing significantly.
- It's worth mentioning that this was probably the first situation that I have encountered where having a tool for automatic parameter tuning (which I don't have) gave a huge advantage. I did consider writing one and/or using existing one. 


## Ablation Study

Results based on the seeds 100000-109999. Output from [psytester](https://github.com/FakePsyho/psytester). I haven't used those seeds during parameter optimization, so there's no overfitting.


```
Tests                  10000    1192    1255    1213    1312    1228    1261    1257    1282    1973    2003    1971    1956    2097    2475    2461    2453    2611
Run                  Overall     C=1     C=2     C=4     C=8    C=16    C=32    C=64   C=128   K=1-2   K=3-4   K=5-6   K=7-8  K=9-10     W=1     W=2     W=3     W=4  Bests  Uniques   Gain
-------------------  -------  ------  ------  ------  ------  ------  ------  ------  ------  ------  ------  ------  ------  ------  ------  ------  ------  ------  -----  -------  -----
final                 93.868  93.193  93.012  93.174  93.218  93.793  94.278  94.768  95.443  91.565  93.741  94.207  94.690  95.073  94.248  93.839  93.611  93.777   2018     1961  0.217
sub_19                90.975  89.985  90.562  89.902  90.627  90.917  91.485  91.844  92.371  89.240  91.099  91.232  91.603  91.661  91.422  91.022  90.803  90.667    940      939  0.239
square_grid           90.309  88.918  88.489  89.234  89.217  89.717  91.258  92.343  93.155  88.352  89.903  90.365  91.135  91.713  90.408  90.209  90.295  90.320   2035     2035  1.232
dig_vertex_last       90.996  88.859  89.381  89.513  89.577  90.970  91.921  93.252  94.323  88.429  90.075  91.343  92.228  92.816  90.759  90.995  90.934  91.281   1448     1398  0.455
dig_vertex_first      88.750  86.167  87.011  86.851  86.755  88.431  90.093  91.433  93.049  86.148  87.659  89.323  89.971  90.564  89.363  88.497  88.382  88.754   1099     1080  0.419
avg_magic_constants   92.096  90.191  91.268  91.949  93.165  93.244  93.561  92.311  90.976  89.688  91.919  92.508  92.853  93.439  92.529  92.169  91.715  91.976   1250     1249  0.329
no_adjust             91.987  93.026  92.760  92.617  92.416  92.204  91.785  90.930  90.253  89.782  91.827  92.398  92.750  93.114  92.262  91.925  91.747  92.009    197      184  0.010
simple_edge_cost      90.897  88.707  89.469  89.311  90.164  90.726  91.862  92.737  93.993  87.998  90.327  91.190  92.111  92.762  91.223  91.133  90.438  90.797   1106     1091  0.260
normal_dig            80.923  86.353  85.692  82.387  81.624  79.553  77.925  76.751  77.458  79.042  80.774  81.357  81.499  81.892  81.373  80.943  80.649  80.737      2        2  0.000
```


- `final`: final submit
- `sub_19`: submit #19 that's linked at the top
- `square_grid`: change the grid to normal "square" grid (one point every 12x12)
- `dig_vertex_last`: in (2) remove (a) tie-breaker
- `dig_vertex_first`: in (2) remove (a) tie-breaker and change (b) to the closest vertex to the source
- `avg_magic_constants`: average all of the magic constants that are split by C (`CLEV`); this shows how much I was able to gain by adjusting parameters for each possible value of C
- `no_adjust`: no adjustment for `dig_lo`/`dig_hi` based on the previous cell; when sequentially digging I try to narrow down the range of possible hardness based on the range of the previous cell
- `simple_edge_cost`: simplified `edge_cost` function to distance * average_hardness, where average_hardness is (`dig_lo` + `dig_hi`) / 2 + (150 if cell is not fully extracted)
- `normal_dig`: when digging sequentially, use the same function that is used for exploration; essentially removes the advantage of the sequential digging

## Other Comments

- I believe the reason why my solution performs very good despite being simple is that (a) I'm wasting very little stamina on exploration (due to heavy optimization of the parameters), (b) sequential digging is probably the best because you can use the result of the previous dig in order to adjust and (c) adjusting parameters per C (look at ablation study to see how much more value you can extract by doing so).
- There are few small tricks that I haven't mentioned, but their overall impact is very low. Probably the most important one is that when my calculated path has a diagonal movement, I make an additional last-minute exploration to see which exact route should I take. IIRC this was a ~1% improvement. 
- I think I'm saying this every contest, but I'd highly recommend having a good local tester. If you're not using one, you're essentially wasting your time and learning at much slower pace. If you're using psytester remember that you should compile rust sources first with `cargo build --release` and then use the generated executables.
