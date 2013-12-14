## Origin profiler

To generate the data run from shell:

```
python src/profiler.py >results/acs_new.txt
```

To generate plots:

```
source("src/plot_profile.R")
```

Results:

![Aggregated skew](https://github.com/bdomokos74/OriginProfiler/raw/master/results/profile_skew.png)

![Aggregated profile](https://github.com/bdomokos74/OriginProfiler/raw/master/results/prf_aggr.png)

![A specific profile](https://github.com/bdomokos74/OriginProfiler/raw/master/results/prf1.png)
