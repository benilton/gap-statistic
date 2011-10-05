## An implementation of the gap statistic algorithm from Tibshirani,
## Walther, and Hastie's "Estimating the number of clusters in a data
## set via the gap statistic".

## Given a matrix `data`, where rows are observations and columns are
## individual dimensions, compute and plot the gap statistic (according
## to a uniform reference distribution).
gap_statistic <- function(data, min_num_clusters=1, max_num_clusters=10, num_reference_bootstraps=10) {
    num_clusters <- min_num_clusters:max_num_clusters
    actual_dispersions <- sapply(num_clusters, function(n) dispersion(data, n))
    ref_dispersions <- do.call(rbind, lapply(num_clusters, function(n) reference_dispersion(data, n, num_reference_bootstraps)))
    mean_ref_dispersions <- ref_dispersions[ , 1]
    stddev_ref_dispersions <- ref_dispersions[ , 2]
    gaps <- mean_ref_dispersions - actual_dispersions
    c(K=num_clusters[which.max(gaps)], gaps=gaps, gap_stddevs=stddev_ref_dispersions)
}

## Calculate log(sum_i(within-cluster_i sum of squares around cluster_i
## mean)).
dispersion <- function(data, num_clusters) {
    ## R's k-means algorithm doesn't work when there is only one cluster.
    if (num_clusters == 1) {
        log(sum(sweep(data, 2, colMeans(data))^2))
    } else {	
        ## Run the k-means algorithm `nstart` times. Each run uses at
        ## most `iter.max` iterations.
        k <- kmeans(data, centers=num_clusters, nstart=10, iter.max=50)
        ## Take the sum, over each cluster, of the within-cluster sum of
        ## squares around the cluster mean. Then take the log. This is
        ## `W_k` in TWH's notation.
        log(sum(k$withinss))
    }
}

## For an appropriate reference distribution (in this case, uniform
## points in the same range as `data`), simulate the mean and standard
## deviation of the dispersion.
reference_dispersion <- function(data, num_clusters, num_reference_bootstraps) {
	dispersions <- sapply(1:num_reference_bootstraps, function(i) dispersion(generate_uniform_points(data), num_clusters))
	mean_dispersion <- mean(dispersions)
        ## the extra factor accounts for simulation error
	stddev_dispersion <- sd(dispersions) / sqrt(1 + 1 / num_reference_bootstraps)
	c(mean_dispersion, stddev_dispersion)
}

## Generate uniform points within the range of `data`.
generate_uniform_points <- function(data) {
    ## Find the min/max values in each dimension, so that we can
    ## generate uniform numbers in these ranges.
    ## row1: min; row2: max.
    rgs <- apply(data, 2, range)
    dims <- dim(data)
    ## For each dimension, generate `num_datapoints` points uniformly in
    ## the min/max range.
    sapply(1:dims[2], function(dim) runif(dims[1], min=rgs[1, dim], max=rgs[2, dim]))
}
