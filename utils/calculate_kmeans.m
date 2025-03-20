function [idx, c, sumd, D] = calculate_kmeans(estimation_2D, cluster_num, distance, rep)
    % K-means clustering
    stream = RandStream('mlfg6331_64');  % Random number stream
    options = statset('UseParallel',1,'UseSubstreams',1,'Streams',stream);
    tic; % Start stopwatch timer

    [idx, c, sumd, D] = kmeans(estimation_2D, cluster_num, 'Options',options,'MaxIter',10000,'Display','final', 'Distance', distance,'Replicates',rep);

end
