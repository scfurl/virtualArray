virtualArrayMergeRecurse <-
function (dfs, by, ...)
{
    if (length(dfs) == 2) {
        merge(dfs[[1]], dfs[[2]], all = FALSE, sort = FALSE, ...)
    }
    else {
        merge(dfs[[1]], Recall(dfs[-1]), all = FALSE, sort = FALSE,
            ...)
    }
}

