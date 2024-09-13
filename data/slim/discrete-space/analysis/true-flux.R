true_flux = function(true_node_states, ts, times, cost.mat, neighbor.mat,
    sample_sets, state_sets)
{
    storage.mode(cost.mat) = "double"
    adjacency_matrix = Matrix::Matrix(neighbor.mat, sparse = TRUE)
    adjacency_matrix = methods::as(adjacency_matrix, "generalMatrix")
    h = .Call(C_treeseq_discrete_mpr_edge_history, ts@treeseq, 
        true_node_states, cost.mat, adjacency_matrix, FALSE)
    H = data.frame(edge_id = rep.int(0:(length(h[[1]]) - 1L), 
        times = h[[3]]), state_id = do.call(c, h[[1]]), time = do.call(c, 
        h[[2]]))
    H = structure(H, node.state = true_node_states, path.offset = c(cumsum(h[[3]]) - 
        h[[3]], nrow(H)) + as.integer(FALSE))
    num_state_sets = as.numeric(max(state_sets))
    num_sample_sets = as.numeric(max(sample_sets))
    num_time_bins = length(times) - 1
    flux = .Call(C_treeseq_discrete_mpr_ancestry_flux, ts@tree, attr(H, 
            "path.offset"), H$state_id, H$time, as.integer(num_state_sets), 
            state_sets - 1L, as.integer(num_sample_sets), sample_sets - 
                1L, times)
    return (flux)
}