require(sf)
require(gaia)
source("proj.R")

data_full = function(chromosome, short=TRUE)
{
    if (short)
    {
        tsfile = sprintf("hgdp_tgp_sgdp_high_cov_ancients_chr%d_p.dated.trees",
            chromosome)
    }
    else
    {
        tsfile = sprintf("hgdp_tgp_sgdp_high_cov_ancients_chr%d_q.dated.trees",
            chromosome)
    }

    ts = treeseq_load(tsfile)
    samples = treeseq_individuals(ts)
    num_samples = nrow(samples)
    nodes = treeseq_nodes(ts)

    # keep georeferenced samples only
    {
        index = logical(num_samples)
        for (i in 1:num_samples)
        {
            s = samples[i, ]
            if (length(s[[2]]))
                index[i] = TRUE
        }
        filtered_samples = samples[index,]
        
        filtered_node_ids = nodes[
            nodes$individual_id %in% unlist(filtered_samples[,1]), "node_id"]

        ts_simpl = treeseq_simplify(ts, filtered_node_ids)

        samples = treeseq_individuals(ts_simpl)
        nodes = treeseq_nodes(ts_simpl)

        sample_ids = unlist(samples[,1])

        # lons,lats of georeferenced samples
        xy = t(apply(samples, 1, "[[", 2))[,2:1]

        # first pass, limit to african and eurasian samples
        filter1 = xy[,1] > -30 
        xy = xy[filter1, ]
        sample_ids = sample_ids[filter1]
        # finish by dropping australasian samples
        filter2 = !(xy[,1] > 110 & xy[,2] < 18) 
        xy = xy[filter2, ]
        sample_ids = sample_ids[filter2]
        sample_idx = match(sample_ids, nodes$individual_id)

        node_ids = nodes[
            c(rbind(
                  sample_idx
                , sample_idx+1L))
            , 1L
        ]

        ts_simpl2 = treeseq_simplify(ts_simpl, node_ids)

        # drop edges to nodes with insane number of parents or children
        edges = treeseq_edges(ts_simpl2)
        etab = table(edges$parent_id)
        ignore_parent = as.integer(names(which(etab > quantile(etab, 0.98))))
        etab = table(
            edges$child_id[edges$child_id >= treeseq_num_samples(ts_simpl2)])
        ignore_child = as.integer(names(which(etab > quantile(etab, 0.98))))

        ts_simpl3 = treeseq_drop_edges(ts_simpl2, ignore_parent, ignore_child)

        samples = treeseq_individuals(ts_simpl3)
        nodes = treeseq_nodes(ts_simpl3)

        sample_ids = unlist(samples[,1])

        xy = t(apply(samples, 1, "[[", 2))[,2:1]

        sample_idx = match(sample_ids, nodes$individual_id)

        node_ids = nodes[
            c(rbind(
                  sample_idx
                , sample_idx+1L))
            , 1L
        ]
    }

    coords = st_transform(st_as_sf(
        data.frame(
            node_id=node_ids,
            lon=rep(xy[,1],each=2),
            lat=rep(xy[,2],each=2)
        )
        , coords=c("lon", "lat")
        , crs=WGS84
    ), crs=st_crs(GRS80))

    landgrid = st_read("landgrid.gpkg", quiet=TRUE)

    sample_locations_hex = cbind(
        node_id=coords$node_id,
        state_id=st_nearest_feature(
            st_transform(coords, st_crs(landgrid))$geometry,
            st_centroid(landgrid$geom)
        )
    )

    sample_coords = cbind(node_id=node_ids,
        st_coordinates(st_transform(coords, crs=st_crs(EEGRS80))))

    list(
        ts=ts_simpl3,
        data=sample_locations_hex,
        sample.coords=sample_coords
    )
}


data_subsample = function(D)
{
    ts = D$ts
    dat = D$data
    nodes = treeseq_nodes(ts)
    indivs = tapply(nodes$individual_id, nodes$population_id, function(i) {
        indiv = unique(i)
        if (length(indiv) > 1)
            return (sample(indiv, 1L))
        else
            return (indiv)
    })
    if (indivs[1L] == -1L)
        indivs = indivs[-1L]
    nodes_to_keep = nodes$node_id[which(nodes$individual_id %in% indivs)]
    idx = which(dat[,1] %in% nodes_to_keep)
    ts2 = treeseq_simplify(ts, nodes_to_keep, node.map=TRUE)
    nodes = treeseq_nodes(ts2)
    node_ids = attr(ts2, 'node.map')[dat[idx, 1]+1]
    dat2 = cbind(node_id=node_ids, state_id=dat[idx, 2])
    coords = cbind(node_id=node_ids, D$sample.coords[idx, 2:3])
    list(
        ts=ts2,
        data=dat2,
        sample.coords=coords
    )
}
