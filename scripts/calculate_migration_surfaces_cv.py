from feems_utils import *

try:
    test_mode = False
    sample_set = snakemake.params['sample_set']
    sample_class = snakemake.params['sample_class']
    flag = snakemake.params['flag']
    buffer = snakemake.params.get('buffer', 2)
    lambdas = snakemake.params.lambdas
    out_cv = snakemake.output.cv
    out_surfaces = snakemake.output.surfaces
except NameError:
    # testing
    test_mode = True
    sample_set = 'pendula'
    sample_class = 'default'
    flag = 'biallelic'
    buffer = 2
    lambdas = np.geomspace(1e-8, 1e2, 11)
    out_cv = 'scratch/cv.png'
    out_surfaces = 'scratch/surfaces.png'

sp_graph, samples, coords = get_sp_graph(sample_set, flag, sample_class=sample_class, buffer=buffer)

projection = get_projection(coords)

# determine cross-validation error
best_lambda = determine_best_lamb(lambdas, sp_graph)

utils.save_fig(out_cv, tight_layout=True, show=test_mode, clear=True)

# fit the graph
sp_graph.fit(lamb=best_lambda)

# draw figure showing the migration surfaces
v = create_plot(sp_graph, projection)

v.draw_edge_colorbar()
v.draw_edges(use_weights=True)
v.draw_obs_nodes(use_ids=False)

utils.save_fig(out_surfaces, tight_layout=True, show=test_mode)
