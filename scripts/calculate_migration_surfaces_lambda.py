"""
Calculate migration surfaces with FEEMS for a certain lambda.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

from feems_utils import *

try:
    test_mode = False
    sample_set = snakemake.params['sample_set']
    sample_class = snakemake.params['sample_class']
    flag = snakemake.params['flag']
    lamb = snakemake.params.get('lamb', None)
    buffer = snakemake.params.get('buffer', 2)
    out_cv_error = snakemake.output.get('cv_error', None)
    out_surfaces = snakemake.output.surfaces
except NameError:
    # testing
    test_mode = True
    sample_set = 'pendula'
    sample_class = 'default'
    flag = 'biallelic'
    lamb = 10.0
    buffer = 2
    out_cv_error = 'scratch/cv_error.txt'
    out_surfaces = 'scratch/surfaces.png'

sp_graph, samples, coords = get_sp_graph(sample_set, flag, sample_class=sample_class, buffer=buffer)

projection = get_projection(coords)

if out_cv_error:
    # determine cross-validation error
    cv_error = determine_cv_error(lamb, sp_graph)[0, 0]

    with open(out_cv_error, 'w') as f:
        f.write(str(cv_error))

# fit the graph
sp_graph.fit(lamb=lamb)

# draw figure showing the migration surfaces
v = create_plot(sp_graph, projection)

v.draw_edge_colorbar()
v.draw_edges(use_weights=True)
v.draw_obs_nodes(use_ids=False)

plt.savefig(out_surfaces, bbox_inches='tight', pad_inches=0)

if test_mode:
    plt.show()
