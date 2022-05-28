import matplotlib.pyplot as plt

import feems_utils

try:
    sample_set = snakemake.params['sample_set']
    sample_class = snakemake.params['sample_class']
    flag = snakemake.params['flag']
    out = snakemake.output[0]
except NameError:
    # testing
    sample_set = 'pendula'
    sample_class = 'default'
    flag = 'biallelic'
    out = 'scratch/feems_locations.svg'

sp_graph, samples, coords = feems_utils.get_sp_graph(sample_set, flag, sample_class=sample_class)

projection = feems_utils.get_projection(coords)

# draw figure showing which nodes the sample snapped to
v = feems_utils.create_plot(sp_graph, projection)
v.draw_samples()
v.draw_obs_nodes(use_ids=False)

plt.savefig(out, bbox_inches='tight', pad_inches=0)
