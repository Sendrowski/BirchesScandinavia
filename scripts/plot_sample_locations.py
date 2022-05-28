import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from feems import viz
from pyproj import Proj

import feems_utils
import utils

try:
    test_mode = False
    sample_set = snakemake.params['sample_set']
    sample_class = snakemake.params['sample_class']

    latitudinal_barriers = snakemake.config['latitudinal_barriers']

    # get latitudinal barrier of given sample set from the config
    if sample_set in latitudinal_barriers:
        latitudinal_barrier = latitudinal_barriers[sample_set]
    else:
        latitudinal_barrier = None

    out = snakemake.output[0]
except NameError:
    # testing
    test_mode = True
    sample_set = 'pendula'
    sample_class = 'default'
    latitudinal_barrier = 64
    out = 'scratch/feems_locations.svg'

# fetch coordinates
coords = utils.get_samples(sample_set, sample_class)[['longitude', 'latitude']].values

# initialize figure and axis
fig = plt.figure()
projection = feems_utils.get_projection(coords)
ax = fig.add_subplot(1, 1, 1, projection=projection)
fig.subplots_adjust(left=0, right=1, bottom=0, top=1, hspace=0)
ax.axis("off")

# add cartopy features
ax.add_feature(cfeature.LAND, facecolor="#f7f7f7", zorder=0)
ax.coastlines("50m", color="#636363", linewidth=0.5, zorder=0)

# project coordinates
project = Proj(projection.proj4_init)
coords_projected = viz.project_coords(coords, project)

# draw coordinates
ax.scatter(
    coords_projected[:, 0],
    coords_projected[:, 1],
    edgecolors="black",
    linewidth=0.5,
    s=150,
    color="#30D5C8",
    marker=".",
    zorder=2
)

# draw line indicating subpopulation boundary
if latitudinal_barrier:
    p1 = project(min(coords[:, 0]), latitudinal_barrier)
    p2 = project(max(coords[:, 0]), latitudinal_barrier)
    ax.plot([p1[0], p2[0]], [p1[1], p2[1]], linestyle='dotted', c='grey', alpha=0.5, linewidth=3)

if test_mode:
    fig.show()

# save plot
plt.savefig(out, bbox_inches='tight', pad_inches=0)
