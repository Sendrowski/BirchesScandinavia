import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import pkg_resources
from feems import Viz, SpatialGraph
from feems.cross_validation import run_cv
from feems.utils import *
from shapely.geometry import MultiPoint, Polygon, Point

import utils


# get the SpatialGraph object for the given sample set
def get_sp_graph(sample_set, flag, sample_class, buffer=2):
    # fetch genotypes
    genotypes, samples, _, _ = utils.get_genotypes(sample_set, flag, sample_class)

    # remove monomorphic sites as they cause division by zero errors
    genotypes = utils.remove_monomorphic(genotypes)

    # setup graph
    coords = samples[['longitude', 'latitude']].values
    grid_path = pkg_resources.resource_filename("feems", "data/") + "grid_100.shp"

    # graph input files
    outer, edges, grid, _ = prepare_graph_inputs(coord=coords, ggrid=grid_path, translated=False, buffer=buffer)

    return SpatialGraph(genotypes, coords, grid, edges, scale_snps=True), samples, coords


# get projection with appropriate central coordinates
def get_projection(coords):
    return ccrs.EquidistantConic(central_latitude=np.median(coords[:, 0]), central_longitude=np.median(coords[:, 1]))


# refine the given set of tiles by filling each polygon with 4 polygons of equal shape
def refine_tiles(tiles):
    polygons = []
    for tile in tiles:
        edges = np.transpose(tile.exterior.xy)[:3]

        for i in range(0, 3):
            new_edges = [edges[i], (edges[i] + edges[(i + 1) % 3]) / 2, (edges[i] + edges[(i + 2) % 3]) / 2]

            polygons.append(Polygon((x, y) for x, y in new_edges))

        new_edges = [(edges[0] + edges[1]) / 2, (edges[1] + edges[2]) / 2, (edges[2] + edges[0]) / 2]

        polygons.append(Polygon((x, y) for x, y in new_edges))

    return polygons


# same function as for feems.utils but with a more refined grid
def prepare_graph_inputs(coord, ggrid, translated, buffer=0, outer=None):
    # no outer so construct with buffer
    if outer is None:
        points = MultiPoint([(x, y) for x, y in coord])
        xy = points.convex_hull.buffer(buffer).exterior.xy
        outer = np.array([xy[0].tolist(), xy[1].tolist()]).T

    if translated:
        outer[:, 0] = outer[:, 0] + 360.0

    tiles = load_tiles(ggrid)

    # intersect outer with discrete global grid
    bpoly = Polygon(outer)
    bpoly_translated = translate(bpoly, xoff=-360.0)
    tiles_restricted = [t for t in tiles if bpoly.intersects(t) or bpoly_translated.intersects(t)]

    # doubly refine tiles
    tiles_refined = refine_tiles(refine_tiles(tiles_restricted))
    pts, rev_pts, e = create_tile_dict(tiles_refined, bpoly)

    # construct grid array
    grid = []
    for i, v in rev_pts.items():
        grid.append((v[0], v[1]))
    grid = np.array(grid)

    assert grid.shape[0] != 0, "grid is empty changing translation"

    # un-translate
    if translated:
        pts = []
        for p in range(len(rev_pts)):
            pts.append(Point(rev_pts[p][0] - 360.0, rev_pts[p][1]))
        grid[:, 0] = grid[:, 0] - 360.0
        outer[:, 0] = outer[:, 0] - 360.0

    # construct edge array
    edges = np.array(list(e))
    ipmap = get_closest_point_to_sample(pts, coord)
    res = (outer, edges, grid, ipmap)
    return res


# visualize the given spatial graph object
def create_plot(sp_graph, projection, opts={}):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=projection)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1, hspace=0)

    opts = {**dict(projection=projection, edge_width=.75, edge_alpha=1, edge_zorder=100, sample_pt_size=50,
                   obs_node_size=6, sample_pt_color="#30D5C8", cbar_font_size=10, cbar_loc='center left',
                   cbar_orientation='vertical', cbar_height="20%", cbar_width="5%", cbar_ticklabelsize=10,
                   coastline_m="50m"), **opts}

    v = Viz(ax, sp_graph, **opts)

    v.draw_map()

    return v


# determine cross-validation error for given lambda
def determine_cv_error(lamb, sp_graph):
    # determine cross-validation error
    cv_err = run_cv(sp_graph, np.array([lamb]))

    # average over folds
    return np.mean(cv_err, axis=0)


# determine best lambda using cross-validation
def determine_best_lamb(lambs, sp_graph):
    # determine cross-validation error
    cv_err = run_cv(sp_graph, np.array(lambs))

    # average over folds
    mean_cv_err = np.mean(cv_err, axis=0)

    return plot_cv_errors(lambs, mean_cv_err)


# plot the cv error for different values of lambda
def plot_cv_errors(lambs, errors):
    # plot cv errors
    plt.plot(lambs, errors, marker='o')
    plt.xlabel("$\lambda$")
    plt.ylabel("L2 CV error")

    # mark lambda with lowest cv error
    best_lamb = lambs[np.argmin(errors)]
    plt.axvline(best_lamb, color="orange")

    plt.gca().set_xscale('log')
    plt.tight_layout()

    return float(best_lamb)
