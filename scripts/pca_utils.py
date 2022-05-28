import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression

import utils

# get the first two principal components of the specified sample set
def get(sample_set, flag, sample_class):
    # fetch genotypes
    genotypes, samples, _, _ = utils.get_genotypes(sample_set, flag, sample_class)

    pca = PCA(n_components=2)
    pc = pca.fit_transform(genotypes)

    return pca, pc, samples

# rotate the given points around specified origin by 'angle' radians
def rotate(origin, points, angle):
    x0, y0 = origin
    x, y = points[:, 0], points[:, 1]

    xr = x0 + np.cos(angle) * (x - x0) - np.sin(angle) * (y - y0)
    yr = y0 + np.sin(angle) * (x - x0) + np.cos(angle) * (y - y0)

    return xr, yr

# get an ellipse encircling the given set of points
def fit_ellipse(points: np.ndarray, buf: float, options={}) -> Ellipse:
    x, y = points[:, 0], points[:, 1]

    c = LinearRegression().fit([[z] for z in x], y).coef_[0]

    angle = 360 - np.rad2deg(np.arctan2(1, c))
    center = ((max(x) + min(x)) / 2, (max(y) + min(y)) / 2)

    x_rotated, y_rotated = rotate(center, points, - angle * 2 * np.pi / 360)

    width = max(x_rotated) - min(x_rotated) + buf / 2
    height = max(y_rotated) - min(y_rotated) + buf / 2

    return Ellipse(center, width, height, angle, **options)


# plot ellipses encircling the given sample set
def encircle_subset(sample_set, subsample_set, sample_class, data, color_index, buf=40):
    samples = utils.get_samples(sample_set, sample_class)
    subsamples = utils.get_samples(subsample_set, sample_class)
    mask = samples.name.isin(subsamples.name).values

    color = cm.viridis(color_index / 2)
    label = utils.get_proper_name(subsample_set)
    opts = dict(lw=1, facecolor=color[:3] + (0.1,), edgecolor=color, label=label)
    ellipse = fit_ellipse(np.array(data[mask]), buf=buf, options=opts)

    plt.gca().add_patch(ellipse)
