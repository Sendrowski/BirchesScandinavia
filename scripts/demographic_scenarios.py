"""
Demographic scenarios for dadi.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

import re
import sys
from abc import abstractmethod
from typing import Dict

import dadi
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import scipy.stats

import sfs_utils
import utils

"""
https://dadi.readthedocs.io/en/latest/user-guide/specifying-a-model/

Units
T: Times in units of 2*N_e generations which runs forward in time from 0 to t_present
nu: population size as fraction of the initial population size N_e
m: migration rate in units of 2*N_e*m_f where m_f is fraction of 
    individuals that migrate to the other population
    
m01 is migration from pendula to pubescens or from south to north
m10 is migration from pubescens to pendula or from north to south
"""


def convert_name(name):
    return name.replace('nu', 'Ne').replace('T', 't')


# convert 'method_name' to 'MethodName'
def to_pascal_case(str):
    return ''.join(x.title() for x in str.split('_'))


# convert 'MethodName' to 'method_name'
def to_snake_case(str):
    return re.sub(r'(?<!^)(?=[A-Z])', '_', str).lower()


def get_name(name):
    return to_pascal_case(name) + 'Scenario'


# instantiate the scenario of given name
def create(name, params: dict):
    return get_class(name)(**params)


def get_class(name):
    return getattr(sys.modules[__name__], get_name(name))


# restore the population scenario of specified name
# with the best likelihood from the given files
def restore(scenario_name, p0, vcf, pops, generation_time, n_pops, N_e, mu, sample_set=None):
    n_proj = 20
    dd = dadi.Misc.make_data_dict_vcf(vcf, pops)
    fs = sfs_utils.get_spectrum_from_dd(dd, n_proj=n_proj, n_pops=n_pops)

    params = dict(fs=fs, dd=dd, n_pops=n_pops, n_proj=n_proj, n_iter=1,
                  N_e=N_e, mu=mu, generation_time=generation_time,
                  sample_set=sample_set)

    scenario = create(scenario_name, params)

    # assign initial value
    scenario.p0 = p0

    # Repeat the simulation with the appropriate initial value
    # for one iteration to re-instantiate the best scenario.
    # This is not ideal as it doesn't exactly produce the same
    # likelihood but it fully restores the scenario.
    scenario.simulate()

    return scenario


class DemographicScenario:
    # Lower and upper bounds for the parameters to be estimated.
    # To be determined at runtime.
    lower = None
    upper = None

    # the parameter's sampling distribution
    dists = None

    # time in years since last ice age
    # i.e. when the ice sheets started to recede
    t_last_ice_age = 16000

    # the parameters names
    params = []

    # array of the parameters to be fixed
    fixed_params = None

    grids = [40, 50, 60]

    # this chunk size yielded a few hundred chunks without
    # any empty entries which was the case for smaller values
    chunk_size_bs = 1000000

    # the standard lower parameter bounds
    lower_bounds_default = {
        'nu': 0.001,
        'T': 0,
        'm': 0,
        's': 0,
        'c': 0.001,
        'e': 0.01
    }

    # the standard upper parameter bounds
    upper_bounds_default = {
        'nu': 1000,
        'T': 1,
        'm': 2,
        's': 1,
        'c': 1000,
        'e': 100
    }

    # the type of distribution to use when sampling the parameters
    param_dists = {
        'nu': 'loguniform',
        'T': 'uniform',
        'm': 'uniform',
        's': 'uniform',
        'c': 'loguniform',
        'e': 'loguniform'
    }

    def __init__(self, fs, dd, n_proj, n_pops, N_e, mu, n_iter, generation_time, p0=None,
                 mode='random_hopping', n_boostraps=100, sample_set=None):
        self.n_proj = n_proj
        self.fs = fs
        self.dd = dd
        self.mode = mode
        self.n_boostraps = n_boostraps
        self.p0 = p0
        self.sample_set = sample_set

        self.N_e = N_e
        self.generation_time = generation_time
        self.mu = mu
        self.n_iter = n_iter

        self.n_pops = n_pops

        self.inferred = None
        self.fs_model = None
        self.lnl = None

        # θ=4*N_e*μ
        self.theta = None

        self.values = None

        self.simulate_fs_extrap = None

        self.stds = None
        self.GIM = None

        self.T_LGM = self.t_last_ice_age / self.infer_scale('T')

        # determine fixed params if not specified explicitly
        if not self.fixed_params:
            self.fixed_params = [None] * len(self.params)
            for key, value in self.get_fixed_params().items():
                if key in self.params:
                    self.fixed_params[self.params.index(key)] = value

        self.determine_default_values()

    # determine upper lower bounds and their sampling distributions
    def determine_default_values(self):
        self.lower = self.determine_values(self.lower, self.lower_bounds_default)
        self.upper = self.determine_values(self.upper, self.upper_bounds_default)
        self.dists = self.determine_values(self.dists, self.param_dists)

    # return a dict of parameters to be fixed
    def get_fixed_params(self):
        scenario = to_snake_case(type(self).__name__)
        fixed_params = {}

        # introduce a little magic to determine
        # the fixed parameters from the class name
        if 'no_migration' in scenario:
            fixed_params['m'] = 0

        if 'since_ice_age' in scenario:
            fixed_params['T'] = self.T_LGM

        return fixed_params

    # replace the non-specified parameters by their default value
    def determine_values(self, bounds, std: Dict):
        if bounds is None:
            bounds = []
            for param in self.params:
                type = self.get_param_type(param)

                # use standard bound
                bounds.append(std[type])

        # assume we have a list instead
        else:
            for i, bound in enumerate(bounds):
                if bound is None:
                    type = self.get_param_type(self.params[i])

                    # use standard bound
                    bounds[i] = std[type]

        return bounds

    # get a random sample of the ith parameter
    def sample_params(self, i):
        if self.dists[i] == 'loguniform':
            return scipy.stats.loguniform.rvs(self.lower[i], self.upper[i])
        elif self.dists[i] == 'uniform':
            return scipy.stats.uniform.rvs(self.lower[i], self.upper[i])
        else:
            raise NotImplementedError(f"Sampling distribution {self.dists[i]} not implemented.")

    @staticmethod
    def get_param_type(name):
        return re.sub('\d', '', name)

    # simulate phi
    @abstractmethod
    def simulate_phi(self, params, phi, xx):
        pass

    # simulate the spectrum being provided with the callback to simulate phi
    def simulate_fs(self, params, ns, pts):
        # obtain default grid
        xx = dadi.Numerics.default_grid(pts)

        # generate equilibrium density
        phi = dadi.PhiManip.phi_1D(xx)

        phi = self.simulate_phi(params, phi, xx)

        return dadi.Spectrum.from_phi(phi, ns, (xx,) * self.n_pops)

    # perform the simulations
    # determining the best fit and additional parameters
    def simulate(self):

        self.infer_best_fit()

        # get SFS scaling factor
        self.theta = dadi.Inference.optimal_sfs_scaling(self.fs_model, self.fs)

    # infer the best fit using the specified optimization mode
    def infer_best_fit(self):
        n_params = len(self.params)
        ns = self.fs.sample_sizes
        self.simulate_fs_extrap = dadi.Numerics.make_extrap_func(self.simulate_fs)

        # run basin hopping algorithm if specified
        # the method of randomly initializing points
        # that are then locally optimized seems
        # to yield better results
        if self.mode == 'basin_hopping':

            print('Performing global optimization using basin hopping.')

            raise BaseException('Basin hopping is deprecated.')

            verbose = 1
            flush_delay = 0.5
            multinom = True
            func_args = []
            func_kwargs = {}
            fixed_params = None
            ll_scale = 1

            args = (self.fs, self.simulate_fs_extrap, self.grids, self.lower, self.upper, verbose,
                    multinom, flush_delay, func_args, func_kwargs, fixed_params,
                    ll_scale)

            if self.p0 is None:
                self.p0 = [self.sample_params(i) for i in range(n_params)]

            res = scipy.optimize.basinhopping(dadi.Inference._object_func_log, self.p0,
                                              minimizer_kwargs={'args': args}, niter=self.n_iter)

            self.inferred = res.x

            self.fs_model = self.simulate_fs_extrap(self.inferred, ns, self.grids)
            self.lnl = dadi.Inference.ll_multinom(self.fs_model, self.fs)

        else:

            if self.mode == 'local_optimization':
                print(f'Performing local optimization using bgfs with '
                      f'{self.n_iter} runs of perturbed initial values.')

            elif self.mode == 'random_hopping':
                print(f'Performing local optimization using bgfs with '
                      f'{self.n_iter} runs of random initializations.')

            else:
                raise NotImplementedError('Specified optimization mode not supported.')

            inferred_vals = [None] * self.n_iter
            likelihoods = [None] * self.n_iter

            for i in range(self.n_iter):

                if i > 0 or self.p0 is None:
                    if self.mode == 'local_optimization':
                        self.p0 = dadi.Misc.perturb_params(self.p0, fold=1,
                                                           lower_bound=self.lower,
                                                           upper_bound=self.upper)
                    else:
                        self.p0 = [self.sample_params(i) for i in range(n_params)]

                # infer values using local optimization
                inferred_vals[i] = dadi.Inference.optimize_log(self.p0, self.fs, self.simulate_fs_extrap, self.grids,
                                                               fixed_params=self.fixed_params,
                                                               lower_bound=self.lower,
                                                               upper_bound=self.upper,
                                                               maxiter=15)

                # obtain likelihood
                fs_model = self.simulate_fs_extrap(inferred_vals[i], ns, self.grids)
                likelihoods[i] = dadi.Inference.ll_multinom(fs_model, self.fs)

                print(f"{i}:\t{likelihoods[i]},\t{str(inferred_vals[i])}", flush=True)

            # choose model with highest likelihood
            self.lnl = max(likelihoods)
            self.inferred = inferred_vals[np.argmax(likelihoods)]
            self.fs_model = self.simulate_fs_extrap(self.inferred, ns, self.grids)

    # obtain standard deviations and Godambe Information Matrix
    # we don't subsample the data so creating all bootstraps at once is fine
    def perform_uncertainty_analysis(self):
        chunks = dadi.Misc.fragment_data_dict(self.dd, self.chunk_size_bs)

        bootstraps = dadi.Misc.bootstraps_from_dd_chunks(chunks, self.n_boostraps,
                                                         sfs_utils.get_pop_ids(self.n_pops),
                                                         [self.n_proj] * self.n_pops)

        try:
            self.stds, self.GIM = dadi.Godambe.GIM_uncert(self.simulate_fs_extrap, self.grids,
                                                          bootstraps, self.inferred, self.fs,
                                                          log=False, multinom=True, return_GIM=True)
        except np.linalg.LinAlgError:
            print('Singular GIM')
            self.stds = np.full(len(self.inferred) + 1, None)

    # plot the trajectory and write to specified file
    @abstractmethod
    def plot_trajectory(self, opts={}):
        pass

    # Get the proper names of the populations.
    # This is only possible if the sample set was specifies.
    def get_proper_pop_names(self):
        return sfs_utils.get_proper_names(self.sample_set) if self.sample_set else None

    # plot the sfs and save to specified file
    def plot_sfs(self, scaling_1d=[1, 1]):
        # plot spectra
        if self.n_pops == 1:
            # plot observed and modelled SFS
            sfs_utils.plot_sfs([self.fs, self.fs_model * self.theta],
                               log_scale=False, labels=['observed', 'modeled'])

            diff = dadi.Inference.Anscombe_Poisson_residual(self.fs, self.fs_model * self.theta)
            diff_cum = np.round(np.mean(np.abs(diff.compressed())), 2)

            plt.plot(1, 0, alpha=0, label='$log(L)=$' + str(np.round(self.lnl, 2)))
            plt.plot(1, 0, alpha=0, label=f'$\Sigma |r_i|/n={diff_cum}$')

            utils.scale_cf(scaling_1d)

            plt.legend()

        elif self.n_pops == 2:

            names = self.get_proper_pop_names()

            dadi.Plotting.plot_2d_comp_multinom(self.fs_model, self.fs, show=False,
                                                vmin=0.001, adjust=True, pop_ids=names)

            plt.gcf().tight_layout(pad=1)

    def extract_subplot(self, n):
        self.plot_sfs()

        fig = plt.gcf()

        # remove unwanted axes
        for i, ax in enumerate(fig.get_axes()):
            if i not in n:
                fig.delaxes(ax)
            else:
                # remove title
                ax.set_title(None)

        # redraw plot
        plt.draw()
        plt.tight_layout()

    def plot_data_2d(self):
        sfs_utils.plot_sfs([self.fs], sample_set=self.sample_set)

    def plot_model_2d(self):
        sfs_utils.plot_sfs([self.fs_model * self.theta], sample_set=self.sample_set)

    def plot_residuals_2d(self):
        self.extract_subplot([3, 4])

        diff = dadi.Inference.Anscombe_Poisson_residual(self.fs, self.fs_model * self.theta)
        diff_cum = np.round(np.mean(np.abs(diff.compressed())), 2)

        text_lnl = '$log(L) =$' + str(np.round(self.lnl, 2))
        text_resid = f'$\Sigma |r_i|/n={diff_cum}$'

        # add cumulative residual
        plt.gca().text(0.6, 17.6, text_lnl + '\n' + text_resid, color='black',
                       bbox=dict(facecolor='none', edgecolor='#c7c7c7',
                                 boxstyle='round,pad=0.2'))

    def plot_residuals_1d(self):

        diff = dadi.Inference.Anscombe_Poisson_residual(self.fs, self.fs_model * self.theta)

        sfs_utils.plot_sfs([diff], log_scale=False)

    # create a dataframe of the inferred parameters
    def to_dataframe(self):
        df = pd.DataFrame(columns=['lnl', 'theta'] + self.params)

        df.loc['params'] = [self.lnl, self.theta] + self.inferred.tolist()

        if self.stds is not None:
            df.loc['std_gim'] = [None, self.stds[-1]] + self.stds[:-1].tolist()

        return df

    # infer the scale from the given parameter name
    def infer_scale(self, name):
        if name.startswith('T'):
            # time in years
            return 2 * self.N_e / self.generation_time
        if name.startswith('nu'):
            # size in effective number of individuals
            return self.N_e

        return 1

    @staticmethod
    def adjust_plot_scaled_pop_size_trajectory(opts={}):
        plt.xlabel('t')
        plt.ylabel('$\\nu$', rotation=0, labelpad=7)
        plt.margins(x=0)

        if 'scaling' not in opts:
            opts['scaling'] = [1, 2 / 3]

        utils.scale_cf(opts['scaling'])

        plt.tight_layout()


class ContinuousPopGrowthScenario(DemographicScenario):
    params = ['T', 'nu0', 'nu1']

    @staticmethod
    def get_nu(T, nu0, nu1):
        pass

    def simulate_phi(self, params, phi, xx):
        T, nu0, nu1 = params
        return dadi.Integration.one_pop(phi, xx, T, self.get_nu(T, nu0, nu1))

    # plot population size trajectory
    def plot_trajectory_from_params(self, opts={}):
        # determine inferred time parameter
        T = self.inferred[self.params.index('T')]

        t = np.arange(0, T, T / 100)
        y = self.get_nu(*self.inferred)(t)

        plt.plot(t, y, linewidth=2)

        self.adjust_plot_scaled_pop_size_trajectory(opts)

    def plot_trajectory(self, opts={}):
        self.plot_trajectory_from_params(opts)


class ContinuousPopSizeScenario(DemographicScenario):
    def get_nu(self):
        pass

    def simulate_phi(self, params, phi, xx):
        T = params[0]
        return dadi.Integration.one_pop(phi, xx, T, self.get_nu(*params))

    def plot_trajectory(self, opts={}):
        t = np.arange(0, self.inferred[0], self.inferred[0] / 10000)
        y = [self.get_nu(*self.inferred)(t_i) for t_i in t]

        plt.plot(t, y, linewidth=2)

        self.adjust_plot_scaled_pop_size_trajectory(opts)


# This model is nested and can be fixed to obtain a constant pop size
class OnePopSizeChangeScenario(ContinuousPopSizeScenario):
    params = ['T', 'nu0', 's', 'nu1']

    def get_nu(self, T, nu0, s, nu1):
        return lambda t: nu0 if (t / T) < s else nu0 * nu1


# This model is nested and can be fixed to obtain a one population size change
class TwoPopSizeChangeScenario(ContinuousPopSizeScenario):
    params = ['T', 'nu0', 's1', 'nu1', 's2', 'nu2']

    def get_nu(self, T, nu0, s1, nu1, s2, nu2):
        def nu(t):
            # fraction of time passed
            tf = t / T

            if tf < s1:
                return nu0
            else:
                # fraction of time passed since s1
                tf2 = tf - s1 / (1 - s1)

                return nu1 if tf2 < s2 else nu2

        return nu

        return lambda t: nu0 if () < s1 else nu0 * nu1


class TwoPopSizeChangeSinceIceAgeScenario(TwoPopSizeChangeScenario):
    pass


class OnePopSizeChangeSinceIceAgeScenario(OnePopSizeChangeScenario):
    pass


# we fix two parameters here to obtain a constant size
class ConstantPopSizeScenario(OnePopSizeChangeScenario):
    def get_fixed_params(self):
        fixed = super().get_fixed_params()
        fixed['nu1'] = 1
        fixed['s'] = 0
        return fixed


class ConstantPopSizeSinceIceAgeScenario(ConstantPopSizeScenario):
    pass


class LinPopGrowthScenario(ContinuousPopGrowthScenario):
    params = ['nu0', 'c', 'T']

    # if c = 1, we have a constant pop size
    def get_nu(self, nu0, c, T):
        return lambda t: nu0 * (1 - (t / T)) + (nu0 * c) * (t / T)

    def simulate_phi(self, params, phi, xx):
        nu0, c, T = params
        return dadi.Integration.one_pop(phi, xx, T, self.get_nu(nu0, c, T))


class LinPopGrowthSinceIceAgeScenario(LinPopGrowthScenario):
    pass


class ExpPopGrowthScenario(ContinuousPopGrowthScenario):
    params = ['nu0', 'e', 'T']

    # if e = 0, we have a constant pop size
    def get_nu(self, nu0, e, T):
        return lambda t: nu0 * e ** (t / T)

    def simulate_phi(self, params, phi, xx):
        nu0, e, T = params
        return dadi.Integration.one_pop(phi, xx, T, self.get_nu(nu0, e, T))


class ExpPopGrowthSinceIceAgeScenario(ExpPopGrowthScenario):
    pass


# abstract class for two-population scenarios
class TwoPopulationScenario(DemographicScenario):

    # simulate phi
    @abstractmethod
    def simulate_phi(self, params, phi, xx):
        pass

    def plot_trajectory(self, opts={}):
        pass

    # simulate the spectrum being provided with the callback to simulate phi
    def simulate_fs(self, params, ns, pts):
        # obtain default grid
        xx = dadi.Numerics.default_grid(pts)

        # generate equilibrium density
        phi = dadi.PhiManip.phi_1D(xx, nu=1)

        # split into two populations
        phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

        phi = self.simulate_phi(params, phi, xx)

        return dadi.Spectrum.from_phi(phi, ns, (xx,) * self.n_pops)

    # population split at time T in the past with fraction
    # s and s-1 going into population 1 and 2, respectively.
    # Both populations then experience exponential population
    # growth until they reach sizes nu1 and nu2 at time 0, respectively.
    @staticmethod
    def get_fractional_exp_pop_size_change(s, nu1, nu2, T):
        return [
            # t ϵ [0, T]
            lambda t: s * (nu1 / s) ** (t / T),
            lambda t: (1 - s) * (nu2 / (1 - s)) ** (t / T)
        ]

    @staticmethod
    def get_nu_one_pop_size_change(T, nu1, s, nu2):
        return lambda t: nu1 if (t / T) < s else nu1 * nu2


# population split at time T in the past with symmetric migration
class SplitSymmetricMigrationScenario(TwoPopulationScenario):
    params = ['T', 'nu1', 'nu2', 'm']

    def simulate_phi(self, params, phi, xx):
        T, nu1, nu2, m = params
        return dadi.Integration.two_pops(phi, xx, T, nu1, nu2, m12=m, m21=m)


# population split at the LGM with no migration
class SplitNoMigrationScenario(SplitSymmetricMigrationScenario):
    pass


# population split since last ice age with symmetric migration
class SplitSymmetricMigrationSinceIceAgeScenario(SplitSymmetricMigrationScenario):
    pass


# population split at the LGM with no migration
class SplitNoMigrationSinceIceAgeScenario(SplitSymmetricMigrationScenario):
    pass


# Continuous asymmetric migration where each population
# experiences a single population size change.
class SplitSymmetricMigrationOnePopSizeChangeScenario(TwoPopulationScenario):
    params = ['T', 'nu11', 'nu12', 's1', 's2', 'nu21', 'nu22', 'm']

    def simulate_phi(self, params, phi, xx):
        T, nu11, nu12, s1, s2, nu21, nu22, m = params

        nu1 = self.get_nu_one_pop_size_change(T, nu11, s1, nu21)
        nu2 = self.get_nu_one_pop_size_change(T, nu12, s2, nu22)

        return dadi.Integration.two_pops(phi, xx, T, nu1, nu2, m12=m, m21=m)


# Continuous asymmetric migration where each population
# experiences a single population size change.
class SplitSymmetricMigrationOnePopSizeChangeSinceIceAgeScenario(SplitSymmetricMigrationOnePopSizeChangeScenario):
    pass


class SplitNoMigrationOnePopSizeChangeScenario(SplitSymmetricMigrationOnePopSizeChangeScenario):
    pass


class SplitNoMigrationOnePopSizeChangeSinceIceAgeScenario(SplitSymmetricMigrationOnePopSizeChangeScenario):
    pass


# Continuous asymmetric migration where each population
# experiences a single population size change.
class SplitAsymmetricMigrationOnePopSizeChangeScenario(TwoPopulationScenario):
    params = ['T', 'nu11', 'nu12', 's1', 's2', 'nu21', 'nu22', 'm01', 'm10']

    def simulate_phi(self, params, phi, xx):
        T, nu11, nu12, s1, s2, nu21, nu22, m01, m10 = params

        nu1 = self.get_nu_one_pop_size_change(T, nu11, s1, nu21)
        nu2 = self.get_nu_one_pop_size_change(T, nu12, s2, nu22)

        return dadi.Integration.two_pops(phi, xx, T, nu1, nu2, m12=m01, m21=m10)


# Continuous asymmetric migration where each population
# experiences a single population size change.
class SplitAsymmetricMigrationOnePopSizeChangeSinceIceAgeScenario(SplitAsymmetricMigrationOnePopSizeChangeScenario):
    pass


# population split with symmetric migration and linear
# migration rate change since last ice age
class SplitVariableSymmetricMigrationSinceIceAgeScenario(TwoPopulationScenario):
    params = ['nu1', 'nu2', 'm0', 'm1']

    def simulate_phi(self, params, phi, xx):
        nu1, nu2, m0, m1 = params

        m = lambda t: m0 * (1 - t / self.T_LGM) + m1 * (t / self.T_LGM)

        return dadi.Integration.two_pops(phi, xx, self.T_LGM, nu1, nu2, m12=m, m21=m)


# population split at time T with asymmetric migration
class SplitAsymmetricMigrationScenario(TwoPopulationScenario):
    params = ['T', 'nu1', 'nu2', 'm01', 'm10']

    def simulate_phi(self, params, phi, xx):
        T, nu1, nu2, m01, m10 = params
        return dadi.Integration.two_pops(phi, xx, T, nu1, nu2, m12=m01, m21=m10)


# population split at time T with asymmetric migration since last ice age
class SplitAsymmetricMigrationSinceIceAgeScenario(SplitAsymmetricMigrationScenario):
    pass


# population split since last ice age with unidirectional migration from pop 2 into pop 1
class SplitUnidirectionalMigrationFromTwoToOneScenario(SplitAsymmetricMigrationScenario):
    def get_fixed_params(self):
        fixed = super().get_fixed_params()
        fixed['m10'] = 0
        return fixed


# population split since last ice age with unidirectional migration from pop 2 into pop 1
class SplitUnidirectionalMigrationFromTwoToOneSinceIceAgeScenario(SplitUnidirectionalMigrationFromTwoToOneScenario):
    pass


# population split since last ice age with unidirectional migration from pop 1 into pop 2
class SplitUnidirectionalMigrationFromOneToTwoScenario(SplitAsymmetricMigrationScenario):
    def get_fixed_params(self):
        fixed = super().get_fixed_params()
        fixed['m01'] = 0
        return fixed


# population split since last ice age with unidirectional migration from pop 1 into pop 2
class SplitUnidirectionalMigrationFromOneToTwoSinceIceAgeScenario(SplitUnidirectionalMigrationFromOneToTwoScenario):
    pass


# a fraction s of pop 1 e.g. pubescens will be derived from pop 2 e.g. pendula
class SplitAdmixtureFromTwoToOneSinceIceAgeScenario(TwoPopulationScenario):
    params = ['s', 'nu1', 'nu2']

    def simulate_phi(self, params, phi, xx):
        s, nu1, nu2 = params
        phi = dadi.PhiManip.phi_2D_admix_2_into_1(phi, s, xx, xx)
        return dadi.Integration.two_pops(phi, xx, self.T_LGM, nu1, nu2)


# a fraction s of pop 2 e.g. pendula will be derived from pop 1 e.g. pubescens
class SplitAdmixtureFromOneToTwoSinceIceAgeScenario(TwoPopulationScenario):
    params = ['s', 'nu1', 'nu2']

    def simulate_phi(self, params, phi, xx):
        s, nu1, nu2 = params
        phi = dadi.PhiManip.phi_2D_admix_1_into_2(phi, s, xx, xx)
        return dadi.Integration.two_pops(phi, xx, self.T_LGM, nu1, nu2)
