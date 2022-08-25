from datetime import datetime
from math import ceil

import numpy
from numpy.linalg import inv, slogdet
from scipy.special import digamma, logsumexp
from sklearn.cluster import KMeans


class DPGMM:


    def __init__(self, components, min_component_weight, x, x_weight, random_number, max_iterations = 500, tolerance = 1e-3):
        self.max_components = components
        self.min_component_weight = min_component_weight
        self.x = x # pca(kmer frequency) + coverage #
        self.x_weight = x_weight
        self.random_number = random_number
        self.beta_prior_alpha = 1.0
        self.gaussian_prior_kappa = 1e-3
        self.max_iterations = max_iterations
        self.tolerance = tolerance


    def initialize_gaussian_prior_mu(self):
        self.gaussian_prior_mu = numpy.average(self.x, axis = 0, weights = self.x_weight)


    def initialize_wishart_prior_df(self):
        self.wishart_prior_df = self.features


    def initialize_wishart_prior_psi(self):
        x_minus_gaussian_prior_mu = self.x - self.gaussian_prior_mu
        self.wishart_prior_psi = x_minus_gaussian_prior_mu.T @ x_minus_gaussian_prior_mu / (self.samples - 1.)


    def initialize_r(self):
        self.r = numpy.zeros(shape = (self.samples, self.components), dtype = numpy.float64)
        kmeans = KMeans(n_clusters = min(self.components, self.samples), n_init = 2, random_state = self.random_number)
        kmeans_labels = kmeans.fit_predict((self.x - self.gaussian_prior_mu) / (numpy.sqrt(numpy.diag(self.wishart_prior_psi)) + numpy.finfo(numpy.float64).eps))
        for label in range(self.components):
            self.r[kmeans_labels == label, label] = 1.0


    def update_beta_gamma(self):
        self.beta_gamma = numpy.empty(shape = (self.components, 2), dtype = numpy.float64)
        self.beta_gamma[ : , 0] = 1. + self.component_weight
        self.beta_gamma[ : , 1] = self.beta_prior_alpha + numpy.hstack((numpy.cumsum(self.component_weight[ : 0 : -1])[ : : -1], 0))


    def update_gaussian_kappa(self):
        self.gaussian_kappa = self.gaussian_prior_kappa + self.component_weight


    def update_gaussian_mu(self):
        self.gaussian_mu = (self.gaussian_prior_kappa * self.gaussian_prior_mu + self.r.T * self.x_weight @ self.x) / self.gaussian_kappa.reshape(self.components, 1)


    def update_wishart_df(self):
        self.wishart_df = self.wishart_prior_df + self.component_weight


    def update_wishart_psi(self):
        self.wishart_psi = numpy.empty(shape = (self.components, self.features, self.features), dtype = numpy.float64)
        for component in range(self.components):
            component_mu_minus_gaussian_prior_mu = self.component_mu[component] - self.gaussian_prior_mu
            self.wishart_psi[component] = (
                self.wishart_prior_psi
                + self.component_weight[component] * (self.component_sigma[component] + self.gaussian_prior_kappa / self.gaussian_kappa[component] * numpy.outer(component_mu_minus_gaussian_prior_mu, component_mu_minus_gaussian_prior_mu))
            )
        self.wishart_psi_inverse = inv(self.wishart_psi)
        self.wishart_psi_log_determinant = slogdet(self.wishart_psi)[1]


    def update_r(self):
        digamma_gamma = digamma(self.beta_gamma) # digamma(γ1), digamma(γ2)
        self.r[ : ] = (
            # component_weight > #
            # = E[ln(v_t)] + E{i:1->t-1}[ln(1-v_i)] #
            # = digamma(γ_t,1) - digamma(γ_t,1 + γ_t,2) + ∑{i:1->t-1}(digamma(γ_i,2) - digamma(γ_i,1 + γ_i,2)) #
            # = digamma(γ_t,1) + ∑{i:1->t-1}(digamma(γ_i,2)) - ∑{i:1->t}(digamma(γ_i,1 + γ_i,2)) #
            digamma_gamma[ : , 0]
            + numpy.hstack((0, numpy.cumsum(digamma_gamma[ : -1, 1])))
            - numpy.cumsum(digamma(numpy.sum(self.beta_gamma, axis = 1)))
            # < component_weight #
            + 0.5 * (
                # ψd(x) = ∑{i: 1 -> d} ψ((df - (i - 1)) / 2) #
                + numpy.sum(digamma((self.wishart_df - numpy.arange(self.features).reshape(-1, 1)) / 2.), axis = 0)
                # d * ln2 #
                # + self.features * numpy.log(2.)
                # ln(|ψ|) #
                - self.wishart_psi_log_determinant
                # d / kappa #
                - self.features / self.gaussian_kappa
            )
        )
        for component in range(self.components):
            x_minus_gaussian_mu = self.x - self.gaussian_mu[component] # (self.samples, self.features) #
            self.r[: , component] += (
                - 0.5 * self.wishart_df[component] * numpy.sum(x_minus_gaussian_mu @ self.wishart_psi_inverse[component] * x_minus_gaussian_mu, axis = 1)
            )
        self.r = numpy.exp(self.r - logsumexp(self.r, axis = 1, keepdims = True))


    def shrink_components(self):
        components = numpy.sum(self.component_weight < self.min_component_weight)
        components = numpy.argsort(self.component_weight)[min(ceil(components * 0.9), self.components - 1) : ]
        self.r = self.r[ : , numpy.sort(components)]
        self.r /= (numpy.sum(self.r, axis = 1, keepdims = True) + numpy.finfo(numpy.float64).eps)
        self.components = components.shape[0]


    def update_component_weight(self):
        self.component_weight = self.x_weight @ self.r + numpy.finfo(numpy.float64).eps


    def update_component_mu(self):
        self.component_mu = (self.r.T * self.x_weight) @ self.x / self.component_weight.reshape(self.components, 1)


    def update_component_sigma(self):
        self.component_sigma = numpy.empty(shape = (self.components, self.features, self.features), dtype = numpy.float64)
        for component in range(self.components):
            x_minus_component_mu = self.x - self.component_mu[component] # (self.samples, self.features) #
            self.component_sigma[component] = (self.r[:, component] * self.x_weight * x_minus_component_mu.T) @ x_minus_component_mu / self.component_weight[component]
            self.component_sigma[component][numpy.diag_indices(self.features)] += 1e-6


    def is_convergence(self):
        return numpy.mean(numpy.sqrt(numpy.sum(numpy.square(self.r_ - self.r), axis = 1))) < self.tolerance


    def main(self):
        self.samples, self.features = self.x.shape
        self.initialize_gaussian_prior_mu() # mean of x
        self.initialize_wishart_prior_df() # features
        self.initialize_wishart_prior_psi() # cov of weighted x
        self.components = self.max_components # The current number of components is from "max_components" to "min_component".
        iterations = 0
        iterations_ = 0
        self.initialize_r()
        self.update_component_weight()
        while True:

            # m step #
            self.update_component_mu()
            self.update_component_sigma()
            self.update_beta_gamma()
            self.update_gaussian_kappa()
            self.update_gaussian_mu()
            self.update_wishart_df()
            self.update_wishart_psi()

            # e step
            self.r_ = numpy.array(self.r)
            self.update_r()
            self.update_component_weight()

            iterations += 1
            iterations_ += 1
            print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', 'DPGMM: {0} iterations have been finished ({1} components).'.format(iterations, self.components), flush = True)
            if (iterations == self.max_iterations) or (self.components == 1):
                break
            elif numpy.min(self.component_weight) >= self.min_component_weight:
                if self.is_convergence() or (iterations_ >= 100):
                    break
            elif (iterations >= 20) and (not (iterations % 10)):
                self.shrink_components()
                self.update_component_weight()
                iterations_ = 0