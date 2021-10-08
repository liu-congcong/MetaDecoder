import numpy
from numpy.linalg import inv, slogdet
from scipy.special import logsumexp


class GMM:

    def __init__(self, components, X, log_prior, prior_weight, small_number = 1e-6, log_likelihood_tolerance = 1e-3, max_iterations = 1000):
        self.__components = components  # number of components
        self.__X = X  # list of samples, each sample is a #bins by #features matrix
        self.__log_prior = log_prior  # log probability
        self.__prior_weight = prior_weight
        self.__small_number = small_number
        self.__log_likelihood_tolerance = log_likelihood_tolerance
        self.__max_iterations = max_iterations


    def __get_x(self):
        '''
        self.__x: samples * features
        '''
        self.__x = numpy.array([numpy.mean(X, axis = 0) for X in self.__X], dtype = numpy.float64)


    def __initialize_log_gaussian_density(self):
        self.__log_gaussian_density = numpy.empty(shape = (self.__samples, self.__components), dtype = numpy.float64)


    def __construct_mapping(self):
        self.__mapping = numpy.array([index for index, X in enumerate(self.__X) for _ in X], dtype = numpy.int64)


    def __flatten_X(self):
        self.__X = numpy.concatenate(self.__X, axis = 0)


    def __initialize_responsibility(self):
        self.__responsibility = numpy.copy(self.__prior)
        self.log_responsibility = numpy.copy(self.__log_prior)


    def __update_mu(self):
        weight = self.__responsibility[self.__mapping]
        self.__mu = weight.T @ self.__X / (numpy.sum(weight, axis = 0, keepdims = True).T)


    def __update_sigma(self):
        self.__sigma = numpy.empty(shape = (self.__components, self.__features, self.__features), dtype = numpy.float64)
        weight = self.__responsibility[self.__mapping]
        for component in range(self.__components):
            temp = self.__X - self.__mu[component]
            self.__sigma[component] = temp.T * weight[:, component] @ temp / numpy.sum(weight[ : , component])
            self.__sigma[component][numpy.diag_indices(self.__features)] += self.__small_number


    def __update_pi(self):
        self.__log_pi = numpy.log(self.__responsibility + self.__prior_weight * self.__prior) - numpy.log(1.0 + self.__prior_weight)
        self.__pi = numpy.exp(self.__log_pi)


    def __calculate_log_gaussian_density(self):
        for component in range(self.__components):
            temp = self.__x - self.__mu[component]
            self.__log_gaussian_density[:, component] = (
                -0.5 * (
                    numpy.sum(temp @ inv(self.__sigma[component]) * temp, axis = 1)
                    + self.__features * numpy.log(2.0 * numpy.pi)
                    + slogdet(self.__sigma[component])[1]
                )
            )


    def __update_responsibility(self):
        self.__calculate_log_gaussian_density()
        # log(pi) + log(gaussian density)
        self.log_responsibility = self.__log_pi + self.__log_gaussian_density
        # log(posterior)
        self.log_responsibility -= logsumexp(self.log_responsibility, axis = 1, keepdims = True)
        # posterior
        self.__responsibility = numpy.exp(self.log_responsibility)
        self.__responsibility += numpy.finfo(numpy.float64).eps
        self.__responsibility /= numpy.sum(self.__responsibility, axis = 1, keepdims = True)


    def __calculate_log_likelihood(self):
        log_likelihood = (
            numpy.sum(logsumexp(self.__log_pi[self.__mapping] + self.__log_gaussian_density[self.__mapping], axis = 1))
            - self.__prior_weight * numpy.sum(self.__prior[self.__mapping] * (self.__log_prior[self.__mapping] - self.__log_pi[self.__mapping]))
        )
        return log_likelihood


    def main(self):
        self.__prior = numpy.exp(self.__log_prior)
        self.__get_x()  # self.__x
        self.__samples, self.__features = self.__x.shape
        self.__initialize_log_gaussian_density()  # self.__log_gaussian_density
        self.__construct_mapping()  # self.__mapping
        self.__flatten_X()
        self.__initialize_responsibility()  # self.__responsibility self.log_responsibility
        lower_log_likelihood = -numpy.inf
        iterations = 0
        while iterations < self.__max_iterations:
            self.__update_mu()  # m step
            self.__update_sigma()  # m step
            self.__update_pi()  # m step
            self.__update_responsibility()  # e step, update posterior distribution
            log_likelihood = self.__calculate_log_likelihood()  # log likelihood
            iterations += 1
            if (log_likelihood - lower_log_likelihood) < self.__log_likelihood_tolerance:
                break
            lower_log_likelihood = log_likelihood