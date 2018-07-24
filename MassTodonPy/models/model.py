class Model(object):
    """A container for storing results of fitting."""
    def fit(self, x, y, **kwds):
        """Fit the model."""
        raise NotImplementedError

    def refit(self, **kwds):
        """Refit the model."""
        raise NotImplementedError

    def __call__(self, x, y, *args, **kwds):
        """Predict the values at the new data points."""
        raise NotImplementedError

    def predict(self, x, y, *args, **kwds):
        """Predict the values at the new data points."""
        self(x, y, *args, **kwds)

    def coef(self):
        """Retrieve spline coefficient.

        Then again, why would you?
        """
        raise NotImplementedError

    def res(self):
        """Get residuals."""
        raise NotImplementedError

    def errors(self):
        return self.res()

    def residuals(self):
        """Get residuals: syntactic sugar for 'res'."""
        return self.res

    def __repr__(self):
        return 'This is the logic for all models.'

    def plot(self, **kwds):
        """Plot results."""
        raise NotImplementedError

    def cv(self, x, y, folds):
        """Run cross-validation."""
        raise NotImplementedError


def predict(model, x, y, *args, **kwds):
    return model(x, y, *args, **kwds)

def fitted(model):
    return model.fitted()

def coef(model):
    return model.coef()

def coefficients(model):
    return model.coef()

def residuals(model):
    return model.res()

def res(model):
    return model.res()

def errors(model):
    return model.res()

def cv(model, **kwds):
    return model.cv(**kwds)

def plot(model, **kwds):
    model.plot(**kwds)