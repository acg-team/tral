# (C) 2014-2015 Elke Schaper

"""
    :synopsis: Configuration Singleton class, ensuring that the configuration
        is created in exactly one instance, even though it is called from
        multiple modules in arbitrary orders.

    .. moduleauthor:: Elke Schaper <elke.schaper@isb-sib.ch>
"""

import configobj
from tral.paths import config_file

P_CONFIG = config_file("config.ini")
P_SPEC = config_file("spec.ini")


class Singleton:

    """
    A non-thread-safe helper class to ease implementing singletons.
    This should be used as a decorator -- not a metaclass -- to the
    class that should be a singleton.

    The decorated class can define one `__init__` function that
    takes only the `self` argument. Other than that, there are
    no restrictions that apply to the decorated class.

    To get the singleton instance, use the `Instance` method. Trying
    to use `__call__` will result in a `TypeError` being raised.

    Limitations: The decorated class cannot be inherited from.

    The credit goes to http://stackoverflow.com/users/627005/paul-manta
    http://stackoverflow.com/questions/42558/python-and-the-singleton-pattern

    """

    def __init__(self, decorated):
        self._decorated = decorated

    def instance(self):
        """
        Returns the singleton instance. Upon its first call, it creates a
        new instance of the decorated class and calls its `__init__` method.
        On all subsequent calls, the already created instance is returned.

        """
        try:
            return self._instance
        except AttributeError:
            self._instance = self._decorated()
            return self._instance

    def __call__(self):
        raise TypeError('Singletons must be accessed through `Instance()`.')

    def __instancecheck__(self, inst):
        return isinstance(inst, self._decorated)


@Singleton
class Configuration:
    """
    A singleton configuration object containing a `ConfigObj` instance.

    Attributes:
        config (configobj.ConfigObj): The `ConfigObj` instance.
    """
    def __init__(self):
        self.config = configobj.ConfigObj(
            P_CONFIG,
            configspec=P_SPEC,
            stringify=True)
