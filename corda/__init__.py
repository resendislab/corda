#  __init__.py
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

from .util import reaction_confidence, test_model
from .corda import CORDA

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
