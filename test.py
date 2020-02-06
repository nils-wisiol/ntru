import doctest
from unittest import TextTestRunner

import ciphers
import polymod

if __name__ == '__main__':
    runner = TextTestRunner()
    for module in [
        ciphers,
        polymod
    ]:
        runner.run(doctest.DocTestSuite(module))
