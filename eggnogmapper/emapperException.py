##
## CPCantalapiedra 2019

import sys

class EmapperException(Exception):
    def __init__(self, *args, **kargs):
        sys.excepthook = lambda exctype,exc,traceback: ""
        super(emapperException, self).__init__(*args, **kargs)

## END
