

class LogMessageFormatter(object):
    ''' 
    can be used to easily format log messages, as in:
    '''

    def __init__(self, *args, **kwargs):
        self.args = args

    def __str__(self):
        msg = ''
        for x in self.args:
            msg += ', %r' % x
        return msg
    
