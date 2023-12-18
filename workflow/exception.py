class BaseException(Exception):

    def __init__(self, message=''):
        Exception.__init__(self, message)


class GenomeNoMarkersFound(BaseException):

    def __init__(self, message=''):
        BaseException.__init__(self, message)
