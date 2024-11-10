
class BaseFile(object):
    def __init__(self, path, mode):
        self._path = path
        self._mode = mode
        self._handle = None

    @property
    def path(self):
        return self._path

    @property
    def mode(self):
        return self._mode

    @property
    def handle(self):
        return self._handle

    def open(self):
        raise NotImplementedError(
            "This method should be overridden in subclasses!")

    def close(self):
        raise NotImplementedError(
            "This method should be overridden in subclasses!")
        
    def fetch(self):
        raise NotImplementedError(
            "This method should be overridden in subclasses!")
    
    def __iter__(self):
        for x in self.fetch():
            yield x

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()
        
