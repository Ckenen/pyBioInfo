
class SyncWrapper(object):
    def __init__(self, index, container):
        self._index = index
        self._container = container
        self._iter = None
        self._end = False
        self._current = None
        
    @property
    def index(self):
        return self._index
    
    def __iter__(self):
        for item in self._container:
            yield item
            
    def _get_next_item(self):
        if self._iter is None:
            self._iter = self.__iter__()
        item = None
        if not self._end:
            try:
                item = next(self._iter)
            except StopIteration:
                self._end = True
        return item
    
    def get_current(self):
        if self._current is None and not self._end:
            self._current = self._get_next_item()
        return self._current
    
    def reset_current(self):
        if self._end:
            self._current = None
        else:
            self._current = self._get_next_item()

class SyncLoader(object):
    def __init__(self, containers):
        self.containers= containers
        self.wrappers = []
        self.valid_wrapper_indexes = []
        for i, container in enumerate(self.containers):
            self.wrappers.append(SyncWrapper(i, container))
            self.valid_wrapper_indexes.append(i)
    
    def __iter__(self):
        while True:
            temp_wrapper = None
            temp_wrapper_index = None
            for index in self.valid_wrapper_indexes:
                wrapper = self.wrappers[index]  
                if wrapper.get_current() is None:
                    self.valid_wrapper_indexes.remove(index)
                    continue
                if temp_wrapper is None:
                    temp_wrapper = wrapper
                    temp_wrapper_index = index
                    continue
                if wrapper.get_current() < temp_wrapper.get_current():
                    temp_wrapper = wrapper
                    temp_wrapper_index = index
            if temp_wrapper is None:
                break
            item = (temp_wrapper_index, temp_wrapper.get_current())
            yield item 
            temp_wrapper.reset_current()
            
