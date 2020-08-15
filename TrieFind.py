class TrieNode:
    def __init__(self, lable='', level=0):
        self.label = lable
        self.childs = []
        self.level = level


    # ########################################## #
    #             searching section              #
    # ########################################## #

    # recursively will find a specific node coresponding to a specific kmer
    def find(self, kmer):
        if len(kmer) == 0:
            return self
        if len(self.childs) == 0:
            return None
        for child in self.childs:
            if child.label[-1] == kmer[0]:
                return child.find(kmer[1:])


    # ########################################## #
    #         saving seen kmers section          #
    # ########################################## #

    def add_frame(self, kmer):

        if len(kmer) == 0:
            # end of the path
            if not hasattr(self, 'count'):
                self.count = 0
            self.count += 1
            return
        
        # searching for proper path
        for child in self.childs:
            if child.label[-1] == kmer[0]:
                return child.add_frame(kmer[1:])

        # proper child (path) dose not exist
        new_child = TrieNode(self.label + kmer[0], self.level+1)
        self.childs += [new_child]
        return new_child.add_frame(kmer[1:])


    def max_count(self):
        self_count = 0
        if hasattr(self, 'count'):
            self_count = self.count
        return max([child.max_count() for child in self.childs] + [self_count])


    def min_count(self):
        self_count = 0
        if hasattr(self, 'count'):
            self_count = self.count
        return min([child.min_count() for child in self.childs] + [self_count])

    def max_min_count(self):
        self_count_min = 1000000000
        self_count_max = -1000000000
        if hasattr(self, 'count'):
            self_count_min = self.count
            self_count_max = self.count
        child_max_min = [child.max_min_count() for child in self.childs]
        return ( max([child_value[0] for child_value in child_max_min] + [self_count_max]) , min([child_value[1] for child_value in child_max_min] + [self_count_min]) )

    def scale_8bit(self, max_min_value):
        if hasattr(self, 'count'):
            return int( ( (self.count - max_min_value[1]) / (max_min_value[0] - max_min_value[1]) ) * 255 )
        return -1

    
    def report(self, max_min_value):
        report = ''
        
        if hasattr(self, 'count'):
            report += 'kmer=%s, count=%d, scale=%d\n'%(self.label, self.count, self.scale_8bit(max_min_value))
        
        for child in self.childs:
            report += child.report(max_min_value)

        return report