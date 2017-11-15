(_ROOT, _DEPTH, _BREADTH) = range(3)


class Node:
    def __init__(self, identifier, cellTypeMembers, const, proj, cellTypeMean):
        self.__identifier = identifier
        self.__children = []
        self.cellTypeMembers = cellTypeMembers
        self.const = const
        self.proj = proj
        self.cellTypeMean = cellTypeMean

    @property
    def identifier(self):
        return self.__identifier

    @property
    def children(self):
        return self.__children

    def add_child(self, identifier):
        self.__children.append(identifier)

    def setChildCellTypeMembers(self, childCellTypeMembers):
        self.childCellTypeMembers = childCellTypeMembers


class Tree:
    def __init__(self):
        self.__nodes = {}
        self.__head = None

    @property
    def nodes(self):
        return self.__nodes

    def head(self):
        return self.__head

    def add_head(self, identifier, cellTypeMembers, const, proj, cellTypeMean):
        node = Node(identifier, cellTypeMembers, const, proj, cellTypeMean)
        self[identifier] = node
        self.__head = node

        return node

    def add_node(self, identifier, cellTypeMembers, const, proj, cellTypeMean,
                 parent):
        node = Node(identifier, cellTypeMembers, const, proj, cellTypeMean)
        self[identifier] = node

        self[parent].add_child(identifier)

        return node

    def display(self, identifier, depth=_ROOT):
        children = self[identifier].children
        if depth == _ROOT:
            print("{0}".format(identifier))
        else:
            print("\t" * depth, "{0}".format(identifier))

        depth += 1
        for child in children:
            self.display(child, depth)  # recursive call

    def traverse(self, identifier, mode=_DEPTH):
        # Python generator. Loosly based on an algorithm from
        # 'Essential LISP' by John R. Anderson, Albert T. Corbett,
        # and Brian J. Reiser, page 239-241
        yield [
            "\ncellType Number: " + str(identifier),
            "cellType Const: " + str(self[identifier].const),
            "cellType Proj: " + str(self[identifier].proj),
            "cellType Mean: " + str(self[identifier].cellTypeMean),
            "Num cellType Members: " +
            str(len(self[identifier].cellTypeMembers)),
            "Num cellType Child Members: " +
            str(len(self[identifier].childCellTypeMembers))
        ]
        #yield self[identifier].cellTypeIterations
        queue = self[identifier].children
        while queue:
            yield [
                "\ncellType Number: " + str(queue[0]),
                "cellType Const: " + str(self[queue[0]].const),
                "cellType Proj: " + str(self[queue[0]].proj),
                "cellType Mean: " + str(self[queue[0]].cellTypeMean),
                "Num cellType Members: " +
                str(len(self[queue[0]].cellTypeMembers)),
                "Num cellType Child Members: " +
                str(len(self[queue[0]].childCellTypeMembers))
            ]
            #yield self[queue[0]].cellTypeIterations
            expansion = self[queue[0]].children
            if mode == _DEPTH:
                queue = expansion + queue[1:]  # depth-first
            elif mode == _BREADTH:
                queue = queue[1:] + expansion  # width-first

    def __getitem__(self, key):
        return self.__nodes[key]

    def __setitem__(self, key, item):
        self.__nodes[key] = item
