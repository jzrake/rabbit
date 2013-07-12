
MAX_DEPTH = 5

class RabbitMesh(object):
    def __init__(self, rank):
        self.volumes = { }
        self.rank = rank

    def add_volume(self, depth, index):
        assert len(index) == self.rank
        self.volumes[(depth, index)] = RabbitVolume(self, depth, index)
        return self.volumes[(depth, index)]

    def del_volume(self, depth, index):
        assert len(index) == self.rank
        del self.volumes[(depth, index)]

    def get_volume(self, depth, index):
        assert len(index) == self.rank
        return self.volumes[(depth, index)]

    def containing_volume(self, depth, index):
        assert len(index) == self.rank
        if depth < 0:
            raise KeyError
        try:
            return self.get_volume(depth, index)
        except KeyError:
            return self.containing_volume(depth-1, tuple([i/2 for i in index]))

    def create_vertices(self):
        self.vertices = set()
        self.faces = set()
        for v in self.volumes.values():
            # for 2D only
            d = v.depth
            f = MAX_DEPTH - d + 1
            vert00 = ((v.index[0] + 0) << f, (v.index[1] + 0) << f)
            vert01 = ((v.index[0] + 0) << f, (v.index[1] + 1) << f)
            vert10 = ((v.index[0] + 1) << f, (v.index[1] + 0) << f)
            vert11 = ((v.index[0] + 1) << f, (v.index[1] + 1) << f)
            self.vertices |= set([vert00, vert01, vert10, vert11])
            self.faces |= set([RabbitFace(self, vert00, vert10),
                               RabbitFace(self, vert00, vert01),
                               RabbitFace(self, vert01, vert11),
                               RabbitFace(self, vert10, vert11)])
        dups = set()
        last_face = None
        for face in sorted(self.faces, key=compare_face):
            if last_face != None:
                if contains_face(last_face, face):
                    dups.add(last_face)
            last_face = face
        self.faces -= dups
        print "there are %d unique faces and %d duplicates" % (
            len(self.faces), len(dups))


class RabbitVolume(object):
    def __init__(self, mesh, depth, index):
        self.mesh = mesh
        self.depth = depth
        self.index = index

    def travel(self, depth, index):
        if depth == 0:
            target_index = self.index
        elif depth < 0:
            target_index = tuple([i >> +depth for i in self.index])
        elif depth > 0:
            target_index = tuple([i << -depth for i in self.index])
        I = tuple([ti + i for ti, i in zip(target_index, index)])
        return self.mesh.volumes[(self.depth + depth, I)]

    def morton(self):
        d = self.depth
        D = MAX_DEPTH
        return tuple([(2 * i + 1) << (D - d) for i in self.index])

    def __repr__(self):
        return "<RabbitVolume @ depth=%d index=%s>" % (self.depth, self.index)


class RabbitFace(object):
    def __init__(self, mesh, vertex0, vertex1):
        di = [b - a for a, b in zip(vertex0, vertex1)]
        assert di.count(0) == len(di) - 1
        axis, size = next(([n,d] for n,d in enumerate(di) if d != 0))
        self.mesh = mesh
        self.vertex0 = vertex0
        self.vertex1 = vertex1
        self.axis = axis
        self.size = size
        self.depth = ({1<<n: n for n in range(MAX_DEPTH + 1)})[self.size]


def compare_face(A):
    """
    Return a tuple which orders the faces depth-first
    """
    if A.axis == 0: other_axis = 1; my_axis = 0
    if A.axis == 1: other_axis = 0; my_axis = 1
    return (A.axis, # orientation (x, y, z) - directed
            A.vertex0[other_axis], # coordinate of other axis
            A.vertex0[my_axis], # left-end point
            A.depth) # face depth


def contains_face(A, B):
    """
    Return True if the face A contains the face B
    """
    if A.axis == 0: other_axis = 1; my_axis = 0
    if A.axis == 1: other_axis = 0; my_axis = 1
    if A.axis != B.axis: # different orientation?
        return False
    if A.vertex0[other_axis] != B.vertex0[other_axis]: # not co-linear?
        return False
    return (A.vertex0[my_axis] <= B.vertex0[my_axis] and
            A.vertex1[my_axis] >= B.vertex1[my_axis])


def preorder_label(depth, index, max_depth=MAX_DEPTH):
    """
    Return the order in which a given node is visited in a preorder traversal of
    a fully fleshed out tree having max_depth
    """
    label = 0
    for d in range(depth):
        if index & (1 << d) == 0:
            label += 1
        else:
            label += 2 << (d + max_depth - depth)
    return label
