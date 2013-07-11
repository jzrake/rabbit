
MAX_DEPTH = 12

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
            return self.containing_volume(depth - 1, tuple([i/2 for i in index]))

    def create_vertices(self):
        self.vertices = set()
        self.faces = set()
        for v in self.volumes.values():
            # for 2D only
            d = v.depth
            f = D -  MAX_DEPTH + 1
            vert00 = ((v.index[0] + 0) << f, (v.index[1] + 0) << f)
            vert01 = ((v.index[0] + 0) << f, (v.index[1] + 1) << f)
            vert10 = ((v.index[0] + 1) << f, (v.index[1] + 0) << f)
            vert11 = ((v.index[0] + 1) << f, (v.index[1] + 1) << f)
            self.vertices |= set([vert00, vert01, vert10, vert11])
            self.faces |= set([((d, 0), (vert00, vert10)),
                               ((d, 1), (vert00, vert01)),
                               ((d, 0), (vert01, vert11)),
                               ((d, 1), (vert10, vert11))])


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


def compare_face(a):
    if a[0][1] == 0: other_axis = 1; my_axis = 0
    if a[0][1] == 1: other_axis = 0; my_axis = 1
    return (a[0][1], # orientation (x, y, z) - directed
            a[1][0][other_axis], # coordinate of other axis
            a[1][0][my_axis], # left-end point
            a[0][0]) # face depth


def contains_face(A, B):
    """
    Returns True if the face A contains the face B
    """
    if A[0][1] == 0: other_axis = 1; my_axis = 0
    if A[0][1] == 1: other_axis = 0; my_axis = 1
    if A[0][1] != B[0][1]: # same orientation
        return False
    if A[1][0][other_axis] != B[1][0][other_axis]: # different line
        return False

    if A[0][1] == 0: my_axis = 0
    if A[0][1] == 1: my_axis = 1
    return (A[1][0][my_axis] <= B[1][0][my_axis] and
            A[1][1][my_axis] >= B[1][1][my_axis])
