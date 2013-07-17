
import numpy as np
from mpi4py import MPI

MAX_LABEL = 1<<32

class RabbitMesh(object):
    def __init__(self, rank, max_depth=5):
        self.volumes = { }
        self.rank = rank
        self.MAX_DEPTH = max_depth
        self.comm = MPI.COMM_WORLD
        self.node_label_range = [0, 0]
        if self.comm.rank == self.comm.size - 1:
            self.node_label_range[1] = MAX_LABEL

    def add_volume(self, depth, index, force=False, ghost=False):
        assert len(index) == self.rank
        assert depth <= self.MAX_DEPTH
        node = RabbitVolume(self, depth, index, ghost=ghost)
        if (self.node_label_range[0] <= node.preorder_label() and
            self.node_label_range[1] >  node.preorder_label()) or (
            force or ghost):
            self.volumes[(depth, index)] = node
            return node
        else:
            return None

    def del_volume(self, depth, index):
        assert len(index) == self.rank
        del self.volumes[(depth, index)]

    def get_volume(self, depth, index):
        assert len(index) == self.rank
        return self.volumes[(depth, index)]

    def containing_volume(self, depth, index):
        assert len(index) == self.rank
        if depth < 0:
            return None
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
            f = self.MAX_DEPTH - d + 1
            vert00 = ((v.index[0] + 0) << f, (v.index[1] + 0) << f)
            vert01 = ((v.index[0] + 0) << f, (v.index[1] + 1) << f)
            vert10 = ((v.index[0] + 1) << f, (v.index[1] + 0) << f)
            vert11 = ((v.index[0] + 1) << f, (v.index[1] + 1) << f)
            self.vertices |= set([vert00, vert01, vert10, vert11])
            self.faces |= set([RabbitFace(self, v, vert00, vert10),
                               RabbitFace(self, v, vert00, vert01),
                               RabbitFace(self, v, vert01, vert11),
                               RabbitFace(self, v, vert10, vert11)])
        dups = set()
        last_face = None
        for face in sorted(self.faces, key=compare_face):
            if last_face != None:
                if contains_face(last_face, face):
                    dups.add(last_face)
                    face.nodes |= last_face.nodes # face eats other face's nodes
            last_face = face
        self.faces -= dups
        print "there are %d unique faces and %d duplicates" % (
            len(self.faces), len(dups))

    def get_volumes(self, include_ghost=False, only_ghost=False):
        if only_ghost:
            return [n for n in self.volumes.values() if n.ghost]
        elif not include_ghost:
            return [n for n in self.volumes.values() if not n.ghost]
        else:
            return [n for n in self.volumes.values()]

    def global_node_count(self):
        recv = np.zeros(1, dtype=int)
        send = np.zeros(1, dtype=int)
        send[0] = len(self.get_volumes())
        self.comm.Allreduce(send, recv, op=MPI.SUM)
        return recv[0]

    def load_balance(self):
        """
        Re-determine the ownership of nodes with respect to processor. The
        global list of nodes is sorted according to its preorder label, then
        partitioned evenly, and finally redistributed.
        """
        # ----------------------------------------------------------------------
        # Determine the number of nodes currently on each processor
        # ----------------------------------------------------------------------
        proc_node_count = np.zeros(self.comm.size, dtype=int)
        proc_node_count[self.comm.rank] = len(self.get_volumes())
        self.comm.Allreduce(MPI.IN_PLACE, proc_node_count)

        # ----------------------------------------------------------------------
        # Assign a global index to each node, indicating its order in a
        # compressed global list of nodes
        # ----------------------------------------------------------------------
        node_index = proc_node_count.cumsum()[self.comm.rank - 1] if (
            self.comm.rank != 0) else 0

        ordered_nodes = sorted(self.get_volumes(),
                               key=lambda node: node.preorder_label())
        nodes_per_proc = self.global_node_count() / self.comm.size

        for node in ordered_nodes:
            node.global_compressed_index = node_index
            node_index += 1
            node.host_proc = node.global_compressed_index / nodes_per_proc
            if node.host_proc >= self.comm.size:
                node.host_proc = self.comm.size - 1

        # ----------------------------------------------------------------------
        # Send each node to its new host processor
        # ----------------------------------------------------------------------
        requests = [ ]

        for proc in range(self.comm.size):
            if proc == self.comm.rank: continue
            nodes_for_proc = [node.pack() for node in ordered_nodes if
                              node.host_proc == proc]
            requests.append(self.comm.isend(nodes_for_proc, proc))

        for proc in range(self.comm.size):
            if proc == self.comm.rank: continue
            nodes_from_proc = self.comm.recv(source=proc)
            for packed_node in nodes_from_proc:
                added_node = self.add_volume(packed_node['depth'],
                                             packed_node['index'], force=True)
                added_node.data = packed_node['data']
                added_node.host_proc = self.comm.rank

        for request in requests:
            request.Wait()

        # ----------------------------------------------------------------------
        # Remove sent nodes from this processor's mesh
        # ----------------------------------------------------------------------
        for key, node in self.volumes.items():
            if node.host_proc != self.comm.rank and not node.ghost:
                del self.volumes[key]

        # ----------------------------------------------------------------------
        # Re-calculate the node label range on all processors
        # ----------------------------------------------------------------------
        new_ordered_nodes = sorted(self.get_volumes(),
                                   key=lambda node: node.preorder_label())
        my_lower_label = 0 if self.comm.rank == 0 else min(
            [node.preorder_label() for node in self.get_volumes()])

        self.node_partition = np.zeros(self.comm.size, dtype=int)
        self.node_partition[self.comm.rank] = my_lower_label
        self.comm.Allreduce(MPI.IN_PLACE, self.node_partition, op=MPI.SUM)
        self.node_label_range[0] = my_lower_label
        self.node_label_range[1] = self.node_partition[self.comm.rank + 1] if (
            self.comm.rank != self.comm.size - 1) else MAX_LABEL

    def locate_ghost_nodes(self):
        for node in self.get_volumes():
            for offset in [[0, -1],
                           [0, +1],
                           [-1, 0],
                           [+1, 0]]:
                target_index = tuple([i + oi for i, oi in
                                      zip(node.index, offset)])
                if self.containing_volume(node.depth, target_index) == None:
                    self.add_volume(node.depth, target_index, ghost=True)

    def synchronize_ghost_nodes(self):
        for node in self.get_volumes(only_ghost=True):
            label = node.preorder_label()
            proc0 = 0
            proc1 = self.comm.size
            while proc1 - proc0 > 1:
                middle = (proc0 + proc1) / 2
                if self.node_partition[middle] < label:
                    proc0 = middle
                else:
                    proc1 = middle
            node.host_proc = proc0


class RabbitVolume(object):
    def __init__(self, mesh, depth, index, ghost=False):
        self.mesh = mesh
        self.depth = depth
        self.index = index
        self.ghost = ghost
        self.data = 0.0

    def travel(self, depth, index, containing=False):
        if depth == 0:
            target_index = self.index
        elif depth < 0:
            target_index = tuple([i >> +depth for i in self.index])
        elif depth > 0:
            target_index = tuple([i << -depth for i in self.index])
        I = tuple([ti + i for ti, i in zip(target_index, index)])
        if containing:
            return self.mesh.containing_volume(self.depth + depth, I)
        else:
            return self.mesh.get_volume(self.depth + depth, I)

    def coordinates(self, normalize=False):
        """
        Return the location of the center of the node's volume in units of leafs
        one below the maximum depth, unless normalize is True in which case node
        positions are mapped to the unit interval
        """
        d = self.depth
        D = self.mesh.MAX_DEPTH
        if normalize:
            return tuple([(2 * i + 1) / (1.0*(2 << d)) for i in self.index])
        else:
            return tuple([(2 * i + 1) << (D - d) for i in self.index])

    def preorder_label(self):
        """
        Return a unique integer label for this node corresponding to the order
        it would be visited in a pre-order traversal of the mesh tree, were it
        densely fleshed out
        """
        zorder = interleave_bits2(*self.index)
        return preorder_label2D(self.depth, zorder, self.mesh.MAX_DEPTH)

    def pack(self):
        """
        Return a dictionary with the necessary data members to rebuild this node
        on a different processor
        """
        return dict(depth=self.depth, index=self.index, data=self.data)

    def __repr__(self):
        return "<node: depth=%d index=%s>" % (self.depth, self.index)


class RabbitFace(object):
    def __init__(self, mesh, node, vertex0, vertex1):
        di = [b - a for a, b in zip(vertex0, vertex1)]
        assert di.count(0) == len(di) - 1
        axis, size = next(([n,d] for n,d in enumerate(di) if d != 0))
        self.mesh = mesh
        self.vertex0 = vertex0
        self.vertex1 = vertex1
        self.axis = axis
        self.size = size
        self.depth = ({1<<n: n for n in range(mesh.MAX_DEPTH + 2)})[self.size]
        self.nodes = set([node])


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


def preorder_label(depth, index, max_depth, m):
    """
    Return the order in which a given node is visited in a preorder traversal of
    a fully fleshed out tree having max_depth and branching ratio m
    """
    assert depth <= max_depth
    label = 0
    for d in range(depth):
        n = depth - d - 1 # bit
        h = max_depth - d - 1 # height
        nb = m**n # number of nodes at the target level below current
        sd = (index / nb) % m
        Md = tree_size(m, h)
        adding = sd * Md + 1
        label += adding
    return label


def preorder_label1D(depth, index, max_depth):
    return preorder_label(depth, index, max_depth, 2)


def preorder_label2D(depth, index, max_depth):
    return preorder_label(depth, index, max_depth, 4)


def preorder_label3D(depth, index, max_depth):
    return preorder_label(depth, index, max_depth, 8)


def tree_size(m, n):
    """ Return the size of a tree of max depth depth n and branching ratio m """
    return (m**(n+1) - 1) / (m - 1)


def interleave_bits2(a, b):
    """
    Create a 64-bit integer by interleaving the bits of the 32-bit integers a
    and b
    """
    assert a < (1 << 32)
    assert b < (1 << 32)
    base_a = [(a >> n) & 1 for n in range(32)]
    base_b = [(b >> n) & 1 for n in range(32)]
    res = [ ]
    for n in range(32):
        res.append(base_a[n])
        res.append(base_b[n])
    return sum([r << n for n,r in enumerate(res)])


def interleave_bits3(a, b, c):
    """
    Create a 64-bit integer by interleaving the bits of the 21-bit integers a,
    b, and c
    """
    assert a < (1 << 21)
    assert b < (1 << 21)
    assert c < (1 << 21)
    base_a = [(a >> n) & 1 for n in range(21)]
    base_b = [(b >> n) & 1 for n in range(21)]
    base_c = [(c >> n) & 1 for n in range(21)]
    res = [ ]
    for n in range(21):
        res.append(base_a[n])
        res.append(base_b[n])
        res.append(base_c[n])
    return sum([r << n for n,r in enumerate(res)])


def plot_mesh(mesh, numbers=True, vertices=True, faces=True,
              volume_args=dict(),
              ghost_args=dict()):
    import matplotlib.pyplot as plt
    Xv = [ ] # vertex coordinates
    Yv = [ ]
    Xn = [ ] # node (volume) coordinates
    Yn = [ ]
    Zn = [ ]
    Xg = [ ] # ghost coordinates
    Yg = [ ]
    Zg = [ ]
    mesh.create_vertices()

    if vertices:
        for coords in mesh.vertices:
            Xv.append(coords[0])
            Yv.append(coords[1])
        plt.scatter(Xv, Yv, c='b')

    if faces:
        for face in mesh.faces:
            plt.plot([face.vertex0[0], face.vertex1[0]],
                     [face.vertex0[1], face.vertex1[1]], c='k')

    for node in mesh.get_volumes(include_ghost=True):
        m = node.coordinates()
        if node.ghost:
            Xg.append(m[0])
            Yg.append(m[1])
            Zg.append(node.data)
        else:
            Xn.append(m[0])
            Yn.append(m[1])
            Zn.append(node.data)
        if numbers:
            plt.text(m[0]+0.1, m[1]+0.1, node.preorder_label())

    plt.scatter(Xn, Yn, c=Zn, **volume_args)
    plt.scatter(Xg, Yg, c=Zg, **ghost_args)
    plt.axis('equal')
    plt.show()


def load_meshes(fnames):
    """
    Load a mesh a collection of pickles, presumably written by different
    processors
    """
    import pickle
    meshes = [ ]

    for fname in fnames:
        infile = open(fname, 'r')
        new_mesh = pickle.load(infile)
        new_mesh.node_label_range = [0, MAX_LABEL]
        meshes.append(new_mesh)
        infile.close()

    for mesh in meshes:
        print "adding %d nodes from mesh" % len(mesh.volumes)
        for node in mesh.volumes.values():
            new_node = meshes[0].add_volume(node.depth, node.index)
            new_node.host_proc = node.host_proc

    meshes[0].comm = MPI.COMM_WORLD
    return meshes[0]
