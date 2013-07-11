#!/usr/bin/env python

MAX_DEPTH = 5


class RabbitMesh(object):
    volumes = { }
    def add_volume(self, depth, index):
        self.volumes[(depth, index)] = RabbitVolume(self, depth, index)
        return self.volumes[(depth, index)]

    def get_volume(self, depth, index):
        return self.volumes[(depth, index)]

    def containing_volume(self, depth, index):
        try:
            return self.get_volume(depth, index)
        except KeyError:
            return self.containing_volume(depth - 1, index / 2)


class RabbitVolume(object):
    def __init__(self, mesh, depth, index):
        self.mesh = mesh
        self.depth = depth
        self.index = index

    def travel(self, depth, index):
        if depth == 0:
            target_index = self.index
        elif depth < 0:
            target_index = self.index >> depth
        elif depth > 0:
            target_index = self.index << -depth
        return self.mesh.volumes[(self.depth + depth, target_index + index)]

    def morton(self):
        i = self.index
        d = self.depth
        return (2 * i + 1) * (1 << MAX_DEPTH) / (1 << d)

    def __repr__(self):
        return "<RabbitVolume @ depth=%d index=%d>" % (self.depth, self.index)


def test1():
    import matplotlib.pyplot as plt

    mesh = RabbitMesh()

    for depth in [0,1,2,3,4]:
        for i in range(1 << depth):
            if depth == 4 and i == 2: continue
            mesh.add_volume(depth, i)

    x = [v.morton() for v in mesh.volumes.values()]
    y = [v.depth for v in mesh.volumes.values()]

    node0 = mesh.get_volume(4, 4)
    node1 = node0.travel(0, -1)

    assert node1.index == 3
    assert node1.depth == 4
    print mesh.containing_volume(4, 2)

    plt.plot(x, y, 'o')
    plt.show()


test1()

