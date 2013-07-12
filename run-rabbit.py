#!/usr/bin/env python

from rabbit.mesh2d import *

def test1():
    mesh = RabbitMesh(2)
    mesh.add_volume(3, (0,2))
    mesh.add_volume(3, (1,2))
   
    n0 = mesh.get_volume(3, (1,2))
    n1 = mesh.containing_volume(4, (2, 4))

    mesh.add_volume(0, (0,0))
    mesh.containing_volume(1, (1,1))

    assert n0 is n1
    assert preorder_label1D(0, 0, max_depth=3) == 0
    assert preorder_label1D(1, 0, max_depth=3) == 1
    assert preorder_label1D(2, 0, max_depth=3) == 2
    assert preorder_label1D(3, 0, max_depth=3) == 3
    assert preorder_label1D(3, 1, max_depth=3) == 4
    assert preorder_label1D(2, 1, max_depth=3) == 5
    assert preorder_label1D(3, 2, max_depth=3) == 6
    assert preorder_label1D(3, 3, max_depth=3) == 7
    assert preorder_label1D(1, 1, max_depth=3) == 8
    assert preorder_label1D(2, 2, max_depth=3) == 9
    assert preorder_label1D(3, 4, max_depth=3) == 10
    assert preorder_label1D(3, 5, max_depth=3) == 11
    assert preorder_label1D(2, 3, max_depth=3) == 12
    assert preorder_label1D(3, 6, max_depth=3) == 13
    assert preorder_label1D(3, 7, max_depth=3) == 14

    assert [preorder_label2D(1, i, max_depth=1) for i in range(4)] == [1,2,3,4]
    assert [preorder_label2D(1, i, max_depth=2) for i in range(4)] == [1,6,11,16]
    assert [preorder_label2D(2, i, max_depth=2) for i in range(16)] == [
        2,3,4,5,7,8,9,10,12,13,14,15,17,18,19,20]


def test2():
    import matplotlib.pyplot as plt

    mesh = RabbitMesh(2, max_depth=2)
    mesh.add_volume(0, (0, 0))
    for i in range(2):
        for j in range(2):
            mesh.add_volume(1, (i, j))
    for i in range(4):
        for j in range(4):
            mesh.add_volume(2, (i, j))

    mesh.create_vertices()

    Xv = [ ]
    Yv = [ ]
    Xn = [ ]
    Yn = [ ]

    for coords in mesh.vertices:
        Xv.append(coords[0])
        Yv.append(coords[1])

    for node in mesh.volumes.values():
        m = node.coordinates()
        Xn.append(m[0])
        Yn.append(m[1])
        plt.text(m[0]+0.1, m[1]+0.1, node.preorder_label())

    for face in mesh.faces:
        plt.plot([face.vertex0[0], face.vertex1[0]],
                 [face.vertex0[1], face.vertex1[1]], c='k')

    plt.scatter(Xn, Yn, c='r')
    plt.scatter(Xv, Yv, c='b')
    plt.axis('equal')
    plt.show()

test1()
test2()
