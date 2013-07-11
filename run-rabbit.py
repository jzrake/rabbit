#!/usr/bin/env python

from rabbit.mesh2d import *

def test1():
    mesh = RabbitMesh(2)
    mesh.add_volume(3, (0,2))
    mesh.add_volume(3, (1,2))
   
    n0 = mesh.get_volume(3, (1,2))
    n1 = mesh.containing_volume(4, (2, 4))

    assert n0 is n1

    mesh.add_volume(0, (0,0))
    mesh.containing_volume(1, (1,1))

    print mesh.containing_volume(1, (1,1))
    print n0.travel(0, (-1,0))

    print mesh.get_volume(0, (0,0)).morton()
    print n1.morton()
    assert mesh.add_volume(5, (31, 31)).morton() == (63, 63)


def test2():
    import matplotlib.pyplot as plt

    mesh = RabbitMesh(2)

    for i in range(2):
        for j in range(2):
            mesh.add_volume(2, (i, j))

    #mesh.add_volume(1, (0, 0))
    mesh.add_volume(1, (0, 1))
    mesh.add_volume(1, (1, 0))
    mesh.add_volume(1, (1, 1))

    mesh.create_vertices()
    Xv = [ ]
    Yv = [ ]
    Xn = [ ]
    Yn = [ ]

    for coords in mesh.vertices:
        Xv.append(coords[0])
        Yv.append(coords[1])

    for node in mesh.volumes.values():
        m = node.morton()
        Xn.append(m[0])
        Yn.append(m[1])

    dups = set()
    last_face = None

    for face in sorted(mesh.faces, key=compare_face):
        if last_face != None:
            if contains_face(last_face, face):
                dups.add(last_face)
        last_face = face

    mesh.faces -= dups

    print "there are %d unique faces and %d duplicates" % (
        len(mesh.faces), len(dups))

    for face in mesh.faces:
        if face[0][0] == 1 or True:
            plt.plot([face[1][0][0], face[1][1][0]],
                     [face[1][0][1], face[1][1][1]], c='k')

    plt.scatter(Xn, Yn, c='r')
    plt.scatter(Xv, Yv, c='b')
    plt.axis('equal')
    plt.show()

test1()
test2()
