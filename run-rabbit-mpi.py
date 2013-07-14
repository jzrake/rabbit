#!/usr/bin/env python

from rabbit.mesh2d_mpi import *

def test1():
    mesh = RabbitMesh(2, max_depth=10)

    #mesh.load_balance()
    print mesh.node_label_range

    for i in range(128):
        mesh.add_volume(7, (i, 0))

    mesh.load_balance()
    #print [node.host_proc for node in mesh.volumes.values()]

test1()
