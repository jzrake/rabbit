#!/usr/bin/env python

from rabbit.mesh2d_mpi import *

def test1():
    mesh = RabbitMesh(2, max_depth=10)

    #mesh.load_balance()
    #print mesh.node_label_range

    for i in range(4):
        for j in range(4):
            mesh.add_volume(2, (i, j))

    mesh.load_balance()
    #print [node.host_proc for node in mesh.volumes.values()]

test1()
