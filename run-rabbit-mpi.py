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

    for rank in range(mesh.comm.size):
        if rank == mesh.comm.rank:
            print rank, mesh.node_label_range, len(mesh.volumes)
        mesh.comm.Barrier()

test1()
