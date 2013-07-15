#!/usr/bin/env python

import pickle
from rabbit.mesh2d_mpi import *

def test1():
    mesh = RabbitMesh(2, max_depth=10)

    for i in range(4):
        for j in range(4):
            if i >= 2 or j >= 2:
                mesh.add_volume(2, (i, j))

    for i in range(4):
        for j in range(4):
            mesh.add_volume(3, (i, j))

    mesh.load_balance()

    for rank in range(mesh.comm.size):
        if rank == mesh.comm.rank:
            print rank, mesh.node_label_range, len(mesh.volumes)
        mesh.comm.Barrier()

    outf = open("mesh.%04d.pickle" % mesh.comm.rank, 'w')
    del mesh.comm
    pickle.dump(mesh, outf)
    outf.close()


def test2():
    mesh = load_meshes(["mesh.%04d.pickle"%i for i in [0,3]])
    for node in mesh.volumes.values():
        node.data = node.host_proc
    plot_mesh(mesh, numbers=False, vertices=False, scatter_args=dict(s=100))


test1()
#test2()
