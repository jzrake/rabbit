#!/usr/bin/env python

import os
import pickle
from rabbit.mesh2d_mpi import *


def build_mesh_parallel():
    mesh = RabbitMesh(2, max_depth=10)
    if True: # uniform-depth mesh
        for i in range(32):
            for j in range(32):
                mesh.add_volume(5, (i, j))
    else: # more interesting mesh
        for i in range(32):
            for j in range(32):
                if i >= 16 or j >= 16:
                    mesh.add_volume(5, (i, j))
        for i in range(8):
            for j in range(8):
                if i < 4 and j < 4:
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


def read_mesh_serial():
    import glob
    mesh = load_meshes(glob.glob("mesh.*.pickle"))
    for node in mesh.volumes.values():
        node.data = node.host_proc
    plot_mesh(mesh, numbers=False, vertices=False, volume_args=dict(s=100))


def test1():
    os.system('rm -f *.pickle')
    build_mesh_parallel()
    MPI.COMM_WORLD.Barrier()
    if MPI.COMM_WORLD.rank == 0:
        read_mesh_serial()
    MPI.COMM_WORLD.Barrier()


def test2():
    mesh = RabbitMesh(2)
    for i in range(8):
        for j in range(8):
            if i < 4 or j < 4:
                node = mesh.add_volume(3, (i, j))

    mesh.locate_ghost_nodes()
    plot_mesh(mesh, numbers=False, ghost_args=dict(marker='x'))


test1()
