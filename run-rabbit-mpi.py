#!/usr/bin/env python

import os
import pickle
from rabbit.mesh2d_mpi import *

def build_mesh_parallel():
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


def read_mesh_serial():
    import glob
    mesh = load_meshes(glob.glob("mesh.*.pickle"))
    for node in mesh.volumes.values():
        node.data = node.host_proc
    plot_mesh(mesh, numbers=False, vertices=False, volume_args=dict(s=100))


os.system('rm -f *.pickle')
build_mesh_parallel()

if MPI.COMM_WORLD.rank == 0:
    read_mesh_serial()

MPI.COMM_WORLD.Barrier()
