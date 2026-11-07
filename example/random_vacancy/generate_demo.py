"""Generate a tiny analytic velocity-wave example, not a material MD prediction."""
import argparse
import json
from pathlib import Path

import numpy as np


def generate(output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    frames = 128
    cell_ids = np.repeat(np.arange(4), 2)
    basis_ids = np.tile(np.arange(2), 4)
    masses = np.array([10.,20.])[basis_ids]
    origins = np.column_stack((np.arange(4), np.zeros((4,2))))
    positions = origins[cell_ids] + .25*basis_ids[:,None]
    keep = np.sort(np.random.RandomState(7).choice(8, 5, replace=False))
    phase = 2*np.pi*(np.arange(frames)[:,None]/frames - cell_ids[None,:]/4)
    phase = phase[:,:,None] + basis_ids[None,:,None]*.17 + np.array([0.,.37,-.2])
    amplitude = np.array([.001,.0007])[basis_ids][None,:,None]*np.array([1.,.5,.25])
    velocities = amplitude*np.cos(phase)  # Angstrom/fs
    trajectory = positions[None,:,:]+amplitude/(2*np.pi/frames)*np.sin(phase)
    rows = np.column_stack((np.arange(1,9),cell_ids+1,basis_ids+1,masses))
    np.savetxt(output/"perfect_basis.in",rows,fmt=["%d","%d","%d","%.12g"],
        header="Complete parent in row order\natom_id cell_id basis_id mass")
    np.savetxt(output/"kept_rows.txt",keep+1,fmt="%d")
    lattice = 'Lattice="4 0 0 0 1 0 0 0 1"'
    species = np.array(["C","Si"])[basis_ids]
    with (output/"perfect.xyz").open("w",encoding="utf-8") as stream:
        stream.write('8\n'+lattice+' Properties=species:S:1:pos:R:3\n')
        for element,position in zip(species,positions):
            stream.write(element+' '+' '.join("%.17g"%value for value in position)+'\n')
    with (output/"dump.xyz").open("w",encoding="utf-8") as stream:
        for t in range(frames):
            stream.write(str(len(keep))+'\n'+lattice+' Properties=species:S:1:pos:R:3:vel:R:3\n')
            for atom in keep:
                values = np.r_[trajectory[t,atom],velocities[t,atom]]
                stream.write(species[atom]+' '+' '.join("%.17g"%v for v in values)+'\n')
    common = """num_atoms = 5
total_num_steps = 128
time_step = 1
output_data_stride = 1
num_splits = 1
compress = 1
compress_block_size = 7
file_format = gpumd
dump_xyz_file = dump.xyz
output_hdf5 = trajectory.hdf5
basis_lattice_file = basis_vacancy.in
perfect_reference_file = perfect_reference.npz
allow_mass_disorder = 0
prim_unitcell = 1 0 0 0 1 0 0 0 1
supercell_dim = 4 1 1
rescale_prim = 0
num_qpaths = 1
q_path_name = GX
q_path = 0 0 0 0.5 0 0
use_parallel = 0
plot_SED = 0
lorentz = 0
"""
    for name,block,prefix in [("input_SED.in",0,"demo"),("input_stream.in",7,"stream")]:
        (output/name).write_text(common+"sed_projection_block_size = %d\nout_files_name = %s\n"%(block,prefix),encoding="utf-8")
    (output/"plot_SED.in").write_text(common.replace("plot_SED = 0","plot_SED = 1")+
        "out_files_name = demo\nplot_cutoff_freq = 30\nplot_interval = 10\nuse_contourf = 1\nif_show_figures = 0\n",encoding="utf-8")
    # Independent legacy-order calculation: FFT each atom, then sum by parent basis.
    qpoints = np.array([[0.,0.,0.],[.25,0.,0.],[.5,0.,0.]])
    atom_fft = np.fft.fft(velocities*100000,axis=0)
    expected = []
    for q in qpoints:
        power = np.zeros(frames)
        for b,mass in enumerate([10.,20.]):
            ids = keep[basis_ids[keep]==b]
            coherent = (atom_fft[:,ids,:]*np.exp(2j*np.pi*origins[cell_ids[ids]]@q)[None,:,None]).sum(axis=1)
            power += mass*1.66054e-27*(np.abs(coherent)**2).sum(axis=1)
        power /= 4*np.pi*128e-15*4
        power[1:64] = (power[1:64]+power[127:64:-1])/2
        expected.append(power[:64])
    np.savez_compressed(output/"expected_SED.npz",SED=np.column_stack(expected),
        qpoints=qpoints,frequency_THz=np.fft.fftfreq(frames,1e-15)[:64]/1e12)
    # Upstream Plot_SED uses the natural logarithm, not log10.
    maximum = int(np.ceil(np.log(np.max(expected))))
    with (output/"plot_SED.in").open("a",encoding="utf-8") as stream:
        stream.write("colorbar_min = %d\ncolorbar_max = %d\n"%(maximum-12,maximum))
    metadata = dict(interpretation="Analytic workflow fixture; labels/masses are placeholders, not a material prediction",
        frames=frames,parent_cells=4,parent_atoms=8,surviving_atoms=5,
        retained_parent_rows=(keep+1).tolist(),injected_frequency_THz=7.8125,
        injected_parent_q_fractional=[.25,0,0])
    (output/"expected.json").write_text(json.dumps(metadata,indent=2),encoding="utf-8")
    return metadata


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=Path("demo_run"))
    args = parser.parse_args()
    print(json.dumps(generate(args.output),indent=2))
