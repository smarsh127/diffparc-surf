import pyvista as pv
from vtk import vtkNIFTIImageReader


class NIFTIReader(pv.BaseReader):
    _class_reader = vtkNIFTIImageReader


reader = NIFTIReader(snakemake.input.nii)
vol = reader.read()
print(vol)
surface = vol.contour([snakemake.params.threshold]).triangulate().strip()
print(surface)
surface.save(snakemake.output.vtk, binary=False)
