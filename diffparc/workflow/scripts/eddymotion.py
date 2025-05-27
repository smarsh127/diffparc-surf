import eddymotion.dmri as dmri
import eddymotion.estimator as estimator
import eddymotion.model as model
import numpy as np
import os

in_dwi = snakemake.input.dwi
in_b0 = snakemake.input.b0
in_bvec = snakemake.input.bvec
in_bval = snakemake.input.bval
in_brainmask = snakemake.input.brainmask

tmp_gradients_file = "temp_gradients.txt"  # this ends up in the shadow dir
out_affine_dir = snakemake.output.affine_dir

# convert bvec and bval into format eddymotion likes
bvec = np.loadtxt(in_bvec).T
bval = np.loadtxt(in_bval).reshape([bvec.shape[0], 1])

gradients = np.hstack([bvec, bval])

np.savetxt(tmp_gradients_file, gradients)

img = dmri.load(
    in_dwi,
    gradients_file=tmp_gradients_file,
    b0_file=in_b0,
    brainmask_file=in_brainmask,
)


est = estimator.EddyMotionEstimator()
est.fit(
    img,
    n_iter=snakemake.params.n_iter,
    model=snakemake.params.model,
    omp_nthreads=snakemake.threads,
)

os.mkdir(out_affine_dir)
for i, aff in enumerate(img.em_affines):
    aff.to_filename(
        os.path.join(out_affine_dir, f"gradvol-{i:03d}_desc-itk_affine.txt"), fmt="itk"
    )
