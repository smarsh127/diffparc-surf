from snakebids import bids


wildcard_constraints:
    shell="[0-9]+",
    dir="[a-zA-Z0-9]+",


if (config["in_prepdwi_dir"] == False) and (config["in_snakedwi_dir"] == False):

    rule import_dwi:
        input:
            nii=[
                re.sub(".nii.gz", ext, input_path["dwi"])
                for ext in [".nii.gz", ".bval", ".bvec", ".json"]
            ],
        output:
            nii=multiext(
                bids(
                    root=root,
                    suffix="dwi",
                    desc="import",
                    datatype="dwi",
                    **input_wildcards["dwi"]
                ),
                ".nii.gz",
                ".bval",
                ".bvec",
                ".json",
            ),
        group:
            "subj"
        run:
            for in_file, out_file in zip(input, output):
                shell("cp -v {in_file} {out_file}")




rule dwidenoise:
    input:
        multiext(
            bids(
                root=root,
                suffix="dwi",
                desc="import",
                datatype="dwi",
                **input_wildcards["dwi"]
            ),
            ".nii.gz",
            ".bvec",
            ".bval",
            ".json",
        ),
    output:
        multiext(
            bids(
                root=root,
                suffix="dwi",
                desc="denoise",
                datatype="dwi",
                **input_wildcards["dwi"]
            ),
            ".nii.gz",
            ".bvec",
            ".bval",
            ".json",
        ),
    container:
        config["singularity"]["diffparc"]
    log:
        bids(root="logs", suffix="denoise.log", **input_wildcards["dwi"]),
    group:
        "subj"
    shell:
        "dwidenoise {input[0]} {output[0]} 2> {log} && "
        "cp {input[1]} {output[1]} && "
        "cp {input[2]} {output[2]} && "
        "cp {input[3]} {output[3]}"


def get_degibbs_inputs(wildcards):
    # if input dwi at least 30 dirs, then grab denoised as input
    # else grab without denoising
    import numpy as np

    in_dwi_bval = re.sub(".nii.gz", ".bval", input_path["dwi"].format(**wildcards))
    bvals = np.loadtxt(in_dwi_bval)
    if bvals.size < 30:
        prefix = bids(
            root=root, suffix="dwi", desc="import", datatype="dwi", **wildcards
        )
    else:
        prefix = bids(
            root=root, suffix="dwi", datatype="dwi", desc="denoise", **wildcards
        )
    return multiext(prefix, ".nii.gz", ".bvec", ".bval", ".json")


rule mrdegibbs:
    input:
        get_degibbs_inputs,
    output:
        multiext(
            bids(
                root=root,
                suffix="dwi",
                datatype="dwi",
                desc="degibbs",
                **input_wildcards["dwi"]
            ),
            ".nii.gz",
            ".bvec",
            ".bval",
            ".json",
        ),
    container:
        config["singularity"]["diffparc"]
    #    log: bids(root='logs',suffix='degibbs.log',**config['input_wildcards']['dwi'])
    group:
        "subj"
    shell:
        "mrdegibbs {input[0]} {output[0]} && "
        "cp {input[1]} {output[1]} && "
        "cp {input[2]} {output[2]} && "
        "cp {input[3]} {output[3]}"
        #2> {log} && ' 


# now have nii with just the b0's, want to create the topup phase-encoding text files for each one:
rule get_phase_encode_txt:
    input:
        bzero_nii=bids(
            root=root,
            suffix="b0.nii.gz",
            datatype="dwi",
            desc="degibbs",
            **input_wildcards["dwi"]
        ),
        json=bids(
            root=root,
            suffix="dwi.json",
            datatype="dwi",
            desc="degibbs",
            **input_wildcards["dwi"]
        ),
    output:
        phenc_txt=bids(
            root=root,
            suffix="phenc.txt",
            datatype="dwi",
            desc="degibbs",
            **input_wildcards["dwi"]
        ),
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    script:
        "../scripts/get_phase_encode_txt.py"


rule concat_phase_encode_txt:
    input:
        phenc_txts=lambda wildcards: expand(
            bids(
                root=root,
                suffix="phenc.txt",
                datatype="dwi",
                desc="degibbs",
                **input_wildcards["dwi"]
            ),
            zip,
            **snakebids.filter_list(input_zip_lists["dwi"], wildcards)
        ),
    output:
        phenc_concat=bids(
            root=root,
            suffix="phenc.txt",
            datatype="dwi",
            desc="degibbs",
            **subj_wildcards
        ),
    group:
        "subj"
    shell:
        "cat {input} > {output}"


# function for either concatenating (if multiple inputs) or copying
def get_concat_or_cp_cmd(wildcards, input, output):
    if len(input) > 1:
        cmd = f"mrcat {input} {output}"
    elif len(input) == 1:
        cmd = f"cp {input} {output}"
    else:
        # no inputs
        cmd = None
    return cmd


rule concat_bzeros:
    input:
        bzero_niis=lambda wildcards: expand(
            bids(
                root=root,
                suffix="b0.nii.gz",
                datatype="dwi",
                desc="degibbs",
                **input_wildcards["dwi"]
            ),
            zip,
            **snakebids.filter_list(input_zip_lists["dwi"], wildcards)
        ),
    params:
        cmd=get_concat_or_cp_cmd,
    output:
        bzero_concat=bids(
            root=root,
            suffix="concatb0.nii.gz",
            datatype="dwi",
            desc="degibbs",
            **subj_wildcards
        ),
    container:
        config["singularity"]["diffparc"]
    log:
        bids(root="logs", suffix="concat_bzeros.log", **subj_wildcards),
    group:
        "subj"
    shell:
        "{params.cmd} 2> {log}"


# this rule should only run if there are multiple images
rule run_topup:
    input:
        bzero_concat=bids(
            root=root,
            suffix="concatb0.nii.gz",
            datatype="dwi",
            desc="degibbs",
            **subj_wildcards
        ),
        phenc_concat=bids(
            root=root,
            suffix="phenc.txt",
            datatype="dwi",
            desc="degibbs",
            **subj_wildcards
        ),
    params:
        out_prefix=bids(root=root, suffix="topup", datatype="dwi", **subj_wildcards),
        config="b02b0.cnf",  #this config sets the multi-res schedule and other params..
    output:
        bzero_corrected=bids(
            root=root,
            suffix="concatb0.nii.gz",
            desc="topup",
            datatype="dwi",
            **subj_wildcards
        ),
        fieldmap=bids(
            root=root,
            suffix="fmap.nii.gz",
            desc="topup",
            datatype="dwi",
            **subj_wildcards
        ),
        topup_fieldcoef=bids(
            root=root,
            suffix="topup_fieldcoef.nii.gz",
            datatype="dwi",
            **subj_wildcards
        ),
        topup_movpar=bids(
            root=root, suffix="topup_movpar.txt", datatype="dwi", **subj_wildcards
        ),
    container:
        config["singularity"]["diffparc"]  #fsl
    log:
        bids(root="logs", suffix="topup.log", **subj_wildcards),
    group:
        "subj"
    shell:
        "topup --imain={input.bzero_concat} --datain={input.phenc_concat} --config={params.config}"
        " --out={params.out_prefix} --iout={output.bzero_corrected} --fout={output.fieldmap} -v 2> {log}"


# this is for equal positive and negative blipped data - method=lsr --unused currently (jac method can be applied more generally)
rule apply_topup_lsr:
    input:
        dwi_niis=lambda wildcards: expand(
            bids(
                root=root,
                suffix="dwi.nii.gz",
                desc="degibbs",
                datatype="dwi",
                **input_wildcards["dwi"]
            ),
            zip,
            **snakebids.filter_list(input_zip_lists["dwi"], wildcards)
        ),
        phenc_concat=bids(
            root=root,
            suffix="phenc.txt",
            desc="degibbs",
            datatype="dwi",
            **subj_wildcards
        ),
        topup_fieldcoef=bids(
            root=root,
            suffix="topup_fieldcoef.nii.gz",
            datatype="dwi",
            **subj_wildcards
        ),
        topup_movpar=bids(
            root=root, suffix="topup_movpar.txt", datatype="dwi", **subj_wildcards
        ),
    params:
        #create comma-seperated list of dwi nii
        imain=lambda wildcards, input: ",".join(input.dwi_niis),
        # create comma-sep list of indices 1-N
        inindex=lambda wildcards, input: ",".join(
            [str(i) for i in range(1, len(input.dwi_niis) + 1)]
        ),
        topup_prefix=bids(root=root, suffix="topup", datatype="dwi", **subj_wildcards),
        out_prefix="dwi_topup",
    output:
        dwi_topup=bids(
            root=root,
            suffix="dwi.nii.gz",
            desc="topup",
            method="lsr",
            datatype="dwi",
            **subj_wildcards
        ),
    container:
        config["singularity"]["diffparc"]  #fsl
    shadow:
        "minimal"
    group:
        "subj"
    shell:
        "applytopup --verbose --datain={input.phenc_concat} --imain={params.imain} --inindex={params.inindex} "
        " -t {params.topup_prefix} -o {params.out_prefix} && "
        " fslmaths {params.out_prefix}.nii.gz {output.dwi_topup}"


def get_applytopup_inindex(wildcards):

    # need to get index of the scan from the subject scans
    # so first filter the ziplist to get the subject(+session)
    subj_filter = {"subject": wildcards.subject}
    if "session" in subj_wildcards.keys():
        subj_filter["session"] = wildcards.session

    zip_list_subj = snakebids.filter_list(input_zip_lists["dwi"], subj_filter)

    # now filter the subj ziplist using all wildcards to get the index of the scan
    indices = snakebids.filter_list(zip_list_subj, wildcards, return_indices_only=True)
    return indices[0] + 1  # get first index, but adjust to start at 1 instead of 0


rule apply_topup_jac:
    input:
        nii=bids(
            root=root,
            suffix="dwi.nii.gz",
            desc="degibbs",
            datatype="dwi",
            **input_wildcards["dwi"]
        ),
        phenc_scan=bids(
            root=root,
            suffix="phenc.txt",
            datatype="dwi",
            desc="degibbs",
            **input_wildcards["dwi"]
        ),
        phenc_concat=bids(
            root=root,
            suffix="phenc.txt",
            desc="degibbs",
            datatype="dwi",
            **subj_wildcards
        ),
        topup_fieldcoef=bids(
            root=root,
            suffix="topup_fieldcoef.nii.gz",
            datatype="dwi",
            **subj_wildcards
        ),
        topup_movpar=bids(
            root=root, suffix="topup_movpar.txt", datatype="dwi", **subj_wildcards
        ),
    params:
        inindex=get_applytopup_inindex,
        topup_prefix=bids(root=root, suffix="topup", datatype="dwi", **subj_wildcards),
    output:
        nii=bids(
            root=root,
            suffix="dwi.nii.gz",
            desc="topup",
            method="jac",
            datatype="dwi",
            **input_wildcards["dwi"]
        ),
    container:
        config["singularity"]["diffparc"]  #fsl
    shadow:
        "minimal"
    group:
        "subj"
    shell:
        " applytopup --verbose --datain={input.phenc_concat} --imain={input.nii} --inindex={params.inindex} "
        " -t {params.topup_prefix} -o dwi_topup --method=jac && mv dwi_topup.nii.gz {output.nii}"


# topup-corrected data is only used for brainmasking..
# here, use the jac method by default (later can decide if lsr approach can be used based on headers)
# with jac approach, the jac images need to be concatenated, then avgshell extracted

"""
rule cp_sidecars_topup_lsr:
    #TODO: BEST WAY TO TO EXEMPLAR DWI? 
    input: multiext(bids(root='work',suffix='dwi',desc='degibbs',datatype='dwi',**config['subj_wildcards'],**dwi_exemplar_dict),\
                '.bvec','.bval','.json')
    output: multiext(bids(root='work',suffix='dwi',desc='topup',method='lsr',datatype='dwi',**config['subj_wildcards']),\
                '.bvec','.bval','.json')
    run:
        for in_file,out_file in zip(input,output):
            shell('cp -v {in_file} {out_file}')
"""


rule cp_sidecars_topup_jac:
    input:
        multiext(
            bids(
                root=root,
                suffix="dwi",
                desc="degibbs",
                datatype="dwi",
                **subj_wildcards
            ),
            ".bvec",
            ".bval",
            ".json",
        ),
    output:
        multiext(
            bids(
                root=root,
                suffix="dwi",
                desc="topup",
                method="jac",
                datatype="dwi",
                **subj_wildcards
            ),
            ".bvec",
            ".bval",
            ".json",
        ),
    group:
        "subj"
    run:
        for in_file, out_file in zip(input, output):
            shell("cp -v {in_file} {out_file}")


rule concat_dwi_topup_jac:
    input:
        dwi_niis=lambda wildcards: expand(
            bids(
                root=root,
                suffix="dwi.nii.gz",
                desc="topup",
                method="jac",
                datatype="dwi",
                **input_wildcards["dwi"]
            ),
            zip,
            **snakebids.filter_list(input_zip_lists["dwi"], wildcards)
        ),
    output:
        dwi_concat=bids(
            root=root,
            suffix="dwi.nii.gz",
            desc="topup",
            method="jac",
            datatype="dwi",
            **subj_wildcards
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        "mrcat {input} {output}"


rule get_eddy_index_txt:
    input:
        dwi_niis=lambda wildcards: expand(
            bids(
                root=root,
                suffix="dwi.nii.gz",
                desc="degibbs",
                datatype="dwi",
                **input_wildcards["dwi"]
            ),
            zip,
            **snakebids.filter_list(input_zip_lists["dwi"], wildcards)
        ),
    output:
        eddy_index_txt=bids(
            root=root,
            suffix="dwi.eddy_index.txt",
            desc="degibbs",
            datatype="dwi",
            **subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    script:
        "../scripts/get_eddy_index_txt.py"


rule concat_degibbs_dwi:
    input:
        dwi_niis=lambda wildcards: expand(
            bids(
                root=root,
                suffix="dwi.nii.gz",
                desc="degibbs",
                datatype="dwi",
                **input_wildcards["dwi"]
            ),
            zip,
            **snakebids.filter_list(input_zip_lists["dwi"], wildcards)
        ),
    params:
        cmd=get_concat_or_cp_cmd,
    output:
        dwi_concat=bids(
            root=root,
            suffix="dwi.nii.gz",
            desc="degibbs",
            datatype="dwi",
            **subj_wildcards
        ),
    container:
        config["singularity"]["diffparc"]
    log:
        bids(root="logs", suffix="concat_degibbs_dwi.log", **subj_wildcards),
    group:
        "subj"
    shell:
        "{params.cmd} 2> {log}"


rule concat_input_dwi:
    """ rule for concatenating input dwi (when skipping dwi preproc) """
    input:
        dwi_niis=lambda wildcards: expand(
            bids(
                root=root,
                suffix="dwi.nii.gz",
                desc="import",
                datatype="dwi",
                **input_wildcards["dwi"]
            ),
            zip,
            **snakebids.filter_list(input_zip_lists["dwi"], wildcards)
        ),
    params:
        cmd=get_concat_or_cp_cmd,
    output:
        dwi_concat=bids(
            root=root,
            suffix="dwi.nii.gz",
            desc="import",
            datatype="dwi",
            **subj_wildcards
        ),
    container:
        config["singularity"]["diffparc"]
    log:
        bids(root="logs", suffix="concat_input_dwi.log", **subj_wildcards),
    group:
        "subj"
    shell:
        "{params.cmd} 2> {log}"


rule concat_runs_bvec:
    input:
        lambda wildcards: expand(
            bids(
                root=root,
                suffix="dwi.bvec",
                desc="{{desc}}",
                datatype="dwi",
                **input_wildcards["dwi"]
            ),
            zip,
            **snakebids.filter_list(input_zip_lists["dwi"], wildcards)
        ),
    output:
        bids(
            root=root,
            suffix="dwi.bvec",
            desc="{desc}",
            datatype="dwi",
            **subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    script:
        "../scripts/concat_bv.py"


rule concat_runs_bval:
    input:
        lambda wildcards: expand(
            bids(
                root=root,
                suffix="dwi.bval",
                desc="{{desc}}",
                datatype="dwi",
                **input_wildcards["dwi"]
            ),
            zip,
            **snakebids.filter_list(input_zip_lists["dwi"], wildcards)
        ),
    output:
        bids(
            root=root,
            suffix="dwi.bval",
            desc="{desc}",
            datatype="dwi",
            **subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    script:
        "../scripts/concat_bv.py"


# combines json files from multiple scans -- for now as a hack just copying first json over..
rule concat_runs_json:
    input:
        lambda wildcards: expand(
            bids(
                root=root,
                suffix="dwi.json",
                desc="{{desc}}",
                datatype="dwi",
                **input_wildcards["dwi"]
            ),
            zip,
            **snakebids.filter_list(input_zip_lists["dwi"], wildcards)
        ),
    output:
        bids(
            root=root,
            suffix="dwi.json",
            desc="{desc}",
            datatype="dwi",
            **subj_wildcards
        ),
    group:
        "subj"
    shell:
        "cp {input[0]} {output}"


#    script: '../scripts/concat_json.py'


rule get_shells_from_bvals:
    input:
        "{dwi_prefix}.bval",
    output:
        "{dwi_prefix}.shells.json",
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    script:
        "../scripts/get_shells_from_bvals.py"


# writes 4d file
rule get_shell_avgs:
    input:
        dwi="{dwi_prefix}.nii.gz",
        shells="{dwi_prefix}.shells.json",
    output:
        avgshells="{dwi_prefix}.avgshells.nii.gz",
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    script:
        "../scripts/get_shell_avgs.py"


# this gets a particular shell (can use to get b0)
rule get_shell_avg:
    input:
        dwi="{dwi_prefix}_dwi.nii.gz",
        shells="{dwi_prefix}_dwi.shells.json",
    params:
        bval="{shell}",
    output:
        avgshell="{dwi_prefix}_b{shell}.nii.gz",
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    script:
        "../scripts/get_shell_avg.py"


def get_dwi_mask():
    if config["masking_method"] == "b0_BET":
        method = "bet_from-b0"
    elif config["masking_method"] == "b0_SyN":
        method = f"b0SyN_from-{config['template']}"

    # then get bids name of file
    return bids(
        root=root,
        suffix="mask.nii.gz",
        desc="brain",
        method=method,
        datatype="dwi",
        **subj_wildcards,
    )


# have multiple brainmasking workflows -- this rule picks the method chosen in the config file
# generate qc snapshot for brain  mask
rule qc_brainmask_for_eddy:
    input:
        img=bids(
            root=root,
            suffix="b0.nii.gz",
            desc="preproc",
            datatype="dwi",
            **subj_wildcards
        ),
        seg=get_dwi_mask(),
    output:
        #        png = bids(root='qc',subject='{subject}',suffix='mask.png',desc='brain'),
        png=report(
            bids(root="qc", suffix="mask.png", desc="brain", **subj_wildcards),
            caption="../report/brainmask_dwi.rst",
            category="Brainmask",
        ),
        html=bids(root="qc", suffix="mask.html", desc="brain", **subj_wildcards),
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    script:
        "../scripts/vis_qc_dseg.py"


# if the --slspec_txt option is not used, use the SliceTiming json, otherwise just get from the file supplied at command-line
if not config["slspec_txt"]:

    rule get_slspec_txt:
        input:
            dwi_jsons=lambda wildcards: expand(
                bids(
                    root=root,
                    suffix="dwi.json",
                    desc="degibbs",
                    datatype="dwi",
                    **input_wildcards["dwi"]
                ),
                zip,
                **snakebids.filter_list(input_zip_lists["dwi"], wildcards)
            ),
        output:
            eddy_slspec_txt=bids(
                root=root,
                suffix="dwi.eddy_slspec.txt",
                desc="degibbs",
                datatype="dwi",
                **subj_wildcards
            ),
        group:
            "subj"
        container:
            config["singularity"]["diffparc"]
        script:
            "../scripts/get_slspec_txt.py"


else:

    rule get_slspec_txt:
        input:
            config["slspec_txt"],
        output:
            eddy_slspec_txt=bids(
                root=root,
                suffix="dwi.eddy_slspec.txt",
                desc="degibbs",
                datatype="dwi",
                **subj_wildcards
            ),
        group:
            "subj"
        shell:
            "cp {input} {output}"
# --- this is where the choice of susceptibility distortion correction (SDC) is made --



#
# get dwi reference for masking, registration -- either topup or degibbs b0
#  currently, if more than one dwi, we use topup, otherwise, no distortion correction (i.e. get degibbs)

#  TODO: implement other fieldmap based correction, and registration-based correction here too


def get_dwi_ref(wildcards):

    # this gets the number of DWI scans for this subject(session)
    filtered = snakebids.filter_list(input_zip_lists["dwi"], wildcards)
    num_scans = len(filtered["subject"])

    if num_scans > 1 and config["use_topup"]:
        return bids(
            root=root,
            suffix="b0.nii.gz",
            desc="topup",
            method="jac",
            datatype="dwi",
            **subj_wildcards,
        )
    else:
        if config["use_eddy"]:
            return bids(
                root=root,
                suffix="b0.nii.gz",
                desc="degibbs",
                datatype="dwi",
                **subj_wildcards,
            )
        else:
            return bids(
                root=root,
                suffix="b0.nii.gz",
                desc="preproc",
                datatype="dwi",
                **subj_wildcards,
            )


def get_eddy_topup_input(wildcards):
    # this gets the number of DWI scans for this subject(session)
    filtered = snakebids.filter_list(input_zip_lists["dwi"], wildcards)
    num_scans = len(filtered["subject"])

    if num_scans > 1 and config["use_topup"]:
        topup_inputs = {
            filename: bids(
                root=root,
                suffix=f"{filename}.nii.gz",
                datatype="dwi",
                **subj_wildcards,
            ).format(**wildcards)
            for filename in ["topup_fieldcoef", "topup_movpar"]
        }
        return topup_inputs
    else:
        return None


def get_eddy_topup_opt(wildcards, input):

    # this gets the number of DWI scans for this subject(session)
    filtered = snakebids.filter_list(input_zip_lists["dwi"], wildcards)
    num_scans = len(filtered["subject"])

    if num_scans > 1 and config["use_topup"]:
        topup_prefix = bids(
            root=root, suffix="topup", datatype="dwi", **subj_wildcards
        ).format(**wildcards)
        return f"--topup={topup_prefix}"
    else:
        return ""


def get_eddy_s2v_opts(wildcards, input):
    options = []
    if config["eddy_no_s2v"]:
        options += [
            f"--{key}={value}"
            for (key, value) in config["eddy"]["without_s2v"].items()
            if value is not None
        ]
    else:
        options += [
            f"--{key}={value}"
            for (key, value) in config["eddy"]["with_s2v"].items()
            if value is not None
        ]
    return " ".join(options)


def get_eddy_slspec_input(wildcards):
    if config["eddy_no_s2v"]:
        return {}
    else:
        return {
            "eddy_slspec_txt": bids(
                root=root,
                suffix="dwi.eddy_slspec.txt",
                desc="degibbs",
                datatype="dwi",
                **subj_wildcards,
            )
        }


def get_eddy_slspec_opt(wildcards, input):
    if config["eddy_no_s2v"]:
        return ""
    else:
        return f"--slspec={input.eddy_slspec_txt}"


def get_eddy_cmd(wildcards):
    if config.get("use_gpu_eddy_container", False):
        return (
            f"singularity exec --nv -e {config['singularity']['diffparc']} eddy_cuda9.1"
        )
    else:
        return "eddy_openmp"


rule run_eddy:
    input:
        unpack(get_eddy_slspec_input),
        dwi_concat=bids(
            root=root,
            suffix="dwi.nii.gz",
            desc="degibbs",
            datatype="dwi",
            **subj_wildcards
        ),
        phenc_concat=bids(
            root=root,
            suffix="phenc.txt",
            desc="degibbs",
            datatype="dwi",
            **subj_wildcards
        ),
        eddy_index_txt=bids(
            root=root,
            suffix="dwi.eddy_index.txt",
            desc="degibbs",
            datatype="dwi",
            **subj_wildcards
        ),
        brainmask=get_dwi_mask(),
        bvals=bids(
            root=root,
            suffix="dwi.bval",
            desc="degibbs",
            datatype="dwi",
            **subj_wildcards
        ),
        bvecs=bids(
            root=root,
            suffix="dwi.bvec",
            desc="degibbs",
            datatype="dwi",
            **subj_wildcards
        ),
    params:
        #set eddy output prefix to 'dwi' inside the output folder
        out_prefix=lambda wildcards, output: os.path.join(output.out_folder, "dwi"),
        flags=" ".join(
            [
                f"--{key}"
                for (key, value) in config["eddy"]["flags"].items()
                if value == True
            ]
        ),
        container=config["singularity"]["diffparc"],
        topup_opt=get_eddy_topup_opt,
        s2v_opts=get_eddy_s2v_opts,
        slspec_opt=get_eddy_slspec_opt,
        eddy_cmd=get_eddy_cmd,
    output:
        #eddy creates many files, so write them to a eddy subfolder instead
        out_folder=directory(
            bids(root=root, suffix="eddy", datatype="dwi", **subj_wildcards)
        ),
        dwi=os.path.join(
            bids(root=root, suffix="eddy", datatype="dwi", **subj_wildcards),
            "dwi.nii.gz",
        ),
        bvec=os.path.join(
            bids(root=root, suffix="eddy", datatype="dwi", **subj_wildcards),
            "dwi.eddy_rotated_bvecs",
        ),
    threads: 16  #this needs to be set in order to avoid multiple gpus from executing
    resources:
        gpus=1,
        time=360,  #6 hours (this is a conservative estimate, may be shorter)
        mem_mb=32000,
    log:
        bids(root="logs", suffix="run_eddy.log", **subj_wildcards),
    group:
        "subj"
    shell:
        "{params.eddy_cmd} "
        " --imain={input.dwi_concat} --mask={input.brainmask} "
        " --acqp={input.phenc_concat} --index={input.eddy_index_txt} "
        " --bvecs={input.bvecs} --bvals={input.bvals} "
        " --out={params.out_prefix} "
        " {params.s2v_opts} "
        " {params.slspec_opt} "
        " {params.topup_opt} "
        " {params.flags}  &> {log}"


rule cp_eddy_outputs:
    input:
        #get nii.gz, bvec, and bval from eddy output
        dwi=os.path.join(
            bids(root=root, suffix="eddy", datatype="dwi", **subj_wildcards),
            "dwi.nii.gz",
        ),
        bvec=os.path.join(
            bids(root=root, suffix="eddy", datatype="dwi", **subj_wildcards),
            "dwi.eddy_rotated_bvecs",
        ),
        bval=bids(
            root=root,
            suffix="dwi.bval",
            desc="degibbs",
            datatype="dwi",
            **subj_wildcards
        ),
        json=bids(
            root=root,
            suffix="dwi.json",
            desc="degibbs",
            datatype="dwi",
            **subj_wildcards
        ),
    output:
        multiext(
            bids(
                root=root, suffix="dwi", desc="eddy", datatype="dwi", **subj_wildcards
            ),
            ".nii.gz",
            ".bvec",
            ".bval",
            ".json",
        ),
    group:
        "subj"
    run:
        for in_file, out_file in zip(input, output):
            shell("cp -v {in_file} {out_file}")


rule eddy_quad:
    input:
        unpack(get_eddy_slspec_input),
        phenc_concat=bids(
            root=root,
            suffix="phenc.txt",
            desc="degibbs",
            datatype="dwi",
            **subj_wildcards
        ),
        eddy_index_txt=bids(
            root=root,
            suffix="dwi.eddy_index.txt",
            desc="degibbs",
            datatype="dwi",
            **subj_wildcards
        ),
        brainmask=get_dwi_mask(),
        bvals=bids(
            root=root,
            suffix="dwi.bval",
            desc="degibbs",
            datatype="dwi",
            **subj_wildcards
        ),
        bvecs=bids(
            root=root,
            suffix="dwi.bvec",
            desc="degibbs",
            datatype="dwi",
            **subj_wildcards
        ),
        eddy_dir=bids(root=root, suffix="eddy", datatype="dwi", **subj_wildcards),
    params:
        eddy_prefix=lambda wildcards, input: os.path.join(input.eddy_dir, "dwi"),
        slspec_opt=get_eddy_slspec_opt,
    output:
        out_dir=directory(
            bids(root=root, suffix="eddy.qc", datatype="dwi", **subj_wildcards)
        ),
        eddy_qc_pdf=bids(
            root=root, suffix="eddy.qc/qc.pdf", datatype="dwi", **subj_wildcards
        ),
    container:
        config["singularity"]["diffparc"]  #fsl
    group:
        "subj"
    shell:
        "rmdir {output.out_dir} && "
        "eddy_quad {params.eddy_prefix} -idx {input.eddy_index_txt} -par {input.phenc_concat} "
        " -m {input.brainmask} -b {input.bvals} -g {input.bvecs} -o {output.out_dir} "
        " {params.slspec_opt} -v"


rule split_eddy_qc_report:
    input:
        eddy_qc_pdf=bids(
            root=root, suffix="eddy.qc/qc.pdf", datatype="dwi", **subj_wildcards
        ),
    output:
        report(
            directory(
                bids(
                    root=root, suffix="eddy.qc_pages", datatype="dwi", **subj_wildcards
                )
            ),
            patterns=["{pagenum}.png"],
            caption="../report/eddy_qc.rst",
            category="eddy_qc",
            subcategory=bids(
                **subj_wildcards, include_subject_dir=False, include_session_dir=False
            ),
        ),
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    script:
        "../scripts/split_pdf.py"


# TODO: gradient correction (optional step -- only if gradient file is provided)..
#  gradient_unwarp.py
#  reg_jacobian
#  convertwarp -> change this to wb_command -convert-warpfield  to get itk transforms
#  applywarp -> change this to antsApplyTransforms


rule eddymotion:
    input:
        dwi=bids(
            root=root,
            suffix="dwi.nii.gz",
            datatype="dwi",
            desc="degibbs",
            **subj_wildcards
        ),
        bvec=bids(
            root=root,
            suffix="dwi.bvec",
            datatype="dwi",
            desc="degibbs",
            **subj_wildcards
        ),
        bval=bids(
            root=root,
            suffix="dwi.bval",
            datatype="dwi",
            desc="degibbs",
            **subj_wildcards
        ),
        b0=bids(
            root=root,
            suffix="concatb0.nii.gz",
            datatype="dwi",
            desc="degibbs",
            **subj_wildcards
        ),
        brainmask=get_dwi_mask(),
    params:
        n_iter=1,
        model="b0",
    output:
        affine_dir=directory(
            bids(
                root=root,
                suffix="transforms",
                desc="eddymotion",
                datatype="dwi",
                **subj_wildcards
            )
        ),
    shadow:
        "minimal"
    threads: 32
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    script:
        "../scripts/eddymotion.py"


# then need a rule to apply the transforms
def get_preproc_dwi_input_dwispace(wildcards):
    if config["use_eddy"]:
        # use eddy dwi as preproc
        return multiext(
            bids(
                root=root, suffix="dwi", desc="eddy", datatype="dwi", **subj_wildcards
            ),
            ".nii.gz",
            ".bvec",
            ".bval",
            ".json",
        )

    else:
        # use moco dwi as preproc
        return multiext(
            bids(
                root=root, suffix="dwi", desc="moco", datatype="dwi", **subj_wildcards
            ),
            ".nii.gz",
            ".bvec",
            ".bval",
            ".json",
        )


rule cp_to_preproc_dwi:
    """ should use config flags to decide what input to use here.
    e.g. degibbs, moco, topup, eddy... """
    input:
        get_preproc_dwi_input_dwispace,
    output:
        expand(
            bids(
                root=root,
                suffix="dwi.{ext}",
                datatype="dwi",
                desc="preproc",
                **subj_wildcards
            ),
            ext=["nii.gz", "bvec", "bval", "json"],
            allow_missing=True,
        ),
    group:
        "subj"
    run:
        for in_file, out_file in zip(input, output):
            shell("cp -v {in_file} {out_file}")
