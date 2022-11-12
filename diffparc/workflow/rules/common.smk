def get_template_prefix(root, subj_wildcards, template):
    """creates prefix for template files, including subject/session wildcards
    so that DAGs for each subject/session are kept independent.
        e.g.: sub-001/tpl-MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym"""

    if "session" in subj_wildcards:
        return (
            "{root}/sub-{subject}/ses-{session}/tpl-{template}/tpl-{template}".format(
                root=root,
                subject=subj_wildcards["subject"],
                session=subj_wildcards["session"],
                template=template,
            )
        )
    else:
        return "{root}/sub-{subject}/tpl-{template}/tpl-{template}".format(
            root=root, subject=subj_wildcards["subject"], template=template
        )
