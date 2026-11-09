"""
Write outputs to a temp file, then rename over the target.

This folder is on OneDrive, which can mark old outputs 'online-only' (dataless). Opening such
a file for writing makes macOS download it first, which can hang and raise
TimeoutError [Errno 60]. Renaming a fresh file over it replaces it without that download.
"""
import json
import os


def _replace(name, write):
    root, ext = os.path.splitext(name)
    tmp = f"{root}.tmp{ext}"
    write(tmp)
    os.replace(tmp, name)


def save_fig(fig, name, **kw):
    _replace(name, lambda p: fig.savefig(p, **kw))


def save_csv(df, name, **kw):
    _replace(name, lambda p: df.to_csv(p, **kw))


def save_json(obj, name, **kw):
    def write(p):
        with open(p, "w") as f:
            json.dump(obj, f, **kw)
    _replace(name, write)
