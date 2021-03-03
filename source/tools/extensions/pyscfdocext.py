from sphinx.util.nodes import split_explicit_title
from docutils import nodes, utils

SOURCE_URI = 'https://github.com/pyscf/pyscf/tree/master/%s'

def source_role(typ, rawtext, text, lineno, inliner, options={}, content=[]):
    has_title, title, target = split_explicit_title(text)
    title = utils.unescape(title)
    target = utils.unescape(target)
    refnode = nodes.reference(title, title, refuri=SOURCE_URI % target)
    return [refnode], []

def setup(app):
    app.add_role('source', source_role)
    return {'version': '0.1', 'parallel_read_safe': True}
