``SDBlock`` has dict-style record access: ``block[name]`` returns the last
value named ``name`` or raises ``KeyError``, ``block[name] = value`` is
``set_record``, ``del block[name]`` removes every such record or raises
``KeyError``, and ``for name in block`` lazily yields the names, duplicates
included, raising ``RuntimeError`` if the records change.
