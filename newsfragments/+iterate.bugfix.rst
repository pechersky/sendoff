``iter_records`` raises ``RuntimeError`` if the block's records change during
iteration, like a dict; it used to end early. Loop over ``records()`` to edit
them.
