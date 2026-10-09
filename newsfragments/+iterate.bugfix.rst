``iter_records`` reads the records as they were when called, so a loop may
edit the block. Editing it mid-iteration used to end the iteration early.
