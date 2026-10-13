"""Separating an SDF into molecule "SDBlock"s of chemical data and metadata."""

import sendoff.native as native

MalformedRecordError = native.MalformedRecordError
SDBlock = native.SDBlock
parse_sdf = native.parse_sdf
write_sdf = native.write_sdf
